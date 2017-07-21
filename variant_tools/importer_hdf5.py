#!/usr/bin/env python
#
# $File: importer.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import os
import sys
import re
import array
import time
from heapq import heappush, heappop, heappushpop
from multiprocessing import Process, Pipe, Value, Lock, Manager,Array
from multiprocessing import Queue as mpQueue
import Queue
if sys.version_info.major == 2:
    from itertools import izip, repeat
else:
    izip = zip
    from itertools import repeat
from collections import defaultdict
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, \
    openFile, DatabaseEngine, hasCommand, \
    downloadFile, env, RefGenome

import numpy as np
from scipy.sparse import csr_matrix,csc_matrix,rand,coo_matrix,hstack
import tables as tb
import HDF5_storage as storage
import glob
from shutil import copyfile


try:
    from variant_tools.cgatools import normalize_variant
except ImportError as e:
    sys.exit('Failed to import module ({})\n'
        'Please verify if you have installed variant tools successfully (using command '
        '"python setup.py install")'.format(e))


# preprocessors
from .preprocessor import *
#



class HDF5GenotypeImportWorker(Process):
    '''This class starts a process, import genotype to a temporary genotype database.'''
    def __init__(self, processor,readQueue,variantIndex, start_sample,end_sample,variant_count, 
        proc_index,geno_info,HDFfileName):
        '''
        variantIndex: a dictionary that returns ID for each variant.
        filelist: files from which variantIndex is created. If the passed filename
            is not in this list, this worker will suicide so that it can be replaced 
            by a worker with more variants.
        encoding, genotypefield, genotype_info, ranges:  parameters to import data
        geno_count:  a shared variable to report number of genotype imported
        status:      an ImportStatus object to monitor the progress
        '''
        Process.__init__(self, name='GenotypeImporter')
        # self.daemon=True
        self.processor=processor
        self.readQueue=readQueue
        # self.indptr_count=indptr_count
        self.variantIndex = variantIndex
        self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample=start_sample
        self.end_sample=end_sample
        self.geno_info=geno_info
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]
        self.colnames=[sample for sample in range(start_sample+1,end_sample+1)]

        self.HDFfileName=HDFfileName
        self.info={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[],[],[],[],[]]
        


    def writeIntoFile(self,chr):
        # print(self.proc_index,self.numVariants,max(self.indptr),min(self.indptr),max(self.indices),min(self.indices))
        # print(self.proc_index,len(self.indptr),self.variant_count.value,max(self.indptr),min(self.indptr))    
        # ids=[x+1 for x in range(len(self.indptr))]    
        # self.indptr_count.value=self.indptr[-1]
        self.variant_count.value = self.variant_count.value + len(self.indptr)
        shape=(self.variant_count.value,self.end_sample-self.start_sample)
        
        #If file doesn't exist, create file
        #if file exists but chromosome does not exist, create chr group
        #if file exists and chromosome exists, append to file
 
        if not os.path.isfile(self.HDFfileName):
            self.indptr=[0]+self.indptr
            storage.store_csr_arrays_into_earray_HDF5(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames,"",chr,self.HDFfileName) 
            if len(self.geno_info)>0:
                for key,value in self.info.items():
                    value[0]=[0]+value[0]
                    # if (key!="AD_geno" and key!="PL_geno"):
                    storage.store_csr_arrays_into_earray_HDF5(value[2],value[1],value[0],shape,value[3],self.colnames,key,chr,self.HDFfileName) 

        else:
            node="/chr"+chr
            f=tb.open_file(self.HDFfileName,"r")
            if node in f:
                group=f.get_node(node)
                currentStart=group.indptr[-1]
                f.close()
                self.indptr=[x+currentStart for x in self.indptr]
                storage.append_csr_arrays_into_earray_HDF5(self.data,self.indices,self.indptr,shape,self.rownames,"",chr,self.HDFfileName) 
                
                if len(self.geno_info)>0:
                    for key,value in self.info.items():
                        f=tb.open_file(self.HDFfileName,"r")
                        group=f.get_node(node+"/"+key)
                        currentStart=group.indptr[-1]
                        f.close()
                        value[0]=[x+currentStart for x in value[0]]
                        storage.append_csr_arrays_into_earray_HDF5(value[2],value[1],value[0],shape,value[3],key,chr,self.HDFfileName) 

            # else:
            #     self.indptr=[0]+self.indptr
            #     f.close()
            #     storage.store_csr_arrays_into_earray_HDF5(self.data,self.indices,self.indptr,shape,self.rownames,"",chr,self.HDFfileName) 
            
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]

        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[0],[],[],[],[]]



    def run(self): 
        genoCount=0
        genoDict={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                genoDict[info.name]=0
        pre_variant_ID=None
        pre_chr=None
        self.start_count = self.variant_count.value
        firstLine=True

        # if self.HDFfileName.split("_")[1]=="1":
        #      print(self.HDFfileName)
        self.processor.reset(import_sample_range=[self.start_sample,self.end_sample])
        while True:
            line=self.readQueue.get()
      
            if line is None:
                self.writeIntoFile(pre_chr)
                break
            for bins,rec in self.processor.process(line):
                # if self.HDFfileName.split("_")[1]=="1":
                #      print(rec)
                variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                if pre_variant_ID==variant_id:
                    continue
                else:
                    pre_variant_ID=variant_id

                if firstLine:
                    pre_chr=rec[0]
                    firstLine=False
                
                if pre_chr!=rec[0]:
                    self.writeIntoFile(pre_chr)
                    pre_chr=rec[0]

                numSamples=self.end_sample-self.start_sample
         
                for idx in range(numSamples):
                    try:
                        pos=4+idx
                        if rec[pos] is not None:
                            if rec[pos]!='0':
                                genoCount=genoCount+1
                                self.indices.append(idx)
                                self.data.append(int(rec[pos]))
                        else:
                            genoCount=genoCount+1
                            self.indices.append(idx)
                            self.data.append(np.nan)
                    except IndexError:
                        env.logger.warning('Incorrect number of genotype fields: {} fields found, {} expected for record {}'.format(
                            len(rec), fld_cols[-1][-1] + 1, rec))
                self.indptr.append(genoCount)
                self.rownames.append(variant_id)


                if len(self.geno_info)>0:
                    for infoIndex,info in enumerate(self.geno_info):
    
                        for idx in range(numSamples):
                            try:
                                pos=4+(infoIndex+1)*numSamples+idx
                                if rec[pos] is not None:
                                    if rec[pos]!='0':
                                        genoDict[info.name]=genoDict[info.name]+1
                                        self.info[info.name][1].append(idx)
                                        try:
                                            self.info[info.name][2].append(int(rec[pos]))
                                        except ValueError:
                                            if rec[pos] != ".":
                                                self.info[info.name][2].append(rec[pos])
                                            else:
                                                self.info[info.name][2].append(0)

                                else:
                                    genoDict[info.name]=genoDict[info.name]+1
                                    self.info[info.name][1].append(idx)
                                    self.info[info.name][2].append(np.nan)
                            except IndexError:
                                env.logger.warning('Incorrect number of genotype fields: {} fields found, {} expected for record {}'.format(
                                    len(rec), fld_cols[-1][-1] + 1, rec))
                        self.info[info.name][0].append(genoDict[info.name])
                        self.info[info.name][3].append(variant_id)








def updateSample(importer,start_sample,end_sample,sample_ids,names,allNames,HDF5fileName):
    cur=importer.db.cursor()
    pos=0
    for id in range(start_sample,end_sample):
        sql="UPDATE sample SET HDF5=?, COLPOS=? WHERE sample_id=? and sample_name=?"
        # task=(HDF5fileName,allNames[names[id]],sample_ids[id],names[id])
        task=(HDF5fileName,pos,sample_ids[id],names[id])
        cur.execute(sql,task)
        pos+=1
    


def manageHDF5(importer):
    # for f in glob.glob("*h5"):
    #     os.remove(f)
    cur=importer.db.cursor()
    sql="ALTER TABLE sample ADD COLUMN HDF5 CHAR(25)"
    try:
        cur.execute(sql)
    except:
        pass
    sql="ALTER TABLE sample ADD COLUMN COLPOS INT"
    try:
        cur.execute(sql)
    except:
        pass
    sql="SELECT sample_name,COLPOS from sample"
    allNames={}
    maxCount=0
    for rec in cur.execute(sql):
        allNames[rec[0]]=int(rec[1])
        if int(rec[1])>maxCount:
            maxCount=int(rec[1])
    return allNames,maxCount+1








def importGenotypesInParallel(importer,num_sample=0):

    # number of genotypes each process have imported
    
    allNames,sample_count=manageHDF5(importer)
    indptr_count=None
    for count, input_filename in enumerate(importer.files):
        # indptr_count=[Value('L', 0) for x in range(importer.jobs)]
        
        env.logger.info('{} variants from {} ({}/{})'.format('Importing', input_filename, count + 1, len(importer.files)))
        importer.importVariant(input_filename)
        env.logger.info('{:,} new variants {}{}{} from {:,} lines are imported.'\
            .format(importer.count[2], "(" if importer.count[2] else '', 
                ', '.join(['{:,} {}'.format(x, y) for x, y in \
                    zip(importer.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'unsupported']) if x > 0]),
                    ")" if importer.count[2] else '', importer.count[0]))
        # genotypes?
        if importer.genotype_field:
            importer.prober.reset()
        # if there are samples?
        sample_ids, genotype_status,names = importer.getSampleIDs(input_filename)

        for name in names:
            if name not in allNames:
                allNames[name]=sample_count
                sample_count+=1
        workload=None
       

        if num_sample>0:
            num_files=int(len(sample_ids)/num_sample)
            workload=[num_sample]*num_files
            if len(sample_ids)%num_sample!=0:
                workload.append(len(sample_ids)%num_sample)
            # if len(workload)<importer.jobs:
            #     workload=workload+[0]*(importer.jobs-len(workload))
        else:
            if len(sample_ids)<importer.jobs:
                workload=[len(sample_ids)]+[0]*(importer.jobs-1)
            else:
                workload = [int(float(len(sample_ids)) / importer.jobs)] * importer.jobs

            # if there are missing ones, spread it across workers ...
            # less than 0 is possible because of the at least 10 policy
            unallocated = max(0, len(sample_ids) - sum(workload))
            for i in range(unallocated):
                workload[i % importer.jobs] += 1


        print(workload)
        numProcess=len(workload)

        # importers = [None] * numProcess
        importers=[None]*importer.jobs
        variant_import_count = [Value('L', -1) for x in range(numProcess)]

        # we should have file line count from importVariant
        num_of_lines = importer.count[0] 

        for i in range(len(importer.count)):
            importer.total_count[i] += importer.count[i]
            importer.count[i] = 0
        if len(sample_ids) == 0:
            continue

        # env.logger.debug('Workload of processes: {}'.format(workload))
        
        start_sample =0
        line_no = 0
        num_lines=0    
        readQueue=[]

        taskQueue=Queue.Queue()

        input_prefix=os.path.basename(input_filename).replace(".vcf","")
        for job in range(numProcess):      
            if workload[job] == 0:
                continue
            readQueue.append(mpQueue())
            end_sample = min(start_sample + workload[job], len(sample_ids))
            print(start_sample,end_sample)     
            if end_sample <= start_sample:
                continue
            HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
            updateSample(importer,start_sample,end_sample,sample_ids,names,allNames,HDFfile_Merge)
            taskQueue.put(HDF5GenotypeImportWorker(importer.processor,readQueue[job], importer.variantIndex, start_sample, end_sample, 
                        variant_import_count[job], job, importer.genotype_info,HDFfile_Merge))
            start_sample = end_sample

        total_genotype_count=importer.total_count[2]*len(sample_ids)
        prog = ProgressBar('Importing genotypes', total_genotype_count,initCount=0)
  
        
        task=None
        
        with openFile(input_filename) as input_file:
            for line in input_file:
                line = line.decode(importer.encoding)
                if line.startswith('#') or "<" in line:
                    continue
                line_no+=1
                num_lines+=1
                for job in range(len(readQueue)):
                    readQueue[job].put(line)
                if (line_no>=2000):
                    print(num_lines)
                    prog.update(num_lines*len(sample_ids))
                    line_no=0
                    for job in range(numProcess):
                        readQueue[job].put(None)
                    while taskQueue.qsize()>0:
                        for i in range(importer.jobs):    
                            if importers[i] is None or not importers[i].is_alive():
                                print(i)
                                task=taskQueue.get()
                                importers[i]=task
                                importers[i].start()           
                                break 
                    for worker in importers:
                        worker.join()
                    start_sample =0
                    for job in range(numProcess):
                        if workload[job] == 0:
                            continue
                        end_sample = min(start_sample + workload[job], len(sample_ids))
                        if end_sample <= start_sample:
                            continue
                        HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
                        taskQueue.put(HDF5GenotypeImportWorker(importer.processor,readQueue[job], importer.variantIndex, start_sample, end_sample, 
                            variant_import_count[job], job, importer.genotype_info,HDFfile_Merge))
                        start_sample = end_sample   


        for job in range(len(readQueue)):
            readQueue[job].put(None)
        while taskQueue.qsize()>0:
            for i in range(importer.jobs):    
                if importers[i] is None or not importers[i].is_alive():
                    print(i)
                    task=taskQueue.get()
                    importers[i]=task
                    importers[i].start()           
                    break 
        for worker in importers:
            worker.join() 
        prog.done()
   