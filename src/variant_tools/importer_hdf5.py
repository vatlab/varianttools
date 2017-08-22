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
#import HDF5_storage as storage
import glob
from shutil import copyfile
from .accessor import *

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
    '''This class starts a process, import genotype to a temporary HDF5 file.
        Args

            processor: a processor to parse the row 
            readQueue: a Queue of vcf rows
            variantIndex: a dictionary that returns ID for each variant.
            start_sample: the sample id of sample in the first column
            end_sample: the sample id of sample in the last column
            sample_ids: a list of sample IDS
            variant_count: count of variants
            proc_index: the index of process
            geno_info: genotype info other than GT
            dbLocation: the HDF5 file name 

    '''
    def __init__(self, processor,readQueue,variantIndex, start_sample,end_sample,sample_ids,variant_count, 
        proc_index,geno_info,dbLocation):

        Process.__init__(self, name='GenotypeImporter')
        # self.daemon=True
        self.processor=processor
        self.readQueue=readQueue
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
        self.colnames=[sample_ids[i] for i in range(start_sample,end_sample)]

        self.dbLocation=dbLocation

        self.info={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[],[],[],[],[]]
        

    def writeIntoFile(self,chr):

        self.variant_count.value = self.variant_count.value + len(self.indptr)
        shape=(self.variant_count.value,self.end_sample-self.start_sample)

        storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
        # make a HMatrix object which is a matrix with rownames and colnames
        hmatrix=HMatrix(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames)
        # write GT into file
        storageEngine.store(hmatrix,chr)
        # write geno info into HDF5 if exists
        if len(self.geno_info)>0:
                for key,value in self.info.items():
                    hmatrix=HMatrix(value[2],value[1],value[0],shape,value[3],self.colnames)
                    storageEngine.store(hmatrix,chr,key) 
        
        # clean up
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]

        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[],[],[],[],[]]



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

        # if self.dbLocation.split("_")[1]=="1":
        #      print(self.dbLocation)
        self.processor.reset(import_sample_range=[self.start_sample,self.end_sample])
        while True:
            line=self.readQueue.get()
            if line is None:
                self.writeIntoFile(pre_chr)
                break
            for bins,rec in self.processor.process(line):
          
                variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                if pre_variant_ID==variant_id:
                    continue
                else:
                    pre_variant_ID=variant_id
                #deal with variant of multiple chromosome in one vcf file
                if firstLine:
                    pre_chr=rec[0]
                    firstLine=False
                if pre_chr!=rec[0]:
                    self.writeIntoFile(pre_chr)
                    pre_chr=rec[0]


                numSamples=self.end_sample-self.start_sample
                
                #Process GT        
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

                #Process other geno_info
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
    
    for id in range(start_sample,end_sample):
        sql="UPDATE sample SET HDF5=? WHERE sample_id=? and sample_name=?"
        # task=(HDF5fileName,allNames[names[id]],sample_ids[id],names[id])
        task=(HDF5fileName,sample_ids[id],names[id])
        cur.execute(sql,task)
        
    


def manageHDF5(importer,allNames={}):

    cur=importer.db.cursor()
    sql="ALTER TABLE sample ADD COLUMN HDF5 CHAR(25)"
    try:
        cur.execute(sql)
    except:
        pass
   
    sql="SELECT sample_name,sample_id from sample"

  
    for rec in cur.execute(sql):
        if rec[0] not in allNames:
            allNames[rec[0]]=int(rec[1])
     
    return allNames




def importGenotypesInParallel(importer,num_sample=0):

    allNames=manageHDF5(importer)
    
    for count, input_filename in enumerate(importer.files):
        
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
        if len(sample_ids) == 0:
            continue
        allNames=manageHDF5(importer,allNames)
        sample_ids=[int(allNames[name]) for name in names]
       
        workload=None
       
        #determine number of samples to be processed in each process
        if num_sample>0:
            #if number of samples in each HDF5 is determined
            num_files=int(len(sample_ids)/num_sample)
            workload=[num_sample]*num_files
            if len(sample_ids)%num_sample!=0:
                workload.append(len(sample_ids)%num_sample)
            if len(workload)<importer.jobs:
                importer.jobs=len(workload)
        else:
            if len(sample_ids)<50:
                importer.jobs=1
                workload=[len(sample_ids)]
            else:
                workload = [int(float(len(sample_ids)) / importer.jobs)] * importer.jobs
          
            unallocated = max(0, len(sample_ids) - sum(workload))
            for i in range(unallocated):
                workload[i % importer.jobs] += 1


        env.logger.debug("work load {}".format(workload))
        numTasks=len(workload)
        # importers = [None] * numProcess
        importers=[None]*importer.jobs
        variant_import_count = [Value('L', 0) for x in range(numTasks)]


        for i in range(len(importer.count)):
            importer.total_count[i] += importer.count[i]
            importer.count[i] = 0
        

        # env.logger.debug('Workload of processes: {}'.format(workload))
        start_sample =0
        line_no = 0
        num_lines=0    
        readQueue=[]

        taskQueue=Queue.Queue()
        input_prefix=os.path.basename(input_filename).replace(".vcf","")
        
        #Put tasks in the queue first
        for job in range(numTasks):      
            if workload[job] == 0:
                continue
            readQueue.append(mpQueue())
            end_sample = min(start_sample + workload[job], len(sample_ids))
            env.logger.debug("{},{}".format(start_sample,end_sample))     
            if end_sample <= start_sample:
                continue
            HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
            updateSample(importer,start_sample,end_sample,sample_ids,names,allNames,HDFfile_Merge)
            taskQueue.put(HDF5GenotypeImportWorker(importer.processor,readQueue[job], importer.variantIndex, start_sample, end_sample, 
                        sample_ids,variant_import_count[job], job, importer.genotype_info,HDFfile_Merge))
            start_sample = end_sample

        total_genotype_count=importer.total_count[2]*len(sample_ids)
        prog = ProgressBar('Importing genotypes', total_genotype_count,initCount=0)
  
        
        task=None
        # This size of cache
        chunckOfLines=20000
        with openFile(input_filename) as input_file:
            for line in input_file:
                line = line.decode(importer.encoding)
                if line.startswith('#') or "<" in line:
                    continue
                line_no+=1
                num_lines+=1
                for job in range(numTasks):
                    readQueue[job].put(line)
                #Read the VCF by chunck
                if (line_no>chunckOfLines):
                    env.logger.debug("number of lines {}".format(num_lines))
                    line_no=0
                    for job in range(numTasks):
                        readQueue[job].put(None)
                    while taskQueue.qsize()>0:
                        for i in range(importer.jobs):    
                            if importers[i] is None or not importers[i].is_alive():
                         
                                task=taskQueue.get()
                                importers[i]=task
                                importers[i].start()           
                                break 
                    for worker in importers:
                        worker.join()
                    prog.update(num_lines*len(sample_ids))
                    start_sample =0
                    for job in range(numTasks):
                        if workload[job] == 0:
                            continue
                        end_sample = min(start_sample + workload[job], len(sample_ids))
                        if end_sample <= start_sample:
                            continue
                        HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
                        taskQueue.put(HDF5GenotypeImportWorker(importer.processor,readQueue[job], importer.variantIndex, start_sample, end_sample, 
                            sample_ids,variant_import_count[job], job, importer.genotype_info,HDFfile_Merge))
                        start_sample = end_sample   

        #Process remaining lines
        for job in range(len(readQueue)):
            readQueue[job].put(None)
        while taskQueue.qsize()>0:
            for i in range(importer.jobs):    
                if importers[i] is None or not importers[i].is_alive():
                    task=taskQueue.get()
                    importers[i]=task
                    importers[i].start()           
                    break 
        for worker in importers:
            worker.join() 
        prog.done()
   
