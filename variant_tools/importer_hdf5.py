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
from multiprocessing import Process, Pipe, Value, Lock, Manager,Queue,Array
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

class mergeHDF5(Process):

    def __init__(self,HDFfileNames,HDFMergeFileName):
        Process.__init__(self)
        self.HDFfileNames=HDFfileNames
        self.HDFMergeFileName=HDFMergeFileName

    def run(self):
        first=True
        prev_end=None
        for HDFfileName in self.HDFfileNames:
            if first:
                if os.path.isfile(self.HDFMergeFileName):
                    os.remove(self.HDFMergeFileName)
                copyfile(HDFfileName,self.HDFMergeFileName)
                first=False
            else:
                f=tb.open_file(self.HDFMergeFileName,"r")
                currentStart=f.root.indptr[-1]
                currentRows=f.root.shape[0]
                currentCols=f.root.shape[1]
                f.close()
                f=tb.open_file(HDFfileName,"r")
                indptr=f.root.indptr

                # new_indptr=[currentStart+indptr[i]-indptr[0] for i in range(1,len(indptr))]
                new_indptr=indptr[1:]
                new_shape=(currentRows+len(new_indptr),currentCols)
                # print(HDFfileName)
                # print(currentStart,new_indptr[0],new_indptr[-1])
                # print(new_indptr,f.root.shape[:])
                storage.append_csr_arrays_into_earray_HDF5(f.root.data,f.root.indices,new_indptr,new_shape,f.root.rownames,self.HDFMergeFileName) 
                f.close()
            os.remove(HDFfileName)



class HDF5GenotypeImportWorker(Process):
    '''This class starts a process, import genotype to a temporary genotype database.'''
    def __init__(self, processor,readQueue,indptr_count,variantIndex, filelist,start_sample,end_sample,ranges, variant_count, 
        proc_index,first,HDFfileName):
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
        self.indptr_count=indptr_count
        self.variantIndex = variantIndex
        self.filelist = filelist
        self.ranges = ranges
        self.variant_count = variant_count
        
        self.proc_index = proc_index
        self.start_sample=start_sample
        self.end_sample=end_sample
        
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.ids=[]
        self.first=first
        self.HDFfileName=HDFfileName
        



    def writeIntoFile(self):
        # print(self.proc_index,self.numVariants,max(self.indptr),min(self.indptr),max(self.indices),min(self.indices))
        # print(self.proc_index,len(self.indptr),self.variant_count.value,max(self.indptr),min(self.indptr))    
        # ids=[x+1 for x in range(len(self.indptr))]    
        self.indptr_count.value=self.indptr[-1]
        self.variant_count.value = self.variant_count.value + len(self.indptr)
        shape=(self.variant_count.value,self.end_sample-self.start_sample)
        
        if self.first:
            if os.path.isfile(self.HDFfileName):
                os.remove(self.HDFfileName)
            storage.store_csr_arrays_into_earray_HDF5(self.data,self.indices,self.indptr,shape,self.ids,"",self.HDFfileName) 
            self.first=False
        else:
            storage.append_csr_arrays_into_earray_HDF5(self.data,self.indices,self.indptr,shape,self.ids,self.HDFfileName) 
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.ids=[]



    def run(self): 
        genoCount=self.indptr_count.value
        if self.first:
            self.indptr.append(0)
        pre_variant_ID=None
        self.start_count = self.variant_count.value
        line_no=0
        self.processor.reset(import_sample_range=[self.start_sample,self.end_sample])
        while True:
            line=self.readQueue.get()
            line_no+=1
            if line is None:
                self.writeIntoFile()
                break
            for bins,rec in self.processor.process(line):
                variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                if (pre_variant_ID==variant_id):
                    continue
                else:
                    pre_variant_ID=variant_id
                   
                for idx in range(self.end_sample-self.start_sample):
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
                self.ids.append(variant_id)



def importGenotypesInParallel(importer):
    '''import files in parallel, by importing variants and genotypes separately, and in their own processes. 
    More specifically, suppose that there are three files

    file1: variant1, sample1.1, sample1.2, sample1.3
    file2: variant2, sample2.1, sample2.2
    file3: variant3, sample3.1, sample3.2

    where variant1, 2, and 3 are three potentially overlapping sets of variants.
    sample1, sample2, sample3 are three groups of samples that are divided by the number of samples
    and number of processes (self.jobs) (say, 3000 samples divided into three groups of 1000 samples).

    '''
    importers = [None] * importer.jobs
    # number of genotypes each process have imported
    variant_import_count = [Value('L', -1) for x in range(importer.jobs)]
    indptr_count=[Value('L', 0) for x in range(importer.jobs)]
   

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
        sample_ids, genotype_status = importer.getSampleIDs(input_filename)
        # we should have file line count from importVariant
        num_of_lines = importer.count[0] 
        print(num_of_lines)

        for i in range(len(importer.count)):
            importer.total_count[i] += importer.count[i]
            importer.count[i] = 0

        if len(sample_ids) == 0:
            continue

        # if len(sample_ids) > 2000:
        #     workload = [max(10, int((1.2*len(sample_ids)/numTempFile)*(1-x/(3.*(numTempFile - 1))))) for x in range(numTempFile)]
        # else:
        workload = [max(10, int(float(len(sample_ids)) / importer.jobs))] * importer.jobs
        # if there are missing ones, spread it across workers ...
        # less than 0 is possible because of the at least 10 policy
        unallocated = max(0, len(sample_ids) - sum(workload))
        for i in range(unallocated):
            workload[i % importer.jobs] += 1
        # 
        # env.logger.debug('Workload of processes: {}'.format(workload))
        
        start_sample = 0
        line_no = 0
        num_lines=0
        first=True
        # indptrDict=Manager().dict({ (i,0)for i in range(self.jobs)})
        readQueue=[]
        inputers=[]

        input_prefix=os.path.basename(input_filename).replace(".vcf","")
        for job in range(importer.jobs):
            readQueue.append(Queue())
            if workload[job] == 0:
                continue
            end_sample = min(start_sample + workload[job], len(sample_ids))
     
            if end_sample <= start_sample:
                continue
            for i in range(importer.jobs):
                if importers[i] is None or not importers[i].is_alive():
                    HDFfile="tmp_"+str(start_sample+1)+"_"+str(end_sample)+"_"+input_prefix+"_genotypes.h5"
                    importers[i] = HDF5GenotypeImportWorker(importer.processor,readQueue[i], indptr_count[i], importer.variantIndex, importer.files[:count+1], start_sample, end_sample, importer.ranges,
                        variant_import_count[i], i, first,HDFfile)
                    importers[i].start()
                    start_sample = end_sample
                    break 
        total_genotype_count=importer.total_count[2]*len(sample_ids)
        prog = ProgressBar('Importing genotypes', total_genotype_count,initCount=0)
        
        with openFile(input_filename) as input_file:
            for line in input_file:
                line = line.decode(importer.encoding)
                if line.startswith('#') or "<" in line:
                    continue
                line_no+=1
                num_lines+=1
                for job in range(importer.jobs):
                    readQueue[job].put(line)
                if (line_no>=10000):
                    prog.update(num_lines*len(sample_ids))
                    first=False
                    line_no=0
                    for job in range(importer.jobs):
                        readQueue[job].put(None)
                    for importer in importers:
                        importer.join()
                    start_sample = 0
                    for job in range(importer.jobs):
                        if workload[job] == 0:
                            continue
                        end_sample = min(start_sample + workload[job], len(sample_ids))
                        if end_sample <= start_sample:
                            continue
                        for i in range(importer.jobs):
                            if importers[i] is None or not importers[i].is_alive():
                                HDFfile="tmp_"+str(start_sample+1)+"_"+str(end_sample)+"_"+input_prefix+"_genotypes.h5"
                                importers[i] = HDF5GenotypeImportWorker(importer.processor,readQueue[i], indptr_count[i], importer.variantIndex, importer.files[:count+1], start_sample, end_sample, importer.ranges,
                                    variant_import_count[i], i,first,HDFfile)
                                importers[i].start()
                                start_sample = end_sample
                                break      

        for job in range(importer.jobs):
            readQueue[job].put(None)

        for worker in importers:
            worker.join() 
        prog.done()
        # add range of id to the file name 
        HDFfileNames=glob.glob("tmp*"+input_prefix+"_genotypes.h5")
        for HDFfileName in HDFfileNames:
            HDFfile=tb.open_file(HDFfileName)
            firstID=HDFfile.root.rownames[0]
            lastID=HDFfile.root.rownames[-1]
            HDFfile.close()
            cols=HDFfileName.split("_")
            newHDFfileName="_".join(cols[:3])+"_"+str(int(firstID))+"_"+str(int(lastID))+"_"+"_".join(cols[3:])
            if os.path.isfile(newHDFfileName):
                os.remove(newHDFfileName)
            os.rename(HDFfileName,newHDFfileName)

    mergeJobs=[]
    if len(importer.files)>1:
        start_sample=0

        for job in range(importer.jobs):
            if workload[job] == 0:
                continue
            end_sample = min(start_sample + workload[job], len(sample_ids))
            if end_sample <= start_sample:
                continue
            HDFfileNames=glob.glob("tmp_"+str(start_sample+1)+"_"+str(end_sample)+"*_genotypes.h5")
            HDFfileNames=sorted(HDFfileNames,key=lambda x : int(os.path.basename(x).split("_")[3]))
            HDFMergeFileName="tmp_"+str(start_sample+1)+"_"+str(end_sample)+"_genotypes.h5"
            p=mergeHDF5(HDFfileNames,HDFMergeFileName)
            p.start()
            mergeJobs.append(p)
            
            start_sample = end_sample
    for job in mergeJobs:
        job.join()