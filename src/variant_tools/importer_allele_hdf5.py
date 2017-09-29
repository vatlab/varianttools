import os
import sys
import re
import array
import time
from heapq import heappush, heappop, heappushpop
from multiprocessing import Process, Pipe, Value, Lock, Manager,Array
from multiprocessing import Queue as mpQueue
import queue
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
from .importer_vcf_to_hdf5 import *
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
    def __init__(self, chunk,variantIndex, start_sample,end_sample,sample_ids,variant_count, 
        proc_index,geno_info,dbLocation,build):

        Process.__init__(self, name='GenotypeImporter')
        # self.daemon=True
        self.chunk=chunk
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
        self.sample_ids=sample_ids
        self.colnames=[sample_ids[i] for i in range(start_sample,end_sample)]
        self.genoCount=0
        self.dbLocation=dbLocation
        self.build=build

        self.info={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[],[],[],[],[]]
        

    # def writeIntoFile(self,chr):

    #     self.variant_count.value = self.variant_count.value + len(self.indptr)
    #     shape=(self.variant_count.value,self.end_sample-self.start_sample)

    #     storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
    #     # make a HMatrix object which is a matrix with rownames and colnames
    #     hmatrix=HMatrix(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames)
    #     # write GT into file
    #     storageEngine.store(hmatrix,chr)
    #     # write geno info into HDF5 if exists
    #     if len(self.geno_info)>0:
    #             for key,value in list(self.info.items()):
    #                 hmatrix=HMatrix(value[2],value[1],value[0],shape,value[3],self.colnames)
    #                 storageEngine.store(hmatrix,chr,key) 
        
    #     # clean up
    #     self.indptr=[]
    #     self.indices=[]
    #     self.data=[]
    #     self.rownames=[]

    #     if len(self.geno_info)>0:
    #         for info in self.geno_info:
    #             #indptr,indices,data,shape,rownames
    #             self.info[info.name]=[[],[],[],[],[]]

    def getGT(self,variant_id,GT,altIndex):
        for idx in range(self.start_sample,self.end_sample):
            if GT[idx] is not None:
                if altIndex==0:
                    if GT[idx]!=0:
                        if GT[idx]!=3 and GT[idx]!=4:
                            self.genoCount=self.genoCount+1
                            self.indices.append(idx-self.start_sample)
                            self.data.append(GT[idx])
                        else:
                            self.genoCount=self.genoCount+1
                            self.indices.append(idx-self.start_sample)
                            self.data.append(np.nan)  
                elif altIndex==1:
                        self.genoCount=self.genoCount+1
                        self.indices.append(idx-self.start_sample)   
                        if GT[idx]==3:
                            self.data.append(1)
                        elif GT[idx]==4:
                            self.data.append(2)
                        else:
                            self.data.append(np.nan)
            else:
                self.genoCount=self.genoCount+1
                self.indices.append(idx-self.start_sample)
                self.data.append(np.nan)
        self.indptr.append(self.genoCount)
        self.rownames.append(variant_id)




    def writeIntoHDF(self):
        storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
        chr=self.chunk["variants/CHROM"][0]
        for i in range(len(self.chunk["variants/ID"])):
            chr=self.chunk["variants/CHROM"][i]
            ref=self.chunk["variants/REF"][i]
            pos=self.chunk["variants/POS"][i]
            GT=self.chunk["calldata/GT"][i].tolist()
        
            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt=self.chunk["variants/ALT"][i][altIndex]
         
                if alt!="":
                    if tuple((chr, ref, alt)) in self.variantIndex:
                        variant_id  = self.variantIndex[tuple((chr, ref, alt))][pos][0]
                        self.getGT(variant_id,GT,altIndex)
                    else:
                        rec=[str(chr),str(pos),ref,alt]  
                        msg=normalize_variant(RefGenome(self.build).crr, rec, 0, 1, 2, 3)
                        if tuple((rec[0], rec[2], rec[3])) in self.variantIndex:
                            variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                            self.getGT(variant_id,GT,altIndex)
         
        shape=(len(self.indptr),len(self.colnames))
        
        # make a HMatrix object which is a matrix with rownames and colnames
        hmatrix=HMatrix(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames)
        # write GT into file
       
        storageEngine.store(hmatrix,"22")
        storageEngine.close()
      


    def run(self): 
        genoCount=0
        genoDict={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                genoDict[info.name]=0

        self.start_count = self.variant_count.value

        # chunk=self.readQueue.get()
        self.writeIntoHDF()

      








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

        line_no = 0
        num_lines=0    
        readQueue=[]

        taskQueue=queue.Queue()
        task=None


        compression='gzip'
        compression_opts=1
        shuffle=False
        overwrite=False
        vlen=True
        fields=None
        types=None
        numbers=None
        alt_number=DEFAULT_ALT_NUMBER
        fills=None
        region=None
        tabix='tabix'
        samples=None
        transformers=None
        buffer_size=DEFAULT_BUFFER_SIZE
        chunk_length=DEFAULT_CHUNK_LENGTH
        chunk_width=DEFAULT_CHUNK_WIDTH
        log=None
       
        _, samples, headers, it = iter_vcf_chunks(
                input_filename, fields=fields, types=types, numbers=numbers, alt_number=alt_number,
                buffer_size=buffer_size, chunk_length=chunk_length, fills=fills, region=region,
                tabix=tabix, samples=samples, transformers=transformers
            )

        #Put tasks in the queue first
        for job in range(numTasks):      
            if workload[job] == 0:
                continue
            # readQueue.append(mpQueue())

        prog = ProgressBar('Importing genotypes', importer.total_count[2],initCount=0)
        
        for chunk, _, _, _ in it:                
            start_sample =0
            for job in range(numTasks):
                # readQueue[job].put(chunk)
                if workload[job] == 0:
                    continue
                end_sample = min(start_sample + workload[job], len(sample_ids))
                if end_sample <= start_sample:
                    continue
                HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
              
                taskQueue.put(HDF5GenotypeImportWorker(chunk, importer.variantIndex, start_sample, end_sample, 
                    sample_ids,variant_import_count[job], job, importer.genotype_info,HDFfile_Merge,importer.build))
                start_sample = end_sample   
            while taskQueue.qsize()>0:
                for i in range(importer.jobs):    
                    if importers[i] is None or not importers[i].is_alive():     
                        task=taskQueue.get()
                        importers[i]=task
                        importers[i].start()           
                        break 
            for worker in importers:
                worker.join()

            prog.update(num_lines)



        prog.done()
   
