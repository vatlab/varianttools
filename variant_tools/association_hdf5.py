#!/usr/bin/env python
#
# $File: association.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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

import sys, os, re
from multiprocessing import Process, Queue, Pipe, Value, Array, Lock
import time
import random
import math
from copy import copy, deepcopy
from .project import Project, Field, AnnoDB, AnnoDBWriter, MaintenanceProcess
from .utils import ProgressBar, consolidateFieldName, DatabaseEngine, delayedAction, \
     env, executeUntilSucceed, ShelfDB, safeMapFloat, PrettyPrinter, flatten, hasGenoInfo
from .phenotype import Sample
from .tester import *
from .rtester import RTest, SKAT

from variant_tools.vt_sqlite3 import OperationalError
import argparse

from .HDF5_accessor import *
import HDF5_storage as storage
import sqlite3
import glob



class GroupHDFGenerator(Process):

        def __init__(self,fileQueue,project,group_names,job):
            Process.__init__(self)
            self.fileQueue=fileQueue
            self.proj=project
            self.group_names=group_names
            self.proc_index=job

        def run(self):
            
            while True:
                HDFfileName=self.fileQueue.get()
                if HDFfileName is None:
                    break
                hdf5_file=tb.open_file(HDFfileName,mode="r")
               
                HDFfileGroupName=HDFfileName.replace(".h5","_multi_genes.h5")
                if (os.path.isfile(HDFfileGroupName)):
                    os.remove(HDFfileGroupName)
                try:
                    rownames=hdf5_file.root.rownames[:]
                    rownames=[int(x.decode("utf-8")) for x in rownames]
                    db = DatabaseEngine()
                    db.connect('{0}.proj'.format(self.proj.name), '__asso_tmp')
                    cur = db.cursor()
            
                    select_group="SELECT {0}, group_concat(variant_id) from __asso_tmp group by {0}".\
                       format(self.group_names[0])
                    for row in cur.execute(select_group):
                        geneSymbol=transformGeneName(row[0])
                        ids=row[1].split(",")
                        # idPos=[int(x) for x in ids]
                        # idPos.sort()     
                        # varPos=[np.where(rownames==varID) for varID in ids]
                        # idPos=[rownames.index(varID) for varID in ids]
                        # idPos.sort()
                        ids=[int(x) for x in ids]
                        ids.sort()

                        minPos=rownames.index(ids[0])
                        maxPos=rownames.index(ids[-1])
                        # if maxPos<minPos:
                        #     temp=minPos
                        #     minPos=maxPos
                        #     maxPos=temp
                        idPos=rownames[minPos:maxPos+1]
                        idPos=[int(x) for x in idPos]
                        idPos.sort()
 
                        # try:
                        #     print(idPos[0]-1,idPos[-1]+1)
                        # except:
                        #     print(ids,idPos,minPos,maxPos,idPos)
                        #     continue
                        sub_indptr=hdf5_file.root.indptr[idPos[0]-1:idPos[-1]+1]
                        sub_indices=hdf5_file.root.indices[min(sub_indptr):max(sub_indptr)]
                        sub_data=hdf5_file.root.data[min(sub_indptr):max(sub_indptr)]
                        sub_indptr=[sub_indptr[i]-sub_indptr[0] for i in range(len(sub_indptr))]
                        # sub_shape=(len(idPos),hdf5_file.root.shape[1])
                        sub_shape=(len(sub_indptr)-1,hdf5_file.root.shape[1])
                        adjust_id=[ idPos[i]-idPos[0] for i in range(len(idPos))]
                        
                        sub_matrix=csr_matrix((sub_data,sub_indices,sub_indptr),shape=sub_shape)    
                        geneGenotypes=sub_matrix[adjust_id,]
                        
                        storage.store_csr_genotype_into_one_HDF5(geneGenotypes,"dataTable_"+geneSymbol,ids,self.proc_index,HDFfileGroupName)        
                        # storage.store_csr_arrays_into_earray_HDF5(sub_data,sub_indices,sub_indptr,sub_shape,ids,"dataTable_"+geneSymbol,HDFfileGroupName) 
                    hdf5_file.close()

                except KeyboardInterrupt as e:
                    hdf5_file.close()
                    pass

def transformGeneName(geneSymbol):
    if ("-" in geneSymbol):
        geneSymbol=geneSymbol.replace("-","_")
    pattern=re.compile(r'\.')
    if pattern.findall(geneSymbol):
        geneSymbol=geneSymbol.replace(".","_")
    return geneSymbol

 
def generateHDFbyGroup(testManager):
        # HDFfileName=self.proj.name+"_genotype.h5"
        # HDFfileGroupName=self.proj.name+"_genotype_multi_genes.h5"

        HDFfileNames=glob.glob("tmp*_genotypes.h5")
        njobs=8
        groupGenerators=[]
        fileQueue=Queue()
        for HDFfileName in HDFfileNames:
            fileQueue.put(HDFfileName)
        for i in range(njobs):
            groupGenerator=GroupHDFGenerator(fileQueue,testManager.proj,testManager.group_names,i)
            groupGenerator.start()
            groupGenerators.append(groupGenerator)
            fileQueue.put(None)
        for groupHDFGenerator in groupGenerators:
            groupHDFGenerator.join()


def getGenotype_HDF5(worker, group):
    '''Get genotype for variants in specified group'''
    # get variant_id

    where_clause = ' AND '.join(['{0}={1}'.format(x, worker.db.PH) for x in worker.group_names])
    cur = worker.db.cursor()
    # variant info
    var_info, variant_id = worker.getVarInfo(group, where_clause)

    # get genotypes / genotype info
    genotype = []
    geno_info = {x:[] for x in worker.geno_info}
    # getting samples locally from my own connection

    geneSymbol=transformGeneName(group[0])

    files=glob.glob("tmp*_genotypes_multi_genes.h5")
    files=sorted(files, key=lambda name: int(name.split("_")[1]))
    for fileName in files:
        startSample=int(fileName.split("_")[1])
        endSample=int(fileName.split("_")[2])
        hdf5db=HDF5Engine_multi("dataTable")
        hdf5db.connect_HDF5(fileName)  
        hdf5db.load_HDF5(geneSymbol)
        snpdict=hdf5db.load_genotype_by_variantID()
        # dbID=self.cached_genes[geneSymbol]-1

        # for ID in self.sample_IDs:
        for ID in range(endSample-startSample+1):
            data=snpdict[ID+1]
            # handle missing values
            gtmp = [data.get(x, [worker.g_na] + [float('NaN')]*len(worker.geno_info)) for x in variant_id]
            # handle -1 coding (double heterozygotes)
            genotype.append([2.0 if x[0] == -1.0 else x[0] for x in gtmp])
            #
            # handle genotype_info
            #
            for idx, key in enumerate(worker.geno_info):
                geno_info[key].append([x[idx+1] if (type(x[idx+1]) in [int, float]) else float('NaN') for x in gtmp])
    gname = ':'.join(list(map(str, group)))
    return worker.filterGenotype(genotype, geno_info, var_info, gname)


                    
# def runAssocaition_HDF5(args,asso,proj,results):
	
#     try:
#         nJobs = max(min(args.jobs, len(asso.groups)), 1)
#         nLoaders = env.associate_num_of_readers
#         # if no env.associate_num_of_readers is set we limit it to a max of 8.
#         if not nLoaders > 0:
#             nLoaders = min(8, nJobs)

#         ready_flags = Array('L', [0]*nLoaders)

#         cached_samples = Array('L', max(asso.sample_IDs) + 1)

     
#         maintenance_flag = Value('L', 1)
#         maintenance = MaintenanceProcess(proj, {'genotype_index': asso.sample_IDs}, maintenance_flag)
#         maintenance.start()

#         grpQueue = Queue()
#         # the result queue is used by workers to return results
#         resQueue = Queue()
#         # see if all workers are ready
#         ready_flags = Array('L', [0]*nJobs)

#         for j in range(nJobs):
#             # the dictionary has the number of temporary database for each sample
#             AssoTestsWorker(asso, grpQueue, resQueue, ready_flags, j, 
#                 {x:y-1 for x,y in enumerate(cached_samples) if y > 0},
#                 results.fields, shelf_lock).start()


#         for grp in asso.groups:
#             grpQueue.put(grp)
#         # the worker will stop once all jobs are finished
#         for j in range(nJobs):
#             grpQueue.put(None)
#         #
#         count = 0
#         s = delayedAction(env.logger.info, "Starting {} association test workers".format(nJobs))
#         while True:
#             if all(ready_flags):
#                 break
#             else:
#                 time.sleep(random.random()*2)
#         del s
#         prog = ProgressBar('Testing for association', len(asso.groups))
#         try:
#             while True:
#                 # if everything is done
#                 if count >= len(asso.groups):
#                     break
#                 # not done? wait from the queue and write to the result recorder
#                 res = resQueue.get()
#                 results.record(res)
#                 # update progress bar
#                 count = results.completed()
#                 prog.update(count, results.failed())
#                 # env.logger.debug('Processed: {}/{}'.format(count, len(asso.groups)))
#         except KeyboardInterrupt as e:
#             env.logger.error('\nAssociation tests stopped by keyboard interruption ({}/{} completed).'.\
#                               format(count, len(asso.groups)))
#             results.done()
#             proj.close()
#             sys.exit(1)
#         # finished
#         prog.done()
#         results.done()
#         # summary
#         env.logger.info('Association tests on {} groups have completed. {} failed.'.\
#                          format(results.completed(), results.failed()))
#         # use the result database in the project
#         if args.to_db:
#             proj.useAnnoDB(AnnoDB(proj, args.to_db, ['chr', 'pos'] if not args.group_by else args.group_by))
#         # tells the maintenance process to stop
#         maintenance_flag.value = 0
#         # wait for the maitenance process to stop
#         s = delayedAction(env.logger.info,
#                           "Maintaining database. This might take a few minutes.", delay=10)
#         maintenance.join()
#         del s
#     except Exception as e:
#         env.logger.error(e)
#         sys.exit(1)
