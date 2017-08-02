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
from multiprocessing import Process, Queue, Pipe, Value, Array, Lock, Manager
from multiprocessing import Queue as mpQueue
import Queue
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


# Maybe useful for single VCF file?
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

                HDFfileGroupName=HDFfileName.replace(".h5","_multi_genes.h5")
                if (os.path.isfile(HDFfileGroupName)):
                    os.remove(HDFfileGroupName)

                hdf5=HDF5Engine_access(HDFfileName)
                hdf5group=HDF5Engine_storage(HDFfileGroupName)
                
                
                try:
                    
                    db = DatabaseEngine()
                    db.connect('{0}.proj'.format(self.proj.name), '__asso_tmp')
                    cur = db.cursor()
            
                    # select_group="SELECT {0}, group_concat(variant_id) from __asso_tmp group by {0}".\
                    #    format(self.group_names[0])
                    select_group="SELECT {0}, group_concat(t.variant_id) from __asso_tmp as t join variant as v on t.variant_id=v.variant_id group by {0} order by v.pos".\
                       format(self.group_names[0])
                    for row in cur.execute(select_group):
                        geneSymbol=transformGeneName(row[0])
                        ids=row[1].split(",")
                        ids=[int(x) for x in ids]
                        # ids.sort()
                        chr= getChr(ids[0],db.cursor())
                        # if (self.proc_index==1):
                        #     print(ids)
                        #     print(hdf5.get_rownames(chr)[:])
                        sub_data,sub_indices,sub_indptr,sub_shape,rownames,colnames=hdf5.get_geno_info_by_variant_IDs(ids,chr)
                        if sub_indices is not None:
                            hdf5group.store_arrays_into_HDF5(sub_data,sub_indices,sub_indptr,sub_shape,rownames,colnames,chr,geneSymbol) 
                    hdf5.close()
                    hdf5group.close()

                except KeyboardInterrupt as e:
                    hdf5.close()
                    hdf5group.close()
                    pass




def generateHDFbyGroup(testManager,njobs):
        # HDFfileName=self.proj.name+"_genotype.h5"
        # HDFfileGroupName=self.proj.name+"_genotype_multi_genes.h5"

        HDFfileNames=glob.glob("tmp*_genotypes.h5")
        groupGenerators=[]
        fileQueue=mpQueue()
        taskQueue=Queue.Queue()
        start=time.time()
        groupGenerators=[None]*njobs
        for HDFfileName in HDFfileNames:
            fileQueue.put(HDFfileName)
        for i in range(len(HDFfileNames)):

            taskQueue.put(GroupHDFGenerator(fileQueue,testManager.proj,testManager.group_names,i))
        while taskQueue.qsize()>0:
            for i in range(njobs):
                if groupGenerators[i] is None or not groupGenerators[i].is_alive():
                    task=taskQueue.get()
                    groupGenerators[i]=task
                    task.start()
                    fileQueue.put(None)
                    break
        for groupHDFGenerator in groupGenerators:
            groupHDFGenerator.join()
        print("group time: ",time.time()-start)


#keep sparase matrix in memory and write onece
class GroupHDFGenerator_memory(Process):

        def __init__(self,geneDict,geneSet,fileQueue,project,group_names,job):
            Process.__init__(self)
            self.geneDict=geneDict
            self.geneSet=geneSet
            self.fileQueue=fileQueue
            self.proj=project
            self.group_names=group_names
            self.proc_index=job

        def run(self):
            
            while True:
                genoDict={}
                HDFfileName=self.fileQueue.get()
                if HDFfileName is None:
                    break
                # hdf5_file=tb.open_file(HDFfileName,mode="r")
                hdf5=HDF5Engine_access(HDFfileName)

                HDFfileGroupName=HDFfileName.replace(".h5","_multi_genes.h5")
                if (os.path.isfile(HDFfileGroupName)):
                    os.remove(HDFfileGroupName)
                hdf5group=HDF5Engine_storage(HDFfileGroupName)
                try:
                    chr="22"
                    rownames=hdf5.get_rownames(chr)
                    colnames=hdf5.get_colnames(chr)
          
                    for idx,id in enumerate(rownames):
                        try:
                            geneNames=self.geneDict[id]
                            for geneName in geneNames:
                                if geneName not in genoDict:
                                    #indptr,indices,data,rownames
                                    genoDict[geneName]=[[0],[],[],[]]
            
                                variant_ID,indices,data=hdf5.get_geno_info_by_row_pos(idx,chr)
                              
                                lastPos=genoDict[geneName][0][-1]
                                if len(indices)==0:
                                    genoDict[geneName][0].append(lastPos)
                                else:
                                    genoDict[geneName][0].append(lastPos+len(indices))
                                    genoDict[geneName][1].extend(indices)
                                    genoDict[geneName][2].extend(data)
                                genoDict[geneName][3].append(variant_ID)
                              
                        except:
                            pass
                    for key,value in genoDict.items():
                        hdf5group.store_arrays_into_HDF5(value[2],value[1],value[0],(len(value[3]),len(colnames)),value[3],colnames,chr,key) 
                    hdf5.close()
                    hdf5group.close()
                except KeyboardInterrupt as e:
                    hdf5.close()
                    hdf5group.close()
                    pass


#read ids one by one and append to gene group HDF5
class GroupHDFGenerator_append(Process):

        def __init__(self,geneDict,geneSet,fileQueue,project,group_names,job):
            Process.__init__(self)
            self.geneDict=geneDict
            self.geneSet=geneSet
            self.fileQueue=fileQueue
            self.proj=project
            self.group_names=group_names
            self.proc_index=job

        def run(self):
            
            while True:
                genoDict={}
                HDFfileName=self.fileQueue.get()
                if HDFfileName is None:
                    break
      
                hdf5=HDF5Engine_access(HDFfileName)
                HDFfileGroupName=HDFfileName.replace(".h5","_multi_genes.h5")
                if (os.path.isfile(HDFfileGroupName)):
                    os.remove(HDFfileGroupName)
                hdf5group=HDF5Engine_storage(HDFfileGroupName)
                chr="22"
                colnames=hdf5.get_colnames(chr)   
                rownames=hdf5.get_rownames(chr)
     
                for geneName in self.geneSet:
                    hdf5group.store_arrays_into_HDF5([],[],[0],(0,0),[],colnames,chr,geneName)
                try:

                    for idx,id in enumerate(rownames):
                        try:
                            geneNames=self.geneDict[id]                           
                            for geneName in geneNames:
      
                                variant_ID,indices,data=hdf5.get_geno_info_by_row_pos(idx,chr)
                                indptr=None
                                if len(indices)>0:
                                    indptr=[len(indices)]
                                shape=(1,len(colnames))
                                hdf5group.append_arrays_into_HDF5(data,indices,indptr,shape,[id],chr,geneName)
                        except KeyError: 
                            pass
                    hdf5.close()
                    hdf5group.close()
                except KeyboardInterrupt as e:
                    hdf5.close()
                    hdf5group.close()
                    pass



def getGroupDict(testManager):
    db = DatabaseEngine()
    db.connect('{0}.proj'.format(testManager.proj.name), '__asso_tmp')
    cur = db.cursor()
    select_group="SELECT {0}, group_concat(variant_id) from __asso_tmp group by {0}".\
       format(testManager.group_names[0])
    geneDict={}
    geneSet=set()
    for row in cur.execute(select_group):
        geneSymbol=transformGeneName(row[0])
        geneSet.add(geneSymbol)
        ids=row[1].split(",")
        ids=[int(id) for id in ids]
        for id in ids:
            if not id in geneDict:
                geneDict[id]=[]
            geneDict[id].append(geneSymbol)
    shareDict=Manager().dict()
    for key,value in geneDict.items():
        shareDict[key]=value

    return shareDict,geneSet



def generateHDFbyGroup_update(testManager,njobs):
        # HDFfileName=self.proj.name+"_genotype.h5"
        # HDFfileGroupName=self.proj.name+"_genotype_multi_genes.h5"

        HDFfileNames=glob.glob("tmp*_genotypes.h5")

        fileQueue=mpQueue()
        taskQueue=Queue.Queue()
        start=time.time()
        groupGenerators=[None]*njobs
        geneDict,geneSet=getGroupDict(testManager)
        for HDFfileName in HDFfileNames:
            fileQueue.put(HDFfileName)
        for i in range(len(HDFfileNames)):
            taskQueue.put(GroupHDFGenerator_memory(geneDict,geneSet,fileQueue,testManager.proj,testManager.group_names,i))
            # taskQueue.put(GroupHDFGenerator_append(geneDict,geneSet,fileQueue,testManager.proj,testManager.group_names,i))
        while taskQueue.qsize()>0:
            for i in range(njobs):
                if groupGenerators[i] is None or not groupGenerators[i].is_alive():
                    task=taskQueue.get()
                    groupGenerators[i]=task
                    fileQueue.put(None)
                    task.start()
                    break
        for groupHDFGenerator in groupGenerators:
            groupHDFGenerator.join()
        print("group time: ",time.time()-start)




def transformGeneName(geneSymbol):
    if ("-" in geneSymbol):
        geneSymbol=geneSymbol.replace("-","_")
    pattern=re.compile(r'\.')
    if pattern.findall(geneSymbol):
        geneSymbol=geneSymbol.replace(".","_")
    return geneSymbol

 


def getChr(variantID,cur):
    find_chr="SELECT chr from variant where variant_id={0}".format(variantID)
    chr= [rec[0] for rec in cur.execute(find_chr)]
    return chr[0]


def getGenotype_HDF5(worker, group):
    '''Get genotype for variants in specified group'''
    # get variant_id

    where_clause = ' AND '.join(['{0}={1}'.format(x, worker.db.PH) for x in worker.group_names])
    cur = worker.db.cursor()
    # variant info
    var_info, variant_id = worker.getVarInfo(group, where_clause)
    
    find_chr="SELECT chr from variant where variant_id={0}".format(variant_id[0])
    chr=getChr(variant_id[0],cur)
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
        
       
        hdf5db=HDF5Engine_access(fileName)
        snpdict=hdf5db.get_geno_info_by_group(geneSymbol,chr)
        # print(snpdict.keys(),startSample,endSample)
        hdf5db.close()
        for ID in range(startSample,endSample+1):
            data=snpdict[ID]
            
            # handle missing values

        # hdf5db=HDF5Engine_multi(fileName)
        # hdf5db.load_HDF5(geneSymbol,chr)
        # snpdict=hdf5db.load_all_GT()
        # for ID in range(startSample,endSample+1):           
        #     data=snpdict[ID]

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

