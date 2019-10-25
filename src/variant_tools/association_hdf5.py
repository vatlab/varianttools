#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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

import glob
import os
import queue
import re
import time
from multiprocessing import Manager, Process
from multiprocessing import Queue as mpQueue

from .accessor import Engine_Access, Engine_Storage
from .utils import DatabaseEngine


class GroupHDFGenerator(Process):
    """ Assume vcf files are sorted and grouped by chromosome positions.
        This class generates temporary HDF5 file which has variants grouped by chromosome and groupName(for example, gene name)
        Each input HDF5 file is queried to get groups of variants. The resut is written into output file.
    """

    def __init__(self, fileQueue, project, group_names, job):
        Process.__init__(self)
        self.fileQueue = fileQueue
        self.proj = project
        self.group_names = group_names
        self.proc_index = job

    def run(self):

        while True:
            HDFfileName = self.fileQueue.get()
            if HDFfileName is None:
                break

            HDFfileGroupName = HDFfileName.replace(".h5", "_multi_genes.h5")
            # if (os.path.isfile(HDFfileGroupName)):
            #     os.remove(HDFfileGroupName)

            accessEngine = Engine_Access.choose_access_engine(HDFfileName)
            storageEngine = Engine_Storage.choose_storage_engine(
                HDFfileGroupName)

            try:

                db = DatabaseEngine()
                db.connect('{0}.proj'.format(self.proj.name), '__asso_tmp')
                time.sleep(self.proc_index)
                cur = db.cursor()

                # select_group="SELECT {0}, group_concat(variant_id) from __asso_tmp group by {0}".\
                #    format(self.group_names[0])
                if self.group_names[0] == "variant_chr" and self.group_names[
                        1] == "variant_pos":
                    select_group = "select t.variant_chr,t.variant_pos,t.variant_id from __asso_tmp as t join variant as v on t.variant_id=v.variant_id  order by v.pos"
                    for row in cur.execute(select_group):
                        chr = row[0]
                        pos = row[1]
                        id = row[2]
                        if not storageEngine.checkGroup(chr, "pos" + str(pos)):
                            updated_rownames, colnames, subMatrix = accessEngine.get_genotype(
                                [id], "", [chr])
                            if subMatrix is not None:
                                storageEngine.store(subMatrix, chr,
                                                    "pos" + str(pos) + "/GT")
                                storageEngine.store(
                                    updated_rownames, chr,
                                    "pos" + str(pos) + "/rownames")
                                if not storageEngine.checkGroup(
                                        chr, "colnames"):
                                    storageEngine.store(colnames, chr,
                                                        "colnames")

                else:
                    select_group="SELECT {0}, group_concat(variant_id) from (select {0},t.variant_id from __asso_tmp as t join variant as v on t.variant_id=v.variant_id  order by v.pos) group by {0}".\
                       format(self.group_names[0])
                    cur.execute(
                        'SELECT value FROM project WHERE name="multiVCF";')
                    multiVCF = cur.fetchall()
                    if len(multiVCF) == 0:
                        multiVCF = 0
                    else:
                        multiVCF = int(multiVCF[0][0])
                    cur.execute(
                        'SELECT DISTINCT(file_id) from sample where HDF5="{0}"'
                        .format(HDFfileName))
                    _ = cur.fetchall()
                    for row in cur.execute(select_group):
                        geneSymbol = transformGeneName(row[0])
                        ids = row[1].split(",")
                        ids = [int(x) for x in ids]
                        # ids.sort()
                        chr = getChr(ids[0], db.cursor())
                        if not storageEngine.checkGroup(chr, geneSymbol):
                            # if int(multiVCF[0][0])==0 or (int(multiVCF[0][0])==1 and int(file_id[0][0])==1):
                            varDict = getChrs(ids, db.cursor())
                            if multiVCF == 0 and len(varDict.keys()) == 1:
                                # updated_rownames,colnames,subMatrix=accessEngine.get_geno_by_variant_IDs(ids,chr)

                                updated_rownames, colnames, subMatrix = accessEngine.get_genotype(
                                    ids, "", [chr])

                                if subMatrix is not None:
                                    storageEngine.store(subMatrix, chr,
                                                        geneSymbol + "/GT")
                                    storageEngine.store(
                                        updated_rownames, chr,
                                        geneSymbol + "/rownames")
                                    if not storageEngine.checkGroup(
                                            chr, "colnames"):
                                        storageEngine.store(
                                            colnames, chr, "colnames")
                            else:
                                for chr, vids in varDict.items():
                                    try:
                                        updated_rownames, colnames, subMatrix = accessEngine.get_geno_by_sep_variant_ids(
                                            vids, chr)
                                        if subMatrix is not None:
                                            storageEngine.store(
                                                subMatrix, chr,
                                                geneSymbol + "/GT")
                                            storageEngine.store(
                                                updated_rownames, chr,
                                                geneSymbol + "/rownames")
                                            if not storageEngine.checkGroup(
                                                    chr, "colnames"):
                                                storageEngine.store(
                                                    colnames, chr, "colnames")
                                        else:
                                            storageEngine.store(
                                                updated_rownames, chr,
                                                geneSymbol + "/rownames")
                                            if not storageEngine.checkGroup(
                                                    chr, "colnames"):
                                                storageEngine.store(
                                                    colnames, chr, "colnames")
                                    except TypeError:
                                        pass

                    # storageEngine.close()
                accessEngine.close()
                storageEngine.close()

            except KeyboardInterrupt:
                accessEngine.close()
                storageEngine.close()
                pass


def generateHDFbyGroup(testManager, njobs):
    """This class manages the processes for generating the temp HDF5 files with variants grouped by gene name.
    """

    HDFfileNames = glob.glob("tmp*_genotypes.h5")
    groupGenerators = []
    fileQueue = mpQueue()
    taskQueue = queue.Queue()

    groupGenerators = [None] * min(njobs, len(HDFfileNames))

    if len(testManager.sample_IDs) > 0:
        HDFfileNames = []
        cur = testManager.db.cursor()
        sampleIDstring = ",".join([str(ID) for ID in testManager.sample_IDs])
        sql = "select distinct HDF5 from sample where sample_id in ({0});".format(
            sampleIDstring)
        cur.execute(sql)
        res = cur.fetchall()
        for filename in res:
            HDFfileNames.append(filename[0])

    for HDFfileName in HDFfileNames:
        fileQueue.put(HDFfileName)
    for i in range(len(HDFfileNames)):

        taskQueue.put(
            GroupHDFGenerator(fileQueue, testManager.proj,
                              testManager.group_names, i))
    while taskQueue.qsize() > 0:
        for i in range(njobs):
            if groupGenerators[i] is None or not groupGenerators[i].is_alive():
                task = taskQueue.get()
                groupGenerators[i] = task
                task.start()
                fileQueue.put(None)
                break
    for groupHDFGenerator in groupGenerators:
        groupHDFGenerator.join()
    # print(("group time: ",time.time()-start))


class GroupHDFGenerator_memory(Process):
    """Assume variant IDs are not ordered by chromosome position. This class keep the group of variants as cache and write together into HDF5.
    """

    def __init__(self, geneDict, geneSet, fileQueue, project, group_names, job):
        Process.__init__(self)
        self.geneDict = geneDict
        self.geneSet = geneSet
        self.fileQueue = fileQueue
        self.proj = project
        self.group_names = group_names
        self.proc_index = job

    def run(self):

        while True:
            genoDict = {}
            HDFfileName = self.fileQueue.get()
            if HDFfileName is None:
                break
            # hdf5_file=tb.open_file(HDFfileName,mode="r")
            accessEngine = Engine_Access.choose_access_engine(HDFfileName)

            HDFfileGroupName = HDFfileName.replace(".h5", "_multi_genes.h5")
            if (os.path.isfile(HDFfileGroupName)):
                os.remove(HDFfileGroupName)
            storageEngine = Engine_Storage.choose_storage_engine(
                HDFfileGroupName)
            try:
                chr = "22"
                rownames = accessEngine.get_rownames(chr)
                colnames = accessEngine.get_colnames(chr)

                for idx, id in enumerate(rownames):
                    try:
                        geneNames = self.geneDict[id]
                        for geneName in geneNames:
                            if geneName not in genoDict:
                                #indptr,indices,data,rownames
                                genoDict[geneName] = [[0], [], [], []]

                            variant_ID, indices, data = accessEngine.get_geno_info_by_row_pos(
                                idx, chr)

                            lastPos = genoDict[geneName][0][-1]
                            if len(indices) == 0:
                                genoDict[geneName][0].append(lastPos)
                            else:
                                genoDict[geneName][0].append(lastPos +
                                                             len(indices))
                                genoDict[geneName][1].extend(indices)
                                genoDict[geneName][2].extend(data)
                            genoDict[geneName][3].append(variant_ID)

                    except:
                        pass
                for key, value in list(genoDict.items()):
                    hdf5matrix = HMatrix(value[2], value[1], value[0],
                                         (len(value[3]), len(colnames)),
                                         value[3], colnames)
                    storageEngine.store(hdf5matrix, chr, key)
                accessEngine.close()
                storageEngine.close()
            except KeyboardInterrupt:
                accessEngine.close()
                storageEngine.close()
                pass


class GroupHDFGenerator_append(Process):
    """Assume variant IDs are not ordered by chromosome position. This class reads variants one by one and appends to the specified group in HDF5.
    """

    def __init__(self, geneDict, geneSet, fileQueue, project, group_names, job):
        Process.__init__(self)
        self.geneDict = geneDict
        self.geneSet = geneSet
        self.fileQueue = fileQueue
        self.proj = project
        self.group_names = group_names
        self.proc_index = job

    def run(self):

        while True:
            HDFfileName = self.fileQueue.get()
            if HDFfileName is None:
                break

            accessEngine = Engine_Access.choose_access_engine(HDFfileName)
            HDFfileGroupName = HDFfileName.replace(".h5", "_multi_genes.h5")
            if (os.path.isfile(HDFfileGroupName)):
                os.remove(HDFfileGroupName)
            storageEngine = Engine_Storage.choose_storage_engine(
                HDFfileGroupName)
            chr = "22"
            colnames = accessEngine.get_colnames(chr)
            rownames = accessEngine.get_rownames(chr)

            for geneName in self.geneSet:
                hdf5matrix = HMatrix([], [], [0], (0, 0), [], colnames)
                storageEngine.store(hdf5matrix, chr, geneName)
            try:

                for idx, id in enumerate(rownames):
                    try:
                        geneNames = self.geneDict[id]
                        for geneName in geneNames:

                            variant_ID, indices, data = accessEngine.get_geno_info_by_row_pos(
                                idx, chr)
                            indptr = None
                            if len(indices) > 0:
                                indptr = [len(indices)]
                            shape = (1, len(colnames))
                            hdf5matrix = HMatrix(data, indices, indptr, shape,
                                                 [id], colnames)
                            storageEngine.store(hdf5matrix, chr, geneName)
                    except KeyError:
                        pass
                accessEngine.close()
                storageEngine.close()
            except KeyboardInterrupt:
                accessEngine.close()
                storageEngine.close()
                pass


def getGroupDict(testManager):
    db = DatabaseEngine()
    db.connect('{0}.proj'.format(testManager.proj.name), '__asso_tmp')
    cur = db.cursor()
    select_group="SELECT {0}, group_concat(variant_id) from __asso_tmp group by {0}".\
       format(testManager.group_names[0])
    geneDict = {}
    geneSet = set()
    for row in cur.execute(select_group):
        geneSymbol = transformGeneName(row[0])
        geneSet.add(geneSymbol)
        ids = row[1].split(",")
        ids = [int(id) for id in ids]
        for id in ids:
            if not id in geneDict:
                geneDict[id] = []
            geneDict[id].append(geneSymbol)
    shareDict = Manager().dict()
    for key, value in list(geneDict.items()):
        shareDict[key] = value

    return shareDict, geneSet


def generateHDFbyGroup_update(testManager, njobs):

    HDFfileNames = glob.glob("tmp*_genotypes.h5")

    fileQueue = mpQueue()
    taskQueue = queue.Queue()
    groupGenerators = [None] * min(njobs, len(HDFfileNames))
    geneDict, geneSet = getGroupDict(testManager)
    for HDFfileName in HDFfileNames:
        fileQueue.put(HDFfileName)
    for i in range(len(HDFfileNames)):
        # taskQueue.put(GroupHDFGenerator_memory(geneDict,geneSet,fileQueue,testManager.proj,testManager.group_names,i))
        taskQueue.put(
            GroupHDFGenerator_append(geneDict, geneSet, fileQueue,
                                     testManager.proj, testManager.group_names,
                                     i))
    while taskQueue.qsize() > 0:
        for i in range(njobs):
            if groupGenerators[i] is None or not groupGenerators[i].is_alive():
                task = taskQueue.get()
                groupGenerators[i] = task
                fileQueue.put(None)
                task.start()
                break
    for groupHDFGenerator in groupGenerators:
        groupHDFGenerator.join()
    # print(("group time: ",time.time()-start))


def transformGeneName(geneSymbol):
    geneSymbol = str(geneSymbol)
    if ("-" in geneSymbol):
        geneSymbol = geneSymbol.replace("-", "_")
    pattern = re.compile(r'\.')
    if pattern.findall(geneSymbol):
        geneSymbol = geneSymbol.replace(".", "_")
    return geneSymbol


def getChr(variantID, cur):
    find_chr = "SELECT chr from variant where variant_id={0}".format(variantID)
    chr = [rec[0] for rec in cur.execute(find_chr)]
    return chr[0]


def getChrs(variantIDs, cur):
    idString = "(" + ",".join([str(variantID) for variantID in variantIDs
                              ]) + ")"
    find_chr = "SELECT chr,variant_id from variant where variant_id in " + idString
    varDict = {}
    for rec in cur.execute(find_chr):
        if rec[0] not in varDict:
            varDict[rec[0]] = []
        varDict[rec[0]].append(rec[1])
    return varDict


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3


def getGenotype_HDF5(worker, group, sample_IDs):
    """This function gets genotype of variants in specified group.
    """
    where_clause = ' AND '.join(
        ['{0}={1}'.format(x, worker.db.PH) for x in worker.group_names])
    cur = worker.db.cursor()
    # variant info
    var_info, variant_ids = worker.getVarInfo(group, where_clause)
    chr = getChr(variant_ids[0], cur)
    chrEnd = getChr(variant_ids[-1], cur)
    if chr != chrEnd:
        varDict = getChrs(variant_ids, cur)
        chrs = list(varDict.keys())
    else:
        varDict = {chr: variant_ids}
        chrs = [chr]

    # get genotypes / genotype info
    genotype = []
    geno_info = {x: [] for x in worker.geno_info}
    # getting samples locally from my own connection
    geneSymbol = transformGeneName(group[0])
    if len(group) == 2:
        geneSymbol = "pos" + str(group[1])
    HDFfileNames = glob.glob("tmp*_genotypes_multi_genes.h5")
    if len(worker.sample_names) > 0:
        HDFfileNames = []
        cur = worker.db.cursor()
        sampleNamestring = ",".join(
            ['"' + str(ID) + '"' for ID in worker.sample_names])
        sql = "select distinct HDF5 from sample where sample_name in ({0});".format(
            sampleNamestring)
        cur.execute(sql)
        res = cur.fetchall()
        for filename in res:
            HDFfileNames.append(filename[0].replace(".h5", "_multi_genes.h5"))

    HDFfileNames = sorted(
        HDFfileNames, key=lambda name: int(name.split("_")[1]))
    for fileName in HDFfileNames:
        accessEngine = Engine_Access.choose_access_engine(fileName)
        colnames = accessEngine.get_colnames(chr)
        colnames = intersection(sample_IDs, colnames.tolist())

        if len(chrs) == 1:
            snpdict = accessEngine.get_geno_by_group(chr, geneSymbol)
            accessEngine.close()

            for ID in colnames:
                data = snpdict[ID]
                gtmp = [
                    data.get(x, [worker.g_na] +
                             [float('NaN')] * len(worker.geno_info))
                    for x in varDict[chr]
                ]

                # handle -1 coding (double heterozygotes)
                genotype.append([2.0 if x[0] == -1.0 else x[0] for x in gtmp])
                #
                # handle genotype_info
                #
                for idx, key in enumerate(worker.geno_info):
                    geno_info[key].append([
                        x[idx + 1] if
                        (type(x[idx + 1]) in [int, float]) else float('NaN')
                        for x in gtmp
                    ])
        else:
            alldict = {}
            for chr in chrs:
                alldict[chr] = accessEngine.get_geno_by_group(chr, geneSymbol)
            for ID in colnames:
                gtmp = []
                for chr in chrs:
                    # snpdict=accessEngine.get_geno_by_group(chr,geneSymbol)
                    snpdict = alldict[chr]
                    data = snpdict[ID]
                    # gtmp.extend([data.get(x, [worker.g_na] + [float('NaN')]*len(worker.geno_info)) for x in varDict[chr]])
                    gtmp.extend([data.get(x) for x in varDict[chr]])
                    # handle -1 coding (double heterozygotes)
                genotype.append([2.0 if x[0] == -1.0 else x[0] for x in gtmp])

                #
                # handle genotype_info
                #
                for idx, key in enumerate(worker.geno_info):
                    geno_info[key].append([
                        x[idx + 1] if
                        (type(x[idx + 1]) in [int, float]) else float('NaN')
                        for x in gtmp
                    ])

            accessEngine.close()
    gname = ':'.join(list(map(str, group)))
    return worker.filterGenotype(genotype, geno_info, var_info, gname)
