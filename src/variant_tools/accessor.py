#!/usr/bin/env python
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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

os.environ['NUMEXPR_MAX_THREADS'] = '8'

import numpy as np
import tables as tb

from .merge_sort_parallel import binarySearch
from .utils import chunks_start_stop, env



class Engine_Storage(object):
    """A factory to make storage engine object
    """

    @staticmethod
    def choose_storage_engine(dbPath):
        """A function to choose which storage engine to start

            Args:

                dbPath: the path to database file
        """
        if dbPath.split(".")[-1] == "h5":
            return HDF5Engine_storage(dbPath)

    # choose_storage_engine=staticmethod(choose_storage_engine)


class Base_Storage(object):
    """An interface for storage APIs
    """

    def __init__(self, dbLocation):
        self.dbPath = dbLocation

    def store(self, data, chr, groupName):
        """The implementation of store API

            Args:

               - data
               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name
        """
        raise NotImplementedError()

    def close(self):
        """This function closes db file
        """

        raise NotImplementedError()


class HDF5Engine_storage(Base_Storage):
    """This is the class for storage of genotype info into HDF5

    """

    def __init__(self, fileName):
        # print("HDF5 engine started")
        Base_Storage.__init__(self, fileName)
        self.file = tb.open_file(self.dbPath, "a")
        self.fileName = fileName

    def store(self, data, chr, groupName=""):

        self.store_genoInfo(data, chr, groupName)

    def store_genoInfo(self, data, chr="", groupName=""):
        filters = tb.Filters(complevel=9, complib='blosc')
        group = self.getGroup(chr)
        if not self.checkGroup(chr, groupName):
            groups = groupName.split("/")
            groupName = groups[-1]
            if len(groups) > 1:
                group = self.getGroup(chr, groups[0])
                # group=self.file.create_group("/chr"+chr+"/"+groups[0],groups[-1])
            if len(data.shape) == 2:
                ds = self.file.create_earray(
                    group,
                    groupName,
                    tb.Atom.from_dtype(data.dtype), (0, data.shape[1]),
                    filters=filters,
                    expectedrows=len(data))
            elif len(data.shape) > 2:
                ds = self.file.create_earray(
                    group,
                    groupName,
                    tb.Atom.from_dtype(data.dtype),
                    (0, data.shape[1], data.shape[2]),
                    filters=filters,
                    expectedrows=len(data))
            else:
                ds = self.file.create_earray(
                    group,
                    groupName,
                    tb.Atom.from_dtype(data.dtype), (0,),
                    filters=filters)
            ds.append(data)

        else:
            if "shape" in groupName:
                group.shape[0] += data[0]
            else:
                self.getGroup(chr, groupName).append(data)

    def geno_fields(self, sampleID):
        fields = []
        chrs = ["X", "Y"]
        chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                for node in self.file.get_node("/chr" + str(chr)):
                    field = node.name
                    if field not in [
                            "Mask", "colnames", "rowmask", "rownames",
                            "samplemask", "shape"
                    ] and field not in fields:
                        if field == "GT":
                            if node[0][0] != -1:
                                fields.append(field.lower())
                        else:
                            fields.append(field.lower())
            except tb.exceptions.NoSuchNodeError:
                pass
        return list(set(fields))

    def remove_variants(self, variantIDs):
        preChr = None
        rownames = None
        group = None
        for res in variantIDs:
            chr = res[1]
            if chr != preChr:
                preChr = chr
                # group=self.file.get_node("/chr"+str(chr)+"/GT/")
                group = self.file.get_node("/chr" + str(chr))
                rownames = group.rownames[:]
            # i=self.rownames.index(variant_id)
            # i=np.where(rownames==res[0])[0][0]
            check = np.where(rownames == res[0])

            if check[0].size != 0:
                group.rowmask[check] = True

    # def recover_variant(self,variant_id,chr,groupName=""):
    #     group=self.file.get_node("/chr"+chr+"/"+groupName)
    #     # i=self.rownames.index(variant_id)
    #     rownames=group.rownames[:]
    #     i=np.where(rownames==variant_id)[0][0]
    #     group.rowmask[i]=False

    def remove_sample(self, sample_id):
        chrs = ["X", "Y"]
        chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                # group=self.file.get_node("/chr"+str(chr)+"/GT/")
                group = self.file.get_node("/chr" + str(chr))
                # i=self.rownames.index(variant_id)
                colnames = group.colnames[:]
                i = np.where(colnames == sample_id)[0][0]

                group.samplemask[i] = True
            except:
                pass

    def remove_genotype(self, cond):
        chrs = []  #["X","Y"]
        chrs.extend(range(1, 20))
        for chr in chrs:
            genoNode = "/chr" + str(chr)
            try:
                node = self.file.get_node(genoNode)
                shape = node.shape[:].tolist()
                chunkPos = chunks_start_stop(shape[0])
                if type(cond) is str:
                    cond = cond.replace("(", "").replace(")", "")
                    # pass
                else:
                    cond = cond[0]
                cond = cond.replace("_geno", "")
                for startPos, endPos in chunkPos:

                    if "GT" in cond and "/chr" + str(chr) + "/GT" in self.file:
                        GT = node.GT[startPos:endPos, :]
                        if "nan" in cond:
                            GT = np.nan_to_num(GT)
                            cond = "GT==0"
                    if "DP" in cond and "/chr" + str(chr) + "/DP" in self.file:
                        DP = node.DP[startPos:endPos, :]
                        DP[DP == -1] = 0
                    if "GQ" in cond and "/chr" + str(chr) + "/GQ" in self.file:
                        GQ = node.GQ[startPos:endPos, :]
                        GQ = np.nan_to_num(GQ)
                    if "GD" in cond and "/chr" + str(chr) + "/GD" in self.file:
                        GD = node.GD[startPos:endPos, :]
                        GD[GD == -1] = 0
                    if "HQ" in cond and "/chr" + str(chr) + "/HQ" in self.file:
                        HQ = node.HQ[startPos:endPos, :]
                        HQ[HQ == -1] = 0
                    if "AD" in cond and "/chr" + str(chr) + "/AD" in self.file:
                        AD = node.AD[startPos:endPos, :]
                        AD[AD == -1] = 0
                    if "PL" in cond and "/chr" + str(chr) + "/PL" in self.file:
                        PL = node.PL[startPos:endPos, :]
                        PL[PL == -1] = 0
                    if "MQ" in cond and "/chr" + str(chr) + "/MQ" in self.file:
                        MQ = node.MQ[startPos:endPos, :]
                        MQ[MQ == -1] = 0
                    if "NS" in cond and "/chr" + str(chr) + "/NS" in self.file:
                        NS = node.NS[startPos:endPos, :]
                        NS[NS == -1] = 0
                    # Mask_geno=group.Mask_geno[:]
                    Mask = np.ones(
                        shape=(endPos - startPos, shape[1]), dtype=np.int8)
                    # node.Mask[startPos:endPos,:]=np.where(eval(cond),np.nan,Mask)
                    node.Mask[startPos:endPos, :] = np.where(
                        eval(cond), -1.0, Mask)
                    startPos = endPos

            except tb.exceptions.NoSuchNodeError as e:
                env.logger.error(f" No such node {e}")

                pass
            except Exception as e:
                # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                env.logger.error(e)
                pass

    def removeNode(self, chrs=[], info=""):
        # for chr in range(1,23):
        if chrs == []:
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                if self.checkGroup(str(chr), info):
                    if len(info) > 0:
                        # print(chr,info,self.fileName)
                        self.file.remove_node("/chr" + str(chr) + "/" + info)
                    else:
                        self.file.remove_node("/chr" + str(chr), recursive=True)
            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                print(e)
                pass

    def remove_genofields(self, items):
        for item in items:
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
            for chr in chrs:
                try:
                    self.file.remove_node("/chr" + str(chr) + "/" + item)

                except tb.exceptions.NoSuchNodeError:
                    pass
                except Exception as e:
                    # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                    print(e)
                    pass

    def checkGroup(self, chr, groupName=""):
        """This function checks the existence of a group

            Args:

               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name

            Returns:

                bool : True if group exists
        """
        if chr.startswith("/chr"):
            chr = chr.replace("/chr", "")
        node = "/chr" + chr
        exist = True if node in self.file else False
        if len(groupName) > 0:
            node = "/chr" + chr + "/" + groupName
            exist = True if node in self.file else False
        return exist

    def getGroup(self, chr, groupName=""):
        """This function gets the node of specified group.

            Args:

               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name

            Returns:

                the node of the group

        """
        if chr.startswith("/chr"):
            chr = chr.replace("/chr", "")
        elif chr.startswith("chr"):
            chr = chr.replace("chr", "")
        node = "/chr" + chr
        if self.checkGroup(chr):
            group = self.file.get_node(node)
        else:
            group = self.file.create_group("/", "chr" + chr, "chromosome")
        if len(groupName) > 0:
            node = "/chr" + chr + "/" + groupName

            if self.checkGroup(chr, groupName):
                group = self.file.get_node(node)
            else:
                group = self.file.create_group("/chr" + chr, groupName)
        return group

    def close(self):
        """This function closes the HDF5 file.

        """
        self.file.close()


class Engine_Access(object):
    """A factory to make data access engine
    """

    @staticmethod
    def choose_access_engine(dbPath, read_only=False):
        """A function to choose which access engine to start

            Args:

                dbPath: the path to database file
        """
        if dbPath.split(".")[-1] == "h5":
            return HDF5Engine_access(dbPath, read_only=read_only)


class Base_Access(object):
    """An interface for data access APIs
    """

    def __init__(self, dbLocation):
        self.dbPath = dbLocation

    def get_geno_info_by_variant_IDs(self, variant_ids, chr, groupName):
        """This function gets a slice of genotype info specified by a list of variant ids.
            **The position of these variants should be consecutive.**

            Args:

                - rowIDs (list): a list of variant ids.
                - chr (string): the chromosome
                - groupName (string): the group name, if not specified, default is genotype

            Returns:

                - A HMatrix object

        """
        raise NotImplementedError()

    def get_geno_info_by_sample_ID(self, sample_id, chr, groupName):
        """This function gets the genotype info of a sample specified by the sample ID in the colnames.

            Args:

                - sampleID : the sample ID stored in colnames
                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict: snpdict[variant_id]=genotype value

        """
        raise NotImplementedError()

    def get_geno_info_by_group(self, groupName, chr):
        """This function gets the genotype info of specified group into a dictionary

            Args:

                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict : snpdict[sample_id][variant_id]=genotype value

        """
        raise NotImplementedError()

    def get_geno_info_by_variant_ID(self, variant_id, chr, groupname):
        """This function gets the genotype info of a variant specified by the variant_id stored in rownames.

            Args:

                - variantID : the variant ID stored in rownames
                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - variant_id (string): the variant_id of the variant
                - indices (list): the position of samples
                - data (list): the genotype info

        """
        raise NotImplementedError()

    def get_geno_info_by_row_pos(self, rowPos, chr, grouopName):
        """This function gets the genotype info of a variant specified by the position of the variant in the matrix.

            Args:
                - rowpos (int): the position of the variant in the matrix
                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - variant_id (string): the variant_id of the variant
                - indices (list): the position of samples
                - data (list): the genotype info

        """

        raise NotImplementedError()

    def close(self):
        """This function closes the HDF5 file.

        """
        raise NotImplementedError()


class HDF5Engine_access(Base_Access):
    """This is the class for access genotype info in HDF5

    """

    def __init__(self, fileName, read_only=False):
        # print("HDF5 engine started")
        Base_Access.__init__(self, fileName)
        self.fileName = fileName
        self.rownames = None
        self.colnames = None
        self.indptr = None
        self.indices = None
        self.data = None
        self.shape = None
        self.chr = None
        self.file = tb.open_file(self.dbPath, "r" if read_only else "a")

    def __load_HDF5_by_chr(self, chr, groupName=""):
        """This function loads matrix, rownames and colnames of specified chromosome into memory

            Args:

               - chr (string): the chromosome


        """
        self.__load_HDF5_by_group(chr, groupName)

    def getGroup(self, chr, groupName=""):
        """This function gets the node of specified group.

            Args:

               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name

            Returns:

                - the node of the group

        """
        if self.chr is None:
            self.chr = chr
        if chr.startswith("/chr"):
            chr = chr.replace("/chr", "")
        # with tb.open_file(self.fileName) as f:
        node = "/chr" + chr
        if node not in self.file:
            self.file.create_group("/", "chr" + chr, "chromosome")
        group = self.file.get_node(node)
        groupName = str(groupName)
        if len(groupName) > 0:
            node = "/chr" + chr + "/" + groupName
            if node not in self.file:
                self.file.create_group("/chr" + chr, groupName)
            group = self.file.get_node(node)
        return group

    def checkGroup(self, chr, groupName=""):
        """This function checks the existence of a group

            Args:

                - chr (string): the chromosome
                - groupName (string): the group name, for example gene name

            Returns:

                - bool : True if group exists

        """
        if chr.startswith("/chr"):
            chr = chr.replace("/chr", "")
        node = "/chr" + chr
        exist = True if node in self.file else False
        if len(groupName) > 0:
            node = "/chr" + chr + "/" + groupName
            exist = True if node in self.file else False
        return exist

    def num_genoinfo(self, sampleID, expr, cond):
        num = 0
        exprCols = expr.split("(")
        method = exprCols[0]
        info = exprCols[1].replace(")", "")
        info = info.replace("_geno", "")
        numVariants = 0
        if method == "min":
            num = 1000000
        for rownames, colnames, sub_all in self.get_all_genotype_genoinfo(
            [sampleID], [], [info]):
            # genotype = sub_all[0]
            genoinfo = sub_all[1]
            if method == "avg":
                localsum = np.nansum(genoinfo).item()
                num += localsum
                numVariants += len(rownames)
            elif method == "min":
                genoinfo[genoinfo == -1] = 0
                localmin = np.nanmin(genoinfo).item()
                if localmin < num:
                    num = localmin
            elif method == "max":
                genoinfo[genoinfo == -1] = 0
                localmax = np.nanmax(genoinfo).item()
                if localmax > num:
                    num = localmax
        if method == "avg":
            num = num / numVariants
        return num

    def num_variants(self, sampleID):
        totalNum = 0
        numCount = {}
        for rownames, colnames, genoinfo in self.get_all_genotype([sampleID]):
            numVariants = len(rownames[:])
            # samplePos=np.where(colnames==sampleID)
            # colPos=np.where(indices==samplePos[0][0])
            colPos = np.where(colnames == sampleID)[0]
            # data=group.data[colPos]
            data = genoinfo[:, colPos]
            numNan = np.where(np.isnan(data))
            # numNone = np.where(data == -1)
            # totalNum+=numVariants-len(numNan[0])-len(numNone[0])
            totalNum += numVariants - len(numNan[0])
            numCount[chr] = numVariants - len(numNan[0])
        return totalNum, numCount

    def num_genotypes(self, sampleID, cond, genotypes):
        totalNum = 0

        for rownames, colnames, genoinfo in self.get_all_genotype_filter(
            [sampleID], genotypes):
            if cond is None:
                numNan = np.where(np.isnan(genoinfo))
                # numNone = np.where(genoinfo == -1)
                totalNum += len(rownames) - len(numNan[0])
            else:
                numNan = np.where(np.isnan(genoinfo))
                if cond == "GT!=0":
                    numCond = np.where(genoinfo != 0)
                    totalNum += len(numCond[0]) - len(numNan[0])
                else:
                    GTcond = cond.split("=")[1]
                    numCond = np.where(genoinfo == int(GTcond))
                    totalNum += len(numCond[0])
        return totalNum

    def sum_genotypes(self, sampleID, cond, genotypes):
        totalGeno = 0
        for rownames, colnames, genoinfo in self.get_all_genotype_filter(
            [sampleID], genotypes):
            numGeno = np.nansum(genoinfo)
            totalGeno += numGeno
        return int(totalGeno)

    def get_geno_by_group(self, chr, groupName):
        group = self.getGroup(chr)
        self.colnames = group.colnames[:].tolist()
        group = self.getGroup(chr, groupName)
        self.rownames = group.rownames[:].tolist()

        if "GT" in group:
            self.GT = group.GT[:]

            snpdict = dict.fromkeys(self.colnames, {})
            for key, value in snpdict.items():
                snpdict[key] = dict.fromkeys(self.rownames, (0,))

            for rowidx, rowID in enumerate(self.rownames):
                for colidx, colID in enumerate(self.colnames):
                    snpdict[colID][rowID] = (self.GT[rowidx, colidx],)

        else:
            snpdict = dict.fromkeys(self.colnames, {})
            for key, value in snpdict.items():
                snpdict[key] = dict.fromkeys(self.rownames, (0,))

            for rowidx, rowID in enumerate(self.rownames):
                for colidx, colID in enumerate(self.colnames):
                    snpdict[colID][rowID] = (None,)
        return snpdict

    def get_geno_by_sep_variant_ids(self, rowIDs, chr, groupName=""):

        group = self.getGroup(chr, groupName)
        # rownames=group.rownames[:].tolist()
        rownames = group.rownames[:]
        colnames = group.colnames[:]
        rowMask = group.rowmask[:]
        sampleMask = group.samplemask[:]

        # rowPos=[rownames.index(id) for id in rowIDs]

        try:
            rowPos = [np.where(rownames == id)[0][0] for id in rowIDs]
            update_rownames = rownames[rowPos]
            sub_geno = group.GT[rowPos, :]
            sub_Mask = group.Mask[rowPos, :]
            sub_Mask = sub_Mask.astype(float)
            sub_Mask[sub_Mask == -1.0] = np.nan
            sub_geno = np.multiply(sub_geno, sub_Mask)
            update_rowMask = rowMask[rowPos]
            rowMasked = np.where(update_rowMask == True)[0]
            sampleMasked = np.where(sampleMask == True)[0]
            if len(rowMasked) > 0:
                update_rownames = np.array(update_rownames)[np.where(
                    update_rowMask == False)[0]]
                sub_geno = np.delete(sub_geno, rowMasked, 0)

            if len(sampleMasked) > 0:
                colnames = colnames[np.where(sampleMask == False)[0]]
                sub_geno = np.delete(sub_geno, sampleMasked, 1)
            return np.array(update_rownames), colnames, np.array(sub_geno)
        except IndexError:
            update_rownames = []
            sub_geno = []
            for id in rowIDs:
                try:
                    rowPos = np.where(rownames == id)[0][0]
                    update_rownames.append(rownames[rowPos])
                    sub_geno.append(group.GT[rowPos, :])
                except IndexError:
                    update_rownames.append(id)
                    sub_geno.append(np.full(len(colnames), np.nan))
            return np.array(update_rownames), colnames, np.array(sub_geno)
            # return np.full(len(rowIDs),np.nan),colnames,None
        except NameError:
            env.logger.error(
                "varaintIDs of this gene are not found on this chromosome {}"
                .format(chr))
        except Exception as e:
            print(e)

    def get_genotype(self, variantIDs, sampleNames, chrs=""):
        updated_rownames = []
        updated_geno = []
        updated_colnames = []
        if chrs == "":
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
        for chr in chrs:
            try:

                node = self.file.get_node("/chr" + str(chr))
                rownames = node.rownames[:].tolist()
                # print(self.fileName,variantIDs,chr,rownames)
                colnames = node.colnames[:].tolist()
                colpos = []
                if (len(sampleNames) > 0):
                    colpos = list(map(lambda x: colnames.index(x), sampleNames))
                else:
                    colpos = list(map(lambda x: colnames.index(x), colnames))

                if len(variantIDs) > 0:
                    for id in variantIDs:
                        try:
                            minPos = rownames.index(id)
                            break
                        except ValueError:
                            continue
                    for id in reversed(variantIDs):
                        try:
                            maxPos = rownames.index(id)
                            maxPos = maxPos + 1
                            break
                        except ValueError:
                            continue
                    try:
                        minPos
                        maxPos
                        # genoinfo=node.GT[minPos:maxPos,colpos]
                        genoinfo = node.GT[minPos:maxPos, :]

                        updated_rownames, updated_colnames, updated_geno = self.filter_removed_genotypes(
                            minPos, maxPos, genoinfo, node, colpos, [], "GT")

                    except NameError:
                        env.logger.error(
                            "varaintIDs of this gene are not found on this chromosome {}"
                            .format(chr))
                else:
                    shape = node.shape[:].tolist()
                    chunkPos = chunks_start_stop(shape[0])
                    for minPos, maxPos in chunkPos:
                        if "/chr" + str(chr) + "/GT" in self.file:
                            genoinfo = node.GT[minPos:maxPos, colpos]
                            sub_rownames, updated_colnames, sub_geno = self.filter_removed_genotypes(
                                minPos, maxPos, genoinfo, node, colpos, [],
                                "GT")
                            updated_rownames.extend(sub_rownames)
                            updated_geno.extend(sub_geno)
            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                print(e)
                pass
        return np.array(updated_rownames), np.array(updated_colnames), np.array(
            updated_geno)

    def get_genotype_by_chunk(self, sampleNames, chrs=""):

        if chrs == "":
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                node = self.file.get_node("/chr" + str(chr))
                colnames = node.colnames[:].tolist()
                colpos = []
                if (len(sampleNames) > 0):
                    colpos = list(map(lambda x: colnames.index(x), sampleNames))
                else:
                    colpos = list(map(lambda x: colnames.index(x), colnames))
                shape = node.shape[:].tolist()
                chunkPos = chunks_start_stop(shape[0])
                for minPos, maxPos in chunkPos:
                    if "/chr" + str(
                            chr) + "/GT" in self.file and minPos != maxPos:
                        genoinfo = node.GT[minPos:maxPos, colpos]
                        sub_rownames, updated_colnames, sub_geno = self.filter_removed_genotypes(
                            minPos, maxPos, genoinfo, node, colpos, [])
                        yield np.array(sub_rownames), np.array(
                            updated_colnames), np.array(sub_geno)

            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                print(e)
                pass

    def get_all_genotype(self, sampleNames, chrs=""):
        return (list(self.get_genotype_by_chunk(sampleNames, chrs)))

    def get_genotype_by_chunk_filter(self, sampleNames, geno_cond, chrs=""):

        if chrs == "":
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                node = self.file.get_node("/chr" + str(chr))
                colnames = node.colnames[:].tolist()
                colpos = []
                if (len(sampleNames) > 0):
                    colpos = list(map(lambda x: colnames.index(x), sampleNames))
                else:
                    colpos = list(map(lambda x: colnames.index(x), colnames))
                shape = node.shape[:].tolist()
                chunkPos = chunks_start_stop(shape[0])
                for minPos, maxPos in chunkPos:
                    if "/chr" + str(
                            chr) + "/GT" in self.file and minPos != maxPos:
                        # genoinfo = node.GT[minPos:maxPos, colpos]
                        sub_rownames, updated_colnames, sub_geno = self.filter_on_genotypes(
                            geno_cond, chr, node, "GT", minPos, maxPos, colpos,
                            [])
                        # sub_rownames,updated_colnames,sub_geno=self.filter_removed_genotypes(minPos,maxPos,genoinfo,node,colpos,[])
                        yield np.array(sub_rownames), np.array(
                            updated_colnames), np.array(sub_geno)

            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                print(e)
                pass

    def get_all_genotype_filter(self, sampleNames, geno_cond, chrs=""):
        return (list(
            self.get_genotype_by_chunk_filter(sampleNames, geno_cond, chrs)))

    def get_genotype_genoinfo_by_chunk(self,
                                       sampleNames,
                                       varIDs,
                                       validGenotypeFields,
                                       cond,
                                       chrs=""):

        if chrs == "":
            chrs = ["X", "Y"]
            chrs.extend(range(1, 23))
        for chr in chrs:
            try:
                node = self.file.get_node("/chr" + str(chr))
                colnames = node.colnames[:].tolist()
                colpos = []
                sampleNames.sort()

                if (len(sampleNames) > 0):
                    colpos = list(map(lambda x: colnames.index(x), sampleNames))
                # else:
                #     colpos=list(map(lambda x:colnames.index(x),colnames))
                shape = node.shape[:].tolist()
                chunkPos = chunks_start_stop(shape[0])
                rowpos = []

                if len(varIDs) > 0:
                    rownames = node.rownames[:].tolist()
                    rowpos = np.array([0] * len(rownames))
                    selectPos = [
                        self.find_element_in_list(id, rownames) for id in varIDs
                    ]
                    selected = [x for x in selectPos if x is not None]
                    rowpos[selected] = 1
                    if np.sum(rowpos) == 0:
                        continue

                for minPos, maxPos in chunkPos:
                    sub_all = []
                    if "/chr" + str(
                            chr) + "/GT" in self.file and minPos != maxPos:

                        sub_rownames, updated_colnames, sub_geno = self.filter_on_genotypes(
                            cond, chr, node, "GT", minPos, maxPos, colpos,
                            rowpos)

                        sub_all.append(np.array(sub_geno))
                        if len(validGenotypeFields) > 0:
                            for pos, field in enumerate(validGenotypeFields):
                                _, _, sub_info = self.filter_on_genotypes(
                                    cond, chr, node, field, minPos, maxPos,
                                    colpos, rowpos)
                                sub_all.append(np.array(sub_info))

                        yield np.array(sub_rownames), np.array(
                            updated_colnames), sub_all

            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                print(e)
                pass

    def get_all_genotype_genoinfo(self,
                                  sampleNames,
                                  varIDs,
                                  validGenotypeFields,
                                  genotypes="",
                                  chrs=""):
        return (list(
            self.get_genotype_genoinfo_by_chunk(sampleNames, varIDs,
                                                validGenotypeFields, genotypes,
                                                chrs)))

    def filter_removed_genotypes(self,
                                 minPos,
                                 maxPos,
                                 genoinfo,
                                 node,
                                 colpos,
                                 rowpos,
                                 field=""):
        #assume rowIDs are sorted by genome position
        rownames = node.rownames[minPos:maxPos]
        colnames = node.colnames[:]
        rowMask = node.rowmask[minPos:maxPos]

        # sampleMask=node.samplemask[:]
        # if field=="GT":
        #     genoinfo[genoinfo==-1]=0
        #     genoinfo=np.nan_to_num(genoinfo)

        try:
            sub_Mask = node.Mask[minPos:maxPos, :]
            sub_Mask = sub_Mask.astype(float)
            sub_Mask[sub_Mask == -1.0] = np.nan
            sub_geno = np.multiply(genoinfo, sub_Mask)

            rowMasked = np.where(rowMask == True)[0]

            # sampleMasked=np.where(sampleMask==True)[0]

            if len(rowpos) > 0:
                selectRows = rowpos[minPos:maxPos]
                rownames = rownames[selectRows == 1]
                sub_geno = sub_geno[selectRows == 1, :]

            colnames = colnames[colpos]
            sub_geno = sub_geno[:, colpos]

            if len(rowpos) == 0 and len(rowMasked) > 0:

                rownames = rownames[np.where(rowMask == False)]
                sub_geno = np.delete(sub_geno, rowMasked, 0)

            # if len(sampleMasked)>0:
            #     # colnames=colnames[np.where(sampleMask==False)]
            #     # sub_geno=np.delete(sub_geno,sampleMasked,1)
            #     colnames=colnames[colpos]
            #     sub_geno=sub_geno[:,colpos]
            return np.array(rownames), colnames, np.array(sub_geno)
        except Exception as e:
            print(e)

    def filter_on_genotypes(self, cond, chr, node, field, startPos, endPos,
                            colpos, rowpos):
        genoinfo = None

        if field in self.file.get_node("/chr" + str(chr)):
            genoinfo = self.file.get_node("/chr" + str(chr) + "/" +
                                          field)[startPos:endPos, :]
            if field != "GT":
                genoinfo[genoinfo == -1] = 0
                genoinfo = np.nan_to_num(genoinfo)
            else:
                if env.treat_missing_as_wildtype:
                    genoinfo = np.nan_to_num(genoinfo)

        # if field=="GT":
        #     genoinfo=node.GT[startPos:endPos,:]
        # if field=="DP" and "/chr"+str(chr)+"/DP" in self.file:
        #     genoinfo=node.DP[startPos:endPos,:]
        #     genoinfo[genoinfo==-1]=0
        # if field=="GQ" and "/chr"+str(chr)+"/GQ" in self.file:
        #     genoinfo=node.GQ[startPos:endPos,:]
        #     genoinfo=np.nan_to_num(genoinfo)

        if len(cond) > 0:
            if type(cond) is str:
                cond = cond.replace("(", "").replace(")", "")
                cond = cond.replace("_geno", "")
            else:
                cond = cond[0]
                cond = cond.replace("_geno", "")
            if "DP" in cond and "/chr" + str(chr) + "/DP" in self.file:
                DP = node.DP[startPos:endPos, :]
                DP[DP == -1] = 0
            if "GQ" in cond and "/chr" + str(chr) + "/GQ" in self.file:
                GQ = node.GQ[startPos:endPos, :]
                GQ = np.nan_to_num(GQ)
            if "GD" in cond and "/chr" + str(chr) + "/GD" in self.file:
                GD = node.GD[startPos:endPos, :]
                GD[GD == -1] = 0
            if "HQ" in cond and "/chr" + str(chr) + "/HQ" in self.file:
                HQ = node.HQ[startPos:endPos, :]
                HQ[HQ == -1] = 0
            if "AD" in cond and "/chr" + str(chr) + "/AD" in self.file:
                AD = node.AD[startPos:endPos, :]
                AD[AD == -1] = 0
            if "PL" in cond and "/chr" + str(chr) + "/PL" in self.file:
                PL = node.PL[startPos:endPos, :]
                PL[PL == -1] = 0
            if "MQ" in cond and "/chr" + str(chr) + "/MQ" in self.file:
                MQ = node.MQ[startPos:endPos, :]
                MQ[MQ == -1] = 0
            if "NS" in cond and "/chr" + str(chr) + "/NS" in self.file:
                NS = node.NS[startPos:endPos, :]
                NS[NS == -1] = 0
            genoinfo = np.where(eval("~(" + cond + ")"), np.nan, genoinfo)
        rownames, colnames, genoinfo = self.filter_removed_genotypes(
            startPos, endPos, genoinfo, node, colpos, rowpos)

        return rownames, colnames, genoinfo

    def find_element_in_list(self, element, list_element):
        try:
            index_element = list_element.index(element)
            return index_element
        except ValueError:
            pass

    def get_geno_by_row_pos(self,
                            rowpos,
                            chr,
                            sortedID,
                            sampleNames,
                            validGenotypeFields,
                            groupName=""):
        try:
            sub_all = []
            node = self.file.get_node("/chr" + str(chr))
            colpos = []
            colnames = node.colnames[:].tolist()
            sampleNames.sort()

            if (len(sampleNames) > 0):
                colpos = list(map(lambda x: colnames.index(x), sampleNames))
            else:
                colpos = list(map(lambda x: colnames.index(x), colnames))
            pos = binarySearch(sortedID, 0, len(sortedID) - 1, rowpos)

            if pos != -1:
                posInNode = sortedID[pos][1]
                _, _, row_geno = self.filter_on_genotypes([], chr, node, "GT",
                                                          posInNode,
                                                          posInNode + 1, colpos,
                                                          [])
                sub_all.append(np.array(row_geno))
                if len(validGenotypeFields) > 0:
                    for pos, field in enumerate(validGenotypeFields):
                        _, _, row_info = self.filter_on_genotypes([], chr, node,
                                                                  field,
                                                                  posInNode,
                                                                  posInNode + 1,
                                                                  colpos, [])
                        # row_geno[0].extend(np.array(row_info))
                        sub_all.append(np.array(row_info))
                return sub_all
            else:
                #return np.zeros(shape=(len(node.GT[1,colpos])),dtype=int)
                # return [None]*len(node.GT[1,colpos])

                sub_all.append([np.full(len(node.GT[1, colpos]), np.nan)])
                if len(validGenotypeFields) > 0:
                    for pos, field in enumerate(validGenotypeFields):
                        sub_all.append(
                            [np.full(len(node.GT[1, colpos]), np.nan)])
                return sub_all
        except tb.exceptions.NoSuchNodeError:
            # return np.zeros(shape=(len(node.GT[1,colpos])),dtype=int)
            # return [None]*len(node.GT[1,colpos])

            sub_all.append(np.full(len(node.GT[1, colpos]), np.nan))
            if len(validGenotypeFields) > 0:
                for pos, field in enumerate(validGenotypeFields):
                    sub_all.append([np.full(len(node.GT[1, colpos]), np.nan)])
            return sub_all
        except Exception as e:
            print("exception", rowpos, chr, pos, sortedID[pos])
            print(e)
            pass

    def get_geno_by_variant_ID(self, variantID, chr, groupName=""):
        pass

    def show_file_node(self):
        """This function prints the nodes in the file.

        """
        for node in self.file:
            print(node)

    def get_rownames(self, chr, groupName=""):
        """This function gets variant IDs of specified group.

            Args:

                - chr (string): the chromosome
                - groupName (string): the group name, for example gene name

            Return:
                - a ndarray of variant IDs

        """
        group = self.getGroup(chr, groupName)
        return group.rownames[:]

    def get_colnames(self, chr, groupName=""):
        """This function gets sample IDs of specified group.

            Args:

                - chr (string): the chromosome
                - groupName (string): the group name, for example gene name

            Return:
                a ndarray of sample IDs

        """
        group = self.getGroup(chr, groupName)
        return group.colnames[:]

    def get_shape(self, chr, groupName=""):
        """This function gets shape of specified group.

            Args:

                - chr (string): the chromosome
                - groupName (string): the group name, for example gene name

            Return:
                shape of matrix

        """
        group = self.getGroup(chr, groupName)
        return group.shape[:]

    def close(self):
        """This function closes db file
        """

        self.file.close()


# class AccessEachHDF5(Process):

#     def __init__(self,fileName,variantID,result,chr,groupName=""):
#         Process.__init__(self)
#         self.hdf5=HDF5Engine_access(fileName)
#         self.variantID=variantID
#         self.chr=chr
#         self.groupName=groupName
#         self.result=result

#     def run(self):
#         variant_ID,indices,data=self.hdf5.get_geno_info_by_variant_ID(self.variantID,self.chr,self.groupName)
#         print((variant_ID,len(indices)))
#         colnames=self.hdf5.get_colnames(self.chr,self.groupName)
#         for idx,col in enumerate(indices):
#             if not math.isnan(data[idx]):
#                 self.result[colnames[col]]=int(data[idx])
#             else:
#                 self.result[colnames[col]]=data[idx]
#         self.hdf5.close()

# class HDF5Engine_access_multi:
#     """This is the class for access genotype info stored in multiple HDF5 (under development).

#     """
#     def __init__(self,files,jobs):
#         self.files=files
#         self.jobs=jobs

#     def get_geno_info_by_variant_ID(self,variantID,chr,groupName=""):
#         """This function gets the genotype info of a variant specified by the variant_id stored in rownames.

#             Args:

#                 - variantID : the variant ID stored in rownames
#                 - chr (string): the chromosome
#                 - groupName (string): the group which variants are in, for example gene name

#             Returns:

#                 right now, just print a dict of dict[sampleID]=genotype info

#         """
#         taskQueue=queue.Queue()
#         importers=[None]*self.jobs
#         result=Manager().dict()
#         for file in self.files:
#             taskQueue.put(AccessEachHDF5(file,variantID,result,chr,groupName))
#         while taskQueue.qsize()>0:
#             for i in range(self.jobs):
#                 if importers[i] is None or not importers[i].is_alive():
#                     task=taskQueue.get()
#                     importers[i]=task
#                     importers[i].start()
#                     break
#         for worker in importers:
#             worker.join()
#         print(result)
