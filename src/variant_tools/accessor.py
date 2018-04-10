#!/usr/bin/env python
from scipy.sparse import csr_matrix
import tables as tb
import numpy as np
import math

import os
import sys
from .utils import env,chunks_start_stop

from multiprocessing import Process,Manager
import queue




class Engine_Storage(object):
    """A factory to make storage engine object
    """
    @staticmethod
    def choose_storage_engine(dbPath):
        """A function to choose which storage engine to start

            Args:

                dbPath: the path to database file
        """
        if dbPath.split(".")[-1]=="h5":
            return HDF5Engine_storage(dbPath)
      

    # choose_storage_engine=staticmethod(choose_storage_engine)




class Base_Storage(object):
    """An interface for storage APIs
    """

    def __init__(self,dbLocation):
        self.dbPath=dbLocation

    def store(self,data,chr,groupName):
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

    def __init__(self,fileName):
        # print("HDF5 engine started")
        Base_Storage.__init__(self,fileName)
        self.file=tb.open_file(self.dbPath,"a")
        self.fileName=fileName


    def store(self,data,chr,groupName=""):

        self.store_genoInfo(data,chr,groupName)
        
       


    def store_genoInfo(self,data,chr="",groupName=""):
        filters = tb.Filters(complevel=9, complib='blosc')       
        group=self.getGroup(chr)
        if not self.checkGroup(chr,groupName):
            groups=groupName.split("/")
            groupName=groups[-1]
            if len(groups)>1:
                group=self.getGroup(chr,groups[0])
                # group=self.file.create_group("/chr"+chr+"/"+groups[0],groups[-1])
            
            if len(data.shape)==2:
                ds = self.file.create_earray(group, groupName, tb.Atom.from_dtype(data.dtype), (0,data.shape[1]),filters=filters,expectedrows=len(data))
            elif len(data.shape)>2:
                ds = self.file.create_earray(group, groupName, tb.Atom.from_dtype(data.dtype), (0,data.shape[1],data.shape[2]),filters=filters,expectedrows=len(data))
            else:
                ds = self.file.create_earray(group, groupName, tb.Atom.from_dtype(data.dtype), (0,),filters=filters)
            ds.append(data)

        else:
            if "shape" in groupName:
                group.shape[0]+=data[0]
            else:
                self.getGroup(chr,groupName).append(data)


               


    def num_variants(self,sampleID):
        totalNum=0
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        numCount={}
        for chr in chrs:
            try:
                group=self.file.get_node("/chr"+str(chr))
                # indices=group.indices[:]
                colnames=group.colnames[:]
                numVariants=len(group.rownames[:])
                # samplePos=np.where(colnames==sampleID)
                # colPos=np.where(indices==samplePos[0][0])
                colPos=np.where(colnames==sampleID)[0]
                # data=group.data[colPos]
                data=group.GT_geno[:,colPos]
                numNan=np.where(np.isnan(data))
                numNone=np.where(data==-1) 
                # totalNum+=numVariants-len(numNan[0])-len(numNone[0])
                totalNum+=numVariants-len(numNan[0])
                numCount[chr]=numVariants-len(numNan[0])
                # totalNum+=numVariants
            except tb.exceptions.NoSuchNodeError:
                pass
        return totalNum,numCount

    def to_csr_matrix(self,group):
        return csr_matrix((group.data[:],group.indices[:],group.indptr[:]),shape=group.shape[:])



    def geno_fields(self,sampleID):
        fields=[]
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
            try:
                for field in ["GT_geno","DP_geno","GQ_geno"]:
                    group=self.file.get_node("/chr"+str(chr)+"/"+field)  
                    if field=="GT_geno":
                        field="GT"     
                    fields.append(field)
            except tb.exceptions.NoSuchNodeError:
                pass
        return list(set(fields))



    def remove_variants(self,variantIDs):
        preChr=None
        rownames=None
        group=None
        for res in variantIDs:
            chr=res[1]
            if chr!=preChr:
                preChr=chr
                # group=self.file.get_node("/chr"+str(chr)+"/GT/")
                group=self.file.get_node("/chr"+str(chr))
                rownames=group.rownames[:]
            # i=self.rownames.index(variant_id)
            i=np.where(rownames==res[0])[0][0]
            group.rowmask[i]=True

    def recover_variant(self,variant_id,chr,groupName=""):
        group=self.file.get_node("/chr"+chr+"/"+groupName)
        # i=self.rownames.index(variant_id)
        rownames=group.rownames[:]
        i=np.where(rownames==variant_id)[0][0]
        group.rowmask[i]=False

    
    def remove_sample(self,sample_id):
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
            try:
                # group=self.file.get_node("/chr"+str(chr)+"/GT/")
                group=self.file.get_node("/chr"+str(chr))
                # i=self.rownames.index(variant_id)
                colnames=group.colnames[:]
                i=np.where(colnames==sample_id)[0][0]

                group.samplemask[i]=True
            except:
                pass
        

    def remove_genotype(self,cond):
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
            genoNode="/chr"+str(chr)
            try:
                node=self.file.get_node(genoNode)
                shape=node.shape[:].tolist()
                chunkPos=chunks_start_stop(shape[0])
                for startPos,endPos in chunkPos:

                    if "GT_geno" in cond and "/chr"+str(chr)+"/GT_geno" in self.file: 
                        GT_geno=node.GT_geno[startPos:endPos,:]
                        if "nan" in cond:
                            GT_geno=np.nan_to_num(GT_geno)
                            cond="GT_geno==0"
                    if "DP_geno" in cond and "/chr"+str(chr)+"/DP_geno" in self.file:    
                        DP_geno=node.DP_geno[startPos:endPos,:]
                        DP_geno[DP_geno==-1]=0
                    if "GQ_geno" in cond and "/chr"+str(chr)+"/GQ_geno" in self.file:    
                        GQ_geno=node.GQ_geno[startPos:endPos,:]
                        GQ_geno=np.nan_to_num(GQ_geno)
                    # Mask_geno=group.Mask_geno[:]
                    Mask_geno=np.ones(shape=(endPos-startPos,shape[1]),dtype=np.int8)
                    node.Mask_geno[startPos:endPos,:]=np.where(eval(cond),np.nan,Mask_geno)
                    startPos=endPos
                # print("mask",group.Mask_geno[:])
                # print(group.Mask_geno[:])
                # print(group.DP_geno[:])
                # print(group.GQ_geno[:])
            except tb.exceptions.NoSuchNodeError:
                pass 
            except Exception as e:
                # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                print(e)
                pass        
  

    def removeNode(self,chr,info=""):
        # for chr in range(1,23):
        try:
            if self.checkGroup(str(chr),info):
                if len(info)>0:
                    self.file.remove_node("/chr"+str(chr)+"/"+info)
                else:
                    self.file.remove_node("/chr"+str(chr),recursive=True)
        except tb.exceptions.NoSuchNodeError:
            pass                              
        except Exception as e:
            print(e)
            pass


    def remove_genofields(self,items):
        for item in items:
            chrs=["X","Y"]
            chrs.extend(range(1,23))
            for chr in chrs:
                genoNode="/chr"+str(chr)
                try:
                    self.file.remove_node("/chr"+str(chr)+"/"+item)
                   
                except tb.exceptions.NoSuchNodeError:
                    pass 
                except Exception as e:
                    # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                    print(e)
                    pass        



    def checkGroup(self,chr,groupName=""):
        """This function checks the existence of a group

            Args:

               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name

            Returns:

                bool : True if group exists
        """
        if chr.startswith("/chr"):
            chr=chr.replace("/chr","")
        node="/chr"+chr
        exist=True if node in self.file else False
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
            exist=True if node in self.file else False
        return exist



    def getGroup(self,chr,groupName=""):
        """This function gets the node of specified group.

            Args:

               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name

            Returns:

                the node of the group

        """
        if chr.startswith("/chr"):
            chr=chr.replace("/chr","")
        elif chr.startswith("chr"):
            chr=chr.replace("chr","")
        node="/chr"+chr
        if self.checkGroup(chr):
            group=self.file.get_node(node)
        else:
            group=self.file.create_group("/","chr"+chr,"chromosome")
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName

            if self.checkGroup(chr,groupName):
                group=self.file.get_node(node)
            else:
                group=self.file.create_group("/chr"+chr,groupName)
        return group

    def close(self):
        """This function closes the HDF5 file.

        """
        self.file.close()



class Engine_Access(object):
    """A factory to make data access engine
    """
    @staticmethod
    def choose_access_engine(dbPath):
        """A function to choose which access engine to start

            Args:

                dbPath: the path to database file
        """
        if dbPath.split(".")[-1]=="h5":
            return HDF5Engine_access(dbPath)





class Base_Access(object):
    """An interface for data access APIs
    """

    def __init__(self,dbLocation):
        self.dbPath=dbLocation

    def get_geno_info_by_variant_IDs(self,variant_ids,chr,groupName):
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

    def get_geno_info_by_sample_ID(self,sample_id,chr,groupName):
        """This function gets the genotype info of a sample specified by the sample ID in the colnames. 

            Args:

                - sampleID : the sample ID stored in colnames
                - chr (string): the chromosome 
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict: snpdict[variant_id]=genotype value

        """
        raise NotImplementedError()

    def get_geno_info_by_group(self,groupName,chr):
        """This function gets the genotype info of specified group into a dictionary

            Args:

                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict : snpdict[sample_id][variant_id]=genotype value

        """
        raise NotImplementedError();

    def get_geno_info_by_variant_ID(self,variant_id,chr,groupname):
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

    def get_geno_info_by_row_pos(self,rowPos,chr,grouopName):
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
    def __init__(self,fileName):
        # print("HDF5 engine started")
        Base_Access.__init__(self,fileName)
        self.fileName=fileName
        self.rownames=None
        self.colnames=None
        self.indptr=None
        self.indices=None
        self.data=None
        self.shape=None
        self.chr=None
        self.file=tb.open_file(self.dbPath,"a")
 

    
    def __load_HDF5_by_chr(self,chr,groupName=""):
        """This function loads matrix, rownames and colnames of specified chromosome into memory

            Args:

               - chr (string): the chromosome 


        """
        self.__load_HDF5_by_group(chr,groupName)


    def getGroup(self,chr,groupName=""):
        """This function gets the node of specified group.

            Args:

               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name

            Returns:

                - the node of the group

        """
        if self.chr is None:
            self.chr=chr
        if chr.startswith("/chr"):
            chr=chr.replace("/chr","")
        # with tb.open_file(self.fileName) as f:
        node="/chr"+chr
        if node not in self.file:
            self.file.create_group("/","chr"+chr,"chromosome")
        group=self.file.get_node(node)

        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName

            if node not in self.file:
                self.file.create_group("/chr"+chr,groupName)
            group=self.file.get_node(node)
        return group

    def checkGroup(self,chr,groupName=""):
        """This function checks the existence of a group

            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Returns:

                - bool : True if group exists

        """
        if chr.startswith("/chr"):
            chr=chr.replace("/chr","")
        node="/chr"+chr
        exist=True if node in self.file else False
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
            exist=True if node in self.file else False
        return exist


    def num_genoinfo(self,sampleID,expr,cond):
        num=0
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
            try:
                group=self.file.get_node("/chr"+str(chr))
                colnames=group.colnames[:]
                colPos=np.where(colnames==sampleID)[0]
                exprCols=expr.split("(")
                method=exprCols[0]
                info=exprCols[1].replace(")","")
                # data=group.get_node(info)[:,colPos]
                
                data=self.file.get_node("/chr"+str(chr)+"/"+info)[:,colPos]
                data[data==-1]=0
                if method=="avg":
                    num=np.average(data).item()
                elif method=="min":
                    num=np.nanmin(data).item()
                elif method=="max":
                    num=np.nanmax(data).item()
            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                print(e)
                pass     
        return num



    def num_genotypes(self,sampleID,cond,genotypes):
        totalNum=0
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
            try:
                group=self.file.get_node("/chr"+str(chr))
                colnames=group.colnames[:]
                numVariants=len(group.rownames[:])
                colPos=np.where(colnames==sampleID)[0]
                
                shape=group.shape[:].tolist()   
                
                chunkPos=chunks_start_stop(shape[0])
    
                for startPos,endPos in chunkPos:
                    # rownames=node.rownames[startPos:endPos].tolist() 
                               
                    if "/chr"+str(chr)+"/GT_geno" in self.file:    
                        rownames,colnames,data=self.filter_on_genotypes(genotypes,chr,group,"GT_geno",startPos,endPos,colPos,"")
                 
                        # data=group.GT_geno[:,colPos]
                        if cond is None:
                            numNan=np.where(np.isnan(data))
                            numNone=np.where(data==-1) 
                            totalNum+=numVariants-len(numNan[0])
                        #     totalNum+=numVariants-len(numNan[0])
                        else:
                            numNan=np.where(np.isnan(data))
                            if cond=="GT!=0":
                                numCond=np.where(data!=0)
                                totalNum+=len(numCond[0])-len(numNan[0])
                            else:
                                GTcond=cond.split("=")[1]
                                numCond=np.where(data==int(GTcond))
                                totalNum+=len(numCond[0])
   
            except tb.exceptions.NoSuchNodeError:
                pass
            except Exception as e:
                # env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))
                print(e)
                pass     
        return totalNum




    def get_geno_by_group(self,chr,groupName):
        group=self.getGroup(chr)
        self.colnames=group.colnames[:].tolist()
        group=self.getGroup(chr,groupName)
        self.rownames=group.rownames[:].tolist()
        self.GT=group.GT_geno[:]
        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.items():
            snpdict[key]=dict.fromkeys(self.rownames,(0,))

        for rowidx,rowID in enumerate(self.rownames):
            for colidx,colID in enumerate(self.colnames):
                snpdict[colID][rowID]=(self.GT[rowidx,colidx],)
        return snpdict

     
    def get_geno_by_sep_variant_ids(self,rowIDs,chr,groupName=""):

     
        group=self.getGroup(chr,groupName)
        # rownames=group.rownames[:].tolist()
        rownames=group.rownames[:]
        colnames=group.colnames[:]
        rowMask=group.rowmask[:]
        sampleMask=group.samplemask[:]

        # rowPos=[rownames.index(id) for id in rowIDs]
        rowPos=[np.where(rownames==id)[0][0] for id in rowIDs]
       
        try:
            update_rownames=rownames[rowPos]
            sub_geno=group.GT_geno[rowPos,:]    
            sub_Mask=group.Mask_geno[rowPos,:]
            sub_geno=np.multiply(sub_geno,sub_Mask)

            update_rowMask=rowMask[rowPos]
            rowMasked=np.where(update_rowMask==True)[0]
            sampleMasked=np.where(sampleMask==True)[0]
            if len(rowMasked)>0:
                update_rownames=np.array(update_rownames)[np.where(update_rowMask==False)[0]]
                sub_geno=np.delete(sub_geno,rowMasked,0)

            if len(sampleMasked)>0:
                colnames=colnames[np.where(sampleMask==False)[0]]
                sub_geno=np.delete(sub_geno,sampleMasked,1)
        
            return np.array(update_rownames),colnames,np.array(sub_geno)
        except NameError:
            env.logger.error("varaintIDs of this gene are not found on this chromosome {}".format(chr))



    def get_genotype(self,variantIDs,sampleName,chrs=""):
        updated_rownames=[]
        updated_geno=[]
        updated_colnames=[]
        if chrs=="":
            chrs=["X","Y"]
            chrs.extend(range(1,23))
        for chr in chrs:
            try:
                node=self.file.get_node("/chr"+str(chr))
                rownames=node.rownames[:].tolist()
                colnames=node.colnames[:].tolist()
                colpos=[]
                if (len(sampleName)>0):
                    colpos=colnames.index(sampleName)
                else:
                    colpos=list(map(lambda x:colnames.index(x),colnames))

                if len(variantIDs)>0:
                    for id in variantIDs:
                        try:
                            minPos=rownames.index(id)
                            break
                        except ValueError:
                            continue
                    for id in reversed(variantIDs):
                        try:
                            maxPos=rownames.index(id)
                            maxPos=maxPos+1
                            break
                        except ValueError:
                            continue
                    try:
                        minPos
                        maxPos
                        genoinfo=node.GT_geno[minPos:maxPos,colpos]           
                        updated_rownames,updated_colnames,updated_geno=self.filter_removed_genotypes(minPos,maxPos,genoinfo,node,colpos,[])
                        
                    except NameError:
                        env.logger.error("varaintIDs of this gene are not found on this chromosome {}".format(chr))
                else:
                    shape=node.shape[:].tolist()
                    chunkPos=chunks_start_stop(shape[0])
                    for minPos,maxPos in chunkPos:
                        if "/chr"+str(chr)+"/"+type in self.file:
                            genoinfo=node.GT_geno[minPos:maxPos,colpos]           
                            sub_rownames,updated_colnames,sub_geno=self.filter_removed_genotypes(minPos,maxPos,genoinfo,node,colpos,[])
                            updated_rownames.append(sub_rownames)
                            updated_geno.append(sub_geno)
     
            except tb.exceptions.NoSuchNodeError:
                pass                              
            except Exception as e:
                print(e)
                pass
        return np.array(updated_rownames),np.array(updated_colnames),np.array(updated_geno)




    def get_geno_by_variant_IDs(self,rowIDs,chr,groupName=""):

        #assume rowIDs are sorted by genome position
     
        group=self.getGroup(chr,groupName)
        rownames=group.rownames[:].tolist()
        colnames=group.colnames[:]
        rowMask=group.rowmask[:]
        sampleMask=group.samplemask[:]
        if len(rowIDs)>0:
            for id in rowIDs:
                try:
                    minPos=rownames.index(id)
                    break
                except ValueError:
                    continue
            for id in reversed(rowIDs):
                try:
                    maxPos=rownames.index(id)
                    break
                except ValueError:
                    continue
        try:
            minPos
            maxPos
            update_rownames=rownames[minPos:maxPos+1]
            
            sub_geno=group.GT_geno[minPos:maxPos+1,:]    
            sub_Mask=group.Mask_geno[minPos:maxPos+1,:]
            sub_geno=np.multiply(sub_geno,sub_Mask)

            update_rowMask=rowMask[minPos:maxPos+1]
            rowMasked=np.where(update_rowMask==True)[0]
            sampleMasked=np.where(sampleMask==True)[0]
            if len(rowMasked)>0:
                update_rownames=np.array(update_rownames)[np.where(update_rowMask==False)[0]]
                sub_geno=np.delete(sub_geno,rowMasked,0)

            if len(sampleMasked)>0:
                colnames=colnames[np.where(sampleMask==False)[0]]
                sub_geno=np.delete(sub_geno,sampleMasked,1)

            
            return np.array(update_rownames),colnames,np.array(sub_geno)
        except NameError:
            env.logger.error("varaintIDs of this gene are not found on this chromosome {}".format(chr))

    
    def get_geno_by_variant_IDs_sample(self,rowIDs,sampleName,chrs):
        NULL_to_0 = env.treat_missing_as_wildtype
        g=set()
        chrs=range(1,23)
        if len(chr)>0:
            chrs=[chr]
        for chr in chrs:
            try:
                node=self.file.get_node("/chr"+str(chr))
                updated_rownames,colnames,subMatrix=self.get_geno_by_variant_IDs(rowIDs,str(chr))
                # print(updated_rownames,colnames,subMatrix)
                colPos=colnames.tolist().index(sampleName)
                GT_geno=subMatrix[:,colPos]
                # print(GT_geno)
                for idx,id in enumerate(updated_rownames):
                    GT=GT_geno[idx]
                    if np.isnan(GT):
                        if NULL_to_0:
                            g.add((id, 0))
                    else:
                        g.add((id, GT))
                # print(g)
            except tb.exceptions.NoSuchNodeError:
                pass                              
            except Exception as e:
                print(e)
                pass
        return g



    def get_geno_by_sample_ID(self,sampleID,type,chrs=range(1,23),groupName=""):
        geno=[]
        for chr in chrs:
            try:
                node=self.file.get_node("/chr"+str(chr))
                colnames=node.colnames[:].tolist()
                colpos=list(map(lambda x:colnames.index(x),colnames))
                shape=node.shape[:].tolist()
                chunkPos=chunks_start_stop(shape[0])
                for startPos,endPos in chunkPos:
                    if "/chr"+str(chr)+"/"+type in self.file:
                        genoinfo=node.GT_geno[startPos:endPos,colpos]           
                        rownames,colnames,genoinfo=self.filter_removed_genotypes(startPos,endPos,genoinfo,node,colpos,[])
                        samplePos=colnames.tolist().index(sampleID)
                        sampleGeno=genoinfo[:,samplePos]
                        for idx,rowname in enumerate(rownames):
                            genotype=sampleGeno[idx]
                            if np.isnan(genotype):
                                genotype=-1
                            geno.append([rowname,genotype])
            except tb.exceptions.NoSuchNodeError:
                pass
        return geno




    def get_noWT_variants(self,samples):
        noWT=[]
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        if len(samples)!=0:
            for chr in chrs:
                try:
                    #haven't dealt with data==-1
                    node=self.file.get_node("/chr"+str(chr))
                    shape=node.shape[:].tolist()
                    samples.sort()
                    colnames=node.colnames[:].tolist()
                    colpos=list(map(lambda x:colnames.index(x),samples))
                    chunkPos=chunks_start_stop(shape[0])
                    for startPos,endPos in chunkPos:
                        if "/chr"+str(chr)+"/GT_geno" in self.file:
                            genoinfo=node.GT_geno[startPos:endPos]       
                            rownames,colnames,genoinfo=self.filter_removed_genotypes(startPos,endPos,genoinfo,node,colpos,[])
                            genoinfo[genoinfo==-1]=0
                            rowsum=np.nansum(genoinfo,axis=1)
                            noWTvariants=rownames[np.where(rowsum>0)].tolist()
                            noWT.extend(noWTvariants)
                        startPos=endPos
                except tb.exceptions.NoSuchNodeError:
                    pass  
                except Exception as e:
                    print(e)
        self.close()
        return noWT



    def filter_removed_genotypes(self,minPos,maxPos,genoinfo,node,colpos,rowpos):
        #assume rowIDs are sorted by genome position
        rownames=node.rownames[minPos:maxPos]
        colnames=node.colnames[:]
        rowMask=node.rowmask[minPos:maxPos]
        # sampleMask=node.samplemask[:]
        try:
            sub_Mask=node.Mask_geno[minPos:maxPos,:]
            sub_geno=np.multiply(genoinfo,sub_Mask)
            rowMasked=np.where(rowMask==True)[0]
            # sampleMasked=np.where(sampleMask==True)[0]  

            
            if len(rowpos)>0:
                rownames=rownames[rowpos]
                sub_geno=sub_geno[rowpos,:]
            colnames=colnames[colpos]
            sub_geno=sub_geno[:,colpos]
    

            if len(rowpos)==0 and len(rowMasked)>0:
                rownames=rownames[np.where(rowMask==False)]
                sub_geno=np.delete(sub_geno,rowMasked,0)

            
            # if len(sampleMasked)>0:
            #     # colnames=colnames[np.where(sampleMask==False)]
            #     # sub_geno=np.delete(sub_geno,sampleMasked,1)
            #     colnames=colnames[colpos]
            #     sub_geno=sub_geno[:,colpos]
    
            return np.array(rownames),colnames,np.array(sub_geno)
        except Exception as e:
            print(e)



    def filter_on_genotypes(self,genotypes,chr,node,field,startPos,endPos,colpos,rowpos):
        genoinfo=None
        if field=="GT_geno":
            genoinfo=node.GT_geno[startPos:endPos,:]
        if field=="DP_geno" and "/chr"+str(chr)+"/DP_geno" in self.file:
            genoinfo=node.DP_geno[startPos:endPos,:]
            genoinfo[genoinfo==-1]=0
        if field=="GQ_geno" and "/chr"+str(chr)+"/GQ_geno" in self.file:
            genoinfo=node.GQ_geno[startPos:endPos,:]
            genoinfo=np.nan_to_num(genoinfo)

        if len(genotypes)>0:
            if type(genotypes) is str:
                genotypes=genotypes.replace("(","").replace(")","")
            else:
                genotypes=genotypes[0]
            if "DP_geno" in genotypes and "/chr"+str(chr)+"/DP_geno" in self.file:
                DP_geno=node.DP_geno[startPos:endPos,:]
                DP_geno[DP_geno==-1]=0
            if "GQ_geno" in genotypes and "/chr"+str(chr)+"/GQ_geno" in self.file:
                GQ_geno=node.GQ_geno[startPos:endPos,:]
                GQ_geno=np.nan_to_num(GQ_geno)
            genoinfo=np.where(eval("~("+genotypes+")"),np.nan,genoinfo)

        rownames,colnames,genoinfo=self.filter_removed_genotypes(startPos,endPos,genoinfo,node,colpos,rowpos)

        return rownames,colnames,genoinfo

    def find_element_in_list(self,element, list_element):
        try:
            index_element = list_element.index(element)
            return index_element
        except ValueError:
            pass


    def get_geno_field_from_HDF5(self,samples,varids,genotypes,fieldSelect,validGenotypeFields,operations):
        vardict={}
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
        # for chr in [1,22]:
            try:
                node=self.file.get_node("/chr"+str(chr))
                shape=node.shape[:].tolist()   
                chunkPos=chunks_start_stop(shape[0])
                samples.sort()
                colnames=node.colnames[:].tolist()
                colpos=list(map(lambda x:colnames.index(x),samples))
                # print(self.fileName,samples,colpos)
                rowpos=[]
                if len(varids)>0:
                    rownames=node.rownames[:].tolist()
                    # rowpos=list(map(lambda x:rownames.index(x),varids))
                    rowpos=[self.find_element_in_list(id,rownames) for id in varids]
                    rowpos=[x for x in rowpos if x is not None]
                    if len(rowpos)==0:
                        continue
                for startPos,endPos in chunkPos:
                    # rownames=node.rownames[startPos:endPos].tolist()              
                    if "/chr"+str(chr)+"/GT_geno" in self.file:    
                        rownames,colnames,genoinfo=self.filter_on_genotypes(genotypes,chr,node,"GT_geno",startPos,endPos,colpos,rowpos)
                    
                        numrow=len(rownames)
                        variants=np.zeros(shape=(numrow,len(validGenotypeFields)+5),dtype=np.int64)   
                        
                        variants[:,3]=np.nansum(~np.isnan(genoinfo),axis=1)
                        variants[:,0]=np.nansum(genoinfo==1,axis=1)
                        variants[:,1]=np.nansum(genoinfo==2,axis=1)
                        variants[:,2]=np.nansum(genoinfo==-1,axis=1) 
                        # print(self.fileName,genoinfo.shape,variants[6,1],np.where(genoinfo[6,:]==2)) 

                        for pos,field in enumerate(validGenotypeFields):
                            if "/chr"+str(chr)+"/"+field in self.file:
                                rownames,colnames,genoinfo=self.filter_on_genotypes(genotypes,chr,node,field,startPos,endPos,colpos,rowpos)
                                # if field=="DP_geno":
                                #     print(self.fileName,genoinfo[0,:100])
                                operation=operations[pos]
                                if operation==0:                      
                                    variants[:,5+pos]=np.nansum(genoinfo,axis=1)
                                if operation==1:
                                    variants[:,5+pos]=np.nansum(genoinfo,axis=1)
                                if operation==2:
                                    variants[:,5+pos]=np.nanmin(genoinfo,axis=1)
                                if operation==3:
                                    variants[:,5+pos]=np.nanmax(genoinfo,axis=1)
                        startPos=endPos    
                        vardict.update(dict(zip(rownames,variants)))
            except tb.exceptions.NoSuchNodeError:
                pass                              
            except Exception as e:
                print(e)
                pass
        return vardict


    def get_genoType_forExport_from_HDF5(self,samples,validGenotypeFields):
        vardict={}
        chrs=["X","Y"]
        chrs.extend(range(1,23))
        for chr in chrs:
        # for chr in [1,22]:
            try:
                node=self.file.get_node("/chr"+str(chr))
                shape=node.shape[:].tolist()          
                chunkPos=chunks_start_stop(shape[0])
                samples.sort()
                colnames=node.colnames[:].tolist()
                colpos=list(map(lambda x:colnames.index(x),samples))
                rownames=node.rownames[:].tolist()
                # print(self.fileName,samples,colpos)
                rowpos=[]
             
                for startPos,endPos in chunkPos:
                    # rownames=node.rownames[startPos:endPos].tolist()              
                    if "/chr"+str(chr)+"/GT_geno" in self.file:    
                        rownames,colnames,genoinfo=self.filter_on_genotypes("",chr,node,"GT_geno",startPos,endPos,colpos,rowpos)
                        numrow,numcol=genoinfo.shape[0],genoinfo.shape[1]
                        info=np.zeros(shape=(numrow,numcol*len(validGenotypeFields)),dtype=float)   
                        for col in range(numcol):
                                info[:,col*len(validGenotypeFields)]=genoinfo[:,col]

                        for pos,field in enumerate(validGenotypeFields):
                            if "/chr"+str(chr)+"/"+field in self.file:
                                rownames,colnames,genoinfo=self.filter_on_genotypes("",chr,node,field,startPos,endPos,colpos,rowpos)
                                for col in range(numcol):
                                    info[:,col*len(validGenotypeFields)+pos]=genoinfo[:,col]
                        startPos=endPos
                   
                        vardict.update(dict(zip(rownames,info)))  
            except tb.exceptions.NoSuchNodeError:
                pass                              
            except Exception as e:
                print(e)
                pass
        return vardict





    


    def get_geno_by_row_pos(self,rowpos,chr,groupName=""):
        pass


    def get_geno_by_variant_ID(self,variantID,chr,groupName=""):
        pass


        
    def compare_HDF5(self,hdf5,chr,groupName=""):
        """This function checks the similarity of two HDF5 files.

            Args:

                - hdf5 (HDF5Engine_access object): another HDF5_access object to compare with
                - chr (string): the chromosome
                - groupName (string): the group name, for example gene name 

        """
        assert sum(self.get_rownames(chr)==hdf5.get_rownames(chr))==len(self.get_rownames(chr)),"rownames not the same"
        assert sum(self.get_data(chr)==hdf5.get_data(chr))==len(self.get_data(chr)),"data not the same" 
        assert sum(self.get_indices(chr)==hdf5.get_indices(chr))==len(self.get_indices(chr)),"indices not the same"             
        assert sum(self.get_indptr(chr)==hdf5.get_indptr(chr))==len(self.get_indptr(chr)),"indptr not the same"       

    
    def show_file_node(self):
        """This function prints the nodes in the file. 

        """
        for node in self.file:
            print(node)



    def get_rownames(self,chr,groupName=""):
        """This function gets variant IDs of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - a ndarray of variant IDs

        """
        group=self.getGroup(chr,groupName)
        return group.rownames[:]

    def get_colnames(self,chr,groupName=""):
        """This function gets sample IDs of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                a ndarray of sample IDs

        """
        group=self.getGroup(chr,groupName)
        return group.colnames[:]

    def get_shape(self,chr,groupName=""):
        """This function gets shape of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                shape of matrix

        """
        group=self.getGroup(chr,groupName)
        return group.shape[:]

   

    def close(self):
        """This function closes db file
        """

        self.file.close()



class AccessEachHDF5(Process):

    def __init__(self,fileName,variantID,result,chr,groupName=""):
        Process.__init__(self)
        self.hdf5=HDF5Engine_access(fileName)
        self.variantID=variantID
        self.chr=chr
        self.groupName=groupName
        self.result=result

    def run(self):
        variant_ID,indices,data=self.hdf5.get_geno_info_by_variant_ID(self.variantID,self.chr,self.groupName)
        print((variant_ID,len(indices)))
        colnames=self.hdf5.get_colnames(self.chr,self.groupName)
        for idx,col in enumerate(indices):
            if not math.isnan(data[idx]):
                self.result[colnames[col]]=int(data[idx])
            else:
                self.result[colnames[col]]=data[idx]
        self.hdf5.close()

class HDF5Engine_access_multi:
    """This is the class for access genotype info stored in multiple HDF5 (under development).

    """
    def __init__(self,files,jobs):
        self.files=files
        self.jobs=jobs


    def get_geno_info_by_variant_ID(self,variantID,chr,groupName=""):
        """This function gets the genotype info of a variant specified by the variant_id stored in rownames. 

            Args:

                - variantID : the variant ID stored in rownames 
                - chr (string): the chromosome  
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                right now, just print a dict of dict[sampleID]=genotype info

        """
        taskQueue=queue.Queue()
        importers=[None]*self.jobs
        result=Manager().dict()
        for file in self.files:
            taskQueue.put(AccessEachHDF5(file,variantID,result,chr,groupName))
        while taskQueue.qsize()>0:
            for i in range(self.jobs):    
                if importers[i] is None or not importers[i].is_alive():
                    task=taskQueue.get()
                    importers[i]=task
                    importers[i].start()           
                    break 
        for worker in importers:
            worker.join() 
        print(result)


if __name__ == '__main__':
    
    # for i in range(1355,1369):
    #     file=HDF5Engine_storage("/Users/jma7/Development/VAT/importTest_13t/tmp_1_2504_genotypes.h5")
    #     file.recover_variant(i,"22","GT")
    #     file.close()

    file=HDF5Engine_storage("/Users/jma7/Development/VAT/genoinfo_allele/tmp_1_491_genotypes.h5")
    file.remove_genotype("(DP<10) & (GQ>10)","22")

    # for i in range(1,3):
    #     file=HDF5Engine_storage("/Users/jma7/Development/VAT/importTest_13t/tmp_1_2504_genotypes.h5")
    #     file.remove_sample(i,"22","GT")
    #     file.close()


    # files=glob.glob("/Users/jma7/Development/VAT/importTest_12t/tmp*genotypes.h5")
    # testHDF=HDF5Engine_access_multi(files,8)
    # testHDF.get_geno_info_by_variant_ID(2,"22")
    # testHDF=HDF5Engine_access("/Users/jma7/Development/VAT/importTest_11t/tmp_1_2504_genotypes.h5")
    # variant_ID,indices,data=testHDF.get_geno_info_by_variant_ID(2,"22")
    # colnames=testHDF.get_colnames("22")
    # result={}
    # for idx,col in enumerate(indices):
    #         if not math.isnan(data[idx]):
    #             result[colnames[col]]=int(data[idx])
    #         else:
    #             result[colnames[col]]=data[idx]
    # print(result)
    
