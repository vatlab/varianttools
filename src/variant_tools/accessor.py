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


# from .importer_vcf_to_hdf5 import _hdf5_setup_datasets,_hdf5_store_chunk,\
# DEFAULT_ALT_NUMBER,DEFAULT_BUFFER_SIZE,DEFAULT_CHUNK_LENGTH,DEFAULT_CHUNK_WIDTH

class HMatrix:
    """This class accepts three arrays which represent the matrix (or sparse matrix,check below for example), the shape of the matrix,
        rownames (variant_ids), colnames (sample names) to create an object to be used as input for storage into HDF5.

            Args:
            
               - data (list): the genotype value for row i are stored in data[indptr[i]:indptr[i+1]]
               - indices (list): the column indices for row i are stored in indices[indptr[i]:indptr[i+1]]
               - indptr (list): the position for each row in the data and indices array
               - shape (tuple): the shape of the matrix
               - rownames (list): variant_id for each variant
               - colnames (list): sample name for each sample

            An example of indptr,indices and data

            >>> indptr = np.array([0, 2, 3, 6])
            >>> indices = np.array([0, 2, 2, 0, 1, 2])
            >>> data = np.array([1, 2, 3, 4, 5, 6])
            >>> csr_matrix((data, indices, indptr), shape=(3, 3)).toarray()
            array([[1, 0, 2],
                   [0, 0, 3],
                   [4, 5, 6]])

    """
    def __init__(self,data,indices,indptr,shape,rownames,colnames):
        self.rownames=rownames
        self.colnames=colnames
        self.indptr=indptr
        self.indices=indices
        self.data=data
        self.shape=shape

    @property
    def rownames(self):
        return self.__rownames

    @property
    def colnames(self):
        return self.__colnames

    @property
    def indptr(self):
        return self.__indptr

    @property
    def indices(self):
        return self.__indices

    @property
    def data(self):
        return self.__data

    @property
    def shape(self):
        return self.__shape

    # @property
    # def rowMask(self):
    #     return self.__rowMask

    @rownames.setter
    def rownames(self,rownames):
        self.__rownames=rownames

    @colnames.setter
    def colnames(self,colnames):
        self.__colnames=colnames

    @indptr.setter
    def indptr(self,indptr):
        self.__indptr=indptr

    @indices.setter
    def indices(self,indices):
        self.__indices=indices

    @data.setter
    def data(self,data):
        self.__data=data

    @shape.setter
    def shape(self,shape):
        self.__shape=shape

    # @rowMask.setter
    # def rowMask(self,rowMask):
    #     self.__rowMask=rowMask




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

class GenoCallData(tb.IsDescription):
    entryMask=tb.BoolCol(pos=1)
    variant_id=tb.Int32Col(dflt=1,pos=2)
    sample_id=tb.Int16Col(dflt=1,pos=3)
    DP=tb.Int8Col(dflt=1,pos=4)
    GQ=tb.Float16Col(dflt=1,pos=5)
    # AD1=tb.Int8Col(dflt=1,pos=6)
    # AD2=tb.Int8Col(dflt=1,pos=7)
    # PL1=tb.Int8Col(dflt=1,pos=8)
    # PL2=tb.Int8Col(dflt=1,pos=9)
    # PL3=tb.Int8Col(dflt=1,pos=10)






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
        
        # if not self.checkGroup(chr,groupName):
        #     self.store_HDF5(data,chr,groupName) 
        # else:
        #     self.append_HDF5(data,chr,groupName)
          

    # def store_table(self,data,tableName,chr="",groupName=""):
    #     filters = tb.Filters(complevel=9, complib='blosc')
    #     if not self.checkGroup(chr,tableName):
    #         group=self.getGroup(chr,tableName)
    #         table=self.file.create_table(group,tableName,GenoCallData,filters=filters)

    #     group=self.getGroup(chr,tableName)
    #     table=group.genoInfo
    #     row=table.row
    #     row["entryMask"]=False
    #     for dataRow in data:
    #         if (len(dataRow)==4):
    #             # for idx,var in enumerate(["variant_id","sample_id","DP","GQ","AD1","AD2","PL1","PL2","PL3"]):
    #             for idx,var in enumerate(["variant_id","sample_id","DP","GQ"]):
    #                 row[var]=dataRow[idx]
    #             row.append()

    #     table.flush()


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
            except tb.exceptions.NoSuchNodeError:
                pass
        return totalNum

    def to_csr_matrix(self,group):
        return csr_matrix((group.data[:],group.indices[:],group.indptr[:]),shape=group.shape[:])

    # def get_sparse_noWT_variants(self,sampleID):
    #     noWT=None
    #     for chr in range(1,23):
    #         try:
    #             #haven't dealt with data==-1
    #             group=self.file.get_node("/chr"+str(chr)+"/GT/")
    #             matrix=self.to_csr_matrix(group)
    #             colnames=group.colnames[:]
    #             rownames=group.rownames[:]
    #             samplePos=np.where(colnames==sampleID)
    #             nonZero=matrix.nonzero()
    #             sampleLoc=np.where(nonZero[1]==samplePos)
    #             rowLoc=nonZero[0][sampleLoc[1]]
    #             noWT=rownames[rowLoc]
    #             # colPos=np.where(indices==samplePos[0][0])
    #             # data=group.data[colPos]
    #             # print(data)
    #             # noNan=rownames[np.where(~np.isnan(data))]
    #             # noNone=rownames[np.where(data!=-1)]
    #             # noWT=np.intersect1d(noNan,noNone)
    #         except Exception as e:
    #             pass
    #     return noWT

   



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
        

    # def remove_genotype(self,cond):
    #     for chr in range(1,23):
    #         genoNode="/chr"+str(chr)+"/genoInfo"
    #         try:
    #             group=self.file.get_node(genoNode)
    #             table=group.genoInfo
    #             for x in table.where(cond):
    #                 x["entryMask"]=True
    #                 x.update()
    #         except:
    #             env.logger.error("The imported VCF file doesn't have DP or GQ value available for chromosome {}.".format(chr))


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
        # group=self.file.get_node("/chr"+chr+"/GT")
        # indptr=group.indptr
        # indices=group.indices
        # rownames=group.rownames[:]
        # data=group.data
        # pre_id=0
        # pos=None
        # sampleList=None
        # for var_id,sample_id in result:
        #     if var_id!=pre_id:
        #         pre_id=var_id
        #         pos=np.where(rownames==var_id)[0][0]
        #         sampleList=indices[indptr[pos]:indptr[pos+1]]
        #         # dataList=data[indptr[pos]:indptr[pos+1]]
        #         # # if len(sampleList)>0:
        #         # #     notNone=np.where(~np.isnan(dataList))[0]
        #         # #     sampleList=notNone
        #         # #     print(var_id,dataList,notNone)
        #     try:
        #         if len(sampleList)>0:
        #             samplePos=np.where(sampleList==sample_id)[0][0]
        #             res=[(x['DP'],x['GQ']) for x in table.where('(variant_id=='+str(var_id)+')&(sample_id=='+str(sample_id)+')')]
        #             print("here",var_id,sample_id,data[indptr[pos]+samplePos],res)
        #         # data[indptr[pos]+samplePos]=np.nan
        #     except Exception as e:
        #             pass

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
        

    def store_HDF5(self,hmatrix,chr,groupName=""):
        """The implementation of store API

            Args:

               - hmatrix (HMatrix): An object contains matrix, rownames,colnames
               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name
        """
        
        filters = tb.Filters(complevel=9, complib='blosc')       
        group=self.getGroup(chr,groupName)
        #check to see if 0 has been added to the start of indptr
        if len(hmatrix.indptr)==len(hmatrix.rownames):
            hmatrix.indptr=[0]+hmatrix.indptr
        
        for par in ('data', 'indices', 'indptr', 'shape',"rownames","colnames","rowMask","sampleMask"):
            arr = None
            atom=tb.Atom.from_dtype(np.dtype(np.int32))
            if (par=='data'):
                arr=np.array(hmatrix.data)
                if groupName=="AD_geno" or groupName=="PL_geno":
                    atom=tb.Atom.from_dtype(np.dtype('S20'))
                else:
                    atom=tb.Atom.from_dtype(np.dtype(np.float64))
            elif (par=='indices'):
                arr=np.array(hmatrix.indices)                
            elif (par=='indptr'):
                arr=np.array(hmatrix.indptr)
            elif (par=='shape'):
                arr=np.array(hmatrix.shape)
            elif(par=="rownames"):
                arr=np.array(hmatrix.rownames)
            elif(par=="colnames"):
                arr=np.array(hmatrix.colnames)
            elif(par=="rowMask"):
                arr=np.zeros(len(hmatrix.rownames),dtype=np.bool)
                atom=tb.Atom.from_dtype(np.dtype(np.bool))
            elif(par=="sampleMask"):
                arr=np.zeros(len(hmatrix.colnames),dtype=np.bool)
                atom=tb.Atom.from_dtype(np.dtype(np.bool))

            ds = self.file.create_earray(group, par, atom, (0,),filters=filters)
            ds.append(arr)


    def append_HDF5(self,hmatrix,chr,groupName=""):
        """This function accepts a HMatrix object and appends it to HDF5 file.
            **The columns of appending matrix should have the exact same samples as exisiting matrix and in the same order.**  

            Args:

               - hmatrix (HMatrix): An object contains matrix, rownames,colnames
               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name

        """
    
        group=self.getGroup(chr,groupName)
        arr = None
        currentStart=group.indptr[-1]

        if hmatrix.indptr is not None:
            hmatrix.indptr=[x+currentStart for x in hmatrix.indptr]
        else:
            hmatrix.indptr=[currentStart]

        if hmatrix.data is not None:
            if groupName=="AD_geno" or groupName=="PL_geno":
                arr=np.array(hmatrix.data,dtype=np.dtype("S20"))
            else:
                arr=np.array(hmatrix.data,dtype=np.dtype(np.float64))
            group.data.append(arr)
        if hmatrix.indices is not None:
            arr=np.array(hmatrix.indices,dtype=np.dtype(np.int32))
            group.indices.append(arr)
        if  hmatrix.indptr is not None:
            arr=np.array(hmatrix.indptr,dtype=np.dtype(np.int32))
            group.indptr.append(arr)
        if  hmatrix.rownames is not None:
            arr=np.array(hmatrix.rownames,dtype=np.dtype(np.int32))
            group.rownames.append(arr)
            arr=np.zeros(len(hmatrix.rownames),dtype=np.dtype(np.bool))
            group.rowMask.append(arr)
        group.shape[1]=hmatrix.shape[1]
        group.shape[0]=len(group.rownames[:])
 



    def __store_arrays_into_HDF5(self,data,indices,indptr,shape,rownames,colnames,chr,groupName=""):
        """This function accepts three arrays which represent the matrix (or sparse matrix,check below for example), the shape of the matrix,
        rownames (variant_ids), colnames (sample names), chromosome and groupName (HDF5 group name)

            Args:
            
               - data (list): the genotype value for row i are stored in data[indptr[i]:indptr[i+1]]
               - indices (list): the column indices for row i are stored in indices[indptr[i]:indptr[i+1]]
               - indptr (list): the position for each row in the data and indices array
               - shape (tuple): the shape of the matrix
               - rownames (list): variant_id for each variant
               - colnames (list): sample name for each sample
               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name
        """
        filters = tb.Filters(complevel=9, complib='blosc')       
        group=self.getGroup(chr,groupName)
        #check to see if 0 has been added to the start of indptr
        if len(indptr)==len(rownames):
            indptr=[0]+indptr
        
        for par in ('data', 'indices', 'indptr', 'shape',"rownames","colnames"):
            arr = None
            atom=tb.Atom.from_dtype(np.dtype(np.int32))
            if (par=='data'):
                arr=np.array(data)
                if groupName=="AD_geno" or groupName=="PL_geno":
                    atom=tb.Atom.from_dtype(np.dtype('S20'))
                else:
                    atom=tb.Atom.from_dtype(np.dtype(np.float64))
            elif (par=='indices'):
                arr=np.array(indices)                
            elif (par=='indptr'):
                arr=np.array(indptr)
            elif (par=='shape'):
                arr=np.array(shape)
            elif(par=="rownames"):
                arr=np.array(rownames)
            elif(par=="colnames"):
                arr=np.array(colnames)

            ds = self.file.create_earray(group, par, atom, (0,),filters=filters)
            ds.append(arr)
           
            # ds = f.create_earray(group, "rownames",  tb.StringAtom(itemsize=200), (0,),filters=filters)
            # ds.append(rownames)
 



    def __append_arrays_into_HDF5(self,data,indices,indptr,shape,rownames,chr,groupName=""):
        """This function appends a new matrix to exisiting matrix stored in HDF5, accepts three arrays as described above,
            shape of matrix, rownames,chromosome and groupName. **The columns of appending matrix should have the exact same samples as 
            exisiting matrix and in the same order.**  

            Args:

               - data (list): the genotype value for row i are stored in data[indptr[i]:indptr[i+1]]
               - indices (list): the column indices for row i are stored in indices[indptr[i]:indptr[i+1]]
               - indptr (list): the position for each row in the data and indices array
               - shape (tuple): the shape of the matrix
               - rownames (list): variant_id for each variant
               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name

        """
    
        group=self.getGroup(chr,groupName)
        arr = None
        currentStart=group.indptr[-1]

        if indptr is not None:
            indptr=[x+currentStart for x in indptr]
        else:
            indptr=[currentStart]

        if data is not None:
            if groupName=="AD_geno" or groupName=="PL_geno":
                arr=np.array(data,dtype=np.dtype("S20"))
            else:
                arr=np.array(data,dtype=np.dtype(np.float64))
            group.data.append(arr)
        if indices is not None:
            arr=np.array(indices,dtype=np.dtype(np.int32))
            group.indices.append(arr)
        if  indptr is not None:
            arr=np.array(indptr,dtype=np.dtype(np.int32))
            group.indptr.append(arr)
        if  rownames is not None:
            arr=np.array(rownames,dtype=np.dtype(np.int32))
            group.rownames.append(arr)
        group.shape[1]=shape[1]
        group.shape[0]=len(group.rownames[:])


    def store_matrix_into_HDF5(self,data,rownames,colnames,chr,groupName=""):
        """This function accepts a matrix and stores the matrix into HDF5 file. 

            Args:

               - data (matrix): a matrix with rows as variants and cols as samples, each cell records the geno info
               - rownames (list): variant_id for each variant
               - colnames (list): sample name for each sample
               - chr (string): the chromosome 
               - groupName (string): the group name, for example gene name

        """

        if (data.__class__ != csr_matrix):
            m=csr_matrix(data.as_matrix())
            rownames=data.index.values
            colnames=data.columns.values
        else:
            m=data            
            if rownames is None:
                raise ValueError("The rownames of sparse matrix can't be None.")
            elif colnames is None:
                raise ValueError("The colnames of sparse matrix can't be None.")
            rownames=np.array(rownames)
            colnames=np.array(colnames)

        msg = "This code only works for csr matrices"
        assert(m.__class__ == csr_matrix), msg
        h5matrix=HMatrix(m.data,m.indices,m.indptr,m.shape,rownames,colnames)
        self.store_HDF5(h5matrix,chr,groupName)


    def append_matrix_into_HDF5(self,data,rownames,chr,groupName=""):
        """This function appends a new matrix to exisiting matrix stored in HDF5, accepts a matrix, rownames,chromosome and groupName. 
            **The columns of appending matrix should have the exact same samples as exisiting matrix and in the same order.**  

            Args:

               - data (list): the genotype value for row i are stored in data[indptr[i]:indptr[i+1]]
               - rownames (list): variant_id for each variant
               - chr (string): the chromosome
               - groupName (string): the group name, for example gene name

        """
        if (data.__class__ != csr_matrix):
            m=csr_matrix(data.as_matrix())
            rownames=data.index.values
        else:
            m=data            
            if rownames is None:
                raise ValueError("The rownames of sparse matrix can't be None.")
            rownames=np.array(rownames)
        h5matrix=HMatrix(m.data,m.indices,m.indptr[1:],m.shape,rownames,[])
        self.append_HDF5(h5matrix,chr,groupName)
 
              


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

        #assume rowIDs are sorted by genome position
     
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

    
    def get_geno_by_variant_IDs_sample(self,rowIDs,sampleName,chr=""):
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
                    rowpos=list(map(lambda x:rownames.index(x),varids))
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
                            # print(self.fileName,colpos,genoinfo)           
                            rownames,colnames,genoinfo=self.filter_removed_genotypes(startPos,endPos,genoinfo,node,colpos,[])
                            # print(self.fileName,samples,colnames)
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


    def get_geno_by_row_pos(self,rowpos,chr,groupName=""):
        pass


    def get_geno_by_variant_ID(self,variantID,chr,groupName=""):
        pass


        
        



    ##used for access sparse GT matrix

    # def __load_HDF5_by_group(self,chr,groupName=""):
    #     """This function loads matrix, rownames and colnames of specified group into memory

    #         Args:

    #            chr (string): the chromosome 
    #            groupName (string): the group name, for example gene name


    #     """
    #     self.chr=chr
    #     group=self.getGroup(chr,groupName)
    #     pars = []
    #     for par in ('data', 'indices', 'indptr', 'shape','rownames',"colnames"):
    #         if hasattr(group, par):
    #             pars.append(getattr(group, par).read())
    #         else:
    #             pars.append([])
    #     # f.close()
    #     self.data=pars[0]
    #     self.indices=pars[1]
    #     self.indptr=pars[2]
    #     self.shape=pars[3]
    #     self.rownames=pars[4]
    #     self.colnames=pars[5]




    # def get_geno_info_by_group(self,groupName,chr=None):
       
    #     if chr is None:
    #         #raise error
    #         pass
    #     if self.chr is None:
    #         self.chr=chr
    #     self.__load_HDF5_by_group(chr,groupName)

    #     snpdict=dict.fromkeys(self.colnames,{})
    #     for key,value in snpdict.items():
    #         snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))
      
    #     for idx,id in enumerate(self.rownames):
    #         variant_ID,indices,data=self.get_geno_info_by_row_pos(idx,chr,groupName)
    #         if len(indices)>0:
    #             for colidx,samplePos in enumerate(indices):
    #                 snpdict[self.colnames[samplePos]][id]=(data[colidx],)
    
    #     return snpdict

    
    # def get_geno_info_by_row_pos(self,rowpos,chr,groupName=""):


    #     group=self.getGroup(chr,groupName)
      
    #     start=group.indptr[rowpos]
    #     end=group.indptr[rowpos+1]
    #     variant_ID=group.rownames[rowpos]
    #     indices=None
    #     data=None
 
    #     if end==start:
    #         indices=[]
    #         data=[]
    #     else:
    #         indices=group.indices[start:end]
    #         data=group.data[start:end]
    #     return variant_ID,indices,data

    # def maskRemovedVariants(self,masked,data,indices,indptr,shape,rownames):
    #     for i in masked:   
    #         n = indptr[i+1] - indptr[i]
    #         # print(n)
    #         if n > 0:
    #             data[indptr[i]:-n] = data[indptr[i+1]:]
    #             data = data[:-n]
    #             indices[indptr[i]:-n] = indices[indptr[i+1]:]
    #             indices = indices[:-n]
    #         indptr[i:-1] = indptr[i+1:]
    #         indptr[i:] -= n
    #         indptr =indptr[:-1]
    #         shape = (shape[0]-1, shape[1])
    #         rownames[i:-1] = rownames[i+1:]
    #         rownames =rownames[:-1]
    #     return data,indices,indptr,shape,rownames

    # def maskRemovedSamples(self,masked,data,indices,indptr,shape,rownames,colnames):
    #     # print(len(data),data)
    #     # print(len(indices),indices)
    #     # print(len(indptr),indptr)
    #     # print(shape)
    #     M=csr_matrix((data,indices,indptr),shape=shape)
    #     C=M.tocoo()
    #     keep = ~np.in1d(C.col, masked)
    #     C.data, C.row, C.col = C.data[keep], C.row[keep], C.col[keep]
    #     C.col -= masked.searchsorted(C.col) 
    #     C._shape = (C.shape[0], C.shape[1] - len(masked))
    #     M=C.tocsr()
    #     for i in masked:
    #         colnames[i:-1] = colnames[i+1:]
    #         colnames =colnames[:-1]
    #     return M.data,M.indices,M.indptr,M.shape,rownames,colnames

    # def maskRemovedGenotypes(self,masked,data,indices,indptr,shape,rownames,colnames):
    #     M=csr_matrix((data,indices,indptr),shape=shape)
    #     C=M.tolil()
    #     rownames=np.array(rownames)
    #     for entry in masked:
    #         rowPos=np.where(rownames==entry[0])
    #         colPos=np.where(colnames==entry[1])
    #         # print(entry,rowPos[0],colPos[0])
    #         if len(rowPos[0]>0) and len(colPos[0])>0:
    #             # print(rowPos[0][0],colPos[0][0])
    #             C[rowPos[0][0],colPos[0][0]]=np.nan
    #     # print(len(masked),count)
    #     M=C.tocsr()
    #     return M.data,M.indices,M.indptr,M.shape


    # def get_geno_info_by_variant_IDs(self,rowIDs,chr,groupName=""):

    #     #assume rowIDs are sorted by genome position
    #     group=self.getGroup(chr,groupName)
    #     rownames=group.rownames[:].tolist()
    #     colnames=group.colnames[:]
    #     rowMask=group.rowMask[:]
    #     sampleMask=group.sampleMask[:]

    #     for id in rowIDs:
    #         try:
    #             minPos=rownames.index(id)
    #             break
    #         except ValueError:
    #             continue
    #     for id in reversed(rowIDs):
    #         try:
    #             maxPos=rownames.index(id)
    #             break
    #         except ValueError:
    #             continue
    #     try:
    #         minPos
    #         maxPos
    #         update_rownames=rownames[minPos:maxPos+1]
    #         sub_indptr=group.indptr[minPos:maxPos+2]

    #         # print(idPos[0],idPos[-1],len(sub_indptr))
    #         sub_indices=sub_data=sub_shape=None
    #         if len(sub_indptr)>0:
    #             sub_indices=group.indices[sub_indptr[0]:sub_indptr[-1]]
    #             sub_data=group.data[sub_indptr[0]:sub_indptr[-1]]
    #             sub_indptr=[sub_indptr[i]-sub_indptr[0] for i in range(0,len(sub_indptr))]
    #             sub_shape=(len(sub_indptr)-1,group.shape[1])

    #         update_rowMask=rowMask[minPos:maxPos+1]
    #         rowMasked=np.where(update_rowMask==True)[0]
    #         sampleMasked=np.where(sampleMask==True)[0]+1
    #         # check removed variants
    #         if len(rowMasked)>0:
    #             sub_data,sub_indices,sub_indptr,sub_shape,update_rownames=self.maskRemovedVariants(rowMasked,sub_data,sub_indices,sub_indptr,sub_shape,update_rownames)
    #             # print("after",len(sub_data),len(sub_indices),len(sub_indptr),sub_shape,len(update_rownames))
    #         #check removed samples
    #         if len(sampleMasked)>0:
    #             # print("before",len(sub_data),len(sub_indices),len(sub_indptr),sub_shape,len(update_rownames),len(colnames))
    #             sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames=self.maskRemovedSamples(sampleMasked,sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames)
    #             # print("after",len(sub_data),len(sub_indices),len(sub_indptr),sub_shape,len(update_rownames),len(colnames))
            
            

    #         return HMatrix(sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames)
    #     except NameError:
    #         env.logger.error("varaintIDs of this gene are not found on this chromosome {}".format(chr))





    # def get_geno_info_by_variant_ID(self,variantID,chr,groupName=""):

    #     group=self.getGroup(chr)
    #     rownames=group.rownames[:].tolist()
    #     try:
    #         pos=rownames.index(variantID)
    #         return self.get_geno_info_by_row_pos(pos,chr,groupName)
    #     except ValueError as e:
    #         env.logger.error("Variant with ID {} not in the specified group.".format(variantID))


    # def get_geno_info_by_sample_ID(self,sampleID,chr,groupName=""):
        
    #     self.__load_HDF5_by_group(chr,groupName)
    #     colnames=self.colnames[:].tolist()
    #     try:
    #         colPos=colnames.index(sampleID)
    #     except ValueError:
    #         env.logger.error("Given sampleID is not in this HDF5 file.")
    #     matrix=csr_matrix((self.data,self.indices,self.indptr),shape=self.shape)
    #     col=matrix.getcol(colPos)
    #     snpdict={}
    #     for idx,value in enumerate(col.toarray()):
    #         genotype=value[0]
    #         if not math.isnan(genotype):
    #             genotype=int(value[0])
    #         snpdict[self.rownames[idx]]=genotype
    #     return snpdict



    # def get_geno_field_from_table(self,sampleID,cond,fields):
    #     for chr in range(1,23):
    #         try:
    #             node=self.file.get_node("/chr"+str(chr))
    #             colnames=node.colnames[:].tolist()
    #             colPos=colnames.index(sampleID)
    #             rownames=node.rownames[:].tolist()
    #             infos=[]

    #             for field in fields:
    #                 if field=="GT":
    #                     field="GT_geno"
    #                 if "/chr"+str(chr)+"/"+field in self.file:
                        
    #                     if field=="GT_geno":
    #                         genoinfo=node.GT_geno
    #                         infos.append(genoinfo)
    #                     elif field=="DP_geno":
    #                         genoinfo=node.DP_geno
    #                         infos.append(genoinfo)
    #                     elif field=="GQ_geno":
    #                         genoinfo=node.GQ_geno
    #                         infos.append(genoinfo)
    #             for rowPos,rowName in enumerate(rownames):
    #                 line=[rowName]
    #                 for info in infos:
    #                     line.append(info[int(rowPos),int(colPos)])
    #                 yield(tuple(line))
        
    #         except Exception as e:
    #             # print(e)
    #             pass

    

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

    # def get_indptr(self,chr,groupName=""):
    #     """This function gets indptr array of specified group.
            
    #         Args:

    #             - chr (string): the chromosome 
    #             - groupName (string): the group name, for example gene name

    #         Return:
    #             - indptr ndarray

    #     """
    #     group=self.getGroup(chr,groupName)
    #     return group.indptr[:]

    # def get_indices(self,chr,groupName=""):
    #     """This function gets indptr array of specified group.
            
    #         Args:

    #             - chr (string): the chromosome 
    #             - groupName (string): the group name, for example gene name

    #         Return:
    #             - indptr ndarray

    #     """
    #     group=self.getGroup(chr,groupName)
    #     return group.indices[:]



    # def get_data(self,chr,groupName=""):
    #     """This function gets data array of specified group.
            
    #         Args:

    #             - chr (string): the chromosome 
    #             - groupName (string): the group name, for example gene name

    #         Return:
    #             - data ndarray

    #     """
    #     group=self.getGroup(chr,groupName)
    #     return group.data[:]


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
    
