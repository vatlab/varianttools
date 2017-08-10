#!/usr/bin/env python
from scipy.sparse import csr_matrix,rand,coo_matrix,csc_matrix
import tables as tb
import numpy as np
import pandas as pd
from numpy import array
import time
import math
import os
import sys

import glob
from multiprocessing import Process,Manager
import queue

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


# class HMatrix(csr_matrix):
#     """This class accept three arrays which represent the matrix (or sparse matrix,check below for example), the shape of the matrix,
#         rownames (variant_ids), colnames (sample names) to create an object.

#             Args:
            
#                - data (list): the genotype value for row i are stored in data[indptr[i]:indptr[i+1]]
#                - indices (list): the column indices for row i are stored in indices[indptr[i]:indptr[i+1]]
#                - indptr (list): the position for each row in the data and indices array
#                - shape (tuple): the shape of the matrix
#                - rownames (list): variant_id for each variant
#                - colnames (list): sample name for each sample

#             An example of indptr,indices and data

#             >>> indptr = np.array([0, 2, 3, 6])
#             >>> indices = np.array([0, 2, 2, 0, 1, 2])
#             >>> data = np.array([1, 2, 3, 4, 5, 6])
#             >>> csr_matrix((data, indices, indptr), shape=(3, 3)).toarray()
#             array([[1, 0, 2],
#                    [0, 0, 3],
#                    [4, 5, 6]])

#     """
#     def __init__(self,data,indices,indptr,shape,rownames,colnames):
#         self.rownames=rownames
#         self.colnames=colnames
#         if len(rownames)==len(indptr):
#             indptr=[0]+indptr
#         csr_matrix.__init__(self,(data,indices,indptr),shape=shape)


#     @property
#     def rownames(self):
#         return self.__rownames

#     @property
#     def colnames(self):
#         return self.__colnames


#     @rownames.setter
#     def rownames(self,rownames):
#         self.__rownames=rownames

#     @colnames.setter
#     def colnames(self,colnames):
#         self.__colnames=colnames


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
        raise NotImplementError()

    def close(self):
        """This function closes db file
        """

        raise NotImplementError()



class HDF5Engine_storage(Base_Storage):
    """This is the class for storage of genotype info into HDF5

    """

    def __init__(self,fileName):
        # print("HDF5 engine started")
        Base_Storage.__init__(self,fileName)
        self.file=tb.open_file(self.dbPath,"a")


    def store(self,hmatrix,chr,groupName=""):

        
        self.file=tb.open_file(self.dbPath,"a")
        if not self.checkGroup(chr,groupName):
            self.store_HDF5(hmatrix,chr,groupName) 
        else:
            self.append_HDF5(hmatrix,chr,groupName)
        self.close()       
 


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
        
        for par in ('data', 'indices', 'indptr', 'shape',"rownames","colnames"):
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
        raise NotImplementError()

    def get_geno_info_by_sample_ID(self,sample_id,chr,groupName):
        """This function gets the genotype info of a sample specified by the sample ID in the colnames. 

            Args:

                - sampleID : the sample ID stored in colnames
                - chr (string): the chromosome 
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict: snpdict[variant_id]=genotype value

        """
        raise NotImplementError()

    def get_geno_info_by_group(self,groupName,chr):
        """This function gets the genotype info of specified group into a dictionary

            Args:

                - chr (string): the chromosome
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict : snpdict[sample_id][variant_id]=genotype value

        """
        raise NotImplementError();

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
        raise NotImplementError()

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

        raise NotImplementError()

    def close(self):
        """This function closes the HDF5 file.

        """
        raise NotImplementError()



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


    def __load_HDF5_by_group(self,chr,groupName=""):
        """This function loads matrix, rownames and colnames of specified group into memory

            Args:

               chr (string): the chromosome 
               groupName (string): the group name, for example gene name


        """
        self.chr=chr
        group=self.getGroup(chr,groupName)
        pars = []
        for par in ('data', 'indices', 'indptr', 'shape','rownames',"colnames"):
            if hasattr(group, par):
                pars.append(getattr(group, par).read())
            else:
                pars.append([])
        # f.close()
        self.data=pars[0]
        self.indices=pars[1]
        self.indptr=pars[2]
        self.shape=pars[3]
        self.rownames=pars[4]
        self.colnames=pars[5]



    def get_geno_info_by_group(self,groupName,chr=None):
       
        if chr is None:
            #raise error
            pass
        if self.chr is None:
            self.chr=chr
            self.__load_HDF5_by_group(chr,groupName)

        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.iteritems():
            snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))
      
        for idx,id in enumerate(self.rownames):
            variant_ID,indices,data=self.get_geno_info_by_row_pos(idx,chr,groupName)
            if len(indices)>0:
                for colidx,samplePos in enumerate(indices):
                    snpdict[self.colnames[samplePos]][id]=(data[colidx],)
    
        return snpdict

    
    def get_geno_info_by_row_pos(self,rowpos,chr,groupName=""):

        # if self.indices is None:
        #     self.load_HDF5_by_chr(chr,groupName)
        # start=self.indptr[rowpos]
        # end=self.indptr[rowpos+1]
        # variant_ID=self.rownames[rowpos]
        # indices=None
        # data=None
 
        # if end==start:
        #     indices=[]
        #     data=[]
        # else:
        #     indices=self.indices[start:end]
        #     data=self.data[start:end]
        # return variant_ID,indices,data

        group=self.getGroup(chr,groupName)
      
        start=group.indptr[rowpos]
        end=group.indptr[rowpos+1]
        variant_ID=group.rownames[rowpos]
        indices=None
        data=None
 
        if end==start:
            indices=[]
            data=[]
        else:
            indices=group.indices[start:end]
            data=group.data[start:end]
        return variant_ID,indices,data


    def get_geno_info_by_variant_IDs(self,rowIDs,chr,groupName=""):
        
        # if not sorted?
        #                  for id in ids:
        #                     try:
        #                         pos=rownames.index(id)
        #                         idPos.append(pos)
        #                     except ValueError:
        #                         continue    
        #                 idPos.sort()

        #assume rowIDs are sorted by genome position
        group=self.getGroup(chr)
        rownames=group.rownames[:].tolist()
        colnames=group.colnames[:]

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
        
        update_rownames=rownames[minPos:maxPos+1]
        sub_indptr=group.indptr[minPos:maxPos+2]
        # print(idPos[0],idPos[-1],len(sub_indptr))
        sub_indices=sub_data=sub_shape=None
        if len(sub_indptr)>0:
            sub_indices=group.indices[sub_indptr[0]:sub_indptr[-1]]
            sub_data=group.data[sub_indptr[0]:sub_indptr[-1]]
            sub_indptr=[sub_indptr[i]-sub_indptr[0] for i in range(1,len(sub_indptr))]
            sub_shape=(len(sub_indptr)-1,group.shape[1])
        return HMatrix(sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames)



    def get_geno_info_by_variant_ID(self,variantID,chr,groupName=""):

        group=self.getGroup(chr)
        rownames=group.rownames[:].tolist()
        try:
            pos=rownames.index(variantID)
            return self.get_geno_info_by_row_pos(pos,chr,groupName)
        except ValueError as e:
            env.logger.error("Variant with ID {} not in the specified group.".format(variantID))


    def get_geno_info_by_sample_ID(self,sampleID,chr,groupName=""):
        
        self.__load_HDF5_by_group(chr,groupName)
        colnames=self.colnames[:].tolist()
        try:
            colPos=colnames.index(sampleID)
        except ValueError:
            env.logger.error("Given sampleID is not in this HDF5 file.")
        matrix=csr_matrix((self.data,self.indices,self.indptr),shape=self.shape)
        col=matrix.getcol(colPos)
        snpdict={}
        for idx,value in enumerate(col.toarray()):
            genotype=value[0]
            if not math.isnan(genotype):
                genotype=int(value[0])
            snpdict[self.rownames[idx]]=genotype
        return snpdict


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

    def get_indptr(self,chr,groupName=""):
        """This function gets indptr array of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - indptr ndarray

        """
        group=self.getGroup(chr,groupName)
        return group.indptr[:]

    def get_indices(self,chr,groupName=""):
        """This function gets indptr array of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - indptr ndarray

        """
        group=self.getGroup(chr,groupName)
        return group.indices[:]



    def get_data(self,chr,groupName=""):
        """This function gets data array of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - data ndarray

        """
        group=self.getGroup(chr,groupName)
        return group.data[:]


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
        print(variant_ID,len(indices))
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
    pass
