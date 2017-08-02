from scipy.sparse import csr_matrix,rand,coo_matrix,csc_matrix
import tables as tb
import numpy as np
import pandas as pd
from numpy import array
import time
import math
import os
import sys





class HDF5Engine_storage:
    """This is the class for storage of genotype info into HDF5

    """

    def __init__(self,fileName):
        # print("HDF5 engine started")
        self.file=tb.open_file(fileName,"a")
 



    def store_arrays_into_HDF5(self,data,indices,indptr,shape,rownames,colnames,chr,groupName=""):
        """This function accept three arrays which represent the matrix (or sparse matrix,check below for example), the shape of matrix,
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



            An example of indptr,indices and data

            >>> indptr = np.array([0, 2, 3, 6])
            >>> indices = np.array([0, 2, 2, 0, 1, 2])
            >>> data = np.array([1, 2, 3, 4, 5, 6])
            >>> csr_matrix((data, indices, indptr), shape=(3, 3)).toarray()
            array([[1, 0, 2],
                   [0, 0, 3],
                   [4, 5, 6]])

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
 



    def append_arrays_into_HDF5(self,data,indices,indptr,shape,rownames,chr,groupName=""):
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


    def store_matrix_into_HDF5(self,data,chr,rownames=None,colnames=None,groupName=""):
        """This function accepts a matrix and store the matrix in to HDF5 file. 

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
        filters = tb.Filters(complevel=9, complib='blosc')
        group=self.getGroup(chr,groupName)
        for par in ('data', 'indices', 'indptr', 'shape'):
            try:
                n = getattr(group, par)
                n._f_remove()
            except AttributeError:
                pass
            arr = np.array(getattr(m, par))
            # atom = tb.Atom.from_dtype(arr.dtype)
            atom=tb.Atom.from_dtype(np.dtype(np.int32))
            if (arr is not None):
                ds = self.file.create_earray(group, par, atom, (0,),filters=filters)
                ds.append(arr)
        ds = self.file.create_earray(group, "rownames", atom, (0,),filters=filters)
        ds.append(rownames)
        ds = self.file.create_earray(group, "colnames", atom, (0,),filters=filters)
        ds.append(colnames)


    def checkGroup(self,chr,groupName=""):
        """This function check the existance of a group

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






class HDF5Engine_access:
    """This is the class for access genotype info in HDF5

    """
    def __init__(self,fileName):
        # print("HDF5 engine started")
        self.fileName=fileName
        self.rownames=None
        self.colnames=None
        self.indptr=None
        self.indices=None
        self.data=None
        self.shape=None
        self.chr=None
        self.file=tb.open_file(fileName,"a")
 

    
    def load_HDF5_by_chr(self,chr,groupName=""):
        """This function loads matrix, rownames and colnames of specified chromosome into memory

            Args:

               - chr (string): the chromosome 


        """
        self.load_HDF5_by_group(chr,groupName)


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
        """This function check the existance of a group

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


    def load_HDF5_by_group(self,chr,groupName=""):
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


    def close(self):
        """This function closes the HDF5 file.

        """
        self.file.close()


    def get_geno_info_by_group(self,groupName,chr=None):
        """This function gets the genotype info of specified group into a dictionay

            Args:

                - chr (string): the chromosome which the group is in 
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - dict : snpdict[sample_id][variant_id]=genotype value

        """
        if chr is None:
            #raise error
            pass
        if self.chr is None:
            self.chr=chr
            self.load_HDF5_by_group(chr,groupName)

        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.iteritems():
            snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))
      
        for idx,id in enumerate(self.rownames):
            variant_ID,indices,data=self.get_geno_info_by_row_pos(idx,chr)
            if len(indices)>0:
                for colidx,samplePos in enumerate(indices):
                    snpdict[self.colnames[samplePos]][id]=(data[colidx],)
    
        return snpdict

    
    def get_geno_info_by_row_pos(self,rowpos,chr,groupName=""):
        """This function gets the genotype info of a row specified by the position of the row in the matrix. 

            Args:

                - chr (string): the chromosome which the group is in 
                - groupName (string): the group which variants are in, for example gene name

            Returns:

                - variant_id (string): the variant_id of the variant
                - indices (list): the position of samples
                - data (list): the genotype info

        """

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
        """This function gets a slice of genotype info specified by a list of variant ids. 
            **The position of these variants should be consecutive.**  

            Args:
                
                - rowIDs (list): a list of variant ids. 
                - chr (string): the chromosome 
                - groupName (string): the group name, if not specified, default is genotype

            Returns:

                - sub_data (list): the data array of sub matrix
                - sub_indices (list): the indices array of sub matrix
                - sub_indptr (list): the indptr array of sub matrix
                - sub_shape (tuple): the shape of sub matrix
                - update_rownames (list): the updated rownames of variants in the matrix
                - colnames (list): the sample names of sub matrix

        """
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

        return sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames



    def get_geno_info_by_variant_ID(self,variantID,chr,groupName=""):
        group=self.getGroup(chr)
        rownames=group.rownames[:].tolist()
        try:
            pos=rownames.index(variantID)
            return self.get_geno_info_by_row_pos(pos,chr,groupName)
        except ValueError as e:
            env.logger.error("Variant with ID {} not in the specified group.".format(variantID))



    def get_rownames(self,chr,groupName=""):
        """This function gets rownames of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - a list of rownames

        """
        group=self.getGroup(chr,groupName)
        return group.rownames[:]

    def get_colnames(self,chr,groupName=""):
        """This function gets colnames of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                a list of colnames

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
                - indptr array

        """
        group=self.getGroup(chr,groupName)
        return group.indptr[:]

    def get_data(self,chr,groupName=""):
        """This function gets data array of specified group.
            
            Args:

                - chr (string): the chromosome 
                - groupName (string): the group name, for example gene name

            Return:
                - data array

        """
        group=self.getGroup(chr,groupName)
        return group.data[:]


    def show_file_node(self):
        for node in self.file:
            print(node)


    def get_number_of_samples(self):
        pass

    def get_sample_ids(self):
        pass

    def get_number_of_variants(self):
        pass

    def get_chromosomes(self):
        pass




class HDF5Engine_multi:
    def __init__(self,fileName):
        # print("HDF5 engine started")
        self.fileName=None
        self.m=None
        self.rownames=None
        self.fileName=fileName

    
    def load_HDF5(self,groupName,chr):
        # start_time = time.time()
        group=None
        with tb.open_file(self.fileName) as f:
            node="/chr"+chr
            if node in f:
                group=f.get_node(node)
            if len(groupName)>0:
                node="/chr"+chr+"/"+groupName
                if node in f:
                    group=f.get_node("/chr"+chr+"/"+groupName)
            pars = []
            for par in ('data', 'indices', 'indptr', 'shape','rownames','colnames'):
                if hasattr(group, par):
                    pars.append(getattr(group, par).read())
                else:
                    pars.append([])
        f.close()
        # print(groupName,pars[3])
        self.m = csr_matrix(tuple(pars[:3]), shape=pars[3])
        self.rownames=pars[4]
        self.colnames=pars[5]

    def load_all_GT(self):      

        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.iteritems():
            snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))

        for idx,id in enumerate(self.rownames):
            value=self.m.getrow(idx)
            cols=value.indices
            genotypes=value.data
            for colindex,samplePos in enumerate(cols): 
                snpdict[self.colnames[samplePos]][id]=(genotypes[colindex],)
               
        return snpdict








class HDF5Engine_csc:
    def __init__(self,dbName):
        # print("HDF5 engine started")
        self.dbName = dbName
        self.fileName=None
        self.m=None


    def connect_HDF5(self,db):
        db = os.path.expanduser(db)
        # self.fileName = "/Users/jma7/Development/VAT/chr22_10t/csc_chr22_test.h5"
        self.fileName = "/Users/jma7/Development/VAT/chr22_50t/csc_chr22_test.h5"
        

    def load_HDF5(self):
        # start_time = time.time()
        with tb.open_file(self.fileName) as f:
            pars = []
            for par in ('data', 'indices', 'indptr', 'shape'):
                pars.append(getattr(f.root, '%s_%s' % (self.dbName, par)).read())
        f.close()
        self.m = csc_matrix(tuple(pars[:3]), shape=pars[3])


    def load_genotype_by_variantID(self,ids,sampleID):
        result=self.m.getcol(int(sampleID)-1).toarray()
        data={}
        for id in ids:
            genotype=result[int(id)-1][0]
            if math.isnan(genotype):
                data[int(id)]=(0,)
            else:
                data[int(id)]=(int(genotype),)
        return data



if __name__ == '__main__':
    pass