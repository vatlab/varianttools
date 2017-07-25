from scipy.sparse import csr_matrix,rand,coo_matrix,csc_matrix
import tables as tb
import numpy as np
import pandas as pd
from numpy import array
import time
import math
import os
import sys
from memory_profiler import profile
import HDF5_storage as storage
from .utils import env





class HDF5Engine_access:
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
 



    def store_arrays_into_HDF5(self,data,indices,indptr,shape,rownames,colnames,chr,groupName=""):
  
        filters = tb.Filters(complevel=9, complib='blosc')       
        group=self.setGroup(chr,groupName)
        
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
    
        group=self.setGroup(chr,groupName)
        arr = None
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



    def store_matrix_into_HDF5(self,data,rownames,chr,groupName=""):

        if (data.__class__ != csr_matrix):
            m=csr_matrix(data.as_matrix())
            rownames=data.index.values
            colnames=data.columns.values
        else:
            m=data
            rownames=np.array(rownames)
        msg = "This code only works for csr matrices"
        assert(m.__class__ == csr_matrix), msg
        filters = tb.Filters(complevel=9, complib='blosc')
        group=self.setGroup(chr,groupName)
        for par in ('data', 'indices', 'indptr', 'shape'):
            try:
                n = getattr(group, par)
                n._f_remove()
            except AttributeError:
                pass
            arr = np.array(getattr(m, par))
            atom = tb.Atom.from_dtype(arr.dtype)
            if (arr is not None):
                ds = f.create_carray(group, par, atom, arr.shape,filters=filters)
                ds[:] = arr
        ds = self.file.create_carray(group, "rownames",  tb.StringAtom(itemsize=200), rownames.shape,filters=filters)
        ds[:] = rownames
        # ds = self.f.create_carray(f.root, name+"_colnames", tb.StringAtom(itemsize=100), colnames.shape,filters=filters)
        # # ds[:] = colnames



    
    def load_HDF5_by_chr(self,chr,groupName=""):
        self.load_HDF5_by_group(chr,groupName)


    def setGroup(self,chr,groupName=""):
        self.chr=chr
        # with tb.open_file(self.fileName) as f:
        node="/chr"+chr
        if node in self.file:
            group=self.file.get_node(node)
        else:
            group=self.file.create_group("/","chr"+chr,"chromosome")
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName

            if node in self.file:
                group=self.file.get_node(node)
            else:
                group=self.file.create_group("/chr"+chr,groupName)
        return group


    def load_HDF5_by_group(self,chr,groupName=""):
        self.chr=chr
        group=self.setGroup(chr,groupName)
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
        self.file.close()


    def get_geno_info_by_group(self,groupName,chr=None):
        if chr is None:
            pass
        if self.chr is None:
            self.chr=chr
            self.load_HDF5_by_group(chr,groupName)

        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.iteritems():
            snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))
      
        for idx,id in enumerate(self.rownames):
            variant_ID,indices,data=self.get_geno_info_by_row_ID(idx)
            if len(indices)>0:
                for colidx,samplePos in enumerate(indices):
                    snpdict[self.colnames[samplePos]][id]=(data[colidx],)
        
        return snpdict

    
    def get_geno_info_by_row_ID(self,rowID,groupName=""):

        if self.indices is None:
            self.load_HDF5_by_chr(self.chr,groupName)
        start=self.indptr[rowID]
        end=self.indptr[rowID+1]
        variant_ID=self.rownames[rowID]
        indices=None
        data=None
 
        if end==start:
            indices=[]
            data=[]
        else:
            indices=self.indices[start:end]
            data=self.data[start:end]
        return variant_ID,indices,data


    def get_GT_by_row_IDs(self,rowIDs,chr,groupName=""):
        # if not sorted?
        #                  for id in ids:
        #                     try:
        #                         pos=rownames.index(id)
        #                         idPos.append(pos)
        #                     except ValueError:
        #                         continue    
        #                 idPos.sort()
        group=self.setGroup(chr)
        rownames=group.rownames[:].tolist()
        colnames=group.colnames[:]

        minPos=rownames.index(rowIDs[0])
        maxPos=rownames.index(rowIDs[-1])
        idPos=rownames[minPos:maxPos+1]
        
        sub_indptr=group.indptr[rowIDs[0]-1:rowIDs[-1]+1]
        sub_indices=group.indices[min(sub_indptr):max(sub_indptr)]
        sub_data=group.data[min(sub_indptr):max(sub_indptr)]
        sub_indptr=[sub_indptr[i]-sub_indptr[0] for i in range(len(sub_indptr))]
        sub_shape=(len(sub_indptr)-1,group.shape[1])
        return sub_data,sub_indices,sub_indptr,sub_shape,rownames,colnames


    def get_rownames(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        group=self.setGroup(chr,groupName)
        return group.rownames

    def get_colnames(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        group=self.setGroup(chr,groupName)
        return group.colnames

    def get_shape(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        group=self.setGroup(chr,groupName)
        return group.shape

    def get_indptr(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        group=self.setGroup(chr,groupName)
        return group.indptr





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





class HDF5Engine:
    def __init__(self,dbName):
        # print("HDF5 engine started")
        self.dbName = dbName
        self.fileName=None
        self.m=None
        self.rownames=None
        self.colnames=None

    def connect_HDF5(self,db):
        db = os.path.expanduser(db)
        self.fileName = db if (db.endswith('.proj') or db.endswith('.h5')) else db.replace(".DB","") + '.h5'
        


    def load_HDF5(self):
        # start_time = time.time()
        with tb.open_file(self.fileName) as f:
            pars = []
            for par in ('data', 'indices', 'indptr', 'shape',"rownames","colnames"):
                pars.append(getattr(f.root, '%s_%s' % (self.dbName, par)).read())
        f.close()
        self.m = csr_matrix(tuple(pars[:3]), shape=pars[3])
        self.rownames=pars[4].astype('int')
        self.colnames=pars[5].astype('str')
        del pars
        # print("--- %s seconds ---" % (time.time() - start_time))


    def binarySearchName(self,nameArray,name):
        start=0
        end=len(nameArray)-1
        found=False
        mid=0
        while (start<=end and not found):
            mid=(end+start)//2
            if (name==nameArray[mid]):
                found=True
            else:
                if (name<nameArray[mid]):
                    end=mid-1
                elif (name>nameArray[mid]):
                    start=mid+1
        if not found:
            mid=-1
        return mid


    def load_genotype_by_variantID(self,ids):
        snpdict={}
        print(len(ids))
        for sample in range(len(self.colnames)):
            snpdict[sample+1]={}
        for id in ids:
            # valindex=self.binarySearchName(self.rownames,id)
            valindex=int(id)-1
            value=self.m.getrow(valindex)
            cols=value.indices
            genotypes=value.data
            for index,pos in enumerate(cols):       
                if math.isnan(genotypes[index]):
                    pass
                else:
                    snpdict[pos+1][int(id)]=(int(genotypes[index]),)
                # print(id,sampleName,value)
        return snpdict

    def representsInt(self,s):
        try: 
            int(s)
            return True
        except ValueError:
            return False

    def load_genotype_by_variantID_samples(self,ids,samples):
        # count=0
        sampleArray=[]
        snpdict={}
        for sample in samples:
            if self.representsInt(sample):
                sample=int(sample)
            sampleArray.append(sample)
            snpdict[sample]={}
        totaltime=0
        for id in ids:
            if self.representsInt(id):
                id=int(id)
            # idIndex=self.binarySearchName(self.rownames,id)
            # env.logger.info(str(idIndex)+" "+str(id))
            idIndex=id-1         
            if idIndex!=-1:
                value=self.m.getrow(idIndex)
                cols=value.indices
                genotypes=value.data
                ## use binary search
                start=time.time()
                for sampleName in sampleArray:
                    samplePos=self.binarySearchName(cols,sampleName-1)  
                    if samplePos!=-1:
                        if math.isnan(genotypes[samplePos]):
                            pass
                        else:
                            genotype=genotypes[samplePos]
                            snpdict[sampleName][id]=(int(genotype),)
                # totaltime=totaltime+time.time()-start           
        # print("get genotype "+str(len(ids))+" "+str(len(samples))+" "+str(totaltime))
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