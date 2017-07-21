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
        self.file=tb.open_file(fileName,"r")
        self.group=None

    
    def load_HDF5_by_chr(self,chr):
        self.load_HDF5_by_group(chr)


    def setGroup(self,chr,groupName=""):
        self.chr=chr
        # with tb.open_file(self.fileName) as f:
        node="/chr"+chr
        if node in self.file:
            self.group=self.file.get_node(node)
        if len(groupName)>0:
            node="/chr"+chr+"/"+groupName
            if node in self.file:
                self.group=self.file.get_node("/chr"+chr+"/"+groupName)
        return self.group


    def load_HDF5_by_group(self,chr,groupName=""):
        self.chr=chr
        self.setGroup(chr,groupName)
        pars = []
        for par in ('data', 'indices', 'indptr', 'shape','rownames',"colnames"):
            if hasattr(self.group, par):
                pars.append(getattr(self.group, par).read())
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


    def load_all_GT_by_group(self,groupName,chr=None):
        if chr is None:
            pass
        if self.chr is None:
            self.chr=chr
            self.load_HDF5_by_group(chr,groupName)

        snpdict=dict.fromkeys(self.colnames,{})
        for key,value in snpdict.iteritems():
            snpdict[key]=dict.fromkeys(self.rownames.tolist(),(0,))
      
        for idx,id in enumerate(self.rownames):
            variant_ID,indices,data=self.get_GT_by_row_ID(idx)
            if len(indices)>0:
                for colidx,samplePos in enumerate(indices):
                    snpdict[self.colnames[samplePos]][id]=(data[colidx],)
        
        return snpdict

    
    def get_GT_by_row_ID(self,rowID):

        if self.indices is None:
            self.load_HDF5_by_chr(self.chr)
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

    def get_rownames(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        self.group=self.setGroup(chr,groupName)
        return self.group.rownames

    def get_colnames(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        self.group=self.setGroup(chr,groupName)
        return self.group.colnames

    def get_shape(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        self.group=self.setGroup(chr,groupName)
        return self.group.shape

    def get_indptr(self,chr,groupName=""):
        if self.chr is None:
            self.chr=chr
        self.group=self.setGroup(chr,groupName)
        return self.group.indptr





    def load_all_geno_info_by_chr(self,type,chr):
        pass

    def load_geno_info_by_variantID(self,type,variant_ID):
        pass

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