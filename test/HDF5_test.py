from cyvcf2 import VCF
from scipy.sparse import csr_matrix,rand,coo_matrix,csc_matrix
import tables as tb
import numpy as np
import pandas as pd
from numpy import array
import time
import math   
import os
import sys

import unittest
import subprocess
from testUtils import ProcessTestCase 
from variant_tools.HDF5_accessor import *


class TestHDF5_storage(ProcessTestCase):

    def __init__(self,methodName='runTest'):
        ProcessTestCase.__init__(self,methodName)
        vcf=VCF('vcf/hdf5_test.vcf')
        # vcf=VCF('../../test/test.vcf')
        # vcf=VCF("/Users/jma7/Development/VAT/importTest_10t/firstchr22.vcf")
        self.lines=[]
        self.variants={}
        variant_id=1
        for line in vcf:
            for alt in line.ALT:
                self.variants[variant_id]=[line.CHROM,line.end,line.REF,alt]
                rec=line.gt_types.tolist()
                variant_id+=1
                # rec.extend(line.gt_types.tolist())
                self.lines.append(rec)
        self.genotype=np.array(self.lines)
        self.genotype_matrix=csr_matrix(self.genotype)
        self.store_matrix=csr_matrix(np.array(self.lines[:500]))
        self.append_matrix=csr_matrix(np.array(self.lines[500:]))
        self.matrixFileName="testHDF5_matrix.h5"
        self.arrayFileName="testHDF5_array.h5"
        self.rownames=list(self.variants.keys())
        self.colnames=[i+1 for i in range(len(vcf.samples))]
        self.chr="22"


    def test_import_matrix(self):
        if os.path.isfile(self.matrixFileName):
            os.remove(self.matrixFileName)
        hdf5=HDF5Engine_storage(self.matrixFileName)
        hdf5.store_matrix_into_HDF5(self.genotype_matrix,self.rownames,self.colnames,self.chr)
        self.assertTrue(os.path.isfile(self.matrixFileName))
        hdf5.close()

        if os.path.isfile(self.matrixFileName):
            os.remove(self.matrixFileName)
        hdf5=HDF5Engine_storage(self.matrixFileName)
        hdf5.store_matrix_into_HDF5(self.store_matrix,self.rownames[:500],self.colnames,self.chr)
        hdf5.append_matrix_into_HDF5(self.append_matrix,self.rownames[500:],self.chr)
        self.assertTrue(os.path.isfile(self.matrixFileName))
        hdf5.close()



    def test_import_matrix_as_arrays(self):
        if os.path.isfile(self.arrayFileName):
            os.remove(self.arrayFileName)
        hdf5=HDF5Engine_storage(self.arrayFileName)
        data=indices=indptr=shape=None
        for par in ('data', 'indices', 'indptr', 'shape'):
            arr = np.array(getattr(self.genotype_matrix, par))
            if par=='data':
                data=arr
            elif par=='indices':
                indices=arr
            elif par=='indptr':
                indptr=arr
            elif par=='shape':
                shape=arr
        hdf5.store_arrays_into_HDF5(data,indices,indptr,shape,self.rownames,self.colnames,self.chr)
        self.assertTrue(os.path.isfile(self.arrayFileName))
        hdf5.close()
        matrixFile=HDF5Engine_access(self.matrixFileName)
        arrayFile=HDF5Engine_access(self.arrayFileName)
        matrixFile.compare_HDF5(arrayFile,self.chr)



class TestHDF5_access(ProcessTestCase):

    def __init__(self,methodName='runTest'):
        ProcessTestCase.__init__(self,methodName)
        self.hdf5=HDF5Engine_access("vcf/hdf5_test.h5")
        self.chr="/chr22"

    def test_load_HDF5(self):
        
        self.assertTrue(self.hdf5.checkGroup(self.chr))
        self.assertTrue(len(self.hdf5.get_data(self.chr))==18676)
        self.assertTrue(len(self.hdf5.get_indptr(self.chr))==1009)
        self.assertTrue(len(self.hdf5.get_rownames(self.chr))==1008)
        self.assertTrue(len(self.hdf5.get_colnames(self.chr))==500) 
        self.hdf5.load_HDF5_by_chr(self.chr)
        self.assertTrue(len(self.hdf5.data)==18676)
        self.assertTrue(len(self.hdf5.indices)==18676)
        self.assertTrue(len(self.hdf5.indptr)==1009)
        self.assertTrue(len(self.hdf5.rownames)==1008)
        self.assertTrue(len(self.hdf5.colnames)==500)

    def test_get_genotype_from_HDF5(self):
        variant_ID,indices,data=self.hdf5.get_geno_info_by_row_pos(6,self.chr)
        self.assertTrue(variant_ID==7)
        self.assertTrue(len(indices)==491)
        self.assertTrue(len(data)==491)
        variant_ID,indices,data=self.hdf5.get_geno_info_by_variant_ID(7,self.chr)
        
        self.assertTrue(len(indices)==491)
        self.assertTrue(len(data)==491)
        variant_IDs=[i for i in range(69,152)]
        sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames=self.hdf5.get_geno_info_by_variant_IDs(variant_IDs,self.chr)
        self.assertTrue(len(sub_data)==1705)
        self.assertTrue(len(sub_indices)==1705)
        self.assertTrue(len(sub_indptr)==83)
        self.assertTrue(len(update_rownames)==83)
        self.assertTrue(len(colnames)==500)
        sub_data,sub_indices,sub_indptr,sub_shape,update_rownames,colnames=self.hdf5.get_geno_info_by_variant_IDs([7],self.chr)
        self.assertTrue(len(sub_data)==491)
        self.assertTrue(len(sub_indices)==491)
        self.assertTrue(len(sub_indptr)==1)
        self.assertTrue(update_rownames[0]==7)
        snpdict=self.hdf5.get_geno_info_by_sample_ID(1,self.chr)
        nonzero=0
        for key,value in snpdict.items():
        	if value!=0:
        		nonzero+=1
        self.assertTrue(nonzero==40)

        

if __name__=='__main__':
    unittest.main()

