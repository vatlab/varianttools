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
import logging
from variant_tools.HDF5_accessor import *


class TestHDF5(ProcessTestCase):

    def __init__(self,methodName='runTest'):
        ProcessTestCase.__init__(self,methodName)
        vcf=VCF('vcf/hdf5_test.vcf')
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
        self.fileName="testHDF5.h5"
        self.rownames=list(self.variants.keys())
        self.colnames=[i+1 for i in range(len(vcf.samples))]


    def test_import_matrix(self):
        if os.path.isfile(self.fileName):
            os.remove(self.fileName)
        hdf5=HDF5Engine_storage(self.fileName)
        hdf5.store_matrix_into_HDF5(self.genotype_matrix,"22",self.rownames,self.colnames)
        self.assertTrue(os.path.isfile(self.fileName))
        hdf5.close()



    def test_import_matrix_as_arrays(self):
        if os.path.isfile(self.fileName):
            os.remove(self.fileName)
        hdf5=HDF5Engine_storage(self.fileName)
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
        hdf5.store_arrays_into_HDF5(data,indices,indptr,shape,self.rownames,self.colnames,"22")
        self.assertTrue(os.path.isfile(self.fileName))

        hdf5.close()
        

        




if __name__=='__main__':
    logging.basicConfig( stream=sys.stderr, level=logging.DEBUG )
    unittest.main()

