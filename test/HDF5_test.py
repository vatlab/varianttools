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



class TestHDF5(ProcessTestCase):

    def testPrepare(self):
        vcf=VCF('vcf/hdf5_test.vcf')
        # vcf=VCF("/Users/jma7/Development/VAT/importTest_10t/firstchr22.vcf")
        # log= logging.getLogger( "TestHDF5.testPrepare" )
        log= logging.getLogger(__name__)
        for line in vcf:
            for alt in line.ALT:
                rec=[line.CHROM,line.end,line.REF,alt]
                rec.extend(line.gt_types.tolist())
                log.debug(rec)




if __name__=='__main__':
    logging.basicConfig( stream=sys.stderr, level=logging.DEBUG )
    unittest.main()

