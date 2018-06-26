#!/usr/bin/env python
#
# $File: test_func.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org)
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
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase

@unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","HDF5 version is not implemented for this test")
class TestFunc(ProcessTestCase):
    def testRefSequence(self):
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertSucc('vtools output variant chr pos ref alt "ref_sequence(chr, pos)"')
        self.assertSucc('vtools output variant chr pos ref alt "ref_sequence(chr, pos-10, pos+10)"')
        self.assertSucc('''vtools output variant chr pos ref alt "ref_sequence('1', pos-10, pos+10)" ''')

    # def testVcfTrack(self):
    #     'Testing vcf track'
    #     self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
    #     self.assertSucc('vtools show track vcf/CEU.vcf.gz')
    #     self.assertSucc('''vtools output variant "track('vcf/CEU.vcf.gz', 'info')" ''')
    #     self.assertSucc('''vtools output variant "track('vcf/CEU.vcf.gz', 'NA12751.GT')" ''')
    #     self.assertSucc('''vtools output variant "track('vcf/CEU.vcf.gz', 'info.AN')" ''')
    #     self.assertSucc('''vtools output variant "track('vcf/CEU.vcf.gz', 'qual')" ''')

    
    def testGenotype(self):
        'Testing function genotype'
        self.runCmd('vtools import vcf/SAMP1.vcf --build hg19')
        self.runCmd('vtools import vcf/SAMP2.vcf')
        self.assertSucc("vtools output variant chr pos 'genotype()'")
        self.assertSucc('''vtools output variant chr pos "genotype('SAMP1')"''')
        self.assertSucc(r'''vtools output variant chr pos "genotype(\"sample_name like 'SAMP%'\", 'd=,')"''')
        self.assertSucc(r'''vtools output variant chr pos "genotype(\"sample_name like 'SAMP%'\", 'd=,&missing=.')"''')
        # create a sample without GT
        self.runCmd('vtools import vcf/SAMP3_complex_variants.vcf --sample_name SAMP3')
        self.assertSucc("vtools output variant chr pos 'genotype()'")

 
    def testSamples(self):
        'Testing function samples'
        self.runCmd('vtools import vcf/SAMP1.vcf --build hg19')
        self.runCmd('vtools import vcf/SAMP2.vcf')
        self.assertSucc("vtools output variant chr pos 'samples()'")
        self.assertSucc(r'''vtools output variant chr pos "samples(\"sample_filter=sample_name like 'SAMP%'&d=,\")"''')
        
if __name__ == '__main__':
    unittest.main()
