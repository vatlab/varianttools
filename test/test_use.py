#!/usr/bin/env python
#
# $File: test_import_vcf.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
from testUtils import ProcessTestCase, runCmd

class TestUse(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf.gz --build hg18')
        runCmd('vtools import_txt input.tsv -c 1 2 4 5 --zero')
        runCmd('vtools import_phenotype phenotype.txt')
        runCmd('vtools import_vcf SAMP1.vcf')
    
    def removeProj(self):
        runCmd('vtools remove project')

    def testUse(self):
        'Test command vtools use'
        self.assertFail('vtools use')
        self.assertSucc('vtools use -h')
        self.assertFail('vtools use ../non_existing_file/testNSFP.ann')
        self.assertFail('vtools use ./non_existing_file.ann')
        self.assertSucc('vtools use ./testNSFP.ann')
        self.assertFail('vtools use ../non_existing_file/testNSFP.DB')
        self.assertSucc('vtools use ./testNSFP.DB')
        runCmd('vtools show fields')

    def testThousandGenomes(self):
        'Test variants in thousand genomes'
        # no hg19
        self.assertFail('vtools use testThousandGenomes.ann --files testThousandGenomes.vcf.head')
        # liftover
        self.assertSucc('vtools liftover hg19')
        # ok now
        self.assertSucc('vtools use testThousandGenomes.ann --files testThousandGenomes.vcf.head')
        # 137 lines, 9 with two alt
        self.assertOutput('vtools execute "SELECT COUNT(1) FROM testThousandGenomes.testThousandGenomes;"', '146\n')
        # test the handling of 10327   CT   C,CA
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=10328;"',
            '1, 10328, T, -\n1, 10328, T, A\n')
        # test the handling of 83934	rs59235392	AG	A,AGAAA	.
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=83935;"',
            '1, 83935, G, -\n')
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=83936;"',
            '1, 83936, -, AAA\n')
        #
        # Now, let us import vcf regularly.
        self.assertSucc('vtools init test --force')
        self.assertSucc('vtools import_vcf testThousandGenomes.vcf.head')
        self.assertSucc('vtools use testThousandGenomes')
        # do we have all the variants matched up?
        self.assertOutput('vtools select variant -c', '146\n')
        self.assertOutput('vtools select variant "testThousandGenomes.chr is not NULL" -c', '146\n')
        
if __name__ == '__main__':
    unittest.main()
