#!/usr/bin/env python
#
# $File: test_init.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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
import shutil
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample

class TestInit(ProcessTestCase):
    def testInit(self):
        'Test command vtools init'
        self.assertFail('vtools init')
        self.assertSucc('vtools init test')
        self.assertFail('vtools init test')
        self.assertSucc('vtools init test -f')
    
    def testParent(self):
        'Test command init --parent'
        try:
            os.mkdir('parent')
        except OSError:
            pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB')
        self.assertSucc('vtools init test --parent parent --variants ceu')
        self.assertEqual(numOfVariant(), 288)
        self.assertEqual(numOfSample(), 61)
        shutil.rmtree('parent')
        
        #The following code is wrotten by Long and will be ended __Long__
    def testSample(self):
        'Test command init --samples'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --samples "filename like \'%CEU%\'"')
        self.assertEqual(numOfVariant(), 577)
        self.assertEqual(numOfSample(), 60)
        shutil.rmtree('parent')
        #ended __Long__

    def testVariantSample(self):
        'Test command init --variants with --samples'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --variants ceu --samples "filename like \'%CEU%\'"')
        self.assertEqual(numOfVariant(), 288)
        self.assertEqual(numOfSample(), 60)
        shutil.rmtree('parent')

    #def testGenotypes(self):
        #'Test command init --genotypes'
        #try:
        #   os.mkdir('parent')
        #except OSError:
        #   pass
        #runCmd('vtools init test -f')
        #runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        #runCmd('vtools import vcf/SAMP1.vcf')
        #runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        #self.assertEqual(numOfVariant('ceu'), 288)
        #shutil.move('test.proj', 'parent/test.proj')
        #shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        #self.assertSucc('vtools init test --parent parent --variants ceu --samples "filename like \'%CEU%\'"')
        #self.assertEqual(numOfVariant(), 288)
        #self.assertEqual(numOfSample(), 60)
        #shutil.rmtree('parent')

if __name__ == '__main__':
    unittest.main()
