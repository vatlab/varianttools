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
from testUtils import ProcessTestCase, runCmd, numOfSample, numOfVariant

class TestImportVCF(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
    def removeProj(self):
        runCmd('vtools remove project')

    def testImportVCF(self):
        'Test command vtools import_vcf'
        self.assertFail('vtools import_vcf')
        self.assertFail('vtools import_vcf non_existing.vcf')
        # help information
        self.assertSucc('vtools import_vcf -h')
        # no build information, fail
        self.assertFail('vtools import_vcf SAMP1.vcf')
        # specify build information
        self.assertSucc('vtools import_vcf SAMP1.vcf --build hg18')
        self.assertEqual(numOfSample(), 1)
        self.assertSucc('vtools import_vcf SAMP2.vcf')
        self.assertEqual(numOfSample(), 2)
        # file will be ignored if re-imported
        self.assertSucc('vtools import_vcf SAMP1.vcf')
        self.assertEqual(numOfSample(), 2)
        # another sample
        self.assertSucc('vtools import_vcf CEU.vcf.gz')
        self.assertFail('vtools import_vcf CEU.vcf.gz --build hg19')
        self.assertEqual(numOfSample(), 62)
        # file will be ignored if re-imported
        self.assertSucc('vtools import_vcf CEU.vcf.gz')
        self.assertEqual(numOfSample(), 62)

    def testMixedBuild(self):
        'Test importing vcf files with different reference genomes'
        self.assertSucc('vtools import_vcf SAMP1.vcf --build hg18')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 289)
        self.assertSucc('vtools import_vcf var_format.vcf --build hg19')
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 289 + 98)
        self.assertSucc('vtools import_vcf var_format.vcf --build hg19')
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 289 + 98)

if __name__ == '__main__':
    unittest.main()
