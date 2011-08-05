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

class TestSelect(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf.gz --build hg18')
        runCmd('vtools import_txt input.tsv -c 1 2 4 5')
        runCmd('vtools import_phenotype phenotype.txt')
        runCmd('vtools import_vcf SAMP1.vcf')
        runCmd('vtools use ./testNSFP.ann')
    def removeProj(self):
        runCmd('vtools remove project')
    def testSelect(self):
        'Test command vtools select'
        self.assertFail('vtools select')
        self.assertSucc('vtools select -h')
        # Variant table non_existing_variant does not exist.
        self.assertFail('vtools select non_existing_variant')
        # Neither --to_table and --output/--count is specified. Nothing to do.
        self.assertFail('vtools select variant \'testNSFP.non_existing_item is not null\'')
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' -t ns')
        # Neither --to_table and --output/--count is specified. Nothing to do.
        self.assertFail('vtools select ns \'genename = "PLEKHN1"\'')
        self.assertSucc('vtools select ns \'genename = "PLEKHN1"\'  -t plekhn1')
        self.assertSucc('vtools select variant "genename = \'PLEKHN1\'" --samples "aff=\'1\' or BMI<20" -t plekhn2')
        self.assertSucc('vtools select variant "genename = \'PLEKHN1\'" --samples "aff=\'1\' and BMI<20" -t plekhn3')
        self.assertSucc('vtools select variant "genename = \'PLEKHN1\'" --samples "aff=\'1\'" "BMI<20" -t plekhn4')
        self.assertSucc('vtools select -c variant')
        self.assertSucc('vtools select variant "genename = \'PLEKHN1\'" --samples \'aff=1\' -t plekhn1_aff')

if __name__ == '__main__':
    unittest.main()