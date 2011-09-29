#!/usr/bin/env python
#
# $File: test_exclude.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, assertOutput, outputOfCmd

class TestExclude(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import_variants --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
        runCmd('vtools import_phenotype phenotype/phenotype.txt')
        runCmd('vtools use ann/testNSFP.ann')
        runCmd('vtools select variant \'testNSFP.chr is not null\' -t ns')
        
    def removeProj(self):
        runCmd('vtools remove project')
        
    def testExclude(self):
        'Test command vtools exclude'
        self.assertFail('vtools exclude')
        self.assertSucc('vtools exclude -h')
        # Variant table non_existing_table does not exist.
        self.assertFail('vtools exclude non_existing_table "aff=\'1\'" -t aff')
        # Cannot overwrite master variant table. Please choose another name for the variant table.
        self.assertFail('vtools exclude ns "sift_score<=0.94" -t variant')
        self.assertSucc('vtools exclude ns "sift_score <= 0.94" -t ns_non_damaging')
        # should have 5 variants
        # select "sift_score >= 0.94" will result in 6 variants
        self.assertOutput('vtools select ns "sift_score > 0.94" -c', 'vtools exclude ns "sift_score <= 0.94" -c')
        self.assertEqual(outputOfCmd('vtools exclude ns "variant_id=604" -c'), '6\n')

if __name__ == '__main__':
    unittest.main()
