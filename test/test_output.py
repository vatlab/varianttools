#!/usr/bin/env python
#
# $File: test_output.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd

class TestOutput(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_variants --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        
    def testOutputContents(self):
        'Test command vtools output'
        self.assertFail('vtools output')
        self.assertSucc('vtools output -h')
        # too few arguments
        self.assertFail('vtools output variant')
        self.assertFail('vtools output variant non_existing_field')
        self.assertSucc('vtools output variant chr pos ref alt')
        # this test now fails. Why?
        runCmd('vtools liftover hg19')
        out1 = outputOfCmd('vtools output variant chr pos alt_pos ref alt')
        out2 = outputOfCmd('cat txt/input.tsv')
        self.assertEqual(out1, out2)
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        self.assertSucc('vtools output CEU chr pos ref alt -l -1')
        self.assertFail('vtools output CEU and variant -l 10')
        
    def testOutputExpression(self):
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools sample_stat variant --num num --hom hom --het het --other other')
        self.assertFail('vtools output variant sum(num)')
        self.assertEqual(outputOfCmd('vtools output variant "sum(num)" -v0'), '6383'+'\n')
        self.assertFail('vtools output variant count(1)')
        self.assertEqual(outputOfCmd('vtools output variant "count(1)"'), '626'+'\n')
        

if __name__ == '__main__':
    unittest.main()
