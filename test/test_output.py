#!/usr/bin/env python
#
# $File: test_output.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd, output2list

class TestOutput(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18')
        
    def testOutputContents(self):
        'Test command vtools output'
        self.assertFail('vtools output')
        self.assertSucc('vtools output -h')
        # too few arguments
        self.assertFail('vtools output variant')
        self.assertFail('vtools output variant non_existing_field')
        self.assertSucc('vtools output variant chr pos ref alt')
        runCmd('vtools liftover hg19')
        out1 = outputOfCmd('vtools output variant chr pos alt_pos ref alt -d"\t"')
        out2 = outputOfCmd('cat txt/input.tsv')
        self.assertEqual(out1, out2)
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        self.assertSucc('vtools output CEU chr pos ref alt -l -1')
        self.assertFail('vtools output CEU and variant -l 10')
        # test output of table with non-ascii name
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t "08#x"')
        self.assertSucc('vtools output "08#x" chr pos ref alt -l -1')


    def testAllOutput(self):
        'Test option --all of command output'
        runCmd('vtools use refGene-hg18_20110909')
        out0 = outputOfCmd('vtools output variant chr ref').split('\n')
        out1 = outputOfCmd('vtools output variant chr ref refGene.name').split('\n')
        out2 = outputOfCmd('vtools output variant chr ref refGene.name --all').split('\n')
        self.assertEqual(len(set(out0)), 6)
        self.assertEqual(len(set(out1)), 33)
        self.assertEqual(len(set(out2)), 38)
        self.assertEqual(len(out0), 339)
        self.assertEqual(len(out1), 339)
        self.assertEqual(len(out2), 407)
        
    def testOutputExpression(self):
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools update variant --from_stat "num=#(alt)" "hom=#(hom)" "het=#(het)" "other=#(other)"')
        self.assertFail('vtools output variant sum(num)')
        self.assertSucc('vtools output variant alt "sum(num)" --group_by alt')
        self.assertSucc('vtools output variant alt "sum(num)" --group_by alt --header alt sum') 
        self.assertSucc('vtools output variant alt "sum(num)" --group_by alt --delimiter ","')  
        self.assertEqual(outputOfCmd('vtools output variant "sum(num)" -v0'), '6383'+'\n')
        self.assertFail('vtools output variant count(1)')
        self.assertEqual(outputOfCmd('vtools output variant "count(1)"'), '626'+'\n')
        self.assertSucc('vtools output variant chr pos ref alt num --order_by num')
        self.assertSucc('vtools output variant chr pos ref alt num --build hg18')
        self.assertSucc('vtools liftover hg19 --flip')
        self.assertSucc('vtools output variant chr pos ref alt num --build hg19')
        self.assertSucc('vtools output variant chr pos ref alt num --header sum_of_num --order_by num')
        self.assertEqual(output2list('vtools output variant num --order_by num')[-10:], ['110', '110', '110', '113', '114', '119', '119', '120', '120', '120'])
        

if __name__ == '__main__':
    unittest.main()
