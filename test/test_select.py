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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd

class TestSelect(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(6)
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        runCmd('vtools sample_stat variant --num num --hom hom --het het --other other --depth depth')
        runCmd('vtools sample_stat CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --het CEU_cases_het')
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
        # Neither --to_table and --output/--count is specified. Nothing to do.
        self.assertFail('vtools select variant \'testNSFP.chr is not null\'')
        self.assertOutput("vtools select variant -c", '915\n')
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' -t ns')
        self.assertOutput("vtools execute 'select count(*) from ns'", '7\n')
        self.assertOutput("vtools select ns -c", '7\n')
        
    def testSelectSample(self):
        self.assertOutput("vtools select variant --samples 'filename like \"%input.tsv\"' -c", '338\n')
        self.assertSucc('vtools select ns --samples "filename like \'%input.tsv\'" -t ns_input')
        self.assertOutput("vtools select ns_input -c", '7\n')
        # Existing table ns_input is renamed to ns_input_Aug06_161348. The command below is equivalent to the former two commands.
        self.assertSucc('vtools select variant "testNSFP.chr is not null" --samples "filename like \'%input.tsv\'" -t ns_input')
        self.assertSucc('vtools select ns \'genename = "PLEKHN1"\'  -t plekhn1')
        self.assertOutput("vtools select plekhn1 -c", '6\n')
        self.assertSucc('vtools select plekhn1 "polyphen2_score>0.9 and sift_score>0.9" -t d_plekhn1')
        self.assertOutput("vtools select d_plekhn1 -c", '5\n')
        self.assertSucc('vtools select variant "testNSFP.chr is not null" --samples "aff=1" -t ns_aff')
        self.assertOutput("vtools select ns_aff -c", '0\n')
        self.assertSucc('vtools select variant --samples "aff=\'1\' or BMI<20" -t ns2')
        self.assertOutput("vtools select ns2 -c", '563\n')
        self.assertOutput("vtools execute 'select count(*) from sample where aff=1 or BMI<20'", '41\n')
        self.assertSucc('vtools select variant --samples "aff=\'1\' and BMI<20" -t ns3')
        self.assertOutput("vtools select ns3 -c", '217\n')
        self.assertOutput("vtools execute 'select count(*) from sample where aff=1 and BMI<20'", '10\n')
        self.assertSucc('vtools select variant --samples "aff=\'1\'" "BMI<20" -t ns3')
        self.assertOutput("vtools select ns3 -c", '217\n')
        self.assertOutput("vtools execute 'select count(*) from sample where aff=1 and BMI<20'", '10\n')
        self.assertSucc('vtools select variant "testNSFP.chr is not null" "genename=\'PLEKHN1\'" "polyphen2_score>0.9 or sift_score>0.9" -t d_plekhn1')
        self.assertSucc('vtools select variant --samples "sample_name like \'NA0%\'" -t NA0')
        self.assertOutput("vtools select NA0 -c", '223\n')
        self.assertOutput("vtools execute 'select count(*) from sample where sample_name like \"NA0%\"'", '9\n')
        self.assertSucc('vtools select CEU -s "BMI<18.5" -t Underweight')
    
    def testSelectFixme(self):
        initTest(2)
        # no sample name for input.tsv
        runCmd('vtools import_txt --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        self.assertOutput('vtools select variant --samples \'filename like "%CEU.vcf.gz" \' -c', '288\n')
        
if __name__ == '__main__':
    unittest.main()
