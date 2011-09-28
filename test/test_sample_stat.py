#!/usr/bin/env python
#
# $File: test_import_sample_stat.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd, output2list

class TestSampleStat(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        
    def removeProj(self):
        runCmd('vtools remove project')
        
    def testSampleStat(self):
        'Test command vtools sample_stat'
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import_phenotype phenotype/phenotype.txt')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        self.assertFail('vtools sample_stat')
        self.assertSucc('vtools sample_stat -h')
        self.assertFail('vtools sample_stat --num --hom --het --other ')
        self.assertFail('vtools sample_stat --num variant')
        self.assertSucc('vtools sample_stat variant --num num --hom hom --het het --other other ')
        total = int(outputOfCmd("vtools execute 'select sum(num) from variant'").split('\n')[0])
        hom = int(outputOfCmd("vtools execute 'select sum(hom) from variant'").split('\n')[0])
        het = int(outputOfCmd("vtools execute 'select sum(het) from variant'").split('\n')[0])
        other = int(outputOfCmd("vtools execute 'select sum(other) from variant'").split('\n')[0])
        self.assertEqual(total, hom*2+het+other)
        self.assertSucc('vtools sample_stat CEU -s "filename like \'%CEU%\'" --num CEU_num')
        self.assertEqual(int(outputOfCmd("vtools execute 'select sum(CEU_num) from CEU'").split('\n')[0]), 6383)
        self.assertSucc('vtools sample_stat CEU --samples "filename like \'%CEU%\'" --num CEU_num --hom CEU_hom --het CEU_het --other CEU_other')
        self.assertSucc('vtools sample_stat CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --het CEU_cases_het')
        self.assertEqual(int(outputOfCmd("vtools execute 'select sum(CEU_cases_het) from CEU'").split('\n')[0]), 1601)
        self.assertSucc('vtools sample_stat CEU -s "filename like \'%CEU%\' and aff=\'1\'" --het CEU_ctrls_het')
        
    def testGenotypeSumstats(self):
        'Test command vtools sample_stat SUMMARYSTATISTIC_FIELD'
        runCmd('vtools import_variants --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        # non-existing field, should fail
        self.assertFail('vtools sample_stat variant --max_GQ max_gq --min_GQ min_gq')
        self.assertSucc('vtools sample_stat variant --max_GQ_INFO max_gq --min_GQ_INFO min_gq --mean_GQ_INFO mean_gq')
        self.assertEqual(output2list('vtools output variant max_gq min_gq mean_gq'), [['NA', 'NA', 'NA'], ['NA', 'NA', 'NA'],  ['NA', 'NA', 'NA'], \
                                                                                                                        ['100', '15', '69.33', ], ['6', '3', '4']])
        self.assertSucc('vtools sample_stat variant --sum_GD_INFO total_dp')
        self.assertEqual(output2list('vtools output variant total_dp'), ['NA', 'NA', 'NA', '60', '7', '4'])
        # ffilter out variants having GQ less than 4,
        # then for each remining variant count the total number of alt genotypes across all samples
        self.assertSucc('vtools sample_stat variant --genotype "GQ_INFO >= 4" --num gq_ge_4')
        self.assertEqual(output2list('vtools output variant gq_ge_4'), ['NA', 'NA', 'NA', '0', '3', '1'])
        
if __name__ == '__main__':
    unittest.main()
