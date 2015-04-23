#!/usr/bin/env python
#
# $File: test_remove.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd, numOfSample

class TestRemove(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(6)
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        runCmd('vtools select variant --samples "aff=\'1\'" -t unaffected')
        runCmd('vtools update CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --from_stat "CEU_cases_num=#(alt)"')
        runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno --var_info DP --build hg18')

    def testRemove(self):
        'Test command vtools remove'
        self.assertFail('vtools remove')
        self.assertSucc('vtools remove -h')

        #remove tables
        self.assertFail('vtools remove table')
        self.assertFail('vtools remove tables')
        self.assertSucc('vtools remove tables unaffected')
        self.assertFail('vtools show table unaffected')
        #
        # remove table with strange names
        runCmd('vtools select variant -t "NAME WITH #$%"')
        self.assertTrue('NAME WITH #$%' in outputOfCmd('vtools show tables'))
        self.assertSucc('vtools remove tables "NAME WITH #$%"')
        self.assertFalse('NAME WITH #$%' in outputOfCmd('vtools show tables'))
        # Removing field CEU_num from variant table CEU
        count1 = len(outputOfCmd('vtools show fields -v0').split('\n'))
        self.assertSucc('vtools remove fields CEU_cases_num')
        count2 = len(outputOfCmd('vtools show fields -v0').split('\n'))
        self.assertEqual(count1-count2, 1)

        #remove annotation
        self.assertFail('vtools remove annotation')
        self.assertFail('vtools remove annotations')
        self.assertSucc('vtools remove annotations testNSFP')
        self.assertFail('vtools show annotation testNSFP')

    def testRemoveSample(self):
        #remove samples
        self.assertFail('vtools remove sample')
        self.assertFail('vtools remove samples')
        self.assertEqual(numOfSample(), 63)
        self.assertSucc('vtools remove samples "sample_name like \'NA070%\'"')
        self.assertEqual(numOfSample(), 60) 
        self.assertSucc('vtools remove samples "BMI > 20"')
        # note that the sample with missing BMI is still there
        self.assertEqual(numOfSample(), 28) 
        self.assertSucc('vtools remove samples "BMI is NULL"')
        self.assertEqual(numOfSample(), 25) 

        #remove variant
    def testRemoveVar(self):
        self.assertFail('vtools remove variant')
        self.assertFail('vtools remove variants')
        out = outputOfCmd('vtools show tables')
        # take only the first two columns (table and #variants)
        str1 = '\n'.join(['\t'.join(x.split()[:2]) for x in out.split('\n')])
        self.assertEqual(str1, '''table\t#variants\nCEU\t288\nunaffected\t552\nvariant\t1,036\n''')
        self.assertSucc('vtools remove variants CEU')
        out = outputOfCmd('vtools show tables')
        str2 = '\n'.join(['\t'.join(x.split()[:2]) for x in out.split('\n')])
        self.assertEqual(str2, '''table\t#variants\nunaffected\t289\nvariant\t748\n''') 

        #remove fields in the variant table
    def testRemoveFields(self):
        #add a field in the variant table
        runCmd('vtools use testNSFP')
        runCmd('vtools update variant --set gene_name=testNSFP.genename')
        self.assertOutput('vtools show table variant', '''Name:                   variant\nDescription:            Master variant table\nCommand:\nFields:                 variant_id, bin, chr, pos, ref, alt, CEU_cases_num,\n                        DP, gene_name\nNumber of variants:     1036\n''', skip=3)
        runCmd('vtools remove fields CEU_cases_num gene_name DP') 
        self.assertOutput('vtools show table variant', '''Name:                   variant\nDescription:            Master variant table\nCommand:\nFields:                 variant_id, bin, chr, pos, ref, alt\nNumber of variants:     1036\n''', skip=3)

    def testRemovePhenotype(self):
        #remove genotype 
        self.assertFail('vtools remove phenotype')
        self.assertFail('vtools remove phenotypes')
        self.assertFail('vtools remove phenotypes sample_name')
        self.assertFail('vtools remove phenotypes filename')
        self.assertFail('vtools remove phenotypes "sex = "F""')
        self.assertSucc('vtools remove phenotypes sex')
        # removing non-existing phenotype should yield just a warning
        self.assertSucc('vtools remove phenotypes non_existing')
        self.assertOutput('vtools show samples', 'sample_name\tfilename      \taff\tBMI',1) 
    
    def testRemoveGenofield(self):
        #runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno --var_info DP--build hg18')
        self.maxDiff=None
        self.assertOutput('vtools show genotypes','''SAMP2      \tvcf/SAMP2.vcf \t288          \tGT,DP_geno\n''',-2)
        self.assertFail('vtools remove geno_fields')
        self.assertFail('vtools remove geno_fields variant_id')
        self.assertFail('vtools remove geno_fields gt')
        self.assertSucc('vtools remove geno_fields DP_geno')
        self.assertOutput('vtools show genotypes', '''SAMP2      \tvcf/SAMP2.vcf \t288          \tGT\n''',-2)
        self.assertFail('vtools remove projects')
        self.assertSucc('vtools remove project')

if __name__ == '__main__':
    unittest.main()
