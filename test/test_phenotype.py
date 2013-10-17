# !/usr/bin/env python
#
# $File: test_phenotype --from_file.py $
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
from testUtils import ProcessTestCase, runCmd, numOfSample, initTest, output2list, PrettyPrinter

class TestPhenotype(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(3)
        runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno')
    def removeProj(self):
        runCmd('vtools remove project')
    def testImportPhenotype(self):
        'Test command phenotype --from_file'
        # too few arguments
        self.assertFail('vtools phenotype --from_file')
        self.assertSucc('vtools phenotype --from_file -h')
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools phenotype --from_file phenotype/phenotype.txt')
        out1 = output2list('vtools show samples -l -1')
        #the output format was changed, so we reorganize the output and compare
        ori_file = open('phenotype/phenotype.txt', 'r')
        prt = PrettyPrinter('new_file')
        prt.cache(ori_file.readline().replace('sample_name', 'sample_name\tfilename').split('\t'))
        for line in ori_file:
            c1,c2,c3,c4 = line.split('\t')
            prt.cache([c1, 'vcf/CEU.vcf.gz', c2, c3, '.' if c4.strip() == 'None' else c4])
        ori_file.close()
        prt.write()
        self.assertOutput('vtools show samples','', 10, 'new_file')
        os.remove('new_file')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype1.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype2.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype3.txt')
        
    def testImportFields(self):
        'Test command phenotype --from_file FIELD'
        # importing only a few fields, not all fields
        runCmd('vtools phenotype --from_file phenotype/phenotype.txt aff')
        out3 = output2list('vtools show samples', space2tab=True)
        #the output format was changed, so we reorganize the output and compare
        ori_file2 = open('phenotype/phenotype.txt', 'r')
        new_file2 = open('new_file2','w')
        new_file2.write('sample_name\tfilename\taff\n')
        ori_file2.readline()
        for line in ori_file2:
            c1,c2,c3,c4 = line.split('\t')
            line = '\t'.join([c1,'vcf/CEU.vcf.gz',c2])
            new_file2.write(line + '\n')
        ori_file2.close()
        new_file2.close()
        with open('new_file2') as inputfile:
            out4 = ['\t'.join((x.strip().split('\t')[:3])) for x in inputfile] 
        self.assertEqual(out3[:20], out4[:20])
        os.remove('new_file2')
        
    def testImportPhenotypeWithFilename(self):
        'Test command phenotype --from_file'
        # too few arguments
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools phenotype --from_file phenotype/pheno_filename.txt')
        out1 = output2list('vtools show samples -l -1')
        #the output format was changed, so we reorganize the output and compare
        ori_file = open('phenotype/pheno_filename.txt', 'r')
        prt = PrettyPrinter('new_file')
        ori_file = open('phenotype/pheno_filename.txt', 'r')
        for line in ori_file:
            c1,c2,c3,c4,c5 = line.strip().split('\t')
            prt.cache([c2, c1, c3, c4, c5 if c5 != 'None' else '.'])
        ori_file.close()
        prt.write()
        self.assertOutput('vtools show samples','', 0, 'new_file')
        os.remove('new_file')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype1.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype2.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype3.txt')
        
    def testImportFieldsWithFilename(self):
        'Test command phenotype --from_file FIELD'
        # importing only a few fields, not all fields
        runCmd('vtools phenotype --from_file phenotype/pheno_filename.txt aff')
        out3 = output2list('vtools show samples', True)
        #the output format was changed, so we reorganize the output and compare
        ori_file2 = open('phenotype/pheno_filename.txt', 'r')
        prt = PrettyPrinter('new_file2')
        for line in ori_file2:
            c1,c2,c3,c4,c5 = line.split('\t')
            prt.cache([c2,c1,c3,c4,c5])
        ori_file2.close()
        prt.write()
        with open('new_file2') as inputfile:
            out4 = ['\t'.join((x.split()[:3])) for x in inputfile] 
        self.assertEqual(out3, out4)
        os.remove('new_file2')
        
    def testSetPhenotype(self):
        'Test command phenotype --set'
        self.assertFail('vtools phenotype --set')
        self.assertSucc('vtools phenotype --set race=1 --samples \'filename like "%CEU%"\'')
        self.assertFail('vtools phenotype --set race="white" --samples \'filename like "%CEU%"\'')
        # have to use quote to pass the test
        self.assertSucc('vtools phenotype --set \'race="white"\' --samples \'filename like "%CEU%"\'')
        # FIXME the following tests pass but have to verify the output. Will do that manually and add assertEqual
        # total genotypes per individual. apply "count" on sample variant tables

    def testPhenotypeFromStat(self):
        'Test command phenotype --from_stat'
        self.assertFail('vtools phenotype --from_stat')
        self.assertSucc('vtools phenotype --from_stat -h')
        self.assertSucc('vtools phenotype --from_stat "numGeno=count(*)"')
        self.assertSucc("vtools phenotype --from_stat 'validGeno=count(*)' --genotypes 'DP_geno>10'")
        # apply some sqlite functions on sample variant tables to provide useful information for genotype qualities
        self.assertSucc('vtools phenotype --from_stat "meanDP=avg(DP_geno)" "minDP=min(DP_geno)" "maxDP=max(DP_geno)"')
        self.assertSucc("vtools phenotype --from_stat 'wildtype=#(wtGT)' 'mutants=#(mutGT)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)'")

    def testPhenotypeOutput(self):
        'Test command phenotype with --output'
        runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        self.assertSucc('vtools phenotype --output sample_name filename aff sex BMI')
        self.assertSucc('vtools phenotype --set race=1 --samples \'filename like "%CEU%"\' --output sex BMI race')
        self.assertSucc('vtools phenotype --from_stat "numGeno=count(*)" --output sample_name sex numGeno --header sample_name sex numGeno') 
        self.assertSucc('vtools phenotype --output sample_name sex numGeno --genotypes "DP_geno>10" --header sample_name sex numGeno')
        #the options of --samples and genotypes could not be used without the 4 primary options (--set, --from_stat, --from_file and --output.

if __name__ == '__main__':
    unittest.main()
