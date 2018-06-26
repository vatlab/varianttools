#!/usr/bin/env python
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
from testUtils import ProcessTestCase


class TestPhenotype(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18 --geno_info DP_geno')
        self.runCmd('vtools import vcf/SAMP1.vcf --geno_info DP_geno')
        self.runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno')

    def testImportPhenotype(self):
        'Test command phenotype --from_file'
        # too few arguments
        self.assertFail('vtools phenotype --from_file')
        self.assertSucc('vtools phenotype --from_file -h')
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools phenotype --from_file phenotype/phenotype.txt')
        self.assertOutput('vtools show samples -l -1', 'output/phenotype_import.txt')
        #the output format was changed, so we reorganize the output and compare
        self.assertSucc('vtools phenotype --from_file phenotype/badphenotype1.txt')
        self.assertSucc('vtools phenotype --from_file phenotype/badphenotype2.txt')
        self.assertSucc('vtools phenotype --from_file phenotype/badphenotype3.txt')
        
    def testImportFields(self):
        'Test command phenotype --from_file FIELD'
        # importing only a few fields, not all fields
        self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt aff')
        self.assertOutput('vtools show samples', 'output/phenotype_fields.txt')
        
    def testImportPhenotypeWithFilename(self):
        'Test command phenotype --from_file'
        # too few arguments
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools phenotype --from_file phenotype/pheno_filename.txt')
        self.assertOutput('vtools show samples -l -1', 'output/phenotype_phenotype_with_filename.txt')
        
    def testImportFieldsWithFilename(self):
        'Test command phenotype --from_file FIELD'
        # importing only a few fields, not all fields
        self.runCmd('vtools phenotype --from_file phenotype/pheno_filename.txt aff')
        self.assertOutput('vtools show samples -l -1', 'output/phenotype_phenotype_with_filename_field.txt')
        
    def testSetPhenotype(self):
        'Test command phenotype --set'
        self.assertFail('vtools phenotype --set')
        self.assertSucc('vtools phenotype --set race=1 --samples \'filename like "%CEU%"\'')
        self.assertFail('vtools phenotype --set race="white" --samples \'filename like "%CEU%"\'')
        # have to use quote to pass the test
        self.assertSucc('vtools phenotype --set \'race="white"\' --samples \'filename like "%CEU%"\'')

    
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
        self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        self.assertSucc('vtools phenotype --output sample_name filename aff sex BMI')
        self.assertSucc('vtools phenotype --set race=1 --samples \'filename like "%CEU%"\' --output sex BMI race')
        self.assertSucc('vtools phenotype --from_stat "numGeno=count(*)" --output sample_name sex numGeno --header sample_name sex numGeno') 
        self.assertSucc('vtools phenotype --output sample_name sex numGeno --genotypes "DP_geno>10" --header sample_name sex numGeno')
        #the options of --samples and genotypes could not be used without the 4 primary options (--set, --from_stat, --from_file and --output.

if __name__ == '__main__':
    unittest.main()
