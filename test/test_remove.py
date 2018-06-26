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
from testUtils import ProcessTestCase


class TestRemove(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        if self.storeMode=="sqlite" and os.path.isfile('TestRemove_sqlite.tar.gz'):
            self.runCmd('vtools admin --load_snapshot TestRemove_sqlite.tar.gz')
        elif self.storeMode=="hdf5" and os.path.isfile("TestRemove_hdf5.tar.gz"):
            self.runCmd('vtools admin --load_snapshot TestRemove_hdf5.tar.gz')
        else:
            self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
            self.runCmd('vtools import vcf/SAMP1.vcf')
            # self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
            self.runCmd('vtools import vcf/input_nogeno.vcf --build hg18 --sample_name input.tsv')
            self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
            self.runCmd('vtools use ann/testNSFP.ann') 
            self.runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
            self.runCmd('vtools select variant --samples "aff=\'1\'" -t unaffected')
            self.runCmd('vtools update CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --from_stat "CEU_cases_num=#(alt)"')
            self.runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno --var_info DP --build hg18')
            if self.storeMode=="sqlite":
                self.runCmd('vtools admin --save_snapshot TestRemove_sqlite.tar.gz "Snapshot for testing command remove"')
            elif self.storeMode=="hdf5":
                self.runCmd('vtools admin --save_snapshot TestRemove_hdf5.tar.gz "Snapshot for testing command remove"')


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
        self.runCmd('vtools select variant -t "NAME WITH #$%"')
        self.assertOutput('vtools show tables', 'NAME WITH #$%', partial=True)
        self.assertSucc('vtools remove tables "NAME WITH #$%"')
        # Removing field CEU_num from variant table CEU
        count1 = len(self.runCmd('vtools show fields -v0', ret='list'))
        self.assertSucc('vtools remove fields CEU_cases_num')
        count2 = len(self.runCmd('vtools show fields -v0', ret='list'))
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
        self.assertProj(numOfSamples= 63)
        self.assertSucc('vtools remove samples "sample_name like \'NA070%\'"')
        self.assertProj(numOfSamples= 60) 
        self.assertSucc('vtools remove samples "BMI > 20"')
        # note that the sample with missing BMI is still there
        self.assertProj(numOfSamples= 28) 
        self.assertSucc('vtools remove samples "BMI is NULL"')
        self.assertProj(numOfSamples= 25) 

        #remove variant
    def testRemoveVar(self):
        self.assertFail('vtools remove variant')
        self.assertFail('vtools remove variants')
        self.assertSucc('vtools remove variants CEU')
        self.assertProj(hasTable='CEU', negate=True)

    def testRemoveFields(self):
        #add a field in the variant table
        self.runCmd('vtools use testNSFP')
        self.runCmd('vtools update variant --set gene_name=testNSFP.genename')
        self.assertOutput('vtools show table variant', 'output/remove_field_before.txt', partial=-3)
        self.runCmd('vtools remove fields CEU_cases_num gene_name DP') 
        self.assertOutput('vtools show table variant', 'output/remove_field_after.txt', partial=-3)

    def testRemovePhenotype(self):
        #remove genotype 
        self.assertFail('vtools remove phenotype')
        self.assertFail('vtools remove phenotypes')
        # return 0 even with incorrect phenotype
        self.assertSucc('vtools remove phenotypes sample_name')
        self.assertSucc('vtools remove phenotypes filename')
        self.assertSucc('vtools remove phenotypes "sex = "F""')
        self.assertSucc('vtools remove phenotypes sex')
        # removing non-existing phenotype should yield just a warning
        self.assertSucc('vtools remove phenotypes non_existing')
        if self.storeMode=="hdf5":
            self.assertOutput('vtools show samples', 'output/remove_phenotype_output.txt')
        elif self.storeMode=="sqlite":
            self.assertOutput('vtools show samples', 'output/remove_phenotype_sqlite.txt')
    
    def testRemoveGenoField(self):
        #runCmd('vtools import vcf/SAMP2.vcf --geno_info DP_geno --var_info DP--build hg18')
        self.maxDiff=None
        if self.storeMode=="sqlite":
            self.assertOutput('vtools show genotypes', 'output/remove_genofield_before_sqlite.txt')
        elif self.storeMode=="hdf5":
            self.assertOutput('vtools show genotypes', 'output/remove_genofield_before_hdf5.txt')
        self.assertFail('vtools remove geno_fields')
        self.assertFail('vtools remove geno_fields variant_id')
        self.assertFail('vtools remove geno_fields gt')
        self.assertSucc('vtools remove geno_fields DP_geno')

        if self.storeMode=="sqlite":
            self.assertOutput('vtools show genotypes', 'output/remove_genofield_after_sqlite.txt')
        elif self.storeMode=="hdf5":
            self.assertOutput('vtools show genotypes', 'output/remove_genofield_after_hdf5.txt')
                
        self.assertFail('vtools remove projects')
        self.assertSucc('vtools remove project')

if __name__ == '__main__':
    unittest.main()
