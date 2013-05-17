#!/usr/bin/env python
#
# $File: test_init.py $
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
import shutil
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample, outputOfCmd

class TestInit(ProcessTestCase):
    def testInit(self):
        'Test command vtools init'
        self.assertFail('vtools init')
        self.assertSucc('vtools init test')
        self.assertFail('vtools init test')
        self.assertSucc('vtools init test -f')
    
    def testParent(self):
        'Test command init --parent'
        try:
            os.mkdir('parent')
        except OSError:
            pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB')
        self.assertSucc('vtools init test --parent parent --variants ceu')
        # non-existing genotype field
        self.assertFail('vtools init test --parent parent --variants ceu --genotypes GD>10')
        self.assertEqual(numOfVariant(), 288)
        self.assertEqual(numOfSample(), 61)
        shutil.rmtree('parent')
        
    def testSample(self):
        'Test command init --samples'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --samples "filename like \'%CEU%\'"')
        self.assertEqual(numOfVariant(), 577)
        self.assertEqual(numOfSample(), 60)
        shutil.rmtree('parent')
        #ended __Long__

    def testVariantSample(self):
        'Test command init --variants with --samples'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --variants ceu --samples "filename like \'%CEU%\'"')
        self.assertEqual(numOfVariant(), 288)
        self.assertEqual(numOfSample(), 60)
        shutil.rmtree('parent')

    def testGenotypes(self):
        'Test command init --genotypes'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        #runCmd('vtools import vcf/SAMP1.vcf')
        #runCmd('vtools import --format fmt/genotypes txt/genotypes.txt --build hg18')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        self.assertEqual(numOfVariant('ceu'), 288)
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --variants variant --genotypes GT=1')
        runCmd('vtools phenotype --from_stat "num=#(GT)" "hom=#(hom)" "het=#(het)"')
        #compare the whole output and result table using "file" option, "output" is null and numOfLines=0 
        self.assertOutput('vtools phenotype --output num hom het', '', 0, 'output/CEU_phynotype_het.txt')
        #compare the first 5 lines among the output and result 
        self.assertOutput('vtools phenotype --output num hom het', '', 5, 'output/CEU_phynotype_het.txt')
        #compare the last 4 lines among the output and result 
        self.assertOutput('vtools phenotype --output num hom het', '', -5, 'output/CEU_phynotype_het.txt')
        #compare the last 2 lines among the output and result 
        self.assertOutput('vtools phenotype --output num hom het', '''63	0	63\n40	0	40\n''', -3)         
        shutil.rmtree('parent')


    def testGenotypes_sample(self):
        'Test command init --genotypes with samples option'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import --format fmt/genotypes txt/genotypes.txt --build hg18')
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --variants variant --samples "filename like \'%geno%\'" --genotypes GT=1') 
        self.assertEqual(numOfSample(),49 )
        runCmd('vtools phenotype --from_stat "num=#(GT)" "hom=#(hom)" "het=#(het)"')
        #compare the first three lines of the output and result using "output" option
        self.assertOutput('vtools phenotype --output num hom het', 
                 '''3	0	3\n7	0	7\n7	0	7''', 3)
        #compare the whole output and result table using "file" option, "output" is null and numOfLines=0 
        self.assertOutput('vtools phenotype --output num hom het', '', 0, 'output/genotype_variant_sample_output.txt')
        shutil.rmtree('parent')

    def testChildren(self):
        'Test command init --children'
        #first project
        try:
           os.mkdir('ceu')
        except OSError:
           pass
        runCmd('vtools init ceu -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools select variant \'ref="A"\' -t refA')
        runCmd('vtools select variant \'ref="G"\' -t refG')
        runCmd('vtools phenotype --set "column_A=sample_name"')
        runCmd('vtools update --set "ref1=ref"')
        self.assertEqual(numOfVariant('variant'), 288)
        self.assertEqual(numOfVariant('refA'), 43)
        self.assertEqual(numOfVariant('refG'), 96)
        self.assertEqual(numOfSample(), 60 )
        shutil.move('ceu.proj', 'ceu/ceu.proj')
        shutil.move('ceu_genotype.DB', 'ceu/ceu_genotype.DB')
        try:
           os.mkdir('sam1')
        except OSError:
           pass
        runCmd('vtools init sam1 -f')
        runCmd('vtools import vcf/SAMP1.vcf --build hg18')
        runCmd('vtools select variant \'ref="A"\' -t refA')
        runCmd('vtools select variant \'ref="C"\' -t refC')
        runCmd('vtools phenotype --set "column-a=sample_name"')
        runCmd('vtools phenotype --set "column_B=sample_name"')
        self.assertEqual(numOfVariant(), 289)
        self.assertEqual(numOfVariant('refA'), 58)
        self.assertEqual(numOfVariant('refC'), 85)
        self.assertEqual(numOfSample(), 1)
        shutil.move('sam1.proj', 'sam1/sam1.proj')
        shutil.move('sam1_genotype.DB', 'sam1/sam1_genotype.DB')
        self.assertSucc('vtools init test --children ceu sam1') 
        self.assertEqual(numOfSample(), 61)
        self.assertEqual(numOfVariant(), 577)
        self.assertEqual(numOfVariant('refA'), 101)
        self.assertEqual(numOfVariant('refC'), 85)
        self.assertEqual(numOfVariant('refG'), 96)
        #
        shutil.rmtree('ceu')
        shutil.rmtree('sam1')

if __name__ == '__main__':
    unittest.main()
