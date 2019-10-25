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
from testUtils import ProcessTestCase

class TestInit(ProcessTestCase):
    def testInit(self):
        'Test command vtools init'
        # Fail because of no project name
        self.assertFail('vtools init')
        # Fail because an project already exists
        self.assertFail('vtools init test')
        # use -f to forcefully create a project
        self.assertSucc('vtools init test -f')
        # can specify build
        self.assertSucc('vtools init test --build hg19 -f')
    

    def testInitFromParentalProject(self):
        'Test command init --parent (create a project from a parent project)'
        try:
            os.mkdir('parent')
        except OSError:
            pass
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('''vtools select variant --samples "sample_name like 'NA12%'" -t na12''')
        self.assertProj(numOfVariants={'variant': 288, 'na12': 280},
            numOfGenotype={1: 287, 2: 287})
        # move the project to parent directory
        shutil.move('test.proj', 'parent/test.proj')
        if self.storeMode=="sqlite":
            shutil.move('test_genotype.DB', 'parent/test_genotype.DB')
        elif self.storeMode=="hdf5":
            hdf5Files=glob.glob("tmp*h5")
            for hdf5File in hdf5Files:
                shutil.move(hdf5File, 'parent/')

        # create a project with parent project parent
        self.assertSucc('vtools init test --parent parent --variants na12 --store '+self.storeMode)
        self.assertProj(numOfVariants={'variant':280, 'na12': 280}, numOfSamples=60,
            numOfGenotype={1: 279, 2: 279})
        # non-existing genotype field
        self.assertFail('vtools init test --parent parent --variants na12 --genotypes GD>10')
        # create project with only homozygous genotype
        self.assertSucc('vtools init test --parent parent --variants na12 --genotypes GT=1 -f --store '+self.storeMode)
        self.assertProj(numOfVariants={'variant':280, 'na12': 280}, numOfSamples=60,
            numOfGenotype={1: 30, 2: 64}, genotype={1: [1]*30, 2: [1]*64})
        # init with selected samples
        self.assertSucc('''vtools init test --parent parent --samples "sample_name like 'NA1%'" -f  --store '''+ self.storeMode)
        self.assertProj(numOfVariants= 288, numOfSamples=51)

    def testGenotypes_sample(self):
        'Test command init --genotypes with samples option'
        try:
           os.mkdir('parent')
        except OSError:
           pass
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools import --format fmt/genotypes txt/genotypes.txt --build hg18')
        shutil.move('test.proj', 'parent/test.proj')
        shutil.move('test_genotype.DB', 'parent/test_genotype.DB') 
        self.assertSucc('vtools init test --parent parent --variants variant --samples "filename like \'%geno%\'" --genotypes GT=1 --store '+self.storeMode) 
        self.assertProj(numOfSamples=49 )
        self.runCmd('vtools phenotype --from_stat "num=#(GT)" "hom=#(hom)" "het=#(het)"')
        #compare the first three lines of the output and result using "output" option
        self.assertOutput('vtools phenotype --output num hom het', 
                 '''3	0	3\n7	0	7\n7	0	7''', 3)
        #compare the whole output and result table using "file" option, "output" is null and numOfLines=0 
        self.assertOutput('vtools phenotype --output num hom het', 'output/genotype_variant_sample_output.txt')
        shutil.rmtree('parent')

    def testChildren(self):
        'Test command init --children'
        #first project
        try:
           os.mkdir('ceu')
        except OSError:
           pass
        self.runCmd('vtools init ceu -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools select variant \'ref="A"\' -t refA')
        self.runCmd('vtools select variant \'ref="G"\' -t refG')
        self.runCmd('vtools phenotype --set "column_A=sample_name"')
        self.runCmd('vtools update variant --set "ref1=ref"')
        self.assertProj(numOfVariants={'variant': 288, 'refA': 43, 'refG': 96}, numOfSamples=60)
        shutil.move('ceu.proj', 'ceu/ceu.proj')
        if self.storeMode=="sqlite":
            shutil.move('ceu_genotype.DB', 'ceu/ceu_genotype.DB')
        try:
           os.mkdir('sam1')
        except OSError:
           pass
        self.runCmd('vtools init sam1 -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/SAMP1.vcf --build hg18')
        self.runCmd('vtools select variant \'ref="A"\' -t refA')
        self.runCmd('vtools select variant \'ref="C"\' -t refC')
        self.runCmd('vtools phenotype --set "column_A=sample_name"')
        self.runCmd('vtools phenotype --set "column_B=sample_name"')
        self.assertProj(numOfVariants={'variant': 289, 'refA': 58, 'refC': 85}, numOfSamples=1)
        shutil.move('sam1.proj', 'sam1/sam1.proj')
        if self.storeMode=="sqlite":
            shutil.move('sam1_genotype.DB', 'sam1/sam1_genotype.DB')
        self.assertSucc('vtools init test --children ceu sam1 --store '+self.storeMode) 
        self.assertProj(numOfVariants={'variant': 577, 'refA': 101, 'refC': 85, 'refG': 96}, numOfSamples=61)
        #
        shutil.rmtree('ceu')
        shutil.rmtree('sam1')

if __name__ == '__main__':
    unittest.main()
