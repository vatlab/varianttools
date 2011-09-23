#!/usr/bin/env python
#
# $File: test_import_variants.py $
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
from testUtils import ProcessTestCase, runCmd, numOfSample, numOfVariant, outputOfCmd

class TestImportVCF(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
    def removeProj(self):
        runCmd('vtools remove project')

    def testAddfield(self):
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import_variants --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        runCmd('vtools import_variants vcf/SAMP1.vcf')
        self.assertEqual(numOfSample(), 61)
        self.assertEqual(numOfVariant(), 915)
        self.assertFail('vtools import_variants --build hg18 --format ../input_fmt/ANNOVAR_output txt/annovar.txt.exonic_variant_function --update')
        self.assertSucc('vtools import_variants --build hg18 --format ../input_fmt/ANNOVAR_output txt/annovar.txt.exonic_variant_function --update variant -v2')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type is not null" -c'), '78'+'\n')
        
    def testUpdate(self):
        runCmd('vtools import_variants vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import_variants vcf/SAMP1.vcf --build hg18')
        runCmd('vtools import_variants --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        # no hg19
        self.assertFail('vtools import_variants --build hg19 --format fmt/dbSNP_hg19validation txt/dbSNP_hg19validation.txt --update variant')
        runCmd('vtools liftover hg19')
        # a variant can have multiple entries -- thus 175 variants in variant table are updated from the 198 imported records
        self.assertSucc('vtools import_variants --build hg19 --format fmt/dbSNP_hg19validation txt/dbSNP_hg19validation.txt --update variant')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type_dbSNP is not null" -c'), '175'+'\n')
        self.assertOutput("vtools select variant alt_pos=753405 -o chr pos mut_type_dbSNP validation", "1\t743268\tuntranslated-5\tby-cluster,by-1000genomes\n") 
        
if __name__ == '__main__':
    unittest.main()
