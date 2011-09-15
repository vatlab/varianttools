#/usr/bin/env python
#
# $File: test_import_txt.py $
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
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample, outputOfCmd

class TestImportTXT(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testImportTXT(self):
        'Test command import_txt'
        self.assertFail('vtools import_txt')
        self.assertFail('vtools import_txt txt/input.tsv')
        # help information
        self.assertSucc('vtools import_txt -h')
        # no format information, fail
        self.assertFail('vtools import_txt txt/input.tsv')
        # no build information, fail
        self.assertFail('vtools import_txt --format ../input_fmt/basic txt/input.tsv')
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/basic txt/input.tsv')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 338)
        self.assertSucc('vtools import_txt --build hg19 --format fmt/basic_hg19 txt/input.tsv')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 338)
        # test downloading fmt file from the website
        self.assertSucc('vtools import_txt --build hg18 --format basic txt/input.tsv')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 338)
        self.assertFail('vtools import_txt --build hg18 --format ../input_fmt/non_existing_fmt txt/input.tsv')
    
    def testANNOVAR(self):
        'Testing the annovar input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR txt/ANNOVAR.txt')
        # one of the variant cannot be imported.
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 11)
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR txt/ANNOVAR.txt --force --sample_name kaiw' )
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 11)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'kaiw'+'\n')
    
    def testCASAVA18_SNP(self):
        'Testing the CASAVA SNP input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/CASAVA18_snps txt/CASAVA18_SNP.txt')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 20)
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/CASAVA18_snps txt/CASAVA18_SNP.txt --force --sample_name casavasnp')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 20)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'casavasnp'+'\n')
        
    def testCASAVA18_INDEL(self):
        'Testing the CASAVA INDEL input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/CASAVA18_indels txt/CASAVA18_INDEL.txt')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 25)
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/CASAVA18_indels txt/CASAVA18_INDEL.txt --force --sample_name casavaindel')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 25)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'casavaindel'+'\n')

    def testPileup_INDEL(self):
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/pileup_indels txt/pileup.indel')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 30)
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/pileup_indels txt/pileup.indel --force --sample_name pileupindel')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 30)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'pileupindel'+'\n')
    
    def testAddfield(self):
        runCmd('vtools import_vcf vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import_txt --build hg18 --format ../input_fmt/basic txt/input.tsv')
        runCmd('vtools import_vcf vcf/SAMP1.vcf')
        self.assertEqual(numOfSample(), 62)
        self.assertEqual(numOfVariant(), 915)
        self.assertFail('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR_output txt/annovar.txt.exonic_variant_function --update')
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR_output txt/annovar.txt.exonic_variant_function --update variant')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type is not null" -c'), '54'+'\n')
        
    def testUpdate(self):
        runCmd('vtools import_vcf vcf/* --build hg18')
        runCmd('vtools import_vcf liftover hg19')
        runCmd('vtools import_txt --build hg19 --format fmt/dbSNP_hg19validation txt/dbSNP_hg19validation.txt --update variant')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type_dbSNP is not null" -c'), '307'+'\n')
        
if __name__ == '__main__':
    unittest.main()
