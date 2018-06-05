#!/usr/bin/env python
#
# $File: test_use.py $
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

class TestUse(ProcessTestCase):

    def testUse(self):
        'Test command vtools use'
        self.assertFail('vtools use')
        self.assertSucc('vtools use -h')
        self.assertFail('vtools use non_existing_file.ann')
        self.assertSucc('vtools use ann/testNSFP.ann')
        self.assertSucc('vtools use ann/testNSFP.ann --files ann/testNSFP.zip')
        self.assertSucc('vtools use ann/testNSFP.ann --files ann/testNSFP.zip --rebuild')
        self.assertFail('vtools use ann/testNSFP.ann --files ann/non_existing_file.zip')
        self.assertSucc('vtools use ann/testNSFP.DB.gz')

    def testThousandGenomes(self):
        'Test variants in thousand genomes'
        if os.path.isfile('TestUse.tar.gz'):
            self.runCmd('vtools admin --load_snapshot TestUse.tar.gz')
        else:
            self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
            self.runCmd('vtools import vcf/SAMP1.vcf')
            self.runCmd('vtools import vcf/input_nogeno.vcf --build hg18 --sample_name input.tsv')
            # self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
            self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
            self.runCmd('vtools admin --save_snapshot TestUse.tar.gz "Snapshot of project to test command use"')
        # no hg19
        self.assertFail('vtools use ann/testThousandGenomes.ann --files ann/testThousandGenomes.vcf.head')
        # liftover
        self.assertSucc('vtools liftover hg19')
        # ok now
        self.assertSucc('vtools use ann/testThousandGenomes.ann --files ann/testThousandGenomes.vcf.head')
        # 137 lines, 9 with two alt
        self.assertOutput('vtools execute "SELECT COUNT(1) FROM testThousandGenomes.testThousandGenomes;"', '146\n')
        # test the handling of 10327   CT   C,CA
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=10328;"',
            '1\t10328\tT\t-\n1\t10328\tT\tA\n')
        # test the handling of 83934	rs59235392	AG	A,AGAAA	.
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=83935;"',
            '1\t83935\tG\t-\n')
        self.assertOutput('vtools execute "SELECT chr, pos, ref, alt FROM testThousandGenomes.testThousandGenomes WHERE pos=83936;"',
            '1\t83936\t-\tAAA\n')
        #
        # Now, let us import vcf regularly.
        self.runCmd('vtools init test --force')
        self.assertSucc('vtools import --format vcf ann/testThousandGenomes.vcf.head --build hg19')
        self.assertSucc('vtools use ann/testThousandGenomes.ann')
        # do we have all the variants matched up?
        self.assertOutput('vtools select variant -c', '145')
        self.assertOutput('vtools select variant "testThousandGenomes.chr is not NULL" -c', '145')

    def testESP(self):
        self.runCmd('vtools import --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        self.assertSucc('vtools use ESP')
        self.assertSucc('vtools show annotation ESP')
        self.assertOutput('vtools execute "select sample_name from sample"', 'WHISP:D967-33\nWHISP:D226958-47\nWHISP:D264508-52\nWHISP:D7476-42\n')
        self.assertOutput('vtools output variant variant_id ref alt DP MQ ANNO SVM --header id ref alt DP MQ ANNO SVM -d"\t"', 'output/evsVariantTest.txt')

    def testNSFP(self):
        'Test variants in dbNSFP'
        if os.path.isfile('TestUse.tar.gz'):
            self.runCmd('vtools admin --load_snapshot TestUse.tar.gz')
        else:
            self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
            self.runCmd('vtools import vcf/SAMP1.vcf')
            self.runCmd('vtools import vcf/input_nogeno.vcf --build hg18 --sample_name input.tsv')
            # self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
            self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
            self.runCmd('vtools admin --save_snapshot TestUse.tar.gz "Snapshot of project to test command use"')
        self.assertSucc('vtools use ann/testNSFP.ann --files ann/testNSFP.zip')
        # see if YRI=10/118 is correctly extracted
        self.assertOutput('''vtools execute "SELECT A_freq, C_freq FROM testNSFP.testNSFP WHERE hg18pos=898186 AND alt='A';"''', '2.0\t0\n')
        self.assertOutput('''vtools execute "SELECT A_freq, C_freq FROM testNSFP.testNSFP WHERE hg18pos=897662 AND alt='C';"''', 'None\tNone\n')

    def testUseRange_1(self):
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #the annotation file of "knownGene" is a range based database
        self.assertSucc('vtools use knownGene --anno_type range --linked_fields chr txStart txEnd')
        self.assertSucc('vtools update variant --set count1=knownGene.exonCount')
        range_out = self.runCmd('vtools execute "select pos, ref, alt, count1 from variant where count1 is not null"', ret='list')
        #using another way to export this table
        self.runCmd('vtools init test -f')
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #this is the default method. the linked_fields have to be in the order in the test below.
        self.assertSucc('vtools use knownGene')
        self.assertSucc('vtools update variant --set count1=knownGene.exonCount')
        def_out = self.runCmd('vtools execute "select pos, ref, alt, count1 from variant where count1 is not null"', ret='list')
        self.assertEqual(range_out, def_out)
        

    def testUseRange_2(self):
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #runCmd('vtools use gwasCatalog')
        #self.assertSucc('vtools show annotation gwasCatalog -v2')
        #the annotation file of "gwasCatalog" is a position-based, but we can use it as range-based as below:
        self.assertSucc('vtools use gwasCatalog --anno_type range --linked_fields chr position-5000 position+5000')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')
        range_out2=len(self.runCmd('vtools execute "select pos, ref, alt, gene_name from variant where gene_name = \'VAMP3\'"', ret='list'))
        self.assertSucc('vtools select variant "gwasCatalog.genes == \'RERE\'" -o variant.chr variant.pos variant.ref variant.alt gwasCatalog.trait gwasCatalog.name gwasCatalog.position gwasCatalog.pValue gwasCatalog.journal gwasCatalog.title gwasCatalog.genes')
        #narrow the range of position
        self.runCmd('vtools init test -f')
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        self.assertSucc('vtools use gwasCatalog --anno_type range --linked_fields chr position-500 position+500')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')
        range_out3=len(self.runCmd('vtools execute "select pos, ref, alt, gene_name from variant where gene_name = \'VAMP3\'"', ret='list'))
        self.assertNotEqual(range_out2, range_out3)
       

    def testUseField_1(self):
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        self.assertFail('vtools use gwasCatalog --anno_type field --linked_fields region')
        #under the option of field in --anno_type, the variable for linked_fields have be in the annotation database that you want to use
        #and --linked_by field in the variant or annotation(already imported) table. 
        #without --linked_by or use the fields that are not in annotation, you will get an error.
        self.assertFail('vtools use gwasCatalog --anno_type field --linked_fields pos')
        self.assertFail('vtools use gwasCatalog --anno_type field --linked_fields pos --linked_by genes')
        self.assertSucc('vtools use gwasCatalog --anno_type field --linked_fields position --linked_by pos')
        #one record was outputed
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')
        #if you use two fields to link the annotation and variant table, the output is different
        #because you added more conditions. The number of linked_fields have to be equal to the number of linked_by
        #nothing was outputed
        self.assertFail('vtools use gwasCatalog --anno_type field --linked_fields chr position --linked_by chr')
        self.assertSucc('vtools use gwasCatalog --anno_type field --linked_fields chr position --linked_by chr pos')
    
    def testUseField_2(self):
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #import the first annotation database
        self.runCmd('vtools use cytoBand')
        self.assertFail('vtools use gwasCatalog --anno_type field --linked_fields region --linked_by genes')
        self.assertSucc('vtools use gwasCatalog --anno_type field --linked_fields region --linked_by cytoBand.name')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')
        self.assertSucc('vtools select variant "gwasCatalog.genes == \'VAMP3\'" -o variant.chr variant.pos variant.ref variant.alt gwasCatalog.trait gwasCatalog.name gwasCatalog.position gwasCatalog.pValue gwasCatalog.journal gwasCatalog.title gwasCatalog.genes')
        #For comparison, if we choise "name" in the linked_fields,nothing will be outputed.  
        self.runCmd('vtools init test -f')
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        self.runCmd('vtools use cytoBand')
        self.assertSucc('vtools use gwasCatalog --anno_type field --linked_fields name --linked_by cytoBand.name')
        self.assertSucc('vtools update variant --set gene_name=1')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')
        self.assertOutput('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"', 'output/use_field.txt')
        
        
    # def testUseVariant(self):
    #     self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
    #     #this is the default method. the linked_fields have to be in the order in the test below.
    #     self.assertSucc('vtools use ESP --anno_type variant --linked_fields chr pos ref alt')
    #     #it is same with
    #     self.assertSucc('vtools use ESP')
    #     # --linked_by option will be ignored if you use the option of --anno_type variant
    #     self.assertFail('vtools use ESP --anno_type variant --linked_fields chr')
    #     self.assertFail('vtools use ESP --anno_type variant --linked_fields chr pos')
    #     self.assertFail('vtools use ESP --anno_type variant --linked_fields chr pos ref')
    #     # because all values are NULL
    #     self.assertFail('vtools update variant --set gene_name=evs.Genes')
    #     self.assertSucc('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"')

    def testUsePosition(self):
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #If we use the option of positon, the linked_fields have to be "chr" and "position"
        #the following command is current, but no variant was linked to the annotation database
        self.assertSucc('vtools use gwasCatalog --anno_type position --linked_fields chr position')
        #the output is none from the command below. in order to get the output, create gene_name in variant first
        self.assertSucc('vtools update variant --set gene_name=1')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        pos_out = self.runCmd('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"', ret='list')
        self.assertEqual(pos_out,['9468354\t-\tA\t1'])
        #using another way to export this table
        self.runCmd('vtools init test -f')
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        self.assertSucc('vtools use gwasCatalog')
        #the output is none from the command too
        self.assertSucc('vtools update variant --set gene_name=1')
        self.assertSucc('vtools update variant --set gene_name=gwasCatalog.genes')
        self.assertOutput('vtools execute "select pos, ref, alt, gene_name from variant where gene_name is not null"', 'output/use_position.txt')
    
    
    def testUseAs(self):
        '''Testing the --as option of command vtools use'''
        self.runCmd('vtools import vcf/SAMP4_complex_variants.vcf --build hg19')
        #this is the default method. the linked_fields have to be in the order in the test below.
        self.assertSucc('vtools use ESP --anno_type variant --linked_fields chr pos ref alt')
        #it is same with
        self.assertSucc('vtools use ESP --as e')
        self.assertSucc('vtools use ESP --as e1')
        self.assertSucc('vtools show annotation e')
        self.assertSucc('vtools output variant chr pos e.chr e1.chr')
        self.assertSucc('vtools select variant "e.chr is not Null" "e1.chr is not Null" --output chr pos e.chr e1.chr')

        
if __name__ == '__main__':
    unittest.main()
