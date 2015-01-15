#/usr/bin/env python
#
# $File: test_import.py $
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
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample, \
    outputOfCmd, getGenotypes, getSamplenames, output2list, getGenotypeInfo

class TestImport(ProcessTestCase):
    
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testInvalidVariant(self):
        'Test importing invalid variants (no position)'
        self.assertSucc('vtools import --build hg18 --format fmt/basic_hg18 txt/invalid.tsv')
        self.assertEqual(numOfSample(), 0)
        # the header should be ignored due ot its POS
        self.assertEqual(numOfVariant(), 10)

    def testImportTXT(self):
        'Test command import'
        self.assertFail('vtools import')
        # Cannot guess input file type from filename
        self.assertFail('vtools import txt/input.tsv')
        # no format information, fail
        self.assertFail('vtools import txt/input.tsv --build hg18')
        # help information
        self.assertSucc('vtools import -h')
        # no build information, fail
        self.assertFail('vtools import --format fmt/basic_hg18 txt/input.tsv')
        # no sample name, a sample with NULL name is created
        self.assertSucc('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        self.assertEqual(numOfSample(), 0)
        self.assertEqual(numOfVariant(), 338)
        self.assertFail('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv --variant chr pos ref not_defined_field --force')
        self.assertFail('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv --variant chr pos --force')
        variants = [x.split() for x in output2list('vtools output variant chr pos ref alt -d"\t"')]
        with open('txt/input.tsv') as inputfile:
            input = [x.split() for x in inputfile.readlines()]
        input = [x[:2] + x[3:] for x in input]
        self.assertEqual(variants, input)
        # test downloading fmt file from the website
        self.assertSucc('vtools import --build hg18 --format ANNOVAR txt/ANNOVAR.txt')
        self.assertFail('vtools import --build hg18 --format ../format/non_existing_fmt txt/input.tsv')
    
    def testGenotypes(self):
        'Test the import of genotypes'
        # use an empty var_info option. The program should not import any var_info for now
        # the following 2 commands are commented out due to running time: very slow to remove the existing sample tables
        #self.assertSucc('vtools import --format fmt/genotypes.fmt txt/genotypes.txt --var_info --build hg18')
        #self.assertFail('vtools output variant chr pos snp_id genet_dist -d"\t"')
        # force re-import the variants using the default info fields now: snp_id and genet_dist are imported 
        self.assertSucc('vtools import --format fmt/genotypes.fmt txt/genotypes.txt --build hg18 --force')
        nsamples = numOfSample()
        nvar = numOfVariant()
        self.assertEqual(nsamples, 49)
        self.assertEqual(nvar, 15)
        # get genotypes for 8 samples
        genotypes = getGenotypes()
        # get samplenames
        samplenames = getSamplenames()
        head = ['#chr','rs','distance','pos','ref','alt'] + samplenames
        variants = [x.split() for x in output2list('vtools output variant chr snp_id genet_dist pos ref alt -d"\t"')]
        with open('txt/genotypes.txt') as inputfile:
            input = [x.split() for x in inputfile.readlines()]
        # test --sample_names
        self.assertEqual(input[0], head)
        input_variants = [x[:6] for x in input[1:]]
        input_genotypes = list(map(list, zip(*[x[6:] for x in input[1:]])))[:8]
        input_genotypes = [list(filter(lambda item: item != '0' and item != '.', item)) for item in input_genotypes]
        # test importing variants 
        self.assertEqual(input_variants, variants)
        # test importing genotypes
        self.assertEqual(input_genotypes, genotypes)

    def testDupGenotype(self):
        'Test importing duplicated genotypes'
        self.assertSucc('vtools import vcf/V1.vcf --sample_name V1 --build hg18')
        self.assertSucc('vtools import vcf/dup_geno.vcf --sample_name DUP --build hg18')
        out = outputOfCmd('vtools show genotypes').split('\n')[1:]
        self.assertEqual(out[0].split()[2], out[1].split()[2])
        #
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertEqual(numOfVariant(), 288)
        g1 = getGenotypes()
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/CEU_dup.vcf.gz --build hg18')
        self.assertEqual(numOfVariant(), 288)
        g2 = getGenotypes()
        self.assertEqual(g1, g2)


    def testANNOVAR(self):
        'Testing the annovar input format'
        self.assertSucc('vtools import --build hg18 --format ../format/ANNOVAR txt/ANNOVAR.txt')
        # one of the variant cannot be imported.
        self.assertEqual(numOfSample(), 0)
        self.assertEqual(numOfVariant(), 11)
        self.assertSucc('vtools import --build hg18 --format ../format/ANNOVAR txt/ANNOVAR.txt --force --sample_name kaiw' )
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 11)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'kaiw\n')
        self.assertSucc('vtools import --build hg18 --format ../format/ANNOVAR_exonic_variant_function txt/annovar.txt.exonic_variant_function' )
        self.assertSucc('vtools output variant mut_type')
        # test for importing user specified var_info
        self.assertSucc('vtools import --build hg18 --format ../format/ANNOVAR_exonic_variant_function txt/annovar.txt.exonic_variant_function --var_info function --force' )
        self.assertEqual(len(output2list('vtools select variant "function is not NULL" -o function')), 78)
        # mut_type should not be imported because it is not specified
        self.assertFail('vtools output variant mut_type')
        
    def testCASAVA18_SNP(self):
        'Testing the CASAVA SNP input format'
        self.assertSucc('vtools import --build hg18 --format ../format/CASAVA18_snps txt/CASAVA18_SNP.txt')
        # 20 new, SNVs, 5 invalid
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 21)
        # sample name should have been scanned from the last line starting with "#"
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'max_gt\n')
        # test for re-naming the sample
        self.assertSucc('vtools import --build hg18 --format ../format/CASAVA18_snps txt/CASAVA18_SNP.txt --force --sample_name casavasnp')
        # both samples exist
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 21)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'max_gt\ncasavasnp\n')
        runCmd('vtools init test -f')
        # test for using user specified genotype information. Have to init a test because of efficiency problem using --force
        self.assertSucc('vtools import --build hg18 --format ../format/CASAVA18_snps txt/CASAVA18_SNP.txt --geno max_gt --geno_info Q_max_gt max_gt_poly_site Q_max_gt_poly_site')
        # now we have 1 genotype field and 3 info field, plus the variant ID: 5 fields in the genotype_x table
        self.assertEqual(len(output2list('vtools execute "PRAGMA table_info(genotype_1)"')), 5)
        # only 1 sample here. Set num=1
        self.assertEqual(getGenotypes(num=1)[0], ['1']*10 + ['2'] + ['1']*10)
        
    def testCASAVA18_INDEL(self):
        'Testing the CASAVA INDEL input format'
        self.assertSucc('vtools import --build hg18 --format ../format/CASAVA18_indels txt/CASAVA18_INDEL.txt')
        # (25 new, 7 insertions, 18 deletions)
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 25)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'max_gtype\n')
        self.assertEqual(getGenotypes(num=1)[0], ['1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '1', '2', '1', '1', '1', '1', '2', '1', '1'])
        
    def testPileup_INDEL(self):
        # this file has one genotype but we do not provide a sample name. Named as "None" when no sample name is specified anywhere
        self.assertSucc('vtools import --build hg18 --format ../format/pileup_indel txt/pileup.indel')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 30)
        # empty sample name
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), '\n')
        # test the MapValue() in the fmt file
        self.assertEqual(getGenotypes(num=1)[0], ['2', '1', '2', '1', '1', '1', '1', '1', '1', '1', '1', '1', \
                                           '1', '2', '1', '1', '1', '1', '1', '2', '1', '1', '2', '1', '1', '1', \
                                           '1', '1', '1', '1'])
    
    def testImportEmpty(self):
        'Test import file without variant'
        self.assertSucc('vtools import vcf/EMPTY.vcf --build hg19')

    def testImportVCF(self):
        'Test command vtools import *.vcf'
        # no build information. Fail
        self.assertFail('vtools import vcf/SAMP1.vcf')
        # use the default vcf format
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'SAMP1\n')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 289)
        self.assertSucc('vtools import vcf/SAMP2.vcf')
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 289+121)
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'SAMP1\nSAMP2\n')
        # geno is empty, i.e, no sample is imported
        self.assertSucc('vtools import vcf/CEU.vcf.gz --geno')
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 698)
        # file will be ignored if re-imported
        self.assertFail('vtools import vcf/SAMP1.vcf')
        self.assertEqual(numOfSample(), 2)
        # force re-import the same file with samples
        self.assertSucc('vtools import vcf/CEU.vcf.gz -f')
        self.assertEqual(numOfSample(), 62)
        self.assertEqual(numOfVariant(), 698)
        # import additional information on variants and on genotypes.
        # DP and DP_geno are fields provided in the default vcf.fmt
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --var_info DP --geno_info DP_geno --build hg18')
        self.assertEqual(numOfSample(), 60)
        self.assertEqual(numOfVariant(), 288)
        self.assertSucc('vtools output variant DP')
        self.assertEqual(len(output2list('vtools execute "PRAGMA table_info(genotype_1)"')), 3)
        self.assertEqual(output2list('vtools execute "select DP_geno from genotype_1"'), ["6","1","1","0","1","0","0","0","3","1","0","0","0","0","0","2","0","0","1","3","0","0","3","0","0","2","2","7","2","4","2","3","2","0","1","6","0","2","1","3","1","5","3","1","5","1","0","1","2","1","4","2","3","7","1","0","1","1","2","0","0","5","1","0","0","1","3","1","0","3","3","1","0","4","0","4","3","2","0","2","2","4","4","0","0","6","2","0","0","1","1","0","0","2","5","1","1","1","2","1","0","1","1","4","0","6","1","3","3","1","3","1","5","0","1","2","1","0","0","1","1","2","1","3","5","1","3","3","1","0","2","0","1","3","5","9","2","4","2","2","1","2","1","2","1","2","0","7","7","9","6","1","1","1","1","1","7","9","3","2","1","1","2","1","1","1","2","2","5","4","1","5","5","3","2","2","0","3","3","0","2","2","2","3","5","1","2","3","1","3","0","8","2","3","6","2","2","0","4","2","7","1","3","0","3","4","7","3","1","3","4","2","1","3","2","1","1","1","7","1","2","2","0","0","1","2","3","1","4","2","1","1","2","4","1","2","4","3","1","5","2","8","8","5","5","3","6","7","8","5","3","2","5","7","0","3","3","3","2","2","5","1","12","1","1","2","2","0","6","1","2","5","3","3","3","1","1","1","0","0","1","2","2","0","1","3","0"])

    def testImportVCFIndel(self):
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18')
        self.assertEqual(numOfSample(), 0)
        self.assertEqual(numOfVariant(), 137)
        self.assertSucc('vtools import vcf/SAMP4_complex_variants.vcf --geno_info')
        self.assertEqual(numOfSample(), 0)
        self.assertEqual(numOfVariant(), 137+12159)
        
    def testMPImport(self):
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j1')
        samples = outputOfCmd('vtools show samples -l -1')
        genotype = outputOfCmd('vtools show genotypes -l -1')
        variants = outputOfCmd('vtools show table variant -l -1')
        genotypes = []
        for i in range(60):
            genotypes.append(outputOfCmd("vtools execute 'select * from genotype.genotype_{}'".format(i+1)))
        #
        # compare results with -j3
        #
        runCmd('vtools init test -f')
        # if more than one reader is used, the order of mutants in some cases will be changed, leading
        # to different variant ids.
        runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j3')
        self.assertEqual(samples, outputOfCmd('vtools show samples -l -1'))
        self.assertEqual(genotype, outputOfCmd('vtools show genotypes -l -1'))
        self.assertEqual(variants, outputOfCmd('vtools show table variant -l -1'))
        for i in range(60):
            self.assertEqual(genotypes[i], outputOfCmd("vtools execute 'select * from genotype.genotype_{}'".format(i+1)))
        #
        # compare results with -j10
        #
        runCmd('vtools init test -f')
        runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j10')
        self.assertEqual(samples, outputOfCmd('vtools show samples -l -1'))
        self.assertEqual(genotype, outputOfCmd('vtools show genotypes -l -1'))
        self.assertEqual(variants, outputOfCmd('vtools show table variant -l -1'))
        for i in range(60):
            self.assertEqual(genotypes[i], outputOfCmd("vtools execute 'select * from genotype.genotype_{}'".format(i+1)))
    
    def testMPImportMultiFiles(self):
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/V1.vcf vcf/V2.vcf vcf/V3.vcf --build hg18 -j1')
        samples = outputOfCmd('vtools show samples -l -1')
        genotype = outputOfCmd('vtools show genotypes -l -1')
        variants = outputOfCmd('vtools show table variant -l -1')
        genotypes = []
        for i in range(3):
            genotypes.append(outputOfCmd("vtools execute 'select * from genotype.genotype_{}'".format(i+1)))
        #
        # compare results with -j3
        #
        runCmd('vtools init test -f')
        runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/V1.vcf vcf/V2.vcf vcf/V3.vcf --build hg18 -j4')
        self.assertEqual(samples, outputOfCmd('vtools show samples -l -1'))
        self.assertEqual(genotype, outputOfCmd('vtools show genotypes -l -1'))
        self.assertEqual(variants, outputOfCmd('vtools show table variant -l -1'))
        for i in range(3):
            self.assertEqual(genotypes[i], outputOfCmd("vtools execute 'select * from genotype.genotype_{}'".format(i+1)))
 
    def testMixedBuild(self):
        'Test importing vcf files with different reference genomes'
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 289)
        # 104 records in SAMP1.vcf failed to map to hg19
        self.assertSucc('vtools import vcf/var_format.vcf --geno safe_GT --build hg19')
        # all records in var_format.vcf are mapped to hg18
        self.assertEqual(numOfSample(), 1+1)
        self.assertEqual(numOfVariant(), 289 + 98)
        # 19 out of 121 records failed to map.
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        self.assertEqual(numOfSample(), 3)
        self.assertEqual(numOfVariant(), 289 + 98 + 121)
        #
        # this is a difficult test to pass, basically, we will create another project
        # in reverse order of reference genomes, using reversed liftover, and see
        # it the output is the same
        out1 = outputOfCmd('vtools output variant bin chr pos alt_bin alt_chr alt_pos')
        self.assertSucc('vtools init test -f')
        self.assertSucc('vtools import vcf/var_format.vcf --geno safe_GT --build hg19')
        self.assertEqual(numOfVariant(), 98)
        self.assertEqual(numOfSample(), 1)
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertEqual(numOfSample(), 2)
        # 101 cannot be mapped. damn.
        self.assertEqual(numOfVariant(), 98 + 289)
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        # 19 out of 121 records failed to map.
        self.assertEqual(numOfSample(), 3)
        self.assertEqual(numOfVariant(), 98 + 289 + 121)
        out2 = outputOfCmd('vtools output variant alt_bin alt_chr alt_pos bin chr pos')
        #
        out1 = '\n'.join([x for x in sorted(out1.split('\n')) if 'NA' not in x]).replace(' ', '')
        out2 = '\n'.join([x for x in sorted(out2.split('\n')) if 'NA' not in x]).replace(' ', '')
        self.assertEqual(out1, out2)
    
    def testImportMyVCF(self):
        'Test a customized vcf import'
        self.assertSucc('vtools import --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        # test importing self defined var_info
        self.assertEqual(output2list('vtools output variant \
                                     DP MQ NS AN AC AF AB LBS_A1 LBS_A2 \
                                     LBS_C1 LBS_C2 LBS_G1 LBS_G2 LBS_T1 LBS_T2 OBS_A1 \
                                     OBS_A2 OBS_C1 OBS_C2 OBS_G1 OBS_G2 OBS_T1 OBS_T2 \
                                     STR STZ CBR CBZ QBR QBZ MBR MSR MBZ \
                                     IOR IOZ IOH IOD AOI AOZ ABE ABZ BCS \
                                     FIC LQR MQ0 MQ10 MQ20 MQ30 ANNO SVM \
                                     -l 1 -d"\t"')[0].split('\t'),
                         ['472', '28', '308', '616', '1', '0.002774', '0.2738', '64', '27', '63', '37', '6558', '1508', \
                          '93', '102', '166', '63', '60', '37', '509234', '225198', '163', '92', '-0.001', '-1.004', '0.002', \
                          '1.354', '-0.032', '-27.065', '0.007', '-0.029', '6.137', '0.52', '-17.939', '0.0', '0.0', '-21.32', '-3.382', \
                          '0.901', '29.686', '-2003.904', '0.948', '0.011', '0.997', '0.997', '1.0', '1.0', \
                          'nonsynonymous:OR4F5:NM_001005484:exon1:c.G26A:p.G9D,', '-1.4352462'])
        # test importing self-defined genotypes with VcfGenotype(default=('0',))
        # code missing genotypes as None and wild-type as '0'
        self.assertEqual(getGenotypes(num=4), [['-1', '-1'], ['0', '2'], ['0', '-1', '-1'], ['0', '-1', '-1']])
        # test importing self-defined geno_info.
        # PL3* are passed into database as a 2X1 "transposed tuple" -- works here.
        # See 'missing_gen.fmt'
        genotypeInfo = getGenotypeInfo(num=4, info=['GT', 'GQ', 'GD', 'PL_1', 'PL_2', 'PL_3', 'PL3_1', 'PL3_2', 'PL3_3'])
        genotypeInfo = ['\t'.join(x) for x in genotypeInfo]
        genotypeInfo = [x.split('\t') for x in genotypeInfo]
        self.assertEqual(genotypeInfo, [['-1', '3', '1', 'None', 'None', 'None', '0', \
                                         '3', '4', '-1', '3', '1', 'None', 'None', 'None', \
                                         '3', '4', '4'], ['0', '93', '27', '0', '81', \
                                        '218', 'None', 'None', 'None', '2', '6', '3', 'None', 'None', \
                                        'None', '43', '6', '0'], ['0', '15', '1', '0', '3', '20', 'None', \
                                        'None', 'None', '-1', '4', '1', 'None', 'None', 'None', '24', '24', '24', \
                                        '-1', '4', '1', 'None', 'None', 'None', '3', '3', '0'], ['0', '100', \
                                        '32', '0', '96', '255', 'None', 'None', 'None', '-1', '3', '2', 'None', \
                                        'None', 'None', '0', '6', '54', '-1', '3', '2', 'None', 'None', 'None', '6', \
                                        '54', '54']])

    def testInsertDelete(self):
        'Testing the number of insertions and deletions'
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18')
        self.assertEqual(numOfVariant(), 137)
        self.assertSucc('vtools select variant \'ref="-"\' --output chr pos ref alt')
        self.assertOutput('vtools select variant \'ref="-"\' --output chr pos ref alt', '',  0, 'output/import_vcf_ref.txt')
        self.assertOutput('vtools select variant \'ref="-"\' --output chr pos ref alt', '''1	10434	-	C\n1	54790	-	T\n1	81963	-	AA\n1	82135	-	AAAAAAAAAAAAAA\n1	83787	-	A''',  5, 'output/import_vcf_ref.txt')
        self.assertOutput('vtools select variant \'ref="-"\' --output chr pos ref alt', '''1	10434	-	C\n1	54790	-	T\n1	81963	-	AA\n1	82135	-	AAAAAAAAAAAAAA\n1	83787	-	A\n''',  -5, 'output/import_vcf_ref.txt')
        self.assertEqual(outputOfCmd('vtools select variant \'ref="-"\' --count'), '74\n') 
        self.assertOutput('vtools select variant \'alt="-"\' --output chr pos ref alt', '',  0, 'output/import_vcf_alt.txt')
        self.assertEqual(outputOfCmd('vtools select variant \'alt="-"\' --count'), '55\n') 

    def testSampleName_single(self):
        'Testing the import of sample names'
        #Testing one sample per file with the default setting in vtools import
        #sample name is given in file
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        self.assertEqual(numOfSample(), 2)
        self.assertEqual(numOfVariant(), 289 + 121)
        self.assertOutput('vtools show genotypes', '',  0, 'output/vcf_single_sampleName_genotype.txt')

    def testNo_SampleName(self):
        #Testing one sample per file with the default setting(NO sample name)
        #sample name was not given in file, then there is no information about sample name and genotypes except you assign one for it.
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18')
        self.assertEqual(numOfSample(), 0)
        self.assertEqual(numOfVariant(), 137)
    
    def testNo_SampleName_assign(self):
        #Assign a sample name if the sample name is not in file
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18 --sample_name vcf_test3')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(), 137) 
        self.assertEqual(outputOfCmd('vtools execute "select sample_name from sample"'), 'vcf_test3\n')
        
    def testSampleName_single_assign(self):
        #Testing one sample per file with the --sample_name option
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18 --sample_name samp_vcf1')
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18 --sample_name samp_vcf2')
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18 --sample_name samp_vcf3')
        self.assertEqual(numOfSample(), 3)
        self.assertEqual(numOfVariant(), 289 + 121 + 137)
        self.assertOutput('vtools show genotypes', '',  0, 'output/vcf_assigned_sample_name_genotype.txt')
        
    def testSampleName_multiple(self):
        #Testing multiple samples in ONE vcf file with default setting
        self.assertSucc('vtools import vcf/500SAMP.vcf --build hg18')
        self.assertEqual(numOfSample(), 501)
        self.assertEqual(numOfVariant(),5)
        self.maxDiff = None
        self.assertOutput('vtools show genotypes', '',  0, 'output/vcf_multiple_samples_genotypes.txt')
       
    def testSampleName_multiple_assign(self):
        #Testing multiple samples in ONE vcf file with --sample_name option
        #Only one sample was generated, no genotype information were imported
        self.assertFail('vtools import vcf/500SAMP.vcf --build hg18 --sample_name output/vcf_multiple_sample_name.txt')
        self.assertFail('vtools import vcf/500SAMP.vcf --build hg18 --sample_name n1 n2 n3 n4 n5')
        self.assertFail('vtools import vcf/SAMP1.vcf --build hg18 --sample_name samp_vcf1 samp_vcf2 samp_vcf3')
        #you can assign sample names as blow
        with open('output/vcf_multiple_sample_name.txt') as inputfile:
              input = [x.strip() for x in inputfile.readlines()]
        self.assertSucc("vtools import vcf/500SAMP.vcf --build hg19 --sample_name {}".format(' '.join(input[1:])))

    def testCsvImport1(self):
        self.assertSucc('vtools import txt/test.csv --format ../format/csv.fmt --build hg18 --sample_name samp_csv')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(),687)
        self.maxDiff = None
        variants = [x.split() for x in output2list('vtools output variant chr pos ref alt')]
        with open('txt/test.csv') as inputfile:
            input = [x.split(',') for x in inputfile.readlines()]
        input = [x[:4] for x in input]
        self.assertEqual(variants, input)

    def testCGAImport(self):
        self.assertSucc('vtools import txt/CGA.tsv.bz2 --format ../format/CGA.fmt --build hg18 --sample_name samp_csv')
        self.assertEqual(numOfSample(), 1)
        self.assertEqual(numOfVariant(),95)
        self.maxDiff = None
        variants = [x.split() for x in output2list('vtools output variant chr pos ref alt')]
        with open('output/CGA_variant.txt') as inputfile:
            input = [x.split(',') for x in inputfile.readlines()]
        input = [x[:4] for x in input]
        self.assertEqual(variants, input) 
        self.assertOutput('vtools show genotypes', '''sample_name\tfilename       \tnum_genotypes\tsample_genotype_fields
samp_csv   \ttxt/CGA.tsv.bz2\t95           \tGT,allele1VarScoreVAF,allele2VarScoreVAF,allele1VarScoreEAF,allele2VarScoreEAF
''')

    def testMultiSamples_1(self):
        #the files are coming from one custmer
        runCmd('vtools init test -f')
        #based on the files format, Dr. Peng modified the format file: new_format.fmt
        #the origianl file contains spaces in lines and the default deliminator is "tab"
        #and then the numbers of columns per line are not equal, and then could not be imported 
        self.assertSucc('vtools import --format fmt/multi_index.fmt txt/sample_chr22.txt  --build hg18')
        self.assertEqual(numOfSample(), 3)
        self.maxDiff=None
        self.assertOutput('vtools show table variant', '''Name:                   variant\nDescription:            Master variant table\nCommand:\nFields:                 variant_id, bin, chr, pos, ref, alt\nNumber of variants:     6\n''', skip=3)
        self.assertOutput('vtools show samples', '''sample_name\tfilename
SMP1       \ttxt/sample_chr22.txt
SMP2       \ttxt/sample_chr22.txt
SMP3       \ttxt/sample_chr22.txt
''')

    def testMultiSamples_2(self):
        #the files are coming from one custmer
        runCmd('vtools init test -f')
        #based on the files format, Dr. Peng modified the format file: new_format.fmt
        self.assertSucc('vtools import --format fmt/multi_index.fmt txt/sample_1_chr22.txt  --build hg18')
        self.assertEqual(numOfSample(), 3)
        self.assertOutput('vtools show table variant', '''Name:                   variant\nDescription:            Master variant table\nCommand:\nFields:                 variant_id, bin, chr, pos, ref, alt\nNumber of variants:     9\n''', skip=3)
        self.assertOutput('vtools show samples', 'sample_name\tfilename\nSMP1       \ttxt/sample_1_chr22.txt\nSMP2       \ttxt/sample_1_chr22.txt\nSMP3       \ttxt/sample_1_chr22.txt\n')
     

if __name__ == '__main__':
    unittest.main()
