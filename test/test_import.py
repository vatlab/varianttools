#!/usr/bin/env python
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
from testUtils import ProcessTestCase 
from variant_tools.accessor import *

class TestImport(ProcessTestCase):
    def testInvalidVariant(self):
        'Test importing invalid variants (<DEL>)'
        self.assertSucc('vtools import --build hg18 --format fmt/basic_hg18 txt/invalid.tsv')
        self.assertProj(numOfSamples= 0, numOfVariants=10)

    def testImportCommand(self):
        'Test command import'
        self.assertFail('vtools import')
        # Cannot guess input file type from filename
        self.assertFail('vtools import txt/input.tsv')
        # no format information, fail
        self.assertFail('vtools import txt/input.tsv --build hg18')

    def testImportTXT(self):
        # help information
        self.assertSucc('vtools import -h')
        # no build information, fail
        self.assertFail('vtools import --format fmt/basic_hg18 txt/input.tsv')
        # no sample name, a sample with NULL name is created
        self.assertSucc('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        self.assertProj(numOfSamples= 0, numOfVariants=1446)
        self.assertFail('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv --variant chr pos ref not_defined_field --force')
        self.assertFail('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv --variant chr pos --force')
        self.assertOutput('vtools output variant chr pos ref alt -d"\t"', 'output/import_txt_1.txt')
        # test downloading fmt file from the website
        self.assertFail('vtools import --build hg18 --format non_existing_fmt txt/input.tsv')
    
    def testGenotypes(self):
        'Test the import of genotypes'
        self.assertSucc('vtools import --format fmt/genotypes.fmt txt/genotypes.txt --build hg18')
        with open('txt/genotypes.txt') as inputfile:
            firstline = inputfile.readline()
        self.assertProj(numOfSamples=49, numOfVariants=15, sampleNames=firstline.strip().split()[6:])
        self.assertOutput('vtools output variant chr snp_id genet_dist pos ref alt -d"\t"', 'output/import_genotype_1.txt')

    def testDupGenotype(self):
        'Test importing duplicated genotypes'
        self.assertSucc('vtools import vcf/V1.vcf --sample_name V1 --build hg18')
        self.assertSucc('vtools import vcf/dup_geno.vcf --sample_name DUP --build hg18')
        self.assertOutput('vtools show genotypes', 'output/import_genotype_2.txt')

    def testDupGenotype1(self):
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertProj(numOfVariants=288)
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import vcf/CEU_dup.vcf.gz --build hg18')
        self.assertProj(numOfVariants=288)


    # @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite3","HDF5 version is not implemented for this test")
    def testANNOVAR(self):
        'Testing the annovar input format'
        self.assertSucc('vtools import --build hg18 --format ../resources/format/ANNOVAR txt/ANNOVAR.txt')
        # one of the variant cannot be imported.
        self.assertProj(numOfSamples= 0, numOfVariants=11)
        self.assertSucc('vtools import --build hg18 --format ../resources/format/ANNOVAR txt/ANNOVAR.txt --force --sample_name kaiw' )
        self.assertProj(numOfSamples= 1, numOfVariants=11)
        self.assertOutput('vtools execute "select sample_name from sample"', 'kaiw\n')
        self.assertSucc('vtools import --build hg18 --format ../resources/format/ANNOVAR_exonic_variant_function txt/annovar.txt.exonic_variant_function' )
        self.assertSucc('vtools output variant mut_type')
        # test for importing user specified var_info
        self.assertSucc('vtools import --build hg18 --format ../resources/format/ANNOVAR_exonic_variant_function txt/annovar.txt.exonic_variant_function --var_info function --force' )
        self.assertSucc('vtools select variant "function is not NULL" -t function')
        self.assertProj(numOfVariants={'function': 78})
        
    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite3","HDF5 version is not implemented for this test")
    def testCASAVA18_SNP(self):
        'Testing the CASAVA SNP input format'
        self.assertSucc('vtools import --build hg18 --format ../resources/format/CASAVA18_snps txt/CASAVA18_SNP.txt')
        # 20 new, SNVs, 5 invalid
        self.assertProj(numOfSamples= 1, numOfVariants=21)
        # sample name should have been scanned from the last line starting with "#"
        self.assertProj(sampleNames=['max_gt'])
        # test for re-naming the sample
        self.assertSucc('vtools import --build hg18 --format ../resources/format/CASAVA18_snps txt/CASAVA18_SNP.txt --force --sample_name casavasnp')
        # both samples exist
        self.assertProj(numOfSamples= 2, numOfVariants=21, sampleNames=['max_gt', 'casavasnp'])
        self.runCmd('vtools init test -f --store '+self.storeMode)
        # test for using user specified genotype information. Have to init a test because of efficiency problem using --force
        self.assertSucc('vtools import --build hg18 --format ../resources/format/CASAVA18_snps txt/CASAVA18_SNP.txt --geno max_gt --geno_info Q_max_gt max_gt_poly_site Q_max_gt_poly_site')
        # now we have 1 genotype field and 3 info field, plus the variant ID: 5 fields in the genotype_x table
        if self.storeMode=="sqlite":
            self.assertProj(numOfColumns={'genotype_1': 5})
        # only 1 sample here. Set num=1
        self.assertProj(genotype={1: ['1']*10 + ['2'] + ['1']*10})
        
    def testCASAVA18_INDEL(self):
        'Testing the CASAVA INDEL input format'
        self.assertSucc('vtools import --build hg18 --format ../resources/format/CASAVA18_indels txt/CASAVA18_INDEL.txt')
        # (25 new, 7 insertions, 18 deletions)
        self.assertProj(numOfSamples= 1, numOfVariants=25, sampleNames=['max_gtype'],
            genotype={1: '1111111111111111121111211'})
        
    def testPileup_INDEL(self):
        # this file has one genotype but we do not provide a sample name. Named as "None" when no sample name is specified anywhere
        self.assertSucc('vtools import --build hg18 --format ../resources/format/pileup_indel txt/pileup.indel')
        self.assertProj(numOfSamples= 1, numOfVariants=30, sampleNames=[''],
            genotype={1: '212111111111121111121121111111'})
    
    def testImportEmpty(self):
        'Test import file without variant'
        self.assertSucc('vtools import vcf/EMPTY.vcf --build hg19')

    def testImportVCF(self):
        'Test command vtools import *.vcf'
        # no build information. Fail
        self.assertFail('vtools import vcf/SAMP1.vcf')
        # use the default vcf format
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertProj(numOfSamples= 1, numOfVariants=289, sampleNames=['SAMP1'])
        self.assertSucc('vtools import vcf/SAMP2.vcf')
        self.assertProj(numOfSamples= 2, numOfVariants=289+121, sampleNames=['SAMP1', 'SAMP2'])
        # geno is empty, i.e, no sample is imported
        self.assertSucc('vtools import vcf/CEU.vcf.gz --geno')
        self.assertProj(numOfSamples= 2, numOfVariants=698)
        # file will be ignored if re-imported
        self.assertFail('vtools import vcf/SAMP1.vcf')
        self.assertProj(numOfSamples= 2)
        # force re-import the same file with samples
        self.assertSucc('vtools import vcf/CEU.vcf.gz -f')
        self.assertProj(numOfSamples= 62, numOfVariants=698)
        # import additional information on variants and on genotypes.
        # DP and DP_geno are fields provided in the default vcf.fmt
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import vcf/CEU.vcf.gz --var_info DP --geno_info DP --build hg18')
        self.assertSucc('vtools output variant DP')
        if self.storeMode=="sqlite":
            self.assertProj(numOfColumns={'genotype_1': 3})
        self.assertProj(numOfSamples= 60, numOfVariants=288)

    def testImportVCFIndel(self):
        'Test importing Indel from VCF files'
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
        self.assertProj(numOfSamples= 0, numOfVariants=134)
        self.assertSucc('vtools import vcf/SAMP4_complex_variants.vcf --geno_info')
        self.assertProj(numOfSamples= 0, numOfVariants=11877)
    
    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","HDF5 version is not implemented for this test")   
    def testMPImport(self):
        'Test multi-processing import'
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j1')
        self.assertOutput('vtools show samples -l -1', 'output/import_mpi_samples.txt')
        self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_genotypes.txt')
        self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_variant.txt', partial=-3)
        self.assertOutput(["vtools execute 'select * from genotype.genotype_{}'".format(i+1) for i in range(20)],
                'output/import_mpi_genotype.txt')
        #
        # compare results with -j3
        #
        self.runCmd('vtools init test -f --store '+self.storeMode)
        # if more than one reader is used, the order of mutants in some cases will be changed, leading
        # to different variant ids.
        self.runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j3')
        self.assertOutput('vtools show samples -l -1', 'output/import_mpi_samples.txt')
        self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_genotypes.txt')
        self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_variant.txt', partial=-3)
        self.assertOutput(["vtools execute 'select * from genotype.genotype_{}'".format(i+1) for i in range(20)],
                'output/import_mpi_genotype.txt')
        #
        # compare results with -j10
        #
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18 -j10')
        self.assertOutput('vtools show samples -l -1', 'output/import_mpi_samples.txt')
        self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_genotypes.txt')
        self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_variant.txt', partial=-3)
        self.assertOutput(["vtools execute 'select * from genotype.genotype_{}'".format(i+1) for i in range(20)],
                'output/import_mpi_genotype.txt')
    
    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","HDF5 version is not implemented for this test")
    def testMPImportMultiFiles_sqlite(self):
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import vcf/V1.vcf vcf/V2.vcf vcf/V3.vcf --build hg18 -j1')
        self.assertOutput('vtools show samples -l -1', 'output/import_mpi_multi_samples.txt')
        self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_multi_genotypes.txt')
        self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_multi_variant.txt', partial=-3)
        self.assertOutput(["vtools execute 'select * from genotype.genotype_{}'".format(i+1) for i in range(3)],
                'output/import_mpi_multi_genotype.txt')
        #
        # compare results with -j3
        #
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools admin --set_runtime_option "import_num_of_readers=0"')
        self.assertSucc('vtools import vcf/V1.vcf vcf/V2.vcf vcf/V3.vcf --build hg18 -j4')
        self.assertOutput('vtools show samples -l -1', 'output/import_mpi_multi_samples.txt')
        self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_multi_genotypes.txt')
        self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_multi_variant.txt', partial=-3)
        self.assertOutput(["vtools execute 'select * from genotype.genotype_{}'".format(i+1) for i in range(3)],
                'output/import_mpi_multi_genotype.txt')

    @unittest.skipUnless(os.getenv("STOREMODE")=="hdf5","sqlite version is not implemented for this test")
    def testMPImportMultiFiles_hdf5(self):
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import vcf/V1.vcf vcf/V2.vcf vcf/V3.vcf --build hg18 -j1')
        # self.assertOutput('vtools show samples -l -1', 'output/import_mpi_multi_samples.txt')
        # self.assertOutput('vtools show genotypes -l -1', 'output/import_mpi_multi_genotypes.txt')
        # self.assertOutput('vtools show table variant -l -1', 'output/import_mpi_multi_variant.txt', partial=-3)
        with open(os.devnull, 'w') as fnull:
            fileResult = subprocess.check_output('vtools execute "SELECT HDF5 FROM sample WHERE sample_id =1"', shell=True,
                            stderr=fnull).decode()
            HDF5FileName=fileResult.rstrip()
            accessEngine=Engine_Access.choose_access_engine(HDF5FileName)
            # geno=accessEngine.get_geno_by_sample_ID(1,"GT_geno")
            geno=[]
            for rownames,colnames,genoinfo in accessEngine.get_all_genotype([1]):
                for idx,rowname in enumerate(rownames):
                    genotype=genoinfo[idx]
                    if np.isnan(genotype):
                        genotype=-1
                    geno.append([rowname,genotype])
            geno=np.array(geno)
            proj_output="".join([str(int(val[0]))+"\t"+str(int(val[1]))+"\n" for val in geno])
            output='output/import_mpi_multi_genotype_hdf5.txt'
            if os.path.isfile(output):
                with open(output, 'r') as cf:
                    output = cf.read()
            self.compare(proj_output, output)


    def testMixedBuild(self):
        'Test importing vcf files with different reference genomes'
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertProj(numOfSamples= 1, numOfVariants=289)
        # 104 records in SAMP1.vcf failed to map to hg19
        self.assertSucc('vtools import vcf/var_format.vcf --geno safe_GT --build hg19')
        # all records in var_format.vcf are mapped to hg18
        self.assertProj(numOfSamples= 1+1, numOfVariants=289 + 98)
        # 19 out of 121 records failed to map.
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        self.assertProj(numOfSamples= 3, numOfVariants=289 + 98 + 121)
        #
        # this is a difficult test to pass, basically, we will create another project
        # in reverse order of reference genomes, using reversed liftover, and see
        # it the output is the same
        self.assertOutput('vtools output variant bin chr pos alt_bin alt_chr alt_pos -d "\t"', 'output/import_mixed_build.txt')
        self.assertSucc('vtools init test -f')
        self.assertSucc('vtools import vcf/var_format.vcf --geno safe_GT --build hg19')
        self.assertProj(numOfVariants=98, numOfSamples= 1)
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        # 101 cannot be mapped. damn.
        self.assertProj(numOfSamples= 2, numOfVariants=98 + 289)
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        # 19 out of 121 records failed to map.
        self.assertProj(numOfSamples= 3, numOfVariants=98 + 289 + 121)
        self.assertOutput('vtools output variant alt_bin alt_chr alt_pos bin chr pos -d "\t"', 'output/import_mixed_build.txt', 
            # compare common variants (without NA) and in sorted order
            lambda x: sorted([i for i in x if 'NA' not in i]))
    
    def testImportMyVCF(self):
        'Test a customized vcf import'
        self.assertSucc('vtools import --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        # test importing self defined var_info
        self.assertOutput('vtools output variant DP MQ NS AN AC AF AB LBS_A1 LBS_A2 LBS_C1 LBS_C2 LBS_G1 LBS_G2 LBS_T1 LBS_T2 OBS_A1 \
            OBS_A2 OBS_C1 OBS_C2 OBS_G1 OBS_G2 OBS_T1 OBS_T2 STR STZ CBR CBZ QBR QBZ MBR MSR MBZ IOR IOZ IOH IOD AOI AOZ ABE ABZ BCS \
            FIC LQR MQ0 MQ10 MQ20 MQ30 ANNO SVM', 'output/import_customized.txt')
        # test importing self-defined genotypes with VcfGenotype(default=('0',))
        # code missing genotypes as None and wild-type as '0'
        if self.storeMode=="sqlite":
            self.assertProj(genotype={1: ['-1', '-1'], 2: ['0', '2'], 3: ['0', '-1', '-1'], 4: ['0', '-1', '-1']},
                genoInfo={(1, 'GT'): ['-1', '-1'], (1, 'GQ'): ['3', '3'], (1, 'GD'): ['1', '1'], 
                    (1, 'PL_1'): ['None', 'None'], (1, 'PL_2'): ['None', 'None'], (1, 'PL3_1'): ['0', '3']})

    def testInsertDelete(self):
        'Testing the number of insertions and deletions'
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
        self.assertOutput('''vtools select variant 'ref="-"' --output chr pos ref alt''', 'output/import_vcf_ref.txt')
        self.assertProj(numOfVariants=134)
        self.assertOutput('''vtools select variant 'ref="-"' --count''', '73\n') 
        self.assertOutput('''vtools select variant 'alt="-"' --output chr pos ref alt''', 'output/import_vcf_alt.txt')
        self.assertOutput('''vtools select variant 'alt="-"' --count''', '53\n') 


    def testSampleName_single(self):
        'Testing the import of sample names'
        #Testing one sample per file with the default setting in vtools import
        #sample name is given in file
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18')
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        self.assertProj(numOfSamples= 2, numOfVariants=289 + 121)
        self.assertOutput('vtools show genotypes', 'output/vcf_single_sampleName_genotype.txt')

    def testNo_SampleName(self):
        #Testing one sample per file with the default setting(NO sample name)
        #sample name was not given in file, then there is no information about sample name and genotypes except you assign one for it.
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
        self.assertProj(numOfSamples= 0, numOfVariants=134)
    
    def testNo_SampleName_assign(self):
        #Assign a sample name if the sample name is not in file
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg19 --sample_name vcf_test3')
        self.assertProj(numOfSamples= 1, numOfVariants=134, sampleNames=['vcf_test3'])
    
    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","hdf5 version is not implemented for this test")    
    def testSampleName_single_assign(self):
        #Testing one sample per file with the --sample_name option
        self.assertSucc('vtools import vcf/SAMP1.vcf --build hg18 --sample_name samp_vcf1')
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18 --sample_name samp_vcf2')
        self.assertSucc('vtools import vcf/SAMP3_complex_variants.vcf --build hg18 --sample_name samp_vcf3')
        self.assertProj(numOfSamples= 3, numOfVariants=545)
        self.assertOutput('vtools show genotypes', 'output/vcf_assigned_sample_name_genotype.txt')
        
    def testSampleName_multiple(self):
        #Testing multiple samples in ONE vcf file with default setting
        self.assertSucc('vtools import vcf/500SAMP.vcf --build hg18')
        self.assertProj(numOfSamples= 501, numOfVariants=5)
        self.assertOutput('vtools show genotypes', 'output/vcf_multiple_samples_genotypes.txt')
       
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
        self.assertSucc('vtools import txt/test.csv --format ../resources/format/csv.fmt --build hg19 --sample_name samp_csv')
        self.assertProj(numOfSamples=1, numOfVariants=687)
        self.assertOutput('vtools output variant chr pos ref alt', 'output/import_csv.txt')

    def testCGAImport(self):
        self.assertSucc('vtools import txt/CGA.tsv.bz2 --format ../resources/format/CGA.fmt --build hg19 --sample_name samp_csv')
        self.assertProj(numOfSamples=1, numOfVariants=95)
        self.assertOutput('vtools output variant chr pos ref alt', 'output/import_cga.txt') 
        if self.storeMode=="sqlite":
            self.assertOutput('vtools show genotypes', 'output/import_cga_phenotype.txt')

    def testMultiSamples_1(self):
        #the files are coming from one custmer
        self.assertSucc('vtools import --format fmt/multi_index.fmt txt/sample_chr22.txt  --build hg19')
        self.assertProj(numOfSamples=3)
        self.assertOutput('vtools show table variant', 'output/import_multi_sample_variant.txt', -4)
        if self.storeMode=="sqlite":
            self.assertOutput('vtools show samples', 'output/import_multi_sample_samples.txt')
        elif self.storeMode=="hdf5":
            self.assertOutput('vtools show samples', 'output/import_multi_sample_samples_hdf5.txt')

    def testMultiSamples_2(self):
        #the files are coming from one custmer
        self.assertSucc('vtools import --format fmt/multi_index.fmt txt/sample_1_chr22.txt  --build hg19')
        self.assertProj(numOfSamples=3)
        self.assertOutput('vtools show table variant', 'output/import_multi_sample2_variant.txt', -4)
        if self.storeMode=="sqlite":
            self.assertOutput('vtools show samples', 'output/import_multi_sample2_samples.txt')
        elif self.storeMode=="hdf5":
            self.assertOutput('vtools show samples', 'output/import_multi_sample2_samples_hdf5.txt')
     

if __name__ == '__main__':
    unittest.main()
