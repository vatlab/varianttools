#!/usr/bin/env python
#
# $File: test_export.py $
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


class TestExport(ProcessTestCase):
    
    def testExportVcfSnv_single(self):
        'Test command export in vcf format'
        #the output format of vcf has several options:
        #                           [--format_string [FORMAT_STRING [FORMAT_STRING ...]]]
        #                           [--geno_info [GENO_INFO [GENO_INFO ...]]]
        #                           [--geno [GENO [GENO ...]]]
        #                           [--var_info [VAR_INFO [VAR_INFO ...]]]
        #                           [--phase_sep [PHASE_SEP [PHASE_SEP ...]]]
        #                           [--wildtype_code [WILDTYPE_CODE [WILDTYPE_CODE ...]]]
        #                           [--id [ID [ID ...]]]
        #                           [--pos [POS [POS ...]]]
        #                           [--ref [REF [REF ...]]]
        #                           [--alt [ALT [ALT ...]]]
        #                           [--qual [QUAL [QUAL ...]]]
        #                           [--filter [FILTER [FILTER ...]]]
        #                           [--info [INFO [INFO ...]]]

        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg19')
        self.runCmd('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
        # self.runCmd('vtools use dbSNP')

        #you can assin any word for the suboption of --format_string
        self.assertSucc('vtools export variant --format vcf --format_string GT')
        self.assertSucc('vtools export variant --format vcf --format_string CA')
        
        #the suboptions you can update using "vtools update variant --from_file" or
        #or using "vtools update variant --set ..."
        #or import an annotation database
        #if you entered the field does not exist, then you will get error message
        
        #"THE FOLLOWING THREE TESTS SHOULD BE FAIL,NEED TO MODIFY THE CODE"
        #self.assertSucc('vtools export variant --format vcf --geno_info test')
        #self.assertSucc('vtools export variant --format vcf --var_info test')
        #self.assertSucc('vtools export variant --format vcf --var_info')
        # self.assertSucc('vtools export variant --format vcf --var_info dbSNP.molType')
        # self.assertSucc('vtools export variant --format vcf --geno_info dbSNP.valid')
        
        # #for those suboptions, you have to assign the fields
        # self.assertSucc('vtools export variant --format vcf --id dbSNP.name')
        # self.assertSucc('vtools export variant --format vcf --pos dbSNP.start')
        # self.assertSucc('vtools export variant --format vcf --ref dbSNP.refNCBI')
        # self.assertSucc('vtools export variant --format vcf --alt dbSNP.alt')
        
        # self.assertSucc('vtools export variant --format vcf --qual dbSNP.valid')
        # self.assertSucc('vtools export variant --format vcf --filter dbSNP.molType')
        # self.assertSucc('vtools export variant --format vcf --info dbSNP.valid')

        #"THE FOLLOWING FOUR TESTS SHOULD BE FAIL, NEED TO MODIFY THE CODE"
        # without name, the cmd can be passed, but no output for the field 
        self.assertSucc('vtools export variant --format vcf --qual')
        self.assertSucc('vtools export variant --format vcf --filter')
        self.assertSucc('vtools export variant --format vcf --info')
        self.assertSucc('vtools export variant --format vcf --geno')

        #--phase_sep should be used with --samples, otherwise useless 
        self.assertSucc('vtools export variant --format vcf --format_string GT ')
        self.assertSucc('vtools export variant --format vcf --samples 1 --format_string GT ')

        #I think --wildtype_code did not work in this way
        self.assertSucc('vtools export variant --format vcf --wildtype_code 0')
        self.assertSucc('vtools export variant --format vcf --samples 1 --format_string GT --wildtype_code 1')
   
    def testExportVcfSnv(self):
        'Test command export in vcf format'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        # test basic vcf output
        self.assertSucc('vtools export variant --output my.vcf')
        variants = self.runCmd('vtools output variant chr pos ref alt -d"\t"', ret='list')
        with open('my.vcf') as infile:
            exported = infile.readlines()
            exported = ['\t'.join([x.split('\t')[0], x.split('\t')[1], x.split('\t')[3], x.split('\t')[4]]) for x in exported]
            # handle multiple alternatives
            for idx, item in enumerate(exported):
                if ',' in item:
                    exported[idx] = '\t'.join(item.split(',')[:-1])
                    exported.insert(idx, exported[idx])
                    tmp = exported[idx].split('\t')[:-1]
                    tmp.extend(item.split(',')[-1])
                    exported[idx+1] = '\t'.join(tmp)
            self.assertEqual(variants, exported)
            # test output with samples
            # FIXME phase information is not properly handled
            self.assertSucc('vtools export variant --format vcf --samples 1 --format_string GT')
            self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --var_info id filter info --geno_info DP_geno')
            self.assertSucc('vtools export variant --format vcf --id id --filter filter --info info')
            self.assertSucc('vtools export variant --format vcf --id id --filter filter --info info --geno_info DP_geno --samples 1 --format_string GT:DP')
            # export selected table with selected info field
            self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --var_info DP')
            self.assertSucc('vtools select variant "DP > 300" -t highDP')
            self.assertSucc('vtools export highDP --format vcf --id id --filter filter --info DP --geno_info DP_geno --samples 1 --format_string GT:DP')

    # def testExportVcfSnv_mul(self):
    #     'Test command export with multiple sub-options'
    #     self.runCmd('vtools import vcf/CEU.vcf.gz --build hg19')
    #     self.runCmd('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
    #     self.runCmd('vtools use dbSNP')
        
    #     # You can assign a file or give the file format using this export "vcf" command 
    #     self.assertSucc('vtools export variant --format vcf --pos pos --ref ref')
    #     self.assertSucc('vtools export variant --format vcf --pos pos --ref ref --alt alt')
    #     self.assertFail('vtools export variant --header CHROM POS ID REF ALT QUAL FILTER INFO')
    #     self.assertSucc('vtools export variant --format vcf --header CHROM POS ID REF ALT QUAL FILTER INFO')
    #     self.assertSucc('vtools export variant -o  my.vcf --header CHROM POS ID REF ALT QUAL FILTER INFO')

    #     #if you assign the postion to the annotation database
    #     #the default is output the chr, pos, ref and alt from the variant table
    #     self.assertSucc('vtools export variant --format vcf --pos pos  --ref ref --alt alt')

    #     #those tests are passed. Please check the fields you want to assign to each suboption, otherwise the output will be messy.
    #     self.assertSucc('vtools export variant --format vcf --pos dbSNP.start --ref ref')
    #     self.assertSucc('vtools export variant --format vcf --pos dbSNP.start --ref dbSNP.refNCBI --alt dbSNP.alt ')
    #     self.assertSucc('vtools export variant -o my.vcf --id dbSNP.name --filter dbSNP.valid --info dbSNP.class')
    #     self.assertSucc('vtools export variant --format vcf --id dbSNP.name --filter dbSNP.valid --info dbSNP.class --samples 1 --format_string GT:DP')

    #     #output in alternative reference genome
    #     self.runCmd('vtools liftover hg18')
    #     self.assertSucc('vtools export variant --format vcf --id dbSNP.name --filter dbSNP.valid --info dbSNP.vlass --samples 1 --format_string GT:DP')

    def testHeader(self):
        'Test option header of vtools export'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertOutput('vtools export variant --format vcf --header something', 'something', 1)
        self.assertOutput('vtools export variant --format vcf --header a1 a2', 'a1\ta2', 1)
        self.assertOutput('vtools export variant --format vcf --samples 1 --header a1 a2 "%(sample_names)s"', 
            'a1\ta2\tNA06985\tNA06986\tNA06994\tNA07000\tNA07037\tNA07051\tNA07346\tNA07347\t'
            'NA07357\tNA10847\tNA10851\tNA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11881\t'
            'NA11894\tNA11918\tNA11919\tNA11920\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995\t'
            'NA12003\tNA12004\tNA12005\tNA12006\tNA12043\tNA12044\tNA12045\tNA12144\tNA12154\t'
            'NA12155\tNA12156\tNA12234\tNA12249\tNA12287\tNA12414\tNA12489\tNA12716\tNA12717\t'
            'NA12749\tNA12750\tNA12751\tNA12760\tNA12761\tNA12762\tNA12763\tNA12776\tNA12812\t'
            'NA12813\tNA12814\tNA12815\tNA12828\tNA12872\tNA12873\tNA12874', 1)
        self.assertOutput('vtools export variant --format vcf --samples \'sample_name like "NA11%"\' --header "chr pos %(sample_names)s"', 
            'chr pos NA11829\tNA11830\tNA11831\tNA11832\tNA11840\tNA11881\tNA11894'
            '\tNA11918\tNA11919\tNA11920\tNA11931\tNA11992\tNA11993\tNA11994\tNA11995', 1)

    def testExportVcfSnv_ind(self):
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools import vcf/SAMP3_complex_variants.vcf --build hg18')
        self.assertSucc('vtools update variant --from_file vcf/SAMP3_complex_variants.vcf --var_info raw_pos raw_ref raw_alt') 
        self.assertSucc('vtools output variant chr pos ref alt raw_pos raw_ref raw_alt')
        self.assertSucc('vtools export variant -o indel.vcf --pos raw_pos --ref raw_ref --alt raw_alt')
        #export indel to the alternative genome
        self.runCmd('vtools liftover hg19')
        self.assertSucc('vtools update variant --from_file vcf/SAMP3_complex_variants.vcf --var_info upstream')
        self.assertSucc('vtools export variant -o indel.vcf --pos \'pos-length(upstream)\' --ref raw_ref --alt raw_alt  --build hg19')

    def testExportANNOVAR(self):
        'Test command export in annovar input format'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg19')
        self.runCmd('vtools import vcf/SAMP3_complex_variants.vcf --build hg19')
        # self.runCmd('vtools use dbSNP')
        #the output format of 'ANNOVAR' is with one options: --comment_string
        self.assertSucc('vtools export variant --format ANNOVAR') 
        # self.assertSucc('vtools export variant --format ANNOVAR --comment_string dbSNP.valid')
        #
        # export table with non-ascii name
        # self.runCmd('vtools select variant -t "8#?"')
        # self.assertSucc('vtools export "8#?" --format ANNOVAR') 

    def testTped(self):
        'Test command export in tped format'
        #the output format of 'tped' is with two options: --name and --style
        #vtools export variant --format tped -- to query
        #
        # we should use hg18 but I do not want to download the hg18 version of dbSNP
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg19')
        self.runCmd('vtools select variant "pos = 533" -t tmp')
        self.assertSucc('vtools export variant --format tped --samples "sample_name like \'NA069%\'"')
        self.assertSucc('vtools export variant --format tped --samples 1 --style numeric')
        self.assertOutput('vtools export tmp --format tped --samples "sample_name like \'NA069%\'"', 
                          '1	.	.	533	G	G	G	G	G	G', 1)
        # self.runCmd('vtools use dbSNP')
        # self.assertSucc('vtools export variant --format tped --samples \'sample_name like "NA069%"\' --name dbSNP.name')

    # def testTpedMissingGen(self):
    #     'Test command export in tped format with missing genotype data'
    #     self.runCmd('vtools import vcf/missing_gen.vcf --build hg19')
    #     self.assertOutput('vtools export variant --format tped --samples 1', 'output/missing_gen.tped') 

if __name__ == '__main__':
    unittest.main()
