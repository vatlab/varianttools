#!/usr/bin/env python
#
# $File: test_update.py $
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
from testUtils import ProcessTestCase, runCmd, numOfSample, numOfVariant, outputOfCmd, initTest, output2list

class TestUpdate(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
    def removeProj(self):
        runCmd('vtools remove project')

    def testAddfield(self):
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        runCmd('vtools import vcf/SAMP1.vcf')
        # no table specified
        self.assertFail('vtools update --format ../format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function')
        #need a format file if you want to add field(s) into the variant table using --from_file
        self.assertFail('vtools update variant --from_file txt/annovar.txt.exonic_variant_function')
        self.assertFail('vtools update variant --format ../format/ANNOVAR_exonic_variant_function')
        self.assertSucc('vtools update variant --format ../format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type is not null" -c'), '78\n')
        #for different version of genome
        self.assertSucc('vtools update variant --format ../format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function --build hg19')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type is not null" -c'), '81\n')

    def testUpdate(self):
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import vcf/SAMP1.vcf --build hg18')
        runCmd('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        # no hg19
        self.assertFail('vtools update variant --format fmt/dbSNP_hg19validation --from_file txt/dbSNP_hg19validation.txt --build hg19')
        runCmd('vtools liftover hg19')
        # a variant can have multiple entries -- thus 175 variants in variant table are updated from the 198 imported records
        self.assertSucc('vtools update variant --format fmt/dbSNP_hg19validation --from_file txt/dbSNP_hg19validation.txt --build hg19')
        #If the file is already imported and you can add field(s) using the orginal file without --format
        self.assertFail('vtools')   
        self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --geno_info DP_geno')
        #you could not use another file which is not loaded into the project to update the current variant table 
        self.assertFail('vtools update variant --from_file vcf/SAMP4_complex_variants.vcf --geno_info DP_geno')
        self.assertEqual(outputOfCmd('vtools select variant "mut_type_dbSNP is not null" -c'), '175\n')
        self.assertOutput("vtools select variant alt_pos=753405 -o chr pos mut_type_dbSNP validation -d'\t'", "1\t743268\tuntranslated-5\tby-cluster,by-1000genomes\n")
        
    def testSampleStat(self):
        'Test command vtools update'
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        self.assertFail('vtools update')
        self.assertSucc('vtools update -h')
        self.assertFail('vtools update "num=#(alt)"')
        self.assertSucc('vtools update variant --from_stat "cnt=#(GT)" "num=#(alt)" "hom=#(hom)" "het=#(het)" "other=#(other)"')
        total = int(outputOfCmd("vtools execute 'select sum(num) from variant'").split('\n')[0])
        hom = int(outputOfCmd("vtools execute 'select sum(hom) from variant'").split('\n')[0])
        het = int(outputOfCmd("vtools execute 'select sum(het) from variant'").split('\n')[0])
        other = int(outputOfCmd("vtools execute 'select sum(other) from variant'").split('\n')[0])
        cnt = int(outputOfCmd("vtools execute 'select cnt from variant LIMIT 1'").split('\n')[0])
        self.assertEqual(cnt, 60)
        self.assertEqual(total, hom*2+het+other)
        self.assertSucc('vtools update CEU --from_stat "CEU_num=#(alt)" -s "filename like \'%CEU%\'"')
        self.assertEqual(int(outputOfCmd("vtools execute 'select sum(CEU_num) from variant'").split('\n')[0]), 6383)
        self.assertSucc('vtools update CEU --from_stat "CEU_num=#(alt)" "CEU_hom=#(hom)" "CEU_het=#(het)" "CEU_other=#(other)"  --samples "filename like \'%CEU%\'"')
        self.assertSucc('vtools update CEU --from_stat "CEU_cases_het=#(het)" --samples "filename like \'%CEU%\' and aff=\'2\'"')
        self.assertEqual(int(outputOfCmd("vtools execute 'select sum(CEU_cases_het) from variant'").split('\n')[0]), 1601)
        self.assertSucc('vtools update CEU --from_stat "CEU_strls_het=#(het)" -s "filename like \'%CEU%\' and aff=\'1\'"')
        
    def testMaf(self):
        'Test computation of MAF'
        # all females
        runCmd('vtools import vcf/chromX.vcf.gz --build hg18')
        runCmd('vtools phenotype --set "sex=2" --samples 1')
        self.assertSucc('vtools update variant --from_stat "cnt0=#(GT)" "num0=#(alt)" "hom0=#(hom)" "het0=#(het)" "other0=#(other)" "maf0=maf()"')
        cnt = outputOfCmd("vtools execute 'select cnt0 from variant LIMIT 10'").split('\n')
        num = outputOfCmd("vtools execute 'select num0 from variant LIMIT 10'").split('\n')
        maf = outputOfCmd("vtools execute 'select maf0 from variant LIMIT 10'").split('\n')
        for i in range(10):
            value = float(num[i])/float(cnt[i])/2.0
            value = value if value < 0.5 else 1 - value
            self.assertAlmostEqual(float(maf[i]), value)
        # all males
        runCmd('vtools phenotype --set "sex=1" --samples 1')
        self.assertSucc('vtools update variant --from_stat "cnt1=#(GT)" "num1=#(alt)" "hom1=#(hom)" "het1=#(het)" "other1=#(other)" "maf1=maf()"')
        cnt = outputOfCmd("vtools execute 'select cnt1 from variant LIMIT 10'").split('\n')
        num = outputOfCmd("vtools execute 'select num1 from variant LIMIT 10'").split('\n')
        maf = outputOfCmd("vtools execute 'select maf1 from variant LIMIT 10'").split('\n')
        for i in range(10):
            value = float(num[i])/float(cnt[i])
            value = value if value < 0.5 else 1 - value
            self.assertAlmostEqual(float(maf[i]), value)
        
       
    def testGenotypeSumStats(self):
        'Test command vtools update min/max/sum/mean_FIELD'
        runCmd('vtools import --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        # non-existing field, should fail
        self.assertFail('vtools update variant --from_stat "max_gq=max(GQ1)" "min_gq=min(GQ)"')
        self.assertSucc("vtools update variant --from_stat 'total=#(GT)' 'num=#(alt)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)' \
                            'minDP=min(GD)' 'maxDP=max(GD)' 'meanDP=avg(GD)' 'minGQv=min(GQ)' 'maxGQv=max(GQ)' 'meanGQv=avg(GQ)'")
        out = output2list('vtools output variant maxGQv minGQv meanGQv -d"\t"')
        self.assertEqual(out[0], '.\t.\t.')
        self.assertTrue(out[3].startswith('100\t15\t69.3333'))
        self.assertTrue(out[4].startswith('6\t3\t4.0'))
        self.assertTrue(out[5].startswith('4\t3\t3.33333'))
        self.assertSucc('vtools update variant --from_stat "total_dp=sum(GD)"')
        self.assertEqual(output2list('vtools output variant total_dp -d"\t"'), ['.', '.', '.', '60', '7', '4'])
        # ffilter out variants having GQ less than 4,
        # then for each remining variant count the total number of alt genotypes across all samples
        self.assertSucc('vtools update variant --from_stat "gq_ge_4=#(alt)" --genotype "GQ >= 4"')
        self.assertEqual(output2list('vtools output variant gq_ge_4 -d"\t"'), ['0', '0', '0', '0', '3', '1'])
        #
        self.assertSucc('vtools update variant --from_stat "missing=#(missing)" "wt=#(wtGT)" "mut=#(mutGT)" "total=#(GT)"')
        self.assertSucc('vtools update variant --set "res=total-wt-hom-het-other"')
        self.assertEqual(output2list('vtools output variant res'), ['0']*6)
        self.assertSucc('vtools update variant --set "res=num-2*hom-het-other"')
        self.assertEqual(output2list('vtools output variant res'), ['0']*6)
        self.assertSucc('vtools update variant --set "res=4-missing-total"')
        self.assertEqual(output2list('vtools output variant res'), ['0']*6)
        self.assertSucc('vtools update variant --set "res=mut-hom-het-other"')
        self.assertEqual(output2list('vtools output variant res'), ['0']*6)
        

    def testSet(self):
        'Testing vtools update --set'
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        runCmd('''vtools select variant "chr='1'" -t chr1''')
        runCmd('vtools update chr1 --from_stat "num=#(alt)" "total=#(GT)"')
        runCmd('vtools update chr1 --set "ratio=num/(total*1.0)"')
        runCmd('vtools select variant "ratio is NULL" -t other')
        lines = output2list('vtools output other chr pos ref alt ratio')
        self.assertEqual(len(lines), 100)
        self.assertSucc('vtools show fields')

    def testGenoAnnoSet(self):
        'Testing command vtools update --set'
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')    
        self.assertSucc("vtools update variant --from_stat 'total=#(GT)' 'num=#(alt)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)' 'minDP=min(GD)' 'maxDP=max(GD)' 'meanDP=avg(GD)' 'minGQv=min(GQ)' 'maxGQv=max(GQ)' 'meanGQv=avg(GQ)'")
        self.assertSucc('vtools update variant --set "maf=num/(total*2.0)"')
        self.assertSucc('vtools output variant chr pos total num maf -l 10')
        #we can set the fields from the annotation file
        self.assertSucc('vtools liftover hg19')
        self.assertSucc('vtools use ESP')
        self.assertSucc('vtools update variant --set evs_gene=ESP.Genes')
        self.assertOutput(('vtools execute "select chr,pos, ref, alt, evs_gene from variant\
                         where evs_gene is not null"'), '22\t49524956\tG\tA\tACR\n')
        self.assertOutput(('vtools select variant "evs_gene is not NULL" -c'), '1\n')

    def testSetMultiValues(self):
        runCmd('vtools import txt/variants.txt --format basic --build hg19')
        runCmd('vtools use refGene')
        runCmd('vtools update variant --set "gname=refGene.name2"')
        lines = output2list('vtools output variant chr pos ref alt gname refGene.name2 -d"\t"')
        # does not count duplicate lines
        line_count = len(set(lines))
        for line in lines:
            values = line.split('\t')
            self.assertEqual(values[-2], values[-1])
        runCmd('vtools update variant --set "vname=refGene.name"')
        lines = output2list('vtools output variant chr pos ref alt vname refGene.name')
        # there should not be more lines
        self.assertEqual(line_count, len(set(lines)))
        
if __name__ == '__main__':
    unittest.main()
