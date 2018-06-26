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
from testUtils import ProcessTestCase

class TestUpdate(ProcessTestCase):
    def testAddField(self):
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        # self.runCmd('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        self.runCmd('vtools import vcf/input_nogeno.vcf --build hg18 --sample_name input.tsv')
        self.runCmd('vtools import vcf/SAMP1.vcf')
        # no table specified
        self.assertFail('vtools update --format ../resources/format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function')
        #need a format file if you want to add field(s) into the variant table using --from_file
        self.assertFail('vtools update variant --from_file txt/annovar.txt.exonic_variant_function')
        self.assertFail('vtools update variant --format ../resources//format/ANNOVAR_exonic_variant_function')
        self.assertSucc('vtools update variant --format ../resources/format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function')
        self.assertOutput('vtools select variant "mut_type is not null" -c', '21')
        #for different version of genome
        self.assertSucc('vtools update variant --format ../resources/format/ANNOVAR_exonic_variant_function --from_file txt/annovar.txt.exonic_variant_function --build hg19')
        self.assertOutput('vtools select variant "mut_type is not null" -c', '24')

    def testUpdate(self):
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools import vcf/SAMP1.vcf --build hg18')
        # self.runCmd('vtools import --build hg18 --format fmt/basic_hg18 txt/input.tsv')
        self.runCmd('vtools import vcf/input_nogeno.vcf --build hg18 --sample_name input.tsv')
        # no hg19, but will automatically add it
        self.assertSucc('vtools update variant --format fmt/dbSNP_hg19validation --from_file txt/dbSNP_hg19validation.txt --build hg19')
        #If the file is already imported and you can add field(s) using the orginal file without --format
        self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --geno_info DP_geno')
        #you could not use another file which is not loaded into the project to update the current variant table 
        # self.assertFail('vtools update variant --from_file vcf/SAMP4_complex_variants.vcf --geno_info DP_geno')
        self.assertOutput('vtools select variant "mut_type_dbSNP is not null" -c', '172')
        self.assertOutput("vtools select variant alt_pos=753405 -o chr pos mut_type_dbSNP validation -d'\t'",
            "1\t743268\tuntranslated-5\tby-cluster,by-1000genomes\n")
     
    
    def testSampleStat(self):
        'Test command vtools update'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        self.runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        self.assertFail('vtools update')
        self.assertSucc('vtools update -h')
        self.assertFail('vtools update "num=#(alt)"')
        self.assertSucc('vtools update variant --from_stat "cnt=#(GT)" "num=#(alt)" "hom=#(hom)" "het=#(het)" "other=#(other)"')
        total = self.runCmd("vtools execute 'select sum(num) from variant'")
        hom = self.runCmd("vtools execute 'select sum(hom) from variant'")
        het = self.runCmd("vtools execute 'select sum(het) from variant'")
        other = self.runCmd("vtools execute 'select sum(other) from variant'")
        cnt = self.runCmd("vtools execute 'select cnt from variant LIMIT 1'")
        for to,ho,he,ot in zip(total.strip().split('\n'), hom.strip().split('\n'), het.strip().split('\n'), other.strip().split('\n')):
            self.assertEqual(int(to), int(ho)*2+int(he)+int(ot))
        self.assertSucc('vtools update CEU --from_stat "CEU_num=#(alt)" -s "filename like \'%CEU%\'"')
        self.assertOutput("vtools execute 'select sum(CEU_num) from variant'", '6383\n', partial=True)
        self.assertSucc('vtools update CEU --from_stat "CEU_num=#(alt)" "CEU_hom=#(hom)" "CEU_het=#(het)" "CEU_other=#(other)"  --samples "filename like \'%CEU%\'"')
        self.assertSucc('vtools update CEU --from_stat "CEU_cases_het=#(het)" --samples "filename like \'%CEU%\' and aff=\'2\'"')
        self.assertOutput("vtools execute 'select sum(CEU_cases_het) from variant'", '1601\n', partial=True)
        self.assertSucc('vtools update CEU --from_stat "CEU_strls_het=#(het)" -s "filename like \'%CEU%\' and aff=\'1\'"')
     
    def testMaf(self):
        'Test computation of MAF'
        # all females
        self.runCmd('vtools import vcf/chromX.vcf.gz --build hg18')
        self.runCmd('vtools phenotype --set "sex=2" --samples 1')
        self.assertSucc('vtools update variant --from_stat "cnt0=#(GT)" "num0=#(alt)" "hom0=#(hom)" "het0=#(het)" "other0=#(other)" "maf0=maf()"')
        cnt = self.runCmd("vtools execute 'select cnt0 from variant LIMIT 10'")
        num = self.runCmd("vtools execute 'select num0 from variant LIMIT 10'")
        maf = self.runCmd("vtools execute 'select maf0 from variant LIMIT 10'")
        for c,n,m in zip(cnt.strip().split('\n'), num.strip().split('\n'), maf.strip().split('\n')):
            value = float(n)/float(c)/2.0
            value = value if value < 0.5 else 1 - value
            self.assertAlmostEqual(float(m), value)
        # all males
        self.runCmd('vtools phenotype --set "sex=1" --samples 1')
        self.assertSucc('vtools update variant --from_stat "cnt1=#(GT)" "num1=#(alt)" "hom1=#(hom)" "het1=#(het)" "other1=#(other)" "maf1=maf()"')
        cnt = self.runCmd("vtools execute 'select cnt1 from variant LIMIT 10'")
        num = self.runCmd("vtools execute 'select num1 from variant LIMIT 10'")
        maf = self.runCmd("vtools execute 'select maf1 from variant LIMIT 10'")
        # for c,n,m in zip(cnt.strip().split('\n'), num.strip().split('\n'), maf.strip().split('\n')):
        #     value = float(n)/float(c)
        #     value = value if value < 0.5 else 1 - value
        #     self.assertAlmostEqual(float(m), value)
        
    
    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","HDF5 version is not implemented for this test")
    def testGenotypeSumStats(self):
        print(os.getenv("STOREMODE"))
        'Test command vtools update min/max/sum/mean_FIELD'
        self.runCmd('vtools import --format fmt/missing_gen vcf/missing_gen.vcf --build hg19')
        # non-existing field, should fail
        # self.assertFail('vtools update variant --from_stat "max_gq=max(GQ1)" "min_gq=min(GQ)"')
        self.assertSucc("vtools update variant --from_stat 'total=#(GT)' 'num=#(alt)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)' \
                        'minDP=min(GD)' 'maxDP=max(GD)' 'meanDP=avg(GD)' 'minGQv=min(GQ)' 'maxGQv=max(GQ)' 'meanGQv=avg(GQ)'")
        self.assertOutput('vtools output variant maxGQv minGQv meanGQv --precision 4', 'output/update_sum_stat.txt')
        self.assertSucc('vtools update variant --from_stat "total_dp=sum(GD)"')
        self.assertProj(info={'total_dp': ['None', 'None', 'None', '60', '7', '4']})
        # ffilter out variants having GQ less than 4,
        # then for each remining variant count the total number of alt genotypes across all samples
        self.assertSucc('vtools update variant --from_stat "gq_ge_4=#(alt)" --genotype "GQ >= 4"')
        self.assertProj(info={'gq_ge_4': ['0', '0', '0', '0', '3', '1']})
        #
        self.assertSucc('vtools update variant --from_stat "missing=#(missing)" "wt=#(wtGT)" "mut=#(mutGT)" "total=#(GT)"')
        self.assertSucc('vtools update variant --set "res=total-wt-hom-het-other"')
        self.assertProj(info={'res': ['0']*6})
        self.assertSucc('vtools update variant --set "res=num-2*hom-het-other"')
        self.assertProj(info={'res': ['0']*6})
        self.assertSucc('vtools update variant --set "res=4-missing-total"')
        self.assertProj(info={'res': ['0']*6})
        self.assertSucc('vtools update variant --set "res=mut-hom-het-other"')
        self.assertProj(info={'res': ['0']*6})


    # def testGenotypeSumStats(self):
    #     'Test command vtools update min/max/sum/mean_FIELD'
    #     self.runCmd('vtools import  vcf/missing_gen_hdf5.vcf --build hg19 --geno_info PL_geno,DP_geno,GQ_geno')
    #     # non-existing field, should fail
    #     # self.assertFail('vtools update variant --from_stat "max_gq=max(GQ1)" "min_gq=min(GQ)"')
    #     self.assertSucc("vtools update variant --from_stat 'total=#(GT)' 'num=#(alt)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)' \
    #                     'minDP=min(DP_geno)' 'maxDP=max(DP_geno)' 'meanDP=avg(DP_geno)' 'minGQv=min(GQ_geno)' 'maxGQv=max(GQ_geno)' 'meanGQv=avg(GQ_geno)'")
    #     self.assertOutput('vtools output variant maxGQv minGQv meanGQv --precision 4', 'output/update_sum_stat.txt')
    #     self.assertSucc('vtools update variant --from_stat "total_dp=sum(DP_geno)"')
    #     self.assertProj(info={'total_dp': ['None', 'None', 'None', '60', '7', '4']})
    #     # ffilter out variants having GQ less than 4,
    #     # then for each remining variant count the total number of alt genotypes across all samples
    #     self.assertSucc('vtools update variant --from_stat "gq_ge_4=#(alt)" --genotype "GQ >= 4"')
    #     self.assertProj(info={'gq_ge_4': ['0', '0', '0', '0', '3', '1']})
    #     #
    #     self.assertSucc('vtools update variant --from_stat "missing=#(missing)" "wt=#(wtGT)" "mut=#(mutGT)" "total=#(GT)"')
    #     self.assertSucc('vtools update variant --set "res=total-wt-hom-het-other"')
    #     self.assertProj(info={'res': ['0']*6})
    #     self.assertSucc('vtools update variant --set "res=num-2*hom-het-other"')
    #     self.assertProj(info={'res': ['0']*6})
    #     self.assertSucc('vtools update variant --set "res=4-missing-total"')
    #     self.assertProj(info={'res': ['0']*6})
    #     self.assertSucc('vtools update variant --set "res=mut-hom-het-other"')
    #     self.assertProj(info={'res': ['0']*6})
        

    def testSet(self):
        'Testing vtools update --set'
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('''vtools select variant "chr='1'" -t chr1''')
        self.runCmd('vtools update chr1 --from_stat "num=#(alt)" "total=#(GT)"')
        self.runCmd('vtools update chr1 --set "ratio=num/(total*1.0)"')
        self.runCmd('vtools select variant "ratio is NULL" -t other')
        self.assertProj(numOfVariants={'other': 100})
        self.assertSucc('vtools show fields')

    def testGenoAnnoSet(self):
        'Testing command vtools update --set'
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')    
        self.assertSucc("vtools update variant --from_stat 'total=#(GT)' 'num=#(alt)' 'het=#(het)' 'hom=#(hom)' 'other=#(other)' 'minDP=min(DP)' 'maxDP=max(DP)' 'meanDP=avg(DP)' 'minGQv=min(GQ)' 'maxGQv=max(GQ)' 'meanGQv=avg(GQ)'")
        self.assertSucc('vtools update variant --set "maf=num/(total*2.0)"')
        self.assertSucc('vtools output variant chr pos total num maf -l 10')
        #we can set the fields from the annotation file
        self.assertSucc('vtools liftover hg19')
        # evs does not exist for some reason
        self.assertSucc('vtools use refGene')
        self.assertSucc('vtools update variant --set evs_gene=refGene.name2')
        self.assertOutput('vtools execute "select chr,pos, ref, alt, evs_gene from variant where evs_gene is not null"',
            '22\t49524956\tG\tA\tACR\n', partial=True)
        self.assertOutput(('vtools select variant "evs_gene is not NULL" -c'), '121\n')

    def testSetMultiValues(self):
        self.runCmd('vtools import txt/variants.txt --format basic --build hg19')
        self.runCmd('vtools use refGene')
        self.runCmd('vtools update variant --set "gname=refGene.name2"')
        lines = self.runCmd('vtools output variant chr pos ref alt gname refGene.name2 -d"\t"', ret='list')
        # does not count duplicate lines
        line_count = len(set(lines))
        for line in lines:
            values = line.split('\t')
            self.assertEqual(values[-2], values[-1])
        self.runCmd('vtools update variant --set "vname=refGene.name"')
        lines = self.runCmd('vtools output variant chr pos ref alt vname refGene.name', ret='list')
        # there should not be more lines
        self.assertEqual(line_count, len(set(lines)))
        
if __name__ == '__main__':
    unittest.main()
