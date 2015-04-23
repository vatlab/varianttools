#!/usr/bin/env python
#
# $File: test_select.py $
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd, output2list

class TestSelect(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(6)
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        runCmd('vtools update variant --from_stat "num=#(alt)" "hom=#(hom)" "het=#(het)" "other=#(other)"')
        runCmd('vtools update CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --from_stat "CEU_cases_het=#(het)"')

    def removeProj(self):
        runCmd('vtools remove project')
    
    def testSelect(self):
        'Test command vtools select'
        self.assertFail('vtools select')
        self.assertSucc('vtools select -h')
        # Variant table non_existing_variant does not exist.
        self.assertFail('vtools select non_existing_variant')
        # Neither --to_table and --output/--count is specified. Nothing to do.
        self.assertFail('vtools select variant \'testNSFP.non_existing_item is not null\'')
        # Neither --to_table and --output/--count is specified. Nothing to do.
        self.assertFail('vtools select variant \'testNSFP.chr is not null\'')
        self.assertOutput("vtools select variant -c", '915\n')
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' -t "ns"')
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' --output chr pos ref alt')
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' --output variant_id testNSFP.CHBJPT_total_lc')
        # Existing table ns_input is renamed to ns_input_Aug06_161348. The command below is equivalent to the former two commands.
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' -t ns')
        #the follosing test should add --output (or -o)
        self.assertSucc('vtools select variant -o \'max(testNSFP.polyphen2_score)\'')
        self.assertSucc('vtools select variant -o \'max(testNSFP.polyphen2_score)\' --header max')
        self.assertOutput("vtools execute 'select count(*) from ns'", '7\n')
        self.assertOutput("vtools select ns -c", '7\n')
        # use strange characters 
        self.assertSucc('vtools select variant \'testNSFP.chr is not null\' -t "* ns@"')
        self.assertSucc('vtools select "* ns@" \'testNSFP.chr is not null\' -t "ns@sub"')
        self.assertTrue('* ns@' in '\n'.join(output2list('vtools show tables')))
        self.assertSucc('vtools show table "* ns@"')
        self.assertSucc('vtools show table "ns@sub"')
        
    def testSelectSample(self):
        self.assertOutput("vtools select variant --samples 'filename like \"%input.tsv\"' -c", '338\n')
        self.assertOutput("vtools select variant --samples 'filename like \"%CEU.vcf.gz\" ' -c", '288\n')
        self.assertSucc('vtools select variant "testNSFP.chr is not null" --samples "filename like \'%input.tsv\'" -t ns_input')
        #nsfp = output2list('vtools execute "select chr, hg18pos, ref, alt from testNSFP"')
        #variantid = [output2list('vtools execute "select variant_id from variant where chr={0} and pos={1} and ref={2} and alt={3}"'.format(x.split()[0], x.split()[1], repr(x.split()[2]), repr(x.split()[3]))) for x in nsfp]
        #print variantid
        self.assertOutput("vtools select ns_input -c", '7\n')
        self.assertSucc('vtools select ns_input \'genename = "PLEKHN1"\'  -t plekhn1')
        self.assertOutput("vtools select plekhn1 -c", '6\n')
        self.assertSucc('vtools select plekhn1 "polyphen2_score>0.9 and sift_score>0.9" -t d_plekhn1')
        self.assertOutput("vtools select d_plekhn1 -c", '5\n')
        self.assertSucc('vtools select variant "testNSFP.chr is not null" --samples "aff=1" -t ns_aff')
        self.assertOutput("vtools select ns_aff -c", '0\n')
        #
        self.assertSucc('vtools select variant --samples "aff=\'1\' and BMI<20" -t ns3')
        namelist = output2list('vtools execute "select sample_id from sample where aff=1 and BMI<20"')
        variantlist = [output2list('vtools execute "select variant_id from genotype_{} where GT <> 0"'.format(x)) for x in namelist]
        variantlist = [x for y in variantlist for x in y]
        lv = str(len(set(variantlist)))        
        self.assertOutput("vtools select ns3 -c", '{}\n'.format(lv))
        self.assertOutput("vtools execute 'select count(*) from sample where aff=1 and BMI<20'", '10\n')
        self.assertSucc('vtools select variant --samples "aff=\'1\'" "BMI<20" -t ns3')
        self.assertOutput("vtools select ns3 -c", '{}\n'.format(lv))
        #
        self.assertSucc('vtools select variant --samples "aff=\'1\' or BMI<20" -t ns2')
        namelist = output2list('vtools execute "select sample_id from sample where aff=1 or BMI<20"')
        variantlist = [output2list('vtools execute "select variant_id from genotype_{} where GT <> 0"'.format(x)) for x in namelist]
        variantlist = [x for y in variantlist for x in y]
        lv = str(len(set(variantlist)))
        self.assertOutput("vtools select ns2 -c", '{}\n'.format(lv))
        self.assertOutput("vtools execute 'select count(*) from sample where aff=1 or BMI<20'", '41\n')
        #
        self.assertSucc('vtools select variant "testNSFP.chr is not null" "genename=\'PLEKHN1\'" "polyphen2_score>0.9 or sift_score>0.9" -t d_plekhn1')
        #
        self.assertSucc('vtools select variant --samples "sample_name like \'NA0%\'" -t NA0')
        namelist = output2list('vtools execute "select sample_id from sample where sample_name like \'NA0%\'"')
        variantlist = [output2list('vtools execute "select variant_id from genotype_{} where gt <> 0"'.format(x)) for x in namelist]
        variantlist = [x for y in variantlist for x in y]
        lv = str(len(set(variantlist)))
        self.assertOutput("vtools select NA0 -c", '{}\n'.format(lv))
        self.assertOutput("vtools execute 'select count(*) from sample where sample_name like \"NA0%\"'", '9\n')
        self.assertSucc('vtools select CEU -s "BMI<18.5" -t Underweight')

    def testSelectSampleWithWildtypeGenotype(self):
        runCmd('vtools import vcf/with_wildtype.vcf --sample_name WT')
        # original 989 variants but some of them have only wildtype genotype
        self.assertOutput('''vtools select variant --samples "sample_name='WT'" -c''', "934\n")


    def testSelectLargeSample(self):
        runCmd('vtools import vcf/500SAMP.vcf')
        self.assertSucc('vtools select variant --samples "sample_name like \'SAM%\'" -c')


    def testFunctionLeast(self):
        'Testing a self-defined sql function least'
        self.assertSucc('vtools update variant --set "nil=NULL"')
        # the least function will ignore value from nil
        self.assertSucc('vtools update variant --set "lst=least_not_null(hom, het, nil)"')
        counts = output2list('vtools output variant hom het lst')
        for line in counts:
            values = [int(x) for x in line.split()]
            self.assertEqual(min(values[0], values[1]), values[2])

    
if __name__ == '__main__':
    unittest.main()
