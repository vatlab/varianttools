#!/usr/bin/env python
#
# $File: test_compare.py $
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
import unittest
from testUtils import ProcessTestCase

class TestCompare(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        if self.storeMode=="sqlite" and os.path.isfile('TestCompare_sqlite.tar.gz'):
            self.runCmd('vtools admin --load_snapshot TestCompare_sqlite.tar.gz')
        elif self.storeMode=="hdf5" and os.path.isfile('TestCompare_hdf5.tar.gz'):
            self.runCmd('vtools admin --load_snapshot TestCompare_hdf5.tar.gz')
        else:
            self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
            self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
            self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
            self.runCmd('vtools use ann/testNSFP.ann')
            self.runCmd('vtools select variant \'testNSFP.chr is not null\' -t ns')
            self.runCmd('vtools select ns \'sift_score > 0.95\' -t ns_damaging')
            self.runCmd('vtools select ns \'genename = "PLEKHN1"\' -t plekhn1')
            self.runCmd('vtools select plekhn1 "polyphen2_score>0.9 or sift_score>0.9" -t d_plekhn1')
            if self.storeMode=="sqlite":
                self.runCmd('vtools admin --save_snapshot TestCompare_sqlite.tar.gz "Snapshot for testing comapre command"')
            elif self.storeMode=="hdf5":
                self.runCmd('vtools admin --save_snapshot TestCompare_hdf5.tar.gz "Snapshot for testing comapre command"')


    def testCompare(self):
        'Test command vtools compare'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertOutput('vtools compare plekhn1 d_plekhn1 ', '431\t0\t996\t1427')
        # error: argument --A_and_B: expected one argument
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection unique')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection common')
        # WARNING: Existing table common is renamed to common_Aug09_170022. No error 
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection common')
        self.assertProj(numOfVariants={'common': 572})
        self.assertSucc('vtools compare d_plekhn1 ns --union AorB')
        self.assertProj(numOfVariants={'AorB': 1446})
        self.assertSucc('vtools compare d_plekhn1 ns --intersection AandB')
        self.assertProj(numOfVariants={'AandB': 996})
        self.assertSucc('vtools compare d_plekhn1 ns --difference AdiffB')
        self.assertProj(numOfVariants={'AdiffB': 0})
        self.assertSucc('vtools compare ns d_plekhn1 --difference BdiffA')
        self.assertProj(numOfVariants={'BdiffA': 450})
        # use both options in one command should be allowed
        self.assertSucc('vtools compare d_plekhn1 plekhn1  --union A_OR_B')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --union "KK@"')
        self.assertProj(hasTable='KK@')
        self.runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare "d@p" ns_damaging --intersection "K&K@"')
        self.assertProj(hasTable='K&K@')
        # handling of wildcard names
        self.assertSucc('vtools compare \'ns*\' --union u')
        self.assertProj(hasTable='u')

    def testCompareExpression1(self):
        'Test command vtools compare (form 1)'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertOutput('vtools compare plekhn1 d_plekhn1 ', ['431\t0\t996\t1427'])
        # error: argument --A_and_B: expected one argument
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "_1 & _2"')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "common=_1 & _2"')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "common=_1 & _2"')
        self.assertProj(numOfVariants={'common': 572})
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AorB=_1 | _2"')
        self.assertProj(numOfVariants={'AorB': 1446})
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AandB=_1 & _2"')
        self.assertProj(numOfVariants={'AandB': 996})
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AdiffB=_1 - _2"')
        self.assertProj(numOfVariants={'AdiffB': 0})
        self.assertSucc('vtools compare ns d_plekhn1 --expression "BdiffA=_1 - _2"')
        self.assertProj(numOfVariants={'BdiffA': 450})
        # use both options in one command should be allowed
        self.assertSucc('vtools compare d_plekhn1 plekhn1  --expression "A_OR_B=_1 | _2"')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "\'KK@\'=_1 | _2"')
        self.assertProj(hasTable='KK@')
        self.runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare "d@p" ns_damaging --expression "\'K&K@\'=_1 & _2"')
        self.assertProj(hasTable='K&K@')

    def testCompareExpression2(self):
        'Test command vtools compare (form 2)'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertOutput('vtools compare plekhn1 d_plekhn1 ', ['431\t0\t996\t1427'], partial=True)
        # error: argument --A_and_B: expected one argument
        self.assertSucc('vtools compare --expression "d_plekhn1 & ns_damaging"')
        self.assertSucc('vtools compare --expression "common=d_plekhn1 & ns_damaging"')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare --expression "common=d_plekhn1 & ns_damaging"')
        self.assertProj(numOfVariants={'common': 572})
        self.assertSucc('vtools compare --expression "AorB=d_plekhn1 | ns"')
        self.assertProj(numOfVariants={'AorB': 1446})
        self.assertSucc('vtools compare --expression "AandB=d_plekhn1 & ns"')
        self.assertProj(numOfVariants={'AandB': 996})
        self.assertSucc('vtools compare --expression "AdiffB=d_plekhn1 -ns"')
        self.assertProj(numOfVariants={'AdiffB': 0})
        self.assertSucc('vtools compare --expression "BdiffA=ns-d_plekhn1"')
        self.assertProj(numOfVariants={'BdiffA': 450})
        # use both options in one command should be allowed
        self.assertSucc('vtools compare --expression "A_OR_B=d_plekhn1|plekhn1"')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare --expression "\'KK@\'=d_plekhn1 |ns_damaging"')
        self.assertProj(hasTable='KK@')
        self.runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare  --expression "\'K&K@\'=\'d@p\' & ns_damaging"')
        self.assertProj(hasTable='K&K@')

    def testSimpleProj(self):
        'Test compare by site'
        # do not use the default project
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.runCmd('vtools import vcf/compare.vcf --build hg18')
        # if self.storeMode=="sqlite":
        #     self.runCmd('vtools execute "DELETE FROM genotype.genotype_1 WHERE GT=0"')
        #     self.runCmd('vtools execute "DELETE FROM genotype.genotype_2 WHERE GT=0"')
        # elif self.storeMode=="hdf5":
        #     self.runCmd('vtools remove genotypes "GT_geno is nan"')
        self.runCmd('vtools select variant "variant_id in (1, 2)" -t T1')
        self.runCmd('vtools select variant "variant_id in (1, 3, 4)" -t T2')
        # variant
        self.assertOutput('vtools compare T1 T2 ', '1\t2\t1\t4\n')
        self.assertOutput('vtools compare T1 T2 --difference', '1')
        self.assertOutput('vtools compare T2 T1 --difference', '2')
        self.assertOutput('vtools compare T1 T2 --intersection', '1')
        self.assertOutput('vtools compare T1 T2 --union', '4')
        # sample
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 ', '1\t1\t1\t3\n')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --difference ', '1')
        self.assertOutput('vtools compare variant --samples SAMP2 SAMP1 --difference ', '1')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --intersection ', '1')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --union', '3')
        # site
        self.assertOutput('vtools compare T1 T2 --mode site ', '0\t1\t3\t4\n')
        self.assertOutput('vtools compare T1 T2 --mode site --difference ', '0')
        self.assertOutput('vtools compare T2 T1 --mode site --difference ', '1')
        self.assertOutput('vtools compare T1 T2 --mode site --intersection ', '3')
        self.assertOutput('vtools compare T1 T2 --mode site --union ', '4')



if __name__ == '__main__':
    unittest.main()
