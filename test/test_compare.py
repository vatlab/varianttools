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
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, initTest, output2list, outputOfCmd

class TestCompare(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(10)

    def removeProj(self):
        runCmd('vtools remove project')
        
    def testCompare(self):
        'Test command vtools compare'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertEqual(output2list('vtools compare plekhn1 d_plekhn1 '), ['0\t0\t6\t6'])
        # error: argument --A_and_B: expected one argument
        self.assertFail('vtools compare d_plekhn1 ns_damaging --intersection')
        self.assertFail('vtools compare d_plekhn1 ns_damaging --intersection unique')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection common')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --intersection common')
        self.assertEqual(output2list("vtools execute 'select * from common'"), ['880', '881', '882', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --union AorB')
        self.assertEqual(output2list("vtools execute 'select * from AorB'"), ['714', '869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --intersection AandB')
        self.assertEqual(output2list("vtools execute 'select * from AandB'"), ['869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --difference AdiffB')
        self.assertEqual(output2list("vtools execute 'select * from AdiffB'"), [])
        self.assertSucc('vtools compare ns d_plekhn1 --difference BdiffA')
        self.assertEqual(output2list("vtools execute 'select * from BdiffA'"), ['714'])
        # use both options in one command should be allowed
        self.assertSucc('vtools compare d_plekhn1 plekhn1  --union A_OR_B')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --union "KK@"')
        self.assertTrue('KK@' in outputOfCmd('vtools show tables'))
        runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare "d@p" ns_damaging --intersection "K&K@"')
        self.assertTrue('K&K@' in outputOfCmd('vtools show tables'))
        #
        # handling of wildcard names
        self.assertSucc('vtools compare \'ns*\' --union u')
        self.assertTrue('u' in outputOfCmd('vtools show tables'))

    def testCompareExpression1(self):
        'Test command vtools compare (form 1)'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertEqual(output2list('vtools compare plekhn1 d_plekhn1 '), ['0\t0\t6\t6'])
        # error: argument --A_and_B: expected one argument
        self.assertFail('vtools compare d_plekhn1 ns_damaging --expression "_1 & _2"')
        self.assertFail('vtools compare d_plekhn1 ns_damaging --expression "unique=_1 & _2"')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "common=_1 & _2"')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "common=_1 & _2"')
        self.assertEqual(output2list("vtools execute 'select * from common'"), ['880', '881', '882', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AorB=_1 | _2"')
        self.assertEqual(output2list("vtools execute 'select * from AorB'"), ['714', '869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AandB=_1 & _2"')
        self.assertEqual(output2list("vtools execute 'select * from AandB'"), ['869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --expression "AdiffB=_1 - _2"')
        self.assertEqual(output2list("vtools execute 'select * from AdiffB'"), [])
        self.assertSucc('vtools compare ns d_plekhn1 --expression "BdiffA=_1 - _2"')
        self.assertEqual(output2list("vtools execute 'select * from BdiffA'"), ['714'])
        # use both options in one command should be allowed
        self.assertSucc('vtools compare d_plekhn1 plekhn1  --expression "A_OR_B=_1 | _2"')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --expression "\'KK@\'=_1 | _2"')
        self.assertTrue('KK@' in outputOfCmd('vtools show tables'))
        runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare "d@p" ns_damaging --expression "\'K&K@\'=_1 & _2"')
        self.assertTrue('K&K@' in outputOfCmd('vtools show tables'))

    def testCompareExpression2(self):
        'Test command vtools compare (form 2)'
        # fix me: the following only test the case with two tables, should also test for one and more than two tables
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertEqual(output2list('vtools compare plekhn1 d_plekhn1 '), ['0\t0\t6\t6'])
        # error: argument --A_and_B: expected one argument
        self.assertFail('vtools compare --expression "d_plekhn1 & ns_damaging"')
        self.assertFail('vtools compare --expression "unique=d_plekhn1 & ns_damaging"')
        self.assertSucc('vtools compare --expression "common=d_plekhn1 & ns_damaging"')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare --expression "common=d_plekhn1 & ns_damaging"')
        self.assertEqual(output2list("vtools execute 'select * from common'"), ['880', '881', '882', '893'])
        self.assertSucc('vtools compare --expression "AorB=d_plekhn1 | ns"')
        self.assertEqual(output2list("vtools execute 'select * from AorB'"), ['714', '869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare --expression "AandB=d_plekhn1 & ns"')
        self.assertEqual(output2list("vtools execute 'select * from AandB'"), ['869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare --expression "AdiffB=d_plekhn1 -ns"')
        self.assertEqual(output2list("vtools execute 'select * from AdiffB'"), [])
        self.assertSucc('vtools compare --expression "BdiffA=ns-d_plekhn1"')
        self.assertEqual(output2list("vtools execute 'select * from BdiffA'"), ['714'])
        # use both options in one command should be allowed
        self.assertSucc('vtools compare --expression "A_OR_B=d_plekhn1|plekhn1"')
        #
        # handling of non-alphanum names
        self.assertSucc('vtools compare --expression "\'KK@\'=d_plekhn1 |ns_damaging"')
        self.assertTrue('KK@' in outputOfCmd('vtools show tables'))
        runCmd('vtools admin --rename_table d_plekhn1 "d@p"')
        self.assertSucc('vtools compare  --expression "\'K&K@\'=\'d@p\' & ns_damaging"')
        self.assertTrue('K&K@' in outputOfCmd('vtools show tables'))

    def testSimpleProj(self):
        'Test compare by site'
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/compare.vcf --build hg18')
        runCmd('vtools execute "DELETE FROM genotype.genotype_1 WHERE GT=0"')
        runCmd('vtools execute "DELETE FROM genotype.genotype_2 WHERE GT=0"')
        runCmd('vtools select variant "variant_id in (1, 2)" -t T1')
        runCmd('vtools select variant "variant_id in (1, 3, 4)" -t T2')
        # variant
        self.assertOutput('vtools compare T1 T2 ', '1\t2\t1\t4\n')
        self.assertOutput('vtools compare T1 T2 --difference', '1\n')
        self.assertOutput('vtools compare T2 T1 --difference', '2\n')
        self.assertOutput('vtools compare T1 T2 --intersection', '1\n')
        self.assertOutput('vtools compare T1 T2 --union', '4\n')
        # sample
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 ', '1\t1\t1\t3\n')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --difference ', '1\n')
        self.assertOutput('vtools compare variant --samples SAMP2 SAMP1 --difference ', '1\n')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --intersection ', '1\n')
        self.assertOutput('vtools compare variant --samples SAMP1 SAMP2 --union', '3\n')
        # site
        self.assertOutput('vtools compare T1 T2 --mode site ', '0\t1\t3\t4\n')
        self.assertOutput('vtools compare T1 T2 --mode site --difference ', '0\n')
        self.assertOutput('vtools compare T2 T1 --mode site --difference ', '1\n')
        self.assertOutput('vtools compare T1 T2 --mode site --intersection ', '3\n')
        self.assertOutput('vtools compare T1 T2 --mode site --union ', '4\n')



if __name__ == '__main__':
    unittest.main()
