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
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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
from testUtils import ProcessTestCase, runCmd, initTest, output2list

class TestCompare(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(10)
    def removeProj(self):
        runCmd('vtools remove project')
        
    def testCompare(self):
        'Test command vtools compare'
        self.assertFail('vtools compare')
        self.assertSucc('vtools compare -h')
        # WARNING: No action parameter is specified. Nothing to do.
        self.assertFail('vtools compare plekhn1 d_plekhn1')
        self.assertSucc('vtools compare plekhn1 d_plekhn1 -c')
        self.assertEqual(output2list('vtools compare plekhn1 d_plekhn1 -c'), ['0\t0\t6\t6'])
        # error: argument --A_and_B: expected one argument
        self.assertFail('vtools compare d_plekhn1 ns_damaging --A_and_B')
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --A_and_B common')
        # WARNING: Existing table common is renamed to common_Aug09_170022.
        self.assertSucc('vtools compare d_plekhn1 ns_damaging --A_and_B common')
        self.assertEqual(output2list('vtools show table common')[1:], ['880', '881', '882', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --A_or_B AorB')
        self.assertEqual(output2list('vtools show table AorB')[1:], ['714', '869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --A_and_B AandB')
        self.assertEqual(output2list('vtools show table AandB')[1:], ['869', '880', '881', '882', '889', '893'])
        self.assertSucc('vtools compare d_plekhn1 ns --A_diff_B AdiffB')
        self.assertEqual(output2list('vtools show table AdiffB')[1:], [])
        self.assertSucc('vtools compare d_plekhn1 ns --B_diff_A BdiffA')
        self.assertEqual(output2list('vtools show table BdiffA')[1:], ['714'])
        # use both options in one command should be allowed
        self.assertSucc('vtools compare d_plekhn1 plekhn1 -c --A_or_B A_OR_B')

if __name__ == '__main__':
    unittest.main()
