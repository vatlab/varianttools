#!/usr/bin/env python
#
# $File: test_remove.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
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
from testUtils import ProcessTestCase, runCmd, initTest, outputOfCmd

class TestRemove(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(6)
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')
        runCmd('vtools select variant --samples "aff=\'1\'" -t unaffected')
        runCmd('vtools sample_stat CEU --samples "filename like \'%CEU%\' and aff=\'2\'" --num CEU_cases_num')
    def testRemove(self):
        'Test command vtools remove'
        self.assertFail('vtools remove')
        self.assertSucc('vtools remove -h')
        self.assertFail('vtools remove tables')
        # Removing table unaffected
        self.assertSucc('vtools remove tables unaffected')
        self.assertFail('vtools show table unaffected')
        # Removing field CEU_num from variant table CEU
        count1 = len(outputOfCmd('vtools show fields').split('\n'))
        self.assertSucc('vtools remove fields CEU_cases_num')
        count2 = len(outputOfCmd('vtools show fields').split('\n'))
        self.assertEqual(count1-count2, 1)
        self.assertSucc('vtools remove annotations testNSFP')
        self.assertFail('vtools show annotation testNSFP')
        self.assertSucc('vtools remove project')

if __name__ == '__main__':
    unittest.main()
