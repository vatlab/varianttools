#!/usr/bin/env python
#
# $File: test_execute.py $
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
from testUtils import ProcessTestCase, runCmd, initTest

class TestExecute(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(5)
    def removeProj(self):
        runCmd('vtools remove project')
    def testExecute(self):
        'Test command vtools execute'
        self.assertFail('vtools execute')
        self.assertSucc('vtools execute -h')
        self.assertFail('vtools execute select non_existing_field from sample -v2')
        self.assertSucc('vtools execute select sample_name from sample -v2')
        self.assertFail('vtools execute select * from sample -v0')
        self.assertSucc('vtools execute \'select * from sample\'')
        self.assertSucc('vtools execute \'delete FROM variant where chr=1\' -v2')
        self.assertSucc('vtools execute \'delete FROM variant where ref="T"\' -v2')
        self.assertSucc('vtools execute \'delete FROM sample where aff=2\' -v2')

if __name__ == '__main__':
    unittest.main()
