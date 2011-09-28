#!/usr/bin/env python
#
# $File: test_show.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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

class TestShow(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(6)
    def removeProj(self):
        runCmd('vtools remove project')
    def testShow(self):
        'Test command vtools show'
        self.assertSucc('vtools show -h')
        # Please specify a table to display
        self.assertFail('vtools show table')
        self.assertSucc('vtools show tables')
        # show project
        self.assertSucc('vtools show')
        # show project
        self.assertSucc('vtools show project')
        self.assertSucc('vtools show table variant')
        self.assertSucc('vtools show table variant -l 20')
        self.assertSucc('vtools show table variant -l -1')
        self.assertSucc('vtools show samples')
        self.assertSucc('vtools show table testNSFP')
        self.assertSucc('vtools show fields')
        self.assertSucc('vtools show formats')
        self.assertSucc('vtools show annotations')

if __name__ == '__main__':
    unittest.main()
