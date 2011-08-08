#!/usr/bin/env python
#
# $File: test_import_vcf.py $
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
from testUtils import ProcessTestCase, runCmd

class TestShow(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf.gz --build hg18')
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
        # Save the first 20 lines to variant.txt
        self.assertSucc('vtools show table variant -l 20 > variant.txt')
        # Save output to variant.txt
        self.assertSucc('vtools show table variant -l -1 > variant.txt')
        self.assertSucc('vtools show sample')
        # Save output to phenotype.txt
        self.assertSucc('vtools show sample > phenotype.txt')
        runCmd('vtools use testNSFP.ann')
        self.assertSucc('vtools show table testNSFP')
        self.assertSucc('vtools show fields')

if __name__ == '__main__':
    unittest.main()