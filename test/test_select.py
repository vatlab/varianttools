#!/usr/bin/env python
#
# $File: test_select.py $
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

class TestSelect(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf --build hg18')
        runCmd('vtools import_vcf CEU.vcf --build hg19')
    def removeProj(self):
        runCmd('vtools remove project')
    def testSelect(self):
        'Test command select'
        self.assertFail('vtools select')
        self.assertFail('vtools select non_existing_variant')
        # help information
        self.assertSucc('vtools select -h')
        # could return tuples by selection table under certain condition 
        self.assertSucc('vtools select variant "genename = "PLEKHN1""')
        # Limit selection from samples      
        self.assertSucc('vtools select variant "genename = "PLEKHN1"" --samples "aff=1"')
        # test with rename table
        self.assertSucc('vtools select variant "genename = "PLEKHN1"" -t plekhn1')
        # count the number of variant
        self.assertSucc('vtools select -c variant')
        # output
        self.assertSucc('vtools output plekhn1 chr pos ref alt genename sift_score --build=hg18')

if __name__ == '__main__':
    unittest.main()
