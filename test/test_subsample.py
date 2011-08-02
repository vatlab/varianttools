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

class TestSubsample(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf --build hg18')
        runCmd('vtools import_txt input.tsv -c 1 2 4 5')
        runCmd('vtools import_phenotype phenotype.txt')
    def removeProj(self):
        runCmd('vtools remove project')
    def testSubsample(self):
        'Test command vtools subsample'
        self.assertFail('vtools subsample')
        self.assertSucc('vtools subsample -h')
        # Cannot overwrite master variant table
        self.assertFail('vtools subsample aff=1 -t variant')
        self.assertSucc('vtools subsample aff=1 -t unaffected1')
        self.assertSucc('vtools subsample "aff=\'1\'" -t unaffected2')
        # Failed to retrieve samples by condition "sex=M"
        self.assertFail('vtools subsample sex=\'M\' -t sexm')
        # Failed to retrieve samples by condition "sex=M"
        self.assertFail('vtools subsample \'sex=M\' -t sexm')
        self.assertSucc('vtools subsample sex=\\\'M\\\' -t sexm')
        self.assertSucc('vtools subsample "sex=\'M\'" -t sexm')
        self.assertSucc('vtools subsample "filename like \'CEU%\'" -t CEU')
        self.assertSucc('vtools subsample "BMI<18.5" -t Underweight')

if __name__ == '__main__':
    unittest.main()