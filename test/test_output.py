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

class TestOutput(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf.gz --build hg18')
        runCmd('vtools import_txt input.tsv -c 1 2 4 5')
        runCmd('vtools import_phenotype phenotype.txt')
        runCmd('vtools import_vcf SAMP1.vcf')
        runCmd('vtools subsample "filename like \'CEU%\'" -t CEU')
        runCmd('vtools subsample "BMI<18.5" -t Underweight')
    def testOutput(self):
        'Test command vtools output'
        self.assertFail('vtools output')
        self.assertSucc('vtools output -h')
        # too few arguments
        self.assertFail('vtools output variant')
        self.assertFail('vtools output variant non_existing_field')
        self.assertSucc('vtools output variant freq chr pos ref alt')
        self.assertSucc('vtools output variant freq num depth other -l 10')
        # too few arguments
        self.assertFail('vtools output variant -s variant.txt')
        self.assertSucc('vtools output variant variant_id freq num depth other -s variant.txt')
        self.assertSucc('vtools output Underweight chr pos ref alt -l -1')
        self.assertFail('vtools output Underweight and variant freq -l 10')
        self.assertFail('vtools output variant sum(num)')
        self.assertSucc('vtools output variant "sum(num)"')
        self.assertFail('vtools output variant count(1)')
        self.assertSucc('vtools output variant "count(1)"')
        self.assertSucc('vtools output variant "count(variant_id)" -v 0')
        self.assertSucc('vtools output variant "avg(depth)''')
        self.assertSucc('vtools output CEU "avg(freq)"')
        self.assertFail('vtools output CEU "avg(bin)"')

if __name__ == '__main__':
    unittest.main()