# !/usr/bin/env python
#
# $File: test_import_phenotype.py $
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
from testUtils import ProcessTestCase, runCmd, numOfSample, initTest, output2list

class TestImportPhenotype(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(3)
        runCmd('vtools import_variants vcf/SAMP2.vcf')
    def removeProj(self):
        runCmd('vtools remove project')
    def testImportPhenotype(self):
        'Test command import_phenotype'
        # too few arguments
        self.assertFail('vtools import_phenotype')
        self.assertSucc('vtools import_phenotype -h')
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools import_phenotype phenotype/phenotype.txt')
        out1 = output2list('vtools show samples')
        out2 = [x.split() for x in file('phenotype/phenotype.txt')]
        self.assertEqual(out1, out2)
        self.assertFail('vtools import_phenotype phenotype/badphenotype1.txt')
        self.assertFail('vtools import_phenotype phenotype/badphenotype2.txt')
        self.assertFail('vtools import_phenotype phenotype/badphenotype3.txt')

if __name__ == '__main__':
    unittest.main()
