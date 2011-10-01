# !/usr/bin/env python
#
# $File: test_phenotype --from_file.py $
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
        runCmd('vtools import vcf/SAMP2.vcf')
    def removeProj(self):
        runCmd('vtools remove project')
    def testImportPhenotype(self):
        'Test command phenotype --from_file'
        # too few arguments
        self.assertFail('vtools phenotype --from_file')
        self.assertSucc('vtools phenotype --from_file -h')
        # opening project project_name. Importing phenotypes into table sample.
        self.assertSucc('vtools phenotype --from_file phenotype/phenotype.txt')
        out1 = output2list('vtools show samples')
        with open('phenotype/phenotype.txt') as inputfile:
            out2 = [x[:-1] for x in inputfile]
        self.assertEqual(out1, out2)
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype1.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype2.txt')
        self.assertFail('vtools phenotype --from_file phenotype/badphenotype3.txt')

if __name__ == '__main__':
    unittest.main()
