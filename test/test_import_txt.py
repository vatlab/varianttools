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
from testUtils import ProcessTestCase, runCmd, numOfVariant

class TestImportTXT(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testImportTXT(self):
        'Test command import_txt'
        self.assertFail('vtools import_txt')
        self.assertFail('vtools import_txt non_existing.txt')
        # help information
        self.assertSucc('vtools import_txt -h')
        # no format information, fail
        self.assertFail('vtools import_txt input.tsv')
        # no build information, fail
        self.assertFail('vtools import_txt --format ../input_fmt/ANNOVAR input.tsv')

    def testANNOVAR(self):
        'Testing the annovar input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR ANNOVAR.txt')
        # one of the variant cannot be imported.
        self.assertEqual(numOfVariant(), 11)

    def testIllumina_SNP(self):
        'Testing the illumina SNP input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/Illumina_SNP Illumina_SNP.txt')
        self.assertEqual(numOfVariant(), 12)

    def testIllumina_INDEL(self):
        'Testing the illumina INDEL input format'
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/Illumina_INDEL Illumina_INDEL.txt')
        self.assertEqual(numOfVariant(), 12)

if __name__ == '__main__':
    unittest.main()
