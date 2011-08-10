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
        # no build information, fail
        self.assertFail('vtools import_txt input.tsv')
        # no columns, fail
        self.assertFail('vtools import_txt input.tsv --build hg18')
        # Four columns are required for each variant (chr, pos, ref, and alt)
        self.assertFail('vtools import_txt input.tsv --build hg18 -c 1 2 3 4 5')
        # import with four columns
        self.assertSucc('vtools import_txt input.tsv -c 1 2 4 5 --build hg18 --zero')
        ## re-import, fail
        self.assertFail('vtools import_txt input.tsv -c 1 2 4 5 --build hg18')
        # import different genome reference
        self.assertFail('vtools import_txt input.tsv -c 1 3 4 5 --build hg19')

if __name__ == '__main__':
    unittest.main()
