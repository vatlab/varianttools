#!/usr/bin/env python
#
# $File: test_import_txt.py $
# $LastChangedDate: 2011-July-29 5:00:00 -0500 (Friday, 29 July 2011) $
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
from testUtils import ProcessTestCase, runCmd, numOfSample

class TestImportTXT(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd(['vtools', 'init', 'test'])

    def testImportVCF(self):
        'Test command import_txt'
        self.assertFail(['vtools', 'import_txt'])
        self.assertFail(['vtools', 'import_txt', 'non_existing.tsv'])
        # no build information, fail
        self.assertFail(['vtools', 'import_txt', 'input.tsv'])
        # specify build information
        self.assertSucc(['vtools', 'import_txt', 'input.tsv', '--build', 'hg18'])
        #file will be ignored if re-imported
        self.assertSucc(['vtools', 'import_txt', 'input.tsv'])


if __name__ == '__main__':
    unittest.main()
