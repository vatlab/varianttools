#!/usr/bin/env python
#
# $File: test_init.py $
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

class TestLiftover(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import_vcf CEU.vcf.gz --build hg18')
    def removeProj(self):
        runCmd('vtools remove project')
    def testLiftover(self):
        'Test command vtools liftover'
        # too few arguments
        self.assertFail('vtools liftover')
        self.assertSucc('vtools liftover -h')
        # too few arguments
        self.assertFail('vtools liftover -v 0')
        # non_existing_build
        self.assertFail('vtools liftover non_existing_build')
        self.assertSucc('vtools liftover hg19')

    #
    # TODO: choose a marker with known hg18 and hg19 coordinates
    # import and liftover, and check the results
    #
if __name__ == '__main__':
    unittest.main()
