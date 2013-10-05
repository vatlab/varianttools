#!/usr/bin/env python
#
# $File: test_liftover.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org)
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
from testUtils import ProcessTestCase, runCmd, outputOfCmd, initTest

class TestLiftover(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(2)

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
        # from hg18 to hg19
        self.assertSucc('vtools liftover hg19')
        out1 = outputOfCmd('vtools output variant bin chr pos alt_bin alt_chr alt_pos -d"\t"')
        out1 = '\n'.join([x for x in out1.split('\n') if '.' not in x])
        #
        # We write in hg19 to a datafile, create a new project, import the
        # data and liftover to hg18, we then compare if coordinates in these
        # projects are the same.
        # 
        data = outputOfCmd('vtools output variant chr pos ref alt --build hg19 -d"\t"') 
        with open('temp_input.txt', 'w') as output:
            output.write(data)
        runCmd('vtools init test -f')
        self.assertSucc('vtools import --build hg19 --format ../format/basic temp_input.txt')
        self.assertSucc('vtools liftover hg18')
        out2 = outputOfCmd('vtools output variant alt_bin alt_chr alt_pos bin chr pos -d"\t"')
        out2 = '\n'.join([x for x in out2.split('\n') if '.' not in x])
        self.assertEqual(out1, out2)
        os.remove('temp_input.txt')

    def testLiftoverFlip(self):
        self.assertFail('vtools liftover hg18')
        self.assertFail('vtools liftover hg18 --flip')
        self.assertSucc('vtools liftover hg19')
        self.assertFail('vtools liftover hg18')
        self.assertSucc('vtools liftover hg19 --flip')
        self.assertSucc('vtools liftover hg18 --flip')   
        var_tab18  = outputOfCmd('vtools output variant bin chr pos')
        var_tab18 = '\n'.join([x for x in var_tab18.split('\n') if 'NA' not in x])
        runCmd('vtools liftover hg19 --flip')   
        var_tab19  = outputOfCmd('vtools output variant alt_bin alt_chr alt_pos')
        var_tab19 = '\n'.join([x for x in var_tab19.split('\n') if 'NA' not in x])
        self.assertEqual(var_tab18, var_tab19)

 
if __name__ == '__main__':
    unittest.main()
