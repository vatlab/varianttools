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
from testUtils import ProcessTestCase

class TestLiftover(ProcessTestCase):

    def testLiftover(self):
        'Test command vtools liftover'
        # too few arguments
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertFail('vtools liftover')
        self.assertSucc('vtools liftover -h')
        # too few arguments
        self.assertFail('vtools liftover -v 0')
        # non_existing_build
        self.assertFail('vtools liftover non_existing_build')
        # from hg18 to hg19
        self.assertSucc('vtools liftover hg19')
        self.assertOutput('vtools output variant bin chr pos alt_bin alt_chr alt_pos -d"\t"',
            'output/liftover_cmp.txt')
        
        #We write in hg19 to a datafile, create a new project, import the
        #data and liftover to hg18, we then compare if coordinates in these
        #projects are the same.
        
        # self.assertOutput('vtools output variant chr pos ref alt --build hg19 -d"\t"', 'output/liftover.txt')
        self.runCmd('vtools init test -f --store '+self.storeMode)
        self.assertSucc('vtools import --build hg19 --format ../resources/format/basic output/liftover.txt')
        self.assertSucc('vtools liftover hg18')
        self.assertOutput('vtools output variant alt_bin alt_chr alt_pos bin chr pos -d"\t"',
            'output/liftover_cmp.txt', lambda x: sorted([i for i in x if '.' not in i]))

    def testLiftoverFlip(self):
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertFail('vtools liftover hg18')
        self.assertFail('vtools liftover hg18 --flip')
        self.assertSucc('vtools liftover hg19')
        self.assertFail('vtools liftover hg18')
        self.assertSucc('vtools liftover hg19 --flip')
        self.assertSucc('vtools liftover hg18 --flip')   
        var_tab18  = self.runCmd('vtools output variant bin chr pos')
        var_tab18 = '\n'.join([x for x in var_tab18.split('\n') if 'NA' not in x])
        self.runCmd('vtools liftover hg19 --flip')   
        var_tab19  = self.runCmd('vtools output variant alt_bin alt_chr alt_pos')
        var_tab19 = '\n'.join([x for x in var_tab19.split('\n') if 'NA' not in x])
        self.assertEqual(var_tab18, var_tab19)

 
if __name__ == '__main__':
    unittest.main()
