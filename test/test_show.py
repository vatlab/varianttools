#!/usr/bin/env python
#
# $File: test_show.py $
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

class TestShow(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
        self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools use ann/testNSFP.ann')
        
    def testShow(self):
        'Test command vtools show'
        self.assertSucc('vtools show -h')
        #Please specify a table to display
        self.assertFail('vtools show table')
        self.assertSucc('vtools show tables')
        # show project
        self.assertSucc('vtools show')
        # show project
        self.assertSucc('vtools show project')
        self.assertSucc('vtools show table variant')
        self.assertSucc('vtools show samples')
        self.assertFail('vtools show field')
        self.assertSucc('vtools show fields')
        self.assertFail('vtools show format')
        self.assertSucc('vtools show formats')
        self.assertFail('vtools show formats ANNOVAR')
        self.assertSucc('vtools show format ANNOVAR')
        self.assertSucc('vtools show genotypes')
        self.assertSucc('vtools show annotations')
        self.assertSucc('vtools show annotations testNSFP')
        self.assertSucc('vtools show annotation testNSFP')
        self.assertFail('vtools show test')
        self.assertSucc('vtools show tests')
        self.assertSucc('vtools show test CollapseBT')
        self.assertFail('vtools show tests CollapseBT')
        # show runtime options
        self.assertSucc('vtools show runtime_options')


if __name__ == '__main__':
    unittest.main()
