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
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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
from testUtils import ProcessTestCase, runCmd, output2list

class TestShow(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')
        runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
        runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        runCmd('vtools use ann/testNSFP.ann')
        
    def removeProj(self):
        runCmd('vtools remove project')
    def testShow(self):
        'Test command vtools show'
        self.assertSucc('vtools show -h')
        # Please specify a table to display
        self.assertFail('vtools show table')
        self.assertSucc('vtools show tables')
        # show project
        self.assertSucc('vtools show')
        # show project
        self.assertSucc('vtools show project')
        self.assertEqual(output2list('vtools show table variant'), ['variant_id, bin, chr, pos, ref, alt', \
                                                         '1, 585, 1, 75927, G, C', \
                                                         '2, 585, 1, 76193, A, G', \
                                                         '3, 585, 1, 77052, G, A', \
                                                         '4, 585, 1, 78178, G, A', \
                                                         '5, 585, 1, 78200, G, A', \
                                                         '6, 585, 1, 81398, G, T', \
                                                         '7, 585, 1, 98172, T, C', \
                                                         '8, 586, 1, 223335, C, G', \
                                                         '9, 586, 1, 224622, A, T', \
                                                         '10, 586, 1, 225791, G, A'])
        self.assertSucc('vtools show table variant -l 20')
        self.assertSucc('vtools show table variant -l -1')
        self.assertSucc('vtools show samples')
        self.assertSucc('vtools show table testNSFP')
        self.assertSucc('vtools show fields')
        self.assertSucc('vtools show formats')
        self.assertSucc('vtools show genotypes')
        self.assertSucc('vtools show annotations')

if __name__ == '__main__':
    unittest.main()
