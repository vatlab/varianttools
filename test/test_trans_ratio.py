#!/usr/bin/env python
#
# $File: test_trans_ratio.py $
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

class TestTransRatio(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools update variant --from_stat "num=#(alt)" "depth=sum(DP_geno)"')

    def testTransRatio(self):
        'Test command vtools_report trans_ratio'
        self.assertFail('vtools_report trans_ratio')
        self.assertSucc('vtools_report trans_ratio -h')
        # no such column: non_existing_column
        self.assertFail('vtools_report trans_ratio non_existing -n non_existing_column')
        self.assertSucc('vtools_report trans_ratio variant -n num')
        # ValueError: invalid literal for int() with base 10: '6292.38333333\n'
        #self.assertFail('vtools_report trans_ratio variant -n depth')
        self.assertSucc('vtools_report trans_ratio variant -n num --group_by num')

if __name__ == '__main__':
    unittest.main()
