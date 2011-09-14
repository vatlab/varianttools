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
from testUtils import ProcessTestCase, runCmd, initTest

class TestSampleStat(ProcessTestCase):
    def setUp(self):
        'Create a project'
        initTest(5)
        runCmd('vtools select variant --samples "filename like \'CEU%\'" -t CEU')
    def removeProj(self):
        runCmd('vtools remove project')
    def testSampleStat(self):
        'Test command vtools sample_stat'
        self.assertFail('vtools sample_stat')
        self.assertSucc('vtools sample_stat -h')
        self.assertFail('vtools sample_stat --num --freq --hom --het --other --depth')
        self.assertFail('vtools sample_stat --num select')
        self.assertSucc('vtools sample_stat variant --num num --freq freq --hom hom --het het --other other --depth depth')
        self.assertSucc('vtools sample_stat CEU -s "filename like \'CEU%\'" --freq CEU_freq')
        self.assertSucc('vtools sample_stat CEU --samples "filename like \'CEU%\'" --num CEU_num --hom CEU_hom --het CEU_het --other CEU_other')
        self.assertSucc('vtools sample_stat CEU --samples "filename like \'CEU%\' and aff=\'2\'" --freq CEU_cases_freq --het CEU_cases_het')
        self.assertSucc('vtools sample_stat CEU -s "filename like \'CEU%\' and aff=\'1\'" --freq CEU_ctrls_freq --het CEU_ctrls_het')

if __name__ == '__main__':
    unittest.main()