#!/usr/bin/env python
#
# $File: test_pipeline.py $
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
import shutil
from testUtils import ProcessTestCase


class TestPipeline(ProcessTestCase):
    def testCheckVariantToolsVersion(self):
        'Test functor CheckVariantToolsVersion'
        self.assertSucc('vtools execute test_pipeline.pipeline checkvtools')
        self.assertFail(
            'vtools execute test_pipeline.pipeline checkvtools_fail')

    def testCheckCommand(self):
        'Test functor CheckCommand'
        self.assertSucc('vtools execute test_pipeline.pipeline checkcommand')
        self.assertFail(
            'vtools execute test_pipeline.pipeline checkcommand_fail')

    def testRunCommand(self):
        self.assertSucc('vtools execute test_pipeline.pipeline runcommand')
        self.assertFail(
            'vtools execute test_pipeline.pipeline runcommand_fail')

    # def testCreatePopulation(self):
    #     self.assertSucc('vtools execute test_pipeline.pipeline createpop')

    # def testEvolvePopulation(self):
    #     self.assertFail('vtools execute test_pipeline.pipeline evolvepop_fail')
    #     self.assertSucc('vtools execute test_pipeline.pipeline evolvepop1')
    #     self.assertSucc('vtools simulate test_pipeline.pipeline evolvepop1')
    #     self.assertSucc('vtools simulate test_pipeline.pipeline evolvepop3')
    #     self.assertSucc('vtools simulate test_pipeline.pipeline evolvepop4')

    # def testDrawCaseCtrlSample(self):
    #     self.assertSucc('vtools simulate test_pipeline.pipeline casectrl1')
    #     self.assertTrue(os.path.isfile('cache/casectrl_sample_12345.pop'))
    #     self.assertSucc('vtools simulate test_pipeline.pipeline casectrl2')
    #     self.assertTrue(os.path.isfile('cache/casectrl_sample_12345_21.pop'))
    #     self.assertTrue(os.path.isfile('cache/casectrl_sample_12345_22.pop'))

    # def testDrawRandomSample(self):
    #     self.assertSucc('vtools simulate test_pipeline.pipeline random1')
    #     self.assertTrue(os.path.isfile('cache/random_sample_12345.pop'))
    #     self.assertSucc('vtools simulate test_pipeline.pipeline random2')
    #     self.assertTrue(os.path.isfile('cache/random_sample_12345_21.pop'))
    #     self.assertTrue(os.path.isfile('cache/random_sample_12345_22.pop'))


if __name__ == '__main__':
    unittest.main()
