#/usr/bin/env python

# $file: test_admin.py $
# $lastChangedDate: 2012-05-14 $

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
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample,\
    getGenotypes, getSamplenames, output2list, getGenotypeInfo, outputOfCmd

class TestAdmin(ProcessTestCase):

    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testMergeSamples(self):
        'Test command vtools admin --merge_samples'
        #
        # FIXME: merge that involves all genotype tables will use a 
        # different algorithm that needs to be tested
        self.assertFail('vtools admin')
        self.assertSucc('vtools admin -h')
        self.assertSucc('vtools import vcf/CEU.vcf.gz --build hg18')
        self.assertEqual(numOfSample(), 60)
        self.assertEqual(numOfVariant(),288)
        self.assertSucc('vtools admin --rename_samples \'sample_name like "%NA069%"\' NA06900')
        #could not merge them together if they are from the same table.
        self.assertFail('vtools admin --merge_samples')
        # Test command vtools admin --merge_samples'
        runCmd('vtools init test -f')
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18') 
        runCmd('vtools import vcf/SAMP1.vcf  --build hg18')
        self.assertEqual(numOfVariant(),577)
        self.assertEqual(numOfSample(), 61)
        runCmd('vtools admin --rename_samples \'filename like "%SAMP1%"\' NA06985')
        self.assertEqual(numOfSample(), 61)
        runCmd('vtools admin --merge_samples')         
        self.assertEqual(numOfVariant(),577)
        self.assertEqual(numOfSample(), 60)
        # Test merge samples with overlapping variants
        runCmd('vtools init test -f')
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg18')
        self.assertSucc('vtools import vcf/SAMP1.vcf  --build hg18')
        self.assertSucc('vtools admin --rename_samples \'filename like "%2%"\' SAMP1')
        self.assertFail('vtools admin --merge_samples')
        #the reason is that the two samepls have some identical variants. If you want to merge them, the samples should have different unique variant information.

    def testRenameSamples(self):
        'Test command vtools admin --rename_samples'
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18') 
        self.assertFail('vtools admin --rename_samples 1')
        # all samples are assigned name NA
        self.assertFail('vtools admin --rename_samples "sample_name like \'NA1\'" NA')
        self.assertFail('vtools admin --rename_samples 1 NA NNA')

    def testRenameTable(self):
        'test rename tables'
        self.assertSucc('vtools select variant -t "%ad name"')
        self.assertTrue('%ad name' in outputOfCmd('vtools show tables'))
        self.assertTrue('vtools show table "%ad name"')
        self.assertFail('vtools admin --rename_table variant not_allowed')
        self.assertSucc('vtools admin --rename_table "%ad name" UNIQUE')
        self.assertTrue('vtools show table UNIQUE')
        self.assertTrue('UNIQUE' in outputOfCmd('vtools show tables'))
        self.assertSucc('vtools admin --rename_table UNIQUE "%ad newname"')
        self.assertTrue('vtools show table "%ad newname"')
        self.assertTrue('%ad newname' in outputOfCmd('vtools show tables'))

    def testDescribeTable(self):
        'test describe tables'
        self.assertSucc('vtools select variant -t "%%" "DESD"')
        self.assertTrue('%%' in outputOfCmd('vtools show tables'))
        self.assertTrue('DESD' in outputOfCmd('vtools show tables'))
        self.assertSucc('vtools admin --describe_table %% "NN NN"')
        self.assertFalse('DESD' in outputOfCmd('vtools show tables'))
        self.assertTrue('NN NN' in outputOfCmd('vtools show tables'))

    def testSaveLoadSnapshot(self):
        'test save/load snapshot'
        self.assertFail('vtools admin --save_snapshot a')
        self.assertSucc('vtools admin --save_snapshot a "some comment"')
        self.assertSucc('vtools admin --load_snapshot a')
        self.assertSucc('vtools admin --save_snapshot a.tar.gz "some comment"')
        self.assertSucc('vtools admin --load_snapshot a.tar.gz')

    def testRuntimeOption(self):
        'test set runtime options'
        # FIXME: test for valid options
        # FIXME: test for value of options -- for whatever value, vtools should be able to load the project
        # FIXME: test for non-exist temp_dir
        # FIXME: test for sqlite_pragma
        pass

    def testResetRuntimeOption(self):
        'test reset runtime options'
        # FIXME: fail no OPT
        # FIXME: fail invalid OPT
        # FIXME: set and rest sqlite_pragma, use vtools show to check
        pass
        
if __name__ == '__main__':
    unittest.main()
