#!/usr/bin/env python

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
import unittest
from testUtils import ProcessTestCase

class TestAdmin(ProcessTestCase):
    def testAdminCommand(self):
        'Test command line options of vtools admin'
        self.assertSucc('vtools admin')
        self.assertSucc('vtools admin -h')

    def testMergeSamples(self):
        # Test command vtools admin --merge_samples'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18') 
        # this should actually be hg19, but we use it anyway for testing.
        self.runCmd('vtools import vcf/SAMP1.vcf  --build hg18')
        self.assertProj(numOfVariants=577, numOfSamples= 61)
        self.runCmd('vtools admin --rename_samples \'filename like "%SAMP1%"\' NA06985')
        self.assertProj(numOfSamples= 61)
        if self.storeMode=="sqlite":
            self.runCmd('vtools admin --merge_samples')         
            self.assertProj(numOfVariants=577, numOfSamples= 60)

    @unittest.skipUnless(os.getenv("STOREMODE")=="sqlite","HDF5 version is not implemented for this test")
    def testMergeWithOverlappingSamples(self):
        # Test merge samples with overlapping variants
        self.assertSucc('vtools import vcf/SAMP2.vcf --build hg19')
        self.assertSucc('vtools import vcf/SAMP1.vcf')
        self.assertSucc('vtools admin --rename_samples \'filename like "%2%"\' SAMP1')
        #the reason is that the two samepls have some identical variants. 
        # If you want to merge them, the samples should have different unique variant information.
        self.assertFail('vtools admin --merge_samples')

    def testRenameSamples(self):
        'Test command vtools admin --rename_samples'
        self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18') 
        self.assertFail('vtools admin --rename_samples 1')
        # all samples are assigned name NA
        self.assertSucc('''vtools admin --rename_samples "sample_name like 'NA1'" NA ''')
        # NA12716 is not renamed because of no match.
        self.assertProj(sampleNames=['NA12716'], partial=True)
        self.assertSucc('''vtools admin --rename_samples "sample_name like 'NA1%'" NA ''')
        # sample no longer exists
        self.assertProj(sampleNames=['NA12716'], partial=True, negate=True)
        self.assertSucc('vtools admin --rename_samples 1 NA NNA')

    def testRenameTable(self):
        'test rename tables'
        self.assertSucc('vtools select variant -t "%ad name"')
        self.assertProj(hasTable='%ad name')
        self.assertTrue('vtools show table "%ad name"')
        self.assertFail('vtools admin --rename_table variant not_allowed')
        self.assertSucc('vtools admin --rename_table "%ad name" UNIQUE')
        self.assertProj(hasTable='UNIQUE')
        self.assertSucc('vtools admin --rename_table UNIQUE "%ad newname"')
        self.assertProj(hasTable='%ad newname')

    def testDescribeTable(self):
        'test describe tables'
        self.runCmd('vtools import vcf/SAMP1.vcf --build hg19')
        self.assertSucc('vtools select variant -t "%%" "DESD"')
        self.assertProj(hasTable='%%', tableDesc={'%%': 'DESD'})
        self.assertSucc('vtools admin --describe_table %% "NN NN"')
        self.assertProj(hasTable='%%', tableDesc={'%%': 'NN NN'})

    def testSaveLoadSnapshot(self):
        'test save/load snapshot'
        self.assertFail('vtools admin --save_snapshot a')
        self.assertSucc('vtools admin --save_snapshot a "some comment"')
        self.assertSucc('vtools admin --load_snapshot a')
        self.assertSucc('vtools admin --save_snapshot a.tar.gz "some comment"')
        self.assertTrue(os.path.isfile('a.tar.gz'))
        self.assertSucc('vtools admin --load_snapshot a.tar.gz')
        os.remove('a.tar.gz')
        
if __name__ == '__main__':
    unittest.main()
