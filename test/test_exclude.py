#!/usr/bin/env python
#
# $File: test_exclude.py $
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


class TestExclude(ProcessTestCase):
    def setUp(self):
        'Create a project'
        ProcessTestCase.setUp(self)
        # self.runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        self.runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
        self.runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
        self.runCmd('vtools use ann/testNSFP.ann')
        self.runCmd('vtools select variant \'testNSFP.chr is not null\' -t ns')
        
    def testExclude(self):
        'Test command vtools exclude'
        self.assertFail('vtools exclude')
        self.assertSucc('vtools exclude -h')
        # Variant table non_existing_table does not exist.
        self.assertFail('vtools exclude non_existing_table "aff=\'1\'" -t aff')
        # Cannot overwrite master variant table. Please choose another name for the variant table.
        self.assertFail('vtools exclude ns "sift_score<=0.94" -t variant')
        self.assertSucc('vtools exclude ns "sift_score <= 0.94" -t ns_non_damaging')
        #should have 5 variants
        #select "sift_score >= 0.94" will result in 6 variants
        self.assertOutput('vtools select ns "sift_score > 0.94" -c', 'output/exclude_sift.txt')
        self.assertOutput('vtools exclude ns "sift_score <= 0.94" -c', 'output/exclude_sift.txt')
        self.assertOutput('vtools exclude ns "variant_id=604" -c', '1445')


    # def testExcludeAnno(self):
    #    self.runCmd('vtools liftover hg19')
    #    # FIXME: this test is bad because it requies downloading a large database
    #    self.runCmd('vtools use dbSNP')
    #    # this will exclude ones called unknown, but will keep the ones that is NULL
    #    # and remove variants with more than one entries
    #    #
    #    #
    #    self.assertSucc(''' vtools exclude variant "dbSNP.valid='unknown'" -t known ''')
    #    self.assertOutput('vtools output known variant_id chr pos ref alt ', 'output/exclude_anno.txt')
    #    self.assertSucc('vtools update variant --set dbSNP_valid=dbSNP.valid ')
    #    self.assertOutput('''vtools select variant 'dbSNP.valid <> "unknown" OR dbSNP.valid is NULL' --output variant_id chr pos ref alt''',
    #        'output/exclude_select_anno.txt')


if __name__ == '__main__':
    unittest.main()
