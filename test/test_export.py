#/usr/bin/env python
#
# $File: test_export.py $
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
#from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample, outputOfCmd, getGenotypes, getSamplenames, output2list, getGenotypeInfo

class TestImport(ProcessTestCase):
    
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def tearDown(self):
        runCmd('vtools remove project')
        try:
            dir = os.getcwd()
            vcf = os.path.join(dir, 'my.vcf')
            os.remove(vcf)
        except Exception as e:
            pass

    def testExportVcfSnv(self):
        'Test command export in vcf format'
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        # test basic vcf output
        self.assertSucc('vtools export variant my.vcf')
        variants = output2list('vtools output variant chr pos ref alt')
        with open('my.vcf') as infile:
            exported = infile.readlines()
            exported = ['\t'.join([x.split('\t')[0], x.split('\t')[1], x.split('\t')[3], x.split('\t')[4]]) for x in exported]
            # handle multiple alternatives
            for idx, item in enumerate(exported):
                if ',' in item:
                    exported[idx] = '\t'.join(item.split(',')[:-1])
                    exported.insert(idx, exported[idx])
                    tmp = exported[idx].split('\t')[:-1]
                    tmp.extend(item.split(',')[-1])
                    exported[idx+1] = '\t'.join(tmp)
            self.assertEqual(variants, exported)
            # test output with samples
            # FIXME phase information is not properly handled
            self.assertSucc('vtools export variant my.vcf --samples 1 --format_string GT')
            self.assertSucc('vtools export variant my.vcf --samples 1 --format_string GT --phase_sep "|"')
            self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --var_info id filter info --geno_info DP_geno')
            self.assertSucc('vtools export variant my.vcf --id id --filter filter --info info')
            self.assertSucc('vtools export variant my.vcf --id id --filter filter --info info --geno_info DP_geno --samples 1 --format_string GT:DP')
            # export selected table with selected info field
            self.assertSucc('vtools update variant --from_file vcf/CEU.vcf.gz --var_info DP')
            self.assertSucc('vtools select variant "DP > 300" -t highDP')
            self.assertSucc('vtools export highDP my.vcf --id id --filter filter --info DP --geno_info DP_geno --samples 1 --format_string GT:DP')
            
    def testExportANNOVAR(self):
        'Test command export in annovar input format'
        pass
   
        
if __name__ == '__main__':
    unittest.main()
