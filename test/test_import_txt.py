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
from testUtils import ProcessTestCase, runCmd, numOfVariant

class TestImportTXT(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testImportTXT(self):
        'Test command import_txt'
        self.assertFail('vtools import_txt')
        self.assertFail('vtools import_txt non_existing.txt')
        # help information
        self.assertSucc('vtools import_txt -h')
        # no format information, fail
        self.assertFail('vtools import_txt input.tsv')
        # no build information, fail
        self.assertFail('vtools import_txt --format ../input_fmt/ANNOVAR input.tsv')

    def testANNOVAR(self):
        'Testing the annovar input format'
        input = open('ANNOVAR.txt', 'w')
        input.write('''\
1 161003087 161003087 C T comments: rs1000050, a SNP in Illumina SNP arrays
1 84647761 84647761 C T comments: rs6576700 or SNP_A-1780419, a SNP in Affymetrix SNP arrays
1 13133880 13133881 TC - comments: rs59770105, a 2-bp deletion
1 11326183 11326183 - AT comments: rs35561142, a 2-bp insertion
1 105293754 105293754 A ATAAA comments: rs10552169, a block substitution
1 67478546 67478546 G A comments: rs11209026 (R381Q), a SNP in IL23R associated with Crohn's disease
2 233848107 233848107 T C comments: rs2241880 (T300A), a SNP in the ATG16L1 associated with Crohn's disease
16 49303427 49303427 C T comments: rs2066844 (R702W), a non-synonymous SNP in NOD2
16 49314041 49314041 G C comments: rs2066845 (G908R), a non-synonymous SNP in NOD2
16 49321279 49321279 - C comments: rs2066847 (c.3016_3017insC), a frameshift SNP in NOD2
13 19661686 19661686 G - comments: rs1801002 (del35G), a frameshift mutation in GJB2, associated with hearing loss
13 19695176 20003944 0 - comments: a 342kb deletion encompassing GJB6, associated with hearing loss
''')
        input.close()
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR ANNOVAR.txt')
        self.assertEqual(numOfVariant(), 12)


if __name__ == '__main__':
    unittest.main()
