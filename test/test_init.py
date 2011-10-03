#!/usr/bin/env python
#
# $File: test_init.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
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
from testUtils import ProcessTestCase, runCmd, numOfVariant, numOfSample

class TestInit(ProcessTestCase):
    def testInit(self):
        'Test command vtools init'
        self.assertFail('vtools init')
        self.assertSucc('vtools init test')
        self.assertFail('vtools init test')
        self.assertSucc('vtools init test -f')
    
    def testParent(self):
        'Test command init --parent'
        try:
            os.mkdir('parent')
        except OSError:
            pass
        
        os.chdir('parent')
        runCmd('vtools init test')
        runCmd('vtools import ../vcf/CEU.vcf.gz --build hg18')
        runCmd('vtools import ../vcf/SAMP1.vcf')
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t ceu')
        os.chdir('..')
        self.assertSucc('vtools init test --parent parent ceu')
        self.assertEqual(numOfVariant(), 288)
        self.assertEqual(numOfSample(), 61)
        for fn in os.listdir('parent'):
            if fn != 'cache': os.remove('parent/{}'.format(fn))
            else:
                os.remove('parent/cache/vcf.fmt')
                os.rmdir('parent/cache')
        os.rmdir('parent')
        


if __name__ == '__main__':
    unittest.main()
