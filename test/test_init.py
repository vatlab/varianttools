#!/usr/bin/env python
#
# $File: test_init $
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
import unittest
import subprocess


class ProcTest(unittest.TestCase):
    def assertOutput(self, cmd, output):
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        self.assertEqual(
            subprocess.check_output(cmd, env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}),
            output)

    def assertSucc(self, cmd):
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        self.assertEqual(subprocess.check_call(cmd, env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}),
            0)

    def assertFail(self, cmd):
        try:
            subprocess.check_call(cmd, env={'PATH': os.pathsep.join(['..', os.environ['PATH']])})
        except subprocess.CalledProcessError:
            return


class TestInit(ProcTest):
    def testInit(self):
        self.assertFail(['vtools', 'init'])
        self.assertSucc(['vtools', 'init', 'test'])
        self.assertFail(['vtools', 'init', 'test1'])


if __name__ == '__main__':
    unittest.main()
