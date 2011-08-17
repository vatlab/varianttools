#!/usr/bin/env python
#
# $File: ProcessTestCase $
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
import shlex, subprocess

class ProcessTestCase(unittest.TestCase):
    'A subclass of unittest.TestCase to handle process output'
# used to be empty; force delete
    def setUp(self):
        'Clear any existing project'
        runCmd('vtools remove project')

    def tearDown(self):
        'Clear any existing project'
        runCmd('vtools remove project')
# compare if the command output is what we want
    def assertOutput(self, cmd, output):
        cmd = shlex.split(cmd)
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        self.assertEqual(
            subprocess.check_output(cmd, stderr=subprocess.PIPE,
                env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}),
            output)

    def assertSucc(self, cmd):
        cmd = shlex.split(cmd)
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        self.assertEqual(subprocess.check_call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}), 0)

    def assertFail(self, cmd):
        cmd = shlex.split(cmd)
        try:
            subprocess.check_call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={'PATH': os.pathsep.join(['..', os.environ['PATH']])})
        except subprocess.CalledProcessError:
            return

def runCmd(cmd):
    cmd = shlex.split(cmd)
    subprocess.call(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            env={'PATH': os.pathsep.join(['..', os.environ['PATH']])})

def numOfSample():
    return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM sample'],
        stderr=subprocess.PIPE,  env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}))

def numOfVariant(table='variant'):
    return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM {}'.format(table)],
        stderr=subprocess.PIPE,  env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}))
