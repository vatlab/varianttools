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
        with open('run_tests.log', 'a') as fnull:
            self.assertEqual(
                subprocess.check_output(cmd, stderr=fnull,
                    env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}).decode(),
                output)

    def assertSucc(self, cmd):
        cmd = shlex.split(cmd)
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        with open('run_tests.log', 'a') as fnull:
            self.assertEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}), 0)

    def assertFail(self, cmd):
        cmd = shlex.split(cmd)
        try:
            with open('run_tests.log', 'a') as fnull:
                subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                    env={'PATH': os.pathsep.join(['..', os.environ['PATH']])})
        except subprocess.CalledProcessError:
            return

def outputOfCmd(cmd):
    cmd = shlex.split(cmd)
    with open('run_tests.log', 'a') as fnull:
        return subprocess.check_output(cmd, stderr=fnull,
            env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}).decode()
        
def output2list(cmd):
    return map(str, ''.join(outputOfCmd(cmd)).split('\n')[:-1])
    
def runCmd(cmd):
    cmd = shlex.split(cmd)
    with open('run_tests.log', 'a') as fnull:
        subprocess.call(cmd, stdout=fnull, stderr=fnull,
            env={'PATH': os.pathsep.join(['..', os.environ['PATH']])})

def numOfSample():
    with open('run_tests.log', 'a') as fnull:
        return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM sample'],
            stderr=fnull,  env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}))

def numOfVariant(table='variant'):
    with open('run_tests.log', 'a') as fnull:
        return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM {}'.format(table)],
            stderr=fnull,  env={'PATH': os.pathsep.join(['..', os.environ['PATH']])}))
    
def getGenotypes(projname='test', num=8):
    nsamples = numOfSample()
    genotypes = []
    # getGenotypes for 8 samples
    for i in range(num):
        genotypes.append(output2list('vtools execute "select variant_type from {}_genotype.sample_variant_{}"'.format(projname, i+1)))
    return genotypes

def getGenotypeInfo(projname='test', num=8, info=['DP_FMT']):
    nsamples = numOfSample()
    genotypeInfo = []
    # getGenotypeInfo for 8 samples
    for i in range(num):
        genotypeInfo.append(output2list('vtools execute "select {} from {}_genotype.sample_variant_{}"'.format(','.join(info), projname, i+1)))
    return genotypeInfo

def getSamplenames():
    return output2list('vtools execute "select sample_name from sample"')
        
def initTest(level):
    i = 1
    while True:
        runCmd('vtools init test -f') #1
        if i == level: break
        else: i += 1
        runCmd('vtools import_vcf vcf/CEU.vcf.gz --build hg18')
        if i == level: break
        else: i += 1
        runCmd('vtools import_vcf vcf/SAMP1.vcf')
        if i == level: break
        else: i += 1
        runCmd('vtools import_txt --build hg18 --format fmt/basic_hg18 txt/input.tsv --sample_name input.tsv')
        if i == level: break
        else: i += 1        
        runCmd('vtools import_phenotype phenotype/phenotype.txt')
        if i == level: break
        else: i += 1        
        runCmd('vtools use ann/testNSFP.ann') #6
        if i == level: break
        else: i += 1        
        runCmd('vtools select variant \'testNSFP.chr is not null\' -t ns')
        if i == level: break
        else: i += 1        
        runCmd('vtools select ns \'sift_score > 0.95\' -t ns_damaging')
        if i == level: break
        else: i += 1        
        runCmd('vtools select ns \'genename = "PLEKHN1"\'  -t plekhn1')
        if i == level: break
        else: i += 1        
        runCmd('vtools select plekhn1 "polyphen2_score>0.9 or sift_score>0.9" -t d_plekhn1')
        if i == level: break
        else: i += 1           
        runCmd('vtools select variant --samples "filename like \'%CEU%\'" -t CEU')  #11
        if i == level: break
        else: i += 1
        runCmd('vtools select variant --samples "filename like \'%input.tsv\'" -t input_tsv')
        break
