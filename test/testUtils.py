#!/usr/bin/env python
#
# $File: ProcessTestCase $
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
import unittest
import shlex, subprocess
import sys

test_env = None
# Sometimes a modified environment for testing is needed.
#
#{'PATH': os.pathsep.join(['/usr/bin', '/usr/local/bin', os.environ['PATH']]),
#   'PYTHONPATH': os.environ.get('PYTHONPATH', '')}

class ProcessTestCase(unittest.TestCase):
    'A subclass of unittest.TestCase to handle process output'
    def setUp(self):
        'Create a new project. This will be called for each new test'
        self.runCmd('vtools init test -f')

    def runCmd(self, cmd):
        'Run a command in shell process. Does not check its output or return value.'
        with open('run_tests.log', 'a') as fcmd:
            fcmd.write(cmd + '\n\n')
        with open(os.devnull, 'w') as fnull:
            subprocess.call(cmd, stdout=fnull, stderr=fnull, shell=True,
                env=test_env)

    def assertOutput(self, cmd, output, partial=None):
        '''Compare the output of cmd to either a string output, or content of output file
        (if output is a filename). If parameter partial is given, texts are split into lines
        and compares the first partial lines if it is a positive number (slice :partial), the
        last few lines if it is a negative number (slice partial:), or the result of it if
        partial is a lambda function.
            
        NOTE: if output is a file (with pattern test/*.txt) and the file does not exist, this
        function will write command output to it with a warning message. This greatly simplies
        the writing of test functions. Make sure to check if the output is correct though.
        '''
        #
        with open('run_tests.log', 'a') as fcmd:
            fcmd.write(cmd + '\n\n')
        with open(os.devnull, 'w') as fnull:
            cmd_output = subprocess.check_output(cmd, stderr=fnull, env=test_env, shell=True).decode()
        if os.path.isfile(output):
            with open(output, 'r') as cf:
                output = cf.read()
        else:
            if output.startswith('test/') and output.endswith('.txt'):
                print('\041[32mWARNING: output file {} does not exist and has just been created.\033[0m'.format(output))
                with open(output, 'w') as cf:
                    cf.write(cmd_output)
        #
        if lines is None:
            self.assertEqual(cmd_output, output)
        else:
            # compare slice?
            if isinstance(partial, int):
                if partial >= 0:
                    self.assertEqual(
                        cmd_output.split('\n')[:partial],
                        output.split('\n')[:partial])
                elif partial < 0:
                    self.assertEqual(
                        cmd_output.split('\n')[partial:],
                        output.split('\n')[partial:])
            elif callable(partial):
                self.assertEqual(
                    partial(cmd_output.split('\n')),
                    partial(output.split('\n')))
            else:
                raise ValueError('Partial can be an integer or a lambda function')

    def assertSucc(self, cmd):
        'Execute cmd and assert its success'
        with open('run_tests.log', 'a') as fcmd:
            fcmd.write(cmd + '\n\n')
        with open(os.devnull, 'w') as fnull:
            self.assertEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull, shell=True,
                env=test_env), 0)

    def assertFail(self, cmd):
        'Execute cmd and assert its failure (return non-zero)'
        with open('run_tests.log', 'a') as fcmd:
            fcmd.write(cmd + '\n\n')
        try:
            with open(os.devnull, 'w') as fnull:
                subprocess.check_call(cmd, stdout=fnull, stderr=fnull, shell=True,
                    env=test_env)
        except subprocess.CalledProcessError:
            return

    def assertProj(self, numOfSamples=None, numOfVariants=None, sampleNames=None):
        '''Check properties of project

        numOfSamples:
            number of samples in the project

        numOfVariants: 
            if a single number is given, assert number of variants in the master variant table.
            Otherwise a dictionary with table name and expected number of variants is checked.
        '''
        if numOfSamples is not None:
            with open(os.devnull, 'w') as fnull:
                proj_num_of_sample = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM sample"', shell=True,
                    stderr=fnull)
            self.assertEqual(int(proj_num_of_sample), numOfSamples)
        if numOfVariants is not None:
            with open(os.devnull, 'w') as fnull:
                if isinstance(numOfVariants, int):
                    proj_num_of_variants = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM variant"', shell=True,
                        stderr=fnull)
                    self.assertEqual(int(proj_num_of_variants), numOfVariants)
                else:
                    for table, number in numOfVariants.items():
                        proj_num_of_variants = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM {}"'.format(table), shell=True,
                            stderr=fnull)
                        self.assertEqual(int(proj_num_of_variants), number)
        if sampleNames is not None:
            with open(os.devnull, 'w') as fnull:
                sample_names = subprocess.check_output('vtools execute "SELECT sample_name FROM sample"', shell=True,
                    stderr=fnull).decode().split('\n')
                self.assertEqual(sorted([x.strip() for x in sample_names]),
                    sorted([x.strip() for x in sampleNames])

          
def getGenotypes(projname='test', num=8):
    nsamples = numOfSample()
    genotypes = []
    # getGenotypes for 8 samples
    for i in range(num):
        genotypes.append(output2list('vtools execute "select GT from {}_genotype.genotype_{}"'.format(projname, i+1)))
    return genotypes

def getGenotypeInfo(projname='test', num=8, info=['DP_geno']):
    nsamples = numOfSample()
    genotypeInfo = []
    # getGenotypeInfo for 8 samples
    for i in range(num):
        genotypeInfo.append(output2list('vtools execute "select {} from {}_genotype.genotype_{}"'.format(','.join(info), projname, i+1)))
    return genotypeInfo

def initTest(level):
    i = 1
    while True:
        runCmd('vtools init test -f') #1
        runCmd('vtools admin --set_runtime_option term_width=78') #1
        if i == level: break
        else: i += 1
        runCmd('vtools import vcf/CEU.vcf.gz --build hg18')
        if i == level: break
        else: i += 1
        runCmd('vtools import vcf/SAMP1.vcf')
        if i == level: break
        else: i += 1
        runCmd('vtools import --format fmt/basic_hg18 txt/input.tsv --build hg18 --sample_name input.tsv')
        if i == level: break
        else: i += 1        
        runCmd('vtools phenotype --from_file phenotype/phenotype.txt')
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
        runCmd('vtools select ns \'genename = "PLEKHN1"\' -t plekhn1')
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
