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
        self.test_command = self.id().split('.')[1]
        with open(self.test_command + '.log', 'a') as fcmd:
            fcmd.write('\n# {}\n# {} \n'.format(self.id().split('.', 1)[-1], 
                '' if self.shortDescription() is None else '\n# '.join(self.shortDescription().split('\n'))))
        self.runCmd('vtools init test -f')

    def compare(self, itemA, itemB, partial=None):
        if isinstance(itemA, str):
            if partial is None:
                self.assertEqual(itemA, itemB)
            elif partial is True:
                self.assertTrue(itemB in itemA)
            elif isinstance(partial, int):
                if partial >= 0:
                    self.assertEqual(itemA.split('\n')[:partial], itemB.split('\n')[:partial])
                elif partial < 0:
                    self.assertEqual(itemA.split('\n')[partial:], itemB.split('\n')[partial:])
            elif callable(partial):
                self.assertEqual(partial(itemA.split('\n')), partial(itemB.split('\n')))
            else:
                raise ValueError('Partial can be an integer, True, or a lambda function')
        else: # a list
            if partial is None:
                self.assertEqual(itemA, itemB)
            elif partial is True:
                if isinstance(itemB, list):
                    self.assertTrue(all([x in itemA for x in itemB]))
                else:
                    self.assertTrue(itemB in itemA)
            elif isinstance(partial, int):
                if partial >= 0:
                    self.assertEqual(itemA[:partial], itemB[:partial])
                elif partial < 0:
                    self.assertEqual(itemA[partial:], itemB[partial:])
            elif callable(partial):
                self.assertEqual(partial(itemA), partial(itemB))
            else:
                raise ValueError('Partial can be an integer, True, or a lambda function')

    def runCmd(self, cmd):
        'Run a command in shell process. Does not check its output or return value.'
        with open(self.test_command + '.log', 'a') as fcmd:
            fcmd.write(cmd + '\n')
        with open(os.devnull, 'w') as fnull:
            return subprocess.check_output(cmd, stderr=fnull, shell=True, env=test_env).decode()

    def assertSucc(self, cmd):
        'Execute cmd and assert its success'
        with open(self.test_command + '.log', 'a') as fcmd:
            fcmd.write('# expect success\n')
            fcmd.write(cmd + '\n')
        with open(os.devnull, 'w') as fnull:
            self.assertEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull, shell=True,
                env=test_env), 0)

    def assertFail(self, cmd):
        'Execute cmd and assert its failure (return non-zero)'
        with open(self.test_command + '.log', 'a') as fcmd:
            fcmd.write('# expect failure\n')
            fcmd.write(cmd + '\n')
        try:
            with open(os.devnull, 'w') as fnull:
                subprocess.check_call(cmd, stdout=fnull, stderr=fnull, shell=True,
                    env=test_env)
        except subprocess.CalledProcessError:
            return


    def assertOutput(self, cmd, output, partial=None):
        '''Compare the output of cmd  to either a string output, or content of output file
        (if output is a filename). cmd can be either a command (string) or a list of commands,
        with output joint together in the latter case. If parameter partial is given,
        a) positive number: texts are split into lines and compares the first partial lines  ([:partial]),
        b) negative number: texts are split into lines and compares the last few lines ([partial:])
        c) a lambda function: texts are split into lines and compare the results returned after
            applying the lambda function
        d) True: if output is a substring of cmd (output in cmd)
            
        NOTE: if output is a file (with pattern output/*.txt) and the file does not exist, this
        function will write command output to it with a warning message. This greatly simplies
        the writing of test functions. Make sure to check if the output is correct though.
        '''
        #
        with open(self.test_command + '.log', 'a') as fcmd:
            if output.startswith('output/') and output.endswith('.txt'):
                fcmd.write('# expect output in {}\n'.format(output))
            else:
                fcmd.write('# expect output string {}\n'.format('\n#'.join(output.split('\n'))))
            if isinstance(cmd, str):
                fcmd.write(cmd + '\n')
            else:
                fcmd.write('\n'.join(cmd) + '\n')
        with open(os.devnull, 'w') as fnull:
            if isinstance(cmd, str):
                cmd_output = subprocess.check_output(cmd, stderr=fnull, env=test_env, shell=True).decode()
            else:
                cmd_output = '\n'.join([subprocess.check_output(c, stderr=fnull, env=test_env, shell=True).decode() for c in cmd])
        if os.path.isfile(output):
            with open(output, 'r') as cf:
                output = cf.read()
        else:
            if output.startswith('output/') and output.endswith('.txt'):
                print('\033[32mWARNING: output file {} does not exist and has just been created.\033[0m'.format(output))
                with open(output, 'w') as cf:
                    cf.write(cmd_output)
                output = cmd_output
        #
        self.compare(cmd_output, output, partial)


    def assertProj(self, numOfSamples=None, numOfVariants=None, sampleNames=None, numOfColumns=None, 
        info=None, genotype=None, genoInfo=None, partial=None):
        '''Check properties of project

        numOfSamples:
            number of samples in the project

        numOfVariants: 
            if a single number is given, assert number of variants in the master variant table.
            Otherwise a dictionary with table name and expected number of variants is checked.

        sampleNames:
            compare sample names with provided list

        numOfColumns:
            Compare number of columns of a table. If a number is given, it is assumed to be the
            'variant' table. Otherwise a dictionary with tablename and expected number of columns
            should be provided.

        info:
            Compare variant info with specified list. This parameter should be a dictionary
            with name of info as key.

        genotype:
            Compare genotype with a provided list. This parameter should be a dictionary with
            sample_id: genotype.
        
        genoInfo:
            Compare genotype info with a provided list. This parameter should be a dictionary
            with key (sample_id, geno_info_name): geno_info.

        partial:
            partial can be True (if specified item is a subset of output list), positive integer
            (compare the first few items), negative number (compare the last few items), or 
            a lambda function (compare result of function call).
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
                    stderr=fnull).decode().strip().split('\n')
                self.compare(sorted([x.strip() for x in sample_names]),
                    sorted([x.strip() for x in sampleNames]), partial=partial)
        if numOfColumns is not None:
            with open(os.devnull, 'w') as fnull:
                if isinstance(numOfColumns, int):
                    proj_num_of_columns = len(subprocess.check_output('vtools execute "PRAGMA table_info(variant)"', shell=True,
                        stderr=fnull).decode().strip().split('\n'))
                    self.assertEqual(proj_num_of_columns, numOfColumns)
                else:
                    for table, number in numOfColumns.items():
                        proj_num_of_columns = len(subprocess.check_output('vtools execute "PRAGMA table_info({})"'.format(table), shell=True,
                            stderr=fnull).decode().strip().split('\n'))
                        self.assertEqual(proj_num_of_columns, number)
        if info is not None:
            with open(os.devnull, 'w') as fnull:
                for field, values in info.items():
                    proj_values = subprocess.check_output('vtools execute "SELECT {} FROM variant"'.format(field), shell=True,
                        stderr=fnull).decode()
                    self.compare([x.strip() for x in proj_values.strip().split('\n')], list(values), partial=partial)
        if genotype is not None:
            with open(os.devnull, 'w') as fnull:
                for table, geno in genotype.items():
                    proj_geno = subprocess.check_output('vtools execute "SELECT GT FROM genotype_{}"'.format(table), shell=True,
                        stderr=fnull).decode()
                    self.compare([x.strip() for x in proj_geno.strip().split('\n')], list(geno), partial=partial)
        if genoInfo is not None:
            with open(os.devnull, 'w') as fnull:
                for table, geno in genoInfo.items():
                    proj_geno = subprocess.check_output('vtools execute "SELECT {} FROM genotype_{}"'.format(table[1], table[0]), shell=True,
                        stderr=fnull).decode()
                    self.compare([x.strip() for x in proj_geno.strip().split('\n')], list(geno), partial=partial)

         
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
