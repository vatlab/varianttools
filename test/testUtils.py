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
# sys.path.append("/Users/jma7/anaconda/envs/VariantTools/lib/python3.6/site-packages/variant_tools-3.0.0.dev0-py3.6-macosx-10.7-x86_64.egg")
from variant_tools.geno_store import *
from variant_tools.accessor import *
# from "variant_tools-3.0.0.dev0-py3.6-macosx-10.7-x86_64.egg".variant_tools.geno_store import *
# from "variant_tools-3.0.0.dev0-py3.6-macosx-10.7-x86_64.egg".variant_tools.accessor import *

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
        if os.environ.get("STOREMODE") is not None:
            self.storeMode=os.getenv("STOREMODE")
            self.runCmd('vtools init test -f --store '+self.storeMode)
        else:
            self.storeMode="hdf5"
            os.environ["STOREMODE"]="hdf5"
            self.runCmd('vtools init test -f --store hdf5')
        if os.environ.get("LOCALRESOURCE") is not None:
            self.local_resource=os.getenv("LOCALRESOURCE")
            self.runCmd('vtools admin --set_runtime_option local_resource='+self.local_resource)

    def compare(self, itemA, itemB, partial=None, negate=None):
        if not isinstance(itemA, list):
            if partial is None:
                if negate:
                    self.assertNotEqual(itemA.strip(), itemB.strip())
                else:
                    self.assertEqual(itemA.strip(), itemB.strip())
            elif partial is True:
                if negate:
                    self.assertFalse(itemB in itemA)
                else:
                    self.assertTrue(itemB in itemA)
            elif isinstance(partial, int):
                if partial >= 0:
                    if negate:
                        self.assertNotEqual(itemA.split('\n')[:partial], itemB.split('\n')[:partial])
                    else:
                        self.assertEqual(itemA.split('\n')[:partial], itemB.split('\n')[:partial])
                elif partial < 0:
                    if negate:
                        self.assertNotEqual(itemA.split('\n')[partial:], itemB.split('\n')[partial:])
                    else:
                        self.assertEqual(itemA.split('\n')[partial:], itemB.split('\n')[partial:])
            elif callable(partial):
                if negate:
                    self.assertNotEqual(partial(itemA.split('\n')), partial(itemB.split('\n')))
                else:
                    self.assertEqual(partial(itemA.split('\n')), partial(itemB.split('\n')))
            else:
                raise ValueError('Partial can be an integer, True, or a lambda function')
        else: # a list
            if partial is None:
                if negate:
                    self.assertNotEqual(itemA, itemB)
                else:
                    self.assertEqual(itemA, itemB)
            elif partial is True:
                if isinstance(itemB, list):
                    if negate:
                        self.assertFalse(any([x in itemA for x in itemB]))
                    else:
                        self.assertTrue(all([x in itemA for x in itemB]))
                else:
                    if negate:
                        self.assertFalse(itemB in itemA)
                    else:
                        self.assertTrue(itemB in itemA)
            elif isinstance(partial, int):
                if partial >= 0:
                    if negate:
                        self.assertNotEqual(itemA[:partial], itemB[:partial])
                    else:
                        self.assertEqual(itemA[:partial], itemB[:partial])
                elif partial < 0:
                    if negate:
                        self.assertNotEqual(itemA[partial:], itemB[partial:])
                    else:
                        self.assertEqual(itemA[partial:], itemB[partial:])
            elif callable(partial):
                if negate:
                    self.assertNotEqual(partial(itemA), partial(itemB))
                else:
                    self.assertEqual(partial(itemA), partial(itemB))
            else:
                raise ValueError('Partial can be an integer, True, or a lambda function')

    def runCmd(self, cmd, ret='string'):
        'Run a command in shell process. Does not check its output or return value.'
        with open(self.test_command + '.log', 'a') as fcmd:
            fcmd.write(cmd + '\n')
        if ret == 'string':
            return subprocess.check_output(cmd, shell=True, env=test_env).decode()
        elif ret == 'list':
            return subprocess.check_output(cmd, shell=True, env=test_env).decode().strip().split('\n')
        else:
            raise ValueError('Unrecognized return type for command runCmd {}'.format(ret))

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
            # we write failed command to .log with # so that these commands are not executed
            # upon examination/rerunning
            fcmd.write('# expect failure\n#')
            fcmd.write(cmd + '\n')
        try:
            with open(os.devnull, 'w') as fnull:
                self.assertNotEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull, shell=True,
                    env=test_env), 0)
        except subprocess.CalledProcessError:
            return


    def assertOutput(self, cmd, output, partial=None, negate=None):
        '''Compare the output of cmd to either a string output, a list of strings, or content of
        output file (if output is a filename). cmd can be either a command (string) or a list of commands,
        with output joint together in the latter case. The output of the command will be converted to
        a list (split by newline) if output is a list.
        
        If parameter partial is given,
        a) positive number: texts are split into lines and compares the first partial lines  ([:partial]),
        b) negative number: texts are split into lines and compares the last few lines ([partial:])
        c) a lambda function: texts are split into lines and compare the results returned after
            applying the lambda function
        d) True: if output is a substring of cmd (output in cmd)
            
        if negate is True, test for negative assertaion (e.g. not equal, not include etc)

        NOTE: if output is a file (with pattern output/*) and the file does not exist, this
        function will write command output to it with a warning message. This greatly simplies
        the writing of test functions. Make sure to check if the output is correct though.
        '''
        #
        with open(os.devnull, 'w') as fnull:
            if isinstance(cmd, str):
                cmd_output = subprocess.check_output(cmd, stderr=fnull, env=test_env, shell=True).decode()
            else:
                cmd_output = '\n'.join([subprocess.check_output(c, stderr=fnull, env=test_env, shell=True).decode() for c in cmd])
            if isinstance(output, list):
                cmd_output = cmd_output.strip().split('\n')
        with open(self.test_command + '.log', 'a') as fcmd:
            if isinstance(output, list):
                fcmd.write('# expect output is a list of {}\n'.format(', '.join(output)))
            elif output.startswith('output/'):
                fcmd.write('# expect output in {} with first 10 lines\n'.format(output))
                if os.path.isfile(output):
                    with open(output, 'r') as cf:
                        output = cf.read()
                        fcmd.write('# ' + '\n# '.join(output.split('\n')[:10]) + '\n')
                else:
                    print('\033[32mWARNING: output file {} does not exist and has just been created.\033[0m'.format(output))
                    with open(output, 'w') as cf:
                        cf.write(cmd_output)
                    output = cmd_output
            else:
                fcmd.write('# expect output {} \n'.format('\n#'.join(output.strip().split('\n'))))
            if isinstance(cmd, str):
                fcmd.write(cmd + '\n')
            else:
                fcmd.write('\n'.join(cmd) + '\n')
        #
        self.compare(cmd_output, output, partial, negate=negate)


    def assertProj(self, numOfSamples=None, numOfVariants=None, numOfGenotype=None, sampleNames=None, numOfColumns=None, 
        info=None, genotype=None, genoInfo=None, hasTable=None, tableDesc=None, partial=None, negate=None):
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

        numOfGenotype:
            Compare number of genotype in specified sample. The parameter should be a dictionary
            with sample_id as keys.

        genotype:
            Compare genotype with a provided list. This parameter should be a dictionary with
            sample_id: genotype.
        
        genoInfo:
            Compare genotype info with a provided list. This parameter should be a dictionary
            with key (sample_id, geno_info_name): geno_info.

        hasTable:
            If the project has the specified table(s). hasTable should be a string or a list of
            strings.

        tableDesc:
            Compare the description of the table. tableDesc should be dictionary with table name
            as key.

        partial:
            partial can be True (if specified item is a subset of output list), positive integer
            (compare the first few items), negative number (compare the last few items), or 
            a lambda function (compare result of function call).
        
        negate:
            If set to True, reverse the test (e.g. assert not equal, not include etc)
        '''
        if numOfSamples is not None:
            with open(os.devnull, 'w') as fnull:
                proj_num_of_sample = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM sample"', shell=True,
                    stderr=fnull).decode()
            if negate:
                self.assertNotEqual(int(proj_num_of_sample), numOfSamples)
            else:
                self.assertEqual(int(proj_num_of_sample), numOfSamples)
        if numOfVariants is not None:
            with open(os.devnull, 'w') as fnull:
                if isinstance(numOfVariants, int):
                    proj_num_of_variants = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM variant"', shell=True,
                        stderr=fnull).decode()
                    if negate:
                        self.assertNotEqual(int(proj_num_of_variants), numOfVariants)
                    else:
                        self.assertEqual(int(proj_num_of_variants), numOfVariants)
                else:
                    for table, number in numOfVariants.items():
                        proj_num_of_variants = subprocess.check_output('vtools execute "SELECT COUNT(1) FROM {}"'.format(table), shell=True,
                            stderr=fnull).decode()
                        if negate:
                            self.assertNotEqual(int(proj_num_of_variants), number)
                        else:
                            self.assertEqual(int(proj_num_of_variants), number)
        if sampleNames is not None:
            with open(os.devnull, 'w') as fnull:
                sample_names = subprocess.check_output('vtools execute "SELECT sample_name FROM sample"', shell=True,
                    stderr=fnull).decode().strip().split('\n')
                self.compare(sorted([x.strip() for x in sample_names]),
                    sorted([x.strip() for x in sampleNames]), partial=partial, negate=negate)
        if numOfColumns is not None:
            with open(os.devnull, 'w') as fnull:
                if isinstance(numOfColumns, int):
                    proj_num_of_columns = len(subprocess.check_output('vtools execute "PRAGMA table_info(variant)"', shell=True,
                        stderr=fnull).decode().strip().split('\n'))
                    if negate:
                        self.assertNotEqual(proj_num_of_columns, numOfColumns)
                    else:
                        self.assertEqual(proj_num_of_columns, numOfColumns)
                else:
                    for table, number in numOfColumns.items():
                        proj_num_of_columns = len(subprocess.check_output('vtools execute "PRAGMA table_info({})"'.format(table), shell=True,
                            stderr=fnull).decode().strip().split('\n'))
                        if negate:
                            self.assertNotEqual(proj_num_of_columns, number)
                        else:
                            self.assertEqual(proj_num_of_columns, number)
        if info is not None:
            with open(os.devnull, 'w') as fnull:
                for field, values in info.items():
                    proj_values = subprocess.check_output('vtools execute "SELECT {} FROM variant"'.format(field), shell=True,
                        stderr=fnull).decode()
                    self.compare([x.strip() for x in proj_values.strip().split('\n')], list(values), partial=partial, negate=negate)
        if numOfGenotype is not None:
            with open(os.devnull, 'w') as fnull:
                for table, numGeno in numOfGenotype.items():
                    print(table,numGeno)
                    if self.storeMode=="sqlite":
                        proj_num_geno = subprocess.check_output('vtools execute "SELECT count(*) FROM genotype_{}"'.format(table), shell=True,
                            stderr=fnull).decode()
                    elif self.storeMode=="hdf5":
                        fileResult = subprocess.check_output('vtools execute "SELECT HDF5 FROM sample WHERE sample_id ={}"'.format(table), shell=True,
                            stderr=fnull).decode()
                        HDF5FileName=fileResult.rstrip()
                        storageEngine=Engine_Storage.choose_storage_engine(HDF5FileName)
                        proj_geno,numCount=storageEngine.num_variants(table)
                        proj_num_geno=proj_geno
                    
                    if negate:
                        self.assertNotEqual(int(proj_num_geno), numGeno)
                    else:
                        self.assertEqual(int(proj_num_geno), numGeno)
        if genotype is not None:
            with open(os.devnull, 'w') as fnull:
                for table, geno in genotype.items():
                    if self.storeMode=="sqlite":
                        proj_geno = subprocess.check_output('vtools execute "SELECT GT FROM genotype_{}"'.format(table), shell=True,
                            stderr=fnull).decode()
                        self.compare([int(x.strip()) for x in proj_geno.strip().split('\n')], list([int(x) for x in geno]), partial=partial, negate=negate)
                    elif self.storeMode=="hdf5":
                        fileResult = subprocess.check_output('vtools execute "SELECT HDF5 FROM sample WHERE sample_id ={}"'.format(table), shell=True,
                            stderr=fnull).decode()
                        HDF5FileName=fileResult.rstrip()
                        accessEngine=Engine_Access.choose_access_engine(HDF5FileName)
                        # proj_geno=accessEngine.get_geno_by_sample_ID(table,"GT_geno")
                        proj_geno=[]
                        for rownames,colnames,genoinfo in accessEngine.get_all_genotype([table]):
                            for idx,rowname in enumerate(rownames):
                                genotype=genoinfo[idx]
                                if np.isnan(genotype):
                                    genotype=-1
                                proj_geno.append([rowname,genotype])
                        proj_geno=np.array(proj_geno)


                        self.compare([int(x[1]) for x in proj_geno], list([int(x) for x in geno]), partial=partial, negate=negate)
                        

        if genoInfo is not None:
            with open(os.devnull, 'w') as fnull:
                for table, geno in genoInfo.items():
                    if self.storeMode=="sqlite":
                        proj_geno = subprocess.check_output('vtools execute "SELECT {} FROM genotype_{}"'.format(table[1], table[0]), shell=True,
                            stderr=fnull).decode()
                        self.compare([x.strip() for x in proj_geno.strip().split('\n')], list(geno), partial=partial, negate=negate)
                    elif self.storeMode=="hdf5":
                        fileResult = subprocess.check_output('vtools execute "SELECT HDF5 FROM sample WHERE sample_id ={}"'.format(table[0]), shell=True,
                            stderr=fnull).decode()
                        HDF5FileName=fileResult.rstrip()
                        accessEngine=Engine_Access.choose_access_engine(HDF5FileName)
                        # proj_geno=accessEngine.get_geno_by_sample_ID(table[0],table[1])
                        proj_geno=[]
                        for rownames,colnames,genoinfo in accessEngine.get_all_genotype([table[0]]):
                            for idx,rowname in enumerate(rownames):
                                genotype=genoinfo[idx]
                                if np.isnan(genotype):
                                    genotype=-1
                                proj_geno.append([rowname,genotype])
                        proj_geno=np.array(proj_geno)

                        self.compare([int(x[1]) for x in proj_geno], list([int(x) for x in geno]), partial=partial, negate=negate)

        if hasTable is not None:
            with open(os.devnull, 'w') as fnull:
                proj_tables = subprocess.check_output('vtools show tables -v0', shell=True).decode().strip().split('\n')
                if isinstance(hasTable, list):
                    for table in hasTable:
                        if negate:
                            self.assertFalse(table in proj_tables)
                        else:
                            self.assertTrue(table in proj_tables)
                else:
                    if negate:
                        self.assertFalse(hasTable in proj_tables)
                    else:
                        self.assertTrue(hasTable in proj_tables)
        if tableDesc is not None:
            with open(os.devnull, 'w') as fnull:
                for table, desc in tableDesc.items():
                    proj_table_desc = subprocess.check_output("vtools show table '{}'".format(table), shell=True, stderr=fnull).decode()
                    proj_table_desc = proj_table_desc.strip().split('\n')[1].split(':', 1)[-1].strip()
                    self.compare(proj_table_desc, desc, partial=partial, negate=negate)
                

        
