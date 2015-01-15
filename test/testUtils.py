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

# we will need to install vtools first, and we put /usr/local/bin
# before any system path that might have path to local vtools ...
test_env = {'PATH': os.pathsep.join(['/usr/bin', '/usr/local/bin', os.environ['PATH']]),
   'PYTHONPATH': os.environ.get('PYTHONPATH', '')}

class PrettyPrinter:
    def __init__(self, out):
        ''' delimiter: use specified field to separate fields
            cache: only use the first cache lines to get column width
        '''
        self.width = []
        self.rows = []
        self.out = out

    def cache(self, d):
        data = [x.strip() for x in d]
        self.rows.append(data)
        if not self.width:
            self.width = [len(str(x)) for x in data]
        else:
            self.width = [max(y, len(str(x))) for y,x in zip(self.width, data)]

    def write(self):
        self.width[-1] = 0
        with open(self.out, 'w') as o:
            o.write('\n'.join([
            '\t'.join(
                [str(col).ljust(width) for col, width in zip(row, self.width)])
            for row in self.rows]) + '\n')

class ProcessTestCase(unittest.TestCase):
    'A subclass of unittest.TestCase to handle process output'
# used to be empty; force delete
    def setUp(self):
        'Clear any existing project'
        runCmd('vtools remove project')

    def tearDown(self):
        'Clear any existing project'
        #runCmd('vtools remove project')

    def fileContent(self, file):
        with open(file) as ifile:
            cont = ifile.read()
        return cont

# compare if the command output is what we want
    def assertOutput(self, cmd, output=None, numOfLines=0, file=None, skip=None):
        cmd = shlex.split(cmd)
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        if skip is not None:
             if numOfLines != 0 or file is not None:
                 raise ValueError('skip is not implemented with other options')
        with open('run_tests.log', 'a') as fnull: 
            if output is None and file is None:
                raise ValueError('Please specify your OUTPUT or give your file name of the output')
            if numOfLines == 0:
                if file is None:        
                    if skip is None:
                        self.assertEqual(
                             subprocess.check_output(cmd, stderr=fnull,
                             env=test_env).decode(), output)
                    else:
                        text = subprocess.check_output(cmd, stderr=fnull,
                             env=test_env).decode().split('\n')
                        text = text[:skip-1] + text[skip:]
                        self.assertEqual('\n'.join(text), output)
                elif file is not None:
                        with open(file) as f:
                             content = f.read()
                        self.assertEqual(
                             subprocess.check_output(cmd, stderr=fnull,
                             env=test_env).decode(),
                             content)
            elif numOfLines > 0:            
                if file is None:            
                    self.assertEqual(
                         '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                         env=test_env).decode().split('\n')[:numOfLines]),
                         output)
                elif file is not None:
                    self.assertEqual(
                         '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                         env=test_env).decode().split('\n')[:numOfLines]),
                         '\n'.join(self.fileContent(file).split('\n')[:numOfLines]))
            elif numOfLines < 0:            
                if file is None:            
                    self.assertEqual(
                         '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                         env=test_env).decode().split('\n')[numOfLines:]),
                         output)
                elif file is not None:
                    self.assertEqual(
                         '\n'.join(subprocess.check_output(cmd, stderr=fnull,
                         env=test_env).decode().split('\n')[numOfLines:]),
                         '\n'.join(self.fileContent(file).split('\n')[numOfLines:]))
    def assertSucc(self, cmd):
        cmd = shlex.split(cmd)
        # '..' is added to $PATH so that command (vtool) that is in the current directory # can be executed.
        with open('run_tests.log', 'a') as fnull:
            self.assertEqual(subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                env=test_env), 0)

    def assertFail(self, cmd):
        cmd = shlex.split(cmd)
        try:
            with open('run_tests.log', 'a') as fnull:
                subprocess.check_call(cmd, stdout=fnull, stderr=fnull,
                    env=test_env)
        except subprocess.CalledProcessError:
            return

def outputOfCmd(cmd):
    cmd = shlex.split(cmd)
    with open('run_tests.log', 'a') as fnull:
        return subprocess.check_output(cmd, stderr=fnull,
            env=test_env).decode()
        

def output2list(cmd, space2tab=False):
    txt = list(map(str, ''.join(outputOfCmd(cmd)).split('\n')[:-1]))
    if space2tab:
        txt = ['\t'.join(t.split())  for t in txt]
    return txt
        
def runCmd(cmd, shell=False):
    if not shell:
        cmd = shlex.split(cmd)
    with open('run_tests.log', 'a') as fnull:
        subprocess.call(cmd, stdout=fnull, stderr=fnull, shell=shell,
            env=test_env)

def numOfSample():
    with open('run_tests.log', 'a') as fnull:
        return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM sample'],
            stderr=fnull,  env=test_env))

def numOfVariant(table='variant'):
    with open('run_tests.log', 'a') as fnull:
        return int(subprocess.check_output(['vtools', 'execute', 'SELECT count(1) FROM {}'.format(table)],
            stderr=fnull,  env=test_env))
    
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

def getSamplenames():
    return output2list('vtools execute "select sample_name from sample"')
        
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
