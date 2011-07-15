#!/usr/bin/env python
#
# $File: runTests.py $
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


import os, sys, argparse, getpass

test_init = '''
# forcefully create a project
../vtools init sample --force
# remove it
../vtools remove -p sample
# should give an error
../vtools init sample
# summarize
../vtools show
'''

testImportVariant = '''
# cleanup
python ../project.py --project sample --remove
# create a project
python ../project.py --project sample --new @connection
# import VCF file
python ../importVariant.py --project sample --build hg19 --input_files SAMP1.vcf
# nothing should be done because SAMP1.vcf has already been imported
python ../importVariant.py --project sample --build hg19 --input_files SAMP1.vcf
# add another VCF file. 
python ../importVariant.py --project sample --build hg19 --input_files SAMP2.vcf
# add alternative reference genome. You should see liftOver is called to add alternative coordinates
python ../importVariant.py --project sample --alt_build=hg18
# remove project
python ../project.py --project sample --remove
# create a project
python ../project.py --project sample --new
# import a file with more samples, should guess refrence genome correctly
python ../importVariant.py --project sample --input_files CEU.vcf 
# summary?
python ../project.py --project sample --summarize
# add alternative ref genome
python ../importVariant.py --project sample --input_files CEU.vcf --alt_build hg19
# cleanup
python ../project.py --project sample --remove
'''

testManageAnnotation = '''
# cleanup
python ../project.py --project sample --remove
# create a new project
python ../project.py --project sample --new @connection
# see available database
python ../manageAnnotation.py --project sample --summarize
# I should see more databases if I set an annotation_dir
python ../manageAnnotation.py --project sample --annotation_dir ../annotation --summarize
# annotation path should be in project summary
python ../project.py --project sample --summarize
# remove existing database testNSFP if exists
python ../manageAnnotation.py --project sample --removeDB testNSFP
# import testNSFP from source
python ../manageAnnotation.py --project sample --importDB testNSFP
# import testNSFP from one or more files
python ../manageAnnotation.py --project sample --importDB testNSFP --source_files testNSFP.zip
# remove existing database testNSFP
python ../manageAnnotation.py --project sample --removeDB testNSFP
# download testNSFP from web, using URL stored in testNSFP.ann
python ../manageAnnotation.py --project sample --downloadDB testNSFP
# download testNSFP from web, using a direct URL
python ../manageAnnotation.py --project sample --downloadDB https://cge.mdanderson.org/~bpeng1/User/annoDB/testNSFP.DB
# cleanup
python ../project.py --project sample --remove
'''

availableTests = [x[5:] for x in dir() if x.startswith('test_')]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='run test scripts.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars='@')
    parser.add_argument('--engine', choices=['mysql', 'sqlite3'], default='sqlite3',
        help='Database engine, can be mysql or sqlite3.')
    parser.add_argument('--host', default='localhost', 
        help='The MySQL server that hosts the project databases. Unused if --engine=sqlite3.')
    parser.add_argument('--user', default=getpass.getuser(),
        help='User name to the MySQL server. Default to current username. Unused if --engine=sqlite3.')
    parser.add_argument('--passwd',
        help='Password to the MySQL server. Unused if --engine=sqlite3.')
    parser.add_argument('-v', '--verbosity', default='1', choices=['0','1','2'],
        help='Produce error and warning (0), info (1) and debug (2) information.')
    parser.add_argument('--pause', action='store_true', 
        help='whether or not stop at each test')
    parser.add_argument('test', default='all', choices=['all'] + availableTests,
        help='tests to run. Default to all.')
    args = parser.parse_args()
    tests = availableTests if args.test == 'all' else [args.test]
    with open('connection', 'w') as conn:
        conn.write('--engine={}\n'.format(args.engine))
        if args.engine == 'mysql':
            conn.write('--host={}\n'.format(args.host))
            conn.write('--user={}\n'.format(args.user))
            conn.write('--passwd={}\n'.format(args.passwd))
    for test in tests:
        script = eval('test_' + test)
        for line in script.split('\n'):
            if line.strip() == '':
                continue
            elif line.startswith('#'):
                print line.strip()
            else:
                cmd = line.strip() 
                cmd += ' -v={}'.format(args.verbosity)
                print '%', cmd
                os.system(cmd)
                if args.pause:
                    print
                    raw_input('Press any key to continue...')
                    print '\n'*5
  
