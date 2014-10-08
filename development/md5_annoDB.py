#!/usr/bin/env python
#
# $File: md5_annoDB.py $
# $LastChangedDate: 2013-01-30 18:29:35 -0600 (Wed, 30 Jan 2013) $
# $Rev: 1663 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
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

#
# This file calculate the md5 signature of decompressed .DB.gz file
# and adds the signature to .ann files.
#

import os
import sys
import argparse
import time
from variant_tools.utils import calculateMD5

def fixMD5(annFile, dbFile):
    #
    md5 = calculateMD5(dbFile, partial=True)
    lines = []
    with open(annFile) as ann:
        for line in ann:
            if ann.startswith('direct_URL'):
                direct_URL = line.lsplit('=')[1].strip()
                if 'md5=' in direct_URL:
                    direct_URL = direct_URL.rsplit(None, 1)[0]
                lines.append('direct_URL={}\tmd5={}'.format(direct_URL, md5))
            else:
                lines.append(line.rstrip())
                
    #
    with open(annFile, 'w') as ann:
        ann.write('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('annFile', 
        help='''A .ann file that will be modified for md5 information.''')
    parser.add_argument('dbFile', nargs='*',
        help='''A corresponding .DB file whose md5 signature will be added to field
            direct_URL of the .ann file. If unspecified, $name.DB corresponding to
            $name.ann will be used.'''
    args = parser.parse_args()
    if not args.annFile.endswith('.ann'):
        raise ValueError('The first file should be an .ann file')
    if not args.dbFile:
        args.dbFile = args.annFile[:-3] + '.DB'
    if not os.path.isfile(args.annFile):
        raise ValueError('Missing .ann file {}'.format(args.annFile))
    if not os.path.isfile(args.dbFile):
        raise ValueError('Missing .DB file {}'.format(args.dbFile))
    #
    fixMD5(args.annFile, args.dbFile)


