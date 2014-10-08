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

def fixMD5(annFile):
    #
    lines = []
    d, f = os.path.split(annFile)
    with open(annFile) as ann:
        for line in ann:
            if line.startswith('direct_url'):
                direct_URL = line.split('=')[1].strip()
                if '\t' in direct_URL:
                    direct_URL = direct_URL.rsplit('\t', 1)[0]
                # file name
                if not direct_URL.endswith('.DB.gz'):
                    print('No .DB.gz file is specified in direct_URL: {}'.format(direct_URL))
                    return
                gzFile = os.path.split(direct_URL)[-1]
                unzippedFile = os.path.join(d, gzFile[:-3])
                if not os.path.isfile(unzippedFile):
                    print('Cannot find unzipped file {}'.format(unzippedFile))
                    return
                md5 = calculateMD5(unzippedFile, partial=True)
                lines.append('direct_url={}\t{}'.format(direct_URL, md5))
                print('md5 signature of {} added to {}'.format(unzippedFile, annFile))
            else:
                lines.append(line.rstrip())
                
    #
    with open(annFile, 'w') as ann:
        ann.write('\n'.join(lines))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('annFiles', nargs='+',
        help='''One or more .ann files that will be modified for md5 information. For 
            each file with name $name.ann, a $name.DB file is expected. ''')
    args = parser.parse_args()
    for ann in args.annFiles:
        if not ann.endswith('.ann'):
            raise ValueError('The first file should be an .ann file')
        if not os.path.isfile(ann):
            raise ValueError('Missing .ann file {}'.format(ann))
        #
        fixMD5(ann)


