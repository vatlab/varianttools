#!/usr/bin/env python
#
# $File: manage_resource.py $
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

import os
import sys
import argparse
import logging
from variant_tools.utils import ResourceManager


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('--generate_local_manifest', nargs='?', default=False,
        help='''Generate a manifest of local resource files. If a directory is not specified,
            $HOME/.variant_tools will be assumed. The manifest will be saved to 
            MANIFEST_local.txt.''')
    #
    args = parser.parse_args()
    logging.basicConfig()
    logger = logging.getLogger()
    if args.generate_local_manifest is not False:
        manager = ResourceManager(logger)
        # --generte_local_manifest without parameter will pass None, which will
        # use ~/.variant_tools.
        manifest = manager.generateLocalManifest('MANIFEST_local.txt', args.generate_local_manifest)
        sys.stderr.write('Local manifest has been saved to MANIFEST_local.txt\n') 


