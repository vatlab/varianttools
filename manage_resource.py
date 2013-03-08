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
import time
from variant_tools.utils import ResourceManager
import base64
from ftplib import FTP

def uploadFile(local_file, remote_file, username, password, logger):
    # upload a local file to houstonbioinformatics.org
    ftp = FTP('www.houstonbioinformatics.org')
    ftp.login(username, password)
    ftp.cwd('vtools')
    d, f = os.path.split(remote_file)
    if d:
        logger.info('CWD {}'.format(d))
        ftp.cwd(d)
    new_f = '{}_{}'.format(f, time.strftime('%b%d', time.gmtime()))
    logger.info('RENAME {} {}'.format(f, new_f))
    ftp.rename(f, new_f)
    logger.info('STOR {}'.format(f))
    ftp.storbinary('STOR {}'.format(f), open(local_file, 'rb'))
    ftp.quit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('--generate_local_manifest', nargs='?', default=False,
        help='''Generate a manifest of local resource files. If a directory is not specified,
            $HOME/.variant_tools will be assumed. The manifest will be saved to 
            MANIFEST_local.txt.''')
    parser.add_argument('--upload', metavar='FILENAME',
        help='''Upload file with name FILENAME to the server. The file should 
            be under the local resource directory ~/.variant_tools.
            A user name and password could be specified via parameters --username
            and --password. The username and password will be saved for future use.''')
    parser.add_argument('--username', nargs='?',
        help='''User name used to connect to vtools.houstonbioinformatics.org.''')
    parser.add_argument('--password', nargs='?',
        help='''Password used to connect to vtools.houstonbioinformatics.org.''')
    #
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('RES')
    if args.generate_local_manifest is not False:
        manager = ResourceManager(logger)
        # --generte_local_manifest without parameter will pass None, which will
        # use ~/.variant_tools.
        manager.scanDirectory(args.generate_local_manifest)
        manager.writeManifest('MANIFEST_local.txt')
        sys.stderr.write('Local manifest has been saved to MANIFEST_local.txt\n') 
    #
    if args.upload:
        if args.username is None or args.password is None:
            try:
                with open(os.path.expanduser('~/.vtools_resource'), 'r') as account:
                    args.username = base64.b64decode(account.readline().decode().strip())
                    args.password = base64.b64decode(account.readline().decode().strip())
            except:
                sys.exit('Please provide username and password.')
            logger.info('Using stored username and password')
        else:
            try:
                with open(os.path.expanduser('~/.vtools_resource'), 'w') as account:
                    account.write('{}\n'.format(base64.b64encode(args.username)))
                    account.write('{}\n'.format(base64.b64encode(args.password)))
            except Exception as e:
                sys.exit('Failed to save username and password: {}'.format(e))

        resource_dir = os.path.expanduser('~/.variant_tools')
        rel_path = os.path.relpath(args.upload, resource_dir)
        manager = ResourceManager(logger)
        # get information about file
        manager.getRemoteManifest()
        manager.addResource(args.upload)
        manager.writeManifest('MANIFEST.tmp')
        uploadFile('MANIFEST.tmp', 'MANIFEST.txt', args.username, args.password, logger)
        uploadFile(args.upload, rel_path, args.username, args.password, logger)


