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

def getOrSaveAuthentication(username, password, logger):
    if username is None or password is None:
        try:
            with open(os.path.expanduser('~/.vtools_resource'), 'r') as account:
                return (base64.b64decode(account.readline().decode().strip()),
                    base64.b64decode(account.readline().decode().strip()))
        except:
            sys.exit('Please provide username and password.')
        logger.info('Using stored username and password')
    else:
        try:
            with open(os.path.expanduser('~/.vtools_resource'), 'w') as account:
                account.write('{}\n'.format(base64.b64encode(username)))
                account.write('{}\n'.format(base64.b64encode(password)))
        except Exception as e:
            sys.exit('Failed to save username and password: {}'.format(e))
        return (username, password)

def deprecateFile(filename, username, password, logger):
    # rename a file to another name if it exists
    ftp = FTP('www.houstonbioinformatics.org')
    ftp.login(username, password)
    ftp.cwd('/vtools')
    # go to directory...
    d, f = os.path.split(filename)
    if d:
        ftp.cwd(d)
        logger.info('CWD {}'.format(ftp.pwd()))
    try:
        new_f = '{}_{}.bak'.format(f, time.strftime('%b%d', time.gmtime()))
        ftp.rename(f, new_f)
        logger.info('RENAME {} {}'.format(f, new_f))
    except:
        # if a new file, we do not need to rename the old file
        logger.info('Do not need to remove {} because it does not exist.'.format(f))
        pass
    ftp.quit()

def uploadFile(local_file, remote_file, username, password, logger):
    # upload a local file to houstonbioinformatics.org
    # This is inefficient for multiple files (repeated login and out), but does not
    # really matter
    ftp = FTP('www.houstonbioinformatics.org')
    ftp.login(username, password)
    ftp.cwd('/vtools')
    # go to directory...
    d, f = os.path.split(remote_file)
    if d:
        ftp.cwd(d)
        logger.info('CWD {}'.format(ftp.pwd()))
    try:
        new_f = '{}_{}.bak'.format(f, time.strftime('%b%d', time.gmtime()))
        ftp.rename(f, new_f)
        logger.info('RENAME {} {}'.format(f, new_f))
    except:
        # if a new file, we do not need to rename the old file
        pass
    logger.info('STOR {}'.format(f))
    ftp.storbinary('STOR {}'.format(f), open(local_file, 'rb'))
    ftp.quit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('--generate_local_manifest', nargs='?', const='~/.variant_tools',
        help='''Generate a manifest of local resource files. If a directory is not specified,
            $HOME/.variant_tools will be assumed. The manifest will be saved to 
            MANIFEST_local.txt.''')
    parser.add_argument('--list', nargs='*',
        help='''List all remote and local resources under ~/.variant_tools and mark them
            as identical, missing, new, and modified. If any argument is given, only
            resources with filename containing specified words are displayed.''')
    parser.add_argument('--upload', metavar='FILE', nargs='+',
        help='''Upload specified files to the server. The file should 
            be under the local resource directory ~/.variant_tools.
            A user name and password could be specified via parameters --username
            and --password. The username and password will be saved for future use.''')
    parser.add_argument('--remove', metavar='FILENAME', nargs='+',
        help='''Remove specified files from the online manifest so that it will no
            longer be listed as part of the resource. The file itself, if exists, will
            be renamed but not removed from the server.''') 
    parser.add_argument('--username', nargs='?',
        help='''User name used to connect to vtools.houstonbioinformatics.org.''')
    parser.add_argument('--password', nargs='?',
        help='''Password used to connect to vtools.houstonbioinformatics.org.''')
    #
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger('RES')
    if args.generate_local_manifest is not None:
        manager = ResourceManager(logger)
        print args.generate_local_manifest
        # --generte_local_manifest without parameter will pass None, which will
        # use ~/.variant_tools.
        manager.scanDirectory(args.generate_local_manifest)
        manager.writeManifest('MANIFEST_local.txt')
        sys.stderr.write('Local manifest has been saved to MANIFEST_local.txt\n') 
    elif args.list is not None:
        manager = ResourceManager(logger)
        manager.scanDirectory('~/.variant_tools', args.list)
        local_manifest = {x:y for x,y in manager.manifest.items()}
        manager.manifest.clear()
        manager.getRemoteManifest()
        remote_manifest = manager.manifest
        #
        # compare manifests
        for f, p in sorted(remote_manifest.iteritems()):
            if args.list and not all([x in f for x in args.list]):
                continue
            if f not in local_manifest:
                print('MISSING   {}'.format(f))
            elif p[0] != local_manifest[f][0] or p[1] != local_manifest[f][1]:
                print('MODIFIED  {}'.format(f))
            else:
                print('IDENTICAL {}'.format(f))
        for f, p in sorted(local_manifest.items()):
            if f not in remote_manifest:
                print('NEW      {}'.format(f))
    elif args.upload:
        username, password = getOrSaveAuthentication(args.username, args.password, logger)
        manager = ResourceManager(logger)
        manager.getRemoteManifest()
        resource_dir = os.path.expanduser('~/.variant_tools')
        # get information about file
        for filename in args.upload:
            rel_path = os.path.relpath(filename, resource_dir)
            manager.addResource(filename)
            uploadFile(filename, rel_path, username, password, logger)
        manager.writeManifest('MANIFEST.tmp')
        uploadFile('MANIFEST.tmp', 'MANIFEST.txt', username, password, logger)
    elif args.remove:
        username, password = getOrSaveAuthentication(args.username, args.password, logger)
        manager = ResourceManager(logger)
        manager.getRemoteManifest()
        removed_count = 0
        for filename in args.remove:
            if filename in manager.manifest:
                manager.manifest.pop(filename)
                deprecateFile(filename, username, password, logger)
                logger.info('Remove {} from manifest'.format(filename))
                removed_count += 1
            else:
                logger.warning('{} does not exist in the manifest'.format(filename))
        # upload manifest
        if removed_count > 0:
            manager.writeManifest('MANIFEST.tmp')
            uploadFile('MANIFEST.tmp', 'MANIFEST.txt', username, password, logger)
    else:
        logger.warning('No option has been provided. Please use -h to get a list of actions.')


