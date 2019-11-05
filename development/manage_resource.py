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
import time
import pexpect
from variant_tools.utils import ResourceManager, env, calculateMD5
import base64

repoURL = 'bpeng1@bcbweb'
repoDir = '/var/www/html/Software/VariantTools/'

dryrun = False

def pexpect_run(cmd, user, passwd):
    p = pexpect.spawn(cmd)
    while True:
        try:
            i = p.expect([
                "(?i)are you sure you want to continue connecting",
                "name:",
                "[pP]assword:",
                pexpect.EOF
            ], timeout=5)

            if i == 0:
                p.sendline('yes')
            elif i == 1:
                p.sendline(user)
            elif i == 2:
                p.sendline(passwd)
            elif i == 3:
                p.close()
                return p.exitstatus
        except Exception as e:
            print(f'Failed to run command {cmd}: {e}')
            return 1
        return 1

def remoteDo(cmd, user, passwd):
    'Run command on the server under directory #repoDir '
    if type(cmd) in (tuple, list):
        cmd = 'ssh {} "cd {}; {}"'.format(repoURL, repoDir, ';'.join(cmd))
    else:
        cmd = 'ssh {} "cd {}; {}"'.format(repoURL, repoDir, cmd)
    env.logger.info(cmd)
    if not dryrun:
        pexpect_run(cmd, user=user, passwd=passwd)

def deprecateFile(filename, repo, user, passwd):
    'Move a file to the deprecated folder with time stamp'
    d, _ = os.path.split(filename)
    remoteDo([
        '[ ! -d deprecated/{0} ] && mkdir -p deprecated/{0}'.format(d),
        '[ -f {0}/{1} ] && mv {0}/{1} deprecated/{1}_{2}'.format(repo, filename, time.strftime('%b%d', time.gmtime()))],
        user=user, passwd=passwd)

def uploadFile(local_file, remote_file, repo, user, passwd):
    d, f = os.path.split(remote_file)
    deprecateFile(remote_file, repo, user=user, passwd=passwd)
    remoteDo('[ ! -d {0}/{1} ] && mkdir -p {0}/{1}'.format(repo, d), user=user, passwd=passwd)
    cmd = 'scp {} {}:{}/{}/{}'.format(local_file, repoURL, repoDir, repo, remote_file)
    env.logger.info(cmd)
    if not dryrun:
        pexpect_run(cmd, user=user, passwd=passwd)
        remoteDo('chmod o+r {}/{}/{}'.format(repoDir, repo, remote_file), user=user, passwd=passwd)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Manage variant tools resources''')
    parser.add_argument('--generate_local_manifest', nargs='?', const=env.local_resource,
        help='''Generate a manifest of local resource files. If a directory is not specified,
            $HOME/.variant_tools will be assumed. The manifest will be saved to
            MANIFEST_local.txt.''')
    parser.add_argument('--list', nargs='*',
        help='''List all remote and local resources under ~/.variant_tools and mark them
            as identical, missing, new, and modified. If any argument is given, only
            resources with filename containing specified words are displayed.''')
    parser.add_argument('--update',  nargs='?', metavar='TYPE',
        const='current',
        choices=['current', 'all', 'existing', 'hg18', 'hg19', 'annotation', 'format',
           'snapshot', 'pipeline'],
        help='''Download resources of specified type, which can be 'current' (latest version
            of all resources excluding snapshots), 'all' (all resources including obsolete
            databases), 'existing' (only update resources that exist locally),
            'hg18' or 'hg19' (all resources for reference genome hg18 or hg19),
            'annotation' (all current annotation databases), 'format' (all formats), and
            'snapshot' (all online snapshots). Identical resources that are available locally
            (under ~/.variant_tools or runtime option $local_resource) are ignored. Note that
            option 'all' will download all versions of annotation databases which can be
            slow and take a lot of disk spaces.''')
    parser.add_argument('--repo', metavar='REPOSITORY', default='repository',
        choices=['repository', 'archive', 'bioinfo', 'devel'],
        help='''What repository the action will be applied to. The default value is
            the main repository. If you would like to move a file from one repository
            to another, you can use options --remove and --upload with different
            --repo.''')
    parser.add_argument('--upload', metavar='FILE', nargs='+',
        help='''Upload specified files to the server. The file should
            be under the local resource directory (default to ~/.variant_tools).''')
    parser.add_argument('--remove', metavar='FILENAME', nargs='+',
        help='''Remove specified files from the online manifest so that it will no
            longer be listed as part of the resource. The file itself, if exists, will
            be renamed but not removed from the server.''')
    parser.add_argument('--local', metavar='LOCAL', default=env.local_resource,
        help='''Local resource directory, default to env.local_resource''')
    parser.add_argument('--user', help='user name')
    parser.add_argument('--passwd', help='password')
    # this set up and use default temporary directory
    env.temp_dir = None
    args = parser.parse_args()
    env.logger = None  # no log file
    if args.generate_local_manifest is not None:
        manager = ResourceManager()
        # --generte_local_manifest without parameter will pass None, which will
        # use ~/.variant_tools.
        manager.scanDirectory(args.generate_local_manifest)
        manager.writeManifest('MANIFEST_local.txt')
        sys.stderr.write('Local manifest has been saved to MANIFEST_local.txt\n')
    elif args.list is not None:
        manager = ResourceManager()
        manager.scanDirectory(args.local, args.list)
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
                print('NEW       {}'.format(f))
    elif args.update:
        res = ResourceManager()
        res.getRemoteManifest()
        res.selectFiles(resource_type=args.update)
        res.excludeExistingLocalFiles()
        env.logger.info('{} files need to be downloaded or updated'.format(len(res.manifest)))
        res.downloadResources()
    elif args.upload:
        manager = ResourceManager()
        manager.getRemoteManifest('http://bioinformatics.mdanderson.org/Software/VariantTools/{}/'.format(args.repo))
        resource_dir = os.path.expanduser(args.local)
        # get information about file
        for filename in args.upload:
            if filename.endswith('.DB'):
                env.logger.info('Ignore uncompressed database file {}'.format(filename))
                continue
            rel_path = os.path.relpath(filename, resource_dir)
            if rel_path in manager.manifest:
                filesize = os.path.getsize(filename)
                md5 = calculateMD5(filename, partial=True)
                refGenome = manager.getRefGenome(filename)
                comment = manager.getComment(filename).replace('\n', ' ').replace('\t', ' ').strip()
                if (filesize, md5, refGenome, comment) == tuple(manager.manifest[rel_path][:4]):
                    env.logger.info('Ignoring identical file {}'.format(rel_path))
                    continue
            manager.addResource(filename)
            uploadFile(filename, rel_path, args.repo, args.user, args.passwd)
        manager.writeManifest('MANIFEST.tmp')
        uploadFile('MANIFEST.tmp', 'MANIFEST.txt', args.repo, args.user, args.passwd)
    elif args.remove:
        manager = ResourceManager()
        manager.getRemoteManifest('http://bioinformatics.mdanderson.org/Software/VariantTools/{}/'.format(args.repo))
        removed_count = 0
        for filename in args.remove:
            if filename in manager.manifest:
                manager.manifest.pop(filename)
                deprecateFile(filename, args.repo, args.user, args.passwd)
                env.logger.info('Remove {} from manifest'.format(filename))
                removed_count += 1
            else:
                env.logger.warning('{} does not exist in the manifest'.format(filename))
        # upload manifest
        if removed_count > 0:
            manager.writeManifest('MANIFEST.tmp')
            uploadFile('MANIFEST.tmp', 'MANIFEST.txt', args.repo, args.user, args.passwd)
    else:
        env.logger.warning('No option has been provided. Please use -h to get a list of actions.')
