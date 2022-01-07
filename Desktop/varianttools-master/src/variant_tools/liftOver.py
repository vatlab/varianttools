#!/usr/bin/env python

#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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
import platform
import stat
import subprocess
import sys

from .project import Project
from .utils import (OS_ENV, ProgressBar, delayedAction, downloadFile, env,
                    getMaxUcscBin, lineCount)


#
class LiftOverTool:
    '''Calling USCS liftover tool to set alternative coordinate for all variants.
    '''

    def __init__(self, proj):
        self.proj = proj
        self.db = proj.db

    def obtainLiftOverTool(self):
        '''Obtain the liftOver tool, download from UCSC website if needed.'''
        try:
            subprocess.Popen(['liftOver'],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             env=OS_ENV)
        except:
            # otherwise download the tool
            env.logger.debug('Failed to execute liftOver -h')
            if 'Darwin' in platform.system():
                # FIXME: Need to test PowerPC for directory macOSX.ppc

                liftOverDir = 'macOSX.x86_64'
            elif platform.system() == 'Linux':
                if platform.architecture()[0] == '64bit':
                    liftOverDir = 'linux.x86_64'
                else:
                    liftOverDir = 'linux.i386'
            else:
                env.logger.error(
                    'You platform does not support USCS liftOver tool. Please use a linux or MacOSX based machine.'
                )
                env.logger.error(
                    'Optionally, you can compile liftOver for your platform and make it available to this script'
                )
                return False
            #
            liftOverURL = 'http://hgdownload.cse.ucsc.edu/admin/exe/{0}/liftOver'.format(
                liftOverDir)
            try:
                env.logger.info('Downloading liftOver tool from UCSC')
                liftOverExe = downloadFile(liftOverURL, env.cache_dir)
                os.chmod(
                    liftOverExe, stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR
                    | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH)
            except Exception as e:
                raise RuntimeError(
                    'Failed to download UCSC liftOver tool from {}: {}'.format(
                        liftOverURL, e))
        return True

    def obtainLiftOverChainFile(self, from_build, to_build):
        '''Obtain liftOver chain file, download from UCSC website if needed.'''
        # select or download the right chain file
        chainFile = '{0}To{1}.over.chain.gz'.format(from_build,
                                                    to_build.title())
        if not os.path.isfile(chainFile):
            try:
                chainFileURL = 'http://hgdownload.cse.ucsc.edu/goldenPath/{0}/liftOver/{1}'.format(
                    from_build, chainFile)
                env.logger.info('Downloading liftOver chain file from UCSC')
                chainFile = downloadFile(chainFileURL, env.cache_dir)
            except Exception as e:
                raise RuntimeError((
                    'Failed to download chain file from {0}'
                    'Please check the URL, change --build and/or --alt_build, and try again'
                    'Optionally, you can download the right chain file and rename it to the one that is needed: {1}'
                ).format(chainFileURL, e))
        return chainFile

    def exportVariantsInBedFormat(self, filename):
        '''Export variants in bed format'''
        cur = self.db.cursor()
        env.logger.info('Exporting variants in BED format')
        cur.execute('SELECT variant_id, chr, pos FROM variant;')
        # output ID so that we can re-insert
        prog = ProgressBar('Exporting variants', self.db.numOfRows("variant"))
        count = 0
        with open(filename, 'w') as var_in:
            for count, rec in enumerate(cur):
                # NOTE: the change from 1-based index to 0-based (assumed by liftOver)
                # NOTE: we add 'chr' to chromosome name (except for those long names) because
                #       liftover uses chr1, chr2 etc
                try:
                    if rec[1]:  # chr and pos of the primary coordinates can be None if they are back lifted from alternative reference genome
                        var_in.write('{0}\t{1}\t{2}\t{3}\n'.format(
                            rec[1] if len(rec[1]) > 2 else 'chr' + rec[1],
                            int(rec[2]) - 1, rec[2], rec[0]))
                except Exception as e:
                    env.logger.debug('Invalid record {}'.format(rec))
                    env.logger.debug(e)
                if count % 10000 == 0:
                    prog.update(count)
        prog.done()
        return count

    def runLiftOver(self, input, chain, output, unmapped):
        '''Executing UCUSC liftOver tool'''
        env.logger.info('Running UCSC liftOver tool')
        proc = subprocess.Popen(['liftOver', input, chain, output, unmapped],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                env=OS_ENV)
        ret = proc.wait()
        if ret != 0:
            err = proc.stderr.read().decode().strip()
            raise RuntimeError(
                'Failed to execute command liftover: {}'.format(err))

    def updateAltCoordinates(self, flip=False):
        '''Download and use the UCSC LiftOver tool to translate coordinates from primary
        to secondary reference genome'''
        if self.proj.build is None:
            raise ValueError('No data is available.')
        if self.proj.alt_build is None:
            raise ValueError('No alternative reference genome is specified.')
        # download liftover
        if not self.obtainLiftOverTool():
            return {}
        # download chain file
        chainFile = self.obtainLiftOverChainFile(self.proj.build,
                                                 self.proj.alt_build)
        if not chainFile:
            return {}
        # Update the master variant table.
        with delayedAction(
                env.logger.info,
                'Adding alternative reference genome {} to the project.'.format(
                    self.proj.alt_build)):
            cur = self.db.cursor()
            headers = self.db.getHeaders('variant')
            for fldName, fldType in [('alt_bin', 'INT'),
                                     ('alt_chr', 'VARCHAR(20)'),
                                     ('alt_pos', 'INT')]:
                if fldName in headers:
                    continue
                self.db.execute('ALTER TABLE variant ADD {} {} NULL;'.format(
                    fldName, fldType))
        # export existing variants to a temporary file
        num_variants = self.exportVariantsInBedFormat(
            os.path.join(env.temp_dir, 'var_in.bed'))
        if num_variants == 0:
            return {}
        self.runLiftOver(
            os.path.join(env.temp_dir, 'var_in.bed'), chainFile,
            os.path.join(env.temp_dir, 'var_out.bed'),
            os.path.join(env.temp_dir, 'unmapped.bed'))
        #
        err_count = 0
        unmapped_file = os.path.join(env.temp_dir, 'unmapped.bed')
        if os.path.isfile(unmapped_file):
            with open(unmapped_file) as var_err:
                for line in var_err:
                    if line.startswith('#'):
                        continue
                    if err_count == 0:
                        env.logger.debug('First 100 unmapped variants:')
                    if err_count < 100:
                        env.logger.debug(line.rstrip())
                    err_count += 1
        if err_count != 0:
            env.logger.info('{0} records failed to map.'.format(err_count))
        #
        #
        mapped_file = os.path.join(env.temp_dir, 'var_out.bed')
        if flip:
            env.logger.info('Flipping primary and alternative reference genome')
            cur.execute(
                'UPDATE variant SET alt_bin=bin, alt_chr=chr, alt_pos=pos;')
            cur.execute('UPDATE variant SET bin=NULL, chr=NULL, pos=NULL')
            query = 'UPDATE variant SET bin={0}, chr={0}, pos={0} WHERE variant_id={0};'.format(
                self.db.PH)
        else:
            query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(
                self.db.PH)
        prog = ProgressBar('Updating table variant', lineCount(mapped_file))
        with open(mapped_file) as var_mapped:
            for count, line in enumerate(var_mapped):
                chr, start, end, id = line.strip().split()
                cur.execute(query, (getMaxUcscBin(int(start), int(end)),
                                    chr[3:] if chr.startswith('chr') else chr,
                                    int(start) + 1, id))
                if count % 10000 == 0:
                    self.db.commit()
                    prog.update(count)
        self.db.commit()
        prog.done()

    def setAltRefGenome(self, alt_build, build_index=True, flip=False):
        if self.proj.build == alt_build:
            raise ValueError(
                'Cannot set alternative build the same as primary build')
        if self.proj.alt_build is not None and self.proj.alt_build != alt_build:
            env.logger.warning(
                'Setting a different alternative reference genome.')
            env.logger.warning(
                'The original alternative genome {} will be overwritten.'
                .format(self.proj.alt_build))
        self.proj.alt_build = alt_build
        self.proj.saveProperty('alt_build', alt_build)
        self.updateAltCoordinates(flip)
        if flip:
            self.proj.alt_build, self.proj.build = self.proj.build, self.proj.alt_build
            self.proj.saveProperty('build', self.proj.build)
            self.proj.saveProperty('alt_build', self.proj.alt_build)
            self.proj.dropIndexOnMasterVariantTable()
        if build_index:
            self.proj.createIndexOnMasterVariantTable()

    def mapCoordinates(self, map_in, from_build, to_build):
        '''Given a input file, run liftover and return output file.
        '''
        # download liftover
        if not self.obtainLiftOverTool():
            raise RuntimeError('Failed to obtain UCSC LiftOver tool')
        # download chain file
        chainFile = self.obtainLiftOverChainFile(from_build, to_build)
        if not chainFile:
            raise RuntimeError(
                'Failed to obtain UCSC chain file {}'.format(chainFile))
        # export existing variants to a temporary file
        self.runLiftOver(map_in, chainFile,
                         os.path.join(env.temp_dir, 'var_out.bed'),
                         os.path.join(env.temp_dir, 'unmapped.bed'))
        #
        err_count = 0
        with open(os.path.join(env.temp_dir, 'unmapped.bed')) as var_err:
            for line in var_err:
                if line.startswith('#'):
                    continue
                if err_count == 0:
                    env.logger.debug('First 100 unmapped variants:')
                if err_count < 100:
                    env.logger.debug(line.rstrip())
                err_count += 1
        #
        return os.path.join(env.temp_dir, 'var_out.bed'), err_count


#
#
# Functions provided by this script
#
#


def liftOverArguments(parser):
    parser.add_argument(
        'build', help='Name of the alternative reference genome'),
    parser.add_argument(
        '--flip',
        action='store_true',
        help='''Flip primary and alternative
	reference genomes so that the specified build will become the primary
        reference genome of the project.''')


def liftOver(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            tool = LiftOverTool(proj)
            tool.setAltRefGenome(args.build, flip=args.flip)
        proj.close()
    except Exception as e:
        sys.exit(e)
