#!/usr/bin/env python
#
# $File: liftOver.py $
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

import sys
import platform
import os
import stat
import subprocess
import tempfile
import shutil
    
from .project import Project
from .utils import ProgressBar, downloadFile, lineCount, getMaxUcscBin

#
class LiftOverTool:
    '''Calling USCS liftover tool to set alternative coordinate for all variants.
    '''
    def __init__(self, proj):
        self.proj = proj
        self.logger = proj.logger
        self.db = proj.db
    
    def obtainLiftOverTool(self):
        '''Obtain the liftOver tool, download from UCSC website if needed.'''
        try:
            subprocess.Popen(['liftOver'], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={'PATH':os.pathsep.join(['.', os.environ['PATH']])})
        except:
            # otherwise download the tool
            self.logger.debug('Failed to execute liftOver -h')
            if 'Darwin' in platform.system():
                # FIXME: Need to test PowerPC for directory macOSX.ppc
                liftOverDir = 'macOSX.i386'
            elif platform.system() == 'Linux':
                if platform.architecture()[0] == '64bit':
                    liftOverDir = 'linux.x86_64'
                else:
                    liftOverDir = 'linux.i386'
            else:
                self.logger.error('You platform does not support USCS liftOver tool. Please use a linux or MacOSX based machine.')
                self.logger.error('Optionally, you can compile liftOver for your platform and make it available to this script')
                return False
            #
            liftOverURL = 'http://hgdownload.cse.ucsc.edu/admin/exe/{0}/liftOver'.format(liftOverDir)
            try:
                self.logger.info('Downloading liftOver tool from UCSC')
                downloadFile(liftOverURL)
                os.chmod('liftOver', stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH) 
            except Exception as e:
                self.logger.warning('Failed to download UCSC liftOver tool from {0}'.format(liftOverURL))
                self.logger.warning('Please check the URL or manually download the file.')
                self.logger.debug(str(e))
                return False
        return True
 
    def obtainLiftOverChainFile(self, from_build, to_build):
        '''Obtain liftOver chain file, download from UCSC website if needed.'''
        # select or download the right chain file
        chainFile = '{0}To{1}.over.chain.gz'.format(from_build, to_build.title())
        if not os.path.isfile(chainFile):
            try:
                chainFileURL = 'http://hgdownload-test.cse.ucsc.edu/goldenPath/{0}/liftOver/{1}'.format(
                    from_build, chainFile)
                self.logger.info('Downloading liftOver chain file from UCSC')
                downloadFile(chainFileURL)
            except Exception as e:
                self.logger.warning('Failed to download chain file from {0}'.format(chainFileURL))
                self.logger.warning('Please check the URL, change --build and/or --alt_build, and try again')
                self.logger.warning('Optionally, you can download the right chain file and rename it to the one that is needed')
                self.logger.debug(e)
                return None
        return chainFile

    def exportVariantsInBedFormat(self, filename):
        '''Export variants in bed format'''
        cur = self.db.cursor()
        self.logger.info('Exporting variants in BED format');
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
                    if rec[1]:   # chr and pos of the primary coordinates can be None if they are back lifted from alternative reference genome
                        var_in.write('{0}\t{1}\t{2}\t{3}\n'.format(rec[1] if len(rec[1]) > 2 else 'chr' + rec[1],
                            int(rec[2]) - 1, rec[2], rec[0]))
                except Exception as e:
                    self.logger.debug('Invalid record {}'.format(rec))
                    self.logger.debug(e)
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        return count

    def runLiftOver(self, input, chain, output, unmapped):
        '''Executing UCUSC liftOver tool'''
        self.logger.info('Running UCSC liftOver tool')
        proc = subprocess.Popen(['liftOver', input, chain, output, unmapped],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                env={'PATH':os.pathsep.join(['.', os.environ['PATH']])})
        proc.wait()
        err = proc.stderr.read().decode().strip()
        if err:
            self.logger.info(err)

    def updateAltCoordinates(self):
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
        chainFile = self.obtainLiftOverChainFile(self.proj.build, self.proj.alt_build)
        if not chainFile:
            return {}
        # export existing variants to a temporary file
        tdir = tempfile.mkdtemp()
        num_variants = self.exportVariantsInBedFormat(os.path.join(tdir, 'var_in.bed'))
        if num_variants == 0:
            return {}
        self.runLiftOver(os.path.join(tdir, 'var_in.bed'), chainFile,
            os.path.join(tdir, 'var_out.bed'), os.path.join(tdir, 'unmapped.bed'))           
        #
        err_count = 0
        with open(os.path.join(tdir, 'unmapped.bed')) as var_err:
            for line in var_err:
                if line.startswith('#'):
                    continue
                if err_count == 0:
                    self.logger.debug('First 100 unmapped variants:')
                if err_count < 100:
                    self.logger.debug(line.rstrip())
                err_count += 1
        if err_count != 0:
            self.logger.info('{0} records failed to map.'.format(err_count))
        #
        # Update the master variant table.
        cur = self.db.cursor()
        headers = self.db.getHeaders('variant')
        if not 'alt_pos' in headers:
            self.logger.info('Adding fields alt_bin, alt_chr and alt_pos to table variant')
            self.db.execute('ALTER TABLE variant ADD alt_bin INT NULL;')
            self.db.execute('ALTER TABLE variant ADD alt_chr VARCHAR(20) NULL;')
            self.db.execute('ALTER TABLE variant ADD alt_pos INT NULL;')
        #
        mapped_file = os.path.join(tdir, 'var_out.bed')
        prog = ProgressBar('Updating table variant', lineCount(mapped_file))
        query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(self.db.PH)
        with open(mapped_file) as var_mapped:
            for count, line in enumerate(var_mapped):
                chr, start, end, id = line.strip().split()
                cur.execute(query, (getMaxUcscBin(int(start), int(end)), chr[3:] if chr.startswith('chr') else chr, int(start) + 1, id))
                if count % self.db.batch == 0:
                    self.db.commit()
                    prog.update(count)
        self.db.commit()
        prog.done()
        #
        # Method 1
        # create indexes separately
        #for field in ['alt_bin', 'alt_chr', 'alt_pos']:
        #    try:
        #        cur.execute('CREATE INDEX variant_{0}_idx ON variant ({0} ASC);'.format(field))
        #    except Exception as e:
        #        # maybe the index already exists
        #        self.logger.debug(e)
        #
        # Method 2
        # Create combined index
        try:
            cur.execute('CREATE INDEX varaint_{}_idx ON variant (alt_bin ASC, alt_chr ASC, alt_pos ASC, alt ASC);'.format(self.proj.alt_build))
        except Exception as e:
            # maybe the index already exists
            self.logger.debug(e)
        #
        self.db.commit()
        shutil.rmtree(tdir)
                
    def setAltRefGenome(self, alt_build):
        if self.proj.build == alt_build:
            raise ValueError('Cannot set alternative build the same as primary build')
        if self.proj.alt_build is not None and self.proj.alt_build != alt_build:
            self.logger.warning('Setting a different alternative reference genome.')
            self.logger.warning('The original alternative genome {} will be overritten.'.format(self.proj.alt_build))
        self.proj.alt_build = alt_build
        self.proj.saveProperty('alt_build', alt_build)
        self.updateAltCoordinates()

    def mapCoordinates(self, coordinates, from_build, to_build):
        '''Given a set of coordinates (chr, pos) in build, from and to build of reference genome,
        return a dictionary that maps (chr, pos) -> (alt_chr, alt_pos)
        '''
        if len(coordinates) == 0:
            return {}
        # download liftover 
        if not self.obtainLiftOverTool():
            raise RuntimeError('Failed to obtain UCSC LiftOver tool')
        # download chain file
        chainFile = self.obtainLiftOverChainFile(from_build, to_build)
        if not chainFile:
            raise RuntimeError('Failed to obtain UCSC chain file {}'.format(chainFile))
        # export existing variants to a temporary file
        tdir = tempfile.mkdtemp()
        with open(os.path.join(tdir, 'var_in.bed'), 'w') as output:
            for cor in sorted(coordinates):
                output.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(cor[0] if len(cor[0]) > 2 else 'chr' + cor[0],
                   cor[1] - 1, cor[1], cor[0], cor[1]))
        #
        self.runLiftOver(os.path.join(tdir, 'var_in.bed'), chainFile,
            os.path.join(tdir, 'var_out.bed'), os.path.join(tdir, 'unmapped.bed'))           
        #
        err_count = 0
        with open(os.path.join(tdir, 'unmapped.bed')) as var_err:
            for line in var_err:
                if line.startswith('#'):
                    continue
                if err_count == 0:
                    self.logger.debug('First 100 unmapped variants:')
                if err_count < 100:
                    self.logger.debug(line.rstrip())
                err_count += 1
        if err_count != 0:
            self.logger.info('{0} out of {1} records failed to map.'.format(err_count, len(coordinates)))
        #
        # create a map
        coordinateMap = {}
        mapped_file = os.path.join(tdir, 'var_out.bed')
        prog = ProgressBar('Reading new coordinates', lineCount(mapped_file))
        with open(mapped_file) as var_mapped:
            for count, line in enumerate(var_mapped.readlines()):
                try:
                    alt_chr, alt_start, alt_end, chr, pos = line.strip().split()
                    coordinateMap[(chr, int(pos))] = (alt_chr[3:] if alt_chr.startswith('chr') else alt_chr, int(alt_start) + 1)
                except Exception as e:
                    self.logger.debug(e)
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        return coordinateMap
    
#
#
# Functions provided by this script
#
#

def liftOverArguments(parser):
    parser.add_argument('build', help='Name of the alternative reference genome')

def liftOver(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            tool = LiftOverTool(proj)
            tool.setAltRefGenome(args.build)
        proj.close()
    except Exception as e:
        sys.exit(e)

