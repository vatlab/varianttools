#!/usr/bin/env python
#
# $File: pipeline.py$
# $LastChangedDate: 2013-04-23 11:58:41 -0500 (Tue, 23 Apr 2013) $
# $Rev: 1855 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2013 Bo Peng (bpeng@mdanderson.org)
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
import subprocess
import glob
import argparse
import logging
import shutil
import tarfile
import copy
import gzip
import bz2
import zipfile
import time
import re
import platform
from collections import namedtuple

from .utils import env, ProgressBar, downloadFile, downloadURL, \
    existAndNewerThan, elapsed_time, run_command, wait_all, TEMP
    
from .project import PipelineDescription, Project

try:
    import pysam
    hasPySam = True
except ImportError as e:
    hasPySam = False


class SequentialAction:
    def __init__(self, actions):
        '''Define an action that calls a list of actions, specified
        by Action1, Action2.'''
        self.actions = []
        for a in actions:
            if hasattr(a, '__call__'):
                self.actions.append(a.__call__)
            else:
                self.actions.append(a)

    def __call__(self, input_files, output_files):
        for a in self.actions:
            if a(input_files, output_files):
                raise RuntimeError('Execution of action failed')
        return 0

class RunCommand:
    def __init__(self, cmd, working_dir=None):
        self.cmd = cmd
        self.working_dir = working_dir

    def __call__(self, input_files, output_files):
        # substitute cmd by input_files and output_files
        pass

class CheckCommand:
    def __init__(self, cmd):
        self.cmd = cmd

    def __call__(self, input_files, output_files):
        if not hasattr(shutil, 'which'):
            # if shutil.which does not exist, use subprocess...
            try:
                subprocess.call(cmd, stdout=open(os.devnull, 'w'), stderr=open(os.devnull, 'w'))
                return 0
            except:
                env.logger.error('Command {} does not exist. Please install it and try again.'
                    .format(cmd))
                return 1
        if shutil.which(cmd) is None:
            env.logger.error('Command {} does not exist. Please install it and try again.'
                .format(cmd))
            return 1

class GetFastqFiles:
    def __init__(self):
        pass

    def __call__(self, input_files, output_files):
        # decompress input files and return a list of output files
        filenames = []
        for filename in input_files:
            if filename.lower().endswith('.bam') or filename.lower().endswith('.sam'):
                filenames.extend(self.picard.bam2fastq(filename))
                continue
            for fastq_file in decompress(filename, env.WORKING_DIR):
                try:
                    with open(fastq_file) as fastq:
                        line = fastq.readline()
                        if not line.startswith('@'):
                            raise ValueError('Wrong FASTA file {}'.foramt(fastq_file))
                    filenames.append(fastq_file)
                except Exception as e:
                    env.logger.error('Ignoring non-fastq file {}: {}'
                        .format(fastq_file, e))
        filenames.sort()
        return filenames        


class Pipeline:
    def __init__(self, name):
        self.pipeline = PipelineDescription(name)
        #
        # resource directory
        #
        self.pipeline_resource = os.path.join(os.path.expanduser(
            env.local_resource), 'var_caller', self.pipeline.name)
        try:
            if not os.path.isdir(self.pipeline_resource):
                sys.stderr.write('Creating pipeline resource directory {}\n'
                    .format(self.pipeline_resource))
                os.makedirs(self.pipeline_resource)
        except:
            raise RuntimeError('Failed to create pipeline resource directory '
                .format(self.pipeline_resource))

    def downloadResource(self):
        '''Download resource'''
        # decompress all .gz files
        saved_dir = os.getcwd()
        os.chdir(self.pipeline_resource)
        skipped = []
        for cnt, URL in enumerate(sorted(self.pipeline.resource)):
            filename = URL.rsplit('/', 1)[-1]
            dest_file = os.path.join(self.pipeline_resource, filename)
            try:
                if os.path.isfile(dest_file):
                    skipped.append(filename)
                else:
                    downloadURL(URL, dest_file, False,
                        message='{}/{} {}'.format(cnt+1, len(self.pipeline.resource), filename))
            except KeyboardInterrupt as e:
                raise e
            except Exception as e:
                env.logger.error('Failed to download {}: {} {}'
                    .format(filename, type(e).__name__, e))
            #
            if filename.endswith('.gz') and not filename.endswith('tar.gz'):
                if not existAndNewerThan(ofiles=filename[:-3], ifiles=filename):
                    decompress(filename, '.')
        os.chdir(saved_dir)
        if skipped:
            env.logger.info('Using {} existing resource files under {}.'
                .format(len(skipped), self.pipeline_resource))
 
    def executeCommand(self):
        pass

    def execute(self, steps, input_files=[], output_files=[]):
        for command in {'init': self.pipeline.init_steps,
                'align': self.pipeline.align_steps,
                'call': self.pipeline.call_steps}[steps]:
            print command.comment
            

###############################################################################
#
# Utility functions
#
###############################################################################


def downloadFile(URL, dest, quiet=False):
    '''Download a file from URL and save to dest'''
    # for some strange reason, passing wget without shell=True can fail silently.
    env.logger.info('Downloading {}'.format(URL))
    if os.path.isfile(dest):
        env.logger.warning('Using existing downloaded file {}.'.format(dest))
        return dest
    p = subprocess.Popen('wget {} -O {} {}'
        .format('-q' if quiet else '', TEMP(dest), URL), shell=True)
    ret = p.wait()
    if ret == 0 and os.path.isfile(TEMP(dest)):
        os.rename(TEMP(dest), dest)
        return dest
    else:
        try:
            os.remove(TEMP(dest))
        except OSError:
            pass
        raise RuntimeError('Failed to download {} using wget'.format(URL))


def fastqVersion(fastq_file):
    '''Detect the version of input fastq file. This can be very inaccurate'''
    #
    # This function assumes each read take 4 lines, and the last line contains
    # quality code. It collects about 1000 quality code and check their range,
    # and use it to determine if it is Illumina 1.3+
    #
    qual_scores = ''
    with open(fastq_file) as fastq:
        while len(qual_scores) < 1000:
            try:
                line = fastq.readline()
            except Exception as e:
                env.logger.error('Failed to read fastq file {}: {}'
                    .format(fastq_file, e))
                sys.exit(1)
            if not line.startswith('@'):
                raise ValueError('Wrong FASTA file {}'.format(fastq_file))
            line = fastq.readline()
            line = fastq.readline()
            if not line.startswith('+'):
                env.logger.warning(
                    'Suspiciout FASTA file {}: third line does not start with "+".'
                    .foramt(fastq_file))
                return 'Unknown'
            line = fastq.readline()
            qual_scores += line.strip()
    #
    min_qual = min([ord(x) for x in qual_scores])
    max_qual = max([ord(x) for x in qual_scores])
    env.logger.debug('FASTA file with quality score ranging {} to {}'
        .format(min_qual, max_qual))
    # Sanger qual score has range Phred+33, so 33, 73 with typical score range 0 - 40
    # Illumina qual scores has range Phred+64, which is 64 - 104 with typical score range 0 - 40
    if min_qual >= 64 or max_qual > 90:
        # option -I is needed for bwa if the input is Illumina 1.3+ read format (quliaty equals ASCII-64).
        return 'Illumina 1.3+'
    else:
        # no option is needed for bwa
        return 'Sanger'

def decompress(filename, dest_dir=None):
    '''If the file ends in .tar.gz, .tar.bz2, .bz2, .gz, .tgz, .tbz2, decompress
    it to dest_dir (current directory if unspecified), and return a list of files.
    Uncompressed files will be returned untouched. If the destination files exist
    and newer, this function will return immediately.'''
    mode = None
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
        mode = 'r:gz'
    elif filename.lower().endswith('.tbz2') or filename.lower().endswith('.tgz'):
        mode = 'r:bz2'
    elif filename.lower().endswith('.tar'):
        mode = 'r'
    elif filename.lower().endswith('.gz'):
        dest_file = os.path.join('.' if dest_dir is None else dest_dir,
            os.path.basename(filename)[:-3])
        if existAndNewerThan(ofiles=dest_file, ifiles=filename):
            env.logger.warning('Using existing decompressed file {}'.format(dest_file))
        else:
            env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
            with gzip.open(filename, 'rb') as gzinput, open(TEMP(dest_file), 'wb') as output:
                content = gzinput.read(10000000)
                while content:
                    output.write(content)
                    content = gzinput.read(10000000)
            # only rename the temporary file to the right one after finishing everything
            # this avoids corrupted files
            os.rename(TEMP(dest_file), dest_file)
        return [dest_file]
    elif filename.lower().endswith('.bz2'):
        dest_file = os.path.join('.' if dest_dir is None else dest_dir, os.path.basename(filename)[:-4])
        if existAndNewerThan(ofiles=dest_file, ifiles=filename):
            env.logger.warning('Using existing decompressed file {}'.format(dest_file))
        else:
            env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
            with bz2.open(filename, 'rb') as bzinput, open(TEMP(dest_file), 'wb') as output:
                content = bzinput.read(10000000)
                while content:
                    output.write(content)
                    content = bzinput.read(10000000)
            # only rename the temporary file to the right one after finishing everything
            # this avoids corrupted files
            os.rename(TEMP(dest_file), dest_file)
        return [dest_file]
    elif filename.lower().endswith('.zip'):
        bundle = zipfile.ZipFile(filename)
        dest_dir = '.' if dest_dir is None else dest_dir
        bundle.extractall(dest_dir)
        env.logger.info('Decompressing {} to {}'.format(filename, dest_dir))
        return [os.path.join(dest_dir, name) for name in bundle.namelist()]
    #
    # if it is a tar file
    if mode is not None:
        env.logger.info('Extracting fastq sequences from tar file {}'
            .format(filename))
        #
        # MOTE: open a compressed tar file can take a long time because it needs to scan
        # the whole file to determine its content. I am therefore creating a manifest
        # file for the tar file in the dest_dir, and avoid re-opening when the tar file
        # is processed again.
        manifest = os.path.join( '.' if dest_dir is None else dest_dir,
            os.path.basename(filename) + '.manifest')
        all_extracted = False
        dest_files = []
        if existAndNewerThan(ofiles=manifest, ifiles=filename):
            all_extracted = True
            for f in [x.strip() for x in open(manifest).readlines()]:
                dest_file = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(f))
                if existAndNewerThan(ofiles=dest_file, ifiles=filename):
                    dest_files.append(dest_file)
                    env.logger.warning('Using existing extracted file {}'.format(dest_file))
                else:
                    all_extracted = False
        #
        if all_extracted:
            return dest_files
        #
        # create a temporary directory to avoid corrupted file due to interrupted decompress
        try:
            os.mkdir('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
        except:
            # directory might already exist
            pass
        #
        dest_files = []
        with tarfile.open(filename, mode) as tar:
            files = tar.getnames()
            # save content to a manifest
            with open(manifest, 'w') as manifest:
                for f in files:
                    manifest.write(f + '\n')
            for f in files:
                # if there is directory structure within tar file, decompress all to the current directory
                dest_file = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(f))
                dest_files.append(dest_file)
                if existAndNewerThan(ofiles=dest_file, ifiles=filename):
                    env.logger.warning('Using existing extracted file {}'.format(dest_file))
                else:
                    env.logger.info('Extracting {} to {}'.format(f, dest_file))
                    tar.extract(f, 'tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
                    # move to the top directory with the right name only after the file has been properly extracted
                    shutil.move(os.path.join('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'), f), dest_file)
            # set dest_files to the same modification time. This is used to
            # mark the right time when the files are created and avoid the use
            # of archieved but should-not-be-used files that might be generated later
            [os.utime(x, None) for x in dest_files]
        return dest_files
    # return source file if 
    return [filename]
   
   
def getReadGroup(fastq_filename, output_bam):
    '''Get read group information from names of fastq files.'''
    # Extract read group information from filename such as
    # GERALD_18-09-2011_p-illumina.8_s_8_1_sequence.txt. The files are named 
    # according to the lane that produced them and whether they
    # are paired or not: Single-end reads s_1_sequence.txt for lane 1;
    # s_2_sequence.txt for lane 2 Paired-end reads s_1_1_sequence.txt 
    # for lane 1, pair 1; s_1_2_sequence.txt for lane 1, pair 2
    #
    # This function return a read group string like '@RG\tID:foo\tSM:bar'
    #
    # ID* Read group identifier. Each @RG line must have a unique ID. The
    # value of ID is used in the RG
    #     tags of alignment records. Must be unique among all read groups
    #     in header section. Read group
    #     IDs may be modifid when merging SAM fies in order to handle collisions.
    # CN Name of sequencing center producing the read.
    # DS Description.
    # DT Date the run was produced (ISO8601 date or date/time).
    # FO Flow order. The array of nucleotide bases that correspond to the
    #     nucleotides used for each
    #     flow of each read. Multi-base flows are encoded in IUPAC format, 
    #     and non-nucleotide flows by
    #     various other characters. Format: /\*|[ACMGRSVTWYHKDBN]+/
    # KS The array of nucleotide bases that correspond to the key sequence
    #     of each read.
    # LB Library.
    # PG Programs used for processing the read group.
    # PI Predicted median insert size.
    # PL Platform/technology used to produce the reads. Valid values: 
    #     CAPILLARY, LS454, ILLUMINA,
    #     SOLID, HELICOS, IONTORRENT and PACBIO.
    # PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for
    #     SOLiD). Unique identifier.
    # SM Sample. Use pool name where a pool is being sequenced.
    #
    filename = os.path.basename(fastq_filename)
    output = os.path.basename(output_bam)
    # sample name is obtained from output filename without file extension
    SM = output.split('.', 1)[0]
    # always assume ILLUMINA for this script and BWA for processing
    PL = 'ILLUMINA'  
    PG = 'BWA'
    #
    # PU is for flowcell and lane information, ID should be unique for each
    #     readgroup
    # ID is temporarily obtained from input filename without exteion
    ID = filename.split('.')[0]
    # try to get lan information from s_x_1/2 pattern
    try:
        PU = re.search('s_([^_]+)_', filename).group(1)
    except AttributeError:
        env.logger.warning('Failed to guess lane information from filename {}'
            .format(filename))
        PU = 'NA'
    # try to get some better ID
    try:
        # match GERALD_18-09-2011_p-illumina.8_s_8_1_sequence.txt
        m = re.match('([^_]+)_([^_]+)_([^_]+)_s_([^_]+)_([^_]+)_sequence.txt', filename)
        ID = '{}.{}'.format(m.group(1), m.group(4))
    except AttributeError as e:
        env.logger.warning('Input fasta filename {} does not match a known'
            ' pattern. ID is directly obtained from filename.'.format(filename))
    #
    rg = r'@RG\tID:{}\tPG:{}\tPL:{}\tPU:{}\tSM:{}'.format(ID, PG, PL, PU, SM)
    env.logger.info('Setting read group tag to {}'.format(rg))
    return rg


def countUnmappedReads(sam_files):
    #
    # count total reads and unmapped reads
    #
    # The previous implementation uses grep and wc -l, but
    # I cannot understand why these commands are so slow...
    #
    targets = ['{}.counts'.format(x) for x in sam_files]
    if not existAndNewerThan(ofiles=targets, ifiles=sam_files):
        for sam_file, target_file in zip(sam_files, targets):
            env.logger.info('Counting unmapped reads in {}'.format(sam_file))
            unmapped_count = 0
            total_count = 0
            with open(sam_file) as sam:
               for line in sam:
                   total_count += 1
                   if 'XT:A:N' in line:
                       unmapped_count += 1
            with open(target_file, 'w') as target:
                target.write('{}\n{}\n'.format(unmapped_count, total_count))
    #
    counts = []
    for count_file in targets:
        with open(count_file) as cnt:
            unmapped = int(cnt.readline())
            total = int(cnt.readline())
            counts.append((unmapped, total))
    return counts

def isBamPairedEnd(input_file):
    # we need pysam for this but pysam does not yet work for Python 3.3.
    if not hasPySam:
        env.logger.error('Cannot detect if input bam file has paired reads (missing pysam). Assuming paired.')
        return True
    bamfile = pysam.Samfile(input_file, 'rb')
    for read in bamfile.fetch():
        if read.is_paired:
            env.logger.info('{} is detected to have paired reads.'.format(input_file))
        else:
            env.logger.info('{} is detected to have single (non-paired) reads.'.format(input_file))
        return read.is_paired


###############################################################################
#
# Wrappers for commands used. Note that
# 
# 1. For each command, it should be initialized with basic information
#    for it to execute (e.g. reference genome) so that they do not have to
#    be supplied later.
#
# 2. The __init__ function should check the existence of the file, even
#    the version of the file, so the initialization of pipeline will check
#    the existence of the commands.
#
# 3. The member functions should wrap around functions of commands. These
#    functions should
#
#    a. accept a list of input files, and essential parameter
#    b. write to output file if specified from parameter (these are generally
#       the last step of the pipeline where the location of results are
#       specified by users.
#    c. these function should make use of existAndNewer and quit if the result
#       files exist and are newer.
#    d. these functions should try to use wait=True to allow parallelized
#       execution.
#    e. these functions should make use of env.WORKING_DIR as temporary output
#    f. extra options should be provided by users through env.options.
#
###############################################################################
class BWA:
    def __init__(self, ref_fasta, version=None):
        '''A BWA wrapper, if version is not None, a specific version is
        required.'''
        checkCmd('bwa')
        self.REF_fasta = ref_fasta
        
    def buildIndex(self):
        '''Create BWA index for reference genome file'''
        # bwa index  -a bwtsw wg.fa
        saved_dir = os.getcwd()
        os.chdir(env.RESOURCE_DIR)
        if os.path.isfile(self.REF_fasta + '.amb'):
            env.logger.warning('Using existing bwa indexed sequence {}.amb'
                .format(self.REF_fasta))
        else:
            run_command('bwa index {}  -a bwtsw {}'
                .format(env.OPT_BWA_INDEX, self.REF_fasta))
        os.chdir(saved_dir)

    def aln(self, fastq_files):
        '''Use bwa aln to process fastq files'''
        for input_file in fastq_files:
            dest_file = '{}/{}.sai'.format(env.WORKING_DIR, os.path.basename(input_file))
            if existAndNewerThan(ofiles=dest_file, ifiles=input_file):
                env.logger.warning('Using existing alignment index file {}'
                    .format(dest_file))
            else:
                # input file should be in fastq format (-t 4 means 4 threads)
                opt = ' -I ' if fastqVersion(input_file) == 'Illumina 1.3+' else ''
                if opt == ' -I ':
                    env.logger.warning('Using -I option for bwa aln command '
                        'because the sequences seem to be in Illumina 1.3+ format.')
                run_command('bwa aln {} {} -t 4 {}/{} {} > {}'
                    .format(opt, env.OPT_BWA_ALN, env.RESOURCE_DIR, 
                        self.REF_fasta, input_file, TEMP(dest_file)),
                    name=os.path.basename(dest_file),
                    upon_succ=(os.rename, TEMP(dest_file), dest_file),
                    wait=False)
        # wait for all bwa aln jobs to be completed
        wait_all()

    def sampe(self, fastq_files, sample_name=None):
        '''Use bwa sampe to generate aligned sam files for paird end reads'''
        sam_files = []
        for idx in range(len(fastq_files)//2):
            f1 = fastq_files[2*idx]
            f2 = fastq_files[2*idx + 1]
            rg = getReadGroup(f1, sample_name)
            sai1 = '{}/{}.sai'.format(env.WORKING_DIR, os.path.basename(f1))
            sai2 = '{}/{}.sai'.format(env.WORKING_DIR, os.path.basename(f2))
            sam_file = '{}/{}_bwa.sam'.format(env.WORKING_DIR, os.path.basename(f1))
            if existAndNewerThan(ofiles=sam_file, ifiles=[f1, f2, sai1, sai2]):
                env.logger.warning('Using existing sam file {}'.format(sam_file))
            else:
                run_command(
                    'bwa sampe {0} -r \'{1}\' {2}/{3} {4} {5} {6} {7} > {8}'
                    .format(
                        env.OPT_BWA_SAMPE, rg, env.RESOURCE_DIR, self.REF_fasta,
                        sai1, sai2, f1, f2, TEMP(sam_file)),
                    name=os.path.basename(sam_file),
                    upon_succ=(os.rename, TEMP(sam_file), sam_file),
                    wait=False)
            sam_files.append(sam_file)
        # wait for all jobs to be completed
        wait_all()
        return sam_files

    def samse(self, fastq_files, sample_name=None):
        '''Use bwa sampe to generate aligned sam files'''
        sam_files = []
        for f in fastq_files:
            sam_file = '{}/{}_bwa.sam'.format(env.WORKING_DIR, os.path.basename(f))
            rg = getReadGroup(f, sample_name)
            sai = '{}/{}.sai'.format(env.WORKING_DIR, os.path.basename(f))
            if existAndNewerThan(ofiles=sam_file, ifiles=[f, sai]):
                env.logger.warning('Using existing sam file {}'.format(sam_file))
            else:
                run_command(
                    'bwa samse {0} -r \'{1}\' {2}/{3} {4} {5} > {6}'
                    .format(
                        env.OPT_BWA_SAMSE, rg, env.RESOURCE_DIR,
                        self.REF_fasta, sai, f, TEMP(sam_file)),
                    name=os.path.basename(sam_file),
                    upon_succ=(os.rename, TEMP(sam_file), sam_file),
                    wait=False)
            sam_files.append(sam_file)
        # wait for all jobs to be completed
        wait_all()
        return sam_files        


class SAMTOOLS:
    def __init__(self, ref_fasta):
        checkCmd('samtools')
        self.REF_fasta = ref_fasta

    def buildIndex(self):
        '''Create index for reference genome used by samtools'''
        saved_dir = os.getcwd()
        os.chdir(env.RESOURCE_DIR)
        if os.path.isfile('{}.fai'.format(self.REF_fasta)):
            env.logger.warning('Using existing samtools sequence index {}.fai'
                .format(self.REF_fasta))
        else:
            run_command('samtools faidx {} {}'
                .format(env.OPT_SAMTOOLS_FAIDX, self.REF_fasta))
        os.chdir(saved_dir)

    def sortSam(self, sam_files):
        '''Convert sam file to sorted bam files.'''
        bam_files = []
        for sam_file in sam_files:
            bam_file = sam_file[:-4] + '.bam'
            if existAndNewerThan(ofiles=bam_file, ifiles=sam_file):
                env.logger.warning('Using existing bam file {}'.format(bam_file))
            else:
                run_command('samtools view {} -bt {}/{}.fai {} > {}'
                    .format(
                        env.OPT_SAMTOOLS_VIEW, env.RESOURCE_DIR,
                        self.REF_fasta, sam_file, TEMP(bam_file)),
                    name=os.path.basename(bam_file),
                    upon_succ=(os.rename, TEMP(bam_file), bam_file),
                    wait=False)
            bam_files.append(bam_file)
        # wait for all sam->bam jobs to be completed
        wait_all()
        #
        # sort bam files
        sorted_bam_files = []
        for bam_file in bam_files:
            sorted_bam_file = bam_file[:-4] + '_sorted.bam'
            if existAndNewerThan(ofiles=sorted_bam_file, ifiles=bam_file):
                env.logger.warning('Using existing sorted bam file {}'
                    .format(sorted_bam_file))
            else:
                run_command('samtools sort {} {} {}'
                    .format(
                        env.OPT_SAMTOOLS_SORT, bam_file,
                        TEMP(sorted_bam_file)),
                    name=os.path.basename(sorted_bam_file),
                    upon_succ=(os.rename, TEMP(sorted_bam_file), sorted_bam_file),
                    wait=False)
            sorted_bam_files.append(sorted_bam_file)
        wait_all()
        return sorted_bam_files

    def indexBAM(self, bam_file):
        '''Index the input bam file'''
        if existAndNewerThan(ofiles='{}.bai'.format(bam_file), ifiles=bam_file):
            env.logger.warning('Using existing bam index {}.bai'.format(bam_file))
        else:
            # samtools index input.bam [output.bam.bai]
            run_command('samtools index {} {} {}'.format(
                env.OPT_SAMTOOLS_INDEX, bam_file, TEMP(bam_file + '.bai')),
                upon_succ=(os.rename, TEMP(bam_file + '.bai'), bam_file + '.bai'))


class PICARD:
    '''Wrapper around command picard'''
    def __init__(self):
        if env.PICARD_PATH:
            if not os.path.isfile(os.path.join(os.path.expanduser(env.PICARD_PATH), 'SortSam.jar')):
                env.logger.error('Specified PICARD_PATH {} does not contain picard jar files.'
                    .format(env.PICARD_PATH))
                sys.exit(1)
        elif 'CLASSPATH' in os.environ:
            if not any([os.path.isfile(os.path.join(os.path.expanduser(x), 'SortSam.jar')) 
                    for x in os.environ['CLASSPATH'].split(':')]):
                env.logger.error('$CLASSPATH ({}) does not contain a path that contain picard jar files.'
                    .format(os.environ['CLASSPATH']))
                sys.exit(1)
            for x in os.environ['CLASSPATH'].split(':'):
                if os.path.isfile(os.path.join(os.path.expanduser(x), 'SortSam.jar')):
                    env.logger.info('Using picard under {}'.format(x))
                    env.PICARD_PATH = os.path.expanduser(x)
                    break
        else:
            env.logger.error('Please either specify path to picard using option '
                '--set PICARD_PATH=path, or set it in environment variable $CLASSPATH.')
            sys.exit(1)

    def sortSam(self, sam_files, output=None):
        '''Convert sam file to sorted bam files using Picard.'''
        # sort bam files
        sorted_bam_files = []
        for sam_file in sam_files:
            sorted_bam_file = sam_file[:-4] + '_sorted.bam'
            if existAndNewerThan(ofiles=sorted_bam_file, ifiles=sam_file):
                env.logger.warning('Using existing sorted bam file {}'
                    .format(sorted_bam_file))
            else:
                run_command(
                    'java {0} -jar {1}/SortSam.jar {2} I={3} O={4} SO=coordinate'
                    .format(
                        env.OPT_JAVA, env.PICARD_PATH, 
                        env.OPT_PICARD_SORTSAM, sam_file,
                        TEMP(sorted_bam_file)),
                    name=os.path.basename(sorted_bam_file),
                    upon_succ=(os.rename, TEMP(sorted_bam_file), sorted_bam_file),
                    wait=False)
            sorted_bam_files.append(sorted_bam_file)
        wait_all()
        return sorted_bam_files

    def mergeBAMs(self, bam_files):
        '''merge sam files'''
        # use Picard merge, not samtools merge: 
        # Picard keeps RG information from all Bam files, whereas samtools uses only 
        # inf from the first bam file
        merged_bam_file = bam_files[0][:-4] + '_merged.bam'
        if existAndNewerThan(ofiles=merged_bam_file, ifiles=bam_files):
            env.logger.warning('Using existing merged bam file {}'
                .format(merged_bam_file))
        else:
            run_command('''java {} -jar {}/MergeSamFiles.jar {} {}
                USE_THREADING=true
                VALIDATION_STRINGENCY=LENIENT
                ASSUME_SORTED=true
                OUTPUT={}'''.format(
                    env.OPT_JAVA, env.PICARD_PATH,
                    env.OPT_PICARD_MERGESAMFILES,
                    ' '.join(['INPUT={}'.format(x) for x in bam_files]),
                    TEMP(merged_bam_file)),
                name=os.path.basename(merged_bam_file),
                upon_succ=(os.rename, TEMP(merged_bam_file), merged_bam_file))
        return merged_bam_file

    def bam2fastq(self, input_file):
        '''This function extracts raw reads from an input BAM file to one or 
        more fastq files.'''
        #
        if isBamPairedEnd(input_file):
            output_files = [os.path.join(env.WORKING_DIR, '{}_{}.fastq'
                .format(os.path.basename(input_file)[:-4], x)) for x in [1,2]]
            if existAndNewerThan(ofiles=output_files, ifiles=input_file):
                env.logger.warning('Using existing sequence files {}'
                    .format(' and '.join(output_files)))
            else:
                run_command('''java {} -jar {}/SamToFastq.jar {} INPUT={}
                    FASTQ={} SECOND_END_FASTQ={}'''
                    .format(env.OPT_JAVA, env.PICARD_PATH,
                        env.OPT_PICARD_SAMTOFASTQ, input_file,
                        TEMP(output_files[0]), output_files[1]),
                    name=os.path.basename(output_files[0]),
                    upon_succ=(os.rename, TEMP(output_files[0]), output_files[0]))
            return output_files
        else:
            output_file = os.path.join(env.WORKING_DIR, '{}.fastq'
                .format(os.path.basename(input_file)[:-4]))
            if existAndNewerThan(ofiles=output_file, ifiles=input_file):
                env.logger.warning('Using existing sequence files {}'
                    .format(output_file))
            else:
                run_command('''java {} -jar {}/SamToFastq.jar {} INPUT={}
                    FASTQ={}'''
                    .format(env.OPT_JAVA, env.PICARD_PATH,
                        env.OPT_PICARD_SAMTOFASTQ, input_file,
                        TEMP(output_file)),
                    name=os.path.basename(output_file),
                    upon_succ=(os.rename, TEMP(output_file), output_file))
            return [output_file]

    def markDuplicates(self, bam_files):
        '''Mark duplicate using picard'''
        dedup_bam_files = []
        for bam_file in bam_files:
            dedup_bam_file = os.path.join(env.WORKING_DIR, os.path.basename(bam_file)[:-4] + '.dedup.bam')
            metrics_file = os.path.join(env.WORKING_DIR, os.path.basename(bam_file)[:-4] + '.metrics')
            if existAndNewerThan(ofiles=dedup_bam_file, ifiles=bam_file):
                env.logger.warning(
                    'Using existing bam files after marking duplicate {}'
                    .format(dedup_bam_file))
            else:
                run_command('''java {0} -jar {1}/MarkDuplicates.jar {2}
                    INPUT={3}
                    OUTPUT={4}
                    METRICS_FILE={5}
                    ASSUME_SORTED=true
                    VALIDATION_STRINGENCY=LENIENT
                    '''.format(
                        env.OPT_JAVA, env.PICARD_PATH,
                        env.OPT_PICARD_MARKDUPLICATES, bam_file,
                        TEMP(dedup_bam_file),
                        metrics_file), 
                    name=os.path.basename(dedup_bam_file),
                    upon_succ=(os.rename, TEMP(dedup_bam_file), dedup_bam_file),
                    wait=False)
            dedup_bam_files.append(dedup_bam_file)
        wait_all()
        return dedup_bam_files


class GATK:
    '''Wrapper around command gatk '''
    def __init__(self, ref_fasta, knownSites):
        '''Check if GATK is available, set GATK_PATH from CLASSPATH if the path
            is specified in CLASSPATH'''
        if env.GATK_PATH:
            if not os.path.isfile(os.path.join(os.path.expanduser(env.GATK_PATH),
                    'GenomeAnalysisTK.jar')):
                env.logger.error('Specified GATK_PATH {} does not contain GATK jar files.'
                    .format(env.GATK_PATH))
                sys.exit(1)
        elif 'CLASSPATH' in os.environ:
            if not any([os.path.isfile(os.path.join(os.path.expanduser(x), 'GenomeAnalysisTK.jar'))
                    for x in os.environ['CLASSPATH'].split(':')]):
                env.logger.error('$CLASSPATH ({}) does not contain a path that contain GATK jar files.'
                    .format(os.environ['CLASSPATH']))
                sys.exit(1)
            else:
                for x in os.environ['CLASSPATH'].split(':'):
                    if os.path.isfile(os.path.join(os.path.expanduser(x), 'GenomeAnalysisTK.jar')):
                        env.logger.info('Using GATK under {}'.format(x))
                        env.GATK_PATH = os.path.expanduser(x)
                        break
        else:
            env.logger.error('Please either specify path to GATK using option '
                '--set GATK_PATH=path, or set it in environment variable $CLASSPATH.')
            sys.exit(1) 
        #
        self.REF_fasta = ref_fasta
        self.knownSites = knownSites


    def realignIndels(self, bam_file):
        '''Create realigner target and realign indels'''
        target = os.path.join(env.WORKING_DIR, os.path.basename(bam_file)[:-4] + '.IndelRealignerTarget.intervals')
        if existAndNewerThan(ofiles=target, ifiles=bam_file):
            env.logger.warning('Using existing realigner target {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                -R {4}/{5}
                -T RealignerTargetCreator
                --mismatchFraction 0.0
                -o {6} {7} '''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_REALIGNERTARGETCREATOR,
                    bam_file, env.RESOURCE_DIR,
                    self.REF_fasta, TEMP(target),
                    ' '.join(['-known {}/{}'.format(env.RESOURCE_DIR, x) for x in self.knownSites])),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        # 
        # realign around known indels
        cleaned_bam_file = os.path.join(env.WORKING_DIR, os.path.basename(bam_file)[:-4] + '.clean.bam')
        if existAndNewerThan(ofiles=cleaned_bam_file, ifiles=target):
            env.logger.warning('Using existing realigner bam file {}'.format(cleaned_bam_file))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                -R {4}/{5}
                -T IndelRealigner 
                --targetIntervals {6}
                --consensusDeterminationModel USE_READS
                -compress 0 -o {7} {8}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_INDELREALIGNER,
                    bam_file, env.RESOURCE_DIR,
                    self.REF_fasta, target, TEMP(cleaned_bam_file),
                    ' '.join(['-known {}/{}'.format(env.RESOURCE_DIR, x) for x in self.knownSites])),
                name=os.path.basename(cleaned_bam_file),
                upon_succ=(os.rename, TEMP(cleaned_bam_file), cleaned_bam_file))
        # 
        return cleaned_bam_file


    def recalibrate(self, bam_file, recal_bam_file):
        '''Create realigner target and realign indels'''
        target = os.path.join(env.WORKING_DIR, os.path.basename(bam_file)[:-4] + '.grp')
        if existAndNewerThan(ofiles=target, ifiles=bam_file):
            env.logger.warning('Using existing base recalibrator target {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} 
                -T BaseRecalibrator
                -I {3} 
                -R {4}/{5}
                -cov ReadGroupCovariate
                -cov QualityScoreCovariate
                -cov CycleCovariate
                -cov ContextCovariate
                -o {6} {7}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_BASERECALIBRATOR, bam_file,
                    env.RESOURCE_DIR, self.REF_fasta, TEMP(target),
                    ' '.join(['-knownSites {}/{}'.format(env.RESOURCE_DIR, x) for x in self.knownSites])),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        #
        # recalibrate
        if existAndNewerThan(ofiles=recal_bam_file, ifiles=target):
            env.logger.warning('Using existing recalibrated bam file {}'.format(recal_bam_file))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2}
                -T PrintReads
                -I {3} 
                -R {4}/{5}
                -BQSR {6}
                -o {7}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_PRINTREADS, bam_file,
                    env.RESOURCE_DIR, self.REF_fasta, target, 
                    TEMP(recal_bam_file)),
                name=os.path.basename(recal_bam_file),
                upon_succ=(os.rename, TEMP(recal_bam_file), recal_bam_file))
        # 
        return recal_bam_file

    def reduceReads(self, input_file):
        target = input_file[:-4] + '_reduced.bam'
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing reduced bam file {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2}
                -T ReduceReads
                -I {3} 
                -R {4}/{5}
                -o {6}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_REDUCEREADS, input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    TEMP(target)),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        # 
        return target

    def unifiedGenotyper(self, input_file):
        target = os.path.join(env.WORKING_DIR,
            os.path.basename(input_file)[:-4] + '_raw.vcf')
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing called variants {}'.format(target))
        else:
            dbSNP_vcf = [x for x in self.knownSites if 'dbsnp' in x][0]
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar
                -T UnifiedGenotyper
                {2}
                -I {3} 
                -R {4}/{5}
                --dbsnp {4}/{7}
                --genotype_likelihoods_model BOTH
                -o {6}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_UNIFIEDGENOTYPER, input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    TEMP(target), dbSNP_vcf),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        # 
        return target

    def haplotypeCall(self, input_file):
        target = os.path.join(env.WORKING_DIR,
            os.path.basename(input_file)[:-4] + '.vcf')
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing called variants {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                -R {4}/{5}
                -T HaplotypeCaller
                -o {6}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_HAPLOTYPECALLER, input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    TEMP(target)),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        # 
        return target


    def variantRecalibration(self, input_file):
        #
        # SNPs error model
        #
        # target is log file
        recal_files = {} 
        target = os.path.join(env.WORKING_DIR,
            os.path.basename(input_file)[:-4] + '.SNPs')
        recal_files['SNP'] = target
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing error model {}'.format(target))
        else:
            dbSNP_vcf = [x for x in self.knownSites if 'dbsnp' in x][0]
            hapmap_vcf = [x for x in self.knownSites if 'hapmap' in x][0]
            kg_vcf = [x for x in self.knownSites if '1000G_omni' in x][0]
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2}
                -input {4}
                -T VariantRecalibrator
                -R {5}/{6}
                {3}
                -resource:hapmap,known=false,training=true,truth=true,prior=15.0 {5}/{7}
                -resource:omni,known=false,training=true,truth=true,prior=12.0	{5}/{8}
                -resource:dbsnp,known=true,training=false,truth=false,prior=2.0	{5}/{9}
                -mode SNP 
                -recalFile {10}.recal
                -tranchesFile {10}.tranches
                -rscriptFile {10}.R
                -log {11}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_VARIANTRECALIBRATOR,
                    env.OPT_GATK_VARIANTRECALIBRATOR_SNV,
                    input_file,
                    env.RESOURCE_DIR,
                    self.REF_fasta,
                    hapmap_vcf,
                    kg_vcf,
                    dbSNP_vcf,
                    target,
                    TEMP(target)),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        #
        # INDELs error model
        #
        # target is log file
        target = os.path.join(env.WORKING_DIR,
            os.path.basename(input_file)[:-4] + '.INDELs')
        recal_files['INDEL'] = target
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing error model {}'.format(target))
        else:
            dbSNP_vcf = [x for x in self.knownSites if 'dbsnp' in x][0]
            indels_vcf = [x for x in self.knownSites if 'gold_standard.indels' in x][0]
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} 
                -input {4} 
                -T VariantRecalibrator
                -R {5}/{6}
                {3}
                -resource:mills,known=false,training=true,truth=true,prior=12.0 {5}/{7}
                -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {5}/{8}
                -mode INDEL 
                -recalFile {9}.recal
                -tranchesFile {9}.tranches
                -rscriptFile {9}.R
                -log {10}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_VARIANTRECALIBRATOR, 
                    env.OPT_GATK_VARIANTRECALIBRATOR_INDEL,
                    input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    indels_vcf,
                    dbSNP_vcf,
                    target, TEMP(target)),
                name=os.path.basename(target),
                upon_succ=(os.rename, TEMP(target), target))
        #
        # SNP recalibration
        #
        target = os.path.join(env.WORKING_DIR,
                              os.path.basename(input_file)[:-4] + '.recal.SNPs.vcf')
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing recalibrated variants {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} {3}
                --input {4} 
                -R {5}/{6}
                -T ApplyRecalibration
                -mode {7} 
                -recalFile {8}.recal
                -tranchesFile {8}.tranches
                -o {9}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_APPLYRECALIBRATION,
                    env.OPT_GATK_APPLYRECALIBRATION_SNV,
                    input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    "SNP", recal_files["SNP"], TEMP(target)),
            name=os.path.basename(target),
            upon_succ=(os.rename, TEMP(target), target))
        #
        # INDEL recalibration
        # input file should be the output from the last step
        #
        input_file = target
        target = os.path.join(env.WORKING_DIR,
                              os.path.basename(input_file)[:-4] + '.INDELs.vcf')
        if existAndNewerThan(ofiles=target, ifiles=input_file):
            env.logger.warning('Using existing recalibrated variants {}'.format(target))
        else:
            run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} {3}
                --input {4} 
                -R {5}/{6}
                -T ApplyRecalibration
                -mode {7} 
                -recalFile {8}.recal
                -tranchesFile {8}.tranches
                -o {9}'''.format(
                    env.OPT_JAVA, env.GATK_PATH,
                    env.OPT_GATK_APPLYRECALIBRATION,
                    env.OPT_GATK_APPLYRECALIBRATION_INDEL,
                    input_file,
                    env.RESOURCE_DIR, self.REF_fasta,
                    "INDEL", recal_files["INDEL"], TEMP(target)),
            name=os.path.basename(target),
            upon_succ=(os.rename, TEMP(target), target))
        return target 


###############################################################################
#
#  Pipelines.
#
#  All pipelines should be subclasses of BasePipeline, and call related functions
#  in BasePipeline before specific actions are taken.
#
#  To add a pipeline, please 
#
#  1. derive a class from BasePipeline or one of its subclasses
#  2. define prepareResourceIfNotExist, align and call if they differ
#     from the function in the parent class.
#
###############################################################################
class BasePipeline:
    '''A vase variant caller that is supposed to provide most of the utility functions
    and common operations that are not specific to any pipeline.
    '''
    def __init__(self):
        pass

    #
    # PREPARE RESOURCE
    #
    # interface
    def prepareResourceIfNotExist(self):
        '''Prepare all resources for the pipeline. '''
        # download test data
        saved_dir = os.getcwd()
        os.chdir(env.RESOURCE_DIR)
        downloadFile('http://vtools.houstonbioinformatics.org/test_data/illumina_test_seq.tar.gz',
            dest='illumina_test_seq.tar.gz')
        if existAndNewerThan(ofiles=['illumina_test_seq_1.fastq', 'illumina_test_seq_2.fastq'],
                ifiles='illumina_test_seq.tar.gz'):
            env.logger.warning('Using existing test sequence files illumina_test_seq_1/2.fastq')
        else:
            decompress('illumina_test_seq.tar.gz', '.')
        os.chdir(saved_dir)

    #
    # align and create bam file
    #
    def align(self, input_files, output):
        '''Align to the reference genome'''
        if not output.endswith('.bam'):
            env.logger.error('Plase specify a .bam file in the --output parameter')
            sys.exit(1)

    def callVariants(self, input_files, ped_file, output):
        '''Call variants from a list of input files'''
        if not output.endswith('.vcf'):
           env.logger.error('Please specify a .vcf file in the --output parameter')
           sys.exit(1)
        for bam_file in input_files:
            if not os.path.isfile(bam_file):
                env.logger.error('Input file {} does not exist'.format(bam_file))
                sys.exit(1)
            if not os.path.isfile(bam_file + '.bai'):
                env.logger.error('Input bam file {} is not indexed.'.format(bam_file))
                sys.exit(1)

class b37_gatk_23(BasePipeline):
    '''A variant caller that uses gatk resource package 2.3 to call variants
    from build b37 of the human genome of the GATK resource bundle'''
    def __init__(self):
        BasePipeline.__init__(self)
        self.pipeline = 'b37_gatk_23'
        self.GATK_resource_url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/b37/*'
        self.REF_fasta = 'human_g1k_v37.fasta'
        self.knownSites = [
            'dbsnp_137.b37.vcf',
            'hapmap_3.3.b37.vcf',
            '1000G_omni2.5.b37.vcf',
            'Mills_and_1000G_gold_standard.indels.b37.vcf',
            '1000G_phase1.indels.b37.vcf',
        ]
        #
        # prepare the commands to be used
        self.bwa = BWA(self.REF_fasta)
        self.samtools = SAMTOOLS(self.REF_fasta)
        self.picard = PICARD()
        self.gatk = GATK(self.REF_fasta, self.knownSites)
        #
        env.RESOURCE_DIR = os.path.join(os.path.expanduser(env.RESOURCE_DIR), self.pipeline)
        if not os.path.isdir(env.RESOURCE_DIR):
            env.logger.info('Creating resource directory {}'.format(env.RESOURCE_DIR))
            os.makedirs(env.RESOURCE_DIR)

    def checkResource(self):
        '''Check if needed resource is available. This pipeline requires
        GATK resource bundle, commands wget, bwa, samtools, picard, and GATK. '''
        saved_dir = os.getcwd()
        os.chdir(env.RESOURCE_DIR)
        files = [self.REF_fasta, self.REF_fasta + '.amb', self.REF_fasta + '.fai']
        if not all([os.path.isfile(x) for x in files]):
            sys.exit('GATK resource bundle does not exist in directory {}. '
                'Please run "call_variants.py prepare_resource" befor you '
                'execute this command.'.format(env.RESOURCE_DIR))
        os.chdir(saved_dir)

    def prepareResourceIfNotExist(self):
        '''This function downloads the UCSC resource boundle for specified build and
        creates bwa and samtools indexed files from the whole genome sequence '''
        # download test data
        BasePipeline.prepareResourceIfNotExist(self)
        # these are pipeline specific
        GATK(self.REF_fasta, self.knownSites).downloadResourceBundle(
            self.GATK_resource_url, files=[self.REF_fasta + '.gz'])
        BWA(self.REF_fasta).buildIndex()
        SAMTOOLS(self.REF_fasta).buildIndex()

    #
    # This is the actual pipeline ....
    #
    def align(self, input_files, output):
        '''Align reads to hg19 reference genome and return a sorted, indexed bam file.'''
        BasePipeline.align(self, input_files, output)
        #
        # the working dir is a directory under output, the internediate files are saved to this
        # directory to avoid name conflict
        if not env._WORKING_DIR:
            env.WORKING_DIR = os.path.join(os.path.split(output)[0], os.path.basename(output) + '_align_cache')
        if not os.path.isdir(env.WORKING_DIR):
            os.makedirs(env.WORKING_DIR)
        #
        # step 1: decompress to get a list of fastq files
        fastq_files = self.getFastqFiles(input_files)
        #
        # step 2: call bwa aln to produce .sai files
        self.bwa.aln(fastq_files)
        #
        # step 3: generate .sam files for each pair of pairend reads, or reach file of unpaired reads
        paired = True
        if len(fastq_files) // 2 * 2 != len(fastq_files):
            env.logger.warning('Odd number of fastq files found, not handled as paired end reads.')
            paired = False
        for idx in range(len(fastq_files)//2):
            f1 = fastq_files[2*idx]
            f2 = fastq_files[2*idx + 1]
            if len(f1) != len(f2):
                env.logger.warning(
                    'Filenames {}, {} are not paired, not handled as paired end reads.'
                    .format(f1, f2))
                paired = False
                break
            diff = [ord(y)-ord(x) for x,y in zip(f1, f2) if x!=y]
            if diff != [1]:
                env.logger.warning(
                    'Filenames {}, {} are not paired, not handled as paired end reads.'
                    .format(f1, f2))
                paired = False
                break
        #
        # sam files?
        if paired:
            sam_files = self.bwa.sampe(fastq_files, output)
        else:
            sam_files = self.bwa.samse(fastq_files, output)
        #
        # check if mapping is successful
        counts = countUnmappedReads(sam_files)
        for f,c in zip(sam_files, counts):
            # more than 40% unmapped
            if c[1] == 0 or c[0]/c[1] > 0.4:
                env.logger.error('{}: {} out of {} reads are unmapped ({:.2f}% mapped)'
                    .format(f, c[0], c[1], 0 if c[1] == 0 else (100*(1 - c[0]/c[1]))))
                sys.exit(1)
            else:
                env.logger.info('{}: {} out of {} reads are unmapped ({:.2f}% mapped)'
                    .format(f, c[0], c[1], 100*(1 - c[0]/c[1])))
        # 
        # step 4: merge per-lane bam files to a single sample bam file
        if len(sam_files) > 1:
            merged_bam_file = self.picard.mergeBAMs(sam_files)
        else:
            merged_bam_file = sam_files[0]

        # step 4: convert sam to sorted bam files
        sorted_bam_file = self.picard.sortSam([merged_bam_file])[0]
        #
        # According to GATK best practice, dedup should be run at the
        # lane level (the documentation is confusing here though)
        #
        # step 5: remove duplicate
        dedup_bam_file = self.picard.markDuplicates([sorted_bam_file])[0]
        #
        #
        # step 7: index the output bam file
        self.samtools.indexBAM(dedup_bam_file)
        #
        # step 8: create indel realignment targets and recall
        cleaned_bam_file = self.gatk.realignIndels(dedup_bam_file)
        self.samtools.indexBAM(cleaned_bam_file)
        #
        # step 9: recalibration
        self.gatk.recalibrate(cleaned_bam_file, output)
        self.samtools.indexBAM(output)
        #
        # step 10: reduce reads
        reduced = self.gatk.reduceReads(output)
        self.samtools.indexBAM(reduced)

    def callVariants(self, input_files, pedfile, output):
        '''Call variants from a list of input files'''
        BasePipeline.callVariants(self, input_files, pedfile, output)
        #
        if not env._WORKING_DIR:
            env.WORKING_DIR = os.path.join(os.path.split(output)[0],
                os.path.basename(output) + '_call_cache')
        if not os.path.isdir(env.WORKING_DIR):
            os.makedirs(env.WORKING_DIR)
        #
        for input_file in input_files:
            #vcf_file = self.haplotypeCall(input_file)
            vcf_file = self.gatk.unifiedGenotyper(input_file)
            #   
            # variant calibrate
            recal_vcf_file = self.gatk.variantRecalibration(vcf_file)
        #
        # copy results to output
        shutil.copy(vcf_file, output)


class hg19_gatk_23(b37_gatk_23):
    '''A variant caller that uses gatk resource package 2.3 to call variants
    from build hg19 of the human genome'''
    def __init__(self):
        # do not call 
        BasePipeline.__init__(self)
        self.pipeline = 'hg19_gatk_23'
        # this piple just uses different resource bundle
        self.GATK_resource_url = 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/hg19/*'
        self.REF_fasta = 'ucsc.hg19.fasta'
        self.knownSites = [
            'dbsnp_137.hg19.vcf',
            'hapmap_3.3.hg19.vcf',
            '1000G_omni2.5.hg19.vcf',
            'Mills_and_1000G_gold_standard.indels.hg19.vcf',
            '1000G_phase1.indels.hg19.vcf',
        ]
        # prepare the commands to be used
        self.bwa = BWA(self.REF_fasta)
        self.picard = PICARD()
        self.samtools = SAMTOOLS(self.REF_fasta)
        self.gatk = GATK(self.REF_fasta, self.knownSites)
        #
        env.RESOURCE_DIR = os.path.join(os.path.expanduser(env.RESOURCE_DIR), self.pipeline)
        if not os.path.isdir(env.RESOURCE_DIR):
            env.logger.info('Creating resource directory {}'.format(env.RESOURCE_DIR))
            os.makedirs(env.RESOURCE_DIR)
        


def alignReadsArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''One or more .txt, .fa, .fastq, .tar, .tar.gz, .tar.bz2, .tbz2, .tgz files
            that contain raw reads of a single sample. Files in sam/bam format are also
            acceptable, in which case raw reads will be extracted and aligned again to 
            generate a new bam file. ''')
    parser.add_argument('-o', '--output',
        help='''Output aligned reads to a sorted, indexed, dedupped, and recalibrated
            BAM file $output.bam.''')
    parser.add_argument('--pipeline', required=True,
        help='Name of the pipeline to be used to call variants.')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs.''')

def alignReads(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            pipeline = Pipeline(args.pipeline)
            # 
            # initialize
            pipeline.downloadResource()
            pipeline.execute('init')
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


def callVariantsArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''One or more BAM files.''')
    parser.add_argument('-o', '--output',
        help='''Output parsered variants to the specified VCF file''')
    parser.add_argument('--pedfile',
        help='''A pedigree file that specifies the relationship between input
            samples, used for multi-sample parsering.''')
    parser.add_argument('--pipeline', required=True,
        help='Name of the pipeline to be used to call variants.')
    parser.add_argument('-j', '--jobs', default=1, type=int,
            help='''Maximum number of concurrent jobs.''')

def callVariants(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            pass
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools
    pass
