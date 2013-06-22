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

from .utils import env, ProgressBar, downloadURL, calculateMD5, \
    existAndNewerThan, TEMP, decompressIfNeeded
    
from .project import PipelineDescription, Project

try:
    import pysam
    hasPySam = True
except ImportError as e:
    hasPySam = False

###################################

class GroupInput:
    '''Select input files of certain types, group them, and send input files
    to action. filetype can be all (all input file types), fastq (check 
    content of files), or one or more file extensions (e.g. ['sam', 'bam']).
    Eligible files are by default sent individually (group_by='single') to 
    action (${INPUT} is a list of a single file), but can also be sent all
    together (group_by='all', ${INPUT} equals to ${INPUT#} where # is the
    index of step) or in pairs (group_by='paired', e.g. filename_1.txt and
    filename_2.txt). Unselected files are by default passed directly as 
    output of a step.'''
    def __init__(self, group_by='single', filetypes='all', pass_unselected=True):
        self.group_by = group_by
        if type(filetypes) == str:
            self.filetypes = [filetypes]
        else:
            self.filetypes = filetypes
        self.pass_unselected = pass_unselected

    def _isFastq(self, filename):
        try:
            with open(filename) as fastq:
                line = fastq.readline()
                if not line.startswith('@'):
                    return False
        except Exception as e:
            return False

    def __call__(self, ifiles):
        selected = []
        unselected = []
        for filename in ifiles:
            for t in self.filetypes:
                if t == 'all':
                    selected.append(filename)
                    break
                if t == 'fastq':
                    if self._isFastq(filename):
                        selected.append(filename)
                        break
                if filename.lower().endswith('.' + t.lstrip('.').lower()):
                    selected.append(filename)
            if self.pass_unselected and (not selected or selected[-1] != filename):
                unselected.append(filename)
        # 
        if self.group_by == 'single':
            return [[x] for x in selected], unselected
        elif self.group_by == 'all':
            return [selected], unselected
        elif self.group_by == 'paired':
            if len(selected) // 2 * 2 != len(selected):
                raise RuntimeError('Cannot emit input files by pairs: odd number of input files {}'
                    .format(', '.join(selected)))
            selected.sort()
            pairs = []
            for idx in range(len(selected)//2):
                f1 = selected[2*idx]
                f2 = selected[2*idx + 1]
                if len(f1) != len(f2):
                    raise RuntimeError('Cannot emit input files by pairs: '
                        'Filenames {}, {} are not paired'
                        .format(f1, f2))
                diff = [ord(y)-ord(x) for x,y in zip(f1, f2) if x!=y]
                if diff != [1]:
                    raise RuntimeError('Cannot emit input files by pairs: '
                        'Filenames {}, {} are not paired'
                        .format(f1, f2))
                pairs.append([f1, f2])
            return pairs, unselected
        

class SkipIfSingle:
    '''Skip action if there is only one input file'''
    def __init__(self, group_by='single'):
        self.group_by = group_by

    def __call__(self, ifiles):
        if len(ifiles) <= 1:
            return [], ifiles
        elif group_by == 'single':
            return [[x] for x in ifiles], []
        elif group_by == 'all':
            return [ifiles], []


###################################

class SequentialActions:
    '''Define an action that calls a list of actions, specified by Action1,
    Action2 etc. The input of the first Action is ${INPUT} (for the first
    action), or the output of the previous action. The output of the last
    action becomes the output of the action sequence, or $OUTPUT if the last
    action does not return any output.'''
    def __init__(self, actions):
        self.actions = []
        for a in actions:
            if hasattr(a, '__call__'):
                self.actions.append(a.__call__)
            else:
                self.actions.append(a)

    def __call__(self, ifiles):
        for a in self.actions:
            # the input of the next action is the output of the
            # previous action.
            ifiles = a(ifiles)
        # return the output of the last action
        return ifiles


class CheckCommands:
    '''Check the existence of specified commands and raise an error if one of
    the commands does not exist.'''
    def __init__(self, cmds):
        if type(cmds) == type(''):
            self.cmd = [cmds]
        else:
            self.cmd = cmds

    def __call__(self, ifiles):
        if not hasattr(shutil, 'which'):
            # if shutil.which does not exist, use subprocess...
            for cmd in self.cmd:
                try:
                    subprocess.call(cmd, stdout=open(os.devnull, 'w'),
                        stderr=open(os.devnull, 'w'))
                    env.logger.info('Command {} is located.'.format(cmd))
                except:
                    raise RuntimeError('Command {} does not exist. Please '
                        'install it and try again.'.format(self.cmd))
            return ifiles
        for cmd in self.cmd:
            if shutil.which(cmd) is None:
                raise RuntimeError('Command {} does not exist. Please install it and try again.'
                    .format(self.cmd))
            else:
                env.logger.info('Command {} is located.'.format(self.cmd))
        return ifiles


class CheckFiles:
    '''Check the existence of specified files and raise an
    error if one of the files does not exist.'''
    def __init__(self, files):
        if type(files) == type(''):
            self.files = [files]
        else:
            self.files = files

    def __call__(self, ifiles):
        for f in self.files:
            if os.path.isfile(f):
                env.logger.info('{} is located.'.format(f))
            else:
                raise RuntimeError('Cannot locate {}.'.format(f))
        return ifiles


class CheckFastqVersion:
    def __init__(self, output):
        self.output = output

    def __call__(self, fastq_file):
        '''Detect the version of input fastq file. This can be very inaccurate'''
        with open(self.output, 'w') as aln_param:
            #
            # This function assumes each read take 4 lines, and the last line contains
            # quality code. It collects about 1000 quality code and check their range,
            # and use it to determine if it is Illumina 1.3+
            #
            qual_scores = ''
            with open(fastq_file[0]) as fastq:
                while len(qual_scores) < 1000:
                    try:
                        line = fastq.readline()
                    except Exception as e:
                        raise RuntimeError('Failed to read fastq file {}: {}'
                            .format(fastq_file, e))
                    if not line.startswith('@'):
                        raise ValueError('Wrong FASTA file {}'.format(fastq_file))
                    line = fastq.readline()
                    line = fastq.readline()
                    if not line.startswith('+'):
                        env.logger.warning(
                            'Suspiciout FASTA file {}: third line does not start with "+".'
                            .foramt(fastq_file))
                        return 
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
                aln_param.write('-I')
        return self.output


class GuessReadGroup:
    def __init__(self, bamfile, rgfile):
        self.output_bam = bamfile
        self.rg_output = rgfile

    def __call__(self, fastq_filename):
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
        with open(self.rg_output, 'w') as rg_output:
            filename = os.path.basename(fastq_filename[0])
            output = os.path.basename(self.output_bam)
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
            rg_output.write(rg)
        return self.rg_output


# NOTE:
#   subprocess.PIPE cannot be used for NGS commands because they tend to send
#   a lot of progress output to stderr, which might block PIPE and cause the
#   command itself to fail, or stall (which is even worse).
#
JOB = namedtuple('JOB', 'proc cmd start_time stdout stderr output')
running_jobs = []
max_running_jobs = 1

def elapsed_time(start):
    '''Return the elapsed time in human readable format since start time'''
    second_elapsed = int(time.time() - start)
    days_elapsed = second_elapsed // 86400
    return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
        time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
 
def run_command(cmd, output=None, wait=True):
    '''Call an external command, raise an error if it fails.
    '''
    global running_jobs
    if not output:
        # subprocess.DEVNULL was introduced in Python 3.3
        proc_out = open(os.devnull, 'w')
        proc_err = open(os.devnull, 'w')
    else:
        proc_out = open(output[0] + '.out_{}'.format(os.getpid()), 'w')
        proc_err = open(output[0] + '.err_{}'.format(os.getpid()), 'w')
    if wait or max_running_jobs == 1:
        try:
            s = time.time()
            env.logger.info('Running {}'.format(cmd))
            proc = subprocess.Popen(cmd, shell=True, stdout=proc_out, stderr=proc_err)
            retcode = proc.wait()
            proc_out.close()
            proc_err.close()
            if retcode < 0:
                raise RuntimeError('Command {} was terminated by signal {} after executing {}'
                    .format(cmd, -retcode, elapsed_time(s)))
            elif retcode > 0:
                if output:
                    with open(output[0] + '.err_{}'.format(os.getpid())) as err:
                        for line in err.read().split('\n')[-20:]:
                            env.logger.error(line)
                raise RuntimeError("Command {} returned {} after executing {}"
                    .format(cmd, retcode, elapsed_time(s)))
            if output:
                with open(output[0] + '.exe_info', 'w') as exe_info:
                    exe_info.write('{}\nStart: {}\nEnd: {}\n'
                        .format(cmd, time.asctime(time.localtime(s)), time.asctime(time.localtime())))
            env.logger.info('Command {} completed successfully in {}'
                .format(cmd, elapsed_time(s)))
        except OSError as e:
            env.logger.error("Execution of command {} failed: {}".format(cmd, e))
            sys.exit(1)
    else:
        # wait for empty slot to run the job
        while True:
            if poll_jobs() >= max_running_jobs:
                time.sleep(5)
            else:
                break
        # there is a slot, start running
        proc = subprocess.Popen(cmd, shell=True, stdout=proc_out, stderr=proc_err)
        running_jobs.append(JOB(proc=proc, cmd=cmd, 
            start_time=time.time(), stdout=proc_out, stderr=proc_err, output=output))
        env.logger.info('Running {}'.format(cmd))

def poll_jobs():
    '''check the number of running jobs.'''
    global running_jobs
    count = 0
    for idx, job in enumerate(running_jobs):
        if job is None:
            continue
        ret = job.proc.poll()
        if ret is None:  # still running
            count += 1
            continue
        #
        # job completed, close redirected stdout and stderr
        job.stdout.close()
        job.stderr.close()
        #
        if ret < 0:
            raise RuntimeError("Command {} was terminated by signal {} after executing {}"
                .format(job.cmd, -ret, elapsed_time(job.start_time)))
        elif ret > 0:
            if job.output:
                with open(job.output[0] + '.err_{}'.format(os.getpid())) as err:
                    for line in err.read().split('\n')[-50:]:
                        env.logger.error(line)
            raise RuntimeError('Execution of command {} failed after {} (return code {}).'
                .format(job.cmd, elapsed_time(job.start_time), ret))
        else:
            if job.output:
                with open(job.output[0] + '.err_{}'.format(os.getpid())) as err:
                    for line in err.read().split('\n')[-10:]:
                        env.logger.info(line)
                with open(job.output[0] + '.exe_info', 'w') as exe_info:
                    exe_info.write('{}\nStart: {}\nEnd:{}\n'
                        .format(job.cmd, job.start_time, time.time()))
            env.logger.info('Command {} completed successfully in {}'
                .format(job.cmd, elapsed_time(job.start_time)))
            #
            running_jobs[idx] = None
    return count

def wait_all():
    '''Wait for all pending jobs to complete'''
    global running_jobs
    while poll_jobs() > 0:
        # sleep ten seconds before checking job status again.
        time.sleep(10)
    running_jobs = []


class RunCommand:
    def __init__(self, cmd='', working_dir=None, output=[]):
        '''This action execute the specified command under the
        specified working directory, and return specified ofiles.
        '''
        # merge mulit-line command into one line and remove extra white spaces
        self.cmd = ' '.join(cmd.split())
        if not self.cmd:
            raise RuntimeError('Invalid command to execute: "{}"'.format(cmd))
        self.working_dir = working_dir
        if type(output) == str:
            self.output = [output]
        else:
            self.output = output

    def __call__(self, ifiles):
        # substitute cmd by input_files and output_files
        if self.output:
            if os.path.isfile(self.output[0] + '.exe_info'):
                with open(self.output[0] + '.exe_info') as exe_info:
                    cmd = exe_info.readline().strip()
                # if the exact command has been used to produce output, and the
                # output files are newer than input file, ignore the step
                if cmd == self.cmd.strip() and existAndNewerThan(ifiles=ifiles,
                        ofiles=self.output, md5=self.output[0] + '.exe_info'):
                    env.logger.info('Reuse existing files {}'.format(', '.join(self.output)))
                    return self.output
        if self.working_dir:
            os.chdir(self.working_dir)
        run_command(self.cmd, output=self.output, wait=False)
        # add md5 signature of input and output files
        if self.output:
            with open(self.output[0] + '.exe_info', 'a') as exe_info:
                for f in ifiles + self.output:
                    # for performance considerations, use partial MD5
                    exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                        calculateMD5(f, partial=True)))
        return self.output


class DecompressFiles:
    '''This action gets a list of fastq files from input file, decompressing
    input files (.tar.gz, .zip, etc) if necessary. Non-fastq files are ignored
    with a warning message. '''
    def __init__(self, dest_dir=None):
        self.dest_dir = dest_dir if dest_dir else '.'

    def _decompress(self, filename):
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
            dest_file = os.path.join(self.dest_dir, os.path.basename(filename)[:-3])
            if existAndNewerThan(ofiles=dest_file, ifiles=filename):
                env.logger.info('Using existing decompressed file {}'.format(dest_file))
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
            dest_file = os.path.join(self.dest_dir, os.path.basename(filename)[:-4])
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
            bundle.extractall(self.dest_dir)
            env.logger.info('Decompressing {} to {}'.format(filename, self.dest_dir))
            return [os.path.join(self.dest_dir, name) for name in bundle.namelist()]
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
            manifest = os.path.join(self.dest_dir, os.path.basename(filename) + '.manifest')
            all_extracted = False
            dest_files = []
            if existAndNewerThan(ofiles=manifest, ifiles=filename):
                all_extracted = True
                for f in [x.strip() for x in open(manifest).readlines()]:
                    dest_file = os.path.join(self.dest_dir, os.path.basename(f))
                    if existAndNewerThan(ofiles=dest_file, ifiles=filename):
                        dest_files.append(dest_file)
                        env.logger.info('Using existing extracted file {}'.format(dest_file))
                    else:
                        all_extracted = False
            #
            if all_extracted:
                return dest_files
            #
            # create a temporary directory to avoid corrupted file due to interrupted decompress
            try:
                os.mkdir(os.path.join(self.dest_dir, 'tmp'))
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
                    dest_file = os.path.join(self.dest_dir, os.path.basename(f))
                    dest_files.append(dest_file)
                    if existAndNewerThan(ofiles=dest_file, ifiles=filename):
                        env.logger.info('Using existing extracted file {}'.format(dest_file))
                    else:
                        env.logger.info('Extracting {} to {}'.format(f, dest_file))
                        tar.extract(f, os.path.join(self.dest_dir, 'tmp'))
                        # move to the top directory with the right name only after the file has been properly extracted
                        shutil.move(os.path.join(self.dest_dir, 'tmp', f), dest_file)
                # set dest_files to the same modification time. This is used to
                # mark the right time when the files are created and avoid the use
                # of archieved but should-not-be-used files that might be generated later
                [os.utime(x, None) for x in dest_files]
            return dest_files
        # return source file if nothing needs to be decompressed
        return [filename]
        
    def __call__(self, ifiles):
        # decompress input files and return a list of output files
        filenames = []
        for filename in ifiles:
            filenames.extend(self._decompress(filename))
        filenames.sort()
        return filenames        

class CountMappedReads:
    '''This action reads the input files in sam format and count the total
    number of reads and number of mapped reads. It raises a RuntimeError
    if the proportion of unmapped reads lower than the specified cutoff value
    (default to 0.8). This action writes a count file to specified output
    files and read from this file directly if the file already exists.'''
    def __init__(self, output, cutoff=0.2):
        self.cutoff = cutoff
        self.output = output

    def __call__(self, sam_file):
        #
        # count total reads and unmapped reads
        #
        # The previous implementation uses grep and wc -l, but
        # I cannot understand why these commands are so slow...
        #
        if not existAndNewerThan(ofiles=self.output, ifiles=sam_file):
            env.logger.info('Counting unmapped reads in {}'.format(sam_file[0]))
            unmapped_count = 0
            total_count = 0
            with open(sam_file[0]) as sam:
               for line in sam:
                   total_count += 1
                   if 'XT:A:N' in line:
                       unmapped_count += 1
            mapped_count = total_count - unmapped_count
            with open(self.output, 'w') as target:
                target.write('{}\n{}\n'.format(mapped_count, total_count))
        else:
            env.logger.info('Using existing counts in {}'.format(self.output))
            #
            with open(self.output) as cnt:
                mapped_count = int(cnt.readline())
                total_count = int(cnt.readline())
        # more than 40% unmapped
        if total_count == 0 or mapped_count * 1.0 / total_count < self.cutoff:
            raise RuntimeError('{}: {} out of {} reads ({:.2f}%) are mapped.'
                .format(sam_file[0], mapped_count, total_count,
                    0 if total_count == 0 else (100.*mapped_count/total_count)))
        else:
            env.logger.info('{}: {} out of {} reads ({:.2f}%) are mapped.'
                .format(sam_file[0], mapped_count, total_count,
                    0 if total_count == 0 else (100.*mapped_count/total_count)))
        return self.output


class Pipeline:
    def __init__(self, name, extra_args=[]):
        self.pipeline = PipelineDescription(name, extra_args)
        #
        # resource directory
        #
        if self.pipeline.resource_dir:
            if os.path.isabs(os.path.expanduser(self.pipeline.resource_dir)):
                self.pipeline_resource = os.path.expanduser(self.pipeline.resource_dir)
            else:
                self.pipeline_resource = os.path.join(os.path.expanduser(
                    env.local_resource), 'pipeline_resource', self.pipeline.resource_dir)
        else:
            self.pipeline_resource = os.path.join(os.path.expanduser(
                    env.local_resource), 'pipeline_resource', self.pipeline.name)
        try:
            if not os.path.isdir(self.pipeline_resource):
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
                    decompressIfNeeded(filename, inplace=False)
        os.chdir(saved_dir)
        if skipped:
            env.logger.info('Using {} existing resource files under {}.'
                .format(len(skipped), self.pipeline_resource))
 
    def substitute(self, text, VARS):
        # now, find ${}
        pieces = re.split('(\${[^}]*})', text)
        for idx, piece in enumerate(pieces):
            if piece.startswith('${') and piece.endswith('}'):
                KEY = piece[2:-1]
                if ':' in KEY:
                    # a lambda function?
                    try:
                        FUNC = eval('lambda {}'.format(KEY))
                    except Exception as e:
                        raise RuntimeError('Failed to substitute variable {}: {}'
                            .format(piece, e))
                    KEY = KEY.split(':', 1)[0].strip()
                    if KEY in VARS:
                        VAL = VARS[KEY]
                    else:
                        raise RuntimeError('Failed to substitute variable {}: key {} not found'
                            .format(piece, KEY))
                    try:
                        pieces[idx] = str(FUNC(VAL))
                    except Exception as e:
                        raise RuntimeError('Failed to substitute variable {}: {}'
                            .format(piece, e))
                else:
                    # if KEY in VARS, replace it
                    if KEY in VARS:
                        if type(VARS[KEY]) == str:
                            pieces[idx] = VARS[KEY]
                        else:
                            pieces[idx] = ', '.join(VARS[KEY])
                    else:
                        raise RuntimeError('Failed to substitute variable {}'
                            .format(piece))
        # now, join the pieces together, but remove all newlines
        return ' '.join(''.join(pieces).split())

    def execute(self, steps, input_files=[], output_files=[], jobs=1):
        global max_running_jobs 
        max_running_jobs = jobs
        VARS = {
            'CMD_INPUT': input_files,
            'CMD_OUTPUT': output_files,
            'RESOURCE_DIR': self.pipeline_resource,
            'TEMP_DIR': env.temp_dir,
            'CACHE_DIR': env.cache_dir,
        }
        ifiles = input_files
        for command in {'init': self.pipeline.init_steps,
                'align': self.pipeline.align_steps,
                'call': self.pipeline.call_steps}[steps]:
            # substitute ${} variables
            if command.input:
                step_input = [x.strip() for x in self.substitute(command.input, VARS).split(',')]
            else:
                step_input = ifiles
            # if there is no input file?
            if not step_input:
                raise RuntimeError('Pipeline stops at step {} of {}: No input file is available.'
                    .format(command.index, steps))
            #
            VARS['INPUT{}'.format(command.index)] = step_input
            env.logger.debug('INPUT of step {} of {}: {}'
                    .format(command.index, steps, step_input))
            # 
            # now, group input files
            if not command.input_emitter:
                emitter = GroupInput()
            else:
                try:
                    emitter = eval(command.input_emitter)
                except Exception as e:
                    raise RuntimeError('Failed to group input files: {}'
                        .format(e))
            #
            saved_dir = os.getcwd()
            igroups, step_output = emitter(step_input)
            try:
                for ig in igroups:
                    if not ig:
                        continue
                    VARS['INPUT'] = ig
                    action = self.substitute(command.action, VARS)
                    env.logger.debug('Input: {} Action: {}' .format(ig, action))
                    action = eval(action)
                    if type(action) == tuple:
                        action = SequentialActions(action)
                    ofiles = action(step_input)
                    if type(ofiles) == str:
                        step_output.append(ofiles)
                    else:
                        step_output.extend(ofiles)
                # wait for all pending jobs to finish
                wait_all()
                VARS['OUTPUT{}'.format(command.index)] = step_output
                env.logger.debug('OUTPUT of step {} of {}: {}'
                    .format(command.index, steps, step_output))
            except Exception as e:
                raise RuntimeError('Failed to execute step {} of {}: {}'
                    .format(command.index, steps, e))
            os.chdir(saved_dir)
            ifiles = step_output
                



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


def alignReadsArguments(parser):
    parser.add_argument('input_files', nargs='+',
        help='''One or more .txt, .fa, .fastq, .tar, .tar.gz, .tar.bz2, .tbz2, .tgz files
            that contain raw reads of a single sample. Files in sam/bam format are also
            acceptable, in which case raw reads will be extracted and aligned again to 
            generate a new bam file. ''')
    parser.add_argument('-o', '--output', nargs='+',
        help='''Output aligned reads to a sorted, indexed, dedupped, and recalibrated
            BAM file $output.bam.''')
    parser.add_argument('--pipeline', required=True,
        help='Name of the pipeline to be used to call variants.')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs.''')

def alignReads(args):
    #try:
        with Project(verbosity=args.verbosity) as proj:
            pipeline = Pipeline(args.pipeline, extra_args=args.unknown_args)
            # 
            # initialize
            pipeline.downloadResource()
            pipeline.execute('init', args.input_files, args.output, args.jobs)
            pipeline.execute('align', args.input_files, args.output, args.jobs)
    #except Exception as e:
    #    env.logger.error(e)
    #    sys.exit(1)


def callVariantsArguments(parser):
    parser.add_argument('input_files', nargs='+',
        help='''One or more BAM files.''')
    parser.add_argument('-o', '--output', nargs='+',
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
