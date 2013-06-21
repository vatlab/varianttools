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

from .utils import env, ProgressBar, downloadURL, \
    existAndNewerThan, TEMP, decompressIfNeeded
    
from .project import PipelineDescription, Project

try:
    import pysam
    hasPySam = True
except ImportError as e:
    hasPySam = False


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


class CheckJavaClasses:
    '''Check the existence of specified java class (.jar files) and raise an
    error if one of the commands does not exist.'''
    def __init__(self, classes):
        if type(classes) == type(''):
            self.java_class = [classes]
        else:
            self.java_class = classes

    def __call__(self, ifiles):
        if 'CLASSPATH' not in os.environ:
            raise RuntimeError('CLASSPATH is not defined.')
        for java_class in self.java_class:
            found = False
            for path in os.environ['CLASSPATH'].split(os.pathsep):
                if os.path.isfile(os.path.join(path, java_class)):
                    found = True
                    env.logger.info('Java class {} is located under {}.'
                        .format(java_class, path))
                    break
            #
            if not found:
                raise RuntimeError('Cannot locate {} from environment variable CLASSPATH.'
                    .format(java_class))
        return ifiles


# NOTE:
#   subprocess.PIPE cannot be used for NGS commands because they tend to send
#   a lot of progress output to stderr, which might block PIPE and cause the
#   command itself to fail, or stall (which is even worse).
#
JOB = namedtuple('JOB', 'proc cmd start_time stdout stderr output')
class RunCommand:
    def __init__(self, cmd='', working_dir=None, output=None):
        '''This action execute the specified command under the
        specified working directory, and return specified ofiles.
        '''
        self.cmd = cmd
        if not self.cmd:
            raise RuntimeError('Invalid command to execute: "{}"'.format(cmd))
        self.working_dir = working_dir
        if type(output) == str:
            self.output = [output]
        else:
            self.output = output
        self.running_jobs = []

    def elapsed_time(self, start):
        '''Return the elapsed time in human readable format since start time'''
        second_elapsed = int(time.time() - start)
        days_elapsed = second_elapsed // 86400
        return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
            time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
     
    #
    # this command duplicate with runCommand, will merge them later on.
    def run_command(self, cmd, output=None, wait=True):
        '''Call an external command, raise an error if it fails.
        '''
        # merge mulit-line command into one line and remove extra white spaces
        cmd = ' '.join(cmd.split())
        if not output:
            # subprocess.DEVNULL was introduced in Python 3.3
            proc_out = open(os.devnull, 'w')
            proc_err = open(os.devnull, 'w')
        else:
            proc_out = open(output[0] + '.out', 'w')
            proc_err = open(output[0] + '.err', 'w')
        if wait or env.jobs == 1:
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
                        with open(output[0] + '.err') as err:
                            for line in err.read().split('\n')[-20:]:
                                env.logger.error(line)
                    raise RuntimeError("Command {} returned {} after executing {}"
                        .format(cmd, retcode, elapsed_time(s)))
                env.logger.info('Command {} completed successfully in {}'
                    .format(cmd, elapsed_time(s)))
            except OSError as e:
                env.logger.error("Execution of command {} failed: {}".format(cmd, e))
                sys.exit(1)
        else:
            # wait for empty slot to run the job
            while True:
                if poll_jobs() >= env.jobs:
                    time.sleep(5)
                else:
                    break
            # there is a slot, start running
            proc = subprocess.Popen(cmd, shell=True, stdout=proc_out, stderr=proc_err)
            self.running_jobs.append(JOB(proc=proc, cmd=cmd, 
                start_time=time.time(), stdout=proc_out, stderr=proc_err, output=output))
            env.logger.info('Running {}'.format(cmd))

    def poll_jobs(self):
        '''check the number of running jobs.'''
        count = 0
        for idx, job in enumerate(self.running_jobs):
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
                    with open(job.output[0] + '.err') as err:
                        for line in err.read().split('\n')[-50:]:
                            env.logger.error(line)
                raise RuntimeError('Execution of command {} failed after {} (return code {}).'
                    .format(job.cmd, elapsed_time(job.start_time), ret))
            else:
                if job.output:
                    with open(job.output[0] + '.err') as err:
                        for line in err.read().split('\n')[-10:]:
                            env.logger.info(line)
                env.logger.info('Command {} completed successfully in {}'
                    .format(job.cmd, elapsed_time(job.start_time)))
                #
                self.running_jobs[idx] = None
        return count

    def wait_all(self):
        '''Wait for all pending jobs to complete'''
        while poll_jobs() > 0:
            # sleep ten seconds before checking job status again.
            time.sleep(10)
        self.running_jobs = []

    def __call__(self, ifiles):
        # substitute cmd by input_files and output_files
        if existAndNewerThan(ifiles=ifiles, ofiles=self.output):
            env.logger.info('Reuse existing files {}'
                .format(', '.join(self.output)))
            return self.output
        if self.working_dir:
            os.chdir(self.working_dir)
        self.run_command(self.cmd, output=self.output)
        return self.output


class GetFastqFiles:
    '''This action gets a list of fastq files from input file, decompressing
    input files (.tar.gz, .zip, etc) if necessary. Non-fastq files are ignored
    with a warning message. '''
    def __init__(self):
        pass

    def _decompress(self, filename, dest_dir):
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
        
    def __call__(self, ifiles):
        # decompress input files and return a list of output files
        filenames = []
        for filename in ifiles:
            for fastq_file in self._decompress(filename, env.cache_dir):
                try:
                    with open(fastq_file) as fastq:
                        line = fastq.readline()
                        if not line.startswith('@'):
                            env.logger.warning('Wrong FASTA file {}'.format(fastq_file))
                            continue
                    filenames.append(fastq_file)
                except Exception as e:
                    env.logger.warning('Ignoring non-fastq file {}: {}'
                        .format(fastq_file, e))
        filenames.sort()
        return filenames        


class CountUnmappedReads:
    '''This action reads the input files in sam format and count the total
    number of reads and number of unmapped reads. It raises a RuntimeError
    if the proportion of unmapped reads exceeds the specified cutoff value
    (default to 0.2). This action writes a count file to specified output
    files and read from this file directly if the file already exists.'''
    def __init__(self, cutoff=0.2):
        self.cutoff = cutoff

    def __call__(self, ifiles):
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


class Pipeline:
    def __init__(self, name):
        self.pipeline = PipelineDescription(name)
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
 
    def substitute(self, action, VARS):
        # FIXME: advanced cases not handled
        for key, value in VARS.iteritems():
            action = action.replace('${' + key + '}', value)
        return action

    def execute(self, steps, input_files=[], output_files=[]):
        VARS = {
            'CMD_INPUT': ', '.join(input_files),
            'CMD_OUTPUT': ', '.join(output_files),
            'RESOURCE_DIR': self.pipeline_resource,
            'TEMP_DIR': env.temp_dir,
            'CACHE_DIR': env.cache_dir,
            'PID': str(os.getpid())
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
            # 
            # now, execute it
            try:
                saved_dir = os.getcwd()
                VARS['INPUT'] = ', '.join(step_input)
                VARS['PID'] = str(os.getpid())
                action = eval(self.substitute(command.action, VARS))
                if type(action) == tuple:
                    action = SequentialActions(action)
                ifiles = action(step_input)
                os.chdir(saved_dir)
            except Exception as e:
                raise RuntimeError('Failed to execute step {} of {}: {}'
                    .format(command.index, steps, e))
                

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
    try:
        with Project(verbosity=args.verbosity) as proj:
            pipeline = Pipeline(args.pipeline)
            # 
            # initialize
            pipeline.downloadResource()
            pipeline.execute('init', args.input_files, args.output)
            pipeline.execute('align', args.input_files, args.output)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


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
