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

from .utils import env, ProgressBar, downloadURL, calculateMD5, delayedAction, \
    existAndNewerThan, TEMP, decompressGzFile
    
from .project import PipelineDescription, Project

try:
    import pysam
    hasPySam = True
except (ImportError, ValueError) as e:
    hasPySam = False

class EmitInput:
    '''Select input files of certain types, group them, and send input files
    to action. Selection criteria can be True (all input file types, default),
    False (select no input file), 'fastq' (check content of files), or one or
    more file extensions (e.g. ['sam', 'bam']).  Eligible files are by default
    sent individually (group_by='single') to action (${INPUT} is a list of a
    single file), but can also be sent altogether (group_by='all', ${INPUT}
    equals to ${INPUT#} where # is the index of step) or in pairs 
    (group_by='paired', e.g. filename_1.txt and filename_2.txt). Unselected
    files are by default passed directly as output of a step.'''
    def __init__(self, group_by='single', select=True, pass_unselected=True):
        self.group_by = group_by
        if type(select) == str:
            self.select = [select]
        else:
            self.select = select
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
            match = False
            if self.select == True:
                match = True
            elif self.select == False:
                pass
            else:   # list of types
                for t in self.select:
                    if t == 'fastq':
                        if self._isFastq(filename):
                            match = True
                            break
                    if filename.lower().endswith('.' + t.lstrip('.').lower()):
                        match = True
                        break
            #
            if match:
                selected.append(filename)
            elif self.pass_unselected:
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
                        'install it and try again.'.format(cmd))
            return ifiles
        for cmd in self.cmd:
            if shutil.which(cmd) is None:
                raise RuntimeError('Command {} does not exist. Please install it and try again.'
                    .format(self.cmd))
            else:
                env.logger.info('Command {} is located.'.format(cmd))
        return ifiles


class CheckOutput:
    '''Check out of of an command, and check if it matches a particular
    pattern. The pipeline will exit if fail is set to True (default).'''
    def __init__(self, cmd, pattern, failIfMismatch=True):
        self.cmd = cmd
        self.pattern = pattern
        self.fail = failIfMismatch

    def __call__(self, ifiles):
        try:
            # do not use subprocess.check_output because I need to get
            # output even when the command returns non-zero return code
            p = subprocess.Popen(self.cmd, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=True)
            odata, edata = p.communicate()
            output = odata.decode() + edata.decode()
            env.logger.debug('Output of command "{}" is "{}"'
                .format(self.cmd, output))
        except Exception as e:
            raise RuntimeError('Failed to execute command "{}": {}'
                .format(self.cmd, e))
        #
        if re.search(self.pattern, output, re.MULTILINE) is None:
            msg = ('Output of command "{}" ("{}") does not ' + 
                    'match specified regular expression "{}".').format(self.cmd,
                        ' '.join(output[:40].split()), self.pattern)
            if self.fail:
                raise RuntimeError(msg)
            else:
                env.logger.warning(msg)
        return ifiles

class CheckFiles:
    '''Check the existence of specified files and raise an
    error if one of the files does not exist.'''
    def __init__(self, files):
        if type(files) == str:
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
 
def run_command(cmd, output=None):
    '''Call an external command, raise an error if it fails.
    '''
    global running_jobs
    if not output:
        # subprocess.DEVNULL was introduced in Python 3.3
        proc_out = None
        proc_err = None
        proc_lck = None
    else:
        proc_out = open(output[0] + '.out', 'w')
        proc_err = open(output[0] + '.err', 'w')
        proc_lck = output[0] + '.lck'
    #
    # wait for empty slot to run the job
    while True:
        if poll_jobs() >= max_running_jobs:
            time.sleep(5)
        else:
            break
    # there is a slot, start running
    if proc_lck:
        if os.path.isfile(proc_lck):
            proc_lck = None   # do not remove lock
            raise RuntimeError('Output of pipeline locked by {}.lck. '
                'Please remove this file and try again if no other '
                'process is writing to this file.'
                .format(output[0])) 
        else:
            open(proc_lck, 'a').close()
    proc = subprocess.Popen(cmd, shell=True, stdout=proc_out, stderr=proc_err)
    running_jobs.append(JOB(proc=proc, cmd=cmd, 
        start_time=time.time(), stdout=proc_out, stderr=proc_err, output=output))
    env.logger.info('Running "{}"'.format(cmd))
    if proc_out is not None:
        env.logger.info('Output redirected to {}.out (and .err)'.format(output[0]))

def poll_jobs():
    '''check the number of running jobs.'''
    global running_jobs
    count = 0
    for idx, job in enumerate(running_jobs):
        if job is None:
            continue
        ret = job.proc.poll()
        if ret is None:  # still running
            # flush so that we can check output easily
            if job.stdout is not None:
                job.stdout.flush()
                job.stderr.flush()
            count += 1
            continue
        #
        # job completed, close redirected stdout and stderr
        if job.stdout is not None:
            job.stdout.close()
            job.stderr.close()
        #
        if ret < 0:
            if job.output:
                try:
                    os.remove(job.output[0] + '.lck')
                except:
                    pass
            raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
                .format(job.cmd, -ret, elapsed_time(job.start_time)))
        elif ret > 0:
            if job.output:
                with open(job.output[0] + '.err') as err:
                    for line in err.read().split('\n')[-50:]:
                        env.logger.error(line)
                try:
                    os.remove(job.output[0] + '.lck')
                except:
                    pass
            raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
                .format(job.cmd, elapsed_time(job.start_time), ret))
        else:
            if job.output:
                with open(job.output[0] + '.err') as err:
                    for line in err.read().split('\n')[-10:]:
                        env.logger.info(line)
                with open(job.output[0] + '.exe_info', 'a') as exe_info:
                    exe_info.write('#End: {}\n'.format(time.asctime(time.localtime())))
                    for f in job.output:
                        if not os.path.isfile(f):
                            raise RuntimeError('Output file {} does not exist after completion of the job.'.format(f))
                        # for performance considerations, use partial MD5
                        exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                            calculateMD5(f, partial=True)))
                try:
                    os.remove(job.output[0] + '.lck')
                except Exception as e:
                    self.logger.warning('Failed to remove lock file {}'
                        .format(job.output[0] + '.lck'))
            env.logger.info('Command {} completed successfully in {}'
                .format(job.cmd, elapsed_time(job.start_time)))
            #
            running_jobs[idx] = None
    return count

def wait_all():
    '''Wait for all pending jobs to complete'''
    global running_jobs
    try:
        while poll_jobs() > 0:
            # sleep ten seconds before checking job status again.
            time.sleep(10)
    except KeyboardInterrupt:
        # clean up lock files
        for job in running_jobs:
            if job is not None and job.stdout:
                try:
                    os.remove(job.output[0] + '.lck')
                except:
                    pass
        env.logger.error('Keyboard interrupted')
        sys.exit(1)
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
        #
        for filename in self.output:
            if filename.lower().rsplit('.', 1)[-1] in ['out', 'err', 'lck']:
                raise RuntimeError('Output file with extension .out, .err and .lck are reserved.')

    def __call__(self, ifiles):
        # substitute cmd by input_files and output_files
        if self.output:
            if os.path.isfile(self.output[0] + '.exe_info'):
                with open(self.output[0] + '.exe_info') as exe_info:
                    cmd = exe_info.readline().strip()
                # if the exact command has been used to produce output, and the
                # output files are newer than input file, ignore the step
                if cmd == self.cmd.strip() and existAndNewerThan(ifiles=ifiles,
                        ofiles=self.output, md5file=self.output[0] + '.exe_info'):
                    env.logger.info('Reuse existing files {}'.format(', '.join(self.output)))
                    return self.output
        if self.working_dir:
            os.chdir(self.working_dir)
        run_command(self.cmd, output=self.output)
        # add md5 signature of input and output files
        if self.output:
            with open(self.output[0] + '.exe_info', 'w') as exe_info:
                exe_info.write('{}\n'.format(self.cmd))
                exe_info.write('#Start: {}\n'.format(time.asctime(time.localtime())))
                for f in ifiles:
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

class LinkToDir:
    '''Create hard links of input files to a specified directory. This is 
    usually used to link input files to a common cache directory so that 
    all operations can be performed on that directory.'''
    def __init__(self, dest):
        self.dest = dest
        if not os.path.isdir(self.dest):
            env.logger.info('Creating directory {}'.format(self.dest))
            try:
                os.mkdirs(self.dest)
            except Exception as e:
                raise RuntimeError('Failed to create directory {}'.format(self.dest))
            if not os.path.isdir(self.dest):
                raise RuntimeError('Failed to create directory {}'.format(self.dest))

    def __call__(self, ifiles):
        ofiles = []
        for filename in ifiles:
            path, basename = os.path.split(filename)
            if not os.path.samefile(filename,  os.path.join(self.dest, basename)):
                env.logger.info('Linking {} to {}'.format(filename, self.dest))
                os.link(filename, os.path.join(self.dest, basename))
            ofiles.append(os.path.join(self.dest, basename))
        return ofiles

        
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


class DownloadResource:
    '''Download resources to specified destination directory. resource_dir can
    be a full path name or a directory relative to 
    $local_resource/pipeline_resource where $local_resource is the local
    resource directory of the project (default to ~/.variant_tools,
    see runtime option local_resource for details). The default pipeline 
    resource directory is $local_resource/pipeline_resource/NAME where NAME
    is the name of the pipeline.'''
    def __init__(self, resource, resource_dir):
        self.resource = resource.split()
        if not resource_dir or type(resource_dir) != str:
            raise ValueError('Invalid resource directory {}'.format(resourece_dir))
        else:
            if os.path.isabs(os.path.expanduser(resource_dir)):
                self.pipeline_resource = os.path.expanduser(resource_dir)
            else:
                self.pipeline_resource = os.path.join(os.path.expanduser(
                    env.local_resource), 'pipeline_resource', resource_dir)
        try:
            if not os.path.isdir(self.pipeline_resource):
                os.makedirs(self.pipeline_resource)
        except:
            raise RuntimeError('Failed to create pipeline resource directory '
                .format(self.pipeline_resource))

    def __call__(self, ifiles):
        saved_dir = os.getcwd()
        os.chdir(self.pipeline_resource)
        lockfile = os.path.join(self.pipeline_resource, '.varianttools.lck')
        while True:
            if os.path.isfile(lockfile):
                env.logger.warning('Wait till file lock {} is release.'.format(lockfile))
                time.sleep(10)
            else:
                break
        # if there is no lock file, proceed
        # create a lock file
        open(lockfile, 'a').close()
        try:
            ofiles, md5files = self._downloadFiles(ifiles)
        finally:
            try:
                os.remove(lockfile)
            except Exception as e:
                env.logger.warning('Failed to remove lock file {}: e'
                    .format(lockfile, e))
        self._validate(md5files)
        os.chdir(saved_dir)
        return ofiles

    def _validate(self, md5_files):
        if md5_files:
            prog = ProgressBar('Validating md5 signature', sum([x[1] for x in md5_files]))
            mismatched_files = []
            for filename, s in md5_files:
                try:
                    downloaded_md5 = open(filename + '.md5').readline().split()[0]
                    calculated_md5 = calculateMD5(filename, partial=False)
                    if downloaded_md5 != calculated_md5:
                        mismatched_files.append(filename)
                except Exception as e:
                    env.logger.warning('Failed to verify md5 signature of {}: {}'
                        .format(filename[:-4], e))
                prog.update(prog.count + s)
            prog.done()
            if mismatched_files:
                env.logger.warning('md5 signature of {} mismatch. '
                      'Please remove {} and try again.'
                      .format(', '.join(mismatched_files),
                      'this file' if len(mismatched_files) == 1 else 'these files'))

    def _downloadFiles(self, ifiles):
        '''Download resource'''
        # decompress all .gz files
        skipped = []
        md5_files = []
        for cnt, URL in enumerate(sorted(self.resource)):
            filename = URL.rsplit('/', 1)[-1]
            dest_file = os.path.join(self.pipeline_resource, filename)
            try:
                if os.path.isfile(dest_file):
                    skipped.append(filename)
                else:
                    downloadURL(URL, dest_file, False,
                        message='{}/{} {}'.format(cnt+1, len(self.resource), filename))
            except KeyboardInterrupt as e:
                raise e
            except Exception as e:
                env.logger.error('Failed to download {}: {} {}'
                    .format(filename, type(e).__name__, e))
            #
            if filename.endswith('.gz') and not filename.endswith('tar.gz'):
                if not existAndNewerThan(ofiles=filename[:-3], ifiles=filename):
                    s = delayedAction(env.logger.info,
                        'Decompressing {}'.format(filename))
                    decompressGzFile(filename, inplace=False, force=True)
                    del s
            #
            if filename.endswith('.md5') and os.path.isfile(filename[:-4]):
                md5_files.append([filename[:-4], os.path.getsize(filename[:-4])])
        if skipped:
            env.logger.info('Using {} existing resource files under {}.'
                .format(len(skipped), self.pipeline_resource))
        return ifiles, md5_files
 

class Pipeline:
    def __init__(self, name, extra_args=[]):
        self.pipeline = PipelineDescription(name, extra_args)

    def substitute(self, text, VARS):
        # now, find ${}
        pieces = re.split('(\${[^}]*})', text)
        for idx, piece in enumerate(pieces):
            if piece.startswith('${') and piece.endswith('}'):
                KEY = piece[2:-1].lower()
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

    def execute(self, pname, input_files=[], output_files=[], jobs=1):
        if pname is None:
            if len(self.pipeline.pipelines) == 1:
                pname = self.pipeline.pipelines.keys()[0]
            else:
                raise ValueError('Name of pipeline should be specified because '
                    '{}.pipeline defines more than one pipelines. '
                    'Available pipelines are: {}.'.format(self.pipeline.name,
                    ', '.join(self.pipeline.pipelines.keys())))
        else:
            if pname not in self.pipeline.pipelines.keys():
                raise ValueError('Pipeline {} is undefined in configuraiton file '
                    '{}. Available pipelines are: {}'.format(pname,
                    self.pipeline.name, ', '.join(self.pipeline.pipelines.keys())))
        psteps = self.pipeline.pipelines[pname]
        #
        global max_running_jobs 
        max_running_jobs = jobs
        VARS = {
            'cmd_input': input_files,
            'cmd_output': output_files,
            'temp_dir': env.temp_dir,
            'cache_dir': env.cache_dir,
            'local_resource': env.local_resource,
        }
        for key, val in self.pipeline.pipeline_vars.items():
            VARS[key.lower()] = self.substitute(val, VARS)
        #
        ifiles = input_files
        for command in psteps:
            env.logger.info('Executing step {}_{} of pipeline {}: {}'
                .format(pname, command.index, self.pipeline.name,
                    command.comment))
            # substitute ${} variables
            if command.input:
                step_input = [x.strip() for x in self.substitute(command.input, VARS).split(',')]
            else:
                step_input = ifiles
            # if there is no input file?
            if not step_input:
                raise RuntimeError('Pipeline stops at step {}_{}: No input file is available.'
                    .format(pname, command.index))
            #
            VARS['input{}'.format(command.index)] = step_input
            env.logger.debug('INPUT of step {}_{}: {}'
                    .format(pname, command.index, step_input))
            # 
            # now, group input files
            if not command.input_emitter:
                emitter = EmitInput()
            else:
                try:
                    # remove ${INPUT} because it is determined by the emitter
                    if 'input' in VARS:
                        VARS.pop('input')
                    # ${CMD_INPUT} etc can be used.
                    emitter = eval(self.substitute(command.input_emitter, VARS))
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
                    VARS['input'] = ig
                    action = self.substitute(command.action, VARS)
                    env.logger.debug('Emitted input of step {}_{}: {}'
                        .format(pname, command.index, ig))
                    env.logger.debug('Action of step {}_{}: {}'
                        .format(pname, command.index, action))
                    action = eval(action)
                    if type(action) == tuple:
                        action = SequentialActions(action)
                    ofiles = action(ig)
                    if type(ofiles) == str:
                        step_output.append(ofiles)
                    else:
                        step_output.extend(ofiles)
                # wait for all pending jobs to finish
                wait_all()
                VARS['output{}'.format(command.index)] = step_output
                env.logger.debug('OUTPUT of step {}_{}: {}'
                    .format(pname, command.index, step_output))
                for f in step_output:
                    if not os.path.isfile(f):
                        raise RuntimeError('Output file {} does not exist after '
                            'completion of step {}_{}'
                            .format(f, pname, command.index))
            except Exception as e:
                raise RuntimeError('Failed to execute step {}_{}: {}'
                    .format(pname, command.index, e))
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


def executeArguments(parser):
    parser.add_argument('pipeline', nargs='+', metavar='PIPELINE/QUERY',
        help='''Name of a pipeline configuration file with optional names of
            pipelines to be executed if the configuration file defines more
            than one pipelines. The configuration file can be identified by
            path to a .pipeline file (with or without extension), or one
            of the online pipelines listed by command "vtools show pipelines".
            To keep backward compatibility, this option also accept a SQL query
            that will be executed, with project genotype database attached
            as "genotype" and annotation databases attached by their names.''')
    parser.add_argument('-i', '--input', nargs='*', metavar='OUTPUT_FILE',
        help='''Input files to the pipeline, which will be passed to the
            pipelines as pipeline variable ${CMD_INPUT}.''')
    parser.add_argument('-o', '--output', nargs='*', metavar='OUTPUT_FILE',
        help='''Names of output files of the pipeline, which will be passed to
            the pipelines as ${CMD_OUTPUT}.''')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs.''')
    parser.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter used to output results of a SQL query.''')

def execute(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # old usage with a SQL query? The pipeline interface should
            # have input files, should have option --output, and should not
            # specify delimiter, and the input file should exist.
            if not args.input or not args.output or args.delimiter != '\t':
                # if there is no output, 
                proj.db.attach('{}_genotype'.format(proj.name), 'genotype')
                # for backward compatibility
                proj.db.attach('{}_genotype'.format(proj.name))
                cur = proj.db.cursor()
                # 
                query = ' '.join(args.pipeline)
                if query.upper().startswith('SELECT'):
                    env.logger.debug('Analyze statement: "{}"'.format(query))
                    cur.execute('EXPLAIN QUERY PLAN ' + query)
                    for rec in cur:
                        env.logger.debug('\t'.join([str(x) for x in rec]))
                # really execute the query
                cur.execute(query)
                proj.db.commit()
                sep = args.delimiter
                for rec in cur:
                    print(sep.join(['{}'.format(x) for x in rec]))
            else:
                pipeline = Pipeline(args.pipeline[0], extra_args=args.unknown_args)
                # unspecified
                if len(args.pipeline) == 1:
                    pipeline.execute(None, args.input, args.output,
                        args.jobs)
                else:
                    for name in args.pipeline[1:]:
                        pipeline.execute(name, args.input, args.output,
                            args.jobs)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)

if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools
    pass
