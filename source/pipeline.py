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
import shlex
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
import csv
import platform
import logging
import random
from collections import namedtuple

from .utils import env, ProgressBar, downloadURL, calculateMD5, delayedAction, \
    existAndNewerThan, TEMP, decompressGzFile, typeOfValues, validFieldName, \
    FileInfo, convertDoubleQuote, openFile, encodeTableName, expandRegions
    
from .project import PipelineDescription, Project

if sys.version_info.major == 2:
    from ucsctools_py2 import showTrack
else:
    from ucsctools_py3 import showTrack

try:
    import pysam
    hasPySam = True
except (ImportError, ValueError) as e:
    hasPySam = False

# define a few functions that can be used in the lambda function 
# of a pipeline to display information about a project
def projInfo(tables=[], samples=[], format_string=''):
    '''Obtain information about one or more tables or samples and output 
    them according to a format string (using python string format syntax for
    a dictionary). The table description contains keys name, num, date,
    cmd, and desc. The sample description contains keys name and num.
    '''
    try:
        with Project(verbosity='0') as proj:
            res = []
            if tables:
                if not format_string:
                    format_string = '{name:20}: {num:12} ({desc})'
                if isinstance(tables, str):
                    tables = [tables]
                for table in tables:
                    table_name = encodeTableName(table)
                    desc, date, cmd = proj.descriptionOfTable(table_name)
                    res.append(format_string.format(name=table, 
                        num=proj.db.numOfRows(table_name),
                        desc=desc, date=date, cmd=cmd))
            if samples:
                proj.db.attach('{}_genotype'.format(proj.name))
                cur = proj.db.cursor()
                # get sample information
                if not format_string:
                    format_string = '{name:20}: {num:12}'
                #
                if isinstance(samples, str):
                    samples = [samples]
                for sample in samples:
                    IDs = proj.selectSampleByPhenotype("sample_name = '{}'".format(sample))
                    if len(IDs) == 0:
                        raise ValueError('No sample with name {} is located'.format(sample))
                    elif len(IDs) > 1:
                        raise ValueError('Name {} matches multiple samples.'.format(sample))
                    cur.execute('SELECT count(*) FROM {}_genotype.genotype_{};'.format(proj.name, IDs[0]))
                    numGenotypes = cur.fetchone()[0] 
                    res.append(format_string.format(name=sample,
                        num=numGenotypes))
            return r'\n'.join(res)
    except Exception as e:
        env.logger.debug('projInfo(tables={},samples={},format_string="{}") failed: {}'.format(tables, samples,format_string, e))
        return ''

class EmitInput:
    '''Select input files of certain types, group them, and send input files
    to action. Selection criteria can be True (all input file types, default),
    'False' (select no input file), 'fastq' (check content of files), or one or
    more file extensions (e.g. ['.sam', '.bam']).  Eligible files are by default
    sent altogether (group_by='all') to action (${INPUT} equals ${INPUT#} where
    # is the index of step, but can also be sent individually (group_by='single',
    ${INPUT} equals to a list of a single file) or in pairs 
    (group_by='paired', e.g. filename_1.txt and filename_2.txt). Unselected
    files are by default passed directly as output of a step.'''
    def __init__(self, group_by='all', select=True, pass_unselected=True):
        self.group_by = group_by
        if type(select) == str:
            if select not in ['fastq', 'bam', 'sam'] and not str(select).startswith('.'):
                raise ValueError("Value to option select can only be True/False, "
                    "'fastq', or a file extension with leading '.': '{}' provided."
                    .format(select))
            self.select = [select]
        elif select in [True, False]:
            self.select = select
        else:
            for s in select:
                if s not in ['fastq', 'bam', 'sam'] and not str(s).startswith('.'):
                    raise ValueError("Value to option select can only be True/False, "
                        "'fastq', or a file extension with leading '.': '{}' provided."
                        .format(s))
            self.select = select
        self.pass_unselected = pass_unselected

    def _isFastq(self, filename):
        try:
            if not os.path.isfile(filename) and not os.path.isfile(filename + '.file_info'):
                raise RuntimeError('File not found')
            fl = FileInfo(filename).firstline()
            if fl is None:
                env.logger.info('Cannot detect the type of file because the {} has been removed.'
                    .format(filename))
                return filename.lower().split('.')[-1] not in ['bam', 'sam', 'gz', 'zip']
            if not fl.startswith('@'):
                return False
            if filename.endswith('.gz'):
                env.logger.warning('{}: compressed fastq file might not be '
                    'acceptable to downstream analysis.'.format(filename))
        except Exception as e:
            env.logger.debug('Input file {} is not in fastq format: {}'.format(filename, e))
            return False
        return True

    def _is_paired(self, f1, f2, at=None):
        if len(f1) != len(f2):
            return False
        if f1 >= f2:
            return False
        diffs = [x != y for x,y in zip(f1, f2)]
        if sum(diffs) != 1:
            return False
        diff_at = diffs.index(True)
        if sorted([f1[diff_at], f2[diff_at]]) != ['1', '2']:
            return False
        if at is not None and diff_at not in at:
            return False
        return True

    def _pairByReadNames(self, selected, unselected):
        # we should pair files by actual read names
        read_map = {}
        for filename in selected:
            read = FileInfo(filename).firstline().strip()
            if read[:-1] in read_map:
                if read.endswith('1'):
                    if read_map[read[:-1]][0] is not None:
                        raise RuntimeError('Fastq file {} has the same first read as {}'
                            .format(filename, read_map[read[:-1]][0]))
                    else:
                        read_map[read[:-1]][0] = filename
                elif read.endswith('2'):
                    if read_map[read[:-1]][1] is not None:
                        raise RuntimeError('Fastq file {} has the same first read as {}'
                            .format(filename, read_map[read[:-1]][1]))
                    else:
                        read_map[read[:-1]][1] = filename
                else:
                    raise RuntimeError('Fastq file {} is not paired because its read name does '
                        'not end with 1 or 2'.format(filename))
            else:
                if read.endswith('1'):
                    read_map[read[:-1]] = [filename, None]
                elif read.endswith('2'):
                    read_map[read[:-1]] = [None, filename]
                else:
                    raise RuntimeError('Fastq file {} is not paired because its read name does '
                        'not end with 1 or 2'.format(filename))
        # now, let us go through files
        pairs = []
        for read, filenames in read_map.items():
            if filenames[0] is None:
                raise RuntimeError('Fastq file {} is not paired (no matching read is found)'
                    .format(filenames[0]))
            elif filenames[1] is None:
                raise RuntimeError('Fastq file {} is not paired (no matching read is found)'
                    .format(filenames[1]))
            else:
                if not self._is_paired(filenames[0], filenames[1]):
                    env.logger.warning('{} and {} contain paired reads but the filenames '
                        'do not follow illumina filename convention'
                        .format(filenames[0], filenames[1]))
                pairs.append(filenames)
        return sorted(pairs), unselected     

    def _pairByFileName(self, selected, unselected):
            #
            # there is a possibility that one name differ at multiple parts
            # with another name. e.g
            #
            #      A1_TAGCTT_L007_R1_001.fastq.gz
            #
            # differ with the following two names by a number
            #
            #      A1_TAGCTT_L007_R1_002.fastq.gz
            #      A1_TAGCTT_L007_R2_001.fastq.gz
            # 
            # the code below tries to find good pairs first, then use matched
            # locations to pair others
            if not selected:
                env.logger.warning('No file matching type "{}" is selected for pairing.'
                    .format(self.select))
                return [], unselected
            all_pairs = [[x,y] for x in selected for y in selected if self._is_paired(x,y)]
            unpaired = [x for x in selected if not any([x in y for y in all_pairs])]
            if unpaired:
                raise ValueError('Failed to pair input filenames: {} is not paired'
                    'with any other names.'.format(', '.join(unpaired)))
            uniquely_paired = [x for x in selected if sum([x in y for y in all_pairs]) == 1]
            # if some filenames are uniquely paired, we can use them to identify
            # index locations.
            if uniquely_paired:
                pairs = [x for x in all_pairs if x[0] in uniquely_paired or x[1] in uniquely_paired]
                if len(pairs) != all_pairs:
                    # find the differentiating index of existing pairs
                    diff_at = set([[i != j for i,j in zip(x[0],x[1])].index(True) for x in pairs])
                    # use the diff_at locations to screen the rest of the pairs
                    pairs.extend([x for x in all_pairs if x not in pairs and self._is_paired(x[0], x[1], diff_at)])
                    #
                    if len(pairs) * 2 != len(selected):
                        unpaired = [x for x in selected if not any([x in y for y in pairs])]
                        raise ValueError('Failed to pair input files because they '
                            'match multiple filenames: {}'.format(', '.join(unpaired)))
                return sorted(pairs), unselected
            else:
                # all filenames match to multiple names, so we try to get all
                # differentiating indexes and see if one of them can pair filenames
                # perfectly. We start from the end because we assume that _1 _2
                # are close to the end of filenames.
                #
                diff_at = set([[i != j for i,j in zip(x[0],x[1])].index(True) for x in all_pairs])
                acceptable_diff_at = []
                for d in diff_at:
                    # try to pair all names at this location.
                    pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], [d])]
                    if len(pairs) * 2 != len(selected):
                        continue
                    # all filename should appear once and only once
                    if not all([sum([x in y for y in pairs]) == 1 for x in selected]):
                        continue
                    acceptable_diff_at.append(d)
                # fortunately, only one perfect pairing is found
                if len(acceptable_diff_at) == 1:
                    pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], acceptable_diff_at)]
                    return sorted(pairs), unselected
                elif len(acceptable_diff_at) > 1:
                    env.logger.warning('There are {} ways to match all filenames '
                        'perfectly. The one using a latter differentiating index '
                        'is used.'.format(len(acceptable_diff_at)))
                    diff_at = sorted(list(acceptable_diff_at))[-1]
                    pairs = [x for x in all_pairs if self._is_paired(x[0], x[1], [diff_at])]
                    return sorted(pairs), unselected
                else:
                    raise ValueError('All filenames match multiple names but no differentiating '
                        'index can pair filenames perfectly.') 

    def __call__(self, ifiles, pipeline=None):
        selected = []
        unselected = []
        for filename in ifiles:
            match = False
            if self.select is True:
                match = True
            elif self.select is False:
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
        # for this special case, the step is skipped
        if self.group_by == 'single':
            return [[x] for x in selected], unselected
        elif self.group_by == 'all':
            return [selected], unselected
        elif self.group_by == 'paired':
            if 'fastq' in self.select:
                try:
                    return self._pairByReadNames(selected, unselected)
                except Exception as e:
                    # if failed to pair by read name, pair by filenames
                    env.logger.warning('Failed to pair fastq files by read names. '
                        'Trying to pair files by filenames: {}'.format(e))
            else:
                # this should not happen becase we do not need to pair non-fastq files 
                # at this point, but I will leave the code here anyway.
                env.logger.warning('It is unsafe to pair input files by names instead of '
                    'their content. Please add option select="fastq" if you need to '
                    'pair input fastq files')
            return self._pairByFileName(selected, unselected)


class SkiptableAction:
    '''Internal actions that can be skipped according to a exe_info file.
    Actions using this class as a base class should define a _execute() function.
    '''
    def __init__(self, cmd, output, ignoreInput=False):
        self.cmd = cmd.strip()
        if isinstance(output, str):
            self.output = [output]
        else:
            self.output = output
        self.ignoreInput = ignoreInput

    def __call__(self, ifiles, pipeline=None):
        exe_info = '{}.exe_info'.format(self.output[0])
        if os.path.isfile(exe_info) and open(exe_info).readline().strip() == self.cmd.strip() \
            and existAndNewerThan(self.output, [] if self.ignoreInput else ifiles,
            exe_info):
            # not that we do not care input file because it might contain different seed
            env.logger.info('Reuse existing {}'.format(self.output[0]))
            return self.output
        with open(exe_info, 'w') as exe_info:
            exe_info.write(self.cmd + '\n')
            exe_info.write('#Start: {}\n'.format(time.asctime(time.localtime())))
            for f in ifiles:
                # for performance considerations, use partial MD5
                exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                    calculateMD5(f, partial=True)))
            self._execute(ifiles, pipeline)
            exe_info.write('#End: {}\n'.format(time.asctime(time.localtime())))
            for f in self.output:
                if not os.path.isfile(f):
                    raise RuntimeError('Output file {} does not exist after completion of the job.'.format(f))
                # for performance considerations, use partial MD5
                exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                    calculateMD5(f, partial=True)))
        return self.output
        
try:
    # some functors are subclasses of SkiptableAction and has to be
    # imported after the definition of that class.
    from simulation import *
except ImportError as e:
    # The simulation functors will not be available if simulation module cannot
    # be loaded.
    env.logger.debug('Failed to import module simulation')



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

    def __call__(self, ifiles, pipeline=None):
        for a in self.actions:
            # the input of the next action is the output of the
            # previous action.
            ifiles = a(ifiles)
        # return the output of the last action
        return ifiles

class OutputText:
    '''Write its input to standard output, or a file if a filename is specified.
    The text can be a list of strings. A new line is added automatically to
    each line of the text.
    '''
    def __init__(self, text='', filename=None, mode='a'):
        if not isinstance(text, str):
            self.text = ''.join([str(x) + '\n' for x in text])
        else:
            self.text = text + '\n'
        self.filename = filename
        self.mode = mode

    def __call__(self, ifiles, pipeline=None):
        if self.filename is not None:
            with open(self.filename, self.mode) as output:
                output.write(self.text)
        else:
            sys.stdout.write(self.text)        
        return ifiles

class CheckCommands:
    '''Check the existence of specified commands and raise an error if one of
    the commands does not exist.'''
    def __init__(self, cmds):
        if type(cmds) == type(''):
            self.cmd = [cmds]
        else:
            self.cmd = cmds

    # the following is copied from shutils.which from Python 3.3
    def _which(self, cmd, mode=os.F_OK | os.X_OK, path=None):
        """Given a command, mode, and a PATH string, return the path which
        conforms to the given mode on the PATH, or None if there is no such
        file.

        `mode` defaults to os.F_OK | os.X_OK. `path` defaults to the result
        of os.environ.get("PATH"), or can be overridden with a custom search
        path.

        """
        # Check that a given file can be accessed with the correct mode.
        # Additionally check that `file` is not a directory, as on Windows
        # directories pass the os.access check.
        def _access_check(fn, mode):
            return (os.path.exists(fn) and os.access(fn, mode)
                    and not os.path.isdir(fn))

        # Short circuit. If we're given a full path which matches the mode
        # and it exists, we're done here.
        if _access_check(cmd, mode):
            return cmd

        path = (path or os.environ.get("PATH", os.defpath)).split(os.pathsep)

        if sys.platform == "win32":
            # The current directory takes precedence on Windows.
            if not os.curdir in path:
                path.insert(0, os.curdir)

            # PATHEXT is necessary to check on Windows.
            pathext = os.environ.get("PATHEXT", "").split(os.pathsep)
            # See if the given file matches any of the expected path extensions.
            # This will allow us to short circuit when given "python.exe".
            matches = [cmd for ext in pathext if cmd.lower().endswith(ext.lower())]
            # If it does match, only test that one, otherwise we have to try
            # others.
            files = [cmd] if matches else [cmd + ext.lower() for ext in pathext]
        else:
            # On other platforms you don't have things like PATHEXT to tell you
            # what file suffixes are executable, so just pass on cmd as-is.
            files = [cmd]

        seen = set()
        for dir in path:
            dir = os.path.normcase(dir)
            if not dir in seen:
                seen.add(dir)
                for thefile in files:
                    name = os.path.join(dir, thefile)
                    if _access_check(name, mode):
                        return name
        return None

    def __call__(self, ifiles, pipeline=None):
        for cmd in self.cmd:
            if self._which(cmd) is None:
                raise RuntimeError('Command {} does not exist. Please install it and try again.'
                    .format(cmd))
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

    def __call__(self, ifiles, pipeline=None):
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
    def __init__(self, files, msg=''):
        if type(files) == str:
            self.files = [files]
        else:
            self.files = files
        self.msg = msg

    def __call__(self, ifiles, pipeline=None):
        for f in self.files:
            if os.path.isfile(f):
                env.logger.info('{} is located.'.format(f))
            else:
                raise RuntimeError('Cannot locate {}. {}'.format(f, self.msg))
        return ifiles

class CheckDirs:
    '''Check the existence of specified directories and raise an
    error if one of the directories does not exist.'''
    def __init__(self, dirs, msg=''):
        if type(dirs) == str:
            self.dirs = [dirs]
        else:
            self.dirs = dirs
        self.msg = msg

    def __call__(self, ifiles, pipeline=None):
        for d in self.dirs:
            if os.path.isdir(d):
                env.logger.info('Directory {} is located.'.format(d))
            else:
                raise RuntimeError('Cannot locate directory {}. {}'.format(d, self.msg))
        return ifiles

class CheckVariantToolsVersion:
    def __init__(self, version=''):
        self.min_version = version

    def __call__(self, ifiles, pipeline=None):
        vtools_version = [int(x) for x in re.sub('\D', ' ', pipeline.VARS['vtools_version']).split()]
        # e.g. minimal 2.2.0, vtools 2.1.1
        if [int(x) for x in re.sub('\D', ' ', self.min_version).split()] > vtools_version:
            raise RuntimeError('Version {} is required to execute this pipeline. '
                'Please upgrade your installation of variant tools (version {})'
                .format(self.min_version, pipeline.VARS['vtools_version']))
        return ifiles
        
class CheckFastqVersion:
    def __init__(self, output):
        self.output = output

    def __call__(self, fastq_file, pipeline=None):
        '''Detect the version of input fastq file. This can be very inaccurate'''
        if not os.path.isfile(fastq_file[0]) and os.path.isfile(fastq_file[0] + '.file_info'):
            if os.path.isfile(self.output):
                return self.output
            else:
                raise RuntimeError('A valid fastq file is needed to check version of fastq: .file_info detected')
        with open(self.output, 'w') as aln_param:
            #
            # This function assumes each read take 4 lines, and the last line contains
            # quality code. It collects about 1000 quality code and check their range,
            # and use it to determine if it is Illumina 1.3+
            #
            qual_scores = ''
            with openFile(fastq_file[0]) as fastq:
                while len(qual_scores) < 1000:
                    try:
                        line = fastq.readline().decode('utf-8')
                    except Exception as e:
                        raise RuntimeError('Failed to read fastq file {}: {}'
                            .format(fastq_file, e))
                    if not line.startswith('@'):
                        raise ValueError('Wrong FASTA file {}'.format(fastq_file))
                    line = fastq.readline().decode('utf-8')
                    line = fastq.readline().decode('utf-8')
                    if not line.startswith('+'):
                        env.logger.warning(
                            'Suspiciout FASTA file {}: third line does not start with "+".'
                            .foramt(fastq_file))
                        return 
                    line = fastq.readline().decode('utf-8')
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


class FieldsFromTextFile:
    '''Read a text file, guess its delimeter, field name (from header)
    and create field descriptions. If a vcf file is encountered, all
    fields will be exported'''
    
    def __init__(self, output):
        self.field_output = output

    def __call__(self, ifiles, pipeline=None):
        try:
            if ifiles[0].endswith('.vcf') or ifiles[0].endswith('.vcf.gz'):
                showTrack(ifiles[0], self.field_output)
            else:
                with open(self.field_output, 'w') as fo:
                    csv_dialect = csv.Sniffer().sniff(open(ifiles[0], 'rU').read(8192))
                    fo.write('delimiter="{}"\n\n'.format(csv_dialect.delimiter.replace('\t', r'\t')))
                    values = []
                    with open(ifiles[0], 'rU') as fi:
                        reader = csv.reader(fi, dialect=csv_dialect)
                        headers = reader.next()
                        values = [[] for x in headers]
                        for line in reader:
                            for idx in range(len(headers)):
                                values[idx].append(line[idx])
                            if len(values[0]) > 100:
                                break
                    #
                    for idx, header in enumerate(headers):
                        fo.write('[{}]\n'.format(validFieldName(header)))
                        fo.write('index={}\n'.format(idx+1))
                        fo.write('type={}\n\n'.format(typeOfValues(values[idx])))
        except Exception as e:
            raise RuntimeError('Failed to guess fields from {}: {}'.format(ifiles[0], e))
        #
        return self.field_output
        
    
class GuessReadGroup:
    def __init__(self, bamfile, rgfile):
        self.output_bam = bamfile
        self.rg_output = rgfile

    def __call__(self, fastq_filename, pipeline=None):
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
JOB = namedtuple('JOB', 'proc cmd start_time stdout stderr working_dir output')
running_jobs = []
max_running_jobs = 1

def elapsed_time(start):
    '''Return the elapsed time in human readable format since start time'''
    second_elapsed = int(time.time() - start)
    days_elapsed = second_elapsed // 86400
    return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
        time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
 
def run_command(cmd, output=None, working_dir=None, max_jobs=1):
    '''Call a list of external command cmd, raise an error if any of them
    fails. '''
    global running_jobs
    if not output:
        # subprocess.DEVNULL was introduced in Python 3.3
        proc_out = None
        proc_err = None
        proc_lck = None
    else:
        proc_out = open(output[0] + '.out_{}'.format(os.getpid()), 'w')
        proc_err = open(output[0] + '.err_{}'.format(os.getpid()), 'w')
        proc_lck = output[0] + '.lck'
    #
    # wait for empty slot to run the job
    while True:
        if poll_jobs() >= min(max_jobs, max_running_jobs):
            time.sleep(5)
        else:
            break
    # there is a slot, start running
    if proc_lck:
        env.lock(proc_lck, str(os.getpid()))
    proc = subprocess.Popen(cmd[0], shell=True, stdout=proc_out, stderr=proc_err,
        cwd=working_dir)
    running_jobs.append(JOB(proc=proc, cmd=cmd,
        start_time=time.time(), stdout=proc_out, stderr=proc_err, 
        working_dir=working_dir, output=output))
    env.logger.info('Running [[{}]]'.format(cmd[0]))
    if proc_out is not None:
        env.logger.debug('Output redirected to {0}.out_{1} and {0}.err_{1} and '
            'will be saved to {0}.exe_info after completion of command.'
            .format(output[0], os.getpid()))

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
            # if a job is running, wait ... but not to the end
            time.sleep(10)
            continue
        #
        # job completed, close redirected stdout and stderr
        if len(job.cmd) == 1 and job.stdout is not None:
            job.stdout.close()
            job.stderr.close()
        #
        if ret < 0:
            if job.output:
                try:
                    env.unlock(job.output[0] + '.lck', str(os.getpid()))
                except:
                    pass
            raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
                .format(job.cmd[0], -ret, elapsed_time(job.start_time)))
        elif ret > 0:
            if job.output:
                with open(job.output[0] + '.err_{}'.format(os.getpid())) as err:
                    for line in err.read().split('\n')[-50:]:
                        env.logger.error(line)
                try:
                    env.unlock(job.output[0] + '.lck', str(os.getpid()))
                except:
                    pass
            raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
                .format(job.cmd[0], elapsed_time(job.start_time), ret))
        else:
            #if job.output:
            #    with open(job.output[0] + '.err_{}'.format(os.getpid())) as err:
            #        for line in err.read().split('\n')[-10:]:
            #            env.logger.info(line)
            env.logger.info('Command "{}" completed successfully in {}'
                .format(job.cmd[0], elapsed_time(job.start_time)))
            # 
            # if there are no more jobs, complete .exe_info
            if len(job.cmd) == 1:
                #
                if job.output:
                    # including output from previous failed runs
                    # this step will fail if the .lck file has been changed by
                    # another process
                    env.unlock(job.output[0] + '.lck', str(os.getpid()))
                    if not os.path.isfile(job.output[0] + '.out_{}'.format(os.getpid())) \
                        or not os.path.isfile(job.output[0] + '.err_{}'.format(os.getpid())):
                        env.logger.warning('Could not locate process-specific output file (id {}), '
                            'which might have been removed by another process that produce '
                            'the same output file {}'.format(os.getpid(), job.output[0]))
                        # try to rerun this step
                        running_jobs[idx] = None
                        raise RewindExecution()
                    with open(job.output[0] + '.exe_info', 'a') as exe_info:
                        exe_info.write('#End: {}\n'.format(time.asctime(time.localtime())))
                        for f in job.output:
                            if not os.path.isfile(f):
                                raise RuntimeError('Output file {} does not exist after completion of the job.'.format(f))
                            # for performance considerations, use partial MD5
                            exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                                calculateMD5(f, partial=True)))
                        # write standard output to exe_info
                        exe_info.write('\n\nSTDOUT\n\n')
                        with open(job.output[0] + '.out_{}'.format(os.getpid())) as stdout:
                            for line in stdout:
                                exe_info.write(line)
                        # write standard error to exe_info
                        exe_info.write('\n\nSTDERR\n\n')
                        with open(job.output[0] + '.err_{}'.format(os.getpid())) as stderr:
                            for line in stderr:
                                exe_info.write(line)
                    # if command succeed, remove all out_ and err_ files, 
                    for filename in glob.glob(job.output[0] + '.out_*') + \
                        glob.glob(job.output[0] + '.err_*'):
                        try:
                            os.remove(filename)
                        except Exception as e:
                            env.logger.warning('Fail to remove {}: {}'
                                .format(filename, e))
                #
                running_jobs[idx] = None
            else:
                # start the next job in the same slot
                proc = subprocess.Popen(job.cmd[1], shell=True, stdout=job.stdout,
                    stderr=job.stderr, cwd=job.working_dir)
                # use the same slot for the next job
                running_jobs[idx] = JOB(proc=proc, cmd=job.cmd[1:],
                    start_time=job.start_time, stdout=job.stdout,
                    stderr=job.stderr, working_dir=job.working_dir, output=job.output)
                env.logger.info('Running [[{}]]'.format(job.cmd[1]))
                # increase the running job count
                count += 1
    return count

def wait_all():
    '''Wait for all pending jobs to complete'''
    global running_jobs
    try:
        while True:
            if not running_jobs or all([x is None for x in running_jobs]):
                break
            # sleep one second, but poll_jobs will wait 10s for each running
            # job. This avoid waiting at least 10s for simple commands
            time.sleep(1)
            if poll_jobs() == 0:
                break
    except KeyboardInterrupt:
        # clean up lock files
        for job in running_jobs:
            if job is not None and job.stdout:
                try:
                    env.unlock(job.output[0] + '.lck')
                except:
                    pass
        # raise an error instead of exit right now to give vtools
        # a chance to close databases
        raise RuntimeError('Keyboard interrupted')
    running_jobs = []


class RewindExecution(Exception):
    pass

class NullAction:
    def __init__(self, output=[], action=''):
        '''A null action that is used to change input, output, or
        execute lambda function of substituted variable.'''
        if type(output) == str:
            self.output = [output]
        else:
            self.output = output

    def __call__(self, ifiles, pipeline=None):
        if self.output:
            return self.output
        else:
            # if no output is specified, pass ifiles through
            return ifiles
        
class RunCommand:
    def __init__(self, cmd='', working_dir=None, output=[], max_jobs=1,
        locking="exclusive"):
        '''This action execute the specified command under the
        specified working directory, and return specified ofiles.
        #
        # parameter locking is unused
        '''
        # merge mulit-line command into one line and remove extra white spaces
        if not cmd:
            self.cmd = ['echo "None command executed."']
        elif type(cmd) == str:
            self.cmd = [' '.join(cmd.split('\n'))]
        else:
            self.cmd = [' '.join(x.split('\n')) for x in cmd]
        self.max_jobs = max_jobs
        self.working_dir = working_dir
        if type(output) == str:
            self.output = [os.path.expanduser(output)]
        else:
            self.output = [os.path.expanduser(x) for x in output]
        #
        for filename in self.output:
            if filename.lower().rsplit('.', 1)[-1] in ['exe_info', 'lck']:
                raise RuntimeError('Output file with extension .exe_info and .lck are reserved.')

    def __call__(self, ifiles, pipeline=None):
        # substitute cmd by input_files and output_files
        if self.output:
            #for ofile in self.output:
            #    nfile = os.path.normpath(os.path.expanduser(ofile))
            #    under_resource = not os.path.relpath(nfile, os.path.normpath(os.path.expanduser(env.local_resource))).startswith('../')
            #    under_cache = not os.path.relpath(nfile, os.path.normpath(os.path.expanduser(env.cache_dir))).startswith('../')
            #    under_temp = not os.path.relpath(nfile, os.path.normpath(os.path.expanduser(env.temp_dir))).startswith('../')
            #    under_proj = not os.path.relpath(nfile, os.path.realpath('.')).startswith('../')
            #    if not (under_resource or under_cache or under_temp or under_proj):
            #        raise RuntimeError('Unable to write to output file ({}) because '
            #            'output file of a pipeline can can be under project, resource, '
            #            'temporary or cache directories.'.format(ofile))
            if os.path.isfile(self.output[0] + '.exe_info'):
                with open(self.output[0] + '.exe_info') as exe_info:
                    cmd = exe_info.readline().strip()
                # if the exact command has been used to produce output, and the
                # output files are newer than input file, ignore the step
                if cmd == '; '.join(self.cmd).strip() and existAndNewerThan(ifiles=ifiles,
                        ofiles=self.output, md5file=self.output[0] + '.exe_info'):
                    env.logger.info('Reuse existing files {}'.format(', '.join(self.output)))
                    return self.output
            # create directory if output directory does not exist
            for dir in [os.path.split(os.path.abspath(x))[0] for x in self.output]:
                if not os.path.isdir(dir):
                    try:
                        os.makedirs(dir)
                    except Exception as e:
                        raise RuntimeError('Failed to create directory {} for output file: {}'.format(dir, e))
                    if not os.path.isdir(dir):
                        raise RuntimeError('Failed to create directory {} for output file: {}'.format(dir, e))
        # now, we cannot ignore this step, but do we have all the input files?
        # the input can be fake .file_info files
        if not all([os.path.isfile(x) for x in ifiles]):
            raise RewindExecution()
        run_command(self.cmd, output=self.output, working_dir=self.working_dir,
            max_jobs=self.max_jobs)
        # add md5 signature of input and output files
        if self.output:
            with open(self.output[0] + '.exe_info', 'w') as exe_info:
                exe_info.write('{}\n'.format('; '.join(self.cmd)))
                exe_info.write('#Start: {}\n'.format(time.asctime(time.localtime())))
                for f in ifiles:
                    # for performance considerations, use partial MD5
                    exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                        calculateMD5(f, partial=True)))
        if self.output:
            return self.output
        else:
            # if no output is specified, pass ifiles through
            return ifiles


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
                with bz2.BZ2File(filename, 'rb') as bzinput, open(TEMP(dest_file), 'wb') as output:
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
                # only extract files
                files = [x.name for x in tar.getmembers() if x.isfile()]
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
        
    def __call__(self, ifiles, pipeline=None):
        # decompress input files and return a list of output files
        filenames = []
        for filename in ifiles:
            filenames.extend(self._decompress(filename))
        filenames.sort()
        return filenames


class RemoveIntermediateFiles:
    def __init__(self, files):
        self.files = files

    def __call__(self, ifiles, pipeline=None):
        env.logger.debug('Remove intermediate files {}'.format(self.files))
        for f in shlex.split(self.files) if isinstance(self.files, str) else self.files:
            if not os.path.isfile(f):
                if os.path.isfile(f + '.file_info'):
                    env.logger.info('Reusing existing {}.file_info.'.format(f))
                else:
                    raise RuntimeError('Failed to create {}.file_info: Missing input file.'
                        .format(f))
            else:
                FileInfo(f).save()
                env.logger.info('Replace {0} with {0}.file_info'.format(f))
                try:
                    os.remove(f)
                except e:
                    env.logger.warning('Failed to remove intermediate file {}'
                        .format(f))
        return ifiles


class LinkToDir:
    '''Create hard links of input files to a specified directory. This is 
    usually used to link input files to a common cache directory so that 
    all operations can be performed on that directory.'''
    def __init__(self, dest_dir):
        self.dest = dest_dir
        if not os.path.isdir(self.dest):
            env.logger.info('Creating directory {}'.format(self.dest))
            try:
                os.makedirs(self.dest)
            except Exception as e:
                raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))
            if not os.path.isdir(self.dest):
                raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))

    def __call__(self, ifiles, pipeline=None):
        ofiles = []
        for filename in ifiles:
            path, basename = os.path.split(filename)
            if not os.path.isfile(filename):
                if os.path.isfile(filename + '.file_info'):
                    dest_file = os.path.join(self.dest, basename) + '.file_info'
                    if os.path.isfile(dest_file):
                        if not os.path.samefile(filename + '.file_info', dest_file):
                            os.remove(dest_file)
                            env.logger.info('Linking {} to {}'.format(filename, self.dest))
                            os.link(filename + '.file_info', os.path.join(self.dest, basename) + '.file_info')
                        else:
                            env.logger.debug('Reusing existing linked file_info file: {}'
                                .format(os.path.join(self.dest, basename) + '.file_info'))
                    else:
                        env.logger.info('Linking {} to {}'.format(filename, self.dest))
                        os.link(filename + '.file_info', os.path.join(self.dest, basename) + '.file_info')
                else:
                    raise RuntimeError('Failed to link {} to directory {}: file does not exist'
                        .format(filename, self.dest))
            else:
                dest_file = os.path.join(self.dest, basename)
                if os.path.isfile(dest_file):
                    if not os.path.samefile(filename, dest_file):
                        os.remove(dest_file)
                        env.logger.info('Linking {} to {}'.format(filename, self.dest))
                        os.link(filename, dest_file)
                    else:
                        env.logger.debug('Reusing existing linked file: {}'.format(dest_file))
                else:
                    env.logger.info('Linking {} to {}'.format(filename, self.dest))
                    os.link(filename, dest_file)
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

    def __call__(self, sam_file, pipeline=None):
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
    '''Download resources to specified destination directory. dest_dir can
    be a full path name or a directory relative to 
    $local_resource/pipeline_resource where $local_resource is the local
    resource directory of the project (default to ~/.variant_tools,
    see runtime option local_resource for details). The default pipeline 
    resource directory is $local_resource/pipeline_resource/NAME where NAME
    is the name of the pipeline.'''
    def __init__(self, resource, dest_dir):
        self.resource = resource.split()
        if not dest_dir or type(dest_dir) != str:
            raise ValueError('Invalid resource directory {}'.format(dest_dir))
        else:
            if os.path.isabs(os.path.expanduser(dest_dir)):
                self.pipeline_resource = os.path.expanduser(dest_dir)
            else:
                self.pipeline_resource = os.path.join(os.path.expanduser(
                    env.local_resource), 'pipeline_resource', dest_dir)
        try:
            if not os.path.isdir(self.pipeline_resource):
                os.makedirs(self.pipeline_resource)
        except:
            raise RuntimeError('Failed to create pipeline resource directory '
                .format(self.pipeline_resource))

    def __call__(self, ifiles, pipeline=None):
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
        env.lock(lockfile)
        try:
            ofiles, md5files = self._downloadFiles(ifiles)
        finally:
            env.unlock(lockfile)
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
                raise RuntimeError('Failed to download {}: {} {}'
                    .format(filename, type(e).__name__, e))
            #
            if filename.endswith('.tar.gz'):
                manifest = filename + '.manifest'
                if not os.path.isfile(manifest):
                    with tarfile.open(filename, 'r:gz') as tar: 
                        s = delayedAction(env.logger.info, 'Extracting {}'.format(filename))
                        tar.extractall(self.pipeline_resource)
                        del s
                        # only extract files
                        files = [x.name for x in tar.getmembers() if x.isfile()]
                        # save content to a manifest
                        with open(manifest, 'w') as manifest:
                            for f in files:
                                manifest.write(f + '\n')
            elif filename.endswith('.gz'):
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
    def __init__(self, proj, name, extra_args=[]):
        self.proj = proj
        self.pipeline = PipelineDescription(name, extra_args)

    def var_expr(self, var):
        if type(var) == str:
            # tries to be clever and quote filenames with space
            if os.path.isfile(var) and ' ' in var:
                return "'{}'".format(var)
            else:
                return var
        elif type(var) == list:
            return ' '.join([self.var_expr(x) for x in var])
        else:
            return str(var)

    def _substitute(self, text, PipelineVars):
        # if text has new line, replace it with space
        text =  ' '.join(text.split())
        # now, find ${}
        pieces = re.split('(\${[^{}]*})', text)
        for idx, piece in enumerate(pieces):
            if piece.startswith('${') and piece.endswith('}'):
                KEY = piece[2:-1].lower()
                if ':' in KEY:
                    # a lambda function?
                    try:
                        FUNC = eval('lambda {}'.format(piece[2:-1]))
                    except Exception as e:
                        env.logger.warning('Failed to interpret {} as a pipeline variable: {}'
                            .format(piece, e))
                        continue
                    KEY = KEY.split(':', 1)[0].strip()
                    try:
                        if not KEY:
                            # if there is no KEY, this is a lamba function without parameter
                            pieces[idx] = self.var_expr(FUNC())
                        elif ',' not in KEY:
                            # single varialbe
                            if KEY in PipelineVars:
                                VAL = PipelineVars[KEY]
                            else:
                                env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
                                    .format(piece, KEY))
                                continue
                            pieces[idx] = self.var_expr(FUNC(VAL))
                        else:
                            # several parameters
                            KEYS = KEY.split(',')
                            VAL = []
                            for KEY in KEYS:
                                # single varialbe
                                if KEY in PipelineVars:
                                    VAL.append(PipelineVars[KEY])
                                else:
                                    env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
                                        .format(piece, KEY))
                                    continue
                            pieces[idx] = self.var_expr(FUNC(*VAL))
                    except Exception as e:
                        env.logger.warning('Failed to interpret {} as a pipeline variable: {}'
                            .format(piece, e))
                        continue
                else:
                    # if KEY in PipelineVars, replace it
                    if KEY in PipelineVars:
                        pieces[idx] = self.var_expr(PipelineVars[KEY])
                    else:
                        env.logger.warning('Failed to interpret {} as a pipeline variable: key "{}" not found'
                            .format(piece, KEY))
                        continue
        # now, join the pieces together, but remove all newlines
        return ' '.join(''.join(pieces).split())

    def substitute(self, text, PipelineVars):
        while True:
            new_text = self._substitute(text, PipelineVars)
            if new_text == text:
                return new_text
            else:
                text = new_text

    def execute(self, pname, input_files=[], output_files=[], jobs=1, **kwargs):
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
        # if there is a output file, write log to .log
        if output_files:
            logfile = output_files[0] + '.log'
            ch = logging.FileHandler(logfile.lstrip('>'), mode = 'a')
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
            env.logger.addHandler(ch)
        #
        global max_running_jobs 
        max_running_jobs = jobs
        self.VARS = {
            'cmd_input': input_files,
            'cmd_output': output_files,
            'temp_dir': env.temp_dir,
            'cache_dir': env.cache_dir,
            'local_resource': env.local_resource,
            'ref_genome_build': '' if self.proj is None or not self.proj.build else  self.proj.build,
            'vtools_version': '' if self.proj is None else self.proj.version,
        }
        self.VARS.update(**kwargs)
        for key, val in self.pipeline.pipeline_vars.items():
            self.VARS[key.lower()] = self.substitute(val, self.VARS)
        #
        ifiles = input_files
        step_index = 0
        rewind_count = 0
        while True:
            # step_index can jump back and forth depending on the 
            # execution status of each step
            command = psteps[step_index]
            env.logger.info('Executing [[{}.{}_{}]]: {}'
                .format(self.pipeline.name, pname, command.index, 
                    ' '.join(command.comment.split())))
            # substitute ${} variables
            if command.input:
                step_input = shlex.split(self.substitute(command.input, self.VARS))
            else:
                step_input = ifiles
            #
            # if there is no input file?
            if not step_input:
                raise RuntimeError('Pipeline stops at step {}_{}: No input file is available.'
                    .format(pname, command.index))
            #
            self.VARS['input{}'.format(command.index)] = step_input
            env.logger.debug('INPUT of step {}_{}: {}'
                    .format(pname, command.index, step_input))
            # 
            # now, group input files
            if not command.input_emitter:
                emitter = EmitInput()
            else:
                try:
                    # remove ${INPUT} because it is determined by the emitter
                    if 'input' in self.VARS:
                        self.VARS.pop('input')
                    # ${CMD_INPUT} etc can be used.
                    emitter = eval(self.substitute(command.input_emitter, self.VARS))
                except Exception as e:
                    raise RuntimeError('Failed to group input files: {}'
                        .format(e))
            #
            saved_dir = os.getcwd()
            # pass Pipeline itself to emitter
            igroups, step_output = emitter(step_input, self)
            try:
                for ig in igroups:
                    if not ig:
                        continue
                    self.VARS['input'] = ig
                    action = self.substitute(command.action, self.VARS)
                    env.logger.debug('Emitted input of step {}_{}: {}'
                        .format(pname, command.index, ig))
                    env.logger.debug('Action of step {}_{}: {}'
                        .format(pname, command.index, action))
                    action = eval(action)
                    if type(action) == tuple:
                        action = SequentialActions(action)
                    # pass the Pipeline object itself to action
                    # this allows the action to have access to pipeline variables
                    # and other options
                    ofiles = action(ig, self)
                    if type(ofiles) == str:
                        step_output.append(ofiles)
                    else:
                        step_output.extend(ofiles)
                # wait for all pending jobs to finish
                wait_all()
                self.VARS['output{}'.format(command.index)] = step_output
                env.logger.debug('OUTPUT of step {}_{}: {}'
                    .format(pname, command.index, step_output))
                for f in step_output:
                    if not (os.path.isfile(f) or os.path.isfile(f + '.file_info')):
                        raise RuntimeError('Output file {} does not exist after '
                            'completion of step {}_{}'
                            .format(f, pname, command.index))
                for key, val in command.pipeline_vars:
                    self.VARS[key.lower()] = self.substitute(val, self.VARS)
                    env.logger.debug('Pipeline variable {} is set to {}'
                        .format(key, self.VARS[key.lower()]))
                #
                ifiles = step_output
                # this step is successful, go to next
                os.chdir(saved_dir)
                step_index += 1
                env.logger.debug('Step {}.{}_{} is executed successfully.'
                    .format(self.pipeline.name, pname, command.index))
                if step_index == len(psteps):
                    break
            except RewindExecution:
                rewind_count += 1
                if rewind_count >= 3:
                    raise RuntimeError('Failed to execute pipeline {}.{}: excessive '
                        'rewind during execution.'.format(self.pipeline.name, pname))
                # unfortunately, a input file has been removed (replaced by .file_info) but
                # a later steps need it. We will have to figure out how to create this 
                # file by looking backward ...
                to_be_regenerated = [x for x in step_input if not os.path.isfile(x)]
                # remove all fony files so that they will be re-generated
                for x in to_be_regenerated:
                    if os.path.isfile(x + '.file_info'):
                        os.remove(x + '.file_info')
                remaining = [x for x in to_be_regenerated]
                while step_index > 0:
                    step_index -= 1
                    command = psteps[step_index]
                    # remove all fony files so that they will be re-generated
                    for x in self.VARS['input{}'.format(command.index)]:
                        if os.path.isfile(x + '.file_info'):
                            os.remove(x + '.file_info')
                    # if any of the input file does not exist, go back further
                    if not all([os.path.isfile(x) for x in self.VARS['input{}'.format(command.index)]]):
                        continue
                    # check if a real file can be generated at this step
                    remaining = [x for x in remaining if x not in self.VARS['output{}'.format(command.index)]]
                    if not remaining:
                        break
                if step_index > 1:
                    ifiles = self.VARS['output{}'.format(psteps[step_index - 1].index)]
                else:
                    ifiles = self.VARS['cmd_input']
                env.logger.warning('Rewinding to [[{}.{}_{}]]: input files {} need to be re-generated.'
                    .format(self.pipeline.name, pname, command.index, ', '.join(to_be_regenerated)))
                os.chdir(saved_dir)
            except Exception as e:
                env.logger.debug('Failed to execute step {}.{}_{}.'
                    .format(self.pipeline.name, pname, command.index))
                raise RuntimeError('Failed to execute step {}_{}: {}'
                    .format(pname, command.index, e))




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
            If no input and output files are specified (options --input and
            --output), values of this option is treated as a SQL query that
            will be executed against the project database, with project genotype
            database attached as "genotype" and annotation databases attached
            by their names.''')
    parser.add_argument('-i', '--input', nargs='*', metavar='INPUT', default=[],
        help='''Input to the pipelines, usually a list of input files, that
            will bepassed to the pipelines as variable ${CMD_INPUT}.''')
    parser.add_argument('-o', '--output', nargs='*', metavar='OUTPUT', default=[],
        help='''Output of the pipelines, usually a list of output files, that
            will be passed to the pipelines as variable ${CMD_OUTPUT}.''')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs to execute, for steps
            of a pipeline that allows multi-processing.''')
    parser.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter used to output results of a SQL query.''')

def execute(args):
    # to keep backward compatibility, the vtools execute command
    # can execute a SQL query and a pipeline
    def executeQuery(proj):
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
        try:
            cur.execute(query)
        except Exception as e:
            raise RuntimeError('Failed to execute SQL query "{}": {}'
                .format(query, e))
        proj.db.commit()
        sep = args.delimiter
        for rec in cur:
            print(sep.join(['{}'.format(x) for x in rec]))
    #
    def executePipeline(proj):                
        pipeline = Pipeline(proj, args.pipeline[0], extra_args=args.unknown_args)
        # unspecified
        if len(args.pipeline) == 1:
            pipeline.execute(None, args.input, args.output,
                args.jobs)
        else:
            for name in args.pipeline[1:]:
                pipeline.execute(name, args.input, args.output,
                    args.jobs)
    # 
    try:
        with Project(verbosity=args.verbosity, mode='ALLOW_NO_PROJ') as proj:
            #
            # definitely a pipeline
            if args.input or args.output or args.unknown_args:
                executePipeline(proj)
            # definitely a sql query
            elif args.delimiter != '\t':
                executeQuery(proj)
            else:
                try:
                    # try to execute as a SQL query
                    executeQuery(proj)
                except RuntimeError as e:
                    env.logger.debug('Failed to execute {} as SQL query: {}'
                        .format(' '.join(args.pipeline), e))
                    executePipeline(proj)
    except Exception as e:
        env.unlock_all()
        env.logger.error(e)
        sys.exit(1)



#
# vtools simulate is implemented using the pipeline execution mechanism. It 
# essentially execute a pipeline that calls various simulation functions 
# to simulate data.
#
def simulateArguments(parser):
    parser.add_argument('model', nargs='+', metavar='MODEL',
        help='''Name of a model to simulate. It should be the name of a model
            description file with optional names for specific models within the
            descrition file to simulate. A list of model-specific parameters could be 
            specified to change the behavior of these models. Please use command
            "vtools show simulations" to get a list all available simulation
            models and "vtools show simulation MODEL" for details of a particular
            simulation model.''')
    parser.add_argument('--seed', type=int,
        help='''Random seed for the simulation. A random seed will be used by
            default but a specific seed could be used to reproduce a previously
            executed simulation.''')
    
def simulate(args):
    try:
        with Project(verbosity=args.verbosity, mode='ALLOW_NO_PROJ') as proj:
            # step 1, create a simulation configuration file.
            model_name = os.path.basename(args.model[0]).split('.', 1)[0]
            if args.seed is None:
                args.seed = random.randint(1, 2**32-1)
            if not os.path.isdir(env.cache_dir):
                os.mkdir(env.cache_dir)
            if len(args.model) == 1:
                cfg_file = '{}/{}_{}.cfg'.format(env.cache_dir, model_name, args.seed)
            else:
                cfg_file = '{}/{}_{}_{}.cfg'.format(env.cache_dir, model_name, args.model[1], args.seed)
            #
            with open(cfg_file, 'w') as cfg:
                cfg.write('model={}\n'.format(' '.join(args.model)))
                cfg.write('seed={}\n'.format(args.seed))
                if '--seed' in sys.argv:
                    # skip the seed option so to stop pipeline from distinguishing the two commands
                    cmd_args = sys.argv[:sys.argv.index('--seed')] + sys.argv[sys.argv.index('--seed') + 2:]
                    cfg.write("command=vtools {}\n".format(subprocess.list2cmdline(cmd_args[1:])))
                else:
                    cfg.write("command={}\n".format(env.command_line))
            #
            env.logger.info('Starting simulation {}'.format(cfg_file))
            pipeline = Pipeline(proj, args.model[0], extra_args=args.unknown_args)
            if len(args.model) == 1:
                pipeline.execute(None, [cfg_file], [], 1, seed=args.seed)
            else:
                for name in args.model[1:]:
                    pipeline.execute(name, [cfg_file], [], 1, seed=args.seed)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


if __name__ == '__main__':
    pass
