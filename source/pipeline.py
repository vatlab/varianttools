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

'''
This module implement Variant Pipeline Tools and VPT-provided
pipeline actions.
'''
import os
import sys
import subprocess
# for parallel execution of steps
import Queue
import threading
#
import glob
import shlex
import argparse
import logging
import shutil
import tarfile
import gzip
import bz2
import zipfile
import time
import re
import csv
import platform
import logging
import random
from multiprocessing import Process
from collections import namedtuple, MutableMapping

from .utils import env, ProgressBar, downloadFile, downloadURL, calculateMD5, delayedAction, \
    existAndNewerThan, TEMP, decompressGzFile, typeOfValues, validFieldName, \
    FileInfo, convertDoubleQuote, openFile, encodeTableName, expandRegions, \
    substituteVars, which
    
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


class SkipIf:
    '''An input emitter that skips the step (does not pass any input to the
    action) if certain condition is met. The input will be passed directly
    as output by default. This emitter is equivalent to 
    ``EmitInput(select=not cond, pass_unselected)`` '''
    def __init__(self, cond=None, pass_unselected=True):
        '''Does not emit input and skip the step if cond is ``True``. 
        
        Parameters:
            cond (boolean): 
                A boolean value True or False. In practice ``cond`` is usually
                a lambda function that checks the existence of a file or value
                of a pipeline variable.

            pass_unselected (boolean):
                Pass input files to output if ``cond`` is not met.

        '''
        self.cond = cond
        self.pass_unselected = pass_unselected

    def __call__(self, ifiles, pipeline=None):
        if self.cond:
            return [[]], ifiles
        else:
            return [ifiles], []

class EmitInput:
    '''An input emitter that emits input files individually, in pairs, or 
    altogether.'''
    def __init__(self, group_by='all', select=True, pass_unselected=True):
        '''Select input files of certain types, group them, and send input files
        to action. Selection criteria can be True (all input file types, default),
        'False' (select no input file), 'fastq' (check content of files), or one or
        more file extensions (e.g. ['.sam', '.bam']).  Eligible files are by default
        sent altogether (group_by='all') to action (${INPUT} equals ${INPUT#} where
        # is the index of step, but can also be sent individually (group_by='single',
        ${INPUT} equals to a list of a single file) or in pairs 
        (group_by='paired', e.g. filename_1.txt and filename_2.txt). Unselected
        files are by default passed directly as output of a step.'''
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


class PipelineAction:
    '''Base class for all pipeline actions. If one or more output files
    are specified, the pipeline will record the runtime signature of
    this action in a file ``$OUTPUT[0].exe_info``, which consists of the
    MD5 signature of input and output files, command used, and additional
    information such as start and end time of execution, standard and error
    outputs. An action will be skipped if the action is re-run with the
    same input, output and command.

    NOTE: The ``__call__`` function or a pipeline action implements the 
    runtime signature feature and calls function ``execute`` for actual work.
    User-defined actions should either override the ``__call__`` function
    (without the runtime signature feature) or function 
    ``_execute(self, ifiles, pipeline)`` (with runtime signature feature).
    User can also define function ``_bypass(self, ifiles, pipeline)`` if the
    step is bypassed due to identical execution signatures.
    '''
    def __init__(self, cmd='', output=[]):
        '''
        Parameters:
            cmd (string or list of strings):
                one or more commands to be executed. It should capture
                the name and all options used by the command.

            output (string or list of strings):
                Output files. If at least one output file is specified,
                the runtime signature of this action will be saved to
                $output[0].exe_info.
        '''
        # multiple command is not allowed.
        if not cmd:
            self.cmd = []
        elif isinstance(cmd, str):
            self.cmd = [' '.join(cmd.split('\n'))]
        else:
            self.cmd = [' '.join(x.split('\n')) for x in cmd]
        #
        if not output:
            self.output = []
        elif isinstance(output, str):
            self.output = [output]
        else:
            self.output = output
        #
        if self.output:
            self.proc_out = '{}.out_{}'.format(self.output[0], os.getpid())
            self.proc_err = '{}.err_{}'.format(self.output[0], os.getpid())
            self.proc_lck = '{}.lck'.format(self.output[0])
            self.proc_info = '{}.exe_info'.format(self.output[0])

    def _bypass(self, ifiles, pipeline=None):
        '''Function called by ``__call__`` if the step is bypassed due to identical
        execution signature. This function can be used, for example, to set pipeline
        variable even when the step is not executed.'''
        return True

    def _execute(self, ifiles, pipeline=None):
        '''Function called by ``__call__`` for actual action performed on ifiles. A user-defined
        action should re-define __call__ or this function. This funciton should return ``True`` if
        the action is completed successfully, ``False`` for pending (signature will be written later,
        and raise an exception for errors. '''
        raise RuntimeError('Please define your own execute function in an derived class of PipelineAction.')

    def _write_info(self):
        if not self.output:
            return
        with open(self.proc_info, 'a') as exe_info:
            exe_info.write('#End: {}\n'.format(time.asctime(time.localtime())))
            for f in self.output:
                if not os.path.isfile(f):
                    raise RuntimeError('Output file {} does not exist after completion of the job.'.format(f))
                # for performance considerations, use partial MD5
                exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                    calculateMD5(f, partial=True)))
            # write standard output to exe_info
            exe_info.write('\n\nSTDOUT\n\n')
            if os.path.isfile(self.proc_out):
                with open(self.proc_out) as stdout:
                    for line in stdout:
                        exe_info.write(line)
            # write standard error to exe_info
            exe_info.write('\n\nSTDERR\n\n')
            if os.path.isfile(self.proc_err):
                with open(self.proc_err) as stderr:
                    for line in stderr:
                        exe_info.write(line)
        # if command succeed, remove all out_ and err_ files, 
        for filename in glob.glob(self.output[0] + '.out_*') + \
            glob.glob(self.output[0] + '.err_*') + glob.glob(self.output[0] + '.done_*'):
            try:
                os.remove(filename)
            except Exception as e:
                env.logger.warning('Fail to remove {}: {}'.format(filename, e))

    def __call__(self, ifiles, pipeline=None):
        '''Execute action with input files ``ifiles`` with runtime information
        stored in ``pipeline``. This function is called by the pipeline and calls
        user-defined ``execute`` function.

        Parameters:

            ifiles (string or list of strings):
                input file names

            pipeline (an pipeline object):
                An Pipeline object for which the action is executed. The action
                can set or retrieve runtime information from a dictionary 
                ``pipeline.VARS``.

        Result:
            An action returns output files (parameter ``output`` of the action)
            if any output is given. Otherwise input files (``ifiles``) are passed
            through and returned.
        '''
        if self.output:
            if os.path.isfile(self.proc_info):
                with open(self.proc_info) as exe_info:
                    cmd = exe_info.readline().strip()
                if cmd == '; '.join(self.cmd).strip() and existAndNewerThan(self.output, ifiles,
                    md5file=self.proc_info):
                    env.logger.info('Reuse existing {}'.format(', '.join(self.output)))
                    self._bypass(ifiles, pipeline)
                    if self.output:
                        return self.output
                    else:
                        return ifiles
            # create directory if output directory does not exist
            for d in [os.path.split(os.path.abspath(x))[0] for x in self.output]:
                if not os.path.isdir(d):
                    try:
                        os.makedirs(d)
                    except Exception as e:
                        raise RuntimeError('Failed to create directory {} for output file: {}'.format(d, e))
        # We cannot ignore this step, but do we have all the input files?
        # If not, we will have to rewind the execution
        for ifile in ifiles:
            if not os.path.isfile(ifile):
                env.logger.warning('Rewind execution because input file {} does not exist.'.format(ifile))
                raise RewindExecution(ifile)
        #
        if self.output:
            with open(self.proc_info, 'w') as exe_info:
                exe_info.write('{}\n'.format('; '.join(self.cmd)))
                exe_info.write('#Start: {}\n'.format(time.asctime(time.localtime())))
                for f in ifiles:
                    if pipeline is not None and f == pipeline.VARS['null_input']:
                        continue
                    # for performance considerations, use partial MD5
                    exe_info.write('{}\t{}\t{}\n'.format(f, os.path.getsize(f),
                        calculateMD5(f, partial=True)))
        # now, run the job, write info if it is successfully finished.
        # Otherwise the job might be forked and it will record the signature by itself.
        ret = self._execute(ifiles, pipeline)
        if ret not in [True, False]:
            env.logger.warning('User defined execute function of a PipelineAction should return True or False')
        if ret:
            self._write_info()#
        #
        if self.output:
            return self.output
        else:
            return ifiles

# for backward compatibility
SkiptableAction=PipelineAction

try:
    from simulation import *
    hasSimuPOP = True
except ImportError as e:
    hasSimuPOP = False

class SequentialActions(PipelineAction):
    '''Define an action that calls a list of actions, specified by Action1,
    Action2 etc. This allows the specification of multiple small tasks in
    a single pipeline step.

    NOTE: this action is automatically applied if a list or tuple of actions
    are specified in the SPEC file (e.g. action=Action1(), Action2()).


    Examples:
        action=CheckCommands('bowtie'), CheckOutput('bowtie --version', '1.1.1')

    '''
    def __init__(self, actions):
        '''
        Parameters:
            actions (a tuple or list of actions):
                A list of actions that will be applied to 
        '''
        self.actions = []
        for a in actions:
            if hasattr(a, '__call__'):
                self.actions.append(a.__call__)
            else:
                self.actions.append(a)

    def __call__(self, ifiles, pipeline=None):
        '''
        Pass ifiles to the first action, take its output and pass it to
        the second action, and so on. Return the output from the last
        action as the result of this ``SequentialAction``.
        '''
        for a in self.actions:
            # the input of the next action is the output of the
            # previous action.
            ifiles = a(ifiles)
        # return the output of the last action
        return ifiles


class CheckVariantToolsVersion(PipelineAction):
    '''Check the version of variant tools and determine if it is
    recent enough to execute the pipeline.

    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Raises:
        Fail if the version of variant tools used to execute the
        pipeline is older than the specified version.

    Examples:
        action=CheckVariantToolsVersion('2.5.0')

    '''
    def __init__(self, version=''):
        '''
        Parameters:
            version (string):
                Oldest version of variant tools that can be used
                to execute this pipeline
        '''
        self.min_version = version
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        vtools_version = [int(x) for x in re.sub('\D', ' ', pipeline.VARS['vtools_version']).split()]
        # e.g. minimal 2.2.0, vtools 2.1.1
        if [int(x) for x in re.sub('\D', ' ', self.min_version).split()] > vtools_version:
            raise RuntimeError('Version {} is required to execute this pipeline. '
                'Please upgrade your installation of variant tools (version {})'
                .format(self.min_version, pipeline.VARS['vtools_version']))
        return ifiles


class ImportModules(PipelineAction):
    '''Import functions and action from a Python module. This action passed input
    files to output and does not change the pipeline.
    
    File Flow: Input passthrough, but import symbols to pipeline.

                     Pipeline
                        ^
            INPUT =============》INPUT

    Raises:
        Raise a RuntimeError if one or more modules can not
        be imported.

    Examples:
        action=ImportModules('DNASeq_tools.py')
        action=ImportModules(['DNASeq_tools.py', 'simuPOP.demography'])
    '''
    def __init__(self, modules=[]):
        '''Import one or more modules to be used by the existing pipeline. 
        
        Parameters:
            module (string or list of strings):
                One or more module, which can be either the name of a system module
                or a .py file. In the latter case, Variant Tools will try to locate
                the file directly (a full path can be given), look for the module in
                the path of the pipeline (if a local pipeline is used), or download
                from the Variant Tools Repository under directory pipeline.


        '''
        if isinstance(modules, str):
            self.modules = [modules]
        else:
            self.modules = modules
        # threads to _monitor the spawned jobs
        self._monitors = Queue.Queue()

    def __call__(self, ifiles, pipeline=None):
        for module in self.modules:
            # this is a path to a .py file
            if module.endswith('.py'):
                if os.path.isfile(module):
                    pyfile = module
                # if the .py file locates in the same directory as the pipeline file
                elif pipeline is not None \
                    and os.path.isfile(os.path.join(os.path.split(pipeline.spec_file)[0], module)):
                    pyfile = os.path.join(os.path.split(pipeline.spec_file)[0], module)
                else:
                    # try to download it from online
                    try:
                        pyfile = downloadFile('simulation/{}'.format(module))
                    except Exception as e:
                        try:
                            pyfile = downloadFile('pipeline/{}'.format(module))
                        except Exception as e:
                            raise ValueError('Failed to download required python module {}: {}'.format(module, e))
                try:
                    p,f = os.path.split(os.path.abspath(os.path.expanduser(pyfile)))
                    sys.path.append(p)
                    local_dict = __import__(f[:-3] if f.endswith('.py') else f, globals(), locals(), module.split('.', 1)[-1:])
                    env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
                    pipeline.GLOBALS.update(local_dict.__dict__)
                except Exception as e:
                    raise RuntimeError('Failed to import module {}: {}'.format(module, e))
            # now a system module
            else:
                try:
                    # allow loading from current directory
                    sys.path.append(os.getcwd())
                    local_dict = __import__(module, globals(), locals(), module.split('.', 1)[-1:])
                    env.logger.info('{} symbols are imported form module {}'.format(len(local_dict.__dict__), module))
                    pipeline.GLOBALS.update(local_dict.__dict__)
                except ImportError as e:
                    raise RuntimeError('Failed to import module {}: {}'.format(module, e))
        return ifiles


class CheckCommands(PipelineAction):
    '''Check the existence of specified commands and raise an error if one of
    the commands does not exist.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT 

    Raises:
        A RuntimeError will be raised if a command is not found.

    Examples:
        action=CheckCommands('java')
        action=CheckCommands(['java', 'tophat2'])
    '''
    def __init__(self, commands):
        '''
        Parameters:
            commands (string or list of strings):
                Name of one of more commands to be checked. No option is allowed.
        '''
        PipelineAction.__init__(self)
        if type(commands) == type(''):
            self.commands= [commands]
        else:
            self.commands= commands

    def __call__(self, ifiles, pipeline=None):
        for cmd in self.commands:
            if which(cmd) is None:
                raise RuntimeError('Command {} does not exist. Please install it and try again.'
                    .format(cmd))
            else:
                env.logger.info('Command {} is located.'.format(cmd))
        return ifiles


class CheckOutput(PipelineAction):
    '''Run a command and check if its output matches at least one of specified
    patterns. The pipeline will be terminated if failIfMismatch is set to True
    (default). Otherwise a warning message will be printed.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Raises:
        Raise a RuntimeError if the output of command does not match
        any of the patterns, if ``failIfMismatch`` is set to ``True``.

    Examples:
        action=CheckOutput('tophat2 --version', ['v2.0.13', 'v2.0.14'])

        # if strict_version is a command line parameter
        action=CheckOutput('samtools', '0.1.19', %(strict_version)s)
    '''
    def __init__(self, command, patterns, failIfMismatch=True):
        '''
        Parameters:
            command (string):
                A command (with or without options)

            patterns (string or list of strings):
                One or more patterns (usually a piece of version string)
                that will be compared to the output of ``command``

            failIfMismatch (boolean):
                If set to ``True`` (default), the action will terminate the
                pipeline if the output of command does not match any of the
                patterns. Otherwise a warning message will be printed when 
                the output of command does not match any of the patterns.
        '''
        self.command = command
        if isinstance(patterns, str):
            self.patterns = [patterns]
        else:
            self.patterns = patterns
        self.fail = failIfMismatch
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        try:
            # do not use subprocess.check_output because I need to get
            # output even when the command returns non-zero return code
            p = subprocess.Popen(self.command, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, shell=True)
            odata, edata = p.communicate()
            output = odata.decode() + edata.decode()
            env.logger.trace('Output of command "{}" is "{}"'
                .format(self.command, output))
        except Exception as e:
            raise RuntimeError('Failed to execute command "{}": {}'
                .format(self.cmd, e))
        #
        if all([re.search(x, output, re.MULTILINE) is None for x in self.patterns]):
            msg = ('Output of command "{}" ("{}") does not ' + 
                    'match specified regular expression {}.').format(self.command,
                        ' '.join(output[:40].split()), ' or '.join(self.patterns))
            if self.fail:
                raise RuntimeError(msg)
            else:
                env.logger.warning(msg)
        return ifiles

class CheckFiles(PipelineAction):
    '''Check the existence of specified files and raise an
    error if one of the files does not exist.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Raises:
        Raise a RuntimeError if any of the files is not found.

    Example:
        # assume gatk_path is a command line argument
        action=CheckFile('%(gatk_path)s/GenomeAnalysisTK.jar',
            'Please point --gatk_path to a directory with GenomeAnalysisTK.jar')
    '''
    def __init__(self, files, message=''):
        '''
        Parameters:
            files (string or list of strings):
                One or more files to check.

            message (string):
                A message when one of the files cannot be found.
        '''
        if type(files) == str:
            self.files = [files]
        else:
            self.files = files
        self.message = message
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        for f in self.files:
            if os.path.isfile(os.path.expanduser(f)):
                env.logger.info('{} is located.'.format(f))
            else:
                raise RuntimeError('Cannot locate {}: {}'.format(f, self.message))
        return ifiles


class CheckDirs(PipelineAction):
    '''Check the existence of specified directories and raise an
    error if one of the directories does not exist.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Raises:
        Raise a RuntimeError if any of the directories is not found.

    Example:
        action=CheckDirs('${CMD_OUTPUT}',
            'Value of parameter --output need to be an existing directory')
    '''
    def __init__(self, dirs, message=''):
        '''
        Parameters:
            files (string or list of strings):
                One or more directories to check.

            message (string):
                A message when one of the directories cannot be found.
        '''

        if type(dirs) == str:
            self.dirs = [dirs]
        else:
            self.dirs = dirs
        self.message = message
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        for d in self.dirs:
            if os.path.isdir(d):
                env.logger.info('Directory {} is located.'.format(d))
            else:
                raise RuntimeError('Cannot locate directory {}. {}'.format(d, self.message))
        return ifiles


class TerminateIf(PipelineAction):
    '''Terminate a pipeline if a condition is not met.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Raises:
        A RuntimeError will be raised to terminate the pipeline if
        the condition is met.

    Examples:
        action=TerminateIf(not '${CMD_OUTPUT}', 'No --output is specified.')
    '''
    def __init__(self, cond, message):
        '''
        Parameters:
            cond (boolean):
                True or False. In practice, ``cond`` is usually
                a lambda function that checks the existence of a file or value
                of a pipeline variable.

            message (string):
                A message to be outputted when the condition is met.
        '''
        self.cond = cond
        self.message = message
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        if self.cond:
            raise RuntimeError(self.message)
        return ifiles


class WarnIf(PipelineAction):
    '''Send a warning message if a condition is not met.

    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Examples:
        action=WarnIf('%(LGD)' == 'NA', 'Default value of parameter --LGD is used.')
    '''
    def __init__(self, cond, message):
        '''
        Parameters:
            cond (boolean):
                True or False. In practice, ``cond`` is usually
                a lambda function that checks the existence of a file or value
                of a pipeline variable.

            message (string):
                A message to be outputted when the condition is met.
        '''
        self.cond = cond
        self.message = message
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        if self.cond:
            env.logger.warning(self.message)
        return ifiles


class OutputText(PipelineAction):
    '''Write specified text to standard output, or a file if a filename is
    specified. The text can be a list of strings. A new line is added 
    automatically to each line of the text.

    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Examples:
        action=OutputText('Hey, the biggest part is done.')
    
    '''
    def __init__(self, text='', output=None, mode='a'):
        '''
        Parameters:
            text (string or list of strings):
                Text to be written to output.

            output (a file name or None):
                Output files. The text will be written to standard output
                if no output is specified.

            mode (string):
                Mode to open file. 'a' for append and 'w' for overwrite.

        '''
        if not isinstance(text, str):
            self.text = ''.join([str(x) + '\n' for x in text])
        else:
            self.text = text + '\n'
        self.filename = filename
        self.mode = mode
        PipelineAction.__init__(self, 'OutputText', filename if filename is not None else '')

    def __call__(self, ifiles, pipeline=None):
        if self.filename is not None:
            with open(self.filename, self.mode) as output:
                output.write(self.text)
        else:
            sys.stdout.write(self.text)        
        return ifiles


class FieldsFromTextFile(PipelineAction):
    '''Read a text file, guess its delimeter, field name (from header)
    and create field descriptions. If a vcf file is encountered, all
    fields will be exported.

    File Flow: extract format of input and output format.

        INPUT ==> Get Format ==> OUTPUT

    Raises:
        Raise a RuntimeError if this action failed to guess format (fields)
        from the input file.

    Examples:
        action=FieldsFromTextFile('format.txt')

    '''
    def __init__(self, output):
        '''
        Parameters:
            output:
                Output file that records the format of the input files.
        '''
        PipelineAction.__init__(self, 'FieldsFromTextFile', output)

    def _execute(self, ifiles, pipeline=None):
        if len(ifiles) > 1:
            env.logger.warning('Only the format of the first input file would be outputted.')
        try:
            if ifiles[0].endswith('.vcf') or ifiles[0].endswith('.vcf.gz'):
                showTrack(ifiles[0], self.output[0])
            else:
                with open(self.output[0], 'w') as fo:
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
        return True
       
class RewindExecution(Exception):
    pass

class NullAction(PipelineAction):
    '''A pipeline action that does nothing. This is usually used when the goal
    of the step is to change input, output, or assign variables to pipelines.
    The action will be assumed if an empty action line is given.

    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Example:
        action=
        action=NullAction()
    '''
    def __init__(self):
        '''A null action that does nothing.'''
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        return ifiles
        
class RunCommand(PipelineAction):
    """This action execute specified commands. If the pipeline is running
    in parallel mode and a submitter is specified, it will use the submitter
    command to execute the commands in a separate job. 
    
    File Flow:

        Input passthrough if no output file is specified.
            INPUT ====> INPUT
        Generate output if one or more output files are specified.
            INPUT ==> CMD ==> OUTPUT
        
    Raises:
        Raises an error if an command fails to execute.

    Examples:
        # simple commands without checking output
        action=RunCommand(cmd='vtools init myproj -f')
        
        action=RunCommand(cmd=[
            '[ -d ${DIR1} ] || mkdir -p ${DIR1}',
            '[ -d ${DIR2} ] || mkdir -p ${DIR2}'
            ])
        
        # multiple commands, change working directory
        # check output
        action=RunCommand([
        	'update_blastdb.pl human_genomic --decompress || true',
        	'update_blastdb.pl nt --decompress || true',
        	],
            working_dir='${NCBI_RESOURCE_DIR}/blast/db',
        	output=['${NCBI_RESOURCE_DIR}/blast/db/human_genomic.nal',
		        '${NCBI_RESOURCE_DIR}/blast/db/nt.nal']
            )
        
        # run command in background, with pipes
        action=RunCommand('''samtools view -h ${INPUT}
           | awk '$6 ~/N/' | awk '{ if ($9 ~ /^-/) {print $1"\t-"} else print $1"\t+" }'
           | sort -T ${TEMP_DIR} -u | wc -l > ${ALIGNMENT_OUT}/junction.count''',
           output='${ALIGNMENT_OUT}/junction.count',
           submitter='sh {} &')

    """
    def __init__(self, cmd='', output=[], working_dir=None, submitter=None, max_jobs=None):
        '''This action accepts one (a string) or more command (a list of strings)
        and executes them in a shell environment, possibly as a separate job. 
        
        Parameters:
            cmd (string or list of strings):
                One or more commands to execute. The commands will be executed in
                shell mode so pipes are allowed.

            output (string or list of strings):
                Expected output files of the action. If specified, the execution
                signature will be created to record the input, output and command
                of the action, and ignore the action if the signature matches of 
                a previous run.

            working_dir (None or string):
                Working directory of the command. Variant Tools will change to
                this directory before executing the commands if a valid directory
                is passed.

            submitter (None or string):
                If a submitter is specified and the pipeline is executed in multi-job
                mode (e.g. --jobs 2), a shell script will be written with the commands
                to be executed. The submitter command will be executed with ``{}`` in
                parameter ``submitter`` replaced by the name of shell script. For
                example, submitter='sh {} &' will run the job as a background job,
                and submitter='qsub -q long < {}' will submit the shell script to the
                long queue of a cluster system. Because the pipeline will be terminated
                if the submitter command fails, `qsub new_job ... && false` can be used
                to replace the running process by start a new job and terminate the
                existing process intentionally.

            max_jobs: (deprecated)
        '''
        self.submitter = submitter
        self.working_dir = working_dir
        if type(output) == str:
            self.output = [os.path.expanduser(output)]
        else:
            self.output = [os.path.expanduser(x) for x in output]
        #
        for filename in self.output:
            if filename.lower().rsplit('.', 1)[-1] in ['exe_info', 'lck']:
                raise RuntimeError('Output file with extension .exe_info and .lck are reserved.')
        #
        if not self.output:
            self.proc_out = None
            self.proc_err = None
            self.proc_lck = None
            self.proc_info = None
        else:
            self.proc_out = '{}.out_{}'.format(self.output[0], os.getpid())
            self.proc_err = '{}.err_{}'.format(self.output[0], os.getpid())
            self.proc_lck = '{}.lck'.format(self.output[0])
            self.proc_info = '{}.exe_info'.format(self.output[0])
        self.start_time = time.time()
        if not cmd:
            cmd = ['echo "None command executed."']
        PipelineAction.__init__(self, cmd=cmd, output=output)

    def _elapsed_time(self):
        '''Return the elapsed time in human readable format since start time'''
        second_elapsed = int(time.time() - self.start_time)
        days_elapsed = second_elapsed // 86400
        return ('{} days '.format(days_elapsed) if days_elapsed else '') + \
            time.strftime('%H:%M:%S', time.gmtime(second_elapsed % 86400))
      
    def _run_command(self):
        '''Call a list of external command cmd, raise an error if any of them
        fails. '''
        if self.proc_lck:
            env.lock(self.proc_lck, str(os.getpid()))
        for cur_cmd in self.cmd:
            env.logger.info('Running [[{}]]'.format(cur_cmd))
            ret = subprocess.call(cur_cmd, shell=True, 
                stdout=None if self.proc_out is None else open(self.proc_out, 'w'),
                stderr=None if self.proc_err is None else open(self.proc_err, 'w'),
                cwd=self.working_dir)
            if ret < 0:
                if self.output:
                    try:
                        env.unlock(self.proc_lck, str(os.getpid()))
                    except:
                        pass
                raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
                    .format(cur_cmd, -ret, self._elapsed_time()))
            elif ret > 0:
                if self.output:
                    with open(self.proc_err) as err:
                        for line in err.read().split('\n')[-50:]:
                            env.logger.error(line)
                    try:
                        env.unlock(self.proc_lck, str(os.getpid()))
                    except:
                        pass
                raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
                    .format(cur_cmd, self._elapsed_time(), ret))

    def _monitor(self):
        while True:
            if os.path.isfile(self.proc_done):
                break
            else:
                time.sleep(10)
        try:
            env.unlock(self.proc_lck, str(os.getpid()))
        except:
            env.logger.warning('Failed to remove lock for file {}'.format(self.output[0]))
            pass
        with open(self.proc_done) as done:
            ret = int(done.read().strip())
        #
        if ret < 0:
            raise RuntimeError("Command '{}' was terminated by signal {} after executing {}"
                .format('; '.join(self.cmd), -ret, self._elapsed_time()))
        elif ret > 0:
            if self.output:
                with open(self.proc_err) as err:
                    for line in err.read().split('\n')[-50:]:
                        env.logger.error(line)
            raise RuntimeError("Execution of command '{}' failed after {} (return code {})."
                .format('; '.join(self.cmd), self._elapsed_time(), ret))
        # remove the .done file
        if not self.output[0] in self.pipeline.THREADS:
            raise RuntimeError('Output is not waited by any threads')
        # DO NOT POP FROM ANOTHER THREAD, this will cause race condition
        # (unless we use thread safe dictionry). In this case, we only need
        # to monitor the status of threads from the master threads.
        #    self.pipeline.THREADS.pop(self.output[0])
        #
        # the thread will end here
        env.logger.trace('Thread for output {} ends.'.format(self.output[0]))
        for filename in glob.glob(self.output[0] + '.done_*'):
            try:
                os.remove(filename)
            except Exception as e:
                env.logger.warning('Fail to remove {}: {}'
                    .format(filename, e))

    def _submit_command(self):
        '''Submit a job and wait for its completion.'''
        # use full path because the command might be submitted to a remote machine
        self.proc_cmd = os.path.abspath(self.output[0] + '.sh')
        self.proc_done = os.path.abspath(self.output[0]) + '.done_{}'.format(os.getpid())
        #
        if os.path.isfile(self.proc_done):
            os.remove(self.proc_done)
        if self.proc_lck:
            env.lock(self.proc_lck, str(os.getpid()))
        #
        # create a batch file for execution
        with open(self.proc_cmd, 'w') as sh_file:
            sh_file.write('#PBS -o {}\n'.format(os.path.abspath(self.proc_out)))
            sh_file.write('#PBS -e {}\n'.format(os.path.abspath(self.proc_err)))
            sh_file.write('#PBS -N {}\n'.format(os.path.basename(self.output[0])))
            #sh_file.write('#PBS -N {}.{}_{}\n'.format(self.proc_err))
            sh_file.write('#PBS -V\n')
            # we try to reproduce the environment as much as possible becaus ehte
            # script might be executed in a different environment
            for k, v in os.environ.items():
                if any([k.startswith('x') for x in ('SSH', 'PBS', '_')]):
                    continue
                sh_file.write('export {}={}\n'.format(k, v))
            #
            sh_file.write('\ncd {}\n'.format(os.path.abspath(os.getcwd())))
            if self.working_dir is not None:
                sh_file.write('[ -d {0} ] || mkdir -p {0}\ncd {0}\n'.format(os.path.abspath(self.working_dir)))
            sh_file.write('\n'.join(self.cmd))
            #
            sh_file.write('\n\nCMD_RET=$?\nif [ $CMD_RET == 0 ]; then vtools admin --record_exe_info {} {}; fi\n'
                .format(os.getpid(), ' '.join(self.output)))
            # a signal to show the successful completion of the job
            sh_file.write('\necho $CMD_RET > {}\n'.format(self.proc_done))
        #
        # try to submit command
        if '{}' in self.submitter:
            submit_cmd = self.submitter.replace('{}', self.proc_cmd)
        else:
            submit_cmd = self.submitter
        #
        env.logger.info('Running job {} with command "{}" from directory {}'.format(
            self.proc_cmd, submit_cmd, os.getcwd()))
        ret = subprocess.call(submit_cmd, shell=True,
            stdout=open(self.proc_out, 'w'), stderr=open(self.proc_err, 'w'),
            cwd=self.working_dir)
        if ret != 0:
            try:
                env.unlock(self.proc_out, str(os.getpid()))
            except:
                pass
        # 
        if ret < 0:
            raise RuntimeError("Failed to submit job {} due to signal {} (submitter='{}')" .format(self.proc_cmd, -ret, self.submitter))
        elif ret > 0:
            if os.path.isfile(self.proc_err):
                 with open(self.proc_err) as err:
                     msg = err.read()
            else:
                msg = ''
            raise RuntimeError("Failed to submit job {} using submiter '{}': {}".format(self.proc_cmd, self.submitter, msg))
        else:
            t = threading.Thread(target=self._monitor)
            t.daemon = True
            t.start()
            if self.output[0] in self.pipeline.THREADS:
                raise RuntimeError('Two spawned jobs have the same self.output[0] file {}'.format(self.output[0]))
            self.pipeline.THREADS[self.output[0]] = t


    def _execute(self, ifiles, pipeline=None):
        # substitute cmd by input_files and output_files
        if pipeline.jobs > 1 and self.submitter is not None and not self.output:
            env.logger.warning('Fail to execute in parallel because no output is specified.')
        #
        self.pipeline = pipeline
        # Submit the job on a cluster system
        # 1. if there is output, otherwise we cannot track the status of the job
        # 2. if a submit command is specified
        # 3. if --jobs with a value greater than 1 is used.
        if self.output and pipeline.jobs > 1 and self.submitter is not None:
            self._submit_command()
            return False
        else:
            self._run_command()
            return True

class DecompressFiles(PipelineAction):
    '''This action gets a list of input files from input file, decompressing
    input files (.tar.gz, .zip, etc) if necessary. The decompressed files
    are returned as output. One particular feature of this action is that
    it records content of large tar or tar.gz files to a manifest file and
    ignores the step if the manifest file exists.

    File Flow: Decompress input files

        INPUT ==> Decompress ==> OUTPUT

    Examples:
        action=DecompressFiles()
    '''
    def __init__(self, dest_dir=None):
        '''
        Parameters:
            dest_dir (None or string):
                Destination directory, default to current directory.
        '''
        self.dest_dir = dest_dir if dest_dir else '.'
        PipelineAction.__init__(self)

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


class RemoveIntermediateFiles(PipelineAction):
    '''This action removes specified files (not the step input files) and replaces
    them with their signature (file size, md5 signature etc). A pipeline can bypass
    completed steps with these files as input or output by checking the signatures.
    In contrast, the steps would have to be re-run if the files are removed from the
    file system. 
    
    File Flow: Input passthrough. Specified files are replaced by their signature.

        INPUT ====> INPUT

    Examples:
        action=RemoveIntermediateFiles('${OUTPUT200}')
        action=RemoveIntermediateFiles('${OUTPUT200} ${OUTPUT330}')
        action=RemoveIntermediateFiles(['${OUTPUT200}', '${OUTPUT330}'])
    '''
    def __init__(self, files):
        '''Replace ``files`` with their signatures. This pipeline passes its 
        input to output and does not change the flow of pipeline.
 
        Parameters:
            files (string or list of strings)
                One or more files to be removed. Multiple files can be specified
                in the same string if they are separated by spaces.

        '''
        if isinstance(files, str):
            self.files_to_remove = [files]
        else:
            self.files_to_remove = files
        PipelineAction.__init__(self)

    def _getFiles(self):
        for name in self.files_to_remove:
            files = shlex.split(name)
            for f in files:
                yield f

    def __call__(self, ifiles, pipeline=None):
        env.logger.trace('Remove intermediate files {}'.format(' '.join(self.files_to_remove)))
        for f in self._getFiles():
            if not os.path.isfile(f):
                if os.path.isfile(f + '.file_info'):
                    env.logger.info('Keeping existing {}.file_info.'.format(f))
                else:
                    raise RuntimeError('Failed to create {}.file_info: Missing input file.'
                        .format(f))
            else:
                FileInfo(f).save()
                env.logger.info('Replace {0} with {0}.file_info'.format(f))
                try:
                    os.remove(f)
                except e:
                    env.logger.warning('Failed to remove intermediate file {}'.format(f))
        return ifiles


class LinkToDir(PipelineAction):
    '''Create hard links of input files to a specified directory. This is 
    usually used to link input files to a common cache directory so that 
    all operations can be performed on that directory.

    File Flow: Link input files to specified destination directory.

        INPUT == LINK ==> DEST_DIR/INPUT

    Examples:
        action=LinkToDir('cache')

    '''
    def __init__(self, dest_dir):
        '''
        Parameters:
            dest_dir (string):
                A directory to which input files will be linked to.
                The directory will be created if it does not exist.
        '''
        self.dest = dest_dir
        if not os.path.isdir(self.dest):
            env.logger.info('Creating directory {}'.format(self.dest))
            try:
                os.makedirs(self.dest)
            except Exception as e:
                raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))
            if not os.path.isdir(self.dest):
                raise RuntimeError('Failed to create directory {}: {}'.format(self.dest, e))
        PipelineAction.__init__(self)

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
                            env.logger.trace('Reusing existing linked file_info file: {}'
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
                        env.logger.trace('Reusing existing linked file: {}'.format(dest_file))
                else:
                    env.logger.info('Linking {} to {}'.format(filename, self.dest))
                    os.link(filename, dest_file)
            ofiles.append(os.path.join(self.dest, basename))
        return ofiles


class DownloadResource(PipelineAction):
    '''Download resources to specified destination directory. dest_dir can
    be a full path name or a directory relative to 
    $local_resource/pipeline_resource where $local_resource is the local
    resource directory of the project (default to ~/.variant_tools,
    see runtime option local_resource for details). The default pipeline 
    resource directory is $local_resource/pipeline_resource/NAME where NAME
    is the name of the pipeline.
    
    File Flow: Input passthrough. 

        INPUT ====> INPUT

    Examples:
        action=DownloadResource(resource='ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Homo_sapiens/UCSC/hg19/Homo_sapiens_UCSC_hg19.tar.gz',
             dest_dir="${LOCAL_RESOURCE}/iGenomes")

        action=DownloadResource(resource='ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz
            ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.sites.vcf.gz.md5',
            dest_dir='${LOCAL_RESOURCE/GATK')
    
    NOTE:
        1. If FILE.md5 file is downloaded, it will be used to validate FILE.
        2. The resources will be automatically decompressed. You would get both 
            FILE and FILE.gz if you downloaded FILE.gz
    '''
    def __init__(self, resource, dest_dir):
        '''Download resources from specified URLs in ``resource``. 

        Parameters:
            dest_dir:
                Directory where the downloaded resources will be placed.
        '''
        self.resource = [x for x in resource.split() if x]
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
        PipelineAction.__init__(self)

    def __call__(self, ifiles, pipeline=None):
        saved_dir = os.getcwd()
        os.chdir(self.pipeline_resource)
        ofiles, md5files = self._downloadFiles(ifiles)
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
 

class _CaseInsensitiveDict(MutableMapping):
    """A case-insensitive ``dict``-like object.
    That limits the type of items to string or list of strings.
    """
    def __init__(self, data=None, **kwargs):
        self._store = dict()
        if data is None:
            data = {}
        self.update(data, **kwargs)

    def __setitem__(self, key, value):
        # Use the uppercased key for lookups, but store the actual
        # key alongside the value.
        if key.upper in self._store and value == self._store[key.upper()][1]:
            env.logger.warning('Changing value of pipeline variable ({} from {} to {}) is strongly discouraged.'
                .format(key, self._store[key.upper()][1], value))
        if not isinstance(value, (str, list, tuple)):
            raise ValueError('Only string or list of strings are allowed for pipeline variables: {} for key {}'.format(value, key))
        if isinstance(value, (list, tuple)) and not all([isinstance(x, str) for x in value]):
            raise ValueError('Only string or list of strings are allowed for pipeline variables: {} for key {}'.format(value, key))
        self._store[key.upper()] = (key, value)

    def __getitem__(self, key):
        return self._store[key.upper()][1]

    def __delitem__(self, key):
        del self._store[key.upper()]

    def __iter__(self):
        return (casedkey for casedkey, mappedvalue in self._store.values())

    def __len__(self):
        return len(self._store)

    def upper_items(self):
        """Like iteritems(), but with all uppercase keys."""
        return (
            (upperkey, keyval[1])
            for (upperkey, keyval)
            in self._store.items()
        )

    def __eq__(self, other):
        if isinstance(other, collections.Mapping):
            other = _CaseInsensitiveDict(other)
        else:
            return NotImplemented
        # Compare insensitively
        return dict(self.upper_items()) == dict(other.upper_items())

    # Copy is required
    def copy(self):
         return _CaseInsensitiveDict(self._store.values())

    def __repr__(self):
        return '%s(%r)' % (self.__class__.__name__, dict(self.items()))


class Pipeline:
    '''The Variant Tools pipeline class. Its instance will be passed to each action
    to provide runtime information. An action should not change any attribute of
    the pipeline, except for setting additional variables through its ``VARS``
    dictionary. Note that VARS is a case-insensitive dictionary but it is generally
    recommended to use CAPTICAL names for pipeline variables. '''
    def __init__(self, name, extra_args=[], pipeline_type='pipeline', verbosity=None, jobs=1):
        self.pipeline = PipelineDescription(name, extra_args, pipeline_type)
        self.spec_file = self.pipeline.spec_file
        self.verbosity = verbosity
        self.jobs = jobs

    def execute(self, pname, input_files=[], output_files=[], **kwargs):
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
        # the project will be opened when needed
        with Project(mode=['ALLOW_NO_PROJ', 'READ_ONLY'], verbosity=self.verbosity) as proj:
            self.VARS = _CaseInsensitiveDict(
                cmd_input=input_files,
                cmd_output=output_files,
                temp_dir=env.temp_dir,
                cache_dir=env.cache_dir,
                local_resource=env.local_resource,
                ref_genome_build=proj.build if proj.build is not None else '',
                pipeline_name=pname,
                spec_file=self.spec_file,
                model_name=pname,
                null_input=env.null_input,
                vtools_version=proj.version)
        self.VARS.update({k:str(v) for k,v in kwargs.items()})
        # we need to put self.pipeline.pipeline_vars in self.VARS because
        # they might refer to each other
        self.VARS.update(self.pipeline.pipeline_vars)
        for key, val in self.pipeline.pipeline_vars.items():
            self.VARS[key.lower()] = substituteVars(val, self.VARS)
        #
        self.GLOBALS = {}
        self.THREADS = {}
        #
        ifiles = input_files
        step_index = 0
        rewind_count = 0
        while True:
            # step_index can jump back and forth depending on the 
            # execution status of each step
            command = psteps[step_index]
            self.VARS['pipeline_step'] = command.index
            env.logger.info('Executing [[{}.{}_{}]]: {}'
                .format(self.pipeline.name, pname, command.index, 
                    ' '.join(command.comment.split())))
            # substitute ${} variables
            if command.input is None:
                step_input = ifiles
            # if this is an empty string
            elif not command.input.strip():
                step_input = []
            else:
                step_input = shlex.split(substituteVars(command.input, self.VARS))
            #
            # if there is no input file?
            if not step_input:
                env.logger.debug('Step ignored because of no input file for step {}_{}.'.format(pname, command.index))
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
                    emitter = eval(substituteVars(command.input_emitter, self.VARS), globals(), self.GLOBALS)
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
                    action = substituteVars(command.action, self.VARS)
                    env.logger.trace('Emitted input of step {}_{}: {}'
                        .format(pname, command.index, ig))
                    env.logger.trace('Action of step {}_{}: {}'
                        .format(pname, command.index, action))
                    # check if the input file is ready. This is used for
                    # parallel execution of the pipeline while the input file
                    # might be worked on by another job
                    for ifile in ig:
                        # is ifile in any of the output files?
                        if ifile in self.THREADS:
                            # wait for the thread to complete
                            env.logger.info('Waiting for the input file {} to be available.'
                                .format(ifile))
                            while self.THREADS[ifile].isAlive():
                                self.THREADS[ifile].join(5)
                            # thread closed, remove from self.THREADS
                            self.THREADS.pop(ifile)
                        if not (os.path.isfile(ifile) or os.path.isfile(ifile + '.file_info')):
                            raise RewindExecution(ifile)
                    #
                    if not action.strip():
                        action = 'NullAction()'
                    action = eval(action, globals(), self.GLOBALS)
                    if isinstance(action, (tuple, list)):
                        action = SequentialActions(action)
                    if not issubclass(action.__class__, PipelineAction):
                        env.logger.warning('Pipeline action {} is not a subclass of PipelineAction'.format(action.__class__))
                    # pass the Pipeline object itself to action
                    # this allows the action to have access to pipeline variables
                    # and other options
                    ofiles = action(ig, self)
                    if type(ofiles) == str:
                        step_output.append(ofiles)
                    else:
                        step_output.extend(ofiles)
                # wait for all pending jobs to finish
                self.VARS['output{}'.format(command.index)] = step_output
                env.logger.debug('OUTPUT of step {}_{}: {}'
                    .format(pname, command.index, step_output))
                for f in step_output:
                    if not (os.path.isfile(f) or os.path.isfile(f + '.file_info') or f in self.THREADS):
                        raise RuntimeError('Output file {} does not exist after '
                            'completion of step {}_{}'
                            .format(f, pname, command.index))
                for key, val in command.pipeline_vars:
                    self.VARS[key.lower()] = substituteVars(val, self.VARS)
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
                # unfortunately, an input file has been removed (replaced by .file_info) but
                # a later steps need it. We will have to figure out how to create this 
                # file by looking backward ...
                to_be_regenerated = [x for x in step_input if not os.path.isfile(x)]
                # we need to check if this file is actually generated at all before
                # otherwise a misspecified input file would cause the whole pipline 
                # to start from step 1 again and again
                all_input_and_output_files = []
                for k,v in self.VARS.items():
                    if ((k.startswith('INPUT') and k[5:].isdigit() and int(k[5:]) < int(command.index)) or \
                        (k.startswith('OUTPUT') and k[6:].isdigit() and int(k[6:]) < int(command.index))) and \
                         isinstance(v, (tuple, list)):
                        all_input_and_output_files.extend(v)
                #
                for x in to_be_regenerated:
                    if x not in all_input_and_output_files:
                        raise RuntimeError('Specified input file "{}" does not exist and is not '
                            'generated by any previous step.'.format(x))
                    # remove all fony files so that they will be re-generated
                    if os.path.isfile(x + '.file_info'):
                        os.remove(x + '.file_info')
                remaining = [x for x in to_be_regenerated]
                env.logger.debug('Missing input file {}'.format(', '.join(remaining)))
                while step_index > 0:
                    step_index -= 1
                    command = psteps[step_index]
                    # remove all fony files so that they will be re-generated
                    for x in self.VARS['input{}'.format(command.index)]:
                        if os.path.isfile(x + '.file_info'):
                            env.logger.debug('Remove file info {}'.format(x + '.file_info'))
                            os.remove(x + '.file_info')
                    # if any of the input file does not exist, go back further
                    if not all([os.path.isfile(x) for x in self.VARS['input{}'.format(command.index)]]):
                        env.logger.debug('Not all input files are available: {}'.format(', '.join(self.VARS['input{}'.format(command.index)])))
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
        #
        # at the end of pipeline wait for all threads to complete
        if self.THREADS:
            for k, v in self.THREADS.items():
                env.logger.trace('Waiting for {} to be completed.'.format(k))
                while v.isAlive():
                    v.join(5)
                # thread closed, remove from self.THREADS
                self.THREADS.pop(ifile)


def executeArguments(parser):
    parser.add_argument('specfile', metavar='SPECFILE',
        help='''Name of a pipeline configuration file, which can be a
            path to a .pipeline file (with or without extension) or one
            of the online pipelines listed by command "vtools show pipelines".
            For backward compatibility, if no input and output files are
            specified (options --input and --output), values of this option 
            is treated as a SQL query that will be executed against the project
            database, with project genotype database attached as "genotype" and
            annotation databases attached by their names.''')
    parser.add_argument('pipelines', nargs='*', metavar='PIPELINES',
        help='''Name of one or more pipelines defined in SPECFILE, which can be
            ignored if the SPECFILE only defines one pipeline. Please use 
            command "vtools show pipeline SPECFILE" for details of available
            pipelines in SPECFILE, including pipeline-specific parameters that
            could be used to change the default behavior of the pipelines.''')
    parser.add_argument('-i', '--input', nargs='*', metavar='INPUT', default=[],
        help='''Input to the pipelines, usually a list of input files, that
            will bepassed to the pipelines as variable ${CMD_INPUT}.''')
    parser.add_argument('-o', '--output', nargs='*', metavar='OUTPUT', default=[],
        help='''Output of the pipelines, usually a list of output files, that
            will be passed to the pipelines as variable ${CMD_OUTPUT}.''')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Execute the pipeline in parallel model if a number other than
            1 is specified. In this mode, the RunCommand action will create
            a shell script and submit the job using a command specified by
            option ``submitter``,  if this parameter is defined.''')
    parser.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter used to output results of a SQL query.''')

def execute(args):
    # to keep backward compatibility, the vtools execute command
    # can execute a SQL query and a pipeline
    def executeQuery():
        with Project(verbosity=args.verbosity) as proj:
            # if there is no output, 
            proj.db.attach('{}_genotype'.format(proj.name), 'genotype')
            # for backward compatibility
            proj.db.attach('{}_genotype'.format(proj.name))
            cur = proj.db.cursor()
            # 
            query = args.specfile + ' ' + ' '.join(args.pipelines)
            if query.upper().startswith('SELECT'):
                env.logger.trace('Analyze statement: "{}"'.format(query))
                cur.execute('EXPLAIN QUERY PLAN ' + query)
                for rec in cur:
                    env.logger.trace('\t'.join([str(x) for x in rec]))
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
    def executePipeline():                
        pipeline = Pipeline(args.specfile, extra_args=args.unknown_args, verbosity=args.verbosity,
            jobs=args.jobs)
        # unspecified
        if not args.pipelines:
            pipeline.execute(None, args.input, args.output)
        else:
            for name in args.pipelines:
                pipeline.execute(name, args.input, args.output)
    # 
    try:
        env.verbosity = args.verbosity
        env.logger = None
        # definitely a pipeline
        if args.specfile.endswith('.pipeline') or args.input or args.output or args.unknown_args:
            executePipeline()
        # definitely a sql query
        elif args.delimiter != '\t':
            executeQuery()
        else:
            try:
                # try to execute as a SQL query
                executeQuery()
            except RuntimeError as e:
                env.logger.debug('Failed to execute {} as SQL query: {}'
                    .format(' '.join(args.pipelines), e))
                executePipeline()
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
    parser.add_argument('specfile', metavar='SPECFILE',
        help='''Name of a model specification file, which can be the name of an
            online specification file, or path to a local .pipeline file. Please
            use command "vtools show simulations" to get a list all available
            simulation models.''')
    parser.add_argument('models', nargs='*', metavar='MODELS',
        help='''Name of one or more simulation models defined in SPECFILE, which
            can be ignored if the SPECFILE only defines one simulation model.
            Please use command "vtools show simulation SPECFILE" for details
            of available models in SPECFILE, including model-specific parameters
            that could be used to change the default behavior of these models.''')
    parser.add_argument('--seed', type=int,
        help='''Random seed for the simulation. A random seed will be used by
            default but a specific seed could be used to reproduce a previously
            executed simulation.''')
    parser.add_argument('--replicates', default=1, type=int,
        help='''Number of consecutive replications to simulate''')
    parser.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs to execute, for steps
            of a pipeline that allows multi-processing.''')

def simulate_replicate(args, rep):
    try:
        env.verbosity = args.verbosity
        env.logger = None
        # step 1, create a simulation configuration file.
        model_name = os.path.basename(args.specfile).split('.', 1)[0]
        if args.seed is None:
            args.seed = random.randint(1, 2**32-1)
        if not os.path.isdir(env.cache_dir):
            os.mkdir(env.cache_dir)

        # set random seed of simulators
        random.seed(args.seed + rep)
        if not args.models:
            cfg_file = '{}/{}_{}.cfg'.format(env.cache_dir, model_name, args.seed + rep)
        else:
            cfg_file = '{}/{}_{}_{}.cfg'.format(env.cache_dir, model_name, '_'.join(args.models), args.seed + rep)
        #
        with open(cfg_file, 'w') as cfg:
            cfg.write('model={} {}\n'.format(args.specfile, ' '.join(args.models)))
            cfg.write('seed={}\n'.format(args.seed + rep))
            if '--seed' in sys.argv:
                # skip the seed option so to stop pipeline from distinguishing the two commands
                cmd_args = sys.argv[:sys.argv.index('--seed')] + sys.argv[sys.argv.index('--seed') + 2:]
                cfg.write("command=vtools {}\n".format(subprocess.list2cmdline(cmd_args[1:])))
            else:
                cfg.write("command={}\n".format(env.command_line))
        #
        env.logger.info('Starting simulation [[{}]]'.format(cfg_file))
        pipeline = Pipeline(args.specfile, extra_args=args.unknown_args,
            pipeline_type='simulation', verbosity=args.verbosity, jobs=args.jobs)
        # using a pool of simulators 
        if not args.models:
            pipeline.execute(None, [cfg_file], [], seed=args.seed+rep)
        else:
            for name in args.models:
                pipeline.execute(name, [cfg_file], [], seed=args.seed+rep)
    except Exception as e:
        env.logger.error('Failed to simulate replicate {} of model {}: {}'.format(rep, model_name, e))
        sys.exit(1)

def simulate(args):
    if not hasSimuPOP:
        raise RuntimeError('Please install simuPOP before running any simulation using Variant Simulation Tools.')
    #
    try:
        ret = 0
        if args.replicates <= 0:
            raise ValueError('No replication studies is requested.')
        #
        for rep in range(args.replicates):
            p = Process(target=simulate_replicate, args=(args, rep))
            p.start()
            p.join()
            if ret == 0:
                ret = p.exitcode
        # return fail if any of the replicates fails
        sys.exit(ret)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)

if __name__ == '__main__':
    pass
