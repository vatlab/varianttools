#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2004 - 2011 Bo Peng (bpeng@mdanderson.org)
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
'''This module provides utility functions and a database engine for sqlite'''

import os
import sys
import glob
import logging
import subprocess
import urllib.request, urllib.parse, urllib.error
import time
import tempfile
import tokenize
import gzip
import traceback
import threading
import re
import shlex
import stat
import random
import shutil
import hashlib
import configparser
from html.parser import HTMLParser
import tarfile
import binascii
from itertools import chain
from shutil import which

# import site_options
from . import site_options

default_user_options = [
    ('sqlite_pragma', '', '''
; separated pragmas for sqlite database that can be used to optimize the
performance of sqlite database operations. Please check the sqlite manual
for more details.'''),
    ('import_num_of_readers', 2, '''
Number of processes to read from input files during multi-process import.
A smaller number might provide better performance for file systems with
slow random I/O access.'''),
    ('associate_num_of_readers', None, '''
Number of processes to read genotype data for association tests. The default
value is the minimum of value of option --jobs and 8. A smaller number
might provide better performance for file system with slow random I/O
access.'''),
    ('local_resource', '~/.variant_tools', '''
A directory to store variant tools related resources such as reference genomes
and annotation database. This option will be ignored if a writable shared
resource directory is specified by the system admin.'''),
    ('temp_dir', None, '''
root of tempory directory to store temporary files (default to system default)
Setting it to a different physical disk than user projects can generally
improve the performance of variant tools. It should also be set to a directory
with at least 500G of free diskspace if the default temporary partition
is small.'''),
    ('search_path',
     'https://bioinformatics.mdanderson.org/Software/VariantTools/repository/;https://bioinformatics.mdanderson.org/Software/VariantTools/archive/',
     '''
A ;-separated list of URL that host the variant tools repository. This option
should only be changed if you have created a local mirror of the variant tools
repository. Adding the URL before the default URL might provide better
downloading performance for your users. Removing the default URL is possible
but not recommended.'''),
    ('shared_resource', None, '''
 A directory for shared resource files. It can be configured as
 1. NO shared resource (default). All users maintain their own resource
   directory ($local_resource, which is usually ~/.variant_tools).

 2. A read-only directory with a mirror of the variant tools repository, with
   .DB.gz files decompressed in the annoDB directory. This is important because
   otherwise each user will have to decompress the files in their local resource
   directory. The system admin can choose to remove outdated databases to reduce
   the use of disk space. This option requires regular update of the resources.

 3. A directory that is writable by all users. The resources will be downloaded
   to this directory by users, and shared by all users. This option is easier
   to implement and requires less maintenance. The system admin can choose to
   mirror the variant tools repository and let the users to keep it up to date.'''
    ),
    ('user_stash', '~/.variant_tools', '''
 A ;-separated list of directories that stores personally-generated file formats,
 annotation databases, pipelines and such. These directories are searched if
 the requested file is not available in any online repository (or their local
 copy). Files under these directories are NOT maintained by variant tools
 (no manifest or md5 signatures are monitors). Currently only local files are
 allowed (no URL to a remote server). By setting this option to ~/.variant_tools
 (default) resource files under that directory will be usable even if they
 are not managed by variant tools''')
]

if os.path.isfile(os.path.expanduser('~/.variant_tools/user_options.py')):
    sys.path.insert(0, os.path.expanduser('~/.variant_tools'))
    try:
        _user_options = __import__('user_options', globals(), locals()).__dict__
    except Exception:
        print(
            'Failed to load user settings from ~/.variant_tools/user_options.py'
        )
    _user_options.update(
        {x: y for x, y, z in default_user_options if x not in _user_options})
else:
    if not os.path.isdir(os.path.expanduser('~/.variant_tools')):
        os.mkdir(os.path.expanduser('~/.variant_tools'))
    with open(os.path.expanduser('~/.variant_tools/user_options.py'),
              'w') as uo:
        uo.write('''#!/usr/bin/env python
# User configuration file that overrides settings from site-wise installation of variant tools.
# This file should be placed in $HOME/.variant_tools/user_options.py. Please
# note that this is a Python file so the strings should be quoted.
#
''' + '\n'.join([
            '#{}\n{}={}\n'.format(
                z.replace('\n', '\n# '), x,
                "'{}'".format(y) if isinstance(y, str) else str(y))
            for x, y, z in default_user_options
        ]))
    _user_options = {x: y for x, y, z in default_user_options}
#
# overriding site option with local setting
for k, v in list(_user_options.items()):
    if not hasattr(site_options, k):
        setattr(site_options, k, v)

try:
    # not all platforms/installations of python support bz2
    import bz2
    bz2_support = True
except:
    bz2_support = False

try:
    import pickle
    from io import StringIO
    import urllib.parse as urlparse
    import variant_tools.vt_sqlite3 as sqlite3
    from variant_tools.cgatools import CrrFile, Location, Range
    from variant_tools.vt_sqlite3 import OperationalError
except ImportError as e:
    sys.exit(
        'Failed to import module ({})\n'
        'Please verify if you have installed variant tools successfully (using command '
        '"python setup.py install")'.format(e))


class ColoredFormatter(logging.Formatter):
    ''' A logging formatter that uses color to differntiate logging messages
    and emphasize texts. Texts that would be empahsized are quoted with
    double backslashes (`` ``).
    '''

    def __init__(self, msg):
        logging.Formatter.__init__(self, msg)
        #
        # color for different logging levels. The current terminal color
        # is used for INFO
        self.LEVEL_COLOR = {
            'TRACE': 'DARK_CYAN',
            'DEBUG': 'BLUE',
            'WARNING': 'PURPLE',
            'ERROR': 'RED',
            'CRITICAL': 'RED_BG',
        }
        self.COLOR_CODE = {
            'ENDC': 0,  # RESET COLOR
            'BOLD': 1,
            'UNDERLINE': 4,
            'BLINK': 5,
            'INVERT': 7,
            'CONCEALD': 8,
            'STRIKE': 9,
            'GREY30': 90,
            'GREY40': 2,
            'GREY65': 37,
            'GREY70': 97,
            'GREY20_BG': 40,
            'GREY33_BG': 100,
            'GREY80_BG': 47,
            'GREY93_BG': 107,
            'DARK_RED': 31,
            'RED': 91,
            'RED_BG': 41,
            'LIGHT_RED_BG': 101,
            'DARK_YELLOW': 33,
            'YELLOW': 93,
            'YELLOW_BG': 43,
            'LIGHT_YELLOW_BG': 103,
            'DARK_BLUE': 34,
            'BLUE': 94,
            'BLUE_BG': 44,
            'LIGHT_BLUE_BG': 104,
            'DARK_MAGENTA': 35,
            'PURPLE': 95,
            'MAGENTA_BG': 45,
            'LIGHT_PURPLE_BG': 105,
            'DARK_CYAN': 36,
            'AUQA': 96,
            'CYAN_BG': 46,
            'LIGHT_AUQA_BG': 106,
            'DARK_GREEN': 32,
            'GREEN': 92,
            'GREEN_BG': 42,
            'LIGHT_GREEN_BG': 102,
            'BLACK': 30,
        }

    def colorstr(self, astr, color):
        return '\033[{}m{}\033[{}m'.format(color, astr, self.COLOR_CODE['ENDC'])

    def emphasize(self, msg, level_color=0):
        # display text within `` and `` in green
        return re.sub(r'``([^`]*)``', '\033[32m\\1\033[{}m'.format(level_color),
                      str(msg))

    def format(self, record):
        level_name = record.levelname
        if level_name in self.LEVEL_COLOR:
            level_color = self.COLOR_CODE[self.LEVEL_COLOR[level_name]]
            record.color_levelname = self.colorstr(level_name, level_color)
            record.color_name = self.colorstr(record.name,
                                              self.COLOR_CODE['BOLD'])
            record.color_msg = self.colorstr(
                self.emphasize(record.msg, level_color), level_color)
        else:
            # for INFO, use default color
            record.color_levelname = record.levelname
            record.color_msg = self.emphasize(record.msg)
        return logging.Formatter.format(self, record)


class RuntimeEnvironments(object):
    # the following make RuntimeEnvironments a singleton class
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            # *args, **kwargs are not passed to avoid
            # DeprecationWarning: object.__new__() takes no parameters
            # cls._instance = super(Singleton, cls).__new__(cls, *args, **kwargs)
            cls._instance = super(RuntimeEnvironments,
                                  cls).__new__(cls)  #, *args, **kwargs)
        return cls._instance

    def __init__(self):
        # these options could be set persistently
        self.hidden_options = {
            'term_width':
                (None,
                 'If set to a fix number, assuming a fixed terminal width for'
                 'outputs. This is used for testing purpose only.')
        }
        self.persistent_options = {
            'logfile_verbosity': (
                '2', 'Verbosity level of the log file, can be 0 for warning '
                'and error only, 1 for general information, or 2 for general and debug information.'
            ),
            'verbosity': (
                '1',
                'Default verbosity level (to the standard output) of the project. '
                'This option can be set during vtools init where the verbosity level set by option'
                ' --verbosity will be set as the project default.'),
            # 'check_update': (True, 'Automatically check update of variant tools and resources.'),
            'sqlite_pragma':
                (site_options.sqlite_pragma,
                 'pragmas for sqlite database that can be used to optimize the '
                 'performance of database operations.'),
            'import_num_of_readers': (
                site_options.import_num_of_readers,
                'variant tools by default uses two processes to read from '
                'input files during multi-process importing (--jobs > 0). You can want to set it '
                'to zero if a single process provides better performance or reduces disk traffic.'
            ),
            # a temporary directory that is used to store temporary files. Will be
            # cleared after project is closed.
            'temp_dir': (
                site_options.temp_dir,
                'Use the specified temporary directory to store temporary files '
                'to improve performance (use separate disks for project and temp files), or '
                'avoid problems due to insufficient disk space.'),
            'treat_missing_as_wildtype': (
                'False', 'Treat missing values as wildtype alleles for '
                'association tests. This option is used when samples are called individuals or '
                'in batch so genotypes for some samples are ignored and treated as missing if '
                'they consist of all wildtype alleles. This option should be used with caution '
                'because it convert real missing genotypes and genotypes removed due to, for '
                'example low quality score, to wildtype genotypes.'),
            'association_timeout': (
                None, 'Cancel associate test and return special values '
                'when a test lasts more than specified time (in seconds). The default '
                'value of this option is None, which stands for no time limit.'
            ),
            'associate_num_of_readers': (
                site_options.associate_num_of_readers,
                'Use specified number of processes to read '
                'genotype data for association tests. The default value is the minimum of value '
                'of option --jobs and 8. Note that a large number of reading processes might '
                'lead to degraded performance or errors due to disk access limits.'
            ),
            'search_path': (
                site_options.search_path, 'A ;-separated list of '
                'directories and URLs that are used to locate annotation database (.ann, .DB), '
                'file format (.fmt) and other files. Reset this option allows alternative '
                'local or online storage of such files. variant tools will append trailing '
                'directories such as annoDB for certain types of data so only root directories '
                'should be listed in this search path.'),
            'local_resource': (
                site_options.local_resource,
                'A directory to store variant tools related '
                'resources such as reference genomes and annotation database. This option will '
                'be ignored if a writable shared resource directory is specified by the system '
                'admin.'),
            'user_stash': (
                site_options.user_stash,
                'A ;-separated list of directories that stores personally-'
                'generated file formats, annotation databases, pipelines and such. These '
                'directories are searched if the requested file is not available in any '
                'online repository (or their local copy). Files under these directories are'
                'NOT maintained by variant tools (no manifest or md5 signatures are monitors).'
                'Currently only local files are allowed (no URL to a remote server). By setting '
                'this option to ~/.variant_tools, resource files under that directory '
                'will be usable even if they are not managed by variant tools')
        }
        # a default value
        self.command_line = ''
        # path to the project cache
        self._cache_dir = os.path.realpath('.vtools_cache')
        #
        self._local_resource = self.persistent_options['local_resource'][0]
        if site_options.shared_resource is None:
            self._shared_resource = self._local_resource
        else:
            self._shared_resource = site_options.shared_resource
            if not os.path.isdir(site_options.shared_resource):
                os.makedirs(site_options.shared_resource)
        #
        self._user_stash = self.persistent_options['user_stash'][0]
        self._logfile_verbosity = self.persistent_options['logfile_verbosity'][
            0]
        self._verbosity = self.persistent_options['verbosity'][0]
        # self._check_update = self.persistent_options['check_update'][0]
        # default sqlite pragma
        self._sqlite_pragma = self.persistent_options['sqlite_pragma'][0]
        # number of processes used for reader under multi-processing mode
        self._import_num_of_readers = self.persistent_options[
            'import_num_of_readers'][0]
        # path to a temporary directory, will be allocated automatically.
        self._temp_dir = self.persistent_options['temp_dir'][0]
        self._proj_temp_dir = self._temp_dir  # per project temp directory
        # how to handle missing data
        self._treat_missing_as_wildtype = self.persistent_options[
            'treat_missing_as_wildtype'][0]
        # association test time out
        self._association_timeout = self.persistent_options[
            'association_timeout'][0]
        # association test number of genotype loaders
        self._associate_num_of_readers = self.persistent_options[
            'associate_num_of_readers'][0]
        # search path
        self._search_path = self.persistent_options['search_path'][0]
        # logger
        self._logger = None
        #
        # a list of lock file that will be removed when the project is killed
        self._lock_files = []
        #
        self._term_width = None
        #
        self._null_input = None

    #
    def lock(self, filename, content=''):
        with open(filename, 'w') as lockfile:
            lockfile.write(content)
        self._lock_files.append(filename)

    def unlock(self, filename, content=''):
        if filename in self._lock_files:
            self._lock_files.remove(filename)
        if not os.path.isfile(filename):
            return
        with open(filename) as lockfile:
            if lockfile.read() != content:
                raise RuntimeError(
                    'Inconsistent lock file. The output file might have been changed by another process.'
                )
        try:
            os.remove(filename)
        except Exception as e:
            self._logger.warning('Failed to remove lock file {}: {}'.format(
                filename, e))

    def unlock_all(self):
        for filename in self._lock_files:
            try:
                os.remove(filename)
            except Exception as e:
                self._logger.warning('Failed to remove lock file {}: {}'.format(
                    filename, e))
        self._lock_files = []

    #
    # attribute check_update
    #
    #def _set_check_update(self, v):
    #    if v in ['1', True, 'T', 'True', 'Y', 'Yes']:
    #        self._check_update = True
    #    else:
    #        self._check_update = False
    #
    #check_update = property(lambda self: self._check_update, _set_check_update)
    #
    # attribute term_width
    def _set_term_width(self, v):
        try:
            self._term_width = int(v)
        except:
            self._term_width = None

    #
    term_width = property(lambda self: self._term_width, _set_term_width)

    #
    # attribute logfile_verbosity
    #
    def _set_logfile_verbosity(self, v):
        if v in ['0', '1', '2']:
            self._logfile_verbosity = v

    #
    logfile_verbosity = property(lambda self: self._logfile_verbosity,
                                 _set_logfile_verbosity)

    #
    #
    # attribute verbosity
    #
    def _set_verbosity(self, v):
        if v in ['0', '1', '2', '3']:
            self._verbosity = v

    #
    verbosity = property(lambda self: self._verbosity, _set_verbosity)

    #
    # attribute pragma
    #
    def _set_sqlite_pragma(self, pragma):
        # 'None' is for backward compatibility
        if pragma in [None, 'None', '']:
            return
        try:
            p = pragma.split(',')
            #
            for item in p:
                if '=' not in str(item):
                    raise ValueError('Invalid pragma {}'.format(item))
            self._sqlite_pragma = pragma
        except:
            sys.stderr.write('Invalid pragma {}\n'.format(pragma))

    #
    sqlite_pragma = property(
        lambda self: self._sqlite_pragma.split(',')
        if self._sqlite_pragma else [], _set_sqlite_pragma)

    #
    # attribute import_num_of_readers
    #
    def _set_import_num_of_readers(self, n):
        try:
            if n is not None:
                int(n)  # test if n is an integer
                self._import_num_of_readers = str(n)
        except:
            sys.stderr.write(
                'Failed to set number of readers to {}\n'.format(n))

    #
    import_num_of_readers = property(
        lambda self: int(self._import_num_of_readers),
        _set_import_num_of_readers)

    #
    # attribute cache_dir, which is not configurable
    #
    def _set_cache_dir(self, path=None):
        if path is not None:
            self._cache_dir = os.path.realpath(os.path.expanduser(path))
        try:
            if not os.path.isdir(self._cache_dir):
                os.makedirs(self._cache_dir)
        except:
            raise RuntimeError('Failed to create cache directory '.format(
                self._cache_dir))

    #
    cache_dir = property(lambda self: self._cache_dir, _set_cache_dir)
    #
    # attribute shared_resource
    #
    shared_resource = property(
        lambda self: os.path.expanduser(self._local_resource)
        if site_options.shared_resource is None else site_options.
        shared_resource, lambda self, path: 0)

    #
    # attribute local_resource
    #
    def _set_local_resource(self, path=None):
        if path is not None:
            self._local_resource = path
        try:
            if not os.path.isdir(os.path.expanduser(self._local_resource)):
                sys.stderr.write(
                    'Creating local resource directory {}\n'.format(
                        self._local_resource))
                os.makedirs(os.path.expanduser(self._local_resource))
        except:
            raise RuntimeError(
                'Failed to create local resource directory '.format(
                    self._local_resource))

    #
    local_resource = property(
        lambda self: os.path.expanduser(self._local_resource),
        _set_local_resource)

    #
    # attribute temp_dir
    #
    def _set_temp_dir(self, path=None):
        # user can explicitly set a path ('None' could be saved by a previous version of vtools)
        if path not in [None, 'None', '']:
            path = os.path.expanduser(path)
            if not os.path.isdir(path):
                raise ValueError(
                    'Temp directory {} does not exist'.format(path))
            if os.path.isdir(path) and (
                (not os.access(path, os.R_OK)) or
                (not os.access(path, os.W_OK)) or
                (os.stat(path).st_mode & stat.S_ISVTX == 512)):
                raise ValueError('Cannot set temporary directory to directory {} because '.format(path) + \
                    'it is not empty or is not writable or deletable. Please clear this directory or use '
                    'command "vtools admin --set_runtime_option temp_dir=DIR" to set it to another path, '
                    'or a random path (empty DIR).')
            self._temp_dir = path
            # create a random subdirectory in this directory
            while True:
                subdir = os.path.join(
                    path, '_tmp_{}'.format(random.randint(1, 1000000)))
                if not os.path.isdir(subdir):
                    if self._proj_temp_dir is not None and os.path.isdir(
                            self._proj_temp_dir):
                        try:
                            shutil.rmtree(env._proj_temp_dir)
                        except:
                            pass
                    self._proj_temp_dir = subdir
                    os.makedirs(subdir)
                    break
        else:
            # the usual case
            if self._temp_dir is None:
                self._proj_temp_dir = tempfile.mkdtemp()
            try:
                if not os.path.isdir(os.path.expanduser(self._proj_temp_dir)):
                    os.makedirs(os.path.expanduser(self._proj_temp_dir))
                while True:
                    subdir = os.path.join(
                        self._proj_temp_dir,
                        '_tmp_{}'.format(random.randint(1, 1000000)))
                    if not os.path.isdir(subdir):
                        if self._proj_temp_dir is not None and os.path.isdir(
                                self._proj_temp_dir):
                            try:
                                shutil.rmtree(env._proj_temp_dir)
                            except:
                                pass
                        self._proj_temp_dir = subdir
                        os.makedirs(subdir)
                        break
            except:
                sys.stderr.write(
                    'Failed to create a temporary directory {}.\n'.format(
                        self._proj_temp_dir))
                self._proj_temp_dir = tempfile.mkdtemp()

    #
    def _get_temp_dir(self):
        if self._proj_temp_dir is None:
            self._set_temp_dir()
        return os.path.expanduser(self._proj_temp_dir)

    #
    temp_dir = property(_get_temp_dir, _set_temp_dir)

    #
    # attribute treat_missing_as_wildtype
    def _set_treat_missing_as_wildtype(self, val):
        if val in [None, 'None', '0', 'False', 'false', 'FALSE']:
            self._treat_missing_as_wildtype = 'False'
        elif val in ['1', 'True', 'TRUE', 'true']:
            self._treat_missing_as_wildtype = 'True'
        else:
            sys.stderr.write(
                'Invalid input ({}) for runtime option treat_missing_as_wildtype\n'
                .format(val))
            self._treat_missing_as_wildtype = 'False'

    #
    treat_missing_as_wildtype = property(
        lambda self: True
        if self._treat_missing_as_wildtype == 'True' else False,
        _set_treat_missing_as_wildtype)

    #
    # attribute association_timeout
    def _set_association_timeout(self, val):
        try:
            if val in ['None', None]:
                self._association_timeout = None
            else:
                # test if val can be converted to int
                int(val)
                self._association_timeout = val
        except:
            pass

    #
    association_timeout = property(
        lambda self: 0 if self._association_timeout is None else int(
            self._association_timeout), _set_association_timeout)

    #
    # attribute associate_num_of_readers
    def _set_associate_num_of_readers(self, val):
        try:
            if val in ['None', None]:
                self._associate_num_of_readers = None
            else:
                # test if val can be converted to int
                int(val)
                self._associate_num_of_readers = val
        except:
            pass

    #
    associate_num_of_readers = property(
        lambda self: 0 if self._associate_num_of_readers is None else int(
            self._associate_num_of_readers), _set_associate_num_of_readers)

    #
    # attribute search_path
    def _set_search_path(self, val):
        if val not in ['None', None]:
            self._search_path = val

    #
    search_path = property(lambda self: self._search_path, _set_search_path)

    #
    # user stash
    def _set_user_stash(self, val):
        if val not in ['None', None]:
            self._user_stash = val

    #
    user_stash = property(lambda self: self._user_stash, _set_user_stash)

    #
    #
    # attribute logger

    def _set_logger(self, logfile=None):
        # create a logger, but shutdown the previous one
        if not hasattr(logging, 'TRACE'):
            logging.TRACE = 5
            logging.addLevelName(logging.TRACE, "TRACE")
        #
        if self._logger is not None:
            self._logger.handlers = []
        self._logger = logging.getLogger()
        self._logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        levels = {
            '0': logging.WARNING,
            '1': logging.INFO,
            '2': logging.DEBUG,
            '3': logging.TRACE,
            None: logging.INFO
        }
        #
        cout.setLevel(levels[self._verbosity])
        cout.setFormatter(
            ColoredFormatter('%(color_levelname)s: %(color_msg)s'))
        self._logger.addHandler(cout)
        self._logger.trace = lambda msg, *args: self._logger._log(
            logging.TRACE, msg, args)
        # output to a log file
        if logfile is not None:
            ch = logging.FileHandler(logfile, mode='a')
            # NOTE: debug informaiton is always written to the log file
            ch.setLevel(levels[self._logfile_verbosity])
            ch.setFormatter(
                logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
            self._logger.addHandler(ch)

    #
    logger = property(lambda self: self._logger, _set_logger)

    def _get_null_input(self):
        if self._null_input is not None:
            return self._null_input
        self._null_input = os.path.join(
            os.path.expanduser(self._local_resource), 'null_input')
        if not os.path.isfile(self._null_input):
            with open(self._null_input, 'w'):
                pass
        return self._null_input

    #
    null_input = property(lambda self: self._get_null_input(), None)


# the singleton object of RuntimeEnvironments
env = RuntimeEnvironments()
# create a default logger without logging to file, this makes sure a logger
# will be usable even when a project is failed to create
env.logger = None
OS_ENV = {
    x: os.pathsep.join(['.', env.cache_dir, os.environ[x]])
    for x in ['PATH', 'LD_LIBRARY_PATH', 'PYTHONPATH', 'PYTHONHOME', 'R_LIBS']
    if x in os.environ
}

SQL_KEYWORDS = set([
    'ADD', 'ALL', 'ALTER', 'ANALYZE', 'AND', 'AS', 'ASC', 'ASENSITIVE',
    'BEFORE', 'BETWEEN', 'BIGINT', 'BINARY', 'BLOB', 'BOTH', 'BY', 'CALL',
    'CASCADE', 'CASE', 'CHANGE', 'CHAR', 'CHARACTER', 'CHECK', 'COLLATE',
    'COLUMN', 'CONDITION', 'CONSTRAINT', 'CONTINUE', 'CONVERT', 'CREATE',
    'CROSS', 'CURRENT_DATE', 'CURRENT_TIME', 'CURRENT_TIMESTAMP',
    'CURRENT_USER', 'CURSOR', 'DATABASE', 'DATABASES', 'DAY_HOUR',
    'DAY_MICROSECOND', 'DAY_MINUTE', 'DAY_SECOND', 'DEC', 'DECIMAL', 'DECLARE',
    'DEFAULT', 'DELAYED', 'DELETE', 'DESC', 'DESCRIBE', 'DETERMINISTIC',
    'DISTINCT', 'DISTINCTROW', 'DIV', 'DOUBLE', 'DROP', 'DUAL', 'EACH', 'ELSE',
    'ELSEIF', 'ENCLOSED', 'ESCAPED', 'EXISTS', 'EXIT', 'EXPLAIN', 'FALSE',
    'FETCH', 'FLOAT', 'FLOAT4', 'FLOAT8', 'FOR', 'FORCE', 'FOREIGN', 'FROM',
    'FULLTEXT', 'GRANT', 'GROUP', 'HAVING', 'HIGH_PRIORITY', 'HOUR_MICROSECOND',
    'HOUR_MINUTE', 'HOUR_SECOND', 'IF', 'IGNORE', 'IN', 'INDEX', 'INFILE',
    'INNER', 'INOUT', 'INSENSITIVE', 'INSERT', 'INT', 'INT1', 'INT2', 'INT3',
    'INT4', 'INT8', 'INTEGER', 'INTERVAL', 'INTO', 'IS', 'ITERATE', 'JOIN',
    'KEY', 'KEYS', 'KILL', 'LEADING', 'LEAVE', 'LEFT', 'LIKE', 'LIMIT', 'LINES',
    'LOAD', 'LOCALTIME', 'LOCALTIMESTAMP', 'LOCK', 'LONG', 'LONGBLOB',
    'LONGTEXT', 'LOOP', 'LOW_PRIORITY', 'MATCH', 'MEDIUMBLOB', 'MEDIUMINT',
    'MEDIUMTEXT', 'MIDDLEINT', 'MINUTE_MICROSECOND', 'MINUTE_SECOND', 'MOD',
    'MODIFIES', 'NATURAL', 'NOT', 'NO_WRITE_TO_BINLOG', 'NULL', 'NUMERIC', 'ON',
    'OPTIMIZE', 'OPTION', 'OPTIONALLY', 'OR', 'ORDER', 'OUT', 'OUTER',
    'OUTFILE', 'PRECISION', 'PRIMARY', 'PROCEDURE', 'PURGE', 'READ', 'READS',
    'REAL', 'REFERENCES', 'REGEXP', 'RELEASE', 'RENAME', 'REPEAT', 'REPLACE',
    'REQUIRE', 'RESTRICT', 'RETURN', 'REVOKE', 'RIGHT', 'RLIKE', 'SCHEMA',
    'SCHEMAS', 'SECOND_MICROSECOND', 'SELECT', 'SENSITIVE', 'SEPARATOR', 'SET',
    'SHOW', 'SMALLINT', 'SONAME', 'SPATIAL', 'SPECIFIC', 'SQL', 'SQLEXCEPTION',
    'SQLSTATE', 'SQLWARNING', 'SQL_BIG_RESULT', 'SQL_CALC_FOUND_ROWS',
    'SQL_SMALL_RESULT', 'SSL', 'STARTING', 'STRAIGHT_JOIN', 'TABLE',
    'TERMINATED', 'THEN', 'TINYBLOB', 'TINYINT', 'TINYTEXT', 'TO', 'TRAILING',
    'TRIGGER', 'TRUE', 'UNDO', 'UNION', 'UNIQUE', 'UNLOCK', 'UNSIGNED',
    'UPDATE', 'USAGE', 'USE', 'USING', 'UTC_DATE', 'UTC_TIME', 'UTC_TIMESTAMP',
    'VALUES', 'VARBINARY', 'VARCHAR', 'VARCHARACTER', 'VARYING', 'WHEN',
    'WHERE', 'WHILE', 'WITH', 'WRITE', 'XOR', 'YEAR_MONTH', 'ZEROFILL',
    'ASENSITIVE', 'CALL', 'CONDITION', 'CONNECTION', 'CONTINUE', 'CURSOR',
    'DECLARE', 'DETERMINISTIC', 'EACH', 'ELSEIF', 'EXIT', 'FETCH', 'GOTO',
    'INOUT', 'INSENSITIVE', 'ITERATE', 'LABEL', 'LEAVE', 'LOOP', 'MODIFIES',
    'OUT', 'READS', 'RELEASE', 'REPEAT', 'RETURN', 'SCHEMA', 'SCHEMAS',
    'SENSITIVE', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE', 'SQLWARNING',
    'TRIGGER', 'UNDO', 'UPGRADE', 'WHILE', 'ABS', 'ACOS', 'ADDDATE', 'ADDTIME',
    'ASCII', 'ASIN', 'ATAN', 'AVG', 'BETWEEN', 'AND', 'BINARY', 'BIN',
    'BIT_AND', 'BIT_OR', 'CASE', 'CAST', 'CEIL', 'CHAR', 'CHARSET', 'CONCAT',
    'CONV', 'COS', 'COT', 'COUNT', 'DATE', 'DAY', 'DIV', 'EXP', 'IS', 'LIKE',
    'MAX', 'MIN', 'MOD', 'MONTH', 'LOG', 'POW', 'SIN', 'SLEEP', 'SORT', 'STD',
    'VALUES', 'SUM'
])


def matchName(pattern, name):
    #
    for char in ('(', ')', '[', ']', '{', '}', '.', '+', '$', '|'):
        pattern = pattern.replace(char, '\\' + char)
    pattern = pattern.replace('?', '.{1}').replace('*', '.*')
    return re.match(pattern, name, re.I)


def validFieldName(name, reserved=[]):
    '''Return a valid field name from a name by converting non-alnum
    characters with _, and add _ if the name starts with a number. If
    the new name is one of reserved, prefix it with _'''
    new_name = re.sub('[\W]+', '_', name.strip())
    if new_name[0].isdigit() or new_name in reserved:
        new_name = '_' + new_name
    return new_name


def decodeTableName(name):
    '''Decode a table to its name that could contain special characters'''
    if name.startswith('_'):
        return binascii.unhexlify(name[1:].encode('utf-8')).decode('utf-8')
    else:
        return name


def encodeTableName(name):
    '''Get a normalized name for variant table. The returned name is a valid
    table name so calling encodeTableName on an encoded name is safe.'''
    # if the table name is not ALPHA + ALPHANUM, use an internal name
    if name.upper() in SQL_KEYWORDS or not name[0].isalpha() \
        or name.startswith('_') or not name.replace('_', '').isalnum():
        return '_' + binascii.hexlify(name.encode('utf-8')).decode('utf-8')
    else:
        return name


def sizeExpr(sz, multiple=1000):
    if sz == 0:
        sz = +0
    SUFFIXES = ["B"
               ] + [i + {
                   1000: "B",
                   1024: "iB"
               }[multiple] for i in "KMGTPEZY"]
    for suffix in SUFFIXES:
        if sz < multiple or suffix == SUFFIXES[-1]:
            if suffix == SUFFIXES[0]:
                return "%d%s" % (sz, suffix)
            else:
                return "%.1f%s" % (sz, suffix)
        else:
            sz /= multiple


#
# Utility functions
#
def lineCount(filename, encoding='UTF-8'):
    '''Estimate the number of lines using file size and line size. This
    function does not attemp to calculate line count exactly because files
    handled by variant tools can be huge. '''
    totalSize = os.path.getsize(filename)
    if totalSize < 500000:
        # small file, read the number of lines directly
        if filename.endswith('.gz'):
            try:
                return len(gzip.open(filename, 'rb').readlines())
            # Python 2.7.4 and 3.3.1 have a regression bug that prevents us from opening
            # certain types of gzip file (http://bugs.python.org/issue17666).
            except TypeError:
                raise RuntimeError(
                    'Failed to open gzipped file {} due to a bug '
                    'in Python 2.7.4 and 3.3.1. Please use a different version '
                    'of Python or decompress this file manually.'.format(
                        filename))
        elif filename.endswith('.bz2'):
            if not bz2_support:
                raise ValueError(
                    'Direct reading of bz2 files is not supported. Please update your python installation or uncompress the file before processing'
                )
            return len(bz2.BZ2File(filename).readlines())
        else:
            return len(open(filename, 'rb').readlines())
    elif filename.endswith('.gz'):
        try:
            input = gzip.open(filename, 'rb')
            input.seek(50000, 0)
            content = input.read(500000).decode(encoding)
            input.close()
            lineCount = len(content.split('\n'))
            input.close()
        # Python 2.7.4 and 3.3.1 have a regression bug that prevents us from opening
        # certain types of gzip file (http://bugs.python.org/issue17666).
        except TypeError:
            raise RuntimeError(
                'Failed to open gzipped file {} due to a bug '
                'in Python 2.7.4 and 3.3.1. Please use a different version '
                'of Python or decompress this file manually.'.format(filename))
        # assuming an arbitrary compression ratio of 5. :-)
        return int(lineCount * (5 * totalSize / 500000.))
    elif filename.endswith('.bz2'):
        if not bz2_support:
            raise ValueError(
                'Direct reading of bz2 files is not supported. Please update your python installation or uncompress the file before processing'
            )
        input = bz2.BZ2File(filename, 'rb')
        input.seek(50000, 0)
        content = input.read(500000).decode(encoding)
        input.close()
        lineCount = len(content.split('\n'))
        input.close()
        # assuming an arbitrary compression ratio of 5. :-)
        return int(lineCount * (5 * totalSize / 500000.))
    else:
        # only binary mode can accomendate end-relative seeks in python 3.
        input = open(filename, 'rb')
        # count from the back because they tend to be records
        # with consistent size
        input.seek(-99000, 2)
        content = input.read()
        input.close()
        lineCount = content.decode(encoding).count('\n')
        input.close()
        return int(lineCount * (totalSize / 99000.))


class PrettyPrinter:

    def __init__(self,
                 delimiter=None,
                 precision=None,
                 na='.',
                 max_width={},
                 cache_size=200):
        ''' delimiter: use specified field to separate fields
            max_width: a dictionary of {col: max_width} to change long
                text to START...END
            cache: only use the first cache lines to get column width
        '''
        self.width = []
        self.rows = []
        self.max_width = max_width
        self.cache_size = cache_size
        self.precision = precision
        if precision is None:
            self.formatter = lambda x: na if x is None else str(x)
        else:
            format_string = '{{:.{}f}}'.format(precision)
            self.formatter = lambda x: na if x is None else (str(
                x) if precision is None else format_string.format(
                    x)) if isinstance(x, float) else str(x)
        # if a delimiter is specified, use it
        if delimiter is not None:
            self.delimiter = delimiter.replace(r'\t', '\t')
            self.write = self.direct_print
            self.write_rest = self.direct_print_rest
        elif max_width:
            self.delimiter = '\t'
            self.write = self.cached_trim_print
            self.write_rest = self.cached_trim_print_rest
        else:
            self.delimiter = '\t'
            self.write = self.cached_print
            self.write_rest = self.cached_print_rest

    #
    # MODE 1: direct print
    #
    def direct_print(self, data):
        '''print data directly using specified delimiter'''
        print((self.delimiter.join([self.formatter(x) for x in data])))

    def direct_print_rest(self):
        '''No cache so do nothing'''
        pass

    #
    # MODE 2: cached, trimmed print
    #
    def cached_trim_print(self, raw_data):
        '''Use cache, figure out column width'''
        data = [self.formatter(x) for x in raw_data]
        trimmed = {}
        for c, m in list(self.max_width.items()):
            if len(data[c]) > m:
                trimmed[c] = data[c][:m //
                                     3] + '...' + data[c][-(m - m // 3 - 3):]
        if trimmed:
            trimmed_data = [x for x in data]
            for c, txt in list(trimmed.items()):
                trimmed_data[c] = txt
        else:
            trimmed_data = data
        #
        self.rows.append(trimmed_data)
        if not self.width:
            self.width = [len(x) for x in trimmed_data]
        else:
            if len(self.width) != len(trimmed_data):
                env.logger.warning(
                    'Inconsistent number of columns ({} and {}), perhaps incorrect header has been specified.'
                    .format(len(self.width), len(trimmed_data)))
                if len(self.width) < len(trimmed_data):
                    self.width.extend([0] *
                                      (len(trimmed_data) - len(self.width)))
            self.width = [
                max(y, len(x)) for y, x in zip(self.width, trimmed_data)
            ]
        # cache size exceeds, use collected width and stop checking
        if len(self.rows) > self.cache_size:
            self.cached_trim_print_rest()
            # change print mode
            self.write = self.uncached_trim_print

    def cached_trim_print_rest(self):
        '''Print and clear cache'''
        if not self.rows:
            return
        # do not ljust the last column. This avoids unnecessary spaces
        # at the end of each line
        self.width[-1] = 0
        # print everything in cache
        print(('\n'.join([
            self.delimiter.join(
                [col.ljust(width)
                 for col, width in zip(row, self.width)])
            for row in self.rows
        ])))
        # clear cache
        self.rows = []

    def uncached_trim_print(self, raw_data):
        data = [self.formatter(x) for x in raw_data]
        trimmed = {}
        for c, m in list(self.max_width.items()):
            if len(data[c]) > m:
                trimmed[c] = data[c][:m //
                                     3] + '...' + data[c][-(m - m // 3 - 3):]
        if trimmed:
            trimmed_data = [x for x in data]
            for c, txt in list(trimmed.items()):
                trimmed_data[c] = txt
        else:
            trimmed_data = data
        #
        print((self.delimiter.join([
            col.ljust(width) for col, width in zip(trimmed_data, self.width)
        ])))

    #
    # MODE 3: cached, untrimmed print
    #
    def cached_print(self, raw_data):
        data = [self.formatter(x) for x in raw_data]
        self.rows.append(data)
        if not self.width:
            self.width = [len(x) for x in data]
        else:
            if len(self.width) != len(data):
                env.logger.warning(
                    'Inconsistent number of columns ({} and {}), perhaps incorrect header has been specified.'
                    .format(len(self.width), len(data)))
                if len(self.width) < len(data):
                    self.width.extend([0] * (len(data) - len(self.width)))
            self.width = [max(y, len(x)) for y, x in zip(self.width, data)]
        # cache size exceeds, use collected width and stop checking
        if len(self.rows) > self.cache_size:
            self.cached_print_rest()
            # change print mode
            self.write = self.uncached_trim_print

    def cached_print_rest(self):
        if not self.rows:
            return
        # do not ljust the last column. This avoids unnecessary spaces
        # at the end of each line
        self.width[-1] = 0
        print(('\n'.join([
            self.delimiter.join(
                [col.ljust(width)
                 for col, width in zip(row, self.width)])
            for row in self.rows
        ])))
        self.rows = []

    def uncached_print(self, data):
        print((self.delimiter.join(
            [col.ljust(width) for col, width in zip(data, self.width)])))


def hasCommand(cmd):
    try:
        fnull = open(os.devnull, 'w')
        result = subprocess.Popen(
            cmd, shell=True, stdout=fnull, stderr=fnull, env=OS_ENV)
        result.terminate()
        fnull.close()
    except OSError:
        # command not found
        return False
    except Exception:
        # other error is OK
        return True
    return True


def runCommand(cmd, instream=None, msg=''):
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    try:
        tc = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=OS_ENV)
        if instream:
            if sys.version_info.major == 3:
                instream = instream.encode(sys.getdefaultencoding())
            out, error = tc.communicate(instream)
        else:
            out, error = tc.communicate()
        if sys.version_info.major == 3:
            out = out.decode(sys.getdefaultencoding())
            error = error.decode(sys.getdefaultencoding())
        if tc.returncode < 0:
            raise ValueError(
                "Command '{0}' was terminated by signal {1}".format(
                    cmd, -tc.returncode))
        elif tc.returncode > 0:
            raise ValueError("{0}".format(error))
        else:
            if error:
                msg = "[WARNING] {0}: {1}".format(msg, error)
                if env.logger is not None:
                    env.logger.debug(msg)
                else:
                    sys.stderr.write(msg + '\n')
    except OSError as e:
        raise OSError("Execution of command '{0}' failed: {1}".format(cmd, e))
    return out


def openFile(filename):
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith(
            '.tgz'):
        raise RuntimeError(
            'Please decompress {} before reading.'.format(filename))
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'rb')
    elif filename.lower().endswith('.bz2'):
        return bz2.BZ2File(filename, 'rb')
    else:
        # text file
        # because readline() from gzip.open will be byte, not string, we should return
        # binary here in order to process them equally in order for things to work
        # correctly under python 3
        return open(filename, 'rb')


def typeOfValues(vals):
    '''Figure out type of values and return INT, FLOAT or VARCHAR(maxLength)'''
    vals = [x for x in vals if x.lower() not in ('na', 'null', 'none', '.', '')]
    if len(vals) == 0:
        # a good default value?
        return 'VARCHAR(255)'
    try:
        list(map(int, vals))
        return 'INT'
    except:
        try:
            list(map(float, vals))
            return 'FLOAT'
        except:
            return 'VARCHAR({})'.format(max([len(x) for x in vals]))


def safeMapFloat(x, nan=True):
    for i, item in enumerate(x):
        try:
            x[i] = float(item)
        except:
            raise
        if not nan and x[i] != x[i]:
            raise
    return x


class delayedAction(object):
    '''Call the passed function with param after a few seconds. It is most often
    used to display certain message only if an action takes a long time.

        action = delayedAction(env.logger.info, 'This might take a while', 5)
        some_action_that_might_take_a_while
        del action

    if the action finishes very quick, the message will not be displayed.
    '''

    def __init__(self, func, param, delay=5):
        self.timer = threading.Timer(delay, func, (param,))

    def __enter__(self):
        self.timer.start()

    def __exit__(self, type, value, traceback):
        self.timer.cancel()


def get_traceback():
    output = StringIO()
    exc_type, exc_value, exc_traceback = sys.exc_info()
    #print "*** print_tb:"
    traceback.print_tb(exc_traceback, limit=1, file=output)
    #print "*** print_exception:"
    try:
        traceback.print_exception(
            exc_type, exc_value, exc_traceback, limit=5, file=output)
    except Exception:
        # the above function call can fail under Python 3.4 for some
        # exception but we do not really care if that happens
        pass
    result = output.getvalue()
    output.close()
    return result


def filesInURL(URL, ext=''):
    '''directory listing of a URL'''
    fh = urllib.request.urlopen(URL)
    files = []
    for line in fh.readlines():
        m = re.search('href="(.*){}"'.format(ext), line.decode())
        if m:
            files.append(m.group(1))
    return files


from array import array
try:
    from fcntl import ioctl
    import termios
except ImportError:
    pass


def getTermWidth():
    if env.term_width is not None:
        return env.term_width
    try:
        h, w = array('h', ioctl(sys.stderr, termios.TIOCGWINSZ, '\0' * 8))[:2]
        return w
    except:
        return 78


class ProgressBar:
    '''A text-based progress bar, it differs from regular progress bar in that
    1. it can start from the middle with init count
    2. it accept update for successful and failed counts
    '''

    def __init__(self, message, totalCount=None, initCount=0,
                 initFailedCount=0):
        if env.verbosity == '0':
            self.update = self.empty
            self.progressBy = self.empty
            self.curlUpdate = self.empty
            self.urllibUpdate = self.empty
            self.sqliteUpdate = self.empty
            self.outputProgress = self.empty
            self.done = self.empty
            self.main = ''
            self.finished = 0
            return
        self.main = message
        self.main_start_time = time.time()
        self.message = self.main
        # get terminal width
        self.handle_resize()
        #
        # It appears that resizing window has caused some threads or processes to
        # stop silently, leading to for example partially imported data. The reason
        # is still unver investigation (might not be related to this) but this feature
        # (resize progress bar dynamically) is temporarily removed.
        #
        #try:
        #    signal.signal(signal.SIGWINCH, self.handle_resize)
        #except:
        #    # signal only works in main thread, so this might not work in all cases
        #    pass
        # total count, including failed ones
        self.count = max(0, initCount)
        self.failed_count = max(0, initFailedCount)
        # total initial count
        self.init_count = initCount
        self.init_failed_count = initFailedCount
        #
        self.finished = 0
        self.reset('', totalCount)

    def handle_resize(self, signum=None, frame=None):
        'Tries to catch resize signals sent from the terminal.'
        self.term_width = getTermWidth()

    def reset(self, msg='', totalCount=None):
        if msg:
            self.message = '{} - {}'.format(self.main, msg)
        self.finished += self.count
        self.count = 0
        self.failed_count = 0
        self.totalCount = totalCount
        self.min_progress_count = None if self.totalCount is None else self.totalCount / 1000
        self.last_progress_count = 0
        self.start_time = None
        self.last_time = None
        self.outputProgress()

    def empty(self, *args, **kwargs):
        return

    def update(self, count, failed_count=0):
        '''completed count jobs, with failed_count failed jobs'''
        if failed_count > count:
            env.logger.warning(
                'Failed count {} greater than completed count {}.'.format(
                    failed_count, count))
            # if there is error ... just give it a number
            failed_count = abs(failed_count)
            count = failed_count
        # do not update if the diferent is less than 0.1% of the total count.
        # this is to avoid excess of calling the time() function
        if self.totalCount is not None and (
                count - self.count) < self.min_progress_count:
            return
        self.count = count
        self.failed_count = failed_count
        self.outputProgress()

    def progressBy(self, count):
        self.last_progress_count += count
        if self.last_progress_count > self.min_progress_count:
            self.count += self.last_progress_count
            self.outputProgress()
            self.last_progress_count = 0

    def curlUpdate(self, total, existing, upload_t, upload_d):
        '''Update called from pycurl'''
        self.count = existing
        self.totalCount = total
        self.outputProgress()

    def urllibUpdate(self, count, blockSize, totalSize):
        '''Update called from urllib'''
        self.count = count * blockSize
        self.totalCount = totalSize
        self.outputProgress()

    def sqliteUpdate(self):
        self.count += 1
        if self.count % 1000 == 0:
            self.outputProgress()

    def outputProgress(self):
        '''Output progress'''
        if not self.start_time:
            self.start_time = time.time()
            self.last_time = self.start_time
        cur_time = time.time()
        # stop update progress bar more than once per second.
        if self.count > 0 and self.count > self.init_count and \
            self.count != self.totalCount and cur_time - self.last_time < 1:
            return
        msg = ['', '', '', '', '', '', '']
        # message
        msg[0] = self.message + ':'
        self.last_time = cur_time
        second_elapsed = cur_time - self.start_time
        if second_elapsed < 0.0001 or self.count == 0:
            msg[4] = ''
        else:
            cps = (self.count - self.init_count) / second_elapsed
            # speed
            if cps > 1000000:
                msg[4] = ' {:.1f}M/s'.format(cps / 1000000)
            elif cps > 1000:
                msg[4] = ' {:.1f}K/s'.format(cps / 1000)
            elif cps > 0.05:
                msg[4] = ' {:.1f}/s'.format(cps)
            elif cps > 1e-6:
                msg[4] = ' {:.1f}s each'.format(1. / cps)
            else:
                msg[4] = ' 0.0/s'
        # estimated time left
        if self.totalCount:
            perc = min(1, float(self.count) / self.totalCount)
            init_perc = min(1, float(self.init_count) / self.totalCount)
            time_left = (second_elapsed / (perc - init_perc) *
                         (1 - perc)) if perc > init_perc else 0
            msg[5] += ' in {}{}'.format(
                '' if time_left < 86400 else '{} day{} '.format(
                    int(time_left / 86400), 's' if time_left > 172800 else ''),
                time.strftime('%H:%M:%S', time.gmtime(time_left)))
        # percentage / progress
        if self.count > 0:
            if self.failed_count == 0:
                # no failed count
                msg[3] = ' {:,}'.format(int(self.count))
                m3Len = len(msg[3])
            else:
                # display failed count in red
                msg[3] = ' {:,}/\033[1;31m{:,}\033[0m'.format(
                    int(self.count), int(self.failed_count))
                m3Len = len(
                    msg[3]
                ) - 11  # the color strings should not be counted as length of message
        else:
            msg[3] = ' '
            m3Len = 1
        if self.totalCount:
            # percentage
            perc = min(1, float(self.count) / self.totalCount)
            failed_perc = min(1, float(self.failed_count) / self.totalCount)
            msg[1] = ' {:5.1f}%'.format(perc * 100)
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(
                msg[4]) - len(msg[5])
            if width > 5:
                front = int((perc - failed_perc) * (width - 4))
                failed_front = int(failed_perc * (width - 4))
                back = width - 4 - front - failed_front
                if failed_front == 0:
                    msg[2] = ' [{}>{}]'.format('=' * front, ' ' * back)
                else:
                    msg[2] = ' [{}\033[1;31m{}\033[0m>{}]'.format(
                        '=' * front, '=' * failed_front, ' ' * back)
        else:
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(
                msg[4])
            msg[6] = ' ' * width
        # use stderr to avoid messing up process output
        sys.stderr.write('\r' + ''.join(msg))

    def done(self, completed=None, failed=None):
        '''Finish, output a new line'''
        if completed is not None:
            self.count = completed
        elif self.totalCount:
            self.count = self.totalCount
        #
        if failed is not None:
            self.falied_count = failed
        #
        msg = ['', '', '', '', '', '']
        # message
        msg[0] = self.main + ':'
        second_elapsed = time.time() - self.main_start_time
        cps = 0 if second_elapsed < 0.0001 else (self.finished +
                                                 self.count) / second_elapsed
        # speed
        if cps > 1000000:
            msg[4] = ' {:.1f}M/s'.format(cps / 1000000)
        elif cps > 1000:
            msg[4] = ' {:.1f}K/s'.format(cps / 1000)
        elif cps > 0.05:
            msg[4] = ' {:.1f}/s'.format(cps)
        elif cps > 1e-6:
            msg[4] = ' {:.1f}s each'.format(1. / cps)
        else:
            msg[4] = ' 0.0/s'
        #
        if self.failed_count == 0:
            msg[3] = ' {:,}'.format(self.finished + self.count)
            m3Len = len(msg[3])
        else:
            msg[3] = ' {:,}/\033[1;31m{:,}\033[0m'.format(
                self.finished + self.count, self.failed_count)
            m3Len = len(msg[3]) - 11
        msg[5] += ' in {}{}'.format(
            '' if second_elapsed < 86400 else '{} day{} '.format(
                int(second_elapsed /
                    86400), 's' if second_elapsed > 172800 else ''),
            time.strftime('%H:%M:%S', time.gmtime(second_elapsed)))
        # percentage / progress
        if self.totalCount:
            # percentage
            msg[1] = ' 100%'
            width = self.term_width - len(msg[0]) - len(msg[1]) - m3Len - len(
                msg[4]) - len(msg[5])
            if width > 4:
                front = int(width - 3)
                if self.count:
                    failed_front = int(
                        float(self.failed_count) / self.count * front)
                else:
                    failed_front = 0
                msg[2] = ' [{}\033[1;31m{}\033[0m]'.format(
                    '=' * (front - failed_front), '=' * failed_front)
        sys.stderr.write('\r' + ''.join(msg) + '\n')
        sys.stderr.flush()


from io import FileIO


class ProgressFileObj(FileIO):
    '''A wrapper of a file object that update a progress bar
    during file read.
    '''

    def __init__(self, prog, *args, **kwargs):
        FileIO.__init__(self, *args, **kwargs)
        self.prog = prog

    def read(self, n, *args):
        self.prog.progressBy(n)
        return FileIO.read(self, n, *args)


def getSnapshotInfo(name):
    '''return meta information for all snapshots'''
    if name.endswith('.tar') or name.endswith('.tar.gz') or name.endswith(
            '.tgz'):
        snapshot_file = name
        mode = 'r' if name.endswith('.tar') else 'r:gz'
    elif name.isalnum():
        snapshot_file = os.path.join(env.cache_dir,
                                     'snapshot_{}.tar'.format(name))
        mode = 'r'
    else:
        raise ValueError(
            'Snapshot name should be a filename with extension .tar, .tgz, or .tar.gz, or a name without any special character.'
        )
    #
    try:
        with tarfile.open(snapshot_file, mode) as snapshot:
            while True:
                tarinfo = snapshot.next()
                if tarinfo.name in ['.snapshot.info', 'README']:
                    readme = snapshot.extractfile(tarinfo)
                    break
                if tarinfo is None:
                    raise ValueError(
                        '{}: cannot find snapshot information'.format(
                            snapshot_file))
            readme.readline().decode()  # header line
            name = readme.readline().decode()[6:].rstrip()  # snapshot name
            date = readme.readline().decode()[6:].rstrip()  # date
            message = ' '.join(
                readme.read().decode()[6:].split('\n'))  # message
            readme.close()
            return (name, date, message)
    except Exception as e:
        env.logger.warning('{}: snapshot read error: {}'.format(
            snapshot_file, e))
        return (None, None, None)


class GenomicRegions(object):
    '''A class to interpret user specified regions. Regions can be
    chr:start-end, chr:end-start, annoDB.field:value,
    chr:start-end,start1-end1, annoDB.field:value1,value2, and their
    union (|), intersection (&), difference(-), and symmetric_difference (^)
    '''

    def __init__(self, regions, zeroBased=False):
        self.raw_regions = regions
        if not self.raw_regions:
            raise ValueError('Empty region is specified.')
        self.proj = None
        self.zeroBased = zeroBased

    def chr_pos_region(self, region):
        # first seprate by ,
        chr, location = region.split(':', 1)
        start, end = location.split('-')
        start = int(start.replace(',', ''))
        end = int(end.replace(',', ''))
        if self.zeroBased:
            start += 1
            end += 1
        if start == 0 or end == 0:
            raise ValueError('0 is not allowed as starting or ending position')
        # start might be after end
        if start > end:
            return (chr[3:] if chr.startswith('chr') else chr, end, start,
                    '(reverse complementary)')
        else:
            return (chr[3:] if chr.startswith('chr') else chr, start, end, '')

    def field_region(self, region):
        if self.proj is None:
            from variant_tools.project import Project
            self.proj = Project()
        # if the regions have been probed before
        regions = self.proj.loadProperty('__region_{}'.format(region), None)
        if regions is not None:
            return eval(regions)
        regions = []
        field, value = region.rsplit(':', 1)
        # what field is this?
        query, fields = consolidateFieldName(self.proj, 'variant', field, False)
        # query should be just one of the fields according to things that are passed
        annoName = query.split('.')[0]
        if query.strip() not in fields:
            # try to run 'vtools use XXXX' if the annotation database if not linked.
            self.proj.close()
            #
            env.logger.info('Linking to {}'.format(annoName))
            try:
                ret = subprocess.call(
                    'vtools use {}'.format(annoName), shell=True)
                if ret:
                    raise RuntimeError(
                        'Could not locate annotation database {} in the project.'
                        .format(annoName))
            except Exception as e:
                raise RuntimeError(
                    'Failed to link to annotation database {}: {}'.format(
                        annoName, e))
            #
            from variant_tools.project import Project
            self.proj = Project()
        # now we have annotation database
        try:
            annoDB = [
                x for x in self.proj.annoDB
                if x.linked_name.lower() == annoName.lower()
            ][0]
        except:
            raise RuntimeError(
                'Could not locate annotation database {} in the project.'
                .format(annoName))
        #
        if annoDB.anno_type != 'range':
            raise ValueError(
                '{} is not linked as a range-based annotation database.'.format(
                    annoDB.linked_name))
        # get the fields?
        chr_field, start_field, end_field = annoDB.build
        #
        # find the regions
        cur = self.proj.db.cursor()
        try:
            cur.execute(
                'SELECT {0},{1},{2} FROM {3}.{4} WHERE {5}="{6}" GROUP BY {0},{1},{2} '
                .format(chr_field, start_field, end_field, annoDB.linked_name,
                        annoDB.name,
                        field.rsplit('.', 1)[-1], value))
        except Exception as e:
            raise ValueError(
                'Failed to search range and comment field: {}'.format(e))
        for idx, (chr, start, end) in enumerate(cur):
            if start > end:
                env.logger.warning(
                    'Ignoring unrecognized region chr{}:{}-{} from {}'.format(
                        chr, start + 1, end + 1, annoDB.linked_name))
                continue
            regions.append(
                (str(chr), int(start), int(end), '{} {}'.format(field,
                                                                idx + 1)))
        if not regions:
            env.logger.warning(
                'No valid chromosomal region is identified for {}'.format(
                    region))
        self.proj.saveProperty('__region_{}'.format(region), str(regions))
        return regions

    def mergeRegions(self, regions):
        while True:
            merged = False
            for i in range(len(regions)):
                for j in range(i + 1, len(regions)):
                    r1 = regions[i]
                    r2 = regions[j]
                    if r1 is None or r2 is None:
                        continue
                    #
                    if r1[0] == r2[0] and r1[2] >= r2[1] and r1[1] <= r2[2]:
                        env.logger.debug(
                            'Merging regions {}:{}-{} ({}) and {}:{}-{} ({})'
                            .format(r2[0], r2[1], r2[2], r2[3], r1[0], r1[1],
                                    r1[2], r1[3]))
                        try:
                            shared_label = [
                                x != y for x, y in zip(r1[3], r2[3])
                            ].index(True)
                        except:
                            # no shared leading string
                            shared_label = 0
                        regions[i] = (r1[0], min(r1[1],
                                                 r2[1]), max(r1[2], r2[2]),
                                      r1[3] + ', ' + r2[3][shared_label:])
                        regions[j] = None
                        merged = True
            if not merged:
                return sorted([x for x in regions if x is not None])

    def eval_regions(self, expr, var_regs):
        #
        # first for all var-regs, we need to expand them to points
        # this can be memory intensive
        for idx in range(len(var_regs)):
            pos = []
            for reg in var_regs[idx]:
                pos.extend([(reg[0], x) for x in range(reg[1], reg[2] + 1)])
            var_regs[idx] = set(pos)
        try:
            positions = eval(expr)
        except Exception as e:
            raise ValueError(
                'Failed to evaluate expression of regions {}: {}'.format(
                    expr, e))
        # re-create regions from positions
        positions = list(positions)
        positions.sort()
        if not positions:
            return []
        regions = [(positions[0][0], positions[0][1], positions[0][1], '')]
        for pos in positions:
            if pos[0] == regions[-1][0] and pos[1] == regions[-1][2] + 1:
                regions[-1][2] += 1
            else:
                regions.append([pos[0], pos[1], pos[1], ''])
        return regions

    def expand(self, proj=None, mergeRegions=True):
        self.proj = proj
        #
        # first, let us identify pieces of the string
        pieces = re.split(
            '(\w+:\d+-\d+(?:,\d+-\d+)*|\w+(-\w+)*\.\w+:[\w.]+(?:,[\w.]+)*)',
            self.raw_regions)
        expr = ''
        var_regs = []
        var_idx = 0
        for piece in [x.strip() for x in pieces if x]:
            if re.match(r'\w+:\d+-\d+(,\d+-\d+)*', piece):
                var_regs.append([])
                chromosome = piece.split(':', 1)[0]
                for reg in piece.split(','):
                    if ':' in reg:
                        var_regs[-1].append(self.chr_pos_region(reg))
                    else:
                        var_regs[-1].append(
                            self.chr_pos_region(chromosome + ':' + reg))
                expr += 'var_regs[{}]'.format(var_idx)
                var_idx += 1
            elif re.match('\w+(-\w+)*\.\w+:[\w.]+(:?,[\w.]+)*', piece):
                var_regs.append([])
                field = piece.split(':', 1)[0]
                for reg in piece.split(','):
                    if ':' in reg:
                        var_regs[-1].extend(self.field_region(reg))
                    else:
                        var_regs[-1].extend(
                            self.field_region(field + ':' + reg))
                expr += 'var_regs[{}]'.format(var_idx)
                var_idx += 1
            elif piece in [',', '|', '&', '^', '-', '(', ')']:
                # treat , as |
                expr += piece.replace(',', '|')
            else:
                raise ValueError(
                    'Incorrect format for regions {}'.format(piece))
        #
        if proj is None and self.proj is not None:
            self.proj.close()
        # if a single expression of regions is specified
        if len(var_regs) == 1:
            regions = var_regs[0]
        else:
            env.logger.debug('Evaluating {}'.format(expr))
            # we have to evaluate an expression to get the regions
            regions = self.eval_regions(expr, var_regs)
        #
        regions = sorted(regions)
        if not regions:
            env.logger.warning('An empty region is specified')
        if regions and mergeRegions:
            regions = self.mergeRegions(regions)
        #env.logger.info('Regions to be simulated ({} bp): {}'.format(
        #    sum([abs(x[2]-x[1])+1 for x in regions]),
        #    ','.join(['{}:{}-{}'.format(x[0], x[1], x[2]) for x in regions])))
        return regions


def expandRegions(regions, proj=None, mergeRegions=True, zeroBased=False):
    return GenomicRegions(regions, zeroBased).expand(proj, mergeRegions)


class ShelfDB:
    '''A sqlite implementation of shelf'''

    def __init__(self, filename, mode='n', lock=None):
        self.filename = filename
        if os.path.isfile(self.filename + '.DB'):
            if mode == 'n':
                os.remove(self.filename + '.DB')
        elif mode == 'r':
            raise ValueError('Temporary database {} does not exist.'.format(
                self.filename))
        self.db = DatabaseEngine()
        self.db.connect(filename, lock=lock)
        self.cur = self.db.cursor()
        self.mode = mode
        if mode == 'n':
            self.cur.execute('CREATE TABLE data (key VARCHAR(255), val TEXT);')
        self.insert_query = 'INSERT INTO data VALUES ({0}, {0});'.format(
            self.db.PH)
        self.select_query = 'SELECT val FROM data WHERE key = {0};'.format(
            self.db.PH)

        if sys.version_info.major >= 3:
            self.add = self._add_py3
            self.get = self._get_py3
        else:
            self.add = self._add_py2
            self.get = self._get_py2

    # python 2 and 3 have slightly different types and methods for pickling.
    def _add_py2(self, key, value):
        # return value from dumps needs to be converted to buffer (bytes)
        self.cur.execute(
            self.insert_query,
            (key,
             memoryview(pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL))))

    def _get_py2(self, key):
        msg = 'Retrieve key {} from ShelfDB'.format(key)
        executeUntilSucceed(self.cur, self.select_query, 5, msg, data=(key,))
        # pickle.loads only accepts string, ...
        return pickle.loads(str(self.cur.fetchone()[0]))

    def _add_py3(self, key, value):
        # return values for dumps is already bytes...
        self.cur.execute(
            self.insert_query,
            (key, pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)))

    def _get_py3(self, key):
        msg = 'Retrieve key {} from ShelfDB'.format(key)
        executeUntilSucceed(self.cur, self.select_query, 5, msg, data=(key,))
        # pickle.loads accepts bytes directly
        return pickle.loads(self.cur.fetchone()[0])

    def close(self):
        if not os.path.isfile(self.filename + '.DB'):
            raise ValueError('Temporary database {} does not exist.'.format(
                self.filename))
        if self.mode == 'n':
            self.db.commit()
            try:
                self.db.execute('CREATE INDEX data_idx ON data (key ASC);')
                self.db.commit()
            except OperationalError as e:
                env.logger.warning('Failed to index temporary database {}: {}. Association tests can still be performed but might be very slow.'.\
                        format(self.filename, e))
                pass
            finally:
                # close the database even if create index failed. In this case
                # shelf retrieval will be slower but still doable
                self.db.close()


def calculateMD5(filename, partial=False):
    filesize = os.path.getsize(filename)
    # calculate md5 for specified file
    md5 = hashlib.md5()
    block_size = 2**20  # buffer of 1M
    try:
        if (not partial) or filesize < 2**26:
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    if not data:
                        break
                    md5.update(data)
        else:
            count = 64
            # otherwise, use the first and last 500M
            with open(filename, 'rb') as f:
                while True:
                    data = f.read(block_size)
                    count -= 1
                    if count == 32:
                        f.seek(-2**25, 2)
                    if not data or count == 0:
                        break
                    md5.update(data)
    except IOError as e:
        sys.exit('Failed to read {}: {}'.format(filename, e))
    return md5.hexdigest()


class ResourceManager:
    '''This class manages a list of resource files managed by
    '''

    def __init__(self):
        self.manifest = {}

    def scanDirectory(self, resource_dir=None, filters=[]):
        '''Returns a manifest for all files under a default or
        specified resource directory. It is used by program manage_resource
        to forcefully generate a manifest file.'''
        if resource_dir is None:
            resource_dir = os.path.expanduser(env.local_resource)
        else:
            resource_dir = os.path.expanduser(resource_dir)
        if not os.path.isdir(resource_dir):
            os.makedirs(resource_dir)
        #
        # go through directories
        filenames = []
        for root, dirs, files in os.walk(resource_dir):
            filenames.extend([(os.path.join(root, x), os.path.getsize(os.path.join(root,x))) for x in files \
                if not filters or all([y in os.path.relpath(os.path.join(root,x), resource_dir) for y in filters])])
        prog = ProgressBar(
            'Scanning {} files under {}'.format(len(filenames), resource_dir),
            sum([min(x[1], 2**26) for x in filenames]))
        total_size = 0
        allowed_directory = [
            'test_data', 'snapshot', 'resource', 'programs', 'pipeline',
            'simulation', 'format', 'annoDB', 'reference'
        ]
        for filename, filesize in filenames:
            if (not any([y in allowed_directory for y in filename.split('/')])) or filename.endswith('.DB') or \
                filename.endswith('.bak') or filename.endswith('.htaccess') or filename.endswith('.log') \
                or filename.endswith('.proj') or '_tmp' in filename:
                env.logger.warning('Ignore {}'.format(filename))
                total_size += min(filesize, 2**26)
            else:
                info = self.addResource(filename, resource_dir)
                total_size += min(info[0], 2**26)
            prog.update(total_size)
        prog.done()

    def writeManifest(self, dest_file=None, URLs=False):
        if dest_file is None:
            if os.path.isdir(env.shared_resource) and \
                os.access(env.shared_resource, os.W_OK):
                dest_file = os.path.join(env.shared_resource,
                                         'MANIFEST_ALL.txt')
            elif os.path.isdir(env.local_resource):
                dest_file = os.path.join(env.local_resource, 'MANIFEST_ALL.txt')
        #
        keys = list(self.manifest.keys())
        keys.sort()
        env.logger.trace('Write manifect to {}'.format(dest_file))
        with open(dest_file, 'w') as manifest:
            for key in keys:
                if URLs:
                    manifest.write(
                        '{0}\t{1[0]}\t{1[1]}\t{1[2]}\t{1[3]}\t{1[4]}\n'.format(
                            key, self.manifest[key]))
                else:
                    manifest.write(
                        '{0}\t{1[0]}\t{1[1]}\t{1[2]}\t{1[3]}\n'.format(
                            key, self.manifest[key]))

    def addResource(self, filename, resource_dir=None):
        if resource_dir is None:
            resource_dir = os.path.expanduser(env.local_resource)
        #
        # if resource_dir is specified, filename
        #
        rel_path = os.path.relpath(filename, resource_dir)
        if rel_path.startswith('.'):
            raise ValueError(
                'Cannot add a resource that is not under the resoure directory {}'
                .format(resource_dir))
        filesize = os.path.getsize(filename)
        md5 = calculateMD5(filename, partial=True)
        refGenome = self.getRefGenome(filename)
        comment = self.getComment(filename).replace('\n',
                                                    ' ').replace('\t',
                                                                 ' ').strip()
        self.manifest[rel_path] = (filesize, md5, refGenome, comment, '')
        return self.manifest[rel_path]

    def getCommentFromConfigFile(self, filename, section, option):
        '''Get comment from annotation description file.'''
        try:
            parser = configparser.ConfigParser()
            parser.read(filename)
            return parser.get(section, option)
        except Exception as e:
            env.logger.warning(
                'Failed to get comment file config file {}: {}'.format(
                    filename, e))
            return ''

    def getRefGenome(self, filename):
        ann_file = filename
        if filename.lower().endswith('.db.gz'):  # annotation database
            if os.path.isfile(filename[:-6] + '.ann'):
                ann_file = filename[:-6] + '.ann'
            else:
                env.logger.warning(
                    'No .ann file could be found for database {}'.format(
                        filename))
                return '*'
        elif filename.endswith('.ann'):
            ann_file = filename
        elif filename.endswith('build36.crr'):
            return 'hg18'
        elif filename.endswith('build37.crr'):
            return 'hg19'
        elif filename.endswith('hg38.crr'):
            return 'hg38'
        elif filename.endswith('.crr'):
            return os.path.basename(filename)[:-4]
        else:
            return '*'
        try:
            parser = configparser.ConfigParser()
            parser.read(ann_file)
            return ','.join([x[0] for x in parser.items('linked fields')])
        except Exception as e:
            env.logger.warning(
                'Failed to get reference genome from .ann file {}: {}'.format(
                    filename, e))
            return '*'

    def getComment(self, filename):
        '''Get the comment from filename according to its type'''
        if filename.lower().endswith('.fmt'):  # file format
            return self.getCommentFromConfigFile(filename, 'format description',
                                                 'description')
        elif filename.lower().endswith('.db.gz'):  # annotation database
            if os.path.isfile(filename[:-6] + '.ann'):
                # we do not have a project .... let us parse .ann file directly
                return self.getCommentFromConfigFile(filename[:-6] + '.ann',
                                                     'data sources',
                                                     'description')
            else:
                return ''
        elif filename.lower().endswith('.ann'):  # annotation
            return self.getCommentFromConfigFile(filename, 'data sources',
                                                 'description')
        elif filename.lower().endswith('.pipeline'):  # pipeline
            return self.getCommentFromConfigFile(filename,
                                                 'pipeline description',
                                                 'description')
        elif filename.lower().endswith('.crr'):  # pipeline
            return 'Reference genome {}'.format(os.path.basename(filename)[:-4])
        elif 'snapshot' in filename and filename.lower().endswith(
                '.tar.gz'):  # snapshot
            (name, date, message) = getSnapshotInfo(filename)
            return '' if message is None else message
        else:  # other files, e.g. crr file
            return ''

    def getLocalManifest(self):
        '''Get a manifest of files from local resource'''
        if os.path.isdir(env.shared_resource) and \
            os.access(env.shared_resource, os.W_OK) and \
            os.path.isfile(os.path.join(env.shared_resource, 'MANIFEST_ALL.txt')):
            # if a local MANIFEST is there in env.local_resource, use it
            manifest_file = os.path.join(env.shared_resource,
                                         'MANIFEST_ALL.txt')
        elif os.path.isdir(env.local_resource) and \
            os.path.isfile(os.path.join(env.local_resource, 'MANIFEST_ALL.txt')):
            # if a local resource is not writable, check user's MANIFEST
            manifest_file = os.path.join(env.local_resource, 'MANIFEST_ALL.txt')
        else:
            # if not, we will have to download a file.
            self.getRemoteManifest()
            return
        #
        self.manifest = {}
        with open(manifest_file, 'r') as manifest:
            for line in manifest:
                filename, sz, md5, refGenome, comment, URLs = line.split(
                    '\t', 5)
                # ref genome might be unsorted
                refGenome = ','.join(sorted(refGenome.split(',')))
                self.manifest[filename] = (int(sz), md5, refGenome,
                                           comment.strip(), eval(URLs))
        # if something wrong with local manifest
        if not self.manifest:
            self.getRemoteManifest()

    def getRemoteManifest(self, URL=None):
        '''Get manifest files from servers (mirrors) and parse it.'''
        try:
            self.manifest = {}
            servers = {}
            for path in env.search_path.split(';'):
                if URL is not None and path != URL:
                    continue
                #
                try:
                    env.logger.trace(
                        'Downloading mirrors.txt from {}'.format(path))
                    if '://' not in path:
                        mirror_file = os.path.join(path, 'MIRRORs.txt')
                    else:
                        mirror_file = downloadURL(
                            path + '/MIRRORS.txt',
                            dest=env.temp_dir,
                            quiet=True)
                    with open(mirror_file) as m:
                        for line in m:
                            if line.startswith('#'):
                                continue
                            fields = line.split()
                            if len(fields) != 2 or (not fields[1].isdigit()
                                                   ) or int(fields[1]) == 0:
                                continue
                            if fields[0] in servers:
                                servers[fields[0]] = max(
                                    servers[fields[0]], int(fields[1]))
                            else:
                                servers[fields[0]] = int(fields[1])
                except Exception as e:
                    env.logger.trace(
                        'Failed to read MIRRORS.txt from {}: {}'.format(
                            path, e))
                    servers[path] = 1
            #
            env.logger.trace('{} mirrors are identified.'.format(len(servers)))
            for server in list(servers.keys()):
                try:
                    # try to download from a local or online repository
                    if '://' not in server:
                        manifest_file = os.path.join(server, 'MANIFEST.txt')
                    else:
                        manifest_file = downloadURL(
                            server + '/MANIFEST.txt',
                            dest=env.temp_dir,
                            quiet=True)
                except Exception as e:
                    env.logger.trace(
                        'Failed to read MANIFEST.txt from {}: {}'.format(
                            server, e))
                    continue
                #
                with open(manifest_file, 'r') as manifest:
                    for line in manifest:
                        filename, sz, md5, refGenome, comment = line.split(
                            '\t', 4)
                        # ref genome might be unsorted
                        refGenome = ','.join(sorted(refGenome.split(',')))
                        if filename not in self.manifest:
                            self.manifest[filename] = [
                                int(sz), md5, refGenome,
                                comment.strip(), [server, servers[server]]
                            ]
                        else:
                            sig = self.manifest[filename]
                            if sig[0] != int(
                                    sz) or sig[1] != md5 or sig[2] != refGenome:
                                env.logger.debug(
                                    '{} in {} have different signature. One of the repositories might need to be updated.'
                                    .format(filename, env.search_path))
                            else:
                                self.manifest[filename][4].extend(
                                    [server, servers[server]])
            self.writeManifest(URLs=True)
        except Exception as e:
            raise RuntimeError(
                'Failed to connect to variant tools resource website: {}'
                .format(e))

    def selectFiles(self, resource_type):
        '''Select files from the remote manifest and see what needs to be downloaded'''
        # if no ceriteria is specified, keep all files
        if resource_type == 'all':
            return
        elif resource_type == 'existing':
            resource_dir = os.path.expanduser(env.shared_resource)
            # go through directories
            filenames = set()
            for root, dirs, files in os.walk(resource_dir):
                filenames |= set([
                    os.path.relpath(os.path.join(root, f), resource_dir)
                    for f in files
                ])
            #
            if not os.access(env.shared_resource, os.W_OK) and \
                env.shared_resource != env.local_resource:
                for root, dirs, files in os.walk(env.local_resource):
                    filenames |= set([
                        os.path.relpath(
                            os.path.join(root, f), env.local_resource)
                        for f in files
                    ])
            #
            self.manifest = {
                x: y for x, y in self.manifest.items() if x in filenames
            }
            return
        elif resource_type == 'format':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('format/')
            }
        elif resource_type == 'snapshot':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('snapshot/')
            }
        elif resource_type in ['annotation', 'annoDB']:
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('annoDB/')
            }
        elif resource_type == 'pipeline':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('pipeline/') and not x.endswith('.py')
            }
        elif resource_type == 'simulation':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('simulation/')
            }
        elif resource_type == 'reference':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if x.startswith('reference/')
            }
        elif resource_type == 'hg18':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if '*' in y[2] or 'hg18' in y[2]
            }
        elif resource_type == 'hg19':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if '*' in y[2] or 'hg19' in y[2]
            }
        elif resource_type == 'hg38':
            self.manifest = {
                x: y
                for x, y in self.manifest.items()
                if '*' in y[2] or 'hg38' in y[2]
            }
        # remove obsolete annotation databases
        if resource_type in ('hg18', 'hg19', 'hg38', 'current', 'annotation'):
            # y[2] is reference genome
            annoDBs = [(x.split('-', 1), y[2])
                       for x, y in self.manifest.items()
                       if x.startswith('annoDB/') and not x.endswith('.ann')]
            # find the latest version of each db
            versions = {}
            for db, refGenome in annoDBs:
                if refGenome not in versions:
                    versions[refGenome] = {}
                if len(db) < 2:
                    # no version
                    versions[refGenome][db[0]] = None
                    continue
                if db[0] in versions[refGenome]:
                    if versions[refGenome][db[0]] < db[1]:
                        versions[refGenome][db[0]] = db[1]
                else:
                    versions[refGenome][db[0]] = db[1]
            # only keep the latest version
            self.manifest = {x:y for x,y in self.manifest.items() if not x.startswith('annoDB/') or x.endswith('.ann') or \
                (x.split('-', 1)[1] if '-' in x else None) == versions[y[2]][x.split('-', 1)[0]]}

    def excludeExistingLocalFiles(self, resource_dir):
        '''Go throughlocal files, check if they are in manifest. If they are
        check if they are identical to remote files'''
        resource_dir = os.path.expanduser(resource_dir)
        # go through directories
        filenames = []
        for root, dirs, files in os.walk(resource_dir):
            for f in files:
                rel_name = os.path.relpath(os.path.join(root, f), resource_dir)
                if rel_name in self.manifest:
                    filenames.append((os.path.join(root, f), rel_name,
                                      os.path.getsize(os.path.join(root, f))))
        prog = ProgressBar(
            'Scanning {} files under {}'.format(len(filenames), resource_dir),
            sum([min(x[2], 2**26) for x in filenames]))
        total_size = 0
        for filename, rel_name, filesize in filenames:
            # if file size are different, will be copied
            if filesize == self.manifest[rel_name][0] and calculateMD5(
                    filename, partial=True) == self.manifest[rel_name][1]:
                self.manifest.pop(rel_name)
            total_size += min(filesize, 2**26)
            prog.update(total_size)
        prog.done()

    def checkUpdate(self, max_updates):
        '''Go through the manifest and download at most max_updates small files
        (.ann, .pipeline etc), and take at most max_updates seconds'''
        changed = []
        for cnt, filename in enumerate(sorted(list(self.manifest.keys()))):
            fileprop = self.manifest[filename]
            dest_dir = os.path.join(env.local_resource,
                                    os.path.split(filename)[0])
            if not os.path.isdir(dest_dir):
                os.makedirs(dest_dir)
            dest_file = os.path.join(env.local_resource, filename)
            #
            if os.path.isfile(dest_file):
                # do not check md5 to increase speed
                if os.path.getsize(dest_file) != fileprop[0]:
                    changed.append(filename)
                continue
            # check all files, but only update small files during casualUpdate
            if filename.rsplit('.', 1)[-1] not in ['ann', 'fmt', 'pipeline']:
                continue
        return changed

    def downloadResources(self, dest_dir=None):
        '''Download resources'''
        for cnt, filename in enumerate(sorted(list(self.manifest.keys()))):
            fileprop = self.manifest[filename]
            #
            try:
                downloaded = downloadFile(
                    filename,
                    dest_dir=dest_dir if dest_dir is None else os.path.join(
                        dest_dir,
                        os.path.split(filename)[0]),
                    checkUpdate=True,
                    quiet=False,
                    message='{}/{} {}'.format(cnt + 1, len(self.manifest),
                                              filename))
                # check md5
                md5 = calculateMD5(downloaded, partial=True)
                if md5 != fileprop[1]:
                    env.logger.error(
                        'Failed to download {}: file signature mismatch.'
                        .format(filename))
            except KeyboardInterrupt as e:
                raise e
            except Exception as e:
                env.logger.error('Failed to download {}: {} {}'.format(
                    filename,
                    type(e).__name__, e))


def compressFile(infile, outfile):
    '''Compress a file from infile to outfile'''
    with open(infile, 'rb') as input, gzip.open(outfile, 'wb') as output:
        buffer = input.read(100000)
        while buffer:
            output.write(buffer)
            buffer = input.read(100000)
    return outfile


def decompressGzFile(filename, inplace=True, force=False, md5=None):
    '''Decompress a file.gz and return file if needed'''
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith(
            '.tar.bz2'):
        dest_files = []
        mode = 'r:gz'
        with tarfile.open(filename, mode) as tar:
            # only extract files
            path = os.path.dirname(filename)
            files = [x.name for x in tar.getmembers() if x.isfile()]
            for f in files:
                dest_file = os.path.join(path, os.path.basename(f))
                dest_files.append(dest_file)
                if not os.path.isfile(dest_file):
                    tar.extract(f, path)
        return dest_files
    elif filename.lower().endswith('.gz'):
        new_filename = filename[:-3]
        # if the decompressed file exists, and is newer than the .gz file, ignore
        if os.path.isfile(new_filename) and not force:
            if md5 is not None and md5 != calculateMD5(
                    new_filename, partial=True):
                env.logger.warning(
                    'MD5 signature mismatch: {} (signature: {} calculated: {})'
                    .format(new_filename, md5,
                            calculateMD5(new_filename, partial=True)))
            else:
                env.logger.debug('Reusing existing decompressed file {}'.format(
                    new_filename))
                return new_filename
        #
        # check if dest_dir is writable
        dest_dir = os.path.dirname(filename)
        if not dest_dir.strip():
            dest_dir = '.'
        if not os.access(dest_dir, os.W_OK):
            # if we are decompressing files from a read-only shared repository
            # write to local_resource
            if os.path.realpath(dest_dir).startswith(
                    os.path.realpath(env.shared_resource)):
                new_filename = '{}/{}'.format(
                    env.local_resource,
                    os.path.realpath(filename)
                    [len(os.path.realpath(env.shared_resource)):-3])
                if os.path.isfile(new_filename) and not force:
                    if md5 is not None and md5 != calculateMD5(
                            new_filename, partial=True):
                        env.logger.warning(
                            'MD5 signature mismatch: {} (signature: {} calculated: {})'
                            .format(new_filename, md5,
                                    calculateMD5(new_filename, partial=True)))
                    else:
                        env.logger.debug(
                            'Reusing existing decompressed file {}'.format(
                                new_filename))
                        return new_filename
            else:
                raise RuntimeError(
                    'Failed to decompress {} to {}: directory not writable'
                    .format(filename, dest_dir))
        # dest_dir can be '' if there is no path for filename
        if dest_dir and not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        #
        env.logger.trace('Decompressing {} to {}'.format(
            filename, new_filename))
        try:
            with gzip.open(filename, 'rb') as input, open(new_filename,
                                                          'wb') as output:
                buffer = input.read(100000)
                while buffer:
                    output.write(buffer)
                    buffer = input.read(100000)
            if md5 is not None and md5 != calculateMD5(
                    new_filename, partial=True):
                env.logger.warning(
                    'MD5 signature mismatch: {} (signature {}, calculated {})'
                    .format(new_filename, md5,
                            calculateMD5(new_filename, partial=True)))
        # Python 2.7.4 and 3.3.1 have a regression bug that prevents us from opening
        # certain types of gzip file (http://bugs.python.org/issue17666).
        except TypeError:
            raise RuntimeError(
                'Failed to open gzipped file {} due to a bug '
                'in Python 2.7.4 and 3.3.1. Please use a different '
                'version of Python or decompress this file manually.'.format(
                    filename))
        #
        if inplace:
            try:
                os.remove(filename)
            except:
                pass
        return new_filename
    else:
        return filename


def TEMP(filename):
    '''Temporary output of filename'''
    # turn path/filename.ext to path/filename_tmp???.ext, where ??? is
    # the process ID to avoid two processes writing to the same temp
    # files. That is to say, if two processes are working on the same step
    # they will produce different temp files, and the final results should
    # still be valid.
    if '.' in os.path.basename(filename):
        return '_tmp{}.'.format(os.getpid()).join(filename.rsplit('.', 1))
    else:
        return '{}_tmp{}'.format(filename, os.getpid())


#
# Well, it is not easy to do reliable download
#
def downloadURL(URL, dest, quiet, message=None):
    # use libcurl? Recommended but not always available
    if 'VTOOLS_ENV' in os.environ and 'NOWEB' in os.environ['VTOOLS_ENV']:
        raise RuntimeError(
            'Failed to download from {}: no internet connection (set by NOWEB in VTOOLS_ENV environment variable)'
            .format(URL))
    env.logger.trace('Download {}'.format(URL))
    filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
    # message during downloading
    if message is None:
        message = filename
    if len(message) > 30:
        message = message[:10] + '...' + message[-16:]
    if os.path.isdir(dest):
        dest = os.path.join(dest, filename)
    #
    try:
        import pycurl
        if not quiet:
            prog = ProgressBar(message)
        dest_tmp = TEMP(dest)
        with open(dest_tmp, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(pycurl.URL, str(URL))
            c.setopt(pycurl.FOLLOWLOCATION, 1)
            c.setopt(pycurl.WRITEFUNCTION, f.write)
            if not quiet:
                c.setopt(pycurl.NOPROGRESS, False)
                c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
            c.perform()
        if not quiet:
            prog.done()
        if c.getinfo(pycurl.HTTP_CODE) == 404:
            try:
                os.remove(dest_tmp)
            except OSError:
                pass
            raise RuntimeError('ERROR 404: Not Found.')
        os.rename(dest_tmp, dest)
        if os.path.isfile(dest):
            return dest
        else:
            raise RuntimeError('Failed to download {} using pycurl'.format(URL))
    except Exception:
        # no pycurl module, or when download fail
        pass
    # use wget? Almost universally available under linux
    if which('wget'):
        # for some strange reason, passing wget without shell=True can fail silently.
        dest_tmp = TEMP(dest)
        p = subprocess.Popen(
            'wget {} -O {} {}'.format('-q' if quiet else '', dest_tmp, URL),
            shell=True)
        ret = p.wait()
        os.rename(dest_tmp, dest)
        if ret == 0 and os.path.isfile(dest):
            return dest
        else:
            try:
                os.remove(dest_tmp)
            except OSError:
                pass
            raise RuntimeError('Failed to download {} using wget'.format(URL))
    else:
        raise RuntimeError(
            'Please install either python module pycurl or program wget for downloading from variant tools repository.'
        )

    # all methods tried
    if os.path.isfile(dest):
        return dest
    # if all failed
    raise RuntimeError('Failed to download {}'.format(URL))


def downloadFile(fileToGet,
                 dest_dir=None,
                 quiet=False,
                 checkUpdate=False,
                 message=None):
    '''Download file from URL to filename.'''
    # two special cases. Move files around to avoid re-download these files.
    if fileToGet == 'reference/hg18.crr':
        if os.path.isfile(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build36.crr')) and \
            not os.path.isfile(os.path.join(env.local_resource, 'reference/hg18.crr')):
            shutil.move(
                os.path.join(
                    env.local_resource,
                    'ftp.completegenomics.com/ReferenceFiles/build36.crr'),
                os.path.join(env.local_resource, 'reference/hg18.crr'))
    elif fileToGet == 'reference/hg19.crr':
        if os.path.isfile(os.path.join(env.local_resource, 'ftp.completegenomics.com/ReferenceFiles/build37.crr')) and \
            not os.path.isfile(os.path.join(env.local_resource, 'reference/hg19.crr')):
            shutil.move(
                os.path.join(
                    env.local_resource,
                    'ftp.completegenomics.com/ReferenceFiles/build37.crr'),
                os.path.join(env.local_resource, 'reference/hg19.crr'))
    #
    # if a complete URL is given, DO NOT download from variant tools repository
    #
    # Downloaded file will look similar to
    #
    # ~/.variant_tools/ftp.completegenomics.com/refgenome/build36.crr
    #
    # unless a specific dest_dir is given. NO md5 check is possible.
    #
    # for backward compatibility, remove http://vtools.houstonbioinformatics.org and
    # use new server
    if fileToGet.startswith('http://vtools.houstonbioinformatics.org/'):
        fileToGet = fileToGet[len('http://vtools.houstonbioinformatics.org/'):]
    #
    if '://' in fileToGet:
        filename = os.path.split(urlparse.urlsplit(fileToGet).path)[-1]
        # get filename from URL
        local_fileToGet = fileToGet.split('://', 1)[1]
        # use root local_resource directory if dest_dir is None
        if dest_dir is not None:
            dest = os.path.join(dest_dir, filename)
            if (not checkUpdate) and os.path.isfile(dest):
                env.logger.trace('Using existing file {}'.format(dest))
                return dest
        else:
            # look for the file in local resource directory
            dest_dir = os.path.join(env.shared_resource,
                                    os.path.split(local_fileToGet)[0])
            dest = os.path.join(env.shared_resource, local_fileToGet)
            # if the file is there, return it directly
            if (not checkUpdate) and os.path.isfile(dest):
                env.logger.trace('Using existing file {}'.format(dest))
                return dest
            # if the shared resource is not writable, write to local_resource
            if not os.access(env.shared_resource, os.W_OK):
                dest_dir = os.path.join(env.local_resource,
                                        os.path.split(local_fileToGet)[0])
                dest = os.path.join(env.local_resource, local_fileToGet)
                # if exists in local user-specific .variant_tools, return it
                if (not checkUpdate) and os.path.isfile(dest):
                    env.logger.trace('Using existing file {}'.format(dest))
                    return dest
        # start to write to it
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        #
        try:
            env.logger.trace('Downloading {} to {}'.format(fileToGet, dest))
            return downloadURL(fileToGet, dest, quiet, message)
        except Exception as e:
            raise ValueError('Failed to download URL {}: {}'.format(
                fileToGet, e))
    #
    # otherwise, download from variant tools repository, but first let us check
    # if the file is in the repository
    #
    filename = os.path.split(fileToGet)[-1]
    local_fileToGet = fileToGet
    resource = ResourceManager()
    resource.getLocalManifest()
    if fileToGet not in resource.manifest:
        # update to the latest manifest and see if we still
        # cannot find the file
        resource.getRemoteManifest()
        if fileToGet not in resource.manifest:
            # look in user stash if avail
            if env.user_stash is not None:
                for us in env.user_stash.split(';'):
                    if not os.path.isdir(os.path.expanduser(us)):
                        env.logger.warning(
                            'Stash directory ({}) does not exist. Check ~/.variant_tools/user_options.py for details.'
                            .format(us))
                    # we usually download file in directories such as annoDB/mydb
                    # so we allow both structure and unstructured stash directory
                    for usf in [
                            os.path.expanduser(os.path.join(us, x))
                            for x in (fileToGet, os.path.basename(fileToGet))
                    ]:
                        if os.path.isfile(usf):
                            return usf
            raise RuntimeError(
                'Failed to download {} because it is not in the variant tools online repository or local stash directories.'
                .format(fileToGet))
    #
    fileSig = resource.manifest[fileToGet]
    if dest_dir is not None:
        dest = os.path.join(dest_dir, os.path.split(filename)[-1])
        # if exists in local user-specific .variant_tools, return it
        if (not checkUpdate) and os.path.isfile(dest):
            env.logger.trace('Using existing file {}'.format(dest))
            if calculateMD5(dest, partial=True) != fileSig[1]:
                env.logger.warning(
                    'MD5 signature mismatch: {} (signature {}, calculated {})'
                    .format(fileToGet, fileSig[1],
                            calculateMD5(dest, partial=True)))
            return dest
    else:
        # look for the file in shared resource directory
        dest_dir = os.path.join(env.shared_resource,
                                os.path.split(local_fileToGet)[0])
        dest = os.path.join(env.shared_resource, local_fileToGet)
        # if the file is there, return it directly
        if (not checkUpdate) and os.path.isfile(dest):
            env.logger.trace('Using existing file {}'.format(dest))
            if calculateMD5(dest, partial=True) != fileSig[1]:
                env.logger.warning(
                    'MD5 signature mismatch: {} (signature {}, calculated {})'
                    .format(fileToGet, fileSig[1],
                            calculateMD5(dest, partial=True)))
            return dest
        # if the share resource is not writable, write to ~/.variant_tools
        if not os.access(env.shared_resource, os.W_OK):
            dest_dir = os.path.join(env.local_resource,
                                    os.path.split(local_fileToGet)[0])
            dest = os.path.join(env.local_resource, local_fileToGet)
            # if exists in local user-specific .variant_tools, return it
            if (not checkUpdate) and os.path.isfile(dest):
                env.logger.trace('Using existing file {}'.format(dest))
                if calculateMD5(dest, partial=True) != fileSig[1]:
                    env.logger.warning(
                        'MD5 signature mismatch: {} (signature {}, calculated {})'
                        .format(fileToGet, fileSig[1],
                                calculateMD5(dest, partial=True)))
                return dest
    #
    if not os.path.isdir(dest_dir):
        os.makedirs(dest_dir)
    #
    # if the file is in the repository, try to find a mirror
    servers = [fileSig[4][2 * i] for i in range(len(fileSig[4]) // 2)]
    weights = [fileSig[4][2 * i + 1] for i in range(len(fileSig[4]) // 2)]
    #
    # if there is a local server, use it regardless of weight
    for server in servers:
        # if the path is a local file repository, do not check for mirrors
        if server.startswith('file://') or '://' not in server:
            source_file = '{}/{}'.format(server, local_fileToGet)
            if source_file.startswith('file://'):
                source_file = source_file[len('file://'):]
            #
            if os.path.isfile(source_file):
                env.logger.trace('Copying {} to {}'.format(source_file, dest))
                shutil.copyfile(source_file, dest)
                if calculateMD5(dest, partial=True) != fileSig[1]:
                    env.logger.warning(
                        'MD5 signature mismatch: {}'.format(dest))
                return dest
            else:
                env.logger.warning(
                    'Cannot locate {} from a local file server {}.'.format(
                        fileToGet, server))
    #
    # no local file server
    while servers:
        if len(servers) == 1:
            idx = 0
        else:
            r = random.random() * sum(weights)
            s = 0
            for idx in range(len(servers)):
                s += weights[idx]
                if s >= r:
                    break
        #
        try:
            env.logger.trace('Download {} from {}'.format(
                fileToGet, servers[idx]))
            downloaded = downloadURL('{}/{}'.format(servers[idx], fileToGet),
                                     dest, quiet, message)
            if calculateMD5(downloaded, partial=True) != fileSig[1]:
                env.logger.warning(
                    'Downloaded file {} is different from remote copy. You might '
                    'want to remove local file and try again to update your local copy.'
                    .format(downloaded))
            return downloaded
        except Exception as e:
            env.logger.warning('Failed to download {} from {}: {}'.format(
                fileToGet, servers[idx], e))
            servers.pop(idx)
            weights.pop(idx)
    # failed to get file
    raise Exception('Failed to download file {}'.format(fileToGet))


class FileInfo:

    def __init__(self, filename):
        self.filename = filename
        self._size = None
        self._md5 = None
        self._mtime = None
        self._ctime = None
        self._first_line = None

    def save(self):
        '''Create a .file_info file with information about the original
        file.'''
        with open(self.filename + '.file_info', 'w') as info:
            info.write('{}\n{}\n{}\n{}\n{}\n'.format(
                os.path.getsize(self.filename),
                calculateMD5(self.filename, partial=True),
                os.path.getctime(self.filename),
                os.path.getmtime(self.filename),
                self._getFirstLine(self.filename)))

    def load(self):
        try:
            with open(self.filename + '.file_info') as info:
                self._size = int(info.readline().strip())
                self._md5 = info.readline().strip()
                self._ctime = float(info.readline().strip())
                self._mtime = float(info.readline().strip())
                try:
                    self._first_line = info.readline().strip()
                except:
                    # early version of file_info does not have first line
                    pass
        except Exception as e:
            raise ValueError('Corrupted file info file {}.file_info: {}'.format(
                self.filename, e))

    def _getFirstLine(self, filename):
        try:
            with openFile(filename) as input:
                return input.readline().decode()
        except:
            return ''

    def md5(self):
        if self._md5 is None:
            if os.path.isfile(self.filename):
                self._md5 = calculateMD5(self.filename, partial=True)
            else:
                self.load()
        return self._md5

    def size(self):
        if self._size is None:
            if os.path.isfile(self.filename):
                self._size = os.path.getsize(self.filename)
            else:
                self.load()
        return self._size

    def mtime(self):
        if self._mtime is None:
            if os.path.isfile(self.filename):
                self._mtime = os.path.getmtime(self.filename)
            else:
                self.load()
        return self._mtime

    def firstline(self):
        if self._first_line is None:
            if os.path.isfile(self.filename):
                self._first_line = self._getFirstLine(self.filename)
            else:
                self.load()
        return self._first_line


def existAndNewerThan(ofiles, ifiles, md5file=None, pipeline=None):
    '''Check if ofiles is newer than ifiles. The oldest timestamp
    of ofiles and newest timestam of ifiles will be used if
    ofiles or ifiles is a list. If a md5file is specified,
    timestamp will be ignored if md5 signature of all ofiles
    and ifiles match.'''
    # if there is no input or output file, ofiles cannot be newer than ifiles.
    if not ofiles or ifiles == ofiles:
        return False
    _ifiles = [ifiles] if not isinstance(ifiles, list) else ifiles
    _ifiles = [x for x in _ifiles if x != env.null_input]
    _ofiles = [ofiles] if not isinstance(ofiles, list) else ofiles
    #
    # file exist?
    for ifile in _ifiles:
        if not (os.path.isfile(ifile) or os.path.isfile(ifile + '.file_info')):
            raise RuntimeError('Input file {} is not found.'.format(ifile))
    # out file does not exist
    if not all([
            os.path.isfile(x) or os.path.isfile(x + '.file_info')
            for x in _ofiles
    ]):
        return False
    #
    # compare timestamp of input and output files
    ifiles_checked = {os.path.realpath(x): False for x in _ifiles}
    md5matched = []
    if md5file:
        nFiles = [0]
        with open(md5file) as md5:
            md5.readline()  # command
            line = md5.readline()
            if not line.startswith('#Start:'):
                env.logger.warning('Invalid exe_info file {}'.format(md5file))
                return False
            for line in md5:
                if line.startswith('#'):
                    if not line.startswith('#End:'):
                        env.logger.warning(
                            'Invalid exe_info file {}'.format(md5file))
                        return False
                    nFiles.append(0)
                    continue
                # stdout and stderr are separated from md5 by newlines
                if not line.strip():
                    break
                try:
                    f_raw, s, m = line.split('\t')
                    f = substituteVars(f_raw, pipeline.VARS, pipeline.GLOBALS)
                    nFiles[-1] += 1
                    s = int(s)
                except Exception as e:
                    env.logger.error('Wrong md5 line {} in {}: {}'.format(
                        line, md5file, e))
                    continue
                # we do not check if f is one of _ifiles or _ofiles because presentation
                # of files might differ
                if not any([
                        os.path.realpath(f) == x
                        for x in list(ifiles_checked.keys())
                ]):
                    if not any([
                            os.path.realpath(f) == os.path.realpath(x)
                            for x in _ofiles
                    ]):
                        env.logger.warning(
                            '{} in exe_info is not an required input or putput file.'
                            .format(f))
                else:
                    ifiles_checked[os.path.realpath(f)] = True
                #
                if not (os.path.isfile(f) or os.path.isfile(f + '.file_info')):
                    env.logger.warning('{} in {} does not exist.'.format(
                        f, md5file))
                    return False
                try:
                    f_info = FileInfo(f)
                    if f_info.size() != s:
                        env.logger.warning(
                            'Size of existing file differ from recorded file: {}'
                            .format(f))
                        return False
                    if f_info.md5() != m.strip():
                        env.logger.warning(
                            'md5 of existing file differ from recorded file: {}'
                            .format(f))
                        return False
                except Exception as e:
                    env.logger.warning(e)
                    return False
                md5matched.append(f)
            #
            if not all(ifiles_checked.values()):
                env.logger.error(
                    'Input or dependent file {} is not recorded in exe_info file.'
                    .format(', '.join(
                        [x for x, y in list(ifiles_checked.items()) if not y])))
                return False
        if len(nFiles) != 2 or nFiles[1] == 0:
            env.logger.warning(
                'Incomplete runtime information for output file {} due to interrupted execution.'
                .format(os.path.basename(md5file).rsplit('.', 1)[0]))
            return False
    #
    def samefile(x, y):
        if x == y:
            return True
        if os.path.isfile(x):
            if os.path.isfile(y):
                return os.path.samefile(x, y)
            elif os.path.isfile(y + '.file_info'):
                return True
            return False
        else:
            if os.path.isfile(y):
                return True
            else:
                return False

    # check if all files have matching signature, do not check timestamp
    if all([any([samefile(x, y) for y in md5matched]) for x in _ifiles]) \
        and all([any([samefile(x, y) for y in md5matched]) for x in _ofiles]):
        return True
    if not _ifiles:
        return True
    # md5 not available
    output_timestamp = min([FileInfo(x).mtime() for x in _ofiles])
    input_timestamp = max([FileInfo(x).mtime() for x in _ifiles])
    if output_timestamp < input_timestamp:
        env.logger.debug('Ignoring older existing output file {}.'.format(
            ', '.join(_ofiles)))
        return False
    else:
        return True


def physicalMemory():
    '''Get the amount of physical memory in the system'''
    # MacOSX?
    import platform
    if platform.platform().startswith('Darwin'):
        # FIXME
        return None
    elif platform.platform().startswith('Linux'):
        try:
            res = subprocess.check_output('free').decode().split('\n')
            return int(res[1].split()[1])
        except Exception:
            return None


def javaXmxCheck(val):
    '''Check if the Xmx option is valid for OPT_JAVA'''
    ram = physicalMemory()
    # cannot check physical memory
    if ram is None:
        return
    # find option matching '-Xmx???'
    m = re.search('-Xmx(\d+)([^\s]*)(?:\s+|$)', val)
    if m is None:  # no -Xmx specified
        return
    try:
        size = int(m.group(1)) * {
            't': 10**9,
            'T': 10**9,
            'g': 10**6,
            'G': 10**6,
            'm': 10**3,
            'M': 10**3,
            '': 1
        }[m.group(2)]
    except:
        sys.exit('Invalid java option {}'.format(val))
    #
    if ram < size:
        sys.exit(
            'Specified -Xms size {} is larger than available physical memory {}'
            .format(size, ram))


def isAnnoDB(annoDB):
    '''Check if a .DB file is a valid annotation database'''
    try:
        db = DatabaseEngine()
        db.connect(os.path.expanduser(annoDB))
        filename = os.path.split(annoDB)[-1]
        if '-' in filename:
            name = filename.split('-')[0]
        else:
            name = os.path.splitext(filename)[0]
        for table in [name, name + '_field', name + '_info']:
            if not db.hasTable(table):
                env.logger.warning(
                    'Skipping invalid annotation database {}. Missing table {}.'
                    .format(annoDB, table))
                return False
        return True
    except:
        return False


#
#
#  reference genome
#
class RefGenome:

    def __init__(self, build):
        if os.path.isfile('{}.crr'.format(build)):
            self.crr = CrrFile('{}.crr'.format(build))
            self.name = build
        elif os.path.isfile(
                os.path.join(env.local_resource, 'reference',
                             '{}.crr'.format(build))):
            self.crr = CrrFile(
                os.path.join(env.local_resource, 'reference',
                             '{}.crr'.format(build)))
            self.name = build
        else:
            try:
                crrFile = downloadFile('reference/{}.crr'.format(build))
                self.crr = CrrFile(crrFile)
                self.name = build
            except Exception as e:
                raise ValueError(
                    'Cannot find reference genome for build {}: {}'.format(
                        build, e))
        #
        self.chrIdx = {}

    def getBase(self, chr, pos):
        try:
            return self.crr.getBase(Location(self.chrIdx[chr], pos - 1))
        except KeyError:
            try:
                self.chrIdx[chr] = self.crr.getChromosomeId(str(chr))
            except Exception as e:
                raise ValueError(
                    'Failed to locate chromosome {} in reference genome {}: e'
                    .format(chr, self.name, e))
            try:
                # ok?
                return self.crr.getBase(Location(self.chrIdx[chr], pos - 1))
            except Exception as e:
                raise ValueError(
                    'Failed to get base {}:{}from reference genome {}: {}'
                    .format(chr, pos, self.name, e))

    def getSequence(self, chr, start, end):
        try:
            return self.crr.getSequence(Range(self.chrIdx[chr], start - 1, end))
        except KeyError:
            try:
                self.chrIdx[chr] = self.crr.getChromosomeId(str(chr))
            except Exception as e:
                raise ValueError(
                    'Failed to locate chromosome {} in reference genome {}: {}'
                    .format(chr, self.name, e))
            # ok?
            try:
                return self.crr.getSequence(
                    Range(self.chrIdx[chr], start - 1, end))
            except Exception as e:
                raise ValueError(
                    'Failed to get sequence {}:{}-{} from reference genome {}: {}'
                    .format(chr, start, end, self.name, e))

    def verify(self, chr, pos, ref):
        try:
            if len(ref) == 1:
                return ref == self.getBase(chr, pos)
            else:
                return ref == self.getSequence(chr, pos, pos + len(ref) - 1)
        except Exception:
            return False


#
#
# Database engine
#


class DatabaseEngine:

    def __init__(self):
        #
        # not connected to any database for now
        self.dbName = None
        self.PH = '?'
        self.database = None

    def describeEngine(self):
        if env.sqlite_pragma == []:
            return 'sqlite (no pragma)'
        else:
            return 'sqlite (with pragma {})'.format(', '.join(
                env.sqlite_pragma))

    #
    # Connection
    #
    def newConnection(self):
        '''Create a new connection from existing configuration'''
        return DatabaseEngine()

    def connect(self, db, readonly=False, lock=None):
        '''Connect to a database'''
        db = os.path.expanduser(db)
        self.dbName = db if (db.endswith('.proj') or
                             db.endswith('.DB')) else db + '.DB'
        self.database = sqlite3.connect(
            self.dbName, check_same_thread=not readonly)
        self.database.enable_load_extension(True)
        #
        # this lock prevents multi process of the same vtools instance, but
        # not for multiple vtools intances...
        if lock is not None:
            lock.acquire()
        while True:
            try:
                self.load_extension(readonly)
                break
            except Exception as e:
                env.logger.warning('Retrying to open database: {}'.format(e))
                time.sleep(10)
        if lock is not None:
            lock.release()

    def load_extension(self, readonly):
        if not readonly:
            # We disable PROGAMA for readonly databases because we often use mutliple readers
            # to read from a readonly database, and applying PRAGMA might cause Operationalerror.
            # We may need to reconsider this though because some pragma applies to
            # readonly databases (e.g. cache_size)
            cur = self.database.cursor()
            for pragma in env.sqlite_pragma:
                # if a pragma is only applicable to certain database, check its name
                if '.' in pragma.split('=')[0] and pragma.split(
                        '.', 1)[0] != self.dbName:
                    continue
                # No error message will be produced for wrong pragma
                # but we may have syntax error.
                try:
                    cur.execute('PRAGMA {}'.format(pragma))
                except Exception as e:
                    # I cannot raise an error because uers need to open the project to reset this value.
                    sys.stderr.write(
                        'Failed to set pragma "{}". Use "vtools admin --set_runtime_option sqlite_pragma=PRAGMA1=VAL,PRAGMA2=VAL" to reset pragmas: {}\n'
                        .format(pragma, e))
                #
                self.database.commit()
        # trying to load extension
        loaded = False
        for path in sys.path:
            ext = glob.glob(os.path.join(path, '_vt_sqlite3_ext*so'))
            if ext:
                cur = self.database.cursor()
                try:
                    cur.execute('SELECT load_extension("{}");'.format(ext[0]))
                except Exception as e:
                    raise SystemError(
                        'Failed to load variant tools sqlite extension from {}: {}'
                        .format(ext[0], e))
                loaded = True
                break
            ext = glob.glob(
                os.path.join(path, 'variant_tools', '_vt_sqlite3_ext*so'))
            if ext:
                cur = self.database.cursor()
                try:
                    cur.execute('SELECT load_extension("{}");'.format(ext[0]))
                except Exception as e:
                    raise SystemError(
                        'Failed to load variant tools sqlite extension from {}: {}'
                        .format(ext[0], e))
                loaded = True
                break
            #
            # pyinstaller bundle this file as 'variant_tools._vt_sqlite3_ext.so'
            ext = glob.glob(
                os.path.join(path, 'variant_tools._vt_sqlite3_ext*so'))
            if ext:
                cur = self.database.cursor()
                try:
                    cur.execute('SELECT load_extension("{}");'.format(ext[0]))
                except Exception as e:
                    raise SystemError(
                        'Failed to load variant tools sqlite extension from {}: {}'
                        .format(ext[0], e))
                loaded = True
                break
        if not loaded:
            env.logger.warning(
                'Failed to load sqlite extension module. No extended SQL functions can be used.'
            )

    def destroy(self):
        self.close()
        try:
            os.remove(self.dbName)
        except:
            self.logger.warning('Failed to remove database file: {}'.format(
                self.dbName))

    def close(self):
        self.database.close()

    def attach(self, db, name=None, lock=None, openExisting=False):
        '''Attach another database to this one. Only needed by sqlite'''
        if db.endswith('.DB') or db.endswith('.proj'):
            db = os.path.expanduser(db)
            if openExisting and not os.path.isfile(db):
                raise ValueError('Database does not exist')
            #if not os.path.isfile(db):
            #    raise RuntimeError('Failed to attach database {}: file does not exist'
            #        .format(db))
            dbName = name if name else os.path.split(db)[-1].split(
                '.')[0].split('-')[0]
            if lock is not None:
                lock.acquire()
            self.execute('''ATTACH DATABASE '{0}' as {1};'''.format(db, dbName))
            for pragma in env.sqlite_pragma:
                if '.' in pragma and pragma.split('.', 1)[0] != dbName:
                    # if pragma is for a specific table with another name, ignore
                    pass
                #
                try:
                    self.execute('PRAGMA {}.{}'.format(dbName, pragma))
                except:
                    pass
            if lock is not None:
                lock.release()
            return dbName
        else:
            db = os.path.expanduser(db)
            if openExisting and not os.path.isfile(db + '.DB'):
                raise ValueError('Database does not exist')
            #if not os.path.isfile(db + '.DB' if db != ':memory:' else db):
            #    raise RuntimeError('Failed to attach database {}: file does not exist'
            #        .format(db))
            dbName = name if name else os.path.split(db)[-1].split(
                '.')[0].split('-')[0]
            if lock is not None:
                lock.acquire()
            self.execute('''ATTACH DATABASE '{0}' as {1};'''.format(
                db + '.DB' if db != ':memory:' else db, dbName))
            for pragma in env.sqlite_pragma:
                # database specific pragma
                if '.' in pragma.split('=')[0]:
                    # pragma for another database
                    if pragma.split('.', 1)[0] != dbName:
                        # if pragma is for a specific table with another name, ignore
                        continue
                    else:
                        # execute it
                        try:
                            self.execute('PRAGMA {}'.format(pragma))
                        except:
                            pass
                else:
                    # apply general pragma to the attached database
                    try:
                        self.execute('PRAGMA {}.{}'.format(dbName, pragma))
                    except:
                        pass
            if lock is not None:
                lock.release()
            return dbName

    def detach(self, db):
        '''Detach a database'''
        self.execute('''DETACH {}'''.format(db))

    def analyze(self):
        '''Analyze a database for better performance'''
        self.execute('analyze;')

    #
    # query
    #
    def execute(self, *args, **kwargs):
        '''A shortcut to get cursor, execute and commit.'''
        cur = self.database.cursor()
        cur.execute(*args, **kwargs)
        return self.database.commit()

    def cursor(self):
        return self.database.cursor()

    def commit(self):
        return self.database.commit()

    #
    # Database
    #
    def hasDatabase(self, db):
        '''Test if a database exists in the current connection'''
        db = os.path.expanduser(db)
        return os.path.isfile(
            db if (db.endswith('.DB') or db.endswith('.proj')) else db + '.DB')

    def removeDatabase(self, db):
        # has to have file extension
        dbFile = db if (db.endswith('.proj') or
                        db.endswith('.DB')) else db + '.DB'
        try:
            os.remove(dbFile)
        except:
            pass
        if os.path.isfile(dbFile):
            sys.exit('Failed to remove database {}'.format(db))

    #
    # Table
    #
    def tables(self, dbName=None):
        '''List all tables in a database'''
        cur = self.database.cursor()
        try:
            if dbName is None:
                cur.execute(
                    "SELECT name FROM sqlite_master WHERE type='table' UNION ALL SELECT name FROM sqlite_temp_master WHERE type='table';"
                )
                return [
                    x[0]
                    for x in cur.fetchall()
                    if not x[0].startswith('sqlite')
                ]
            else:
                cur.execute(
                    "SELECT name FROM {0}.sqlite_master WHERE type='table' UNION ALL SELECT name FROM sqlite_temp_master WHERE type='table';"
                    .format(dbName))
                return [
                    x[0]
                    for x in cur.fetchall()
                    if not x[0].startswith('sqlite')
                ]
        except:
            return []

    def hasTable(self, table):
        '''Test if a table exists in the current database '''
        if '.' not in table:
            return table.lower() in [x.lower() for x in self.tables()]
        else:
            dbName, tableName = table.split('.', 1)
            return tableName.lower() in [x.lower() for x in self.tables(dbName)]

    def hasIndex(self, index):
        '''Test if an index exists in the current database '''
        cur = self.database.cursor()
        if '.' in index:
            db, idx = index.split('.', 1)
            cur.execute(
                "SELECT count(name) FROM {}.sqlite_master WHERE type='index' AND name={};"
                .format(db, self.PH), (idx,))
            return cur.fetchone()[0] > 0
        else:
            cur.execute(
                "SELECT count(name) FROM sqlite_master WHERE type='index' AND name={0} UNION ALL SELECT name FROM sqlite_temp_master WHERE type='index' AND name={0};"
                .format(self.PH), (index, index))
            return cur.fetchone()[0] > 0

    def dropIndex(self, index, table):
        self.execute('DROP INDEX {};'.format(index))

    def removeTable(self, table):
        '''Remove specified table'''
        cur = self.database.cursor()
        cur.execute('DROP TABLE {};'.format(table))
        # FIXME: should we automatically do VACUUM, this can be slow when the table is deletec
        # but can help performance for the creation of new tables.
        # NOTE: It seems that re-generating a table can be VERY slow without vacuum.
        #    cur.execute('VACUUM;'

        self.database.commit()

    def truncateTable(self, table):
        '''Clear all record in a table'''
        cur = self.database.cursor()
        cur.execute('DELETE FROM {};'.format(table))
        self.database.commit()

    def renameTable(self, fromTable, toTable):
        '''Rename a table from fromTable to toTable'''
        cur = self.database.cursor()
        cur.execute('ALTER TABLE {} RENAME TO {};'.format(fromTable, toTable))
        self.database.commit()

    def backupTable(self, table):
        '''Backup a table to table_timestamp'''
        while True:
            new_table = encodeTableName('{}_{}'.format(
                decodeTableName(table),
                time.strftime('%b%d_%H%M%S', time.gmtime())))
            if not self.hasTable(new_table):
                self.renameTable(table, new_table)
                return new_table
            time.sleep(1)

    def fieldsOfTable(self, table):
        '''Get the name and type of fields in a table'''
        cur = self.database.cursor()
        if '.' not in table:
            cur.execute(
                'SELECT sql FROM sqlite_master WHERE UPPER(name) = "{}";'
                .format(table.upper()))
        else:
            db, tbl = table.rsplit('.', 1)
            cur.execute(
                'SELECT sql FROM {}.sqlite_master WHERE UPPER(name) = "{}";'
                .format(db, tbl.upper()))
        try:
            schema = cur.fetchone()[0]
        except:
            raise ValueError('Could not get schema of table {}'.format(table))
        fields = [x.strip() for x in schema.split(',')]
        fields[0] = fields[0].split('(')[1].strip()
        fields[-1] = fields[-1].rsplit(')', 1)[0].strip()
        return [x.split(None, 1) for x in fields]

    def binningRanges(self, build, keys, anno_name):
        '''Create a binning table for a range based annotation
        database for specified build and keys'''
        cur = self.database.cursor()
        tbl = '__rng_' + anno_name + '_' + encodeTableName(
            '_'.join([build] + keys))
        if self.hasTable(tbl):
            return
        cur.execute('SELECT rowid, {} FROM {}'.format(','.join(keys),
                                                      anno_name))
        ranges = cur.fetchall()
        cur.execute(
            'CREATE TABLE {} (bin INT, chr VARCHAR(255), start INT, end INT, range_id INT)'
            .format(tbl))
        insert_query = 'INSERT INTO {0} VALUES ({1}, {1}, {1}, {1}, {1});'.format(
            tbl, self.PH)
        prog = ProgressBar('Binning ranges', len(ranges))
        for idx, (rowid, chr, start, end) in enumerate(ranges):
            try:
                if start > end:
                    if start > end + 1:
                        # start == end in the original database, start > end after adjusting start position to 1-based.
                        env.logger.warnning(
                            'Start position {} greater than ending position {} in database {}'
                            .format(start, end, anno_name))
                    sbin = getMaxUcscBin(start - 1, start)
                    ebin = sbin
                else:
                    sbin = getMaxUcscBin(start - 1, start)
                    ebin = getMaxUcscBin(end - 1, end)
                if sbin > ebin:
                    raise SystemError('Start bin greater than end bin...')
                cur.executemany(insert_query,
                                [(bin, chr, start, end, rowid)
                                 for bin in range(sbin, ebin + 1)])
            except Exception:
                env.logger.warning(
                    'Failed to create bins for range {} - {}'.format(
                        start, end))
            if idx % 100 == 99:
                prog.update(idx + 1)
        prog.done()
        cur.execute(
            'CREATE INDEX {0}_idx ON {0} (bin ASC, chr ASC, range_id ASC);'
            .format(tbl))
        self.database.commit()

    def removeFields(self, table, cols):
        '''Remove fields from a table'''
        if len(cols) == 0:
            return
        cur = self.database.cursor()
        if '.' not in table:
            # for my sqlite, we have to create a new table
            fields = self.fieldsOfTable(table)
            new_fields = [
                '{} {}'.format(x, y)
                for x, y in fields
                if x.lower() not in [z.lower() for z in cols]
            ]
            if len(fields) == len(new_fields):
                raise ValueError(
                    'No field could be removed from table {}'.format(table))
            # rename existing table
            cur.execute('ALTER TABLE {0} RENAME TO _{0}_tmp_;'.format(table))
            # create a new table
            cur.execute('CREATE TABLE {} ('.format(table) +
                        ',\n'.join(new_fields) + ');')
            # insert data back
            cur.execute('INSERT INTO {0} SELECT {1} FROM _{0}_tmp_;'.format(
                table, ','.join([x.split()[0] for x in new_fields])))
            # remove old table
            cur.execute('DROP TABLE _{}_tmp_;'.format(table))
        else:
            db, tbl = table.rsplit('.', 1)
            fields = self.fieldsOfTable(table)
            new_fields = [
                '{} {}'.format(x, y)
                for x, y in fields
                if x.lower() not in [z.lower() for z in cols]
            ]
            if len(fields) == len(new_fields):
                raise ValueError(
                    'No field could be removed from table {}'.format(table))
            # rename existing table
            cur.execute('ALTER TABLE {1}.{0} RENAME TO _{0}_tmp_;'.format(
                tbl, db))
            # create a new table
            cur.execute('CREATE TABLE {1}.{0} ('.format(tbl, db) +
                        ',\n'.join(new_fields) + ');')
            # insert data back
            cur.execute(
                'INSERT INTO {2}.{0} SELECT {1} FROM {2}._{0}_tmp_;'.format(
                    tbl, ','.join([x.split()[0] for x in new_fields]), db))
            # remove old table
            cur.execute('DROP TABLE {1}._{0}_tmp_;'.format(tbl, db))

    def typeOfColumn(self, table, col, simplify=False):
        '''Return type of col in table'''
        fields = self.fieldsOfTable(table)
        for n, t in fields:
            if n.lower() == col.lower():
                if simplify:
                    return 'char' if 'CHAR' in t.upper() else (
                        'int' if 'INT' in t.upper() else 'float')
                else:
                    return t
        raise ValueError('No column called {} in table {}'.format(col, table))

    def numOfRows(self, table, exact=True):
        cur = self.database.cursor()
        if not exact:
            # this is much faster if we do not need exact count
            if '.' in table:
                db, tbl = table.rsplit('.', 1)
                cur.execute(
                    'SELECT seq FROM {}.sqlite_sequence WHERE name = {};'
                    .format(db, self.PH), (tbl,))
            else:
                cur.execute(
                    'SELECT seq FROM sqlite_sequence WHERE name = {};'.format(
                        self.PH), (table,))
            res = cur.fetchone()
            if res is not None:
                return res[0]
        cur.execute('SELECT count(*) FROM {};'.format(table))
        return cur.fetchone()[0]

    def startProgress(self, text):
        self.prog = ProgressBar(text)
        self.database.set_progress_handler(self.prog.sqliteUpdate, 10000)

    def stopProgress(self):
        self.prog.done()
        self.database.set_progress_handler(None, 10000)

    def getHeaders(self, table):
        '''Obtain field names of a table'''
        cur = self.database.cursor()
        try:
            cur.execute('SELECT * FROM {} LIMIT 1;'.format(table))
            return [x[0] for x in cur.description]
        except:
            return None


import token


def consolidateFieldName(proj, table, clause_or_list, alt_build=False):
    '''For input sift_score > 0.5, this function expand it to
    dbNSFP.sift_score > 0.5 and return a list of fields (dbNSFP.sift_score
    in this case). It also change pos to alt_pos if alt_build is true.
    We are using a Python tokenizer here so the result might be wrong.

    If clause is passed as a list of fields, they will be connected by ','.
    However, the list can potentially be changed to reflect for example
    wildcard character expansion in functions such as track('d*.vcf').
    It is therefore highly recommended that you pass a list instead of a
    joint fields.
    '''
    if isinstance(clause_or_list, list):
        clause = ', '.join(clause_or_list)
    else:
        clause = clause_or_list
    tokens = [x for x in tokenize.generate_tokens(StringIO(clause).readline)]
    res = []
    fields = []
    has_ref_query = False
    ref_tokens = {
        'ref_sequence': '__REFSEQ__',
        'mut_sequence': '__MUTSEQ__',
        'vcf_variant': '__VCFVARIANT__',
    }
    has_track_query = False
    track_tokens = {'track': '__TRACK__'}
    has_samples_query = False
    samples_tokens = {'samples': '__SAMPLES__', 'genotype': '__GENOTYPE__'}
    has_in_table_query = False
    in_table_tokens = {'in_table': '__IN_TABLE__'}
    #
    for i in range(len(tokens)):
        before_dot = (i + 1 != len(tokens)) and tokens[i + 1][1] == '.'
        after_dot = i > 1 and tokens[i - 1][1] == '.'
        #
        toktype, toval, _, _, _ = tokens[i]
        # replace chr by alt_chr if using an alternative reference genome.
        if alt_build and toval in ['chr', 'pos'] and not before_dot:
            toval = 'alt_' + toval
        #
        if toval.lower() in ref_tokens:
            toval = ref_tokens[toval.lower()]
            has_ref_query = True
        elif toval.lower() in track_tokens:
            toval = track_tokens[toval.lower()]
            has_track_query = True
        elif toval.lower() in samples_tokens:
            toval = samples_tokens[toval.lower()]
            has_samples_query = True
        elif toval.lower() in in_table_tokens:
            toval = in_table_tokens[toval.lower()]
            has_in_table_query = True
        #
        if toktype == token.NAME and toval.upper() not in SQL_KEYWORDS:
            if before_dot:
                # A.B, does not try to expand A
                res.append((toktype, toval))
            elif after_dot:
                # try to get fields:
                try:
                    for info in proj.linkFieldToTable(
                            '{}.{}'.format(tokens[i - 2][1], toval), table):
                        fields.append(info.field)
                    # if variant field (e.g. variant.chr, do not expand
                    # if annotation field (e.g. linked_name.table_name.field, put in table_name
                    if info.field.count('.') == 2:
                        res.append((toktype, info.field.split('.', 1)[1]))
                    else:
                        # A.B, do not expand
                        res.append((toktype, toval))
                except ValueError as e:
                    res.append((toktype, toval))
                    env.logger.debug(e)
            elif toval.startswith('__'):
                res.append((toktype, toval))
            else:
                # A: try to expand A and identify fields
                try:
                    for info in proj.linkFieldToTable(toval, table):
                        fields.append(info.field)
                    # use expanded field, ONLY the last one should have the expanded fieldname
                    res.append((toktype, info.field))
                except ValueError as e:
                    env.logger.debug(e)
                    res.append((toktype, toval))
        else:
            # fasttrack for symbols or function names
            res.append((toktype, toval))
    # a quick fix for a.b parsed to a .b. :-(
    #
    query = tokenize.untokenize(res).replace(' .', '.')
    # for function ref_base and ref_sequence, we need to add a parameter for
    # reference genome file
    if has_ref_query:
        build = proj.alt_build if alt_build else proj.build
        if os.path.isfile('{}.crr'.format(build)):
            crrFile = '{}.crr'.format(build)
        elif os.path.isfile(
                os.path.join(env.local_resource, 'reference',
                             '{}.crr'.format(build))):
            crrFile = os.path.join(env.local_resource, 'reference',
                                   '{}.crr'.format(build))
        else:
            try:
                crrFile = downloadFile('reference/{}.crr'.format(build))
            except Exception as e:
                raise ValueError(
                    'Cannot find reference genome for build {}: {}'.format(
                        build, e))
        for k, v in list(ref_tokens.items()):
            if v == '__MUTSEQ__':
                query = re.sub(
                    r'{}\s*\('.format(v),
                    " {}('{}', chr, pos, ref, alt, ".format(k, crrFile), query)
            else:
                query = re.sub(r'{}\s*\('.format(v),
                               " {}('{}', ".format(k, crrFile), query)

        # chr and pos will be passed to ref_sequence
        if alt_build:
            fields.append('variant.alt_chr')
            fields.append('variant.alt_pos')
        else:
            fields.append('variant.chr')
            fields.append('variant.pos')
    if has_track_query:
        build = proj.alt_build if alt_build else proj.build

        def handleTrackParams(matchObj):
            try:
                filename = eval(matchObj.group(1).strip())
            except Exception as e:
                raise ValueError(
                    'A filename (quoted string) is needed for the first parameter of function track(), "{}" provided: {}'
                    .format(filename, e))
            param = matchObj.group(2)
            if os.path.isfile(filename) or '://' in filename:
                filenames = [filename]
            else:
                filenames = glob.glob(filename)
            #
            if len(filenames) > 1:
                env.logger.info('Filename "{}" matches {} files: {}'.format(
                    filename, len(filenames), ', '.join(filenames)))
                # this will expand clause_or_list ...
                if isinstance(clause_or_list, list):
                    for idx, cl in enumerate(clause_or_list):
                        if isinstance(cl, str) and filename in cl:
                            clause_or_list[idx] = [
                                clause_or_list[idx].replace(filename, x)
                                for x in filenames
                            ]
            # if the string has option 'matched', we need reference genome information
            if os.path.isfile('{}.crr'.format(build)):
                crrFile = '{}.crr'.format(build)
            elif os.path.isfile(
                    os.path.join(env.local_resource, 'reference',
                                 '{}.crr'.format(build))):
                crrFile = os.path.join(env.local_resource, 'reference',
                                       '{}.crr'.format(build))
            else:
                try:
                    crrFile = downloadFile('reference/{}.crr'.format(build))
                except Exception as e:
                    raise ValueError(
                        'Cannot find reference genome for build {}: {}'.format(
                            build, e))
            if alt_build:
                return ', '.join([
                    "track(variant.alt_chr, variant.alt_pos, variant.ref, variant.alt, '{}', '{}' {})"
                    .format(crrFile, x, param) for x in filenames
                ])
            else:
                return ', '.join([
                    "track(variant.chr, variant.pos, variant.ref, variant.alt, '{}', '{}' {})"
                    .format(crrFile, x, param) for x in filenames
                ])

        # these fields are needed for the track function to execute
        if alt_build:
            fields.append('variant.alt_chr')
            fields.append('variant.alt_pos')
        else:
            fields.append('variant.chr')
            fields.append('variant.pos')
        fields.append('variant.ref')
        fields.append('variant.alt')
        query = re.sub("__TRACK__\s*\(\s*([^,\)]+)([^\)]*)\)",
                       handleTrackParams, query)
    if has_samples_query:
        #
        idListFiles = {}

        def writeIDList(cond=''):
            #
            # cond is a string (e.g. "'1'") and need to be evaluated as a string
            cond = str(eval(cond)) if cond else '1'
            if not cond:
                cond = '1'
            # return a file for condition
            if cond in idListFiles:
                return idListFiles[cond]
            cur = proj.db.cursor()
            # first get all sample names
            try:
                cur.execute(
                    'SELECT sample_id FROM sample WHERE sample_name = {} LIMIT 0,1'
                    .format(proj.db.PH), (cond,))
                res = cur.fetchone()
                return res[0]
            except Exception:
                try:
                    cur.execute(
                        'SELECT sample_id FROM sample, filename '
                        'WHERE sample.file_id = filename.file_id AND ({});'
                        .format(cond))
                except:
                    raise ValueError(
                        'Failed to identify a sample using name or condition "{}"'
                        .format(cond))
                filename = '{}/_sample_id_list_{}.txt'.format(
                    env.cache_dir,
                    binascii.hexlify(cond.encode('utf-8')).decode('utf-8'))
                with open(filename, 'w') as idList:
                    for rec in cur:
                        idList.write('{}\n'.format(rec[0]))
                return filename

        #
        idNameFiles = {}

        def writeSampleIdMap(cond=''):
            # cond is a string (e.g. "'1'") and need to be evaluated as a string
            cond = cond if cond else '1'
            if not cond:
                cond = '1'
            if cond in idNameFiles:
                return idNameFiles[cond]
            cur = proj.db.cursor()
            if cond == '1':
                filename = '{}/_sample_id_map.txt'.format(env.cache_dir)
            else:
                filename = '{}/_sample_id_map_{}.txt'.format(
                    env.cache_dir,
                    binascii.hexlify(cond.encode('utf-8')).decode('utf-8'))
            cur.execute(
                'SELECT sample_id, sample_name FROM sample, filename '
                'WHERE sample.file_id = filename.file_id AND ({});'.format(
                    cond))
            with open(filename, 'w') as idMap:
                for rec in cur:
                    idMap.write('{}\t{}\n'.format(
                        rec[0], rec[1] if rec[1].strip() else '.'))
            return filename

        #
        def handleGenotypeParams(matchObj):
            # optional parameters
            #    1: sample name or filter
            #    2: field to output
            params = matchObj.group(1).strip().split(',', 1)
            # default, all samples
            ret = writeIDList(params[0])
            if len(params) == 0:
                # a filename of IDs
                return ("genotype('{}_genotype.DB', variant.variant_id, '{}')"
                        .format(proj.name, ret))
            elif len(params) == 1:
                if type(ret) == str:
                    # a filename of IDs
                    return (
                        "genotype('{}_genotype.DB', variant.variant_id, '{}')"
                        .format(proj.name, ret))
                else:
                    # a single ID
                    return ("genotype('{}_genotype.DB', variant.variant_id, {})"
                            .format(proj.name, ret))
            elif len(params) == 2:
                if type(ret) == str:
                    # a filename of IDs
                    return (
                        "genotype('{}_genotype.DB', variant.variant_id, '{}', {})"
                        .format(proj.name, ret, params[1]))
                else:
                    # a single ID
                    return (
                        "genotype('{}_genotype.DB', variant.variant_id, {}, {})"
                        .format(proj.name, ret, params[1]))

        #
        def handleSamplesParams(matchObj):
            # optional parameters
            #    1: sample filter
            #    2: genotype filter
            try:
                par_string = matchObj.group(1).strip()
                if par_string:
                    params = eval(par_string).split('&')
                else:
                    params = []
            except Exception as e:
                raise ValueError(
                    'Invalid parameter to function samples: {}, {}'.format(
                        matchObj.group(1), e))
            sample_filter = ''
            for p in params:
                if p.startswith('sample_filter='):
                    sample_filter = p[14:]
            # default, all samples
            ret = writeSampleIdMap(sample_filter)
            # put the rest of the parameter together
            params = '&'.join(
                [x for x in params if not x.startswith('sample_filter=')])
            if not params:
                # no parameter
                return ("samples('{}_genotype.DB', variant.variant_id, '{}')"
                        .format(proj.name, ret))
            else:
                return (
                    "samples('{}_genotype.DB', variant.variant_id, '{}', '{}')"
                    .format(proj.name, ret, params))

        #
        query = re.sub("__SAMPLES__\s*\(([^\)]*)\)", handleSamplesParams, query)
        query = re.sub("__GENOTYPE__\s*\(([^\)]*)\)", handleGenotypeParams,
                       query)
        # variant_id will be passed to samples() and genotype() function and is
        # therefore needed.
        fields.append('variant.variant_id')
    if has_in_table_query:
        query = re.sub("__IN_TABLE__\s*\(", 'in_table(variant_id, ', query)
    #
    return query, fields


def hasGenoInfo(proj, IDs, geno_info):
    '''return a vector of true/false for each geno_info.
    all "true" only if all sample tables has every genotype info'''
    if len(geno_info) == 0:
        return []
    result = [True] * len(geno_info)
    proj.db.attach(proj.name + '_genotype')
    try:
        for table in [
                '{}_genotype.genotype_{}'.format(proj.name, id) for id in IDs
        ]:
            header = [x.lower() for x in proj.db.getHeaders(table)]
            for x in geno_info:
                if x not in header:
                    raise ValueError(x)
    except ValueError as e:
        result[geno_info.index('{}'.format(e))] = False
    return result


def extractField(field):
    '''Extract pos from strings such as pos + 100'''
    if field.isalnum():
        return field
    tokens = [x for x in tokenize.generate_tokens(StringIO(field).readline)]
    for i in range(len(tokens)):
        toktype, toval, _, _, _ = tokens[i]
        if toktype == token.NAME:
            return toval
    raise ValueError('Invalid field name: {}'.format(field))


def splitField(clause):
    '''split chr,alt into ['chr', 'alt'], but keep func(chr, pos) + 1 together.'''
    if not clause.strip():
        return []
    if ',' not in clause:
        return [clause]
    tokens = [x for x in tokenize.generate_tokens(StringIO(clause).readline)]
    fields = []
    cache = []
    level = 0
    for i in range(len(tokens)):
        toktype, tokval, _, _, _ = tokens[i]
        if toktype == token.OP:
            if tokval == '(':
                level += 1
            elif tokval == ')':
                level -= 1
            elif tokval == ',':
                if level == 0:
                    fields.append(''.join(cache))
                    cache = []
                    continue
        cache.append(tokval)
    if cache:
        fields.append(''.join(cache))
    # remove empty fields
    return [x for x in fields if x]


#
# Utility function to calculate bins.
#
# This function implements a hashing scheme that UCSC uses (developed by Jim Kent) to
# take in a genomic coordinate range and return a set of genomic "bins" that your range
# intersects.  I found a Java implementation on-line (I need to find the URL) and I
# simply manually converted the Java code into Python code.

# IMPORTANT: Because this is UCSC code the start coordinates are 0-based and the end
# coordinates are 1-based!!!!!!

# BINRANGE_MAXEND_512M = 512 * 1024 * 1024
# binOffsetOldToExtended = 4681; #  (4096 + 512 + 64 + 8 + 1 + 0)

_BINOFFSETS = (
    512 + 64 + 8 + 1,  # = 585, min val for level 0 bins (128kb binsize)
    64 + 8 + 1,  # =  73, min val for level 1 bins (1Mb binsize)
    8 + 1,  # =   9, min val for level 2 bins (8Mb binsize)
    1,  # =   1, min val for level 3 bins (64Mb binsize)
    0)  # =   0, only val for level 4 bin (512Mb binsize)

#    1:   0000 0000 0000 0001    1<<0
#    8:   0000 0000 0000 1000    1<<3
#   64:   0000 0000 0100 0000    1<<6
#  512:   0000 0010 0000 0000    1<<9

_BINFIRSTSHIFT = 17
# How much to shift to get to finest bin.
_BINNEXTSHIFT = 3
# How much to shift to get to next larger bin.
_BINLEVELS = len(_BINOFFSETS)


#
# IMPORTANT: the start coordinate is 0-based and the end coordinate is 1-based.
#
def getUcscBins(start, end):
    bins = []
    startBin = start >> _BINFIRSTSHIFT
    endBin = (end - 1) >> _BINFIRSTSHIFT
    for i in range(_BINLEVELS):
        offset = _BINOFFSETS[i]
        if startBin == endBin:
            bins.append(startBin + offset)
        else:
            for bin in range(startBin + offset, endBin + offset):
                bins.append(bin)
        startBin >>= _BINNEXTSHIFT
        endBin >>= _BINNEXTSHIFT
    return bins


def getMaxUcscBin(start, end):
    bin = 0
    startBin = start >> _BINFIRSTSHIFT
    endBin = (end - 1) >> _BINFIRSTSHIFT
    for i in range(_BINLEVELS):
        offset = _BINOFFSETS[i]
        if startBin == endBin:
            if startBin + offset > bin:
                bin = startBin + offset
        else:
            for i in range(startBin + offset, endBin + offset):
                if i > bin:
                    bin = i
        startBin >>= _BINNEXTSHIFT
        endBin >>= _BINNEXTSHIFT
    return bin


def executeUntilSucceed(cur, query, attempts, operation_msg, data=None):
    '''try to execute queries a few times before it fails'''
    for attempt in range(attempts):
        try:
            if data:
                cur.execute(query, data)
            else:
                cur.execute(query)
            if attempt != 0:
                env.logger.debug(
                    'Operation "' + operation_msg +
                    '" succeeded after {} attempts'.format(attempt + 1))
            break
        except Exception as e:
            if attempt == attempts - 1:
                env.logger.error('Operation "' + operation_msg + '" failed after {} attempts: {}'.\
                                 format(attempt + 1, e))
                raise
            else:
                time.sleep(1 + attempt + random.random() * 10)


def parenthetic_contents(string):
    """Generate parenthesized contents in string as pairs (level, contents)."""
    stack = []
    for i, c in enumerate(string):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start = stack.pop()
            yield (len(stack), string[start + 1:i])


def longest_parenthetic_content(string):
    """Generate longest parenthesized contents in string"""
    stack = []
    for i, c in enumerate(string):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            start = stack.pop()
            if len(stack) == 0:
                return string[start + 1:i]


def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError:
        if os.path.isdir(path):
            pass
        else:
            raise


def installRPackage(libloc, package):
    try:
        mkdir_p(libloc)
        sys.stderr.write('Installing {0} to {1} ...\n'.format(package, libloc))
        runCommand(['R', '-e', 'install.packages("{0}", lib="{1}", repos="http://cran.stat.ucla.edu")'.\
                    format(package, libloc)])
        runCommand([
            "R", "-e", "library('{1}', lib.loc='{0}')".format(libloc, package)
        ])
    except Exception as e:
        raise ValueError("Cannot auto-install / load R library {1}: {0}".format(
            e, package))
    return libloc


def whereisRPackage(package):
    libloc = 'NULL'
    for item in [
            None,
            os.path.join(env.shared_resource, 'Rlib'),
            os.path.join(env.local_resource, 'Rlib')
    ]:
        try:
            library = "'{}'{}".format(
                package, '' if item is None else ", lib.loc='{0}'".format(item))
            runCommand(["R", "-e", "library({})".format(library)])
            libloc = item
            break
        except:
            continue
    if libloc == 'NULL':
        libloc = installRPackage(
            os.path.join(env.local_resource, 'Rlib'), package)
    return libloc


def flatten(listOfLists):
    "Flatten one level of nesting"
    return list(chain.from_iterable(listOfLists))


def pairwise(x, y):
    return flatten([[(i, j) for j in y] for i in x])


def convertDoubleQuote(x):
    return '"{}"'.format(x.replace('"', "'"))


def make_unique(lst):
    extension = 1
    used = set()
    result = []
    for xname in lst:
        if xname in used:
            xname = '{}V{}'.format(xname, extension)
            extension += 1
        result.append(xname)
        used.add(xname)
    return result


class VariableSubstitutor:

    def __init__(self, text, asString):
        self.text = text
        self.asString = asString

    def var_expr(self, var):
        if type(var) == str:
            # tries to be clever and quote filenames with space
            if os.path.isfile(var) and ' ' in var:
                return "'{}'".format(var)
            else:
                return var
        elif type(var) == list:
            if self.asString:
                return ' '.join([self.var_expr(x) for x in var])
            else:
                return [self.var_expr(x) for x in var]
        else:
            env.logger.debug(
                'Return value of pipeline variable is not string or list of strings: {}'
                .format(var))
            return str(var)

    def _substitute(self, text, PipelineVars, PipelineGlobals):
        if float(PipelineVars['pipeline_format']) <= 1.0:
            # for the first version of pipeline specification file, newlines are
            # replaced with ' '. The newer version (1.0+) keeps newline to faciliate the
            # inclusion of multi-line scripts etc.
            text = ' '.join(text.split())
        # now, find ${}, excluding simple {}, {{}} etc
        pieces = re.split('(\${(?:[^{}]|{[^{}]*{[^}]*}[^{}]*}|{[^}]*})*})',
                          text)
        for idx, piece in enumerate(pieces):
            if piece.startswith('${') and piece.endswith('}'):
                expr = piece[2:-1].strip()
                # remove the "name:" part for the old lambda expression syntax
                has_name_piece = re.match('^([\w\d_, ]*:).*', expr)
                if has_name_piece:
                    expr = expr[len(has_name_piece.group(1)):].strip()
                #
                # evaluate the expression in the dictionary of VARS and GLOBAL
                try:
                    # FIXME: # we should try not to contaminating globals()
                    globals().update(PipelineGlobals)
                    ret_value = eval(expr, globals(), PipelineVars)
                    pieces[idx] = self.var_expr(ret_value)
                except Exception as e:
                    env.logger.warning(
                        'Failed to evaluate expression "{}": {}'.format(
                            expr, e))
        #
        if self.asString:
            if float(PipelineVars['pipeline_format']) <= 1.0:
                # now, join the pieces together, but remove all newlines
                return ' '.join(''.join(pieces).split())
            else:
                return ''.join(pieces)
        else:
            pieces = [x for x in pieces if x]
            if not pieces:
                return ''
            elif len(pieces) == 1:
                return pieces[0]
            else:
                if all([isinstance(x, str) or len(x) <= 1 for x in pieces]):
                    return ''.join([
                        x if isinstance(x, str) else
                        (x[0] if len(x) == 1 else '') for x in pieces
                    ])
                else:
                    raise ValueError(
                        'Variables must be string type (or list of length 1) in variable assignment if text and expressions are mixed: {} evalulated as {}'
                        .format(text, pieces))

    def substituteWith(self, PipelineVars, PipelineGlobals):
        count = 1
        while count < 10:
            if isinstance(self.text, str):
                new_text = self._substitute(self.text, PipelineVars,
                                            PipelineGlobals)
            else:
                new_text = [
                    self._substitute(x, PipelineVars, PipelineGlobals)
                    for x in self.text
                ]
            if new_text == self.text:
                return new_text
            else:
                self.text = new_text
            count += 1
        raise ValueError(
            'Failed to evaluate pipeline varialbe {}. Perhpas the variable is nested.'
            .format(self.text))


def substituteVars(text, PipelineVars, PipelineGlobals, asString=True):
    # if asString is to, the return value is forced to be string
    # Otherwise, the evaluate values are returned.
    return VariableSubstitutor(text,
                               asString).substituteWith(PipelineVars,
                                                        PipelineGlobals)


def determineSexOfSamples(proj, sample_IDs=None):
    '''Determine the sex of samples by checking phenotype sex or gender.
    Only selected samples will be checked if sample_IDs is provided.
    This function returns a dictionary with {ID:sex} where sex is 1 for male
    and 2 for female.
    '''
    # find sex information
    sample_fields = [x[0].lower() for x in proj.db.fieldsOfTable('sample')]
    if 'sex' in sample_fields:
        sex_field = 'sex'
    elif 'gender' in sample_fields:
        sex_field = 'gender'
    else:
        raise ValueError(
            'Failed to determine sex of samples: '
            'A phenotype named sex or gender with values 1/2, M/F or '
            'Male/Female is required.')
    sex_dict = {
        'M': 1,
        1: 1,
        'Male': 1,
        'MALE': 1,
        'F': 2,
        2: 2,
        'Female': 2,
        'FEMALE': 2,
        None: None
    }
    try:
        cur = proj.db.cursor()
        # get sex
        cur.execute('SELECT sample_id, {} FROM sample;'.format(sex_field))
        if sample_IDs is not None:
            return {
                x[0]: sex_dict[x[1]]
                for x in cur.fetchall()
                if x[0] in sample_IDs
            }
        else:
            return {x[0]: sex_dict[x[1]] for x in cur.fetchall()}
    except KeyError:
        raise ValueError(
            'Invalid or missing value detected for field {}. Allowed '
            'values are M/F, 1/2, Male/Female.'.format(sex_field))


class PsudoAutoRegion:

    def __init__(self, chrom, build):
        if build == ['hg18', 'build36'] and chrom.lower() in ['x', '23']:
            self.check = self.checkChrX_hg18
        elif build in ['hg18', 'build36'] and chrom.lower() in ['y', '24']:
            self.check = self.checkChrY_hg18
        elif build in ['hg19', 'build37'] and chrom.lower() in ['x', '23']:
            self.check = self.checkChrX_hg19
        elif build in ['hg19', 'build37'] and chrom.lower() in ['y', '24']:
            self.check = self.checkChrY_hg19
        else:
            env.logger.warning(
                'Checking psudo-autosomal regions for build {} on chromosome {} is not supported'
                .format(build, chrom))
            self.check = self.notWithinRegion

    def checkChrX_hg18(self, pos):
        return (pos >= 1 and pos <= 2709520) or \
            (pos >= 154584238 and pos <= 154913754)

    def checkChrY_hg18(self, pos):
        return (pos >= 1 and pos <= 2709520) or \
            (pos >= 57443438 and pos <= 57772954)

    def checkChrX_hg19(self, pos):
        return (pos >= 60001 and pos <= 2699520) or \
            (pos >= 154931044 and pos <= 155270560)

    def checkChrY_hg19(self, pos):
        return (pos >= 10001 and pos <= 2649520) or \
            (pos >= 59034050 and pos <= 59373566)

    def notWithinRegion(self, pos):
        return False


def withinPseudoAutoRegion(chrom, pos, build):
    # return True if position is in autosomal or pseudo-autosomal regions
    # the position information are based on personal communication with
    # Dr. Bert Overduin from 1000 genomes
    return PsudoAutoRegion(chrom, build).check(pos)


# 1000 genomes record pseduo-autosomal regions on chromosome X, and
# record genotypes as homozygotes if they appear on PAR1 and PAR2 of
# both regions. Anyway, the following code removes variants within
# these regions and treat them as autosome variants.


# 23 -- X
# 24 -- Y
# 25 -- MT
# XY -- pseduo-autosomal region
#
# NOTE: Some pipelines use 24 for XY... this can be a mess
def getVariantsOnChromosomeX(proj, variant_table='variant'):
    cur = proj.db.cursor()
    if variant_table == 'variant':
        cur.execute(
            "SELECT variant_id, pos FROM variant WHERE chr in ('X', 'x', '23')")
    else:
        cur.execute(
            "SELECT {0}.variant_id, pos FROM {0}, variant WHERE {0}.variant_id "
            "= variant.variant_id AND variant.chr in ('X', 'x', '23')".format(
                variant_table))
    #
    var_chrX = set(cur.fetchall())
    nPrev = len(var_chrX)
    if nPrev > 0:
        paRegion = PsudoAutoRegion('X', proj.build)
        var_chrX = [x for x in var_chrX if not paRegion.check(x[1])]
        if nPrev > len(var_chrX):
            env.logger.info(
                '{} variants in pseudo-autosomal regions on '
                'chromosome X are treated as autosomal variants.'.format(
                    nPrev - len(var_chrX)))
    return set([x[0] for x in var_chrX])


def getVariantsOnChromosomeY(proj, variant_table='variant'):
    cur = proj.db.cursor()
    if variant_table == 'variant':
        cur.execute(
            "SELECT variant_id, pos FROM variant WHERE chr in ('Y', 'y', '24')")
    else:
        cur.execute(
            "SELECT {0}.variant_id, pos FROM {0}, variant WHERE {0}.variant_id "
            "= variant.variant_id AND chr in ('Y', 'y', '24')".format(
                variant_table))
    #
    var_chrY = set(cur.fetchall())
    nPrev = len(var_chrY)
    if nPrev > 0:
        paRegion = PsudoAutoRegion('Y', proj.build)
        var_chrY = [x for x in var_chrY if not paRegion.check(x[1])]
        if nPrev > len(var_chrY):
            env.logger.info(
                '{} variants in pseudo-autosomal regions on '
                'chromosome Y are treated as autosomal variants.'.format(
                    nPrev - len(var_chrY)))
    return set([x[0] for x in var_chrY])


def getVariantsOnManifolds(proj, variant_table='variant'):
    cur = proj.db.cursor()
    if variant_table == 'variant':
        cur.execute(
            "SELECT variant_id FROM variant WHERE chr NOT IN ('1', '2', "
            "'3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', "
            "'14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y',"
            "'x', 'y', 'XY', 'xy', '23', '24')")
    else:
        cur.execute(
            "SELECT {0}.variant_id FROM {0}, variant WHERE {0}.variant_id "
            "= variant.variant_id AND chr NOT IN ('1', '2', "
            "'3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', "
            "'14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y',"
            "'x', 'y', 'XY', 'xy', '23', '24')".format(variant_table))
    var_chrOther = set([x[0] for x in cur.fetchall()])
    #
    env.logger.debug('{} variants on other chromosomes are identifield'.format(
        len(var_chrOther)))
    return var_chrOther


def call_sex(dat):
    # make sex calls based on very apparent information
    # dat = [chr, pos, gt]
    # this assumes coding for homozygouse on chrX is '2'
    sex = 'unknown'
    xhomo = False
    for locus in dat:
        # FIXME: build is fixed to hg19
        if withinPseudoAutoRegion(locus[0], int(locus[1]), 'hg19'):
            continue
        # call 'M' if '1' on Y chromosome is observed
        # FIXME: will be problematic if XY is coded 24 instead
        if locus[0].lower() in ['y', '24']:
            sex = 'M'
            break
        if locus[0].lower() in ['x', '23']:
            if locus[2] == '2':
                xhomo = True
    # call 'F' if both '0' and '2' are observed in gt and '1' on Y chromosome not observed
    if xhomo and sex == 'unknown':
        sex = 'F'
    return sex


codon_table = {
    'TTT': 'F',
    'TTC': 'F',
    'TCT': 'S',
    'TCC': 'S',
    'TAT': 'Y',
    'TAC': 'Y',
    'TGT': 'C',
    'TGC': 'C',
    'TTA': 'L',
    'TCA': 'S',
    'TAA': '*',
    'TGA': '*',
    'TTG': 'L',
    'TCG': 'S',
    'TAG': '*',
    'TGG': 'W',
    'CTT': 'L',
    'CTC': 'L',
    'CCT': 'P',
    'CCC': 'P',
    'CAT': 'H',
    'CAC': 'H',
    'CGT': 'R',
    'CGC': 'R',
    'CTA': 'L',
    'CTG': 'L',
    'CCA': 'P',
    'CCG': 'P',
    'CAA': 'Q',
    'CAG': 'Q',
    'CGA': 'R',
    'CGG': 'R',
    'ATT': 'I',
    'ATC': 'I',
    'ACT': 'T',
    'ACC': 'T',
    'AAT': 'N',
    'AAC': 'N',
    'AGT': 'S',
    'AGC': 'S',
    'ATA': 'I',
    'ACA': 'T',
    'AAA': 'K',
    'AGA': 'R',
    'ATG': 'M',
    'ACG': 'T',
    'AAG': 'K',
    'AGG': 'R',
    'GTT': 'V',
    'GTC': 'V',
    'GCT': 'A',
    'GCC': 'A',
    'GAT': 'D',
    'GAC': 'D',
    'GGT': 'G',
    'GGC': 'G',
    'GTA': 'V',
    'GTG': 'V',
    'GCA': 'A',
    'GCG': 'A',
    'GAA': 'E',
    'GAG': 'E',
    'GGA': 'G',
    'GGG': 'G'
}

# codon on the reverse strand, but read on the forward-strand in forward direction
codon_table_reverse_complement = {
    'AAG': 'L',
    'CTA': '*',
    'TGT': 'T',
    'TTT': 'K',
    'GAT': 'I',
    'GTT': 'N',
    'TAT': 'I',
    'CCT': 'R',
    'AGG': 'P',
    'AGT': 'T',
    'GCT': 'S',
    'CTT': 'K',
    'TCT': 'R',
    'ATG': 'H',
    'ATT': 'N',
    'AAT': 'I',
    'CAG': 'L',
    'TAG': 'L',
    'GAG': 'L',
    'GTG': 'H',
    'CCA': 'W',
    'CGG': 'P',
    'ACT': 'S',
    'TGG': 'P',
    'TTG': 'Q',
    'GGG': 'P',
    'ATA': 'Y',
    'ACC': 'G',
    'ACA': 'C',
    'TCG': 'R',
    'CTG': 'Q',
    'AGA': 'S',
    'ATC': 'D',
    'CCG': 'R',
    'AAA': 'F',
    'GCA': 'C',
    'CCC': 'G',
    'TCA': '*',
    'TCC': 'G',
    'TTA': '*',
    'CGT': 'T',
    'GTA': 'Y',
    'GAA': 'F',
    'CGA': 'S',
    'TAA': 'L',
    'CAA': 'L',
    'GGA': 'S',
    'GGT': 'T',
    'TGA': 'S',
    'TGC': 'A',
    'TAC': 'V',
    'GGC': 'A',
    'GAC': 'V',
    'GCC': 'G',
    'CGC': 'A',
    'CAC': 'V',
    'CTC': 'E',
    'AAC': 'V',
    'AGC': 'A',
    'GTC': 'D',
    'ACG': 'R',
    'TTC': 'E',
    'CAT': 'M',
    'GCG': 'R'
}

complement_table = {
    'A': 'T',
    'T': 'A',
    'a': 't',
    't': 'a',
    'G': 'C',
    'C': 'G',
    'g': 'c',
    'c': 'g',
    'N': 'N',
    'n': 'n'
}


def genesInRegions(regions, proj):
    '''Locate isoforms (refGene.name) that overlap with any part of specified
    regions. '''
    #
    if 'refgene' not in [x.linked_name.lower() for x in proj.annoDB]:
        # try to run 'vtools use XXXX' if the annotation database if refGene not linked.
        proj.close()
        #
        env.logger.info('Linking to refGene')
        try:
            ret = subprocess.call('vtools use refGene', shell=True)
            if ret:
                raise RuntimeError(
                    'Could not locate annotation database refGene in the project.'
                )
        except Exception as e:
            raise RuntimeError(
                'Failed to link to annotation database refGene: {}'.format(e))
        #
        from variant_tools.project import Project
        proj = Project()
    #
    cur = proj.db.cursor()
    # if a record has been processed
    isoforms = []
    for (chr, start, end, comment) in regions:
        try:
            cur.execute(
                'SELECT name FROM refGene WHERE chr = ? AND txStart <= ? AND txEnd >= ?',
                (chr, end, start))
        except Exception as e:
            raise RuntimeError(
                'Failed to retrieve gene information from refGene database. '
                'Please use "vtools admin --update_resource existing" to update to the latest version '
                'and command "vtools use refGene" to link the refGene database: {}'
                .format(e))
        isoforms.extend([x[0] for x in cur.fetchall()])
    return sorted(list(set(isoforms)))


def dissectGene(gene, proj):
    if 'refgene' not in [x.linked_name.lower() for x in proj.annoDB]:
        # try to run 'vtools use XXXX' if the annotation database if refGene not linked.
        proj.close()
        #
        env.logger.info('Linking to refGene')
        try:
            ret = subprocess.call('vtools use refGene', shell=True)
            if ret:
                raise RuntimeError(
                    'Could not locate annotation database refGene in the project.'
                )
        except Exception as e:
            raise RuntimeError(
                'Failed to link to annotation database refGene: {}'.format(e))
        #
        from variant_tools.project import Project
        proj = Project()
    cur = proj.db.cursor()
    try:
        cur.execute(
            'SELECT name, chr, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, name2 '
            'FROM refGene WHERE name=?', (gene,))
    except Exception as e:
        raise RuntimeError(
            'Failed to retrieve gene information from refGene database. '
            'Please use "vtools admin --update_resource existing" to update to the latest version '
            'and command "vtools use refGene" to link the refGene database: {}'
            .format(e))
    rec = cur.fetchall()
    if len(rec) > 1:
        env.logger.warning(
            'Multiple records was found for isoform {}. Only one of them is used.'
            .format(gene))
    name, chr, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds, name2 = rec[
        0]
    #
    exonStarts = [int(x) for x in exonStarts.split(',') if x]
    exonEnds = [int(x) for x in exonEnds.split(',') if x]
    if cdsStart >= cdsEnd:
        # if no coding ...
        env.logger.debug('{} is non-coding.'.format(gene))
    #
    upstream = []
    intron = []
    exon = []
    coding = []
    downstream = []
    for s, e in zip(exonStarts, exonEnds):
        if exon:
            intron.append((chr, exon[-1][1] + 1, s - 1))
        exon.append((chr, s, e))
        if e >= cdsStart and s <= cdsEnd:
            cs = max(s, cdsStart)
            ce = min(e, cdsEnd)
            coding.append((chr, cs, ce))
        if cdsStart > s:
            upstream.append((chr, s, cdsStart - 1))
        if cdsEnd < e:
            downstream.append((chr, cdsEnd + 1, e))
    return {
        'name': gene,
        'strand': strand,
        'intron': intron,
        'exon': exon,
        'coding': coding,
        'upstream': upstream,
        'downstream': downstream,
        'build': proj.build
    }


def getRNASequence(structure, mutants=[], TtoU=True):
    '''Get protein sequence, mark locations of mutants if a list of variants
    are given (as a list of (chr, pos, ref, alt))'''
    ref = RefGenome(structure['build'])
    seq = ''
    for reg in structure['exon']:
        seq += ref.getSequence(reg[0], reg[1], reg[2])
    if mutants:
        loc_map = {}
        index = 0
        for reg in structure['exon']:
            for pos in range(reg[1], reg[2] + 1):
                loc_map[(reg[0], pos)] = index
                index += 1
        for m in mutants:
            if len(m) > 2:
                if m[2] == '-' or m[3] == '-' or len(m[2]) != 1 or len(
                        m[3]) != 1:
                    env.logger.warning(
                        'Get RNA sequence does not support indels yet.')
                    continue
            if m[0].startswith('chr'):
                loc = (m[0][3:], int(m[1]))
            else:
                loc = (m[0], int(m[1]))
            if loc in loc_map:
                idx = loc_map[loc]
                if len(m) > 2:
                    if m[2] != seq[idx]:
                        env.logger.warning(
                            'Reference alleles mismatch (chr {}, pos {}, ref {}, mutant ref {})'
                            .format(m[0], m[1], seq[idx], m[2]))
                    seq = seq[:idx] + m[3].lower() + seq[idx + 1:]
                else:
                    seq = seq[:idx] + seq[idx].lower() + seq[idx + 1:]
            #else:
            #    env.logger.debug('Failed to mark mutant {}'.format(loc))
    #
    # if len(seq) // 3 * 3 != len(seq):
    #    raise ValueError('Transcribed sequence should have length that is multiple of 3')
    if structure['strand'] == '-':
        # if len(seq) == 9, range(0, 9, 3) ==> 0, 3, 6
        seq = [complement_table[x] for x in reversed(seq)]
    TtoU_map = {
        'a': 'a',
        'A': 'A',
        'c': 'c',
        'C': 'C',
        'g': 'g',
        'G': 'G',
        't': 'u',
        'T': 'U'
    }
    if TtoU:
        return ''.join([TtoU_map[x] for x in seq])
    else:
        return ''.join(seq)


def getProteinSequence(structure, mutants=[]):
    '''Get protein sequence, mark locations of mutants if a list of variants
    are given (as a list of (chr, pos, ref, alt))'''
    ref = RefGenome(structure['build'])
    seq = ''
    for reg in structure['coding']:
        seq += ref.getSequence(reg[0], reg[1], reg[2])
    mut_idx = []
    if mutants:
        loc_map = {}
        index = 0
        for reg in structure['coding']:
            for pos in range(reg[1], reg[2] + 1):
                loc_map[(reg[0], pos)] = index
                index += 1
        for m in mutants:
            if len(m) > 2:
                if m[2] == '-' or m[3] == '-' or len(m[2]) != 1 or len(
                        m[3]) != 1:
                    env.logger.warning(
                        'Get protein sequence does not support indels yet.')
                    continue
            if m[0].startswith('chr'):
                loc = (m[0][3:], int(m[1]))
            else:
                loc = (m[0], int(m[1]))
            if loc in loc_map:
                idx = loc_map[loc]
                if len(m) > 2:
                    if m[2] != seq[idx]:
                        env.logger.warning(
                            'Reference alleles mismatch (chr {}, pos {}, ref {}, mutant ref {})'
                            .format(m[0], m[1], seq[idx], m[2]))
                    seq = seq[:idx] + m[3] + seq[idx + 1:]
                mut_idx.append(idx)
            #else:
            #    env.logger.debug('Failed to mark mutant {}'.format(loc))
    #
    if len(seq) // 3 * 3 != len(seq):
        raise ValueError(
            'Translated sequence should have length that is multiple of 3')
    if structure['strand'] == '+':
        # if len(seq) == 9, range(0, 9, 3) ==> 0, 3, 6
        pseq = [codon_table[seq[i:i + 3]] for i in range(0, len(seq), 3)]
    else:
        # if len(seq) == 9, range(6, -1, -3) ==> 6, 3, 0
        pseq = [
            codon_table_reverse_complement[seq[i:i + 3]]
            for i in range(len(seq) - 3, -1, -3)
        ]
    if mut_idx:
        if structure['strand'] == '+':
            for i in mut_idx:
                pseq[i // 3] = pseq[i // 3].lower()
        else:
            for i in mut_idx:
                pseq[-(1 + i // 3)] = pseq[-(1 + i // 3)].lower()
    return ''.join(pseq)


class _DeHTMLParser(HTMLParser):

    def __init__(self):
        HTMLParser.__init__(self)
        self.__text = []

    def handle_data(self, data):
        text = data.strip()
        if len(text) > 0:
            text = re.sub('[ \t\r\n]+', ' ', text)
            self.__text.append(text + ' ')

    def handle_starttag(self, tag, attrs):
        if tag == 'p':
            self.__text.append('\n\n\n\n')
        elif tag == 'br':
            self.__text.append('\n\n')
        elif tag == 'ul':
            self.__text.append('')
        elif tag == 'li':
            self.__text.append('\n\n  * ')

    def handle_endtag(self, tag):
        if tag == 'ul':
            self.__text.append('\n\n')
        if tag == 'li':
            self.__text.append('\n\n')

    def handle_startendtag(self, tag, attrs):
        if tag == 'br':
            self.__text.append('\n\n')

    def text(self):
        return ''.join(self.__text).strip()


def dehtml(text):
    try:
        parser = _DeHTMLParser()
        parser.feed(text)
        parser.close()
        return parser.text()
    except Exception as e:
        env.logger.warning('Failed to dehtml text: {}'.format(e))
        return text


def chunks_start_stop(length, rows=50):
    cycle = int(length / rows)
    startPos = 0
    endPos = 0
    for i in range(cycle + 1):
        startPos = endPos
        if i == cycle:
            endPos = length
        else:
            endPos = startPos + rows
        yield startPos, endPos


def chunks(data, rows=200):
    for i in range(0, len(data), rows):
        yield data[i:i + rows]


class RuntimeFiles:

    def __init__(self, output_files=[], pid=None):
        if not output_files:
            self.sig_file = None
            self.proc_out = None
            self.proc_err = None
            self.proc_lck = None
            self.proc_info = None
            self.proc_cmd = None
            self.proc_prog = None
            self.proc_done = None
            self.manifest = None
        else:
            if isinstance(output_files, list):
                output_file = output_files[0]
            elif isinstance(output_files, str):
                output_file = output_files
            else:
                raise ValueError('Invalid output file specification: {}'.format(
                    output_files))
            #
            # what is the relative
            # The parental directory of cache?
            cache_parent = os.path.dirname(env.cache_dir.rstrip(os.sep))
            #
            # is the file relative to this cache_parent?
            rel_path = os.path.relpath(
                os.path.realpath(os.path.expanduser(output_file)), cache_parent)
            # if this file is not relative to cache, use global signature file
            if rel_path.startswith('../'):
                self.sig_file = os.path.join(
                    env.local_resource, '.runtime',
                    os.path.realpath(os.path.expanduser(output_file)).lstrip(
                        os.sep))
            else:
                # if this file is relative to cache, use cache to store signature file
                self.sig_file = os.path.join(env.cache_dir, '.runtime',
                                             rel_path)
            # path to file
            sig_path = os.path.split(self.sig_file)[0]
            if not os.path.isdir(sig_path):
                try:
                    os.makedirs(sig_path)
                except Exception as e:
                    raise RuntimeError(
                        'Failed to create runtime directory {}: {}'.format(
                            sig_path, e))
            env.logger.trace('Using signature file {} for output {}'.format(
                self.sig_file, output_file))
            if pid is None:
                self.pid = os.getpid()
            else:
                self.pid = pid
            self.proc_out = '{}.out_{}'.format(self.sig_file, self.pid)
            self.proc_err = '{}.err_{}'.format(self.sig_file, self.pid)
            self.proc_lck = '{}.lck'.format(self.sig_file)
            self.proc_info = '{}.exe_info'.format(self.sig_file)
            self.proc_cmd = '{}.cmd'.format(self.sig_file)
            self.proc_done = '{}.done_{}'.format(self.sig_file, self.pid)
            self.proc_prog = '{}.working_{}'.format(self.sig_file, self.pid)
            self.manifest = '{}.manifest'.format(self.sig_file)
            #
            # now if there is an old signature file, let us move it to the new location
            if os.path.isfile('{}.exe_info'.format(output_file)):
                env.logger.info(
                    'Moving {}.exe_info from older version of variant tools to local cache'
                    .format(output_file))
                shutil.move('{}.exe_info'.format(output_file), self.proc_info)

    def clear(self, types=['out', 'err', 'done']):
        if self.sig_file is None:
            return
        for filename in sum(
            [glob.glob(self.sig_file + '.{}_*'.format(x)) for x in types], []):
            try:
                os.remove(filename)
            except Exception as e:
                env.logger.warning('Fail to remove {}: {}'.format(filename, e))
