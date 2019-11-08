#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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

import argparse
import io
import os
import re
import signal
import time
from collections import OrderedDict
from configparser import SafeConfigParser as RCP

from .project import Field
from .tester import ExternTest, freq
from .utils import (downloadFile, env, flatten, make_unique, mkdir_p, pairwise,
                    runCommand)

StringType = str
_mybytes = lambda s:bytes(s, 'utf8') #'ascii')
_mystr = lambda s:str(s, 'utf8')


def NoneStr(obj):
    return 'NA'

def BoolStr(obj):
    return obj and 'TRUE' or 'FALSE'

def ReprStr(obj):
    return repr(obj)

def LongStr(obj):
    rtn = repr(obj)
    if rtn[-1] == 'L': rtn = rtn[:-1]
    return rtn

def ComplexStr(obj):
    return repr(obj).replace('j', 'i')

def SeqStr(obj, method='c(', tail=')'):
    if not obj:
        return method + tail
    # detect types
    if isinstance(obj, set):
        obj = list(obj)
    obj_not_none = [x for x in obj if x is not None]
    if len(obj_not_none) == 0:
        # all is None
        return (method + ','.join(list(map(Str4R, obj))) + tail)
    obj0 = obj_not_none[0]
    tp0 = type(obj0)
    simple_types = [str, bool, int, int, float, complex]
    num_types = [int, int, float, complex]
    is_int = tp0 in (int, int)
    if tp0 not in simple_types and method == 'c(':
        method = 'list('
    else:
        # tps = isinstance(obj0, str) or num_types
        tps = isinstance(obj0, str) and [StringType] or num_types
        for i in obj_not_none[1:]:
            tp = type(i)
            if tp not in tps:
                if method == 'c(':
                    method = 'list('
                is_int = False
                break
            elif is_int and tp not in (int, int):
                is_int = False
    # convert

    return (is_int and 'as.integer(' or '') + \
      method + ','.join(list(map(Str4R, obj))) + tail + \
      (is_int and ')' or '')

def DictStr(obj):
    return 'list(' + ','.join(['%s=%s' % (Str4R(a[0]), Str4R(a[1])) for a in list(obj.items())]) + ')'

def OtherStr(obj):
    if hasattr(obj, '__iter__'): # for iterators
        if hasattr(obj, '__len__') and len(obj) <= 10000:
            return SeqStr(list(obj))
        else: # waiting for better solution for huge-size containers
            return SeqStr(list(obj))
    return repr(obj)

base_tps = [type(None), bool, int, int, float, complex, str, str, list, tuple, set, frozenset, dict] # use type(None) instead of NoneType since the latter cannot be found in the types module in Python 3
base_tps.reverse()
str_func = {type(None):NoneStr, bool:BoolStr, int:repr, float:repr, complex:ComplexStr,
        str:repr, list:SeqStr, tuple:SeqStr, set:SeqStr, frozenset:SeqStr, dict:DictStr}

def Str4R(obj, method = 'c'):
    '''
    convert a Python basic object into an R object in the form of string.
    '''
    if method != 'c':
        return SeqStr(obj, method=method + '(')
    if type(obj) in str_func:
        return str_func[type(obj)](obj)
    for tp in base_tps:
        if isinstance(obj, tp):
            return str_func[tp](obj)
    return OtherStr(obj)

class TimeoutException(Exception):
    pass

def timectrl(timeout_time, default):
    def timeout_function(f):
        def f2(*args):
            def timeout_handler(signum, frame):
                raise TimeoutException()

            old_handler = signal.signal(signal.SIGALRM, timeout_handler)
            signal.alarm(timeout_time) # triger alarm in timeout_time seconds
            try:
                retval = f(*args)
            except TimeoutException:
                return default
            finally:
                signal.signal(signal.SIGALRM, old_handler)
            signal.alarm(0)
            return retval
        return f2
    return timeout_function

class RConfig:
    def __init__(self, rlist):
        self.vars = OrderedDict()
        self.conf = RCP(allow_no_value=True)
        configstr = self.load(rlist)
        if configstr is None:
            return
        self.conf.readfp(io.StringIO(configstr))
        for section in self.conf.sections():
            try:
                self.parse(section)
            except Exception as e:
                raise ValueError("Invalid setting in configuration section [{0}]:\n\t{1}\n"
                                 "For help, please read "
                                 "http://varianttools.sf.net/Association/RTest".\
                                 format(section, e))
    def getvars(self):
        return self.vars

    def load(self, rlist):
        '''obtain raw configuration strings from R script'''
        res = None
        started = False
        ended = False
        for item in rlist:
            item = item.strip()
            if not item.startswith("#"):
                continue
            item = item[1:].strip()
            if item.upper().startswith("BEGINCONF"):
                started = True
                res = ''
            elif item.upper().startswith("ENDCONF"):
                ended = True
                break
            else:
                if res is not None:
                    res += item + '\n'
        if started and not ended:
            raise ValueError("Cannot find keyword '#ENDCONF'")
        if res is None or res == '':
            env.logger.warning("No output configuration specified in the R program. "
                               "As a result nothing will be written to standard output or result database.")
        return res

    def parse(self, section):
        '''self.vars[section] = {type: .., comment: .., name: .., [column_name:] ..}'''
        def typemapper(i):
            i = re.sub(r'"|\'', '', i)
            try:
                if i.lower() == 'integer':
                    return "INT"
                elif i.lower() == 'float':
                    return "FLOAT"
                elif i.lower() == 'string':
                    return "VARCHAR(255)"
                else:
                    raise
            except:
                raise ValueError("Unsupported type speficiation '{0}'. "
                                 "Should be one of 'integer', 'float' or 'string'".format(i))
        #
        if len(section.split()) > 1:
            raise ValueError ("White space is not allowed in section name '{0}'".\
                              format(section))
        self.vars[section] = {}
        options = self.conf.options(section)
        #
        for option in options:
            if option not in ['name','type','comment','n','columns','column_name']:
                raise ValueError("Unknown configuration parameter '{0}'".\
                                 format(option, section))
        #
        if "columns" in options and not "n" in options:
            raise ValueError("Please specify the length of output (n=?)".format(section))
        #
        if "n" in options:
            n = int(self.conf.get(section, "n"))
            self.vars[section]['n'] = n
            #
            if "name" in options:
                # 1. array length n
                names = [x.strip() for x in self.conf.get(section, "name").split(',')]
                if len(names) != n:
                    raise ValueError("length of 'name' does not equal specified "
                                     "{0}".format(n))
                else:
                    self.vars[section]['name'] = names
            else:
                self.vars[section]['name'] = ["{}{}".format(section, i+1) for i in range(n)]
            #
            if not "columns" in options or int(self.conf.get(section, "columns")) == 1:
                self.vars[section]['p'] = 1
                if "type" in options:
                    self.vars[section]['type'] = [typemapper(self.conf.get(section, 'type')) for i in range(n)]
                else:
                    self.vars[section]['type'] = ['FLOAT' for i in range(n)]
                if "column_name" in options:
                    raise ValueError ("Cannot set column_name without setting "
                                      "number of columns (columns=?)")
                if "comment" in options:
                    self.vars[section]['comment'] = \
                      [self.conf.get(section, "comment").rstrip('.') + ', {0}'.format(x) \
                       for x in self.vars[section]['name']]
                else:
                    self.vars[section]['comment'] = None
            else:
                # 2. is a n by p data frame or matrix
                p = int(self.conf.get(section, "columns"))
                self.vars[section]['p'] = p
                if "column_name" in options:
                    cnames = [x.strip() for x in self.conf.get(section, "column_name").split(',')]
                    if len(cnames) != p:
                        raise ValueError("length of 'column_name' does not equal specified "
                                         "{0}".format(p))
                    else:
                        self.vars[section]['column_name'] = cnames
                else:
                    self.vars[section]['column_name'] = [i+1 for i in range(p)]
                if "type" in options:
                    tp = [x.strip() for x in self.conf.get(section, "type").split(',')]
                    if len(tp) != p:
                        raise ValueError("length of 'type' does not equal specified "
                                         "{0}".format(p))
                    else:
                        self.vars[section]['type'] = \
                          flatten([[typemapper(x.strip()) for i in range(n)] for x in tp])
                else:
                    self.vars[section]['type'] = ['FLOAT' for i in range(n * p)]
                if "comment" in options:
                    self.vars[section]['comment'] = \
                      [self.conf.get(section, "comment").rstrip('.') + ', {0}_{1}'.format(x, y) \
                       for x, y in pairwise(self.vars[section]['name'], self.vars[section]['column_name'])]
                else:
                    self.vars[section]['comment'] = None
        else:
            self.vars[section]['n'] = 1
            self.vars[section]['p'] = 1
            # 3. is just a number
            if "type" in options:
                self.vars[section]['type'] = typemapper(self.conf.get(section, 'type'))
            else:
                self.vars[section]['type'] = 'FLOAT'
            if "name" in options:
                self.vars[section]['name'] = self.conf.get(section, "name")
            else:
                self.vars[section]['name'] = section
            if len(self.vars[section]['name'].split()) > 1:
                raise ValueError ("White space is not allowed in name string '{0}' for section '{1}'".\
                                  format(self.vars[section]['name'], section))
            if "comment" in options:
                self.vars[section]['comment'] = self.conf.get(section, "comment")
            else:
                self.vars[section]['comment'] = None
        # finally merge name and column name
        if "column_name" in self.vars[section]:
            self.vars[section]['name'] = \
               ['{}_{}'.format(i,j) for i,j \
                in pairwise(self.vars[section]['name'], self.vars[section]['column_name'])]
            del self.vars[section]['column_name']
        return

#
# R tests
# Wraps R packages via simple piping
#

class RTest(ExternTest):
    '''A general framework for association analysis using R programs'''
    def __init__(self, ncovariates, *method_args):
        '''Properly parse command arguments to parts of R script'''
        ExternTest.__init__(self, *method_args)
        self.ncovariates = ncovariates
        # step 1: load customer script
        self.source()
        # step 2: parse function input, combine with additional command method_args
        self.parseInput()
        # step 3: parse output, prepare output database field. field type have to be explicitly specified.
        self.parseOutput()
        # step 4: self.calculate
        # 4.1  format python data into S4 object string that matches the input: self.loadData()
        # 4.2 finish script for the output part: self.formatOutput
        # 4.3 pipe and run
        # 4.4 write output

    def source(self):
        '''find the R script and load it'''
        # _determine_algorithm() is empty for the base RTest
        basename = self._determine_algorithm()
        # No pre-defined R script available. Have to look for it
        if basename is None:
            # check filename
            basename, ext = os.path.splitext(self.script)
            if ext.upper() != '.R':
                raise ValueError("Improper R program filename '{0}' (should be '{1}.R')".\
                                 format(self.script, basename))
            # look for file
            found = False
            for path in [None, '.', os.path.join(env.shared_resource, 'programs'),
                         os.path.join(env.local_resource, 'programs'), env.cache_dir]:
                if path is None:
                    script = os.path.expanduser(self.script)
                else:
                    script = os.path.join(path, self.script)
                if os.path.exists(script):
                    found = True
                    basename = os.path.splitext(script)[0]
                    env.logger.info("Loading R script '{0}'".format(script))
                    with open(script, 'r') as f:
                        self.Rscript = [i.strip() for i in f.readlines()]
                    break
            if not found:
                raise ValueError("Cannot find R program '{0}'".format(self.script))
        self.basename = os.path.split(basename)[-1]

    def parseInput(self):
        '''input argument passed to the Rscript'''
        foostatement = ''
        start = None
        end = None
        for idx, item in enumerate(self.Rscript):
            # locate the main function
            if item.strip().startswith(self.basename):
                foostatement = item
                start = idx
            elif foostatement:
                foostatement += item
            if '{' in item and foostatement:
                end = idx
                break
        foo = re.match(r'{0}(\s*)(=|<-)(\s*)function(\s*)\((\s*)(?P<a>.+?)(\s*)\)(\s*){{'.\
                        format(self.basename), foostatement, re.DOTALL)
        if foo:
            if self.params is not None:
                updated_params = self._determine_params(foo.group('a'))
                self.Rscript = [i for j, i in enumerate(self.Rscript) if j < start or j > end]
                self.Rscript.insert(start, re.sub(r'function(\s*)\((.*?)\)(\s*){', r'function ({0}) {{'.\
                                                  format(updated_params), foostatement, re.DOTALL)
                    )
        else:
            raise ValueError("Cannot find main function '{0}' in the R program {1}".\
                             format(self.basename, self.script))
        return

    def parseOutput(self):
        '''analyze output and prepare database fields'''
        # remember to replace dot by underscore in field names
        self.outconf = RConfig(self.Rscript).getvars()
        if not self.outconf:
            self.fields.append(Field(name='name', index=None, type="VARCHAR(255)",
                                     adj=None, fmt=None, comment=''))
            self.outvars = []
            return
        #
        flds = []; tps = []; comments = []
        for var in self.outconf:
            for k, item in list(self.outconf[var].items()):
                if k == 'name':
                    fld = item
                elif k == 'type':
                    tp = item
                elif k == 'comment':
                    comment = item
                else:
                    continue
            if isinstance(fld, list):
                flds.extend(fld)
                tps.extend(tp)
                comments.extend(comment if comment else ['' for i in range(len(tps))])
            else:
                flds.append(fld)
                tps.append(tp)
                comments.append(comment if comment else '')
        #
        for fld, tp, comment in zip(flds, tps, comments):
            self.fields.append(Field(name='{0}'.format(fld.replace('.', '_')),
                                 index=None,
                                 type='{0}'.format(tp),
                                 adj=None,
                                 fmt=None,
                                 comment='{0}'.format(comment)))
        self.outvars = [(x, y) for x, y in zip(flds, tps)]
        return

    def loadData(self):
        '''data to S4 object string'''
        if not self.pydata:
            raise ValueError("Empty data set")
        self.Rdata = []
        self.Rdata.append('''setClass(Class="VATData",
        representation=representation(
        name="character",
        X="data.frame",
        Y="data.frame",
        V="data.frame",
        G="list"
        )
        )''')
        self.Rdata.append('''{0} = new(Class="VATData")'''.format(self.datvar))
        self.nvar = len(self.pydata['genotype'][0])
        self.nsample = len(self.pydata['genotype'])
        # handle duplicate coordinate names
        self.pydata['coordinate'] = make_unique(self.pydata['coordinate'])
        # Group name
        self.Rdata.append('''{0}@name = "{1}"'''.format(self.datvar, self.pydata['name']))
        # X matrix
        self.Rdata.append('''{0}@X = as.data.frame({1})'''.\
                            format(self.datvar, Str4R(self.pydata['genotype'], method = 'rbind')))
        self.Rdata.append('''colnames({0}@X) = {1}'''.\
                            format(self.datvar, Str4R([':'.join(x) for x in self.pydata['coordinate']])))
        self.Rdata.append('''rownames({0}@X) = {1}'''.\
                            format(self.datvar, Str4R(self.pydata['sample_name'])))
        # Y matrix
        self.Rdata.append('''{0}@Y = as.data.frame({1})'''.\
                            format(self.datvar, Str4R(self.pydata['covariates'] + [self.pydata['phenotype']], method = 'cbind')))
        self.Rdata.append('''colnames({0}@Y) = {1}'''.\
                            format(self.datvar, Str4R(self.pydata['covariate_name'] + self.pydata['phenotype_name'])))
        self.Rdata.append('''rownames({0}@Y) = {1}'''.\
                            format(self.datvar, Str4R(self.pydata['sample_name'])))

        # Variant information matrix
        if len(self.pydata['var_info']):
                self.Rdata.append('''{0}@V = as.data.frame({1})'''.\
                                    format(self.datvar, Str4R(list(zip(*self.pydata['var_info'])), method = 'cbind')))
                self.Rdata.append('''rownames({0}@V) = {1}'''.\
                                    format(self.datvar, Str4R([':'.join(x) for x in self.pydata['coordinate']])))
                self.Rdata.append('''colnames({0}@V) = {1}'''.\
                                    format(self.datvar, Str4R(self.pydata['var_info_header'])))

        # Genotype information matrix
        if len(self.pydata['geno_info']):
                for idx, item in enumerate(self.pydata['geno_info_header']):
                        self.Rdata.append('''{0}@G${2} = as.data.frame({1})'''.\
                        format(self.datvar, Str4R([[j[idx] for j in x] for x in self.pydata['geno_info']], method = 'rbind'), item))
                        self.Rdata.append('''colnames({0}@G${2}) = {1}'''.\
                            format(self.datvar, Str4R([':'.join(x) for x in self.pydata['coordinate']]), item))
                        self.Rdata.append('''rownames({0}@G${2}) = {1}'''.\
                            format(self.datvar, Str4R(self.pydata['sample_name']), item))
        # run calculation
        self.Rdata.append('''tryCatch({0}, error = function(err) q("no",1,FALSE))'''.format(self.basename))
        self.Rdata.append('''tryCatch({{result = {0}({1})}}, error = function(err) {{
        write(err$message, stderr()); q("no", 1, FALSE)}})'''.format(self.basename, self.datvar))
        return

    def formatOutput(self):
        '''some R code to convert the output to a plain string
            that matches exactly with the self.outvars'''
        # check if any return value needs to be collected
        script = '\n'
        if len(self.outvars) == 0:
            return script
        output = list(self.outconf.keys())
        for item in output:
            n = self.outconf[item]['n']
            p = self.outconf[item]['p']
            script += '''
            # check if exist
            if (! "{0}" %in% names(result)) {{
                write("[FATAL ERROR] Cannot find variable {0} in output", stderr())
                q("no", 1, FALSE)
            }}
            '''.format(item)
            if p == 1:
                script += '''
                # check if length matches for arrays
                if (is.null(result[["{0}"]])) {{
                    result[["{0}"]] = rep(NA, {1} )
                }} else {{
                    if (!is.null(dim(result[["{0}"]]))) result[["{0}"]] = as.vector(result[["{0}"]])
                    if (length(result[["{0}"]]) != {1}) {{
                        write("[WARNING] Length of output {0} does not match specified length {1}", stderr())
                        if (length(result[["{0}"]]) < {1}) result[["{0}"]][{1}] = NA
                        else result[["{0}"]] = result[["{0}"]][1:{1}]
                    }}
                }}
                '''.format(item, n)
            else:
                script += '''
                # check if dimension matches for data frames
                if (is.null(result[["{0}"]])) {{
                    result[["{0}"]] = rep(NA, {1} * {2})
                }} else {{
                    if (sum(dim(result[["{0}"]]) == c({1},{2})) != 2) {{
                        write("[FATAL ERROR] Dimension of output {0} does not match specified ({1}, {2})",
                         stderr())
                        q("no", 1, FALSE)
                    }}
                    # convert data frames / matrices to arrays
                    result[["{0}"]] = as.vector(t(result[["{0}"]]))
                }}
                '''.format(item, n, p)
        # write output
        script += 'write("BEGIN-VATOUTPUT", stdout());' + \
          'write(c({0}), stdout());'.format(', '.join(['result[["{0}"]]'.format(x) for x in output])) + \
          'write("END-VATOUTPUT", stdout())\n'
        return script

    def calculate(self, tmout):
        '''run command and return output results'''
        def forcefloat(i):
            try:
                i = float(i)
            except:
                i = float('nan')
            return i
        #
        def forceint(i):
            try:
                i = int(i)
            except:
                i = float('nan')
            return i
        #
        def typemapper(i):
            if i == "INT":
                return forceint
            elif i == 'FLOAT':
                return forcefloat
            else:
                return str
        #
        @timectrl(tmout, 'time expired')
        def timed_command(cmd, instream = None, msg = ''):
            return runCommand(cmd, instream, msg)
        #
        self.loadData()
        # write data and R script to log file for this group
        if self.data_cache_counter < self.data_cache or self.data_cache < 0:
            try:
                with open(os.path.join(env.cache_dir, '{0}.dat.R'.format(self.gname)), 'w') as f:
                    f.write('\n'.join(self.Rscript + self.Rdata))
                self.data_cache_counter += 1
            except:
                pass
        # print self.formatOutput()
        cmd = "R --slave --no-save --no-restore"
        cmd_logger = env.logger.debug
        na_value = float('nan')
        try:
            out = timed_command(cmd, "{0}".format('\n'.join(self.Rscript + self.Rdata) + self.formatOutput()),
                              "R message for {0}".format(self.pydata['name']))
            if out == 'time expired':
                cmd_logger = env.logger.warning
                na_value = -9
                raise ValueError("running time {}s has expired".format(tmout))
            else:
                out = out.split()
            if len(self.outvars):
                out = out[out.index("BEGIN-VATOUTPUT") + 1:out.index("END-VATOUTPUT")]
                res = [typemapper(x[1])(y) for x, y in zip(self.outvars, out)]
            else:
                res = [self.pydata['name']]
        except Exception as e:
            cmd_logger("Association test {} failed while processing '{}': {}".\
                              format(self.name, self.gname, e))
            res = [na_value]*len(self.outvars)
        return res

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''R test''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('script', type=str,
            help='R program to be loaded. The R script format has to follow the convention '
            'documented at http://varianttools.sf.net/Association/RTest')
        parser.add_argument('--name', default='RTest',
            help='''Name of the test that will be appended to names of output fields.''')
        parser.add_argument('--data_cache', '--data-cache', metavar='N', type=int, default = 0,
            help='''Name of R data sets to be written into cache folder, for debug purpose.''')
        # incorporate args to this class
        args, unknown = parser.parse_known_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        # handle unknown arguments
        self.params = self._determine_params(unknown, 'parse')
        self.data_cache_counter = 0

    def _determine_params(self, params, option = 'replace'):
        if option == 'parse':
            # format input argument list to input argument string dictionary,
            # e.g., ["--kernel", "linear"] --> {"kernel":"linear"}.
            rparams = {}
            current_key = None
            for p in params:
                if p.startswith('-'):
                    if not current_key:
                        current_key = re.sub(r"^\-+", '', p)
                    else:
                        raise ValueError("Invalid input '{0}' to Argument '{1}'".format(p, current_key))
                else:
                    rparams[current_key] = p
                    current_key = None
        else:
            self.datvar = None
            rparams = []
            defaults = []
            # input parameter is a string, from the original R function
            # e.g., kernel = 'linear', df = 2
            # have to update information of input string and from self.params
            params_list = [[j.strip() for j in x.split('=')] for x in self._param_split(params)]
            for idx, item in enumerate(params_list):
                if len(item) == 1 and idx == 0:
                    self.datvar = item[0]
                    env.logger.debug("Data set variable name is '{0}'".format(self.datvar))
                    rparams.append(self.datvar)
                elif len(item) == 2 and idx == 0:
                    raise ValueError("[R script error] Missing required argument: the input data set variable.")
                elif idx != 0:
                    if len(item) == 1:
                        if item[0] not in self.params:
                            raise ValueError("[R script error] Input value required for argument '{0}'".\
                                             format(item[0]))
                        else:
                            item.append('')
                    defaults.append(item[0])
                    if item[0] in self.params:
                        s = self.params[item[0]].strip()
                        # s should be a raw string, determined from the function,
                        # but was not inputted as raw
                        if re.match(r'(.*[^\s])\((.*)\)$', re.sub(r'\'|"', '', s).strip()):
                            # input is something like "c()", "rbind()" ...
                            # will remove quote marks around them
                            s = re.sub(r'^\'|\'$|^"|"$', '', s)
                        if re.match(r'^\'|\'$|^"|"$', item[1].strip()) \
                          and (not re.match(r'^\'|\'$|^"|"$', s)):
                          item[1] = '"{0}"'.format(s)
                        else:
                          item[1] = '{0}'.format(s)
                    rparams.append(' = '.join(item))
                else:
                    raise ValueError("[R script error] Invalid R argument '{0}'".format(' = '.join(item)))
            # check if input data variable is specified
            if self.datvar is None:
                raise ValueError("[R script error] Invalid parameter specification: "
                                 "no variable name for data object is found in '{0}'".format(params))
            # check if all self.params belong to default parameters
            for k in list(self.params.keys()):
                if k is None:
                    raise ValueError("Broken input argument {0}".format(self.params[k]))
                if k not in defaults:
                    raise ValueError("[R script error] Input parameter '{0}' is not defined in R function '{1}'".\
                                     format(k, self.basename))
            rparams = ', '.join(rparams)
        return rparams

    def _determine_algorithm(self):
        return None

    def _param_split(self, pm):
        # quotes match
        qmatch = True
        # () match
        pmatch = True
        res = []
        start = 0
        pm = list(pm)
        for idx, item in enumerate(pm):
            # handle "" ''
            if item in ['"', "'"]:
                qmatch = not qmatch
            # handle ()
            if item in ['(', ')']:
                pmatch = not pmatch
            # handle "()"
            if not pmatch and not qmatch:
                pmatch = True
            if item == ',' and pmatch and qmatch:
                # collect contents before this comma
                res.append(''.join(pm[start:idx]))
                start = idx + 1
        # finally collect the last one
        res.append(''.join(pm[start:]))
        return res

class SKAT(RTest):
    '''SKAT (Wu et al 2011) and SKAT-O (Lee et al 2012)'''
    def __init__(self, ncovariates, *method_args):
        self.skat_dir = self._check_skat()
        RTest.__init__(self, ncovariates, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''SNP-set (Sequence) Kernel Association Test (Wu et al 2011; Lee at all 2012).
            This test adapts the R package "SKAT" implemented & maintained by Dr. Seunggeun Lee, with a less
            complete interface and collection of features. Please refer to the SKAT documentation in R
            http://cran.r-project.org/web/packages/SKAT/
            for details of usage. Installation of SKAT R package v0.82 or higher is required to use this test.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('trait_type', type=str, choices = ['quantitative','disease'], default='quantitative',
            help='''Phenotype is quantitative trait or disease trait (0/1 coding).
            Default set to "quantitative"''')
        parser.add_argument('--name', default='SKAT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by SKAT
        #
        parser.add_argument('-k','--kernel', type=str, choices = ["linear", "linear.weighted", "IBS", "IBS.weighted", "quadratic", "2wayIX"],
            default="linear.weighted",
            help='''A type of kernel. Default set to "linear.weighted". Please refer to SKAT documentation for details.''')
        parser.add_argument('--beta_param', '--beta-param', nargs=2, type=float, default = [1,25],
            help='''Parameters for beta weights. It is only used with weighted kernels. Default set to (1,25).
                Please refer to SKAT documentation for details.''')
        parser.add_argument('--p_method', '--p-method', type=str, choices = ['davies','liu','liu.mod','optimal','optimal.adj'], default='davies',
            help='''A method to compute the p-value. The "optimal.adj" refers to the SKAT-O method. Default set to "davies". Please refer to SKAT documentation for details.''')
        parser.add_argument('-i','--impute', type=str, choices = ['fixed','random'], default='fixed',
            help='''A method to impute missing genotypes. Default set to "fixed". Please refer to SKAT documentation for details.''')
        parser.add_argument('--logistic_weights', '--logistic-weights', nargs=2, type=float, metavar='PARAM',
            help='''This option, if specified, will get the logistic weights from genotype matrix Z and apply this weight to SKAT. It requires
                two input parameters par1 and par2. To use the SKAT default setting, type `--logistic_weights 0.07 150'.
                Please refer to SKAT documentation for details.''')
        parser.add_argument('-r','--corr', nargs='*', type=float, default=[0],
            help='''The pho parameter of SKAT test. Default is 0. Please refer to SKAT documentation for details.''')
        parser.add_argument('--missing_cutoff', '--missing-cutoff', metavar='C', type=freq, default=0.15,
            help='''a cutoff of the missing rates of SNPs. Any SNPs with missing
                rates higher than cutoff will be excluded from the analysis. Default set to 0.15''')
        parser.add_argument('--resampling', metavar='N', type=int, default=0,
            help='''Number of resampling using bootstrap method. Set it to '0' if you do not want to apply resampling.''')
        parser.add_argument('--small_sample', '--small-sample', metavar='N', type=int, default=0,
            help='''Apply small sample model "SKAT_Null_Model_MomentAdjust" for binary traits having sample size smaller than N. Please refer to SKAT documentation for details.''')
        parser.add_argument('--resampling_kurtosis', '--resampling-kurtosis', metavar='N', type=int, default=0,
            help='''Number of resampling to estimate kurtosis, for small sample size adjustment.
            Set it to '0' if you do not wnat to apply the adjustment. The SKAT default setting is 10000. Please
            refer to SKAT documentation for details.''')
        parser.add_argument('--data_cache', '--data-cache', metavar='N', type=int, default = 0,
            help=argparse.SUPPRESS)
        # incorporate args to this class
        args = parser.parse_args(method_args)
        self.__dict__.update(vars(args))
        if self.trait_type == 'quantitative' and self.small_sample:
            env.logger.warning('Small sample-size adjustment does not apply to quantitative traits.')
            self.small_sample = False
        if self.resampling_kurtosis and not self.small_sample:
            env.logger.warning('Kurtosis adjustment is only effective with "--small_sample" option.')
            self.resampling_kurtosis = 0
        self.params = None
        self.data_cache_counter = 0

    def _determine_algorithm(self):
        '''Generate SKAT R script'''
        conf = '''
        # BEGINCONF
        #[sample.size]
        #type=integer
        #comment=Sample size of the test unit
        #[Q.stats]
        #comment=Q statistic for SKAT
        #[pvalue]
        #comment=p-value{0}
        {1}# ENDCONF
        '''.format(' (from resampled outcome)' if self.resampling else '',
                   '#[pvalue.noadj]\n#comment=The p-value of SKAT without small sample adjustment, '
                        'when small sample null model (SKAT_Null_Model_MomentAdjust) is applied for '
                        'association units having sample size smaller than {0}\n'.format(self.small_sample) if self.small_sample else '')
        self.Rscript = [x.strip() for x in conf.split('\n')]
        if self.skat_dir is None:
            self.Rscript.append('suppressMessages(library("SKAT"))')
        else:
            self.Rscript.append("suppressMessages(library('SKAT', lib.loc='{0}'))".format(self.skat_dir))
        self.datvar = 'dat'
        basename = 'SKAT{0}'.format(time.strftime('%b%d%H%M%S', time.gmtime()))
        # write header
        self.Rscript.append(basename + ' <- function(dat) {')
        self.Rscript.append('m <- ncol(dat@Y)')
        self.Rscript.append('Z <- as.matrix(dat@X)')
        # self.Rscript.append('Z[which(is.na(Z))] <- 9')
        # get parameters and residuals from the H0 model
        y_type = 'out_type="{0}", '.format('C' if self.trait_type == 'quantitative' else 'D')
        adjust_arg = 'Adjustment=TRUE'
        moment_adjust_arg = 'is_kurtosis_adj={0}, n.Resampling.kurtosis={1}'.\
            format('TRUE' if self.resampling_kurtosis else 'FALSE', self.resampling_kurtosis)
        # small sample null model
        self.Rscript.append('small.sample <- F')
        h_null = 'obj <- SKAT_Null_Model(dat@Y[, m]~{0}, {1}n.Resampling={2}, type.Resampling="bootstrap",{3})'.\
                format(1 if not self.ncovariates else 'as.matrix(dat@Y[, -m])',
                       y_type, self.resampling, adjust_arg)
        if self.small_sample:
            h_null = '''
            if (nrow(dat@Y) < {0}) {{
            small.sample <- T
            obj <- SKAT_Null_Model_MomentAdjust(dat@Y[, m]~{1}, n.Resampling={2},
                 type.Resampling="bootstrap", {3})
            }} else {{
            {4}
            }}
            '''.format(self.small_sample, 1 if not self.ncovariates else 'as.matrix(dat@Y[, -m])',
                       self.resampling, moment_adjust_arg, h_null)
        self.Rscript.extend([x.strip() for x in h_null.split('\n')])
        # get logistic weights
        if self.logistic_weights:
            self.Rscript.append('weights <- Get_Logistic_Weights(Z, par1={0}, par2={1})'.\
                     format(self.logistic_weights[0], self.logistic_weights[1]))
        # SKAT test
        skat_test = '''re <- SKAT(Z, obj, kernel="{0}", method="{1}", weights.beta=c({2}, {3}), {4}
                impute.method="{5}", r.corr={6}, is_check_genotype={7}, is_dosage={8}, missing_cutoff={9})'''.\
                        format(self.kernel, self.p_method, self.beta_param[0],
                               self.beta_param[1], 'weights=weights, ' if self.logistic_weights else '',
                               self.impute,
                               self.corr[0] if len(self.corr) == 1 else 'c('+','.join(map(str, self.corr))+')',
                               'FALSE' if self.logistic_weights else 'TRUE',
                               'TRUE' if self.logistic_weights else 'FALSE',
                               self.missing_cutoff)
        self.Rscript.append(skat_test)
        # get results
        self.Rscript.extend(['p <- re[["p.value"]]', 'stat <- ifelse (is.null(dim(re$Q)), NA, re$Q[1,1])'])
        if self.resampling:
            self.Rscript.append('p <- Get_Resampling_Pvalue(re)[["p.value"]]')
        # return
        self.Rscript.append('return(list(sample.size = nrow(dat@Y), '
                 'Q.stats = stat, '
                 '{0} '
                 'pvalue = p))'.\
                 format('pvalue.noadj = ifelse(small.sample,re[["p.value.noadj"]],NA),'))
        self.Rscript.append('}')
        # write
        self.script = os.path.join(env.cache_dir, basename + ".R")
        with open(self.script, "w") as f:
                  f.write('\n'.join(self.Rscript))
        return os.path.splitext(self.script)[0]

    def _install_skat(self, skat_dir):
        try:
            version = '0.82'
            skat_url = 'programs/SKAT_{0}.tar.gz'.format(version)
            env.logger.info('Downloading SKAT (V{0}) ...'.format(version))
            skat_tar = downloadFile(skat_url, env.temp_dir)
            mkdir_p(skat_dir)
            env.logger.info('Installing SKAT (V{0}) to {1} ...'.format(version, skat_dir))
            runCommand(['R', '-e', 'install.packages("{0}", repos = NULL, type = "source", lib="{1}")'.format(skat_tar, skat_dir)])
            runCommand(["R", "-e", "library('SKAT', lib.loc='{0}')".format(skat_dir)])
        except Exception as e:
            raise ValueError("Cannot auto-install / load R library SKAT: {0}".format(e))
        return skat_dir

    def _check_skat(self):
        # Check for R/SKAT installation
        skat_dir = 'NULL'
        for item in [None, os.path.join(env.shared_resource, 'Rlib'),
                         os.path.join(env.local_resource, 'Rlib')]:
            try:
                library = "'SKAT'{}".format('' if item is None else ", lib.loc='{0}'".format(item))
                runCommand(["R", "-e", "library({})".format(library)])
                skat_dir = item
                break
            except:
                continue
        if skat_dir == 'NULL':
            skat_dir = self._install_skat(os.path.join(env.local_resource, 'Rlib'))
        return skat_dir
