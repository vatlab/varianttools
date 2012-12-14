#!/usr/bin/env python
#
# $File: importer.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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
import re
import array
import time
from heapq import heappush, heappop, heappushpop
from multiprocessing import Process, Pipe, Value, Lock, Manager
if sys.version_info.major == 2:
    from itertools import izip, repeat
else:
    izip = zip
    from itertools import repeat
from collections import defaultdict
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, \
    normalizeVariant, openFile, DatabaseEngine, hasCommand, consolidateFieldName, \
    downloadFile, runOptions

# preprocessors
from .preprocessor import *

#
#
# Functors to process input 
#
#
# Extractors to extract value from a field
class ExtractField:
    def __init__(self, index, sep=';', default=None):
        '''Define an extractor that returns the index-th (1-based) field of the fields
        separated by specified delimiter. Return default if unsuccessful.'''
        self.index = index - 1
        self.sep = sep
        self.default = default
    
    def __call__(self, item):
        try:
            return item.split(self.sep, self.index + 1)[self.index]
        except:
            return self.default

class CheckSplit:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. Return default if unsuccessful. It differs from
        SplitField in that it will return the item itself (instead of a tuple
        of one element) when there is only one element. The item will then
        will copy if multiple items exist.'''
        self.sep = sep
    
    def __call__(self, item):
        return item if self.sep not in item else tuple(item.split(self.sep))
    
class SplitField:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        return tuple(item.split(self.sep))

class ExtractFlag:
    def __init__(self, name, sep=';'):
        '''Define an extractor that returns 1 is item contains name as one of the fields,
        and 0 otherwise. No space is allowed between delimiter and flag.'''
        self.n = name
        self.s = name + sep
        self.e = sep + name
        self.m = sep + name + sep
    
    def __call__(self, item):
        # this should be faster than
        #
        #     if self.name in item.split(self.sep):
        # 
        # because we do not have to split the whole string.
        #
        if self.n not in item:
            return '0'
        # put the most common case first
        if self.m in item or item.startswith(self.s) or item.endswith(self.e) or item == self.n:
            return '1'
        else:
            return '0'

class CommonLeading:
    '''Find the common leading piece of two input strings (ref and alt).
    '''
    def __init__(self):
        pass

    def _commonLeading(self, ref, alt):
        common_leading = 0
        for i in range(min(len(ref), len(alt))):
            if ref[i] == alt[i]:
                common_leading += 1
        return ref[:common_leading]

    def __call__(self, item):
        if ',' in item[1]:
            return tuple([self._commonLeading(item[0], alt) for alt in item[1].split(',')])
        else:
            return self._commonLeading(item[0], item[1])

class CommonEnding:
    '''Find the common ending piece of two input strings (ref and alt).
    '''
    def __init__(self):
        pass
    
    def _commonEnding(self, ref, alt):
        common_leading = 0
        for i in range(min(len(ref), len(alt))):
            if ref[i] == alt[i]:
                common_leading += 1
        if common_leading > 0:
            ref = ref[common_leading:]
            alt = alt[common_leading:]
        common_ending = 0
        for i in range(-1, - min(len(ref), len(alt)) - 1, -1):
            if ref[i] == alt[i]:
                common_ending -= 1
            else:
                break
        if common_ending < 0:
            return ref[common_ending:]
        else:
            return ''
    
    def __call__(self, item):
        if ',' in item[1]:
            return tuple([self._commonEnding(item[0], alt) for alt in item[1].split(',')])
        else:
            return self._commonEnding(item[0], item[1])


class __FieldFromFormat:
    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that return the value of a field according 
        to a format string. This is used to extract stuff from the format
        string of vcf files.
        '''
        self.name = name
        self.sep = sep
        self.default = default
        self.factory = defaultdict(dict)
        self.index = {}

    def __call__(self, item):
        try:
            # first try to get from a global factory
            return self.factory[item[0]][item[1]]
        except:
            fmt, val = item
            try:
                # now split .... assuming the format has been handled before.
                # this should be the case most of the time
                res = val.split(self.sep)[self.index[fmt]]
                # we assume that the most common ones has been added...
                # and we do not want to add all sorts of rare values forever
                if len(self.factory[fmt]) < 10000:
                    self.factory[fmt][val] = res
                return res
            except:
                # if the format has not been handled before.
                if fmt not in self.index:
                    fields = fmt.split(self.sep)
                    if self.name in fields:
                        self.index[fmt] = fields.index(self.name)
                    else:
                        self.index[fmt] = None
                # try again
                try:
                    res = val.split(self.sep)[self.index[fmt]]
                    if len(self.factory[fmt]) < 10000:
                        self.factory[fmt][val] = res
                    return res
                # if still error
                except:
                    self.factory[fmt][val] = self.default
                    return self.default

__all_field_from_format = {}

def FieldFromFormat(name, sep=';', default=None):
    # this is a factory of __FieldFromFormat class
    #
    global __all_field_from_format
    if (name, sep, default) in __all_field_from_format:
        return __all_field_from_format[(name, sep, default)]
    else:
        obj = __FieldFromFormat(name, sep, default)
        __all_field_from_format[(name, sep, default)] = obj
        return obj

class VcfGenotype:
    def __init__(self, default=None):
        '''Define an extractor that extract genotype from a .vcf file'''
        self.default = default
        self.map = {'0/0': default, '0|0': default,
            '0/1': ('1',), '1/0': ('1',), '0|1': ('1',), '1|0': ('1',),
            '1/1': ('2',), '1|1': ('2',),
            '0/2': ('0', '1'), '2/0': ('0', '1'), '0|2': ('0', '1'), '2|0': ('0', '1'), 
            '1/2': ('-1', '-1'), '2/1': ('-1', '-1'), '1|2': ('-1', '-1'), '2|1': ('-1', '-1'),
            '2/2': ('0', '2'), '2|2': ('0', '2'),
            '0': default, '1': ('1',)}

    def __call__(self, item):
        # the most common and correct case...
        try:
            return self.map[item.partition(':')[0]]
        except KeyError:
            return None

class VcfGenoFromFormat:
    def __init__(self, default=None):
        '''Define an extractor that return genotype according to a format string.
        This is used to extract genotype from the format string of vcf files.
        '''
        self.fmt = '\t'
        self.idx = None
        self.default = default
        self.map = {'0/0': default, '0|0': default,
            '0/1': ('1',), '1/0': ('1',), '0|1': ('1',), '1|0': ('1',),
            '1/1': ('2',), '1|1': ('2',),
            '0/2': ('0', '1'), '2/0': ('0', '1'), '0|2': ('0', '1'), '2|0': ('0', '1'), 
            '1/2': ('-1', '-1'), '2/1': ('-1', '-1'), '1|2': ('-1', '-1'), '2|1': ('-1', '-1'),
            '2/2': ('0', '2'), '2|2': ('0', '2'),
            '0': default, '1': ('1',)}

    def __call__(self, item):
        # the most common and correct case...
        try:
            if item[0][:2] == 'GT':
                return self.map[item[1].partition(':')[0]]
            elif item[0] != self.fmt:
                fmt, val = item
                self.fmt = fmt
                fields = fmt.split(':')
                if 'GT' in fields:
                    self.idx = fields.index('GT')
                    return self.map[val.split(':')[self.idx]]
                else:
                    self.idx = None
                    return self.default
            return self.map[item[1].split(':', self.idx + 1)[self.idx]] if self.idx is not None else self.default
        except KeyError:
            return None
        
class ExtractValue:
    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that returns the value after name in one of the fields,
        and a default value if no such field is found. No space is allowed between 
        delimiter and the name.'''
        self.name = name
        self.sep = sep
        #self.pos = len(name)
        self.default = default

    def __call__(self, item):
        if self.name not in item:
            return self.default
        #
        # Using two partisions seems to be a tiny bit faster than 
        # split and startswith
        #
        #for field in item.split(self.sep):
        #    if field.startswith(self.name):
        #        return field[self.pos:]
        ss = item.partition(self.name)
        return ss[2].partition(self.sep)[0] if ss[2] is not None else self.default

class IncreaseBy:
    '''Increase passed value by a given number, will convert input to integer'''
    def __init__(self, inc=1):
        '''Adjust position'''
        self.inc = inc

    def __call__(self, item):
        return str(int(item) + self.inc) if item.isdigit() else None

class MapValue:
    '''Map value to another one, return default if unmapped'''
    def __init__(self, map, default=None):
        self.map = map
        self.default = default

    def __call__(self, item):
        try:
            return self.map[item]
        except:
            return self.default
        
class RemoveLeading:
    '''Remove specified leading string if the input string starts with it. Used
    for example to remove chr from inputs such as chr15'''
    def __init__(self, val):
        self.val = val
        self.vlen = len(val)

    def __call__(self, item):
        return item[self.vlen:] if item.startswith(self.val) else item

class EncodeGenotype:
    '''Encode 1/1, 1/2 etc to variant tools code'''
    def __init__(self, default=None):
        self.map = {'0/0': default, '0|0': default,
            '0/1': ('1',), '1/0': ('1',), '0|1': ('1',), '1|0': ('1',),
            '1/1': ('2',), '1|1': ('2',),
            '0/2': ('0', '1'), '2/0': ('0', '1'), '0|2': ('0', '1'), '2|0': ('0', '1'), 
            '1/2': ('-1', '-1'), '2/1': ('-1', '-1'), '1|2': ('-1', '-1'), '2|1': ('-1', '-1'),
            '2/2': ('0', '2'), '2|2': ('0', '2'),
            '0': default, '1': ('1',)}

    def __call__(self, item):
        return self.map[item]
        
class Nullify:
    '''Change specified input value to NULL '''
    def __init__(self, val):
        self.val = val
        if type(self.val) == str:
            self.__call__ = self.nullify_single
        else:
            self.__call__ = self.nullify_multiple

    def nullify_single(self, item):
        return None if item == self.val else item

    def nullify_multiple(self, item):
        return None if item in self.val else item

class IgnoredRecord(Exception):
    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        return repr(self.value)

class DiscardRecord:
    def __init__(self, val):
        self.val = val
        if type(self.val) == str:
            self.__call__ = self.discard_single
        else:
            self.__call__ = self.discard_multiple

    def discard_single(self, item):
        if item == self.val:
            raise IgnoredRecord()
        return item

    def discard_multiple(self, item):
        if item in self.val:
            raise IgnoredRecord()
        return item
    
__databases = {}

class _DatabaseQuerier:
    '''This query a field from an annotation database'''
    def __init__(self, cursor, name, res_field, cond_fields, default=None):
        '''Supose res_field is alt, cond_fields are chr,pos, this querier
        will get alt using query
          SELECT dbSNP.alt FROM dbSNP WHERE chr=VAL1 AND pos=VAL2
        '''
        self.default = default
        self.cur = cursor
        self.query = 'SELECT {} FROM {} WHERE {}'.format(res_field,
            name, ' AND '.join(['{}=?'.format(x) for x in cond_fields]))

    def __call__(self, item):
        self.cur.execute(self.query, item)
        res = self.cur.fetchall()
        if len(res) == 1:
            return res[0][0]
        elif len(res) > 1:
            return tuple([x[0] for x in res])
        else:
            return self.default

def DatabaseQuerier(dbfile, res_field, cond_fields, default=None):
    global __databases
    if dbfile not in __databases:
        db = DatabaseEngine()
        if not os.path.isfile(dbfile):
            raise ValueError('Database file {} does not exist'.format(dbfile))
        db.connect(dbfile)
        cur = db.cursor()
        tables = db.tables()
        try:
            name = [x for x in tables if x.endswith('_info')][0][:-5]
        except Exception as e:
            raise ValueError('Incorrect database (missing info table): {}'.format(e))
        if not name in tables:
            raise ValueError('Incorrect database (missing table {})'.format(name))
        if not name + '_field':
            raise ValueError('Incorrect database (missing field table)')
        for fld in cond_fields.split(','):
            if not db.hasIndex('{}_idx'.format(fld)):
                cur.execute('CREATE INDEX {0}_idx ON {1} ({0} ASC);'.format(fld, name))
        __databases[dbfile] = (cur, name)
    return _DatabaseQuerier(__databases[dbfile][0], __databases[dbfile][1], 
        res_field, cond_fields.split(','), default)

class SequenceExtractor:
    '''This sequence extractor extract subsequence from a pre-specified sequence'''
    def __init__(self, filename):
        if not os.path.isfile(filename):
            filename = downloadFile(filename)
        if not os.path.isfile(filename):
            raise valueError('Failed to obtain sequence file {}'.format(filename))
        # a dictionary for seq for each chromosome
        self.seq = {}
        # we assume that the input file has format
        #
        # >chr9
        # seq
        # >chr10
        # seq
        #
        current_chr = None
        # openFile can open .gz file directly
        cnt = lineCount(filename)
        prog = ProgressBar('Reading ref genome sequences', cnt)
        with openFile(filename) as input:   
            # for each chromosome? need to fix it here
            for idx, line in enumerate(input):
                line = line.decode()
                if line.startswith('>'):
                    chr = line[1:].split()[0]
                    if chr.startswith('chr'):
                        chr = chr[3:] 
                    self.seq[chr] = array.array('b', [])
                else:
                    self.seq[chr].fromstring(line.rstrip())
                if idx % 10000 == 0:
                    prog.update(idx)
        # use another key with 'chr' to point to the same item so that the dictionary 
        # works with both ['2'] and ['chr2']
        for key in self.seq:
            self.seq['chr' + key] = self.seq[key]
        prog.done()

    def __call__(self, item):
        return self.seq[item[0]][item[1]]

# this is a dictionary to save extractors for each file used
g_SeqExtractor = {}
def SeqAtLoc(filename):
    # return the same object for multiple instances of SeqAtLoc because
    # we do not want to read the fasta file multiple times
    if filename not in g_SeqExtractor:
        g_SeqExtractor[filename] = SequenceExtractor(filename)
    return g_SeqExtractor[filename]
    
class SequentialExtractor:
    def __init__(self, extractors):
        '''Define an extractor that calls a list of extractors. The string extracted from
        the first extractor will be passed to the second, and so on.'''
        self.extractors = []
        for e in extractors:
            if hasattr(e, '__call__'):
                self.extractors.append(e.__call__)
            else:
                self.extractors.append(e)

    def __call__(self, item):
        for e in self.extractors:
            # if multiple records are returned, apply to each of them
            if type(item) is tuple:
                if type(item[0]) is tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = tuple(e(x) for x in item)
            # if item is None or ''
            elif not item:
                return item
            else:
                item = e(item)
        return item

#
#
# Process each line using the above functors
#
#
class LineProcessor:
    '''An interpreter that read a record (a line), process it and return processed records.'''
    def __init__(self, fields, build, delimiter, ranges, logger):
        '''
        fields: a list of fields with index, adj (other items are not used)
        builds: index(es) of position, reference allele and alternative alleles. If 
            positions are available, UCSC bins are prepended to the records. If reference
            and alternative alleles are available, the records are processed for correct
            format of ref and alt alleles.
        delimiter: how to split line
        ranges: range of fields (var, var_info, GT, GT_info), used to determine
            var_info and GT fields when subsets of samples are imported.
        '''
        self.logger = logger
        self.build = build
        self.raw_fields = fields
        self.fields = []
        self.delimiter = delimiter
        # column range contains the range [start, end) of output for each
        # raw field. If it is None, no output is available. This field tells
        # the user how to split and handle output fields
        self.columnRange = [None] * len(self.raw_fields)
        self.ranges = ranges
        #
        self.first_time = True
        self.nColumns = 0     # number of columns 
        self.import_var_info = True
        self.import_sample_range = None  # genotype fields might be disabled
        self.maxInputCol = 0      # sometimes processor do not have to split input all the way through

    def reset(self, import_var_info = True, import_sample_range = None):
        '''Reset the processor for its internal state, which is used when the processor
        is used to process a new batch of data (with different format), or with
        modified behavior.
        
        import_var_info: if set to False, variant info will not be imported.
        import_sample_range: if set to a range, only selected samples are handled
        '''
        self.first_time = True
        self.fields = []
        self.nColumns = 0
        #
        self.import_var_info = import_var_info
        self.import_sample_range = import_sample_range

    def process(self, tokens):
        if not self.first_time:
            # if not first time, we know the maximum fields we need, if
            # self.maxInputCol is n, we need to have n+1 fields to guarantee that
            # the piece n is a properly split one
            if type(tokens) is not list:
                tokens = [x.strip() for x in tokens.split(self.delimiter, self.maxInputCol)][:self.maxInputCol]
        else:
            if type(tokens) is not list:
                tokens = [x.strip() for x in tokens.split(self.delimiter)]
            self.nColumns = len(tokens)
            cIdx = 0             # column index
            num_sample = -1      # number of samples ...
            for fIdx, field in enumerate(self.raw_fields):
                if not self.import_var_info:
                    # if sample range is not None, we do not import variant information either
                    if fIdx >= self.ranges[1] and fIdx < self.ranges[2]:
                        continue
                # if do not import any sample
                if self.import_sample_range is not None and \
                    self.import_sample_range[0] == self.import_sample_range[1] and \
                    fIdx >= self.ranges[2]:
                    continue
                try:
                    # get an instance of an extractor, or a function
                    e = eval(field.adj) if field.adj else None
                    # 1. Not all passed object has __call__ (user can define a lambda function)
                    # 2. Althoug obj(arg) is equivalent to obj.__call__(arg), saving obj.__call__ to 
                    #    e will improve performance because __call__ does not have to be looked up each time.
                    # 3. Passing object directly has an unexpected side effect on performance because comparing
                    #    obj to 1 and 'c' later are very slow because python will look for __cmp__ of the object.
                    if hasattr(e, '__iter__'):
                        # if there are multiple functors, use a sequential extractor to handle them
                        e = SequentialExtractor(e)
                    if hasattr(e, '__call__'):
                        e = e.__call__
                    indexes = []
                    # indexes gives the real indexes, for example
                    # 8,8::2 might yield
                    #
                    # indexes = [8, [8,10,12,14]]
                    #
                    for x in field.index.split(','):
                        if ':' in x:
                            # a slice
                            if x.count(':') == 1:
                                start,end = map(str.strip, x.split(':'))
                                step = None
                            else:
                                start,end,step = map(str.strip, x.split(':'))
                                step = int(step) if step else None
                            start = int(start) - 1 if start else None
                            if end.strip():
                                if int(end) >= 0:   # position index, shift by 1
                                    end = int(end) - 1
                                else:               # relative to the back, do not move
                                    end = int(end)
                            else:
                                end = None
                            indexes.append(slice(start, end, step))
                        else:
                            # easy, an integer
                            indexes.append(int(x) - 1)
                    #
                    if ':' not in field.index:
                        # case of 'index=10'
                        if len(indexes) == 1:
                            # int, True means 'not a tuple'
                            self.fields.append((indexes[0], True, e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            self.maxInputCol = max(self.maxInputCol, indexes[0])
                            cIdx += 1
                        # case of index=7,8,9
                        else:
                            # a tuple
                            self.fields.append((tuple(indexes), False, e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            self.maxInputCol = max(self.maxInputCol, max(indexes))
                            cIdx += 1
                    # if there is only one slice
                    # case of index=8::2
                    elif len(indexes) == 1:
                        # single slice
                        cols = range(len(tokens))[indexes[0]]
                        if self.import_sample_range is not None:
                            # limiting the columns to import
                            if self.import_sample_range[0] >= len(cols) or self.import_sample_range[1] > len(cols):
                                raise ValueError('ERROR PROCESSING subset of samples.')
                            cols = cols[self.import_sample_range[0]:self.import_sample_range[1]]
                        for c in cols:
                            self.fields.append((c, True, e))
                            self.maxInputCol = max(self.maxInputCol, c)
                        if num_sample == -1:
                            num_sample = len(cols)
                        elif num_sample != len(cols):
                            sys.exit('The first line of input has inconsistent number of fields for samples, perhaps due to incorrect use of delimiters.')
                        self.columnRange[fIdx] = (cIdx, cIdx + len(cols))
                        cIdx += len(cols)
                    else:
                        # we need to worry about mixing integer and slice
                        expanded_indexes = [repeat(s, len(tokens)) if type(s) == int else range(len(tokens))[s] for s in indexes]
                        count = 0
                        for idx, c in enumerate(izip(*expanded_indexes)):
                            if self.import_sample_range is not None:
                                if idx < self.import_sample_range[0] or idx >= self.import_sample_range[1]:
                                    continue
                            count += 1
                            self.fields.append((tuple(c), False, e))
                            self.maxInputCol = max(self.maxInputCol, max(c))
                        if num_sample == -1:
                            num_sample = count
                        elif num_sample != count:
                            sys.exit('The first line of input has inconsistent number of fields for samples, perhaps due to incorrect use of delimiters.')
                        self.columnRange[fIdx] = (cIdx, cIdx + count)
                        cIdx += count
                except Exception as e:
                    sys.exit('Incorrect value adjustment functor or function {}: {}'.format(field.adj, e))
                    raise ValueError('Incorrect value adjustment functor or function {}: {}'.format(field.adj, e))
            self.first_time = False
            self.maxInputCol += 1    # use 1-indexes maxInputCol
        #
        try:
            # we first trust that nothing can go wrong and use a quicker method
            records = [(tokens[col] if t else [tokens[x] for x in col]) if adj is None else \
                (adj(tokens[col]) if t else adj([tokens[x] for x in col])) for col,t,adj in self.fields]
        except IgnoredRecord as e:
            return
        except Exception:
            # If anything wrong happends, process one by one to get a more proper error message (and None values)
            records = []
            for col, t, adj in self.fields:
                try:
                    item = tokens[col] if t else [tokens[x] for x in col]
                except IndexError:
                    raise ValueError('Cannot get column {} of the input line, which has only {} columns (others have {} columns).'.format(\
                        col + 1 if type(col) is int else [x+1 for x in col], len(tokens), self.nColumns))
                if adj is not None:
                    try:
                        item = adj(item)
                    except Exception as e:
                        self.logger.debug('Failed to process field {}: {}'.format(item, e))
                        # missing ....
                        item = None
                records.append(item)
        #
        num_records = max([len(item) if type(item) is tuple else 1 for item in records]) if records else 1
        # handle records
        if not self.build:
            # there is no build information, this is 'field' annotation, nothing to worry about
            if num_records == 1:
                yield [], [x[0] if type(x) is tuple else x for x in records]
            else:
                for i in range(num_records):
                    yield [], [(x[i] if i < len(x) else None) if type(x) is tuple else x for x in records]
        elif len(self.build[0]) == 1:
            for i in range(num_records):
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) is tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None) if type(x) is tuple else x for x in records]
                bins = [getMaxUcscBin(int(rec[pos_idx]) - 1, int(rec[pos_idx])) if rec[pos_idx] else None for pos_idx, in self.build]
                yield bins, rec
        else:
            for i in range(num_records):
                bins = []
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) is tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None) if type(x) is tuple else x for x in records]
                for pos_idx, ref_idx, alt_idx in self.build:
                    bin, pos, ref, alt = normalizeVariant(int(rec[pos_idx]) if rec[pos_idx] else None, rec[ref_idx], rec[alt_idx])
                    bins.append(bin)
                    rec[pos_idx] = pos
                    rec[ref_idx] = ref
                    rec[alt_idx] = alt
                yield bins, rec


#
# Read record from disk file
#

def TextReader(processor, input, varIdx, getNew, jobs, encoding, logger):
    '''
    input: input file
    varIdx: variant index, if specified, only matching variants will be returned
        used by updater, or only new variants will be returned by an variant reader
    getNew: if getNew is true, only return variants that are NOT in varIdx. Otherwise,
        return variants that are in varIdx.
    jobs: number of jobs
    encoding: file encoding
    '''
    if jobs == 0:
        return EmbeddedTextReader(processor, input, varIdx, getNew, encoding, logger)
    elif jobs == 1:
        return StandaloneTextReader(processor, input, varIdx, getNew, encoding, logger)
    else:
        return MultiTextReader(processor, input, varIdx, getNew, jobs, encoding, logger)

class EmbeddedTextReader:
    #
    # This reader read the file from the main process. No separate process is spawned.
    #
    def __init__(self, processor, input, varIdx, getNew, encoding, logger):
        self.num_records = 0
        self.skipped_lines = 0
        self.processor = processor
        self.input = input
        self.varIdx = varIdx
        self.getNew = getNew
        self.encoding = encoding
        self.logger = logger

    def records(self): 
        first = True
        line_no = 0
        with openFile(self.input) as input_file:
            for line in input_file:
                line_no += 1
                line = line.decode(self.encoding)
                try:
                    if line.startswith('#'):
                        continue
                    for bins,rec in self.processor.process(line):
                        if first:
                            self.columnRange = self.processor.columnRange
                            first = False
                        self.num_records += 1
                        if self.varIdx is not None:
                            var_key = (rec[0], rec[2], rec[3])
                            if self.getNew:
                                # only need new variant, so continue if variant in varIdx.
                                if var_key in self.varIdx and rec[1] in self.varIdx[var_key]:
                                    continue
                            else:
                                # only need existing variant, continue if variant not in varIdx
                                if var_key not in self.varIdx or rec[1] not in self.varIdx[var_key]:
                                    continue
                        yield (line_no, bins, rec)
                except Exception as e:
                    self.logger.debug('Failed to process line "{}...": {}'.format(line[:20].strip(), e))
                    self.skipped_lines += 1


class ReaderWorker(Process):
    '''
    This class starts a process and use passed LineProcessor
    to process input line. If multiple works are started,
    they read lines while skipping lines (e.g. 1, 3, 5, 7, ...)
    '''
    def __init__(self, processor, input, varIdx, getNew, output, step, index, encoding, logger):
        '''
        processor:  line processor
        input:      input filename
        output:     a pipe to write output
        varIdx:     a dictionary of variants, used to only return matching lines if specified
        step:       step between processing lines because multiple workers might read the same file
        index:      index of this worker in the worker group
        encoding:   file encoding
        logger:     logger object to write error messages when things go wrong
        '''
        Process.__init__(self, name='FileReader')
        self.processor = processor
        self.input = input
        self.output = output
        self.step = step
        self.varIdx = varIdx
        self.getNew = getNew
        self.index = index
        self.encoding = encoding
        self.logger = logger

    def run(self): 
        first = True
        num_records = 0
        skipped_lines = 0
        line_no = 0
        with openFile(self.input) as input_file:
            for line in input_file:
                line_no += 1
                if line_no % self.step != self.index:
                    continue
                line = line.decode(self.encoding)
                try:
                    if line.startswith('#'):
                        continue
                    for bins,rec in self.processor.process(line):
                        if first:
                            self.output.send(self.processor.columnRange)
                            first = False
                        num_records += 1
                        if self.varIdx is not None:
                            var_key = (rec[0], rec[2], rec[3])
                            if self.getNew:
                                # only need new variant, so continue if variant in varIdx.
                                if var_key in self.varIdx and rec[1] in self.varIdx[var_key]:
                                    continue
                            else:
                                # only need existing variant, continue if variant not in varIdx
                                if var_key not in self.varIdx or rec[1] not in self.varIdx[var_key]:
                                    continue
                        self.output.send((line_no, bins, rec))
                except Exception as e:
                    self.logger.debug('Failed to process line "{}...": {}'.format(line[:20].strip(), e))
                    skipped_lines += 1
        # if still first (this thread has not read anything), still send the columnRange stuff
        if first:
            self.output.send(self.processor.columnRange)
        # everything is done, stop the pipe
        self.output.send(None)
        # and send the summary statistics
        self.output.send((num_records, skipped_lines))
        self.output.close()


class StandaloneTextReader:
    ''' This processor fire up 1 worker to read an input file
    and gather their outputs
    '''
    def __init__(self, processor, input, varIdx, getNew, encoding, logger):
        self.num_records = 0
        self.skipped_lines = 0
        #
        self.reader, w = Pipe(False)
        self.worker = ReaderWorker(processor, input, varIdx, getNew, w, 1, 0, encoding, logger)
        self.worker.start()
        # the send value is columnRange
        self.columnRange = self.reader.recv()
        
    def records(self):
        while True:
            val = self.reader.recv()
            if val is None:
                self.num_records, self.skipped_lines = self.reader.recv()
                break
            else:
                yield val
        self.worker.join()

class MultiTextReader:
    '''This processor fire up num workers to read an input file
    and gather their outputs
    '''
    def __init__(self, processor, input, varIdx, getNew, jobs, encoding, logger):
        self.readers = []
        self.workers = []
        self.num_records = 0
        self.skipped_lines = 0
        for i in range(jobs):
            r, w = Pipe(False)
            p = ReaderWorker(processor, input, varIdx, getNew, w, jobs, i, encoding, logger)
            self.readers.append(r)
            self.workers.append(p)
            p.start()
        # the send value is columnRange
        for reader in self.readers:
            self.columnRange = reader.recv()
        
    def records(self):
        all_workers = len(self.readers)
        still_working = len(self.readers)
        #
        # we need a heap to keep records read from multiple processes in order
        # we can not really guarantee this if there are large trunks of ignored
        # records but a heap size = 4 * number of readers should work in most cases
        #
        heap = []
        filled = False
        while True:
            for idx, reader in enumerate(self.readers):
                val = reader.recv()
                if val is None:
                    # one thread died ...
                    still_working -= 1
                    nr, sl = reader.recv()
                    self.num_records += nr
                    self.skipped_lines += sl
                    self.readers[idx] = None
                elif filled:
                    yield heappushpop(heap, val)
                else:
                    heappush(heap, val)
                    filled = len(heap) == 4 * len(self.readers)
            if still_working < all_workers:
                self.readers = [x for x in self.readers if x is not None]
                all_workers = len(self.readers)
                if all_workers == 0:
                    # return all things in the heap
                    for i in range(len(heap)):
                        yield heappop(heap)
                    break
        for p in self.workers:
            p.join()

#
#
# Write genotype to disk
# 
#

class GenotypeWriter:
    '''This class write genotypes to a genotype database, which does
    not have to be the project genotype database.'''
    def __init__(self, geno_db, geno_info, genotype_status, sample_ids, logger=None):
        '''geno_db:   genotype database
        geno_info:    genotype information fields
        sample_ids:   ID of samples that will be written to this database
        '''
        #
        self.db = DatabaseEngine()
        self.db.connect(geno_db)
        if genotype_status == 2:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'\
                .format(','.join([self.db.PH] * (2 + len(geno_info))))
        else:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'.format(self.db.PH)
        self.cur = self.db.cursor()
        if logger:
            s = delayedAction(logger.info, 'Creating {} genotype tables'.format(len(sample_ids)))
        for idx, sid in enumerate(sample_ids):
            # create table
            self.createNewSampleVariantTable(self.cur,
                'genotype_{0}'.format(sid), genotype_status == 2, geno_info)
        self.db.commit()
        if logger:
            del s
        #
        # number of genotype batches that has been written
        self.count = 0
        # cache of genotypes. This class will accumulate 1000 genotypes before
        # it writes to the disk using 'executemany', which will be faster than
        # calling 1000 'execute'.
        self.cache = {}

    def createNewSampleVariantTable(self, cur, table, genotype=True, fields=[]):
        '''Create a table ``genotype_??`` to store genotype data'''
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS {0} (
                variant_id INT NOT NULL
            '''.format(table) + 
            (', GT INT' + ''.join([', {} {}'.format(f.name, f.type) for f in fields]) if genotype else '')
            + ');'
         )
     
    def write(self, id, rec):
        try:
            if len(self.cache[id]) == 1000:
                # this will fail if id does not exist, so we do not need 
                # extra test if id is valid
                self.cur.executemany(self.query.format(id), self.cache[id])
                self.cache[id] = [rec]
                self.count += 1
            else:
                self.cache[id].append(rec)
        except KeyError:
            # if this is a new id
            self.cache[id] = [rec]
        if self.count % 1000 == 0:
            self.db.commit()
    
    def close(self):
        for id, val in self.cache.iteritems():
            if len(val) > 0:
                self.cur.executemany(self.query.format(id), val)
        self.db.commit()
        self.db.close()

#
#
# utility function to get sample name
#
def probeSampleName(filename, prober, encoding, logger):
    '''Probe text file for sample name. Essentially speaking
    this function will go to the last comment line, break it in pieces
    and see if we can grab some headers'''
    header_line = None
    count = 0
    with openFile(filename) as input:
        for line in input:
            line = line.decode(encoding)
            # the last # line
            if line.startswith('#'):
                header_line = line
            else:
                try:
                    for bins, rec in prober.process(line):
                        if header_line is None:
                            return len(rec), []
                        elif len(rec) == 0:
                            return 0, []
                        else:
                            cols = [x[0] for x in prober.fields]
                            if type(cols[0]) is tuple:
                                fixed = False
                                # mutiple ones, need to figure out the moving one
                                for i,idx in enumerate(prober.raw_fields[0].index.split(',')):
                                    if ':' in idx:
                                        cols = [x[i] for x in cols]
                                        fixed = True
                                        break
                                if not fixed:
                                    cols = [x[-1] for x in cols]
                            header = [x.strip() for x in header_line.split()] # #prober.delimiter)]
                            if max(cols) - min(cols)  < len(header) and len(header) > max(cols):
                                return len(rec), [header[len(header) - prober.nColumns + x] for x in cols]
                            else:
                                return len(rec), []
                except IgnoredRecord:
                    continue
                except Exception as e:
                    # perhaps not start with #, if we have no header, use it anyway
                    if header_line is None:
                        header_line = line
                    count += 1
                    if count == 100:
                        raise ValueError('No genotype column could be determined after 1000 lines.')
                    logger.debug(e)



#
#  A status object to control the import process. I cannot use a simple
#  queue solution because different types of importers are needed to handle
#  different inputs.
#
class ImportStatus:
    def __init__(self):
        self.manager = Manager()
        self.tasks = self.manager.dict()
        self.lock = Lock()
        self.total_sample_count = 0
        self.total_genotype_count = 0
        self.all_done = Value('i', 0)

    def add(self, item, num_lines):
        '''Add a job, each item has
        input_filename:     source filename
        tmp_file:           temporary genotype file
        start_sample, end_sample: index of samples
        sample_ids:         actual sample ids

        Status of the items can be

        (0,0):  initial input, not imported
        (1,x):  being imported
        (2,x):  imported
        (3,1):  being copied
        (4,1):  copied (done)

        x is 0, or 1 indicating if the genotype tables have been created.
        '''
        if item in self.tasks:
            raise RuntimeError('Item already added to ImportStatus')
        self.lock.acquire()
        self.tasks[item] = (0, 0)
        self.total_sample_count += len(item[-1])
        self.total_genotype_count += num_lines * len(item[-1])
        self.lock.release()

    def itemToImport(self, filelist):
        '''Try to find some item that can be imported. filelist is a list of files
        that the importer can work on
        '''
        ret = None
        self.lock.acquire()
        for item in sorted(self.tasks.keys()):
            # only check for status 0 items
            if self.tasks[item][0] == 0 and item[0] in filelist:
                self.tasks[item] = (1, self.tasks[item][1])
                ret = item
                break
        self.lock.release()
        return ret

    def pendingImport(self):
        # return waiting (queued, to be imported) and ongoing import
        return sum([x[0] == 0 for x in self.tasks.values()]), sum([x[0] == 1 for x in self.tasks.values()])

    def pendingCopy(self):
        # return waiting 
        return sum([x[0] == 2 for x in self.tasks.values()])

    def numItems(self):
        return len(self.tasks)

    def itemToCopy(self):
        '''Try to find some item that can be copied. Another field
        indicating whether or not the tables have been created.'''
        ret = None
        self.lock.acquire()
        for item in sorted(self.tasks.keys()):
            # only check for status 2 items
            if self.tasks[item][1] == 0:
                self.tasks[item] = (self.tasks[item][0], 1)
                ret = (item, 0)
                break
            if self.tasks[item] == (2, 1):
                self.tasks[item] = (3, 1)
                ret = (item, 1)
                break
        self.lock.release()
        return ret

    def set(self, item, status):
        self.lock.acquire()
        self.tasks[item] = (status, self.tasks[item][1])
        self.lock.release()

#
#   A working process to import genotype from a file, or part of a file
#   and write to a temporary genotype database.
#
# 
class GenotypeImportWorker(Process):
    '''This class starts a process, import genotype to a temporary genotype database.'''
    def __init__(self, variantIndex, filelist, processor, encoding, genotype_field, genotype_info, ranges, 
        geno_count, proc_index, status, logger):
        '''
        variantIndex: a dictionary that returns ID for each variant.
        filelist: files from which variantIndex is created. If the passed filename
            is not in this list, this worker will suicide so that it can be replaced 
            by a worker with more variants.
        encoding, genotypefield, genotype_info, ranges:  parameters to import data
        geno_count:  a shared variable to report number of genotype imported
        status:      an ImportStatus object to monitor the progress
        '''
        Process.__init__(self, name='GenotypeImporter')
        self.variantIndex = variantIndex
        self.filelist = filelist
        self.encoding = encoding
        self.processor = processor
        #
        self.genotype_field = genotype_field
        self.genotype_info = genotype_info
        self.ranges = ranges
        #
        self.geno_count = geno_count
        self.proc_index = proc_index
        #
        self.status = status
        self.logger = logger

    def run(self): 
        self.logger.debug('Importer {} starts with variants from {} files'.format(self.proc_index, len(self.filelist)))
        while True:
            item = self.status.itemToImport(self.filelist)
            if item is None:
                # wait a second, make sure there is no job
                time.sleep(1)
                item = self.status.itemToImport(self.filelist)
                if item is None:
                    self.logger.debug('Importer {} exits normally'.format(self.proc_index))
                    break
            # get parameters
            self.input_filename, self.genotype_file, self.genotype_status, start_sample, end_sample, self.sample_ids = item
            self.processor.reset(import_var_info=False, 
                import_sample_range=[0,0] if self.genotype_status == 1 else [start_sample, end_sample])
            self.count = [0, 0]
            # use the last value as start
            self.start_count = self.geno_count.value
            start_import_time = time.time()
            self._importData()
            # set the status to be imported (2)
            self.status.set(item, 2)
            end_import_time = time.time()
            todo, going = self.status.pendingImport()
            self.logger.debug('Importing {} samples ({} - {}) to {} took importer {} {:.1f} seconds, {} onging, {} to go.'.format(
                len(self.sample_ids), min(self.sample_ids), max(self.sample_ids), os.path.basename(self.genotype_file),
                self.proc_index, end_import_time - start_import_time, going, todo))

    def _importData(self):
        self.logger.debug('Importer {} starts importing genotypes for {} samples ({} - {})'.format(self.proc_index, len(self.sample_ids),
            min(self.sample_ids), max(self.sample_ids)))
        reader = TextReader(self.processor, self.input_filename, None, False, 0, self.encoding, self.logger)
        writer = GenotypeWriter(self.genotype_file, self.genotype_info, self.genotype_status, self.sample_ids)
        fld_cols = None
        last_count = 0
        for self.count[0], bins, rec in reader.records():
            try:
                variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
            except KeyError:
                self.logger.debug('Variant {} {} {} {} not found'.format(rec[0], rec[1], rec[2], rec[3]))
                continue
            # if there is genotype 
            if self.genotype_status == 2:
                if fld_cols is None:
                    col_rngs = [reader.columnRange[x] for x in range(self.ranges[2], self.ranges[4])]
                    fld_cols = []
                    for idx in range(len(self.sample_ids)):
                        fld_cols.append([sc + (0 if sc + 1 == ec else idx) for sc,ec in col_rngs])
                    if col_rngs[0][1] - col_rngs[0][0] != len(self.sample_ids):
                        self.logger.error('Number of genotypes ({}) does not match number of samples ({})'.format(
                            col_rngs[0][1] - col_rngs[0][0], len(self.sample_ids)))
                for idx, id in enumerate(self.sample_ids):
                    try:
                        # variant info is not read, ranges[1] should be used because ranges[2] is the index after
                        # variant index
                        if rec[self.ranges[1] + idx] is not None:
                            self.count[1] += 1
                            writer.write(id, [variant_id] + [rec[c] for c in fld_cols[idx]])
                    except IndexError:
                        self.logger.warning('Incorrect number of genotype fields: {} fields found, {} expected for record {}'.format(
                            len(rec), fld_cols[-1][-1] + 1, rec))
            else:
                # should have only one sample
                for id in self.sample_ids:
                    writer.write(id, [variant_id])
            if self.count[0] - last_count > 100:
                self.geno_count.value = self.start_count + self.count[0] * len(self.sample_ids)
                last_count = self.count[0]
        writer.close()


class GenotypeCopier(Process):
    def __init__(self, main_genotype_file, genotype_info, copied_samples, status, logger):
        '''copied_samples is a shared variable that should be increased with
        each sample copy.
        '''
        Process.__init__(self)
        self.main_genotype_file = main_genotype_file
        self.genotype_info = genotype_info
        self.copied_samples = copied_samples
        self.status = status
        self.logger = logger

    def run(self):
        self.db = DatabaseEngine()
        self.db.connect(self.main_genotype_file)
        file_count = 0
        while True:
            item = self.status.itemToCopy()
            if item is None:
                if self.status.all_done.value == 1:
                    self.db.close()
                    self.logger.debug('Genotype of {} samples are copied from {} files'.format(
                        self.copied_samples.value, file_count))
                    break
                time.sleep(1)
            elif item[1] == 0:   # create table only
                infile, self.genotype_file, self.genotype_status, ss, es, self.sample_ids = item[0]
                start_create_time = time.time()
                self._createTables()
                self.db.commit()
                end_create_time = time.time()
                self.logger.debug('Creating genotype tables for {} samples ({} - {}) took {:.1f} seconds.'.format(
                    len(self.sample_ids), min(self.sample_ids), max(self.sample_ids),
                    end_create_time - start_create_time))
            else:
                infile, self.genotype_file, gs, ss, es, self.sample_ids = item[0]
                start_copy_time = time.time()
                # get parameters
                self._copySamples()
                # set the status of the item to be copied (4)
                self.status.set(item[0], 4)
                end_copy_time = time.time()
                self.logger.debug('Copying {} samples ({} - {}) from {} took {:.1f} seconds, {} to go.'.format(
                    len(self.sample_ids), min(self.sample_ids), max(self.sample_ids),
                    os.path.basename(self.genotype_file), end_copy_time - start_copy_time, self.status.pendingCopy()))
                file_count += 1


    def createNewSampleVariantTable(self, cur, table, genotype=True, fields=[]):
        '''Create a table ``genotype_??`` to store genotype data'''
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS {0} (
                variant_id INT NOT NULL
            '''.format(table) + 
            (', GT INT' + ''.join([', {} {}'.format(f.name, f.type) for f in fields]) if genotype else '')
            + ');'
         )

    def _createTables(self):
        cur = self.db.cursor()
        for ID in self.sample_ids:
            self.createNewSampleVariantTable(cur, 'genotype_{}'.format(ID), self.genotype_status == 2, self.genotype_info)
        
    def _copySamples(self):
        self.db.attach(self.genotype_file, '__from')
        # start copying genotype
        # copy genotype table
        cur = self.db.cursor()
        for count, ID in enumerate(self.sample_ids):
            query = 'INSERT INTO genotype_{0} SELECT * FROM __from.genotype_{0};'.format(ID)
            cur.execute(query)
            self.db.commit()
            # update progress
            self.copied_samples.value += 1
        self.db.detach('__from')
        # remove the temporary file
        try:
            os.remove(self.genotype_file)
        except:
            pass
        
#
#
#  Command import
#
#

class Importer:
    '''A general class for importing variants'''
    def __init__(self, proj, files, build, format, sample_name=None, force=False, jobs=1, fmt_args=[]):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        self.sample_in_file = []
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
        #
        # for #record, #genotype (new or updated), #new variant, SNV, insertion, deletion, complex variants, invalid record, updated record
        self.count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.total_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.import_alt_build = False
        #
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                self.logger.info('Using primary reference genome {} of the project.'.format(self.build))
            else:
                raise ValueError('Please specify a reference genome using parameter --build')
        else:
            self.build = build
        #
        if self.proj.build is None:
            self.proj.setRefGenome(self.build)
        elif self.build == self.proj.build:
            # perfect case
            pass
        elif self.build == self.proj.alt_build:
            # troublesome
            self.import_alt_build = True
        elif self.proj.alt_build is None:
            # even more troublesome
            self.logger.warning('The new files uses a different reference genome ({}) from the primary reference genome ({}) of the project.'.format(self.build, self.proj.build))
            self.logger.info('Adding an alternative reference genome ({}) to the project.'.format(self.build))
            tool = LiftOverTool(self.proj)
            # in case of insert, the indexes will be dropped soon so do not build
            # index now
            tool.setAltRefGenome(self.build, build_index=False)
            self.import_alt_build = True
        #
        self.jobs = max(1, jobs)
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        # try to guess file type
        if not format:
            filename = files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            fmt = fileFMT(format, fmt_args, logger=self.logger)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        if fmt.preprocessor is not None:
            self.logger.info('Preprocessing files {} to generate intermediate input files for import'.format(', '.join(files)))
            # if this is the case, only one input stream will be allowed.
            # process command line
            command = fmt.preprocessor
            # replace command with other stuff, if applicable
            command = command.replace('$build', "'{}'".format(self.build))
            #
            # create temp files
            temp_files = [os.path.join(runOptions.cache_dir, os.path.basename(x) + '.' + fmt.name) for x in files]
            try:
                processor = eval(command)
                # intermediate files will be named as "cache_dir/$inputfilename.$(fmt.name)"
                processor.convert(files, temp_files, self.logger)
                for output in temp_files:
                    if not os.path.isfile(output):
                        raise ValueError("Preprocessed file {} does not exist.".format(output))
            except Exception as e:
                raise ValueError("Failed to execute preprocessor '{}': {}".format(command, e))
            #
            # we record file as cache files
            files = temp_files
        #
        self.files = []
        cur = self.db.cursor()
        cur.execute('SELECT filename from filename;')
        existing_files = [x[0] for x in cur.fetchall()]
        for filename in files:
            if filename in existing_files:
                if force:
                    self.logger.info('Re-importing imported file {}'.format(filename))
                    IDs = proj.selectSampleByPhenotype('filename = "{}"'.format(filename))
                    self.proj.db.attach(self.proj.name + '_genotype')
                    proj.removeSamples(IDs)
                    self.proj.db.detach(self.proj.name + '_genotype')
                    # remove file record
                    cur = self.db.cursor()
                    cur.execute('DELETE FROM filename WHERE filename={};'.format(self.db.PH), (filename,))
                    self.db.commit()
                    self.files.append(filename)
                else:
                    self.logger.info('Ignoring imported file {}'.format(filename))
            elif not os.path.isfile(filename):
                raise ValueError('File {} does not exist'.format(filename))
            else:
                self.files.append(filename)
        #
        if len(self.files) == 0:
            raise ValueError('No file to import')
        #
        self.sample_name = sample_name
        #
        # how to split processed records
        #
        self.ranges = fmt.ranges
        self.variant_fields = [x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]]
        self.variant_info = [x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]]
        #
        # NOTE: for a given format, number of genotype_fields and associated genotype info
        # are determined by the format file, and are fixed, but the number of samples in
        # a file can vary. In particular, there can be no sample, and in this case, if a
        # sample name is given, a sample will be created with only variant ids.
        #
        self.genotype_field = [x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]]
        self.genotype_info = [x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]]
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for pos, ref, alt are 1, 2, 3 in fields.
            self.processor = LineProcessor(fmt.fields, [(1, 2, 3)], fmt.delimiter, self.ranges, self.logger)
        else:  # position or range type
            raise ValueError('Can only import data with full variant information (chr, pos, ref, alt)')
        # probe number of samples
        if self.genotype_field:
            self.prober = LineProcessor([fmt.fields[fmt.ranges[2]]], [], fmt.delimiter, self.ranges, self.logger)
        # there are variant_info
        if self.variant_info:
            cur = self.db.cursor()
            headers = self.db.getHeaders('variant')
            for f in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]:
                # either insert or update, the fields must be in the master variant table
                self.proj.checkFieldName(f.name, exclude='variant')
                if f.name not in headers:
                    s = delayedAction(self.logger.info, 'Adding column {}'.format(f.name))
                    cur.execute('ALTER TABLE variant ADD {} {};'.format(f.name, f.type))
                    del s
        #
        if fmt.input_type != 'variant':
            self.logger.info('Only variant input types that specifies fields for chr, pos, ref, alt could be imported.')
        #
        self.input_type = fmt.input_type
        self.encoding = fmt.encoding
        fbin, fchr, fpos = ('alt_bin', 'alt_chr', 'alt_pos') if self.import_alt_build else ('bin', 'chr', 'pos')
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH)
        self.variant_insert_query = 'INSERT INTO variant ({0}, {1}, {2}, ref, alt {3}) VALUES ({4});'\
            .format(fbin, fchr, fpos, ' '.join([', ' + x for x in self.variant_info]), ', '.join([self.db.PH]*(len(self.variant_info) + 5)))
        #
        self.variantIndex = self.proj.createVariantMap('variant', self.import_alt_build)
        # drop index here after all possible exceptions have been raised.
        self.proj.dropIndexOnMasterVariantTable()


    def __del__(self):
        # remove existing indexes, which will be created when the project is open
        # by a non-import command
        self.proj.dropIndexOnMasterVariantTable()

    def recordFileAndSample(self, filename, sampleNames):
        cur = self.db.cursor()
        # get header of file
        header = ''
        with openFile(filename) as input:
            for line in input:
                line = line.decode(self.encoding)
                if line.startswith('#'):
                    header += line
                else:
                    break
        cur.execute("INSERT INTO filename (filename, header) VALUES ({0}, {0});".format(self.db.PH), (filename, header))
        filenameID = cur.lastrowid
        sample_ids = []
        s = delayedAction(self.logger.info, 'Creating {} genotype tables'.format(len(sampleNames)))
        #
        for samplename in sampleNames:
            cur.execute('INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'.format(self.db.PH),
                (filenameID, samplename))
            sid = cur.lastrowid
            sample_ids.append(sid)
        del s
        return sample_ids

    def addVariant(self, cur, rec):
        #
        if rec[4] == '-':
            self.count[5] += 1
        elif rec[3] == '-':
            self.count[4] += 1
        elif len(rec[4]) == 1 and len(rec[3]) == 1:
            self.count[3] += 1
        else:
            self.count[6] += 1
        #
        var_key = tuple((rec[1], rec[3], rec[4]))
        if var_key in self.variantIndex and rec[2] in self.variantIndex[var_key]:
            variant_id = self.variantIndex[var_key][rec[2]][0]
            if len(rec) > 5:
                self.count[8] += 1
                cur.execute(self.update_variant_query, rec[5:] + [variant_id])
            return variant_id
        else:
            # new variant!
            # alt_chr and alt_pos are updated if adding by alternative reference genome
            self.count[2] += 1
            cur.execute(self.variant_insert_query, rec)
            variant_id = cur.lastrowid
            # one for new variant
            if var_key in self.variantIndex:
                self.variantIndex[var_key][rec[2]] = (variant_id, 1)
            else:
                self.variantIndex[var_key] = {rec[2]: (variant_id, 1)}
            return variant_id

    def getSampleIDs(self, input_filename):
        '''Return a list of sample ids, and a status field indicating
        0:  no sample. sample id has to be empty
        1:  no genotype but with sample (a list of variant ids, no genotype), there can be only one sample id.
        2:  genotype with samples.
        '''
        if not self.sample_name:
            # if no sample name is specified
            if not self.genotype_field:
                self.logger.warning('Sample information is not recorded for a file without genotype and sample name.')
                self.sample_in_file = []
                return ([], 0)
            else:
                try:
                    numSample, names = probeSampleName(input_filename, self.prober, self.encoding, self.logger)
                    if not names:
                        if numSample == 1:
                            self.logger.debug('Missing sample name (name None is used)'.format(numSample))
                            self.sample_in_file = [None]
                            return (self.recordFileAndSample(input_filename, [None]), 2)
                        elif numSample == 0:
                            self.logger.debug('No genotype column exists in the input file so no sample will be recorded.')
                            self.sample_in_file = []
                            return ([], 0)
                        else:
                            raise ValueError('Failed to guess sample name. Please specify sample names for {} samples using parameter --sample_name, or add a proper header to your input file. See "vtools import_variants -h" for details.'.format(numSample))
                    else:
                        self.sample_in_file = [x for x in names]
                        return (self.recordFileAndSample(input_filename, names), 2)
                except ValueError:
                    # cannot find any genotype column, perhaps no genotype is defined in the file (which is allowed)
                    self.logger.warning('No genotype column could be found from the input file. Assuming no genotype.')
                    self.sample_in_file = []
                    return ([], 0)
        else:
            self.sample_in_file = [x for x in self.sample_name]
            if not self.genotype_field:
                # if no genotype, but a sample name is given
                self.logger.debug('Input file does not contain any genotype. Only the variant ownership information is recorded.')
                return (self.recordFileAndSample(input_filename, self.sample_name), 1)
            else:
                try:
                    numSample, names = probeSampleName(input_filename, self.prober, self.encoding, self.logger)
                except ValueError as e:
                    self.logger.debug(e)
                    numSample = 0
                #
                if numSample == 0:
                    self.logger.warning('No genotype column could be found from the input file. Assuming no genotype.')
                    # remove genotype field from processor
                    self.processor.reset(import_var_info=True, import_sample_range=[0,0])
                    if len(self.sample_name) > 1:
                        raise ValueError("When there is no sample genotype, only one sample name is allowed.")
                elif len(self.sample_name) != numSample:
                    raise ValueError('{} sample detected but only {} sample names are specified'.format(numSample, len(self.sample_name)))                        
                return (self.recordFileAndSample(input_filename, self.sample_name), 1 if numSample == 0 else 2)
 
    def importVariantAndGenotype(self, input_filename):
        '''Input variant and genotype at the same time, appropriate for cases with
        no or one sample in a file'''
        self.processor.reset()
        if self.genotype_field:
            self.prober.reset()
        #
        sample_ids, genotype_status = self.getSampleIDs(input_filename)
        #
        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc//200, 100), 100000)
        # one process is for the main program, the
        # other threads will handle input
        reader = TextReader(self.processor, input_filename, None, False, self.jobs - 1, self.encoding, self.logger)
        if genotype_status != 0:
            writer = GenotypeWriter(
                # write directly to the genotype table
                '{}_genotype'.format(self.proj.name),
                self.genotype_info, genotype_status,
                sample_ids, self.logger)
        # preprocess data
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0
        fld_cols = None
        for self.count[0], bins, rec in reader.records():
            variant_id = self.addVariant(cur, bins + rec[0:self.ranges[2]])
            if genotype_status == 2:
                if fld_cols is None:
                    col_rngs = [reader.columnRange[x] for x in range(self.ranges[2], self.ranges[4])]
                    fld_cols = []
                    for idx in range(len(sample_ids)):
                        fld_cols.append([sc + (0 if sc + 1 == ec else idx) for sc,ec in col_rngs])
                    if col_rngs[0][1] - col_rngs[0][0] != len(sample_ids):
                        self.logger.error('Number of genotypes ({}) does not match number of samples ({})'.format(
                            col_rngs[0][1] - col_rngs[0][0], len(sample_ids)))
                for idx, id in enumerate(sample_ids):
                    try:
                        if rec[self.ranges[2] + idx] is not None:
                            self.count[1] += 1
                            writer.write(id, [variant_id] + [rec[c] for c in fld_cols[idx]])
                    except IndexError:
                        self.logger.warning('Incorrect number of genotype fields: {} fields found, {} expected for record {}'.format(
                            len(rec), fld_cols[-1][-1] + 1, rec))
            elif genotype_status == 1:
                # should have only one sample
                for id in sample_ids:
                    writer.write(id, [variant_id])
            if (last_count == 0 and self.count[0] > 200) or (self.count[0] - last_count > update_after):
                self.db.commit()
                last_count = self.count[0]
                prog.update(self.count[0])
        prog.done()
        self.count[7] = reader.skipped_lines
        # stop writers
        if genotype_status != 0:
            writer.close()
        self.db.commit()
       
    def importVariant(self, input_filename):
        # reset text processor to allow the input of files with different number of columns
        self.processor.reset(import_var_info=True, import_sample_range=[0,0])
        #
        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc//200, 100), 100000)
        # one process is for the main program, the
        # other threads will handle input
        # getNew=True so the reader only read variants not in variantIndex if no additional
        # variant info is imported
        if self.variant_info:
            reader = TextReader(self.processor, input_filename, None, True, 
                runOptions.import_num_of_readers, self.encoding, self.logger)
        else:
            reader = TextReader(self.processor, input_filename, self.variantIndex, True,
                runOptions.import_num_of_readers, self.encoding, self.logger)
        # preprocess data
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0
        fld_cols = None
        for self.count[0], bins, rec in reader.records():
            variant_id = self.addVariant(cur, bins + rec[0:self.ranges[2]])
            if (last_count == 0 and self.count[0] > 200) or (self.count[0] - last_count > update_after):
                self.db.commit()
                last_count = self.count[0]
                prog.update(self.count[0])
        prog.done()
        self.count[7] = reader.skipped_lines
        self.db.commit()

    def finalize(self):
        # this function will only be called from import
        cur = self.db.cursor()
        total_new = sum(self.total_count[3:7])
        if total_new > 0:
            # analyze project to get correct number of rows for the master variant table
            self.proj.analyze(force=True)
        if total_new == 0 or self.proj.alt_build is None:
            # if no new variant, or no alternative reference genome, do nothing
            return
        # we need to run lift over to convert coordinates before importing data.
        tool = LiftOverTool(self.proj)
        to_be_mapped = os.path.join(runOptions.temp_dir, 'var_in.bed')
        loci_count = 0
        with open(to_be_mapped, 'w') as output:
            for key in self.variantIndex:
                for pos, status in self.variantIndex[key].iteritems():
                    if status[1] == 1:
                        output.write('{0}\t{1}\t{2}\t{3}/{4}/{5}\n'.format(key[0] if len(key[0]) > 2 else 'chr' + key[0],
                           pos - 1, pos, key[1], key[2], status[0]))
                        loci_count += 1
        # free some RAM
        self.variantIndex.clear()
        #
        if self.import_alt_build:
            self.logger.info('Mapping new variants at {} loci from {} to {} reference genome'.format(loci_count, self.proj.alt_build, self.proj.build))
            query = 'UPDATE variant SET bin={0}, chr={0}, pos={0} WHERE variant_id={0};'.format(self.db.PH)
            mapped_file, err_count = tool.mapCoordinates(to_be_mapped, self.proj.alt_build, self.proj.build)
        else:
            self.logger.info('Mapping new variants at {} loci from {} to {} reference genome'.format(loci_count, self.proj.build, self.proj.alt_build))
            query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(self.db.PH)
            # this should not really happen, but people (like me) might manually mess up with the database
            s = delayedAction(self.logger.info, 'Adding alternative reference genome {} to the project.'.format(self.proj.alt_build))
            headers = self.db.getHeaders('variant')
            for fldName, fldType in [('alt_bin', 'INT'), ('alt_chr', 'VARCHAR(20)'), ('alt_pos', 'INT')]:
                if fldName in headers:
                    continue
                self.db.execute('ALTER TABLE variant ADD {} {} NULL;'.format(fldName, fldType))
            del s
            mapped_file, err_count = tool.mapCoordinates(to_be_mapped, self.proj.build, self.proj.alt_build)
        # update records
        prog = ProgressBar('Updating coordinates', total_new)
        # 1: succ mapped
        count = 0
        with open(mapped_file) as var_mapped:
            for line in var_mapped.readlines():
                try:
                    chr, start, end, name = line.strip().split()
                    ref, alt, var_id = name.split('/')
                    if chr.startswith('chr'):
                        chr = chr[3:]
                    pos = int(start) + 1
                    var_id = int(var_id)
                except:
                    continue
                cur.execute(query, (getMaxUcscBin(pos - 1, pos), chr, pos, var_id))
                count += 1
                if count % self.db.batch == 0:
                    self.db.commit()
                    prog.update(count)
        self.db.commit()
        prog.done()
        self.logger.info('Coordinates of {} ({} total, {} failed to map) new variants are updated.'\
            .format(count, total_new, err_count))
            
    def importFilesSequentially(self):
        '''import files one by one, adding variants along the way'''
        sample_in_files = []
        for count,f in enumerate(self.files):
            self.logger.info('{} variants and genotypes from {} ({}/{})'.format('Importing', f, count + 1, len(self.files)))
            self.importVariantAndGenotype(f)
            total_var = sum(self.count[3:7])
            self.logger.info('{:,} variants ({:,} new{}) from {:,} lines are imported, {}.'\
                .format(total_var, self.count[2],
                    ''.join([', {:,} {}'.format(x, y) for x, y in \
                        zip(self.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                    self.count[0],
                    'no sample is created' if len(self.sample_in_file) == 0 else 'with a total of {:,} genotypes from {}'.format(
                        self.count[1], 'sample {}'.format(self.sample_in_file[0]) if len(self.sample_in_file) == 1 else '{:,} samples'.format(len(self.sample_in_file)))))
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
            sample_in_files.extend(self.sample_in_file)
        if len(self.files) > 1:
            total_var = sum(self.total_count[3:7])
            self.logger.info('{:,} variants ({:,} new{}) from {:,} lines are imported, {}.'\
                .format(total_var, self.total_count[2],
                    ''.join([', {:,} {}'.format(x, y) for x, y in \
                        zip(self.total_count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                    self.total_count[0],
                    'no sample is created' if len(sample_in_files) == 0 else 'with a total of {:,} genotypes from {}'.format(
                        self.total_count[1], 'sample {}'.format(sample_in_files[0]) if len(sample_in_files) == 1 else '{:,} samples'.format(len(sample_in_files)))))

    def importFilesInParallel(self):
        '''import files in parallel, by importing variants and genotypes separately, and in their own processes. 
        More specifically, suppose that there are three files

        file1: variant1, sample1.1, sample1.2, sample1.3
        file2: variant2, sample2.1, sample2.2
        file3: variant3, sample3.1, sample3.2

        where variant1, 2, and 3 are three potentially overlapping sets of variants.
        sample1, sample2, sample3 are three groups of samples that are divided by the number of samples
        and number of processes (self.jobs) (say, 3000 samples divided into three groups of 1000 samples).

        Then, there are 

        A: 2 TextReader processes <--> main process read variant1, 2, 3 --> master variant table
        B: self.jobs GenotypeImportWorker reads sample 1.1, 1.2, 1.3, 2.1 etc ... --> temporary genotype tables
        C: 1 GenotypeCopier -> copy temporary genotype tables to the master genotype table

        The processes are organized so that 
        1. sample A.x if read after variant A is read
        2. genotypes are copied after temporary tables are completed

        Because of the processes work together, although there are multiple progress bars, they
        might start from the middle because some work might have already been done in the previous step.
        '''
        importers = [None] * self.jobs
        # number of genotypes each process have imported
        genotype_import_count = [Value('i', 0) for x in range(self.jobs)]
        # number of sample copied
        sample_copy_count = Value('i', 0)
        # import queue that accept jobs sample 1.1, 1.2, etc
        status = ImportStatus()
        # start copier
        copier = GenotypeCopier('{}_genotype.DB'.format(self.proj.name), 
            self.genotype_info, sample_copy_count, status, self.logger)
        copier.start()
        #
        # The logic of importer is complex here. Because an importer needs to know variantIndex to 
        # write genotype tables, and because importers starts after each file is read, an importer
        # created earlier (eg. from file 1) cannot be used to import genotype for file 2. This is
        # why we
        #
        # 1. create importer only after a file is processed.
        # 2. If all workers are busy, no more new worker will be created.
        # 3. When an old importer finds that it cannot process new genotype file, it will commit succide.
        # 4. The master process will create new importers with updated variantIndex when
        #    there are empty slots.
        #
        # We do not pass variantIndex to importers because variantIndex is huge and it is slow
        # to pass it with other parameters.
        #
        # process each file
        for count, input_filename in enumerate(self.files):
            self.logger.info('{} variants from {} ({}/{})'.format('Importing', input_filename, count + 1, len(self.files)))
            self.importVariant(input_filename)
            self.logger.info('{:,} new variants ({}) from {:,} lines are imported.'\
                .format(self.count[2],
                    ', '.join(['{:,} {}'.format(x, y) for x, y in \
                        zip(self.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                    self.count[0]))
            # genotypes?
            if self.genotype_field:
                self.prober.reset()
            # if there are samples?
            sample_ids, genotype_status = self.getSampleIDs(input_filename)
            #
            # we should have file line count from importVariant
            num_of_lines = self.count[0]
            #
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
            #
            if len(sample_ids) == 0:
                continue
            #
            # determine workload:
            # from our benchmark, if we split jobs evenly, the last trunk will
            # take about double time to finish because the extra time to reach the end
            # if the lines are long. Therefore, we are using an algorithm that the last
            # piece will handle half of the samples of the first one.
            #
            # n -- len(sample_ids)
            # m -- number of processes
            # k -- workload of the first one
            # k - k/2(m-1) -- workload of the second one
            # ...
            # k/2 -- workload of the last one
            # 
            # k =  4n/3m
            #
            # each process handle at least 10 samples
            workload = [max(10, int(4*len(sample_ids)/(3.*self.jobs)*(1-x/(2.*(self.jobs - 1))))) for x in range(self.jobs)]
            # if there are missing ones, give it to workload
            workload[0] += max(0, len(sample_ids) - sum(workload))
            start_sample = 0
            for job in range(self.jobs):
                if workload[job] == 0:
                    continue
                end_sample = min(start_sample + workload[job], len(sample_ids))
                if end_sample <= start_sample:
                    continue
                # tell the processor do not import variant info, import part of the sample
                tmp_file = os.path.join(runOptions.temp_dir, 'tmp_{}_{}_genotype.DB'.format(count, job))
                if os.path.isfile(tmp_file):
                    os.remove(tmp_file)
                if os.path.isfile(tmp_file):
                    raise RuntimeError('Failed to remove existing temporary database {}. Remove clean your cache directory'.format(tmp_file))
                # send a import job to the importer workers
                status.add((input_filename, tmp_file, genotype_status, start_sample, end_sample, 
                    tuple(sample_ids[start_sample : end_sample])), num_of_lines)
                # start an importer if needed
                for i in range(self.jobs):
                    if importers[i] is None or not importers[i].is_alive():
                        importers[i] = GenotypeImportWorker(self.variantIndex, self.files[:count+1], 
                            self.processor, self.encoding, self.genotype_field, self.genotype_info, self.ranges,
                            genotype_import_count[i], i, status, self.logger)
                        importers[i].start()
                        break
                start_sample = end_sample
        # 
        # monitor the import of genotypes
        prog = ProgressBar('Importing genotypes', status.total_genotype_count, initCount=sum([x.value for x in genotype_import_count]))
        while True:
            # each process update their passed shared value
            # the master process add them and get the total number of genotypes imported
            prog.update(sum([x.value for x in genotype_import_count]))
            # 
            queued, importing = status.pendingImport()
            if queued + importing == 0:
                # the importer will kill themselves after there is no pending job
                prog.done()
                break
            # create importers if any of the importer is gone. This might be a waste of resource
            # but an importer that is pending does not cost much
            if queued > 0:
                new_count = 0
                for i in range(self.jobs):
                    if importers[i] is None or not importers[i].is_alive():
                        importer = GenotypeImportWorker(self.variantIndex, self.files, self.processor, self.encoding,
                            self.genotype_field, self.genotype_info, self.ranges,
                            genotype_import_count[i], i, status, self.logger)
                        importers[i] = importer
                        importer.start()
                        new_count += 1
                    if new_count >= queued:
                        break
            time.sleep(2)
        # monitor the copy of genotypes
        prog = ProgressBar('Copying genotype', status.total_sample_count, initCount=sample_copy_count.value)
        while True:
            prog.update(sample_copy_count.value)
            if sample_copy_count.value == status.total_sample_count:
                prog.done()
                status.all_done.value = 1
                copier.join()
                break
            time.sleep(1)
        # final status line
        if len(self.files) > 1:
            self.logger.info('{:,} new variants ({}) from {:,} lines ({:,} samples) are imported.'\
                .format(self.total_count[2],
                    ', '.join(['{:,} {}'.format(x, y) for x, y in \
                        zip(self.total_count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                    self.total_count[0], status.total_sample_count))


def importVariantsArguments(parser):
    parser.add_argument('input_files', nargs='+',
        help='''A list of files that will be imported. The file should be delimiter
            separated with format described by parameter --format. Gzipped files are
            acceptable. If a preprocessor is defined in the format, processor specific
            parameters will be passed as $1, $2 etc to generate a single input stream.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, it is assumed to be the primary reference genome of the project.
            If a reference genome that is different from the primary reference genome of the
            project is specified, it will become the alternative reference genome of the
            project. The UCSC liftover tool will be automatically called to map input
            coordinates between the primary and alternative reference genomes. If you
            are uncertain about the reference genome used in the data, use a recent
            standard reference genome (e.g. hg19) and validate it after the data is
            imported (c.f. "vtools admin --validate_build").''')
    parser.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as VCF (cf. 'vtools show formats'), or a 
            local format specification file (with extension .fmt). If unspecified,
            variant tools will try to guess format from file extension. Some file
            formats accept parameters (cf. 'vtools show format FMT') and allow you
            to import additional or alternative fields defined for the format. ''')
    parser.add_argument('--sample_name', nargs='*', default=[],
        help='''Name of the samples imported by the input files. The same names will be
            used for all files if multiple files are imported. If unspecified, headers
            of the genotype columns of the last comment line (line starts with #) of the
            input files will be used (and thus allow different sample names for input files).
            If sample names are specified for input files without genotype, samples
            will be created without genotype. If sample names cannot be determined from
            input file and their is no ambiguity (only one sample is imported), a sample
            with NULL sample name will be created.''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.'''),
    parser.add_argument('-j', '--jobs', metavar='N', default=4, type=int,
        help='''Number of processes to import input file. Variant tools by default
            uses four processes to import variants and samples genotypes in 
            parallel, and you can use more or less processes by adjusting this
            parameter. Due to the overhead of inter-process communication, more
            jobs do not automatically lead to better performance.''')

def importVariants(args):
    try:
        # the project is opened with verify=False so index on the master
        # variant table will not be created if it does not exist (because the
        # last command was a vtools import command)
        with Project(verbosity=args.verbosity, verify=False) as proj:
            importer = Importer(proj=proj, files=args.input_files,
                build=args.build, format=args.format, sample_name=args.sample_name,
                force=args.force, jobs=args.jobs, fmt_args=args.unknown_args)
            if args.jobs == 1:
                # if jobs == 1, use the old algorithm that insert variant and
                # genotype together ...
                importer.importFilesSequentially()
            else:
                importer.importFilesInParallel()
            importer.finalize()
        proj.close()
    except Exception as e:
        sys.exit(e)
