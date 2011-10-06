#!/usr/bin/env python
#
# $File: importer.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
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
import gzip
import re
from itertools import izip, repeat
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, normalizeVariant

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
            return item.split(self.sep)[self.index]
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


class FieldFromFormat:
    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that return the value of a field according 
        to a format string. This is used to extract stuff from the format
        string of vcf files.
        '''
        self.name = name
        self.sep = sep
        self.fmt = '\t'
        self.idx = None
        self.default = default

    def __call__(self, item):
        if not item[0] == self.fmt:
            fmt, val = item
            self.fmt = fmt
            fields = fmt.split(self.sep)
            if self.name in fields:
                self.idx = fields.index(self.name)
                return val.split(self.sep)[self.idx]
            else:
                self.idx = None
                return self.default
        try:
            return item[1].split(self.sep)[self.idx]
        except:
            # in the case that self.idx is None, or if item[1].split(self.sep) does not have enough items
            return self.default

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
            # if item is None or ''
            if not item:
                return item
            # if multiple records are returned, apply to each of them
            if type(item) is tuple:
                if type(item[0]) is tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item


class LineImporter:
    '''An intepreter that read a record, process it and return processed records.'''
    def __init__(self, fields, build, delimiter, logger):
        '''Fields: a list of fields with index, adj (other items are not used)
        builds: index(es) of position, reference allele and alternative alleles. If 
            positions are available, UCSC bins are prepended to the records. If reference
            and alternative alleles are available, the records are processed for correct
            format of ref and alt alleles.
        '''
        self.logger = logger
        self.build = build
        self.raw_fields = fields
        self.fields = []
        self.delimiter = delimiter
        self.columnRange = [None] * len(self.raw_fields)
        self.first_time = True
        self.valid_till = None  # genotype fields might be disabled

    def reset(self, validTill=None):
        self.first_time = True
        self.fields = []
        self.nColumns = 0
        self.valid_till = validTill

    def process(self, line):
        tokens = [x.strip() for x in line.split(self.delimiter)]
        if self.first_time:
            self.nColumns = len(tokens)
            cIdx = 0
            for fIdx, field in enumerate(self.raw_fields):
                if self.valid_till is not None and fIdx >= self.valid_till:
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
                    for x in field.index.split(','):
                        if ':' in x:
                            # a slice
                            if x.count(':') == 1:
                                start,end = map(str.strip, x.split(':'))
                                step = None
                            else:
                                start,end,step = map(str,strip, x.split(':'))
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
                        if len(indexes) == 1:
                            # int, True means 'not a tuple'
                            self.fields.append((indexes[0], True, e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            cIdx += 1
                        else:
                            # a tuple
                            self.fields.append((tuple(indexes), False, e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            cIdx += 1
                    elif len(indexes) == 1:
                        # single slice
                        cols = range(len(tokens))[indexes[0]]
                        for c in cols:
                            self.fields.append((c, True, e))
                        self.columnRange[fIdx] = (cIdx, cIdx + len(cols))
                        cIdx += len(cols)
                    else:
                        # we need to worry about mixing integer and slice
                        indexes = [repeat(s, len(tokens)) if type(s) == int else range(len(tokens))[s] for s in indexes]
                        count = 0
                        for c in izip(*indexes):
                            count += 1
                            self.fields.append((tuple(c), False, e))
                        self.columnRange[fIdx] = (cIdx, cIdx + count)
                        cIdx += count
                except Exception as e:
                    self.logger.debug(e)
                    raise ValueError('Incorrect value adjustment functor or function: {}'.format(field.adj))
            self.first_time = False
        #
        try:
            # we first trust that nothing can go wrong and use a quicker method
            records = [(tokens[col] if t else [tokens[x] for x in col]) if adj is None else \
                (adj(tokens[col]) if t else adj([tokens[x] for x in col])) for col,t,adj in self.fields]
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

class BaseImporter:
    '''A general class for importing variants'''
    def __init__(self, proj, files, build, force, mode='insert'):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        self.mode = mode
        self.sample_in_file = []
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
        #
        if mode == 'insert':
            self.files = []
            cur = self.db.cursor()
            cur.execute('SELECT filename from filename;')
            existing_files = [x[0] for x in cur.fetchall()]
            for filename in files:
                if filename in existing_files:
                    if force:
                        self.logger.info('Re-importing imported file {}'.format(filename))
                        IDs = proj.selectSampleByPhenotype('filename = "{}"'.format(filename))
                        proj.removeSamples(IDs)
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
        else:
            for filename in files:
                if not os.path.isfile(filename):
                    raise ValueError('File {} does not exist'.format(filename))
            self.files = files
        # for #record, #sample variant (new or updated), #variant, new SNV, insertion, deletion, complex variants, invalid record, updated record
        self.count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.total_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.import_alt_build = False
        if len(self.files) == 0:
            raise ValueError('No file to import')
        #
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                self.logger.info('Using primary reference genome {} of the project.'.format(self.build))
            else:
                self.build = self.guessBuild(self.files[0])
                if self.build is None:
                    raise ValueError('Cannot determine a reference genome from files provided. Please specify it using parameter --build')
                else:
                    self.logger.info('Reference genome is determined to be {}'.format(self.build))
        else:
            self.build = build
        #
        if self.proj.build is None:
            if mode == 'insert':
                self.proj.setRefGenome(self.build)
            else:
                raise ValueError('Cannot update variants of a project without variants.')
        elif self.build == self.proj.build:
            # perfect case
            pass
        elif self.build == self.proj.alt_build:
            # troublesome
            self.import_alt_build = True
        elif self.proj.alt_build is None:
            # even more troublesome
            self.logger.warning('The new files uses a different refrence genome ({}) from the primary reference genome ({}) of the project.'.format(self.build, self.proj.build))
            self.logger.info('Adding an alternative reference genome ({}) to the project.'.format(self.build))
            tool = LiftOverTool(self.proj)
            tool.setAltRefGenome(self.build, build_index= (mode == 'insert'))
            self.import_alt_build = True
        #

    def __del__(self):
        if self.mode == 'insert':
            self.proj.createIndexOnMasterVariantTable()

    def guessBuild(self, file):
        # by default, reference genome cannot be determined from file
        return None

    def openFile(self, filename):
        if filename.lower().endswith('.gz'):
            return gzip.open(filename, 'rb')
        else:
            # text file
            # because readline() from gzip.open will be byte, not string, we should return
            # binary here in order to process them equally in order for things to work
            # correctly under python 3 
            return open(filename, 'rb')

    def createLocalVariantIndex(self):
        '''Create index on variant (chr, pos, ref, alt) -> variant_id'''
        self.variantIndex = {}
        cur = self.db.cursor()
        numVariants = self.db.numOfRows('variant')
        if numVariants == 0:
            return
        self.logger.debug('Creating local indexes for {:,} variants'.format(numVariants));
        if self.import_alt_build:
            cur.execute('SELECT variant_id, alt_chr, alt_pos, ref, alt FROM variant;')
        else:
            cur.execute('SELECT variant_id, chr, pos, ref, alt FROM variant;')
        prog = ProgressBar('Getting existing variants', numVariants)
        for count, rec in enumerate(cur):
            # zero for existing loci
            self.variantIndex[(rec[1], rec[2], rec[3], rec[4])] = (rec[0], 0)
            if count % self.db.batch == 0:
                prog.update(count)
        prog.done()

    def getSampleName(self, filename, prober):
        '''Prove text file for sample name'''
        header_line = None
        count = 0
        with self.openFile(filename) as input:
            for line in input:
                line = line.decode()
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
                                if max(cols) - min(cols)  < len(header):
                                    return len(rec), [header[len(header) - prober.nColumns + x] for x in cols]
                                else:
                                    return len(rec), []
                    except Exception as e:
                        # perhaps not start with #, if we have no header, use it anyway
                        if header_line is None:
                            header_line = line
                        count += 1
                        if count == 100:
                            raise ValueError('No genotype column could be determined after 1000 lines.')
                        self.logger.debug(e)

    def recordFileAndSample(self, filename, sampleNames, hasGenotype=True, sampleFields = []):
        cur = self.db.cursor()
        # get header of file
        header = ''
        with self.openFile(filename) as input:
            for line in input:
                line = line.decode()
                if line.startswith('#'):
                    header += line
                else:
                    break
        cur.execute("INSERT INTO filename (filename, header) VALUES ({0}, {0});".format(self.db.PH), (filename, header))
        filenameID = cur.lastrowid
        sample_ids = []
        s = delayedAction(self.logger.info, 'Creating {} sample variant tables'.format(len(sampleNames)))
        for samplename in sampleNames:
            cur.execute('INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'.format(self.db.PH),
                (filenameID, samplename))
            sample_ids.append(cur.lastrowid)
            self.proj.createNewSampleVariantTable('{0}_genotype.genotype_{1}'.format(self.proj.name, cur.lastrowid),
                hasGenotype, sampleFields)
        del s
        return sample_ids

    def importData(self):
        '''Start importing'''
        sample_in_files = []
        for count,f in enumerate(self.files):
            self.logger.info('{} variants from {} ({}/{})'.format('Importing' if self.mode == 'insert' else 'Updating', f, count + 1, len(self.files)))
            self.importFromFile(f)
            if self.mode == 'insert':
                new_var = sum(self.count[3:7])
                self.logger.info('{:,} variants ({:,} new{}) from {:,} records are imported, {}.'\
                    .format(self.count[2], new_var, 
                        ''.join([', {:,} {}'.format(x, y) for x, y in \
                            zip(self.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                        self.count[0],
                        'no sample is created' if len(self.sample_in_file) == 0 else 'with a total of {:,} genotypes from {}'.format(
                            self.count[1], 'sample {}'.format(self.sample_in_file[0]) if len(self.sample_in_file) == 1 else '{:,} samples'.format(len(self.sample_in_file)))))
            else:
                self.logger.info('Field{} {} of {:,} variants{} are updated'.format('' if len(self.variant_info) == 1 else 's', ', '.join(self.variant_info), self.count[8],
                    '' if self.count[1] == 0 else ' and geno fields of {:,} genotypes'.format(self.count[1])))
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
            sample_in_files.extend(self.sample_in_file)
        if len(self.files) > 1:
            if self.mode == 'insert':
                new_var = sum(self.total_count[3:7])
                self.logger.info('{:,} variants ({:,} new{}) from {:,} records are imported, {}.'\
                    .format(self.total_count[2], new_var, 
                        ''.join([', {:,} {}'.format(x, y) for x, y in \
                            zip(self.total_count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'invalid']) if x > 0]),
                        self.total_count[0],
                        'no sample is created' if len(sample_in_files) == 0 else 'with a total of {:,} genotypes from {}'.format(
                            self.total_count[1], 'sample {}'.format(sample_in_files[0]) if len(sample_in_files) == 1 else '{:,} samples'.format(len(sample_in_files)))))
            else:
                self.logger.info('Field{} {} of {:,} variants{} are updated'.format('' if len(self.variant_info) == 1 else 's', ', '.join(self.variant_info), self.total_count[8],
                    '' if self.total_count[1] == 0 else ' and geno fields of {:,} genotypes'.format(self.total_count[1])))
        if self.mode == 'insert' and sum(self.total_count[3:7]) > 0 and self.proj.alt_build is not None:
            coordinates = set([(x[0], x[1]) for x,y in self.variantIndex.iteritems() if y[1] == 1])
            # we need to run lift over to convert coordinates before importing data.
            tool = LiftOverTool(self.proj)
            if self.import_alt_build:
                self.logger.info('Mapping new variants at {} loci from {} to {} reference genome'.format(len(coordinates), self.proj.alt_build, self.proj.build))
                coordinateMap = tool.mapCoordinates(coordinates, self.proj.alt_build, self.proj.build)
                query = 'UPDATE variant SET bin={0}, chr={0}, pos={0} WHERE variant_id={0};'.format(self.db.PH)
            else:
                self.logger.info('Mapping new variants at {} loci from {} to {} reference genome'.format(len(coordinates), self.proj.build, self.proj.alt_build))
                coordinateMap = tool.mapCoordinates(coordinates, self.proj.build, self.proj.alt_build)
                query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(self.db.PH)
                # this should not really happen, but people (like me) might manually mess up with the database
                headers = self.db.getHeaders('variant')
                if not 'alt_pos' in headers:
                    self.logger.info('Adding fields alt_bin, alt_chr and alt_pos to table variant')
                    self.db.execute('ALTER TABLE variant ADD alt_bin INT NULL;')
                    self.db.execute('ALTER TABLE variant ADD alt_chr VARCHAR(20) NULL;')
                    self.db.execute('ALTER TABLE variant ADD alt_pos INT NULL;')
            # update records
            prog = ProgressBar('Updating coordinates', len(self.variantIndex))
            count = 0
            cur = self.db.cursor()
            for k,v in self.variantIndex.iteritems():
                if v[1] == 1:
                    count += 1
                    try:
                        (chr, pos) = coordinateMap[(k[0], k[1])]
                    except Exception as e:
                        self.logger.debug(e)
                        # unmapped
                        continue
                    cur.execute(query, (getMaxUcscBin(pos - 1, pos), chr, pos, v[0]))
                    if count % self.db.batch == 0:
                        self.db.commit()
                        prog.update(count)
            self.db.commit()
            prog.done()
            

class TextImporter(BaseImporter):
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, files, build, format, sample_name=None, 
        force=False, fmt_args=[]):
        BaseImporter.__init__(self, proj, files, build, force, mode='insert')
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        # try to guess file type
        if not format:
            filename = self.files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            fmt = fileFMT(format, fmt_args)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        self.sample_name = sample_name
        #
        # how to split processed records
        self.ranges = fmt.ranges
        self.variant_fields = [x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]]
        self.variant_info = [x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]]
        self.genotype_field = [x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]]
        self.genotype_info = [x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]]
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for pos, ref, alt are 1, 2, 3 in fields.
            self.processor = LineImporter(fmt.fields, [(1, 2, 3)], fmt.delimiter, self.logger)
        else:  # position or range type
            raise ValueError('Can only import data with full variant information (chr, pos, ref, alt)')
        # probe number of sample
        if self.genotype_field:
            self.prober = LineImporter([fmt.fields[fmt.ranges[2]]], [], fmt.delimiter, self.logger)
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
        fbin, fchr, fpos = ('alt_bin', 'alt_chr', 'alt_pos') if self.import_alt_build else ('bin', 'chr', 'pos')
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH)
        self.variant_insert_query = 'INSERT INTO variant ({0}, {1}, {2}, ref, alt {3}) VALUES ({4});'\
            .format(fbin, fchr, fpos, ' '.join([', ' + x for x in self.variant_info]), ', '.join([self.db.PH]*(len(self.variant_info) + 5)))
        #
        self.createLocalVariantIndex()
        # drop index here after all possible exceptions have been raised.
        self.proj.dropIndexOnMasterVariantTable()
                    
    def addVariant(self, cur, rec):
        #
        var_key = tuple(rec[1:5])
        if var_key in self.variantIndex:
            variant_id = self.variantIndex[var_key][0]
            if len(rec) > 5:
                self.count[8] += 1
                cur.execute(self.update_variant_query, rec[5:] + [variant_id])
            return variant_id
        else:
            # new varaint!
            if rec[4] == '-':
                self.count[5] += 1
            elif rec[3] == '-':
                self.count[4] += 1
            elif len(rec[4]) == 1 and len(rec[3]) == 1:
                self.count[3] += 1
            else:
                self.count[6] += 1
            # alt_chr and alt_pos are updated if adding by alternative reference genome
            cur.execute(self.variant_insert_query, rec)
            variant_id = cur.lastrowid
            # one for new variant
            self.variantIndex[var_key] = (variant_id, 1)
            return variant_id

    def importFromFile(self, input_filename):
        '''Import a TSV file to sample_variant'''
        # reset text processor to allow the input of files with different number of columns
        self.processor.reset()
        if self.genotype_field:
            self.prober.reset()
        #
        if not self.sample_name:
            # if no sample name is specified
            if not self.genotype_field:
                self.logger.warning('Sample information is not recorded for a file without genotype and sample name.')
                sample_ids = []
                self.sample_in_file = []
            else:
                try:
                    numSample, names = self.getSampleName(input_filename, self.prober)
                    if not names:
                        if numSample == 1:
                            self.logger.debug('Missing sample name (name None is used)'.format(numSample))
                            sample_ids = self.recordFileAndSample(input_filename, [None], True,
                                self.genotype_info)
                            self.sample_in_file = [None]
                        elif numSample == 0:
                            self.logger.debug('No genotype column exists in the input file so no sample will be recorded.')
                            sample_ids = []
                            self.sample_in_file = []
                        else:
                            raise ValueError('Failed to guess sample name. Please specify sample names for {} samples using parameter --sample_name, or add a proper header to your input file. See "vtools import_variants -h" for details.'.format(numSample))
                    else:
                        sample_ids = self.recordFileAndSample(input_filename, names, True,
                            self.genotype_info)
                        self.sample_in_file = [x for x in names]
                except ValueError:
                    # cannot find any genotype column, perhaps no genotype is defined in the file (which is allowed)
                    self.logger.warning('No genotype column could be found from the input file. Assuming no genotype.')
                    self.genotype_field = []
                    sample_ids = []
                    self.sample_in_file = []
        else:
            self.sample_in_file = [x for x in self.sample_name]
            if not self.genotype_field:
                # if no genotype, but a sample name is given
                self.logger.debug('Input file does not contain any genotype. Only the variant ownership information is recorded.')
                sample_ids = self.recordFileAndSample(input_filename, self.sample_name, False, self.genotype_info)
            else:
                try:
                    numSample, names = self.getSampleName(input_filename, self.prober)
                    if len(self.sample_name) != numSample:
                        raise ValueError('{} sample detected but only {} names are specified'.format(numSample, len(self.sample_name)))                        
                except ValueError:
                    self.logger.warning('No genotype column could be found from the input file. Assuming no genotype.')
                    self.genotype_field = []
                    self.genotype_info = []
                    # remove genotype field from processor
                    self.processor.reset(validTill=self.ranges[2])
                sample_ids = self.recordFileAndSample(input_filename, self.sample_name, len(self.genotype_field) > 0, self.genotype_info)
        #
        cur = self.db.cursor()
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        if sample_ids:
            genotype_insert_query = {id: 'INSERT INTO {0}_genotype.genotype_{1} VALUES ({2});'\
                .format(self.proj.name, id, ','.join([self.db.PH] * (1 + len(self.genotype_field) + len(self.genotype_info))))
                for id in sample_ids}
        fld_cols = None
        # cache genotype status
        if len(self.genotype_field) > 0:
            # has genotype
            genotype_status = 1
        elif len(sample_ids) > 0:
            # no genotype but with sample
            genotype_status = 2
        else:
            # no genotype no sample
            genotype_status = 0
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    line = line.decode()
                    if line.startswith('#'):
                        continue
                    for bins, rec in self.processor.process(line):
                        self.count[2] += 1
                        variant_id = self.addVariant(cur, bins + rec[0:self.ranges[2]])
                        if genotype_status == 1:
                            if fld_cols is None:
                                col_rngs = [self.processor.columnRange[x] for x in range(self.ranges[2], self.ranges[4])]
                                fld_cols = []
                                for idx in range(len(sample_ids)):
                                    fld_cols.append([sc + (0 if sc + 1 == ec else idx) for sc,ec in col_rngs])
                                if col_rngs[0][1] - col_rngs[0][0] != len(sample_ids):
                                    self.logger.error('Number of genotypes ({}) does not match number of samples ({})'.format(
                                        col_rngs[0][1] - col_rngs[0][0], len(sample_ids)))
                            for idx, id in enumerate(sample_ids):
                                if rec[self.ranges[2] + idx] is not None:
                                    self.count[1] += 1
                                    cur.execute(genotype_insert_query[id], [variant_id] + [rec[c] for c in fld_cols[idx]])
                        elif genotype_status == 2:
                            # should have only one sample
                            for id in sample_ids:
                                cur.execute(genotype_insert_query[id], [variant_id])
                    self.count[0] += 1
                except Exception as e:
                    self.logger.debug('Failed to process line "{}...": {}'.format(line[:20].strip(), e))
                    self.count[7] += 1
                if self.count[0] % self.db.batch == 0:
                    self.db.commit()
                    prog.update(self.count[0])
            self.db.commit()
            prog.done()

class TextUpdater(BaseImporter):
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, table, files, build, format, fmt_args=[]):
        # if update is None, recreate index
        BaseImporter.__init__(self, proj, files, build, True, mode='update')
        #
        if not proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        # try to guess file type
        if not format:
            filename = self.files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            fmt = fileFMT(format, fmt_args)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        # how to split processed records
        self.ranges = fmt.ranges
        self.variant_fields = [x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]]
        self.variant_info = [x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]]
        self.genotype_field = [x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]]
        self.genotype_info = [x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]]
        #
        if not self.variant_info and not self.genotype_info:
            raise ValueError('No variant or genotype info needs to be updated')
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for pos, ref, alt are 1, 2, 3 in fields.
            self.processor = LineImporter(fmt.fields, [(1, 2, 3)], fmt.delimiter, self.logger)
        else:  # position or range type
            self.processor = LineImporter(fmt.fields, [(1,)], fmt.delimiter, self.logger)
        # probe number of sample
        if self.genotype_field and self.genotype_info:
            self.prober = LineImporter([fmt.fields[fmt.ranges[2]]], [], fmt.delimiter, self.logger)
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
        #if len(self.variant_info) == 0 and len(self.genotype_info == 0:
        #    raise ValueError('No field could be updated using this input file')
        #
        self.input_type = fmt.input_type
        fbin, fchr, fpos = ('alt_bin', 'alt_chr', 'alt_pos') if self.import_alt_build else ('bin', 'chr', 'pos')
        from_table = 'AND variant.variant_id IN (SELECT variant_id FROM {})'.format(table) if table != 'variant' else ''
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1} {2};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH, from_table)
        self.update_position_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} = {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        self.update_range_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} >= {0} AND variant.{4} <= {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        #
        self.createLocalVariantIndex()

    def updateVariant(self, cur, bins, rec):
        if self.input_type == 'variant':
            var_key = tuple(rec[0:4])
            if var_key in self.variantIndex:
                variant_id = self.variantIndex[var_key][0]
                # update by variant_id, do not need bins
                if len(rec) > 4:
                    cur.execute(self.update_variant_query, rec[4:] + [variant_id])
                    self.count[8] += cur.rowcount
                return variant_id
        elif self.input_type == 'position':
            cur.execute(self.update_position_query, rec[2:] + bins + [rec[0], rec[1]])
            self.count[8] += cur.rowcount
        else:  # range based
            cur.execute(self.update_range_query, rec[3:] + bins + [rec[0], rec[1], rec[2]])
            self.count[8] += cur.rowcount
        return None

    def getSampleIDs(self, filename):
        if not self.genotype_field:
            # no genotype_field, good, do not have to worry about genotype
            return []
        # has the file been imported before?
        cur = self.db.cursor()
        cur.execute('SELECT filename from filename;')
        existing_files = [x[0] for x in cur.fetchall()]
        if filename not in existing_files:
            return []
        #
        # what are the samples related to this file?
        cur.execute('SELECT sample_id, sample_name FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id WHERE filename.filename = {}'\
            .format(self.db.PH), (filename,))
        sample_ids = []
        sample_names = []
        for rec in cur:
            sample_ids.append(rec[0])
            sample_names.append(rec[1])
        # what is the sample names get from this file?
        nSample, names = self.getSampleName(filename, self.prober)
        if nSample != len(sample_ids):
            self.logger.warning('Number of samples mismatch. Cannot update genotype')
            return []
        if nSample == 1:
            # if only one sample, update it regardless of sample name.
            return sample_ids
        if sample_names == names:
            # if sample name matches, get sample_ids
            return sample_ids
        else:
            self.logger.warning('Sample names mismatch. Cannot update genotype.')
            return []
        
    def importFromFile(self, input_filename):
        '''Import a TSV file to sample_variant'''
        self.processor.reset()
        if self.genotype_field and self.genotype_info:
            self.prober.reset()
        #
        # do we handle genotype info as well?
        sample_ids = self.getSampleIDs(input_filename) if self.genotype_info else []
        #
        # do we need to add extra columns to the genotype tables
        if sample_ids:
            s = delayedAction(self.logger.info, 'Preparing genotype tables (adding needed fields)...')
            cur = self.db.cursor()
            for id in sample_ids:
                headers = [x.upper() for x in self.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                if 'GT' not in headers:  # for genotype
                    self.logger.debug('Adding column GT to table genotype_{}'.format(id))
                    cur.execute('ALTER TABLE {}_genotype.genotype_{} ADD {} {};'.format(self.proj.name, id, 'GT', 'INT'))
                for field in self.genotype_info:
                    if field.name.upper() not in headers:
                        self.logger.debug('Adding column {} to table genotype_{}'.format(field, id))
                        cur.execute('ALTER TABLE {}_genotype.genotype_{} ADD {} {};'.format(self.proj.name, id, field.name, field.type))
            del s
        else:
            # do not import genotype even if the input file has them
            self.genotype_field = []
            self.genotype_info = []
            self.processor.reset(validTill=self.ranges[2])
        #
        cur = self.db.cursor()
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        if sample_ids:
            genotype_update_query = {id: 'UPDATE {0}_genotype.genotype_{1} SET {2} WHERE variant_id = {3};'\
                .format(self.proj.name, id,
                ', '.join(['{}={}'.format(x, self.db.PH) for x in [y.name for y in self.genotype_info]]),
                self.db.PH)
                for id in sample_ids}
        fld_cols = None
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    line = line.decode()
                    if line.startswith('#'):
                        continue
                    for bins, rec in self.processor.process(line):
                        variant_id = self.updateVariant(cur, bins, rec[0:self.ranges[2]])
                        # variant might not exist
                        if variant_id is not None and sample_ids:
                            if fld_cols is None:
                                col_rngs = [self.processor.columnRange[x] for x in range(self.ranges[3], self.ranges[4])]
                                fld_cols = []
                                for idx in range(len(sample_ids)):
                                    fld_cols.append([sc + (0 if sc + 1 == ec else idx) for sc,ec in col_rngs])
                                if col_rngs[0][1] - col_rngs[0][0] != len(sample_ids):
                                    self.logger.error('Number of genotypes ({}) does not match number of samples ({})'.format(
                                        col_rngs[0][1] - col_rngs[0][0], len(sample_ids)))
                            for idx, id in enumerate(sample_ids):
                                if rec[self.ranges[2] + idx] is not None:
                                    cur.execute(genotype_update_query[id], [rec[c] for c in fld_cols[idx]] + [variant_id])
                                    self.count[1] += 1
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                    self.count[7] += 1
                if self.count[0] % self.db.batch == 0:
                    self.db.commit()
                    prog.update(self.count[0])
            self.db.commit()
            prog.done()

#
#
# Functions provided by this script
#
#

def importVariantsArguments(parser):
    parser.add_argument('input_files', nargs='+',
        help='''A list of files that will be imported. The file should be delimiter
            separated with format described by parameter --format. Gzipped files are
            acceptable.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, it is assumed to be the primary reference genome of the project.
            If a reference genome that is different from the primary reference genome of the
            project is specified, it will become the alternative reference genome of the
            project. The UCSC liftover tool will be automatically called to map input
            coordinates between the primary and alternative reference genomes.''')
    parser.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as VCF (c.f. 'vtools show formats'), or a 
            local format specification file (with extension .fmt). If unspecified,
            variant tools will try to guess format from file extension. Some file
            formats accept parameters (c.f. 'vtools show format FMT') and allow you
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
            not remove wrongfully imported variants from the master variant table.''')

def importVariants(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = TextImporter(proj=proj, files=args.input_files,
                build=args.build, format=args.format, sample_name=args.sample_name,
                force=args.force, fmt_args=args.unknown_args)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)


def updateArguments(parser):
    parser.add_argument('table', help='''variants to be updated.''')
    parser.add_argument('input_files', nargs='+',
        help='''A list of files that will be used to add or update existing fields of
            variants. The file should be delimiter separated with format described by
            parameter --format. Gzipped files are acceptable. If input files contains
            genotype information, have been inputted before, and can be linked to the
            samples they created without ambiguity (e.g. single sample, or samples with
            detectable sample names), genotypes and their information will also be
            updated.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input files,
            which should be the primary (used by default) or alternative (if available)
            reference genome of the project. An alternative reference genome will be
            added to the project if needed.''')
    parser.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as ANNOVAR_mut_type (c.f. 'vtools show formats'),
            or a local format specification file (with extension .fmt). Some formats 
            accept parameters (c.f. 'vtools show format FMT') and allow you to update
            additional or alternative fields from the input file.''')

def update(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = TextUpdater(proj=proj, table=args.table, files=args.input_files,
                build=args.build, format=args.format, fmt_args=args.unknown_args)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

