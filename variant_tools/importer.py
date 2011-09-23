#!/usr/bin/env python
#
# $File: importer.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
        return item[1].split(self.sep)[self.idx] if self.idx is not None else self.default

class VcfGenotype:
    def __init__(self, default=None):
        '''Define an extractor that extract genotype from a .vcf file'''
        self.default = default
        self.map = {'0/0': default, '0|0': default, '0': default,
            '0/1': ('1',), '1/0': ('1',), '0|1': ('1',), '1|0': ('1',),
            '1/1': ('2',), '1|1': ('2',),
            '0/2': ('0', '1'), '2/0': ('0', '1'), '0|2': ('0', '1'), '2|0': ('0', '1'), 
            '1/2': ('-1', '-1'), '2/1': ('-1', '-1'), '1|2': ('-1', '-1'), '2|1': ('-1', '-1'),
            '2/2': ('0', '2'), '2|2': ('0', '2'),
            '0': default, '1': ('1',)}

    def __call__(self, item):
        # the most common and correct case...
        return self.map[item.partition(':')[0]]

class VcfGenoFromFormat:
    def __init__(self, default=None):
        '''Define an extractor that return genotype according to a format string.
        This is used to extract genotype from the format string of vcf files.
        '''
        self.fmt = '\t'
        self.idx = None
        self.default = default
        self.map = {'0/0': default, '0|0': default, '0': default,
            '0/1': ('1',), '1/0': ('1',), '0|1': ('1',), '1|0': ('1',),
            '1/1': ('2',), '1|1': ('2',),
            '0/2': ('0', '1'), '2/0': ('0', '1'), '0|2': ('0', '1'), '2|0': ('0', '1'), 
            '1/2': ('-1', '-1'), '2/1': ('-1', '-1'), '1|2': ('-1', '-1'), '2|1': ('-1', '-1'),
            '2/2': ('0', '2'), '2|2': ('0', '2'),
            '0': default, '1': ('1',)}

    def __call__(self, item):
        # the most common and correct case...
        if item[0][:2] == 'GT':
            return self.map[item[1].partition(':')[0]]
        elif item[0] != self.fmt:
            fmt, val = item
            self.fmt = fmt
            fields = fmt.split(self.sep)
            if self.name in fields:
                self.idx = fields.index(self.name)
                return self.map[val.split(self.sep)[self.idx]]
            else:
                self.idx = None
                return self.default
        return self.map[item[1].split(self.sep, self.idx + 1)[self.idx]] if self.idx is not None else self.default

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
            if type(item) == tuple:
                if type(item[0]) == tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item


class TextProcessor:
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

    def process(self, line):
        tokens = [x.strip() for x in line.split(self.delimiter)]
        if self.first_time:
            cIdx = 0
            for fIdx, field in enumerate(self.raw_fields):
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
                            # int
                            self.fields.append((indexes[0], e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            cIdx += 1
                        else:
                            # a tuple
                            self.fields.append((tuple(indexes), e))
                            self.columnRange[fIdx] = (cIdx, cIdx+1)
                            cIdx += 1
                    elif len(indexes) == 1:
                        # single slice
                        cols = range(len(tokens))[indexes[0]]
                        for c in cols:
                            self.fields.append((c, e))
                        self.columnRange[fIdx] = (cIdx, cIdx + len(cols))
                        cIdx += len(cols)
                    else:
                        # we need to worry about mixing integer and slice
                        indexes = [repeat(s, len(tokens)) if type(s) == int else range(len(tokens))[s] for s in indexes]
                        count = 0
                        for c in izip(*indexes):
                            count += 1
                            self.fields.append((tuple(c), e))
                        self.columnRange[fIdx] = (cIdx, cIdx + count)
                        cIdx += count
                except Exception as e:
                    self.logger.debug(e)
                    raise ValueError('Incorrect value adjustment functor or function: {}'.format(field.adj))
            self.first_time = False
        #        
        num_records = 1
        records = [None]*len(self.fields)
        for idx, (col, adj) in enumerate(self.fields):
            item = tokens[col] if type(col) == int else [tokens[x] for x in col]
            if adj is not None:
                try:
                    item = adj(item)
                    if type(item) == tuple and len(item) > num_records:
                        num_records = len(item)
                except Exception as e:
                    self.logger.debug(e)
                    # missing ....
                    item = None
            records[idx] = item
        # handle records
        if not self.build:
            # there is no build information, this is 'field' annotation, nothing to worry about
            if num_records == 1:
                yield [], [x[0] if type(x) == tuple else x for x in records]
            else:
                for i in range(num_records):
                    yield [], [(x[i] if i < len(x) else None) if type(x) == tuple else x for x in records]
        elif len(self.build[0]) == 1:
            for i in range(num_records):
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) == tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None) if type(x) == tuple else x for x in records]
                bins = [getMaxUcscBin(int(rec[pos_idx]) - 1, int(rec[pos_idx])) if rec[pos_idx] else None for pos_idx, in self.build]
                yield bins, rec
        else:
            for i in range(num_records):
                bins = []
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) == tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None) if type(x) == tuple else x for x in records]
                for pos_idx, ref_idx, alt_idx in self.build:
                    bin, pos, ref, alt = normalizeVariant(int(rec[pos_idx]) if rec[pos_idx] else None, rec[ref_idx], rec[alt_idx])
                    bins.append(bin)
                    rec[pos_idx] = pos
                    rec[ref_idx] = ref
                    rec[alt_idx] = alt
                yield bins, rec



class Importer:
    '''A general class for importing variants'''
    def __init__(self, proj, files, build, force, recreateIndex=True):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
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
                    for ID in IDs:
                        proj.removeSample(ID)
                    # remove file record
                    cur = self.db.cursor()
                    cur.execute('DELETE FROM filename WHERE filename={};'.format(self.db.PH), (filename,))
                    self.db.commit()
                    self.files.append(filename)
                else:
                    self.logger.info('Ignoring imported file {}'.format(filename))
            else:
                self.files.append(filename)
        # for all record, new SNV, insertion, deletion, complex variants, invalid record, and updated record
        self.count = [0, 0, 0, 0, 0, 0, 0]
        self.total_count = [0, 0, 0, 0, 0, 0, 0]
        self.import_alt_build = False
        self.recreateIndex = False
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
            self.proj.setRefGenome(self.build)
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
            # we will drop indexes soon so do not build index
            tool.setAltRefGenome(self.build, build_index=False)
            self.import_alt_build = True
        else:
            raise ValueError('Specified build {} does not match either the primary '.format(self.build) + \
                ' {} or the alternative reference genome of the project.'.format(self.proj.build, self.proj.alt_build))
        #
        self.recreateIndex = recreateIndex
        if self.recreateIndex:
            self.proj.dropIndexOnMasterVariantTable()
        #
        self.createLocalVariantIndex()

    def __del__(self):
        if self.recreateIndex:
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

    def recordFileAndSample(self, filename, sampleNames, hasGenotype=True, sampleFields = []):
        cur = self.db.cursor()
        cur.execute("INSERT INTO filename (filename) VALUES ({});".format(self.db.PH), (filename,))
        filenameID = cur.lastrowid
        sample_ids = []
        s = delayedAction(self.logger.info, 'Creating {} sample variant tables'.format(len(sampleNames)))
        for samplename in sampleNames:
            cur.execute('INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'.format(self.db.PH),
                (filenameID, samplename))
            sample_ids.append(cur.lastrowid)
            self.proj.createNewSampleVariantTable('{0}_genotype.sample_variant_{1}'.format(self.proj.name, cur.lastrowid),
                hasGenotype, sampleFields)
        del s
        return sample_ids
        
    def importData(self):
        '''Start importing'''
        for count,f in enumerate(self.files):
            self.logger.info('Importing genotype from {} ({}/{})'.format(f, count + 1, len(self.files)))
            self.importFromFile(f)
            self.logger.info('{:,} new variants from {:,} records are imported, with {:,} SNVs, {:,} insertions, {:,} deletions, and {:,} complex variants.{}'\
                .format(sum(self.count[1:5]), self.count[0], self.count[1], self.count[2], self.count[3], self.count[4],
                ' {} invalid records are ignored'.format(self.count[5]) if self.count[5] > 0 else ''))
            if self.count[6] > 0:
                self.logger.info('{:,} exiting variants are updated'.format(self.count[6]))
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
        if len(self.files) > 1:
            self.logger.info('{:,} new variants from {:,} records in {} files are imported, with {:,} SNVs, {:,} insertions, {:,} deletions, and {:,} complex variants.{}'\
                .format(sum(self.total_count[1:5]), self.total_count[0], len(self.files), self.total_count[1], self.total_count[2], self.total_count[3], self.total_count[4],
                ' {} invalid records are ignored'.format(self.total_count[5]) if self.total_count[5] > 0 else ''))
            if self.total_count[6] > 0:
                self.logger.info('{:,} exiting variants are updated'.format(self.total_count[6]))
        if sum(self.total_count[1:5]) > 0 and self.proj.alt_build is not None:
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
            prog = ProgressBar('Updating coordinates', self.total_count[0])
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
            


class vcfImporter(Importer):
    '''A vcf importer to import genotype from one or more vcf files.
    In case of vcf file, it records the type of variant for this sample, 
    which can be 1 for genotype 0/1, heterozygous, 2 for genotype 1/1,
    homozygous of alternative alleles, and -1 for gentoype 1/2, which 
    consists of two different alternative alleles. In the last case,
    an individual will have two variants with different alternative
    alleles, each with a type -1. That it to say, allele frequency should
    be calculated as sum (abs(type))/ (2*num_of_sample). '''
    def __init__(self, proj, files, build=None, variant_only=False, info=[], force=False):
        '''see importVariant.py -h for details about parameters. Additional
        keyword paramters such as user, passwd and host are passed to
        MySQLdb.connect.
        '''
        Importer.__init__(self, proj, files, build, force)
        # vcf tools only support DP for now
        self.variant_only = variant_only
        self.import_depth = 'DP' in info
        # importing data to alternative reference genome
        if self.import_alt_build:
            self.variant_insert_query = 'INSERT INTO variant (alt_bin, alt_chr, alt_pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        else:
            self.variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)

    def guessBuild(self, filename):
        '''Called by the initializer to determine reference genome
        '''
        with self.openFile(filename) as input:
            for line in input:
                line = line.decode()
                if line.startswith('##reference'):
                    if 'NCBI36' in line.upper() or 'HG18' in line.upper() or 'HUMAN_B36' in line.upper():
                        return 'hg18'
                    elif 'GRCH37' in line.upper() or 'HG19' in line.upper() or 'HUMAN_B37' in line.upper():
                        return 'hg19'
                if not line.startswith('#'):
                    return None

    def getMetaInfo(self, filename):
        '''Probe vcf files for additional information to be put to the variant_meta table.
        '''
        samples = []
        with self.openFile(filename) as input:
            line = input.readline().decode()
            if not line.startswith('##fileformat=VCF'):
                self.logger.error('Invalid vcf file. file not started with line ##fileformat')
                raise ValueError('Invalid vcf file')
            if not line.strip().endswith('4.0') and not line.strip().endswith('4.1'):
                raise ValueError('This importer tool only supports VCF format v4.0 and v4.1. \
                    Please use vcftools to convert your vcf file to a supported format')
            #
            for line in input:
                line = line.decode()
                if line.startswith('#CHR'):
                    samples = line.split()[9:]
                if not line.startswith('#'):
                    break
        return samples

    def addVariant(self, cur, chr, pos, ref, alt):
        # if chr, pos are from alternative reference genome
        bin, pos, ref, alt = normalizeVariant(pos, ref, alt)
        var_key = (chr, pos, ref, alt)
        #
        try:
            return self.variantIndex[var_key][0]
        except:
            # new varaint!
            if alt == '-':
                self.count[3] += 1
            elif ref == '-':
                self.count[2] += 1
            elif len(alt) == 1 and len(ref) == 1:
                self.count[1] += 1
            else:
                self.count[4] += 1
            # alt_chr and alt_pos are updated if adding by alternative reference genome
            cur.execute(self.variant_insert_query, (bin, chr, pos, ref, alt))
            variant_id = cur.lastrowid
            # one for new variant
            self.variantIndex[var_key] = (variant_id, 1)
            return variant_id

    def importFromFile(self, input_filename):
        '''Import a VCF file to sample_variant'''
        #
        # handle meta information and get sample names
        sampleNames = self.getMetaInfo(input_filename)
        if len(sampleNames) == 0:
            self.logger.warning('No sample is created for this file because it does not contain any genotype.')
        elif self.variant_only:
            self.logger.warning('No sample is created for this file because no genotype is imported.')
        #
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        no_sample = self.variant_only or len(sampleNames) == 0
        if not no_sample:
            sample_ids = self.recordFileAndSample(input_filename, [] if no_sample else sampleNames, True,
                ['DP'] if self.import_depth else [])   # record individual depth, total depth is divided by number of sample in a file
            #
            nSample = len(sample_ids)
        #
        DP_pattern = re.compile('.*DP=(\d+)')
        #
        cur = self.db.cursor()
        # sample variants are inserted into different tables in a separate database.
        if not no_sample:
            genotype_insert_query = {x: 'INSERT INTO {1}_genotype.sample_variant_{3} VALUES ({0}, {0} {2});'\
                .format(self.db.PH, self.proj.name, ',' + self.db.PH if self.import_depth else '', x) for x in sample_ids}
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    line = line.decode()
                    if line.startswith('#'):
                        continue
                    self.count[0] += 1
                    # get data
                    tokens = [x.strip() for x in line.split('\t')]
                    chr = tokens[0][3:] if tokens[0].startswith('chr') else tokens[0]
                    pos = int(tokens[1])
                    ref = tokens[3]
                    # we only extract info get depth.
                    if not no_sample and self.import_depth:
                        m = DP_pattern.match(tokens[7])
                        DP = [None if m is None else float(m.group(1))/nSample]
                    else:
                        DP = []
                    # is GT the first field?
                    try:
                        # the format field of a vcf can be empty ('.', does not contain GT) if there is no sample
                        GT_idx = 0 if no_sample or tokens[8].startswith('GT') else tokens[8].split(':').index('GT')
                    except Exception as e:
                        self.logger.debug(e)
                        raise ValueError('The genotype format field does not have GT')
                    # for efficiency, we separte out this most common case ...
                    if len(ref) == 1 and len(tokens[4]) == 1:   
                        # the easy case: there is only one alternative allele,
                        # all genotypes should be 0/0, 0/1, 1/0, or 1/1.
                        alt = tokens[4][0]
                        variant_id = self.addVariant(cur, chr, pos, ref, alt)
                        #
                        if not no_sample:
                            variants = [x.split(':')[GT_idx].count('1') for x in tokens[-len(sample_ids):]]
                            for var_idx, var in enumerate(variants):
                                if var != 0:  # genotype 0|0 are ignored
                                    cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id, var] + DP)
                    else:
                        # now, this is the common case with insertion, deletion, and multiple alternative variants
                        alts = tokens[4].split(',')
                        variant_id = [0] * len(alts)
                        for altidx, alt in enumerate(alts):
                            variant_id[altidx] = self.addVariant(cur, chr, pos, ref, alt)
                        if not no_sample:
                            # process variants
                            for var_idx, var in enumerate([x.split(':')[GT_idx] for x in tokens[-len(sample_ids):]]):
                                if len(var) == 3:  # regular
                                    gt = var[0] + var[2]  # GT can be separated by / or |
                                    if gt in ['01', '10']:
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif gt in ['02', '20']:
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
                                    elif gt == '11':
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[0], 2] + DP)
                                    elif gt in ['12', '21']:
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[0], -1] + DP)
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[1], -1] + DP)
                                    elif gt == '22':
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[1], 2] + DP)
                                    elif gt == '00':
                                        pass
                                    else:
                                        raise ValueError('I do not know how to process genotype {}'.format(var))
                                else: # should have length 1
                                    if var == '1':
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif var == '2':
                                        cur.execute(genotype_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
                                    elif var == '0':
                                        pass
                                    else:
                                        raise ValueError('I do not know how to process genotype {}'.format(var))
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                    self.count[5] += 1
                if self.count[0] % self.db.batch == 0:
                    self.db.commit()
                    prog.update(self.count[0])
            self.db.commit()
            prog.done()


class txtImporter(Importer):
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, files, build, format, sample_name=None, variant_info=None,
        genotype_info=[], update=None, force=False):
        # if update is None, recreate index
        Importer.__init__(self, proj, files, build, force, update is None)
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
            fmt = fileFMT(format, variant_info=None, genotype_info=None)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Input file format {} is not currently supported by variant tools'.format(format))
        #
        self.sample_name = sample_name
        #
        # how to split processed records
        self.ranges = fmt.ranges
        self.variant_fields = [x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]]
        self.variant_info = [x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]]
        self.genotype_field = [x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]]
        self.genotype_info = [x.name for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]]
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for pos, ref, alt are 1, 2, 3 in fields.
            self.processor = TextProcessor(fmt.fields, [(1, 2, 3)], fmt.delimiter, self.logger)
        else:  # position or range type
            self.processor = TextProcessor(fmt.fields, [(1,)], fmt.delimiter, self.logger)
        # probe number of sample
        if self.genotype_field:
            self.prober = TextProcessor([fmt.fields[fmt.ranges[2]]], [], fmt.delimiter, self.logger)
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
        self.update = update
        if fmt.input_type != 'variant':
            if not self.update:
                self.logger.info('Parameter --update variant is assumed for input file with type {}'.format(fmt.input_type))
                self.update = 'variant'
        if self.update and len(self.variant_info) == 0:
            raise ValueError('No field could be updated using this input file')
        #
        self.input_type = fmt.input_type
        fbin, fchr, fpos = ('alt_bin', 'alt_chr', 'alt_pos') if self.import_alt_build else ('bin', 'chr', 'pos')
        from_table = 'AND variant.variant_id IN (SELECT variant_id FROM {})'.format(self.update) if self.update != 'variant' else ''
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1} {2};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH, from_table)
        self.update_position_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} = {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        self.update_range_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} >= {0} AND variant.{4} <= {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        self.variant_insert_query = 'INSERT INTO variant ({0}, {1}, {2}, ref, alt {3}) VALUES ({4});'\
            .format(fbin, fchr, fpos, ' '.join([', ' + x for x in self.variant_info]), ', '.join([self.db.PH]*(len(self.variant_info) + 5)))

    def getSampleName(self, filename):
        '''Prove text file for sample name'''
        header = None
        count = 0
        with self.openFile(filename) as input:
            for line in input:
                line = line.decode()
                # the last # line
                if line.startswith('#'):
                    header = line
                else:
                    try:
                        for bins, rec in self.prober.process(line):
                            if header is None:
                                return len(rec), []
                            elif len(rec) == 0:
                                return 0, []
                            else:
                                self.logger.info('HERE')
                                cols = [x[0] for x in self.prober.fields]
                                if type(cols[0]) == tuple:
                                    fixed = False
                                    # mutiple ones, need to figure out the moving one
                                    for i,idx in enumerate(self.prober.raw_fields[0].index.split(',')):
                                        if ':' in idx:
                                            cols = [x[i] for x in cols]
                                            fixed = True
                                            break
                                    if not fixed:
                                        cols = [x[-1] for x in cols]
                                header = header.split(self.prober.delimiter)
                                if max(cols) < len(header):
                                    return len(rec), [header[x] for x in cols]
                                else:
                                    return len(rec), []
                    except Exception as e:
                        count += 1
                        if count == 100:
                            raise ValueError('Cannot determine header of file')
                        self.logger.debug(e)
                    
    def addVariant(self, cur, rec):
        #
        var_key = tuple(rec[1:5])
        if var_key in self.variantIndex:
            variant_id = self.variantIndex[var_key][0]
            if len(rec) > 5:
                self.count[6] += 1
                cur.execute(self.update_variant_query, rec[5:] + [variant_id])
            return variant_id
        else:
            # new varaint!
            if rec[4] == '-':
                self.count[3] += 1
            elif rec[3] == '-':
                self.count[2] += 1
            elif len(rec[4]) == 1 and len(rec[3]) == 1:
                self.count[1] += 1
            else:
                self.count[4] += 1
            # alt_chr and alt_pos are updated if adding by alternative reference genome
            cur.execute(self.variant_insert_query, rec)
            variant_id = cur.lastrowid
            # one for new variant
            self.variantIndex[var_key] = (variant_id, 1)
            return variant_id

    def updateVariant(self, cur, bins, rec):
        if self.input_type == 'variant':
            var_key = tuple(rec[0:4])
            if var_key in self.variantIndex:
                variant_id = self.variantIndex[var_key][0]
                # update by variant_id, do not need bins
                cur.execute(self.update_variant_query, rec[4:] + [variant_id])
                self.count[6] += cur.rowcount
        elif self.input_type == 'position':
            cur.execute(self.update_position_query, rec[2:] + bins + [rec[0], rec[1]])
            self.count[6] += cur.rowcount
        else:  # range based
            cur.execute(self.update_range_query, rec[3:] + bins + [rec[0], rec[1], rec[2]])
            self.count[6] += cur.rowcount

    def importFromFile(self, input_filename):
        '''Import a TSV file to sample_variant'''
        if self.update:
            sample_ids = []
        else:
            if not self.sample_name:
                # if no sample name is specified
                if not self.genotype_field:
                    self.logger.warning('Sample information is not recorded for a file without genotype and sample name.')
                    sample_ids = []
                else:
                    numSample, names = self.getSampleName(input_filename)
                    if not names:
                        if numSample == 1:
                            self.logger.warning('Missing sample name (name None is used)'.format(numSample))
                            sample_ids = self.recordFileAndSample(input_filename, [None], True,
                                self.genotype_info)
                        elif numSample == 0:
                            self.logger.info('No genotype column exists in the input file so no sample will be recorded.')
                            sample_ids = []
                        else:
                            raise ValueError('Failed to guess sample name. Please specify sample names for {} samples using parameter --sample_name'.format(numSample))
                    else:
                        self.logger.info('Using sample name {}'.format(', '.join(names)))
                        sample_ids = self.recordFileAndSample(input_filename, names, True,
                            self.genotype_info)
            else:
                if not self.genotype_field:
                    # if no genotype, but a sample name is given
                    self.logger.info('Input file does not contain any genotype. Only the variant ownership information is recorded.')
                    sample_ids = self.recordFileAndSample(input_filename, self.sample_name, False, self.genotype_info)
                else:
                    numSample, names = self.getSampleName(input_filename)
                    if len(self.sample_name) != numSample:
                        raise ValueError('{} sample detected but only {} names are specified'.format(numSample, len(self.sample_name)))                        
                    sample_ids = self.recordFileAndSample(input_filename, self.sample_name, True, self.genotype_info)
        #
        cur = self.db.cursor()
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        if sample_ids:
            genotype_insert_query = {id: 'INSERT INTO {0}_genotype.sample_variant_{1} VALUES ({2});'\
                .format(self.proj.name, id, ','.join([self.db.PH] * (1 + len(self.genotype_field) + len(self.genotype_info))))
                for id in sample_ids}
        rngs = None
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    line = line.decode()
                    if line.startswith('#'):
                        continue
                    for bins, rec in self.processor.process(line):
                        if self.update:
                            self.updateVariant(cur, bins, rec[0:self.ranges[2]])
                        else:
                            variant_id = self.addVariant(cur, bins + rec[0:self.ranges[2]])
                            if not rngs:
                                rngs = [self.processor.columnRange[x] for x in range(self.ranges[2], self.ranges[4])]
                                if rngs[0][1] - rngs[0][0] != len(sample_ids):
                                    raise ValueError('Number of genotypes ({}) does not match number of samples ({})'.format(rngs[0][1] - rngs[0][0], len(sample_ids)))
                            for idx, id in enumerate(sample_ids):
                                if rec[self.ranges[2] + idx]:
                                    cur.execute(genotype_insert_query[id], [variant_id] + [rec[sc + (0 if sc + 1 == ec else idx)] for sc,ec in rngs])
                            self.count[0] += 1
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                    self.count[5] += 1
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

def importVCFArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''A list of files that will be imported. The file should be in 
            VCF 4.0 format and can be compressed in gzip format.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, variant tools will try to determine the reference genome from the
            header of VCF files, and use the primary reference genome of the project if a
            reference genome cannot be determined. If a reference genome that is different
            from the primary reference genome of the project is determined or specified, it
            will become the alternative referenge genome of the project. The UCSC liftover
            tool will be automatically called to map input coordinates to the primary
            reference genome.''')
    parser.add_argument('--variant_only', action='store_true',
        help='''Import only variants. No sample will be created and all sample variants will
            be ignored.''')
    parser.add_argument('--variant_info', nargs='*', default=[],
        help='''Variant information fields to import, which should be one or more of the INFO
            fields specified in the header of VCF file, QUAL (the quality score listed in the
            sixth column of a VCF file), or FILTER (the seventh column). variant tools by default
            does not import any variant field. ''')
    parser.add_argument('--info', nargs='*', default=['DP'],
        help='''Genotype information fields to import, which should be one or more of the FORMAT
            fields specified in the header of VCF file, QUAL (the quality score listed in the
            sixth column of a VCF file), or FILTER (the seventh column). If read depth (DP) is
            specified but no DP field is found in genotype FORMAT, the average variant depth
            (variant info DP divided by the number of samples) will be used as individual read
            depth. variant tools by default does not import any genotype field.''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.''')


def importVCF(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = vcfImporter(proj=proj, files=args.input_files, build=args.build,
                variant_only=args.variant_only, info=[] if args.variant_only else args.info,
                force=args.force)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)


def importTxtArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''A list of files that will be imported. The file should be delimiter
            separated with format described by parameter --format. Gzipped files are
            acceptable. If parameter --update is specified, the input files will be used
            to update fields of existing variants. Otherwise, this command will import
            variants and optional genotypes of one sample from input files with fields
            chr, pos, ref and alt.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, it is assumed to be the primary reference genome of the project.
            If a reference genome that is different from the primary reference genome of the
            project is specified, it will become the alternative referenge genome of the
            project. The UCSC liftover tool will be automatically called to map input
            coordinates to the primary reference genome.''')
    parser.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types (use 'vtools show formats' to list them, or 
            'vtools show format FMT' for details about a specific format), or a local
            format specification file (with extension .fmt,
            see http://varianttools.sourceforge.net/Format/New for details). If
            unspecified, variant tools will try to guess format from file extension.
        ''')
    parser.add_argument('--sample_name', nargs='*', default=[],
        help='''Name of the samples imported by the text file. If samples are specified
            for input files without genotype, samples will be created with NULL genotype.
            If no sample name is specified for input files with genotype, a sample with
            NULL sample name will be created.''')
    parser.add_argument('--variant_info', nargs='*', default=None,
        help='''Variant information fields to be imported. This option overrides variant
            info fields of the format used and can be used to import fewer or alternative
            fields defined in a format.'''),
    parser.add_argument('--genotype_info', nargs='*', default=None,
        help='''Genotype information fields to be imported. This option overrides genotype
            info fields of the format used and can be used to import fewer or alternative
            fields defined in a format.'''),
    parser.add_argument('--update', 
        help='''Add or update fields of existing variants of specified variant table
            instead of adding new variants to the master variant table. This option
            is required for position, range or field input files because these input files
            can only be used to update existing variants.''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.''')

def importTxt(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = txtImporter(proj=proj, files=args.input_files,
                build=args.build, format=args.format, sample_name=args.sample_name,
                variant_info=args.variant_info, genotype_info=args.genotype_info,
                update=args.update, force=args.force)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

