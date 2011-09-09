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

class SplitField:
    def __init__(self, sep=',', default=[]):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. Return default if unsuccessful. These items
        will lead to multiple records in the database.'''
        self.sep = sep
        self.default = default
    
    def __call__(self, item):
        try:
            return item.split(self.sep)
        except:
            return self.default

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

class ExtractValue:
    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that returns the value after name in one of the fields,
        and a default value if no such field is found. No space is allowed between 
        delimiter and the name.'''
        self.name = name
        self.sep = sep
        self.pos = len(name)
        self.default = default

    def __call__(self, item):
        if self.name not in item:
            return self.default
        for field in item.split(self.sep):
            if field.startswith(self.name):
                return field[self.pos:]
        return self.default

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

class Nullify:
    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return None if item == self.val else item

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
            if type(item) == list:
                if type(item[0]) == list:
                    raise ValueError('Nested vector extracted is not allowed')
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item


class TextProcessor:
    '''An intepreter that read a record, process it and return processed records.'''
    def __init__(self, fields, build, logger):
        '''Fields: a list of fields with index, adj (other items are not used)
        builds: index(es) of position, reference allele and alternative alleles. If 
            positions are available, UCSC bins are prepended to the records. If reference
            and alternative alleles are available, the records are processed for correct
            format of ref and alt alleles.
        '''
        self.logger = logger
        self.build = build
        self.fields = []
        for field in fields:
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
                idx = [int(x) - 1 for x in field.index.split()]
                if len(idx) == 1:
                    self.fields.append((idx[0], e))
                else:
                    self.fields.append((tuple(idx), e))
            except Exception as e:
                self.logger.debug(e)
                raise ValueError('Incorrect value adjustment functor or function: {}'.format(field.adj))

    def process(self, line):
        tokens = [x.strip() for x in line.split('\t')]
        records = []
        num_records = 1
        #
        for col, adj in self.fields:
            item = tokens[col] if type(col) == int else '\t'.join([tokens[x] for x in col])
            if adj is not None:
                try:
                    item = adj(item)
                    if type(item) == list:
                        if len(item) == 1:
                            # trivial case
                            item = item[0]
                        elif num_records == 1:
                            # these records will be handled separately.
                            num_records = len(item)
                        elif num_records != len(item):
                            raise ValueError('Fields in a record should generate the same number of annotations.')
                except Exception as e:
                    self.logger.debug(e)
                    # missing ....
                    item = None
            #
            records.append(item)
        # handle records
        if not self.build:
            # there is no build information, this is 'field' annotation, nothing to worry about
            if num_records == 1:
                yield records
            else:
                for i in range(num_records):
                    yield [x[i] if type(x) == list else x for x in records]
        elif len(self.build[0]) == 1:
            if num_records == 1:
                # there is no ref and alt, easy
                bins = []
                for pos_idx, in self.build:
                    try:
                        # zero-based: v, v+1 (adj=1)
                        # one-based:  v-1, v (adj=0)
                        bin = getMaxUcscBin(int(records[pos_idx]) - 1, int(records[pos_idx]))
                    except:
                        # position might be None (e.g. dbNSFP has complete hg18 coordinates,
                        # but incomplete hg19 coordinates)
                        bin = None
                    bins.append(bin)
                yield bins + records
            else:
                for i in range(num_records):
                    # get the i-th record 
                    rec = [x[i] if type(x) == list else x for x in records]
                    bins = []
                    for pos_idx, in self.build:
                        try:
                            bin = getMaxUcscBin(int(rec[pos_idx]) - 1, int(rec[pos_idx]))
                        except:
                            bin = None
                        bins.append(bin)
                    yield bins + rec
        elif num_records == 1:
            # variant, single record, although alt can lead to multiple records
            #
            # We assume that the fields for ref and alt are the same for multiple reference genomes.
            # Otherwise we cannot do this.
            all_alt = records[self.build[0][2]].split(',')
            # support multiple alternative alleles
            for input_alt in all_alt:
                bins = []
                # if there are multiple alternative alleles, we need to use a copy of records because
                # some of them will be adjusted. Otherwise, we can change the record itself. Note that
                # rec = records does not cost much in Python.
                rec = [x for x in records] if len(all_alt) > 0 else records
                for pos_idx, ref_idx, alt_idx in self.build:
                    bin, pos, ref, alt = normalizeVariant(int(rec[pos_idx]) if rec[pos_idx] else None, rec[ref_idx], input_alt)
                    # these differ build by build
                    bins.append(bin)
                    rec[pos_idx] = pos
                    rec[ref_idx] = ref
                    rec[alt_idx] = alt
                yield bins + rec
        else:
            all_alt = records[self.build[0][2]].split(',')
            for input_alt in all_alt:
                for i in range(num_records):
                    bins = []
                    rec = [x[i] if type(x) == list else x for x in records]
                    for pos_idx, ref_idx, alt_idx in self.build:
                        bin, pos, ref, alt = normalizeVariant(int(rec[pos_idx]) if rec[pos_idx] else None, rec[ref_idx], input_alt)
                        bins.append(bin)
                        rec[pos_idx] = pos
                        rec[ref_idx] = ref
                        rec[alt_idx] = alt
                    yield bins + rec



class Importer:
    '''A general class for importing variants'''
    def __init__(self, proj, files, build, force):
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
        for f in files:
            filename = os.path.split(f)[-1]
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
                    self.files.append(f)
                else:
                    self.logger.info('Ignoring imported file {}'.format(filename))
            else:
                self.files.append(f)
        # for all record, new SNV, insertion, deletion, complex variants, and invalid record
        self.count = [0, 0, 0, 0, 0, 0]
        self.total_count = [0, 0, 0, 0, 0, 0]
        if len(self.files) == 0:
            return
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
        self.import_alt_build = False
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
            self.logger.info('Adding an alternative referenge genome ({}) to the project.'.format(self.build))
            tool = LiftOverTool(self.proj)
            tool.setAltRefGenome(self.build)
            self.import_alt_build = True
        else:
            raise ValueError('Specified build {} does not match either the primary '.format(self.build) + \
                ' {} or the alternative reference genome of the project.'.format(self.proj.build, self.proj.alt_build))
        #
        # importing data to alternative reference genome
        if self.import_alt_build:
            self.variant_insert_query = 'INSERT INTO variant (alt_bin, alt_chr, alt_pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        else:
            self.variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        #
        self.proj.dropIndexOnMasterVariantTable()
        #
        self.createLocalVariantIndex()

    def __del__(self):
        self.proj.createIndexOnMasterVariantTable()

    def guessBuild(self, file):
        # by default, reference genome cannot be determined from file
        return None

    def openFile(self, filename):
        if filename.lower().endswith('.gz'):
            return gzip.open(filename, 'rb')
        else:
            # text file
            return open(filename, 'r')

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

    def recordFileAndSample(self, filename, sampleNames, sampleFields = []):
        cur = self.db.cursor()
        cur.execute("INSERT INTO filename (filename) VALUES ({});".format(self.db.PH), (filename,))
        filenameID = cur.lastrowid
        self.proj.createSampleTableIfNeeded()
        sample_ids = []
        s = delayedAction(self.logger.info, 'Creating {} sample variant tables'.format(len(sampleNames)))
        for samplename in sampleNames:
            cur.execute('INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'.format(self.db.PH),
                (filenameID, samplename))
            sample_ids.append(cur.lastrowid)
            self.proj.createNewSampleVariantTable('{0}_genotype.sample_variant_{1}'.format(self.proj.name, cur.lastrowid),
                sampleFields)
        del s
        return sample_ids
        
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

    def importData(self):
        '''Start importing'''
        for count,f in enumerate(self.files):
            self.logger.info('Importing genotype from {} ({}/{})'.format(f, count + 1, len(self.files)))
            self.importFromFile(f)
            self.logger.info('{:,} new variants from {:,} records are imported, with {:,} SNVs, {:,} insertions, {:,} deletions, and {:,} complex variants.{}'\
                .format(sum(self.count[1:-1]), self.count[0], self.count[1], self.count[2], self.count[3], self.count[4],
                ' {} invalid records are ignored'.format(self.count[5]) if self.count[5] > 0 else ''))
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
        if len(self.files) > 1:
            self.logger.info('{:,} new variants from {:,} records in {} files are imported, with {:,} SNVs, {:,} insertions, {:,} deletions, and {:,} complex variants.{}'\
                .format(sum(self.total_count[1:-1]), self.total_count[0], len(self.files), self.total_count[1], self.total_count[2], self.total_count[3], self.total_count[4],
                ' {} invalid records are ignored'.format(self.total_count[5]) if self.total_count[5] > 0 else ''))
        if self.total_count[0] > 0 and self.proj.alt_build is not None:
            coordinates = set([(x[0], x[1]) for x,y in self.variantIndex.iteritems() if y[1] == 1])
            # we need to run lift over to convert coordinates before importing data.
            tool = LiftOverTool(self.proj)
            if self.import_alt_build:
                self.logger.info('Mapping {} new variants from {} to {} reference genome'.format(self.total_count[0], self.proj.alt_build, self.proj.build))
                coordinateMap = tool.mapCoordinates(coordinates, self.proj.alt_build, self.proj.build)
                query = 'UPDATE variant SET bin={0}, chr={0}, pos={0} WHERE variant_id={0};'.format(self.db.PH)
            else:
                self.logger.info('Mapping {} new variants from {} to {} reference genome'.format(self.total_count[0], self.proj.build, self.proj.alt_build))
                coordinateMap = tool.mapCoordinates(coordinates, self.proj.build, self.proj.alt_build)
                query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(self.db.PH)
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

    def guessBuild(self, filename):
        '''Called by the initializer to determine reference genome
        '''
        with self.openFile(filename) as input:
            for line in input:
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
            line = input.readline()
            if not line.startswith('##fileformat=VCF'):
                self.logger.error('Invalid vcf file. file not started with line ##fileformat')
                raise ValueError('Invalid vcf file')
            if not line.strip().endswith('4.0') and not line.strip().endswith('4.1'):
                raise ValueError('This importer tool only supports VCF format v4.0 and v4.1. \
                    Please use vcftools to convert your vcf file to a supported format')
            #
            for line in input:
                if line.startswith('#CHR'):
                    samples = line.split()[9:]
                if not line.startswith('#'):
                    break
        return samples

    def importFromFile(self, input_filename):
        '''Import a VCF file to sample_variant'''
        #
        # handle meta information and get sample names
        sampleNames = self.getMetaInfo(input_filename)
        #
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        no_sample = self.variant_only or len(sampleNames) == 0
        sample_ids = self.recordFileAndSample(os.path.split(input_filename)[-1], [None] if no_sample else sampleNames, 
            ['DP'] if self.import_depth else [])   # record individual depth, total depth is divided by number of sample in a file
        #
        nSample = len(sample_ids)
        #
        DP_pattern = re.compile('.*DP=(\d+)')
        #
        cur = self.db.cursor()
        # sample variants are inserted into different tables in a separate database.
        sample_variant_insert_query = {x: 'INSERT INTO {1}_genotype.sample_variant_{3} VALUES ({0}, {0} {2});'\
            .format(self.db.PH, self.proj.name, ',' + self.db.PH if self.import_depth else '', x) for x in sample_ids}
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    if line.startswith('#'):
                        continue
                    self.count[0] += 1
                    # get data
                    tokens = [x.strip() for x in line.split('\t')]
                    chr = tokens[0][3:] if tokens[0].startswith('chr') else tokens[0]
                    pos = int(tokens[1])
                    ref = tokens[3]
                    # we only extract info get depth.
                    if self.import_depth:
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
                        if no_sample:
                            cur.execute(sample_variant_insert_query[sample_ids[0]], [variant_id, 1] + DP)
                        else:
                            variants = [x.split(':')[GT_idx].count('1') for x in tokens[-len(sample_ids):]]
                            for var_idx, var in enumerate(variants):
                                if var != 0:  # genotype 0|0 are ignored
                                    cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id, var] + DP)
                    else:
                        # now, this is the common case with insertion, deletion, and multiple alternative variants
                        alts = tokens[4].split(',')
                        variant_id = [0] * len(alts)
                        for altidx, alt in enumerate(alts):
                            variant_id[altidx] = self.addVariant(cur, chr, pos, ref, alt)
                        if no_sample:
                            for i in range(len(alts)):
                                cur.execute(sample_variant_insert_query[sample_ids[0]], [variant_id[i], 1] + DP)
                        else:
                            # process variants
                            for var_idx, var in enumerate([x.split(':')[GT_idx] for x in tokens[-len(sample_ids):]]):
                                if len(var) == 3:  # regular
                                    gt = var[0] + var[2]  # GT can be separated by / or |
                                    if gt in ['01', '10']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif gt in ['02', '20']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
                                    elif gt == '11':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 2] + DP)
                                    elif gt in ['12', '21']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], -1] + DP)
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], -1] + DP)
                                    elif gt == '22':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 2] + DP)
                                    elif gt == '00':
                                        pass
                                    else:
                                        raise ValueError('I do not know how to process genotype {}'.format(var))
                                else: # should have length 1
                                    if var == '1':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif var == '2':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
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
    def __init__(self, proj, files, build, format, force):
        Importer.__init__(self, proj, files, build, force)
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        try:
            fmt = fileFMT(format)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Input file format {} is not currently supported by variant tools'.format(item))
        #
        pos_idx = [i for i,x in enumerate(fmt.fields) if x.name == fmt.variant_fields[1]][0]
        ref_idx = [i for i,x in enumerate(fmt.fields) if x.name == fmt.variant_fields[2]][0]
        alt_idx = [i for i,x in enumerate(fmt.fields) if x.name == fmt.variant_fields[3]][0]
        self.processor = TextProcessor(fmt.fields, [(pos_idx, ref_idx, alt_idx)], self.logger)

    def importFromFile(self, input_filename):
        '''Import a TSV file to sample_variant'''
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        filename = os.path.split(input_filename)[-1]
        # assuming one sample for each file
        sample_id = self.recordFileAndSample(filename, [None])[0]
        #
        cur = self.db.cursor()
        sample_variant_insert_query = 'INSERT INTO {0}_genotype.sample_variant_{1} VALUES ({2}, {2});'\
            .format(self.proj.name, sample_id, self.db.PH)
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    if line.startswith('#'):
                        continue
                    for rec in processor.process(line):
                        self.count[0] += 1
                        cur.execute(insert_query, rec)
                    #
                    variant_id = self.addVariant(cur, chr, pos, ref, alt)
                    # sample variant, the variant type is always hetero???
                    cur.execute(sample_variant_insert_query, (variant_id, 1))
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
        help='''Import only variants. Sample variants will be ignored but a sample with
            NULL sample name will still be created to trace the source of variants.''')
    parser.add_argument('--info', nargs='*', default=['DP'],
        help='''Variant information fields to import. This command only support
            'DP' (total depth). When 'DP' is listed (default), vtools will look
            for total depth (DP=) in the INFO field of each variant and set average
            depth to each individual (DP/numSample in the vcf file) in field 'DP'.''')
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
        help='''A list of files that will be imported. The file should be in 
            tab or command separated value format. Gzipped files are acceptable.''')
    grp = parser.add_argument_group('Description of input files')
    grp.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, it is assumed to be the primary reference genome of the project.
            If a reference genome that is different from the primary reference genome of the
            project is specified, it will become the alternative referenge genome of the
            project. The UCSC liftover tool will be automatically called to map input
            coordinates to the primary reference genome.''')
    grp.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types (use 'vtools show formats' to list them, or 
            'vtools show format FMT' for details about a specific format), or a local
            format specification file (with extension .fmt,
            see http://varianttools.sourceforge.net/Format/New for details).''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.''')

def importTxt(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = txtImporter(proj=proj, files=args.input_files,
                build=args.build, format=args.format, force=args.force)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

