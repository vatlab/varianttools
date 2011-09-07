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
from .project import Project
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
            self.logger.info('Reading coordinates in preparation for importing from alternative reference genome')
            # we need to run lift over to convert coordinates before importing data.
            alt_coordinates = set()
            for f in self.files:
                alt_coordinates |= set(self.getCoordinates(f))
            tool = LiftOverTool(self.proj)
            self.coordinateMap = tool.mapCoordinates(alt_coordinates, self.build, self.proj.build)
            self.logger.info('{:,} coordinates are mapped from {} to {}, variants in {} unmapped records will have NULL chr and pos.'.format(len(self.coordinateMap),
                self.build, self.proj.build, len(alt_coordinates) - len(self.coordinateMap)))
            self.variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt, alt_bin, alt_chr, alt_pos) VALUES ({0}, {0}, {0}, {0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
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
            cur.execute('SELECT variant_id, chr, pos, ref, alt, alt_chr, alt_pos FROM variant;')
        else:
            cur.execute('SELECT variant_id, chr, pos, ref, alt FROM variant;')
        prog = ProgressBar('Getting existing variants', numVariants)
        for count, rec in enumerate(cur):
            if self.import_alt_build:
                self.variantIndex[(rec[1], rec[2], rec[3], rec[4], rec[5], rec[6])] = rec[0]
            else:
                self.variantIndex[(rec[1], rec[2], rec[3], rec[4])] = rec[0]
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
        if self.import_alt_build:
            alt_chr, input_pos = chr, pos 
            try:
                # find chr and pos in primary reference genome
                chr, pos = self.coordinateMap[(alt_chr, input_pos)]
            except:
                # Coordinate cannot be mapped to the primary reference genome.
                bin = chr = pos = None
            # handling coordinate in alternative reference genome
            alt_bin, alt_pos, ref, alt = normalizeVariant(input_pos, ref, alt)
            # if there is a valid coordinate in the primary reference genome
            if pos:
                # if alt_pos is shifted, pos will also be shifted
                pos += alt_pos - input_pos
                bin = getMaxUcscBin(pos - 1, pos)
            var_key = (chr, pos, ref, alt, alt_chr, alt_pos)
        else:
            bin, pos, ref, alt = normalizeVariant(pos, ref, alt)
            var_key = (chr, pos, ref, alt)
        #
        try:
            return self.variantIndex[var_key]
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
            if self.import_alt_build:
                cur.execute(self.variant_insert_query, (bin, chr, pos, ref, alt, alt_bin, alt_chr, alt_pos))
            else:
                cur.execute(self.variant_insert_query, (bin, chr, pos, ref, alt))
            variant_id = cur.lastrowid
            self.variantIndex[var_key] = variant_id
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

    def getCoordinates(self, input_filename):
        '''Get (chr, pos) of all variants, for backward liftOver'''
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        count = 0
        coordinates = []
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    if line.startswith('#'):
                        continue
                    count += 1
                    tokens = [x.strip() for x in line.split('\t')]
                    chr = tokens[0][3:] if tokens[0].startswith('chr') else tokens[0]
                    pos = int(tokens[1])
                    coordinates.append((chr, pos))
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                if count % self.db.batch == 0:
                    self.db.commit()
                    prog.update(count)
            self.db.commit()
            prog.done()
        return coordinates
        
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
    def __init__(self, proj, files, col, build, delimiter, zero, force):
        Importer.__init__(self, proj, files, build, force)
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        self.col = [x - 1 for x in col]
        if len(self.col) != 4:
            raise ValueError('Four columns are required for each variant (chr, pos, ref, and alt)')
        self.delimiter = delimiter
        self.zero = zero

    def getCoordinates(self, input_filename):
        '''Read coordinates in prepration for reverse liftOver'''
        count = 0
        coordinates = []
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    if line.startswith('#'):
                        continue
                    count += 1
                    # get data
                    tokens = [x.strip() for x in line.split(self.delimiter)]
                    chr, pos = tokens[self.col[0]], tokens[self.col[1]]
                    if chr.startswith('chr'):
                        chr = chr[3:]
                    pos = int(pos) + 1 if self.zero else int(pos)
                    coordinates.append((chr, pos))
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                if count % self.db.batch == 0:
                    self.db.commit()
                    prog.update(count)
            self.db.commit()
            prog.done()
        return coordinates

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
                    self.count[0] += 1
                    # get data
                    tokens = [x.strip() for x in line.split(self.delimiter)]
                    chr, pos, ref, alt = [tokens[x] for x in self.col]
                    if chr.startswith('chr'):
                        chr = chr[3:]
                    pos = int(pos) + 1 if self.zero else int(pos)
                    if len(ref) != 1:
                        raise ValueError('Incorrect reference allele: {}'.format(ref))
                    if len(alt) != 1:
                        raise ValueError('Incorrect alternative allele: {}'.format(alt))
                    # variant
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
    grp.add_argument('-c', '--columns', default=[1,2,3,4], nargs='+', type=int,
        help='Columns for chromosome, position, reference and alternative alleles.')
    grp.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter, default to tab, a popular alternative is ',' for csv output''')
    grp.add_argument('-z', '--zero', action='store_true',
        help='''Whether or not specified file uses zero-based index. If unspecified, the
            position column is assumed to be 1-based.''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.''')

def importTxt(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = txtImporter(proj=proj, files=args.input_files,
                col=args.columns, build=args.build,
                delimiter=args.delimiter, zero=args.zero, force=args.force)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

