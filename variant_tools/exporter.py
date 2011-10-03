#!/usr/bin/env python
#
# $File: exporter.py $
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


class TextExporter:
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

class Exporter:
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
        # for #record, #sample variant, #variant, new SNV, insertion, deletion, complex variants, invalid record, updated record
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
        if mode == 'insert':
            self.proj.dropIndexOnMasterVariantTable()
        #
        self.createLocalVariantIndex()

        #Importer.__init__(self, proj, files, build, force, mode='insert')
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
            self.processor = LineImporter(fmt.fields, [(1,)], fmt.delimiter, self.logger)
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
            self.proj.createNewSampleVariantTable('{0}_genotype.sample_variant_{1}'.format(self.proj.name, cur.lastrowid),
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
                self.logger.info('{:,} variants are updated'.format(self.count[8]))
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
                self.logger.info('{:,} variants are updated'.format(self.total_count[8]))
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
            
    def getSampleName(self, filename):
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
                        for bins, rec in self.prober.process(line):
                            if header_line is None:
                                return len(rec), []
                            elif len(rec) == 0:
                                return 0, []
                            else:
                                cols = [x[0] for x in self.prober.fields]
                                if type(cols[0]) is tuple:
                                    fixed = False
                                    # mutiple ones, need to figure out the moving one
                                    for i,idx in enumerate(self.prober.raw_fields[0].index.split(',')):
                                        if ':' in idx:
                                            cols = [x[i] for x in cols]
                                            fixed = True
                                            break
                                    if not fixed:
                                        cols = [x[-1] for x in cols]
                                header = [x.strip() for x in header_line.split()] # #self.prober.delimiter)]
                                if max(cols) - min(cols)  < len(header):
                                    return len(rec), [header[len(header) - self.prober.nColumns + x] for x in cols]
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
                    numSample, names = self.getSampleName(input_filename)
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
                    numSample, names = self.getSampleName(input_filename)
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
            genotype_insert_query = {id: 'INSERT INTO {0}_genotype.sample_variant_{1} VALUES ({2});'\
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

# Functions provided by this script
#
#

def exportArguments(parser):
    parser.add_argument('table', help='''A variant table whose variants will be exported.'''),
    parser.add_argument('filename', help='''Name of output file.'''),
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Samples that will be exported, specified by conditions such as 'aff=1'
            and 'filename like "MG%%"'. Multiple samples could be exported to a
            file if the output format allows. No sample will be outputted if this
            parameter is ignored.''')
    parser.add_argument('--format',
        help='''Format of the exported file. It can be one of the variant tools
            supported file types such as VCF (c.f. 'vtools show formats') or a local
            format specification file (with extension .fmt). If unspecified, variant
            tools will try to guess format from file extension. Some formats accept
            additional parameters (c.f. 'vtools show format FMT') and allows you to
            export additional or alternative fields.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the exported data. It
            can only be one of the primary (default) of alternative (if exists) reference
            genome of the project.'''),

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            exporter = txtExporter(proj=proj, table=args.table, filename=args.filename,
                samples=args.samples, format=args.format, build=args.build, 
                fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        sys.exit(e)


