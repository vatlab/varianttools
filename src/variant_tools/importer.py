#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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
import re
import sys
from itertools import repeat

from .geno_store import GenoStore
from .liftOver import LiftOverTool
from .preprocessor import SequentialExtractor, IgnoredRecord
from .preprocessor import *
from .project import Project, fileFMT
from .text_reader import TextReader
from .utils import (ProgressBar, RefGenome, delayedAction, env, getMaxUcscBin,
                    lineCount, openFile)

try:
    from variant_tools.cgatools import normalize_variant
except ImportError as e:
    sys.exit(
        'Failed to import module ({})\n'
        'Please verify if you have installed variant tools successfully (using command '
        '"python setup.py install")'.format(e))



#
#
# Process each line using the above functors
#
#
class LineProcessor:
    '''An interpreter that read a record (a line), process it and return processed records.'''

    def __init__(self, fields, build_info, delimiter, ranges):
        '''
        fields: a list of fields with index, adj (other items are not used)
        build_info: reference genome and index(es) of chromosome, position, reference allele and alternative alleles. If
            positions are available, UCSC bins are prepended to the records. If reference
            and alternative alleles are available, the records are processed for correct
            format of ref and alt alleles.
        delimiter: how to split line
        ranges: range of fields (var, var_info, GT, GT_info), used to determine
            var_info and GT fields when subsets of samples are imported.
        '''
        self.build_info = build_info
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
        self.nColumns = 0  # number of columns
        self.import_var_info = True
        self.import_sample_range = None  # genotype fields might be disabled
        self.maxInputCol = 0  # sometimes processor do not have to split input all the way through

    def reset(self, import_var_info=True, import_sample_range=None):
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
                tokens = [
                    x.strip()
                    for x in tokens.split(self.delimiter, self.maxInputCol)
                ][:self.maxInputCol]
        else:
            if type(tokens) is not list:
                tokens = [x.strip() for x in tokens.split(self.delimiter)]
            self.nColumns = len(tokens)
            cIdx = 0  # column index
            num_sample = -1  # number of samples ...
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
                    if not field.index:
                        raise ValueError(
                            'Invalid index ({}) for field {}. Perhaps it is not defined in the .fmt file.'
                            .format(field.index, field.name))
                    for x in field.index.split(','):
                        if ':' in x:
                            # a slice
                            if x.count(':') == 1:
                                start, end = list(map(str.strip, x.split(':')))
                                step = None
                            else:
                                start, end, step = list(
                                    map(str.strip, x.split(':')))
                                step = int(step) if step else None
                            start = int(start) - 1 if start else None
                            if end.strip():
                                if int(end) >= 0:  # position index, shift by 1
                                    end = int(end) - 1
                                else:  # relative to the back, do not move
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
                            self.columnRange[fIdx] = (cIdx, cIdx + 1)
                            self.maxInputCol = max(self.maxInputCol, indexes[0])
                            cIdx += 1
                        # case of index=7,8,9
                        else:
                            # a tuple
                            self.fields.append((tuple(indexes), False, e))
                            self.columnRange[fIdx] = (cIdx, cIdx + 1)
                            self.maxInputCol = max(self.maxInputCol,
                                                   max(indexes))
                            cIdx += 1
                    # if there is only one slice
                    # case of index=8::2
                    elif len(indexes) == 1:
                        # single slice
                        cols = list(range(len(tokens)))[indexes[0]]
                        if self.import_sample_range is not None:
                            # limiting the columns to import
                            if self.import_sample_range[0] >= len(
                                    cols) or self.import_sample_range[1] > len(
                                        cols):
                                raise ValueError(
                                    'ERROR PROCESSING subset of samples.')
                            cols = cols[self.import_sample_range[0]:self
                                        .import_sample_range[1]]
                        for c in cols:
                            self.fields.append((c, True, e))
                            self.maxInputCol = max(self.maxInputCol, c)
                        if num_sample == -1:
                            num_sample = len(cols)
                        elif num_sample != len(cols):
                            sys.exit(
                                'The first line of input has inconsistent number of fields for samples, perhaps due to incorrect use of delimiters.'
                            )
                        self.columnRange[fIdx] = (cIdx, cIdx + len(cols))
                        cIdx += len(cols)
                    else:
                        # we need to worry about mixing integer and slice
                        expanded_indexes = [
                            repeat(s, len(tokens)) if type(s) == int else list(
                                range(len(tokens)))[s] for s in indexes
                        ]
                        count = 0
                        for idx, c in enumerate(zip(*expanded_indexes)):
                            if self.import_sample_range is not None:
                                if idx < self.import_sample_range[
                                        0] or idx >= self.import_sample_range[1]:
                                    continue
                            count += 1
                            self.fields.append((tuple(c), False, e))
                            self.maxInputCol = max(self.maxInputCol, max(c))
                        if num_sample == -1:
                            num_sample = count
                        elif num_sample != count:
                            sys.exit(
                                'The first line of input has inconsistent number of fields for samples, perhaps due to incorrect use of delimiters.'
                            )
                        self.columnRange[fIdx] = (cIdx, cIdx + count)
                        cIdx += count
                except Exception as e:
                    env.logger.error(
                        'Incorrect value adjustment functor or function {}: {}'
                        .format(field.adj, e))
                    sys.exit(1)
            self.first_time = False
            self.maxInputCol += 1  # use 1-indexes maxInputCol
        #
        try:
            # we first trust that nothing can go wrong and use a quicker method
            records = [(tokens[col] if t else [tokens[x] for x in col]) if adj is None else \
                (adj(tokens[col]) if t else adj([tokens[x] for x in col])) for col,t,adj in self.fields]
        except IgnoredRecord:
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
                        env.logger.debug(
                            'Failed to process field {}: {}'.format(item, e))
                        # missing ....
                        item = None
                records.append(item)
        #
        num_records = max([
            len(item) if type(item) is tuple else 1 for item in records
        ]) if records else 1
        # handle records
        if not self.build_info:
            # there is no build information, this is 'field' annotation, nothing to worry about
            if num_records == 1:
                yield [], [x[0] if type(x) is tuple else x for x in records]
            else:
                for i in range(num_records):
                    yield [], [(x[i] if i < len(x) else None)
                               if type(x) is tuple else x for x in records]
        elif len(self.build_info[0]) == 1:
            # this is a range format
            for i in range(num_records):
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) is tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None)
                           if type(x) is tuple else x for x in records]
                bins = [
                    getMaxUcscBin(int(rec[pos_idx]) - 1, int(rec[pos_idx]))
                    if rec[pos_idx] else None for pos_idx, in self.build_info
                ]
                yield bins, rec
        else:
            # this is a variant format
            for i in range(num_records):
                bins = []
                if i == 0:  # try to optimize a little bit because most of the time we only have one record
                    rec = [x[0] if type(x) is tuple else x for x in records]
                else:
                    rec = [(x[i] if i < len(x) else None)
                           if type(x) is tuple else x for x in records]
                for ref_genome, chr_idx, pos_idx, ref_idx, alt_idx in self.build_info:
                    # bin, pos, ref, alt = normalizeVariant(int(rec[pos_idx]) if rec[pos_idx] else None, rec[ref_idx], rec[alt_idx])
                    #env.logger.error('PRE {} {} {} {}'.format(rec[chr_idx], rec[pos_idx], rec[ref_idx], rec[alt_idx]))
                    msg = normalize_variant(ref_genome.crr, rec, chr_idx,
                                            pos_idx, ref_idx, alt_idx)
                    #env.logger.error('POST {} {} {} {}'.format(rec[chr_idx], rec[pos_idx], rec[ref_idx], rec[alt_idx]))
                    #
                    if msg:
                        # if the message says 'Unrecognized allele', the variant will be ignored.
                        if msg[0] == 'U':
                            raise ValueError(msg)
                        else:
                            env.logger.warning('{}: {}'.format(
                                ref_genome.name, msg))
                    # normalization will convert rec[pos_idx] to int if necessary
                    bins.append(getMaxUcscBin(rec[pos_idx] - 1, rec[pos_idx]))
                yield bins, rec


#
#
# utility function to get sample name
#
def probeSampleName(filename, prober, encoding):
    '''Probe text file for sample name. Essentially speaking
    this function will go to the last comment line, break it in pieces
    and see if we can grab some headers'''
    header_line = None
    count = 0
    with openFile(filename) as input:
        try:
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
                                    for i, idx in enumerate(
                                            prober.raw_fields[0].index.split(
                                                ',')):
                                        if ':' in idx:
                                            cols = [x[i] for x in cols]
                                            fixed = True
                                            break
                                    if not fixed:
                                        cols = [x[-1] for x in cols]
                                header = [
                                    x.strip() for x in header_line.split()
                                ]  # #prober.delimiter)]
                                if max(cols) - min(cols) < len(header) and len(
                                        header) > max(cols):
                                    return len(rec), [
                                        header[len(header) - prober.nColumns +
                                               x] for x in cols
                                    ]
                                else:
                                    header = [
                                        x.strip() for x in header_line.split(
                                            prober.delimiter)
                                    ]
                                    if max(cols) - min(cols) < len(
                                            header) and len(header) > max(cols):
                                        return len(rec), [
                                            header[len(header) -
                                                   prober.nColumns + x]
                                            for x in cols
                                        ]
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
                            raise ValueError(
                                'No genotype column could be determined after 1000 lines.'
                            )
                        env.logger.debug(e)
        except TypeError as e:
            # Python 2.7.4 and 3.3.1 have a regression bug that prevents us from opening
            # certain types of gzip file (http://bugs.python.org/issue17666).
            if filename.lower().endswith('.gz'):
                raise RuntimeError(
                    'Failed to open gzipped file {} due to a bug '
                    'in Python 2.7.4 and 3.3.1. Please use a different version '
                    'of Python or decompress this file manually.'.format(
                        filename))
            else:
                raise e
        # if there is no data line, we cannot use cols to determine sample name
        return 0, []


#
#
#  Command import
#
#


class Importer:
    '''A general class for importing variants'''

    def __init__(self,
                 proj,
                 files,
                 build,
                 format,
                 sample_name=None,
                 force=False,
                 jobs=1,
                 fmt_args=[],
                 sort=False):
        self.proj = proj
        self.db = proj.db
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
                env.logger.info(
                    'Using primary reference genome {} of the project.'.format(
                        self.build))
            else:
                raise ValueError(
                    'Please specify a reference genome using parameter --build')
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
            env.logger.warning(
                'The new files uses a different reference genome ({}) from the primary reference genome ({}) of the project.'
                .format(self.build, self.proj.build))
            env.logger.info(
                'Adding an alternative reference genome ({}) to the project.'
                .format(self.build))
            tool = LiftOverTool(self.proj)
            # in case of insert, the indexes will be dropped soon so do not build
            # index now
            tool.setAltRefGenome(self.build, build_index=False)
            self.import_alt_build = True
        #
        self.jobs = max(1, jobs)
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError(
                'Please specify the reference genome of the input data.')
        #
        # try to guess file type
        self.format = format
        if not format:
            filename = files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
                self.format = format
            else:
                raise ValueError(
                    'Cannot guess input file format from filename "{}"'.format(
                        files[0]))

        try:
            fmt = fileFMT(format, fmt_args)
        except Exception as e:
            env.logger.debug(e)
            raise IndexError(
                'Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '
                .format(e, format))
        # if there are any invalid field, quite with an error message
        # (field idx+1 can be used for output, but not for input)
        for fld in fmt.fields:
            if fld.index is None:
                raise ValueError(
                    'Cannot import field {} from input file.'.format(fld.name))
        #
        if fmt.preprocessor is not None:
            env.logger.info(
                'Preprocessing data [{}] to generate intermediate input files for import'
                .format(', '.join(files)))
            # if this is the case, only one input stream will be allowed.
            # process command line
            command = fmt.preprocessor
            # replace command with other stuff, if applicable
            command = command.replace('$build', "'{}'".format(self.build))
            #
            # create temp files
            temp_files = [
                os.path.join(env.cache_dir,
                             os.path.basename(x) + '.' + fmt.name)
                for x in files
            ]
            try:
                processor = eval(command)
                # intermediate files will be named as "cache_dir/$inputfilename.$(fmt.name)"
                processor.convert(files, temp_files)
                for output in temp_files:
                    if not os.path.isfile(output):
                        raise ValueError(
                            "Preprocessed file {} does not exist.".format(
                                output))
            except Exception as e:
                raise ValueError("Failed to execute preprocessor '{}': {}".\
                                 format(re.sub(r'\((.*)\)', '', command), e))
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
                    env.logger.warning(
                        'Re-importing imported file {}, duplicated samples may occur.'
                        .format(filename))
                    self.files.append(filename)
                else:
                    env.logger.info(
                        'Ignoring imported file {}'.format(filename))
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
        self.variant_fields = [
            x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]
        ]
        self.variant_info = [
            x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]
        ]
        #
        # NOTE: for a given format, number of genotype_fields and associated genotype info
        # are determined by the format file, and are fixed, but the number of samples in
        # a file can vary. In particular, there can be no sample, and in this case, if a
        # sample name is given, a sample will be created with only variant ids.
        #
        self.genotype_field = [
            x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]
        ]
        self.genotype_info = [
            x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]
        ]

        if 'GT' in self.genotype_info:
            raise ValueError(
                'GT (genotype) field should not be explicitly specified.')
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for chr, pos, ref, alt are 0, 1, 2, 3 in fields.
            self.processor = LineProcessor(
                fmt.fields, [(RefGenome(self.build), 0, 1, 2, 3)],
                fmt.delimiter, self.ranges)
        else:  # position or range type
            raise ValueError(
                'Can only import data with full variant information (chr, pos, ref, alt)'
            )
        # probe number of samples
        if self.genotype_field:
            self.prober = LineProcessor([fmt.fields[fmt.ranges[2]]], [],
                                        fmt.delimiter, self.ranges)
        # there are variant_info
        if self.variant_info:
            cur = self.db.cursor()
            headers = self.db.getHeaders('variant')
            for f in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]:
                # either insert or update, the fields must be in the master variant table
                self.proj.checkFieldName(f.name, exclude='variant')
                if f.name not in headers:
                    with delayedAction(env.logger.info,
                                       'Adding column {}'.format(f.name)):
                        cur.execute('ALTER TABLE variant ADD {} {};'.format(
                            f.name, f.type))
        #
        if fmt.input_type != 'variant':
            env.logger.info(
                'Only variant input types that specifies fields for chr, pos, ref, alt could be imported.'
            )
        #
        self.input_type = fmt.input_type
        self.encoding = fmt.encoding
        self.header = fmt.header
        fbin, fchr, fpos = ('alt_bin', 'alt_chr',
                            'alt_pos') if self.import_alt_build else ('bin',
                                                                      'chr',
                                                                      'pos')
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH)
        self.variant_insert_query = 'INSERT INTO variant ({0}, {1}, {2}, ref, alt {3}) VALUES ({4});'\
            .format(fbin, fchr, fpos, ' '.join([', ' + x for x in self.variant_info]), ', '.join([self.db.PH]*(len(self.variant_info) + 5)))
        #
        self.variantIndex = self.proj.createVariantMap('variant',
                                                       self.import_alt_build)
        # drop index here after all possible exceptions have been raised.
        self.proj.dropIndexOnMasterVariantTable()
        self.sort = sort

    def __del__(self):
        # remove existing indexes, which will be created when the project is open
        # by a non-import command
        self.proj.dropIndexOnMasterVariantTable()

    def recordFileAndSample(self, filename, sampleNames):
        cur = self.db.cursor()
        cur.execute(
            'SELECT file_id FROM filename WHERE filename={0}'.format(
                self.db.PH), (filename,))
        rec = cur.fetchall()
        if rec:
            filenameID = rec[0][0]
        else:
            cur.execute(
                "INSERT INTO filename (filename) VALUES ({0});".format(
                    self.db.PH), (filename,))
            filenameID = cur.lastrowid
        sample_ids = []
        with delayedAction(
                env.logger.info,
                'Creating {} genotype tables'.format(len(sampleNames))):
            #
            cur.execute(
                'SELECT sample_id FROM sample where sample_name={0}'.format(
                    self.db.PH), (sampleNames[0],))
            rec = cur.fetchall()
            cur.execute('SELECT count(*) FROM sample'.format(self.db.PH))
            count = cur.fetchall()
            if len(rec) == 0 and int(count[0][0]) > 0:
                cur.execute(
                    'UPDATE project SET value={0} WHERE name={0};'.format(
                        self.db.PH), (1, "multiVCF"))

            for samplename in sampleNames:
                cur.execute(
                    'INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'
                    .format(self.db.PH),
                    (filenameID, '' if samplename is None else samplename))
                sid = cur.lastrowid
                sample_ids.append(sid)
        return sample_ids

    def addVariant(self, cur, rec):
        #
        try:
            # try to see if rec[1] can be translated to a proper integer
            int(rec[2])
        except:
            env.logger.debug(
                'Variant without valid position is ignored: {}'.format(rec))
            return None
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
            #env.logger.info('Overlapping variant {}:{}-{}/{}'.format(rec[1], rec[2], rec[3], rec[4]))
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
                env.logger.warning(
                    'Sample information is not recorded for a file without genotype and sample name.'
                )
                self.sample_in_file = []
                return ([], 0, [])
            else:
                try:
                    numSample, names = probeSampleName(input_filename,
                                                       self.prober,
                                                       self.encoding)
                    if not names:
                        if numSample == 1:
                            env.logger.debug(
                                'Missing sample name (name None is used)'
                                .format(numSample))
                            self.sample_in_file = [None]
                            return (self.recordFileAndSample(
                                input_filename, ['']), 2, names)
                        elif numSample == 0:
                            env.logger.debug(
                                'No genotype column exists in the input file so no sample will be recorded.'
                            )
                            self.sample_in_file = []
                            return ([], 0, [])
                        else:
                            raise ValueError(
                                'Failed to guess sample name. Please specify sample names for {} samples using parameter --sample-name, or add a proper header to your input file that matches the columns of samples. See "vtools import -h" for details.'
                                .format(numSample))
                    else:
                        self.sample_in_file = [x for x in names]
                        return (self.recordFileAndSample(input_filename,
                                                         names), 2, names)
                except ValueError as e:
                    # cannot find any genotype column, perhaps no genotype is defined in the file (which is allowed)
                    env.logger.warning(
                        'Cannot import genotype from the input file: {0}'
                        .format(e))
                    self.sample_in_file = []
                    return ([], 0, [])
        else:
            self.sample_in_file = [x for x in self.sample_name]

            if not self.genotype_field:
                # if no genotype, but a sample name is given
                env.logger.debug(
                    'Input file does not contain any genotype. Only the variant ownership information is recorded.'
                )
                return (self.recordFileAndSample(input_filename,
                                                 self.sample_name), 1, [])
            else:
                try:
                    numSample, names = probeSampleName(input_filename,
                                                       self.prober,
                                                       self.encoding)

                except ValueError as e:
                    env.logger.debug(e)
                    numSample = 0
                #
                if numSample == 0:
                    env.logger.warning(
                        'No genotype column could be found from the input file. Assuming no genotype.'
                    )
                    # remove genotype field from processor
                    self.processor.reset(
                        import_var_info=True, import_sample_range=[0, 0])
                    if len(self.sample_name) > 1:
                        raise ValueError(
                            "When there is no sample genotype, only one sample name is allowed."
                        )
                elif len(self.sample_name) != numSample:
                    raise ValueError(
                        '{} sample detected but only {} sample names are specified'
                        .format(numSample, len(self.sample_name)))
                return (self.recordFileAndSample(input_filename,
                                                 self.sample_name),
                        1 if numSample == 0 else 2, self.sample_in_file)

    def importVariant(self, input_filename):
        # reset text processor to allow the input of files with different number of columns
        self.processor.reset(import_var_info=True, import_sample_range=[0, 0])
        #
        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc // 200, 100), 100000)
        # one process is for the main program, the
        # other threads will handle input
        # getNew=True so the reader only read variants not in variantIndex if no additional
        # variant info is imported
        # if self.variant_info:
        #     reader = TextReader(self.processor, input_filename, None, True,
        #         env.import_num_of_readers, self.encoding, self.header, quiet=False)
        # else:
        #     reader = TextReader(self.processor, input_filename, self.variantIndex, True,
        #         env.import_num_of_readers, self.encoding, self.header, quiet=False)
        if self.variant_info:
            reader = TextReader(
                self.processor,
                input_filename,
                None,
                True,
                1,
                self.encoding,
                self.header,
                quiet=False)
        else:
            reader = TextReader(
                self.processor,
                input_filename,
                self.variantIndex,
                True,
                1,
                self.encoding,
                self.header,
                quiet=False)
        # preprocess data
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0
        for self.count[0], bins, rec in reader.records():
            variant_id = self.addVariant(cur, bins + rec[0:self.ranges[2]])
            if variant_id is None:
                continue
            if (last_count == 0 and self.count[0] > 200) or (
                    self.count[0] - last_count > update_after):
                self.db.commit()
                last_count = self.count[0]
                prog.update(self.count[0])
        prog.done()
        self.count[7] = reader.unprocessable_lines
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
        to_be_mapped = os.path.join(env.temp_dir, 'var_in.bed')
        loci_count = 0
        with open(to_be_mapped, 'w') as output:
            for key in self.variantIndex:
                for pos, status in self.variantIndex[key].items():
                    if status[1] == 1:
                        output.write('{0}\t{1}\t{2}\t{3}/{4}/{5}\n'.format(
                            key[0] if len(key[0]) > 2 else 'chr' + key[0],
                            pos - 1, pos, key[1], key[2], status[0]))
                        loci_count += 1
        # free some RAM
        self.variantIndex.clear()
        #
        if self.import_alt_build:
            env.logger.info(
                'Mapping new variants at {} loci from {} to {} reference genome'
                .format(loci_count, self.proj.alt_build, self.proj.build))
            query = 'UPDATE variant SET bin={0}, chr={0}, pos={0} WHERE variant_id={0};'.format(
                self.db.PH)
            mapped_file, err_count = tool.mapCoordinates(
                to_be_mapped, self.proj.alt_build, self.proj.build)
        else:
            env.logger.info(
                'Mapping new variants at {} loci from {} to {} reference genome'
                .format(loci_count, self.proj.build, self.proj.alt_build))
            query = 'UPDATE variant SET alt_bin={0}, alt_chr={0}, alt_pos={0} WHERE variant_id={0};'.format(
                self.db.PH)
            # this should not really happen, but people (like me) might manually mess up with the database
            with delayedAction(
                    env.logger.info,
                    'Adding alternative reference genome {} to the project.'
                    .format(self.proj.alt_build)):
                headers = self.db.getHeaders('variant')
                for fldName, fldType in [('alt_bin', 'INT'),
                                         ('alt_chr', 'VARCHAR(20)'),
                                         ('alt_pos', 'INT')]:
                    if fldName in headers:
                        continue
                    self.db.execute(
                        'ALTER TABLE variant ADD {} {} NULL;'.format(
                            fldName, fldType))
            mapped_file, err_count = tool.mapCoordinates(
                to_be_mapped, self.proj.build, self.proj.alt_build)
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
                cur.execute(query,
                            (getMaxUcscBin(pos - 1, pos), chr, pos, var_id))
                count += 1
                if count % 10000 == 0:
                    self.db.commit()
                    prog.update(count)
        self.db.commit()
        prog.done()
        env.logger.info('Coordinates of {} ({} total, {} failed to map) new variants are updated.'\
            .format(count, total_new, err_count))


def importVariantsArguments(parser):
    parser.add_argument(
        'input_files',
        nargs='+',
        help='''A list of files that will be imported. The file should be delimiter
            separated with format described by parameter --format. Gzipped files are
            acceptable. If a preprocessor is defined in the format, input files will
            be processed by the preprocessor before they are imported.''')
    parser.add_argument(
        '--build',
        help='''Build version of the reference genome (e.g. hg18) of the input data. If
            unspecified, it is assumed to be the primary reference genome of the project.
            If a reference genome that is different from the primary reference genome of the
            project is specified, it will become the alternative reference genome of the
            project. The UCSC liftover tool will be automatically called to map input
            coordinates between the primary and alternative reference genomes. If you
            are uncertain about the reference genome used in the data, use a recent
            standard reference genome (e.g. hg19) and validate it after the data is
            imported (c.f. "vtools admin --validate_build").''')
    parser.add_argument(
        '--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as VCF (cf. 'vtools show formats'), or a
            local format specification file (with extension .fmt). If unspecified,
            variant tools will try to guess format from file extension. Some file
            formats accept parameters (cf. 'vtools show format FMT') and allow you
            to import additional or alternative fields defined for the format. '''
    )
    parser.add_argument(
        '--sample_name',
        '--sample-name',
        nargs='*',
        default=[],
        help='''Name of the samples imported by the input files. The same names will be
            used for all files if multiple files are imported. If unspecified, headers
            of the genotype columns of the last comment line (line starts with #) of the
            input files will be used (and thus allow different sample names for input files).
            If sample names are specified for input files without genotype, samples
            will be created without genotype. If sample names cannot be determined from
            input file and their is no ambiguity (only one sample is imported), a sample
            with empty sample name will be created.''')
    parser.add_argument(
        '-f',
        '--force',
        action='store_true',
        help='''Import files even if the files have been imported before. This option
            can be used to import from updated file or continue disrupted import, but will
            not remove wrongfully imported variants from the master variant table.'''
    ),
    parser.add_argument(
        '-j',
        '--jobs',
        metavar='N',
        default=4,
        type=int,
        help='''Number of processes to import input file. Variant tools by default
            uses four processes to import variants and samples genotypes in
            parallel, and you can use more or less processes by adjusting this
            parameter. Due to the overhead of inter-process communication, more
            jobs do not automatically lead to better performance.'''),
    parser.add_argument(
        '--sort', action='store_true', help='''Import another VCF file.''')


def importVariants(args):

    try:
        # the project is opened with verify=False so index on the master
        # variant table will not be created if it does not exist (because the
        # last command was a vtools import command)
        with Project(
                verbosity=args.verbosity, mode='SKIP_VERIFICATION') as proj:
            importer = Importer(
                proj=proj,
                files=args.input_files,
                build=args.build,
                format=args.format,
                sample_name=args.sample_name,
                force=args.force,
                jobs=args.jobs,
                fmt_args=args.unknown_args,
                sort=args.sort)
            store = GenoStore(proj, importer)
            store.importGenotypes(importer)
            importer.finalize()
            #transform genotype stored in sqlite to hdf5
            if proj.store == "hdf5" and importer.format != "vcf":
                importer.format = "vcf"
                store = GenoStore(proj, importer)
                if os.path.isfile(proj.name +
                                  "_genotype.DB") and os.path.getsize(
                                      proj.name + "_genotype.DB") > 0:
                    store.load_Genotype_From_SQLite(
                        [proj.name + "_genotype.DB"], proj, importer)
        proj.close()
    except Exception as e:
        env.logger.error(e, exc_info=True)
        sys.exit(1)
