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
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, normalizeVariant, \
    consolidateFieldName


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
    def __init__(self, proj, table, filename, samples, format, build, header, fmt_args):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        #
        # table
        if not self.proj.isVariantTable(table):
            raise ValueError('{} is not a valid variant table'.format(table))
        self.table = table
        #
        # filename
        self.filename = filename
        #
        # samples
        self.IDs = self.proj.selectSampleByPhenotype(samples) if samples else []
        # 
        # build
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                if self.proj.alt_build is not None:
                    self.logger.info('Using primary reference genome {} of the project.'.format(self.build))
            else:
                raise ValueError('This project does not have any data to export.')
            self.export_alt_build = False
        elif build in [self.proj.build, self.proj.alt_build]:
            self.build = build
            self.export_alt_build = build == self.proj.alt_build
        else:
            raise ValueError('Invalid parameter for --build, which should be either the primary ({}) or the alternative ({}) reference genome of the project.'\
                .format(self.proj.build, self.proj.alt_build))
        #
        # format
        if not format:
            filename = self.files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            self.format = fileFMT(format, fmt_args)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        if not self.format.columns:
            raise ValueError('Cannot output in format {} because no output column is defined for this format.'.format(format)) 
        #
        # header
        self.header = self.getHeader(header)
    
    def getHeader(self, filename):
        if not filename:
            return ''
        #
        if os.path.isfile(filename):
            if filename.lower().endswith('.gz'):
                input = gzip.open(filename, 'rb')
            else:
                input = open(filename, 'rb')
            header = ''
            for line in input:
                line = line.decode()
                if line.startswith('#'):
                    header += line
                else:
                    return header
        # if file does not exist, try the filename table
        cur = self.db.cursor()
        try:
            cur.execute('SELECT header from filename WHERE filename = {};'.format(self.db.PH), (filename,))
            return cur.fetchone()[0]
        except Exception as e:
            self.logger.debug(e)
            return ''
        # nothing works
        self.logger.warning('Cannot get header from filename {}. Please check if this file exists in disk or in the sample table (vtools show samples)'.format(filename))
        return ''

    def exportData(self):
        '''Export data in specified format'''
        #
        # get all fields
        fields = []
        for col in self.format.columns:
            fields.extend([x.strip() for x in col.field.split(',')] if col.field else [])
        output_fields = list(set(fields))
        #
        # fields for each column
        col_fields = []
        for col in self.format.columns:
            f = [x.strip() for x in col.field.split(',')] if col.field else []
            col_fields.append([output_fields.index(x) for x in f])
        #
        # fields
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(output_fields), self.export_alt_build)
        #
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        # WHERE clause
        where_clause = ''
        # GROUP BY clause
        group_clause = ''
        #if args.group_by:
        #    group_fields, tmp = consolidateFieldName(proj, table, ','.join(args.group_by))
        #    group_clause = ' GROUP BY {}'.format(group_fields)
        query = 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, group_clause)
        self.logger.debug('Running query {}'.format(query))
        #
        # get variant and their info
        cur = self.db.cursor()
        cur.execute(query)
        with open(self.filename, 'w') as output:
            # write header
            if self.header:
                print >> output, self.header.rstrip()
            for rec in cur:
                print>> output, ','.join([' '.join([str(rec[x]) for x in col]) for col in col_fields])
            return


# Functions provided by this script
#
#

def exportArguments(parser):
    parser.add_argument('table', help='''A variant table whose variants will be exported.
        If parameter --samples is specified, only variants belong to one or more of the
        samples will be exported.'''),
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
    parser.add_argument('--header',
        help='''If specified, the header (leading comment lines starting with #) of this
            file will become the header of the exported file. Because vtools saves header
            of imported files, the original header will be used even if the input files
            have been removed or renamed.''')

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            exporter = Exporter(proj=proj, table=args.table, filename=args.filename,
                samples=args.samples, format=args.format, build=args.build, header=args.header,
                fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        sys.exit(e)


