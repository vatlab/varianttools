#!/usr/bin/env python
#
# $File: exporter.py $
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
from threading
from itertools import izip, repeat
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, normalizeVariant, \
    consolidateFieldName, DatabaseEngine

 
class JoinFields:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        try:
            if type(item) == str:
                return item
            else:
                return self.sep.join(item)
        except:
            return str(item)

class IfMulti:
    def __init__(self, ifFunc=None, elseFunc=None):
        if hasattr(ifFunc, '__call__'):
            self.ifFunc = ifFunc.__call__
        else:
            self.ifFunc = ifFunc
        if hasattr(elseFunc, '__call__'):
            self.elseFunc = elseFunc.__call__
        else:
            self.elseFunc = elseFunc

    def __call__(self, item):
        if type(item) == tuple:
            return item[0] if self.ifFunc is None else self.ifFunc(item)
        else:
            return item if self.elseFunc is None else self.elseFunc(item)
        
class JoinRecords:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        try:
            if type(item) == str:
                return item
            else:
                return self.sep.join([str(x) for x in item])
        except:
            return str(item)

class ValueOfNull:
    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return self.val if item in ('', None) else item


class Formatter:
    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, item):
        try:
            return self.fmt.format(item)
        except:
            return str(item)

class Constant:
    def __init__(self, val=''):
        self.val = val

    def __call__(self, item):
        return self.val

class SequentialCollector:
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
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item

MAX_COLUMN = 63
def VariantReader(proj, table, export_by_fields, var_fields, geno_fields,
        export_alt_build, IDs, jobs):
    if jobs == 0 and len(IDs) < MAX_COLUMN:
        # using a single thread
        return EmbeddedVariantReader(proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    elif jobs > 0 and len(IDs) < MAX_COLUMN:
        # using a standalone process to read things and
        # pass information using a pipe
        return StandaloneVariantReader(proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    else:
        # using multiple process to handle more than 1500 samples
        if len(IDs) // MAX_COLUMN + 2 > jobs:
            proj.logger.info('Using {} threads to handle {} samples'.format(len(IDs) // MAX_COLUMN + 2, len(IDs)))
        return MultiVariantReader(proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, max(jobs, len(IDs) // MAX_COLUMN + 2))

class BaseVariantReader:
    def __init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        self.proj = proj
        self.table = table
        self.export_by_fields = export_by_fields
        self.var_fields = var_fields
        self.geno_fields = geno_fields
        self.export_alt_build = export_alt_build
        self.IDs = IDs

    def getQuery(self):
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(self.var_fields), self.export_alt_build)
        if self.geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in self.geno_fields:
                    if fld.lower() in header:
                        select_clause += ', {}_genotype.genotype_{}.{}'.format(self.proj.name, id, fld)
                    else:
                        select_clause += ', NULL'
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        if self.geno_fields:
            for id in self.IDs:
                from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} ON {0}_genotype.genotype_{1}.variant_id = {2}.variant_id '\
                    .format(self.proj.name, id, table)
        # WHERE clause
        where_clause = ''
        # GROUP BY clause
        if self.export_by_fields:
            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.export_by_fields)
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ''
        return 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)

    def getVariantQuery(self):
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(['variant_id'] + self.var_fields), self.export_alt_build)
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
        if self.export_by_fields:
            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.export_by_fields + ',variant_id')
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ' ORDER BY {}.variant_id'.format(self.table)
        return 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)

    def getSampleQuery(self, IDs):
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            'variant_id', False)
        for id in IDs:
            header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
            for fld in self.geno_fields:
                if fld.lower() in header:
                    select_clause += ', {}_genotype.genotype_{}.{}'.format(self.proj.name, id, fld)
                else:
                    select_clause += ', NULL'
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        processed = set()
        for id in IDs:
            from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} ON {0}_genotype.genotype_{1}.variant_id = {2}.variant_id '\
                .format(self.proj.name, id, self.table)
        # WHERE clause
        where_clause = ''
        # GROUP BY clause
        order_clause = ' ORDER BY {}.variant_id'.format(self.table)
        return 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)


class EmbeddedVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        self.proj = proj
        self.logger = proj.logger
        self.var_fields = var_fields
        BaseVariantReader.__init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)

    def records(self):
        self.logger.debug('Running query {}'.format(self.getQuery()))
        cur = self.proj.db.cursor()
        try:
            cur.execute(self.getQuery())
        except Exception as e:
            raise ValueError('Failed to generate output: {}\nIf your project misses one of the following fields {}, you might want to add them to the project (vtools update TABLE INPUT_FILE --var_info FIELDS) or stop exporting them using format parameters (if allowed).'\
                .format(e, ', '.join(self.var_fields)))
        for rec in cur:
            yield rec


class StandaloneVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        BaseVariantReader.__init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)
        self.proj = proj
        self.var_fields = var_fields
        self.reader, w = Pipe(False)
        self.worker = VariantWorker(proj.name, self.getQuery(), w, proj.logger)
        self.worker.start()

    def records(self):
        while True:
            rec = self.reader.recv()
            if rec is None:
                break
            else:
                yield rec

class MultiVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, jobs):
        BaseVariantReader.__init__(self, proj, table, export_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)
        self.proj = proj
        self.var_fields = var_fields
        # the first job for variants
        r, w = Pipe(False)
        block = len(IDs) // (jobs-1) + 1
        p = VariantWorker(proj.name, self.getVariantQuery(), w, proj.logger)
        self.workers = [p]
        self.readers = [r]
        IDs = list(IDs)
        IDs.sort()
        for i in range(jobs - 1):
            r, w = Pipe(False)
            subIDs = IDs[(block*i):(block *(i + 1))]
            p = VariantWorker(proj.name, self.getSampleQuery(subIDs), w, proj.logger)
            self.workers.append(p)
            self.readers.append(r)
        for w in self.workers:
            w.start()

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
        rec = []
        id = None
        last = len(self.readers) - 1
        while True:
            try:
                for idx, reader in enumerate(self.readers):
                    val = reader.recv()
                    if val is None:
                        break
                    if idx == 0:
                        id = val[0]
                    elif id != val[0]:
                        raise ValueError('Read different IDs from multiple processes')
                    rec.extend(val[1:])
                    if idx == last:
                        yield rec
                        rec = []
            except Exception as e:
                print e
        for p in self.workers:
            p.join()

class VariantWorker(threading):
    # this class starts a process and used passed query to read variants
    def __init__(self, dbname, query, output, logger):
        self.dbname = dbname
        self.query = query
        self.output = output
        self.logger = logger
        Process.__init__(self)

    def run(self): 
        db = DatabaseEngine()
        db.connect(self.dbname + '.proj', readonly=True)
        if '{}_genotype.'.format(self.dbname) in self.query:
            db.attach('{}_genotype.DB'.format(self.dbname), '{}_genotype'.format(self.dbname))
        cur = db.cursor()
        cur.execute(self.query)
        for rec in cur:
            self.output.send(rec)
        self.output.send(None)
        db.close()


class Exporter:
    '''A general class for importing variants'''
    def __init__(self, proj, table, filename, samples, format, build, header, jobs, fmt_args):
        self.proj = proj
        self.db = proj.db
        self.jobs = jobs
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
        if samples:
            self.logger.info('{} samples are selected'.format(len(self.IDs)))
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
            filename = self.filename.lower()
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
        self.header = self.getHeader(header, self.IDs)
    
    def getHeader(self, filename, IDs):
        if filename is '':
            return ''
        elif filename == []:
            if IDs:
                # try to get filename for the sample
                cur = self.db.cursor()
                cur.execute('SELECT filename.header FROM filename LEFT JOIN sample ON filename.file_id = sample.file_id WHERE sample.sample_id = {};'.format(self.db.PH),
                    (list(IDs)[0],))
                return cur.fetchone()[0]
            else:
                return ''
        #
        filename = filename[0]
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
        #
        # if file does not exist, try the filename table
        cur = self.db.cursor()
        try:
            cur.execute('SELECT header from filename WHERE filename = {};'.format(self.db.PH), (filename,))
            return cur.fetchone()[0]
        except Exception as e:
            self.logger.debug('Failed to get header: {}'.format(e))
            return ''
        # nothing works
        self.logger.warning('Cannot get header from filename {}. Please check if this file exists in disk or in the sample table (vtools show samples)'.format(filename))
        return ''

    def getAdjFunc(self, code):
        if not code:
            return None
        e = eval(code)
        if hasattr(e, '__iter__'):
            # if there are multiple functors, use a sequence functor
            e = SequentialCollector(e)
        if hasattr(e, '__call__'):
            e = e.__call__
        return e

    def getFormatters(self, indexes, fields):
        #
        # indexes: indexes of data from the input data (returned by sql)
        # fields: field names
        #
        # formatter... : key: value on how to format field(s)
        #
        # Return:
        #    formatters to process all values
        #
        # get first keys for self.formatter
        fmt_first_keys = [x.lower() if ',' in x else x.split(',')[0].lower() for x in self.format.formatter.keys()]
        fmt_keys = [x.lower() for x in self.format.formatter.keys()]
        #
        formatters = []
        i = 0
        while i < len(indexes):
            if fields[i].lower() in fmt_keys:
                # use adj to handle value at index[i]
                formatters.append((self.getAdjFunc(self.format.formatter[fields[i].lower()]), indexes[i]))
                i += 1
            elif fields[i].lower() in fmt_first_keys:
                # start of a multi-field entry
                found = False
                for key,length in [(x,len(x)) for x in fmt_keys if type(x) != str]:
                    if ','.join(fields[i:i+length]).lower() == key.lower():
                        found = True
                        formatters.append(self.getAdjFunc(self.format.formatter[fields[i]]),
                            [indexes[j] for j in range(i, i+length)])
                        i += length
                if not found:
                    formatters.append((None, indexes[i]))
                    i += 1
            else:
                formatters.append((None, indexes[i]))
                i += 1
        return formatters

    def exportData(self):
        '''Export data in specified format'''
        #
        # get all fields
        var_fields = [x.strip() for x in self.format.export_by_fields.split(',')] if self.format.export_by_fields else []
        geno_fields = []
        for col in self.format.columns:
            col_fields = [x.strip() for x in col.field.split(',') if x] if col.field.strip() else []
            if 'GT' in col_fields:
                for fld in col_fields:
                    if fld not in geno_fields:
                        geno_fields.append(fld)
            else:
                for fld in col_fields:
                    if fld not in var_fields:
                        var_fields.append(fld)
        #
        # get indexes of fields
        #
        field_indexes = {}
        for idx,fld in enumerate(var_fields):
            # var_info, field = idx
            field_indexes[(-1, fld.lower())] = idx
        #
        idx = len(var_fields)
        if geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in geno_fields:
                    field_indexes[(id, fld.lower())] = idx
                    idx += 1
        # 
        # how to process each column
        sep = '\t' if self.format.delimiter is None else self.format.delimiter
        formatters = [] # formatters that will be used to produce strings from values
        col_adj = []        # adjust functions to combine values to one column.
        #
        col_idx = 0  # index of things after formatter.
        for col in self.format.columns:
            # indexes to get values for each column
            fields = [x.strip() for x in col.field.split(',') if x] if col.field else []
            if 'GT' in fields:
                for id in self.IDs:
                    col_indexes = []
                    indexes = [field_indexes[(id, x.lower())] for x in fields]
                    fmt = self.getFormatters(indexes, fields)
                    formatters.extend(fmt if fmt else [(None, None)])
                    col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else range(col_idx, col_idx +  len(fmt))])
                    col_idx += max(1, len(fmt))
                    if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                        raise ValueError('Columns with multiple fields must have an adjust function to merge values')
            else:
                indexes = [field_indexes[(-1, x.lower())] for x in fields]
                fmt = self.getFormatters(indexes, fields)
                formatters.extend(fmt if fmt else [(None, None)])
                col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else range(col_idx, col_idx + len(fmt))])
                col_idx += max(1, len(fmt))
                if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                    raise ValueError('Columns with multiple fields must have an adjust function to merge values')

        # needs fmt and adj
        count = 0
        prog = ProgressBar(self.filename, self.db.numOfRows(self.table))
        rec_stack = []
        nFieldBy = len(self.format.export_by_fields.split(','))
        #
        reader = VariantReader(self.proj, self.table, self.format.export_by_fields,
            var_fields, geno_fields, self.export_alt_build, self.IDs, max(self.jobs - 1, 0))
        with open(self.filename, 'w') as output:
            # write header
            if self.header:
                print >> output, self.header.rstrip()
            for idx, raw_rec in enumerate(reader.records()):
                multi_records = False
                try:
                    if nFieldBy != 4:
                        if not rec_stack:
                            rec_stack.append(raw_rec)
                            continue
                        # if the same, wait for the next record
                        elif rec_stack[-1][:nFieldBy] == raw_rec[:nFieldBy]:
                            rec_stack.append(raw_rec)
                            continue
                        elif len(rec_stack) == 1:
                            rec = rec_stack[0]
                            rec_stack = [raw_rec]
                        else:
                            n = len(rec_stack)
                            rec = [tuple([rec_stack[i][x] for i in range(n)]) for x in range(len(raw_rec))]
                            multi_records = True
                            rec_stack = [raw_rec]
                    else:
                        rec = raw_rec
                    # step one: apply formatters 
                    # if there is no fmt, the item must be either empty or a single item
                    #
                    # fmt: single or list
                    # no fmt: single or None
                    #
                    # this is extremely ugly but are we getting any performance gain?
                    if multi_records:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col][0] is None) else str(rec[col][0])) for fmt, col in formatters]
                    else:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col] is None) else str(rec[col])) for fmt, col in formatters]
                    # step two: apply adjusters
                    #
                    # adj: single or list
                    # no adj: must be single
                    columns = [adj(fields[col] if type(col) is int else [fields[x] for x in col]) if adj else fields[col] for adj, col in col_adj]
                    # step three: output columns
                    print >> output, sep.join(columns)
                    count += 1
                except Exception as e:
                    self.logger.debug('Failed to process record {}: {}'.format(rec, e))
                if idx % 10000 == 0:
                    prog.update(idx)
            # the last block
            if rec_stack:
                try:
                    n = len(rec_stack)
                    if n == 1:
                        rec = rec_stack[0]
                    else:
                        rec = [[rec_stack[i][x] for i in range(n)] for x in range(len(raw_rec))]
                    # step one: apply formatters
                    fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) if fmt else ('' if (col is None or rec[col] is None) else str(rec[col])) for fmt, col in formatters]
                    # step two: apply adjusters
                    columns = [adj(fields[col] if type(col) is int else [fields[x] for x in col]) if adj else fields[col] for adj, col in col_adj]
                    # step three: output columns
                    print >> output, sep.join(columns)
                    count += 1
                except Exception as e:
                    self.logger.debug('Failed to process record {}: {}'.format(rec, e))
        prog.done()
        self.logger.info('{} lines are exported'.format(count))



# Functions provided by this script
#
#

def exportArguments(parser):
    parser.add_argument('table', help='''A variant table whose variants will be exported.
        If parameter --samples is specified, only variants belong to one or more of the
        samples will be exported.'''),
    parser.add_argument('filename', help='''Name of output file.'''),
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Samples that will be exported, specified by conditions such as 'aff=1'
            and 'filename like "MG%%"'. Multiple samples could be exported to a
            file if the output format allows. No sample will be outputted if this
            parameter is ignored.''')
    parser.add_argument('--format',
        help='''Format of the exported file. It can be one of the variant tools
            supported file types such as VCF (cf. 'vtools show formats') or a local
            format specification file (with extension .fmt). If unspecified, variant
            tools will try to guess format from file extension. Some formats accept
            additional parameters (cf. 'vtools show format FMT') and allows you to
            export additional or alternative fields.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the exported data. It
            can only be one of the primary (default) of alternative (if exists) reference
            genome of the project.'''),
    parser.add_argument('--header', nargs='*', default='',
        help='''If a valid file is specified, the header (leading comment lines starting with
            #) of this file will become the header of the exported file. If the file does not
            exist, variant tools will retrieve saved header of this file (if filename exists
            in the sample table) or the file from which the first sample is imported (if
            filename is empty and --samples are specified). No header will be used for
            the exported file if this parameter is left unspecified.''')
    parser.add_argument('-j', '--jobs', type=int, default=1,
        help='''Number of processes to export data. Multiple threads will be automatically
            used if there are a large number of samples.''')

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            exporter = Exporter(proj=proj, table=args.table, filename=args.filename,
                samples=' AND '.join(['({})'.format(x) for x in args.samples]), format=args.format,
                build=args.build, header=args.header, jobs=args.jobs,
                fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        sys.exit(e)


