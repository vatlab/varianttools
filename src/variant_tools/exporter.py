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

import datetime
import os
import re
import sys
from argparse import SUPPRESS

from .exporter_reader import VariantReader
from .preprocessor import PlainFormatter, SequentialCollector
from .project import Project, fileFMT
from .utils import (ProgressBar, decodeTableName, encodeTableName, env,
                    splitField)

MAX_COLUMN = 62



class Exporter:
    '''A general class for importing variants'''
    def __init__(self, proj, table, filename, samples, format, build, header,
        jobs, fmt_args):
        self.proj = proj
        self.db = proj.db
        self.jobs = jobs
        #
        # table
        self.table = encodeTableName(table)
        if not self.proj.isVariantTable(self.table):
            raise ValueError('{} is not a valid variant table'.format(decodeTableName(self.table)))
        #
        # filename
        self.filename = filename
        #
        # samples
        self.IDs = self.proj.selectSampleByPhenotype(samples) if samples else []

        self.samples = []
        if samples:
            env.logger.info('Genotypes of {} samples are exported.'.format(len(self.IDs)))
            cur = self.db.cursor()
            for ID in self.IDs:
                cur.execute('SELECT filename, sample_name FROM sample, filename WHERE sample.file_id = filename.file_id AND sample.sample_id = {};'\
                    .format(self.db.PH), (ID,))
                for rec in cur:
                    self.samples.append('{}'.format(rec[1]))
                    env.logger.debug('\t'.join(['{}'.format(x) for x in rec]))
        #
        # build
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                if self.proj.alt_build is not None:
                    env.logger.info('Using primary reference genome {} of the project.'.format(self.build))
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
            filename = self.filename.lower() if self.filename is not None else ''
            if filename.lower().endswith('.vcf'):
                format = 'vcf'
            elif filename.lower().endswith('.csv'):
                format = 'csv'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            self.format = fileFMT(format, fmt_args)
        except Exception as e:
            env.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        if not self.format.columns:
            raise ValueError('Cannot output in format {} because no output column is defined for this format.'.format(format))
        #
        # header
        if header is None:
            self.header = ''
        elif header == ['-']:
            # read from standard input
            env.logger.info('Reading header from standard input')
            self.header = sys.stdin.read()
        else:
            self.header = self.format.delimiter.join(header)

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
        fmt_first_keys = [x.lower() if ',' in x else x.split(',')[0].lower() for x in list(self.format.formatter.keys())]
        fmt_keys = [x.lower() for x in list(self.format.formatter.keys())]
        #
        formatters = []
        i = 0
        while i < len(indexes):
            fld_fmt = [f.fmt for f in self.format.fields if f.name.lower() == fields[i].lower()]
            if not fld_fmt:
                fld_fmt = [f.fmt for f in self.format.other_fields if f.name.lower() == fields[i].lower()]
            # if fmt is defined for the field
            if fld_fmt and fld_fmt[0] is not None:
                formatters.append((self.getAdjFunc(fld_fmt[0]), indexes[i]))
                i += 1
            # if fmt_ is defined in the [formatters] section (deprecated for single field)
            elif fields[i].lower() in fmt_keys:
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

    def exportTfam(self, fname):
        fname = fname.replace('$table', self.table)
        fname += '.tfam' if not fname.endswith('.tfam') else ''
        if os.path.exists(fname):
            os.remove(fname)
        env.logger.info('Sample names are exported to {}'.format(fname))
        # a tfam file is the first 6 columns of a ped file
        # FID, ID, paternal ID, maternal ID, sex, phenotype
        # Will output ID only; other fields will be populated with placeholders
        with open(fname, 'a') as f:
            for i in range(len(self.samples)):
                f.write('\t'.join([str(i), self.samples[i], '0', '0', '0', '-9']) + '\n')
        return True

    def exportData(self):
        '''Export data in specified format'''
        if self.format.additional_exports:
            additionalFiles = [x.strip() for x in self.format.additional_exports.split(',')]
            for item in additionalFiles:
                base, ext = os.path.splitext(item)
                if len(ext) == 0 or len(base) == 0:
                    raise ValueError('Invalid filename specification "{}".'.format(item))
                try:
                    eval('{}'.format('self.export{}(base)'.format(ext[1:].capitalize())))
                except Exception:
                    raise ValueError('Additional export to file *{} is not supported'.format(ext))
        else:
            pass
        # get all fields
        var_fields = [x.strip() for x in self.format.export_by_fields.split(',')] if self.format.export_by_fields else []



        geno_fields = []
        for col in self.format.columns:
            col_fields = splitField(col.field)

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
                if self.proj.store=="sqlite":
                    header = [x.lower() for x in self.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in geno_fields:
                    field_indexes[(id, fld.lower())] = idx
                    idx += 1


        #
        # how to process each column
        sep = '\t' if self.format.delimiter is None else self.format.delimiter
        formatters = [] # formatters that will be used to produce strings from values
        col_adj = []        # adjust functions to combine values to one column.

        default_formatter = PlainFormatter().__call__ if '*' not in list(self.format.formatter.keys()) else self.getAdjFunc(self.format.formatter['*'])
        #
        col_idx = 0  # index of things after formatter.

        for col in self.format.columns:
            # indexes to get values for each column
            fields = splitField(col.field)
            if 'GT' in fields:
                for id in self.IDs:
                    indexes = [field_indexes[(id, x.lower())] for x in fields]
                    fmt = self.getFormatters(indexes, fields)
                    formatters.extend(fmt if fmt else [(None, None)])
                    col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else list(range(col_idx, col_idx +  len(fmt)))])
                    col_idx += max(1, len(fmt))
                    if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                        raise ValueError('Columns with multiple fields must have an adjust function to merge values')
            else:
                indexes = [field_indexes[(-1, x.lower())] for x in fields]
                fmt = self.getFormatters(indexes, fields)
                formatters.extend(fmt if fmt else [(None, None)])
                col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else list(range(col_idx, col_idx + len(fmt)))])
                col_idx += max(1, len(fmt))
                if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                    raise ValueError('Columns with multiple fields must have an adjust function to merge values')

        # needs fmt and adj
        count = 0
        failed_count = 0
        nr = self.db.numOfRows(self.table)
        last_count = 0
        update_after = min(max(100, nr/100), 10000)
        rec_stack = []
        # if export_by_fields is empty
        nFieldBy = len([x for x in self.format.export_by_fields.split(',') if x])
        #
        # cur=self.db.cursor()
        # cur.execute("select group_concat(sample_id) from sample group by HDF5;")
        # res=cur.fetchall()
        # for ids in res:
        #     # IDlist=[]
        #     # groupIDs=[int(id) for id in ids[0].split(",")]
        #     # if len(self.IDs)>0:

        #     #     IDlist=list(set(self.IDs).intersection(set(groupIDs)))
        #     # else:
        #     #     IDlist=groupIDs
        #     self.IDs=ids[0]
        #     if len(self.IDs)>0:

        reader = VariantReader(self.proj, self.table, self.format.export_by_fields, self.format.order_by_fields,
            var_fields, geno_fields, self.export_alt_build, self.IDs, max(self.jobs - 1, 0))
        reader.start()


        prog = ProgressBar(self.filename if self.filename else 'Writing', nr)
        #
        # we cannot use a with statement here because we cannot close sys.stdout
        # perhaps a better solution is available.
        output = open(self.filename, 'w') if self.filename else sys.stdout
        # write header
        if self.header:
            # we can support more variables later
            header_keywords = ['sample_names', 'datetime', 'command']
            # interpolate %(VAR)s with values
            for m in re.finditer('%\((\w+)\)s', self.header):
                if m.group(1) not in header_keywords:
                    env.logger.warning('Unknown variable {} in header specification. '
                        'Only {} are supported'.format(m.group(1), ', '.join(header_keywords)))
            #
            # sample_names
            header = self.header.replace('%(sample_names)s',
                        self.format.delimiter.join(self.samples))
            # datetime
            header = header.replace('%(datetime)s',
                datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
            # command
            header = header.replace('%(command)s',
                env.command_line)
            print(header.rstrip(), file=output)
        global rec_alleles
        for idx, raw_rec in enumerate(reader.records()):
            multi_records = False
            try:
                if nFieldBy != 0:
                    if not rec_stack:
                        rec_stack.append(raw_rec)
                        continue
                    # if the same, wait for the next record
                    elif rec_stack[-1][2:(nFieldBy+2)] == raw_rec[2:(2+nFieldBy)]:
                        # the record needs to be merged even with the same export_by_fields (e.g. chr pos ref with two alts)
                        if rec_stack[-1] != raw_rec:
                            rec_stack.append(raw_rec)
                        continue
                    elif len(rec_stack) == 1:
                        rec_alleles[0] = rec_stack[0][0]
                        rec_alleles[1] = rec_stack[0][1]
                        rec = rec_stack[0][2:]
                        rec_stack = [raw_rec]
                    else:
                        n = len(rec_stack)
                        rec_alleles[0] = [rec_stack[i][0] for i in range(n)]
                        rec_alleles[1] = [rec_stack[i][1] for i in range(n)]
                        rec = [tuple([rec_stack[i][x] for i in range(n)]) for x in range(2, len(raw_rec))]
                        multi_records = True
                        rec_stack = [raw_rec]
                else:
                    rec_alleles[0] = raw_rec[0]
                    rec_alleles[1] = raw_rec[1]
                    rec = raw_rec[2:]

                # step one: apply formatters
                # if there is no fmt, the item must be either empty or a single item
                #
                # fmt: single or list
                # no fmt: single or None
                #
                # this is extremely ugly but are we getting any performance gain?
                if multi_records:
                    try:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col][0] is None) else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'.format(
                                    rec[col] if type(col) is int else [rec[x] for x in col], col, e))
                else:
                    try:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col] is None) else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'.format(
                                    rec[col] if type(col) is int else [rec[x] for x in col], col, e))
                # step two: apply adjusters
                try:
                    #
                    # adj: single or list
                    # no adj: must be single
                    columns = [adj(fields[col] if type(col) is int else [fields[x] for x in col]) if adj else fields[col] for adj, col in col_adj]
                except:
                    # this part just try to get a more precise error message
                    for adj, col in col_adj:
                        try:
                            if type(col) is int:
                                adj(fields[col])
                            elif adj:
                                [fields[x] for x in col]
                            else:
                                fields[col]
                        except:
                            raise ValueError('Failed to apply adjust functor to column {}'.format(col))
                # step three: output columns
                # print(rec,columns)

                line = sep.join(columns)
                output.write(line + '\n')
                count += 1

            except Exception as e:
                env.logger.debug('Failed to process record {}: {}'.format(rec, e))
                failed_count += 1
            if idx - last_count > update_after:
                last_count = idx
                prog.update(idx)
        # the last block
        if rec_stack:
            try:
                n = len(rec_stack)
                multi_records = n > 1
                if n == 1:
                    rec_alleles[0] = rec_stack[0][0]
                    rec_alleles[1] = rec_stack[0][1]
                    rec = rec_stack[0][2:]
                else:
                    rec_alleles[0] = [rec_stack[i][0] for i in range(n)]
                    rec_alleles[1] = [rec_stack[i][1] for i in range(n)]
                    rec = [tuple([rec_stack[i][x] for i in range(n)]) for x in range(2, len(raw_rec))]
                # step one: apply formatters
                if multi_records:
                    try:
                        fields = [fmt(None if col is None else (rec[col] \
                            if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col][0] is None)\
                            else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] \
                                    if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'
                                    .format(rec[col] if type(col) is int else \
                                        [rec[x] for x in col], col, e))
                else:
                    try:
                        fields = [fmt(None if col is None else (rec[col] \
                            if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col] is None)\
                            else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] \
                                    if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'
                                    .format(rec[col] if type(col) is int else \
                                        [rec[x] for x in col], col, e))
                # step two: apply adjusters
                columns = [adj(fields[col] if type(col) is int else \
                    [fields[x] for x in col]) if adj else fields[col] \
                    for adj, col in col_adj]
                # step three: output columns
                line = sep.join(columns)
                output.write(line + '\n')
                count += 1
            except Exception as e:
                env.logger.debug('Failed to process record {}: {}'.format(rec, e))
                failed_count += 1
        if self.filename is not None:
            output.close()
        prog.done()
        env.logger.info('{} lines are exported from variant table {} {}'
            .format(count, decodeTableName(self.table), '' if failed_count == 0 else \
                'with {} failed records'.format(failed_count)))



# Functions provided by this script
#
#

def exportArguments(parser):
    parser.add_argument('table', help='''A variant table whose variants will be exported.
        If parameter --samples is specified, only variants belong to one or more of the
        samples will be exported.''')
    parser.add_argument('filename', nargs='?', help=SUPPRESS)
    parser.add_argument('-o', '--output', nargs='?', help='''Name of output file.
        Output will be written to the standard output if this parameter is left
        unspecified.''')
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Samples that will be exported, specified by conditions such as 'aff=1'
            and 'filename like "MG%%"'. Multiple samples could be exported to a
            file if the output format allows. No sample will be outputted if this
            parameter is ignored.''')
    parser.add_argument('--format',
        help='''Format of the exported file. It can be one of the variant tools
            supported file types such as VCF (cf. 'vtools show formats') or a local
            format specification file (with extension .fmt). Some formats accept
            additional parameters (cf. 'vtools show format FMT') and allows you to
            export additional or alternative fields.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the exported data. It
            can only be one of the primary (default) of alternative (if exists) reference
            genome of the project.''')
    parser.add_argument('--header', nargs='*',
        help='''A complete header or a list of names that will be joined by a
            delimiter specified by the file format to form a header. If a special
            name - is specified, the header will be read from the standard input,
            which is the preferred way to specify large multi-line headers (e.g.
            cat myheader | vtools export --header -). Strings in the form of
            %%(VAR)s will be interpolated to values of variable VAR, which can be
            "sample_names" for list of sample names, "datetime" for current date
            and time, and "command" for the command used to create output.''')
    parser.add_argument('-j', '--jobs', type=int, default=1,
        help='''Number of processes to export data. Multiple threads will be automatically
            used if there are a large number of samples.''')

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if proj.store=="sqlite":
                proj.db.attach(proj.name + '_genotype')
            if args.filename:
                raise ValueError('Specifying output filename without --output is deprecated.')
            exporter = Exporter(proj=proj, table=args.table, filename=args.output,
                samples=' AND '.join(['({})'.format(x) for x in args.samples]),
                format=args.format, build=args.build, header=args.header,
                jobs=args.jobs, fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
