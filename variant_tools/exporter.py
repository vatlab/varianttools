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
from itertools import izip, repeat
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, normalizeVariant, \
    consolidateFieldName

 
class JoinFields:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        try:
            return self.sep.join(item)
        except:
            return str(item)

class ValueOfNull:
    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return self.val if item is None else item


class Constant:
    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return self.val

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
        if filename is None:
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

    def exportData(self):
        '''Export data in specified format'''
        #
        # get all fields
        var_fields = []
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
        # We are going to query:
        #
        #    var_fields   geno_fields * IDs (only if sample_variant has the geno_field)
        #
        field_indexes = {}
        for idx,fld in enumerate(var_fields):
            # var_info, field = idx
            field_indexes[(-1, fld.lower())] = idx
        #
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(var_fields), self.export_alt_build)
        #
        idx = len(var_fields)
        if geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.db.getHeaders('{}_genotype.sample_variant_{}'.format(self.proj.name, id))]
                for fld in geno_fields:
                    if fld.lower() in header:
                        select_clause += ', {}_genotype.sample_variant_{}.{}'.format(self.proj.name, id, fld)
                    else:
                        select_clause += ', NULL'
                    field_indexes[(id, fld.lower())] = idx
                    idx += 1
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        if geno_fields:
            for id in self.IDs:
                from_clause += ' LEFT OUTER JOIN {0}_genotype.sample_variant_{1} ON {0}_genotype.sample_variant_{1}.variant_id = {2}.variant_id '\
                    .format(self.proj.name, id, self.table)
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
        # how to process each column
        sep = '\t' if self.format.delimiter is None else self.format.delimiter
        col_formatters = [] # formatters that will be used to produce strings from values
        col_adj = []        # adjust functions to combine values to one column.
        #
        # get first keys for self.formatter
        fmt_first_keys = [x.lower() if ',' in x else x.split(',')[0].lower() for x in self.format.formatter.keys()]
        fmt_keys = [x.lower() for x in self.format.formatter.keys()]
        #
        col_idx = 0  # index of things after formatter.
        for col in self.format.columns:
            #
            # indexes to get values for each column
            fields = [x.strip() for x in col.field.split(',') if x] if col.field else []
            if 'GT' in fields:
                for id in self.IDs:
                    col_indexes = []
                    indexes = [field_indexes[(id, x.lower())] for x in fields]
                    i = 0
                    while i < len(indexes):
                        if fields[i].lower() in fmt_keys:
                            # use adj to handle value at index[i]
                            col_formatters.append((self.getAdjFunc(self.format.formatter[fields[i].lower()]), [indexes[i]]))
                            col_indexes.append(col_idx)
                            col_idx += 1
                            i += 1
                        elif fields[i].lower() in fmt_first_keys:
                            # start of a multi-field entry
                            found = False
                            for key,length in [(x,len(x)) for x in fmt_keys if type(x) != str]:
                                if ','.join(fields[i:i+length]).lower() == key.lower():
                                    found = True
                                    col_formatters.append(self.getAdjFunc(self.format.formatter[fields[i]]),
                                        [indexes[j] for j in range(i, i+length)])
                                    i += length
                                    col_idx += 1
                                    col_indexes.append(col_idx)
                            if not found:
                                col_formatters.append((None, [indexes[i]]))
                                col_indexes.append(col_idx)
                                col_idx += 1
                                i += 1
                        else:
                            col_formatters.append((None, [indexes[i]]))
                            col_indexes.append(col_idx)
                            col_idx += 1
                            i += 1
                    # adjust function/functors
                    col_adj.append([self.getAdjFunc(col.adj), col_indexes])
            else:
                indexes = [field_indexes[(-1, x.lower())] for x in fields]
                col_indexes = []
                i = 0
                while i < len(indexes):
                    if fields[i].lower() in self.format.formatter:
                        # use adj to handle value at index[i]
                        col_formatters.append((self.getAdjFunc(self.format.formatter[fields[i]]), [indexes[i]]))
                        col_indexes.append(col_idx)
                        col_idx += 1
                        i += 1
                    elif fields[i].lower() in fmt_keys:
                        # start of a multi-field entry
                        found = False
                        for key,length in [(x,len(x)) for x in fmt_keys if type(x) != str]:
                            if ','.join(fields[i:i+length]).lower() == key.lower():
                                found = True
                                col_formatters.append(self.getAdjFunc(self.format.formatter[fields[i]]),
                                    [indexes[j] for j in range(i, i+length)])
                                i += length
                                col_idx += 1
                                col_indexes.append(col_idx)
                        if not found:
                            col_formatters.append((None, [indexes[i]]))
                            col_indexes.append(col_idx)
                            col_idx += 1
                            i += 1
                    else:
                        col_formatters.append((None, [indexes[i]]))
                        col_indexes.append(col_idx)
                        col_idx += 1
                        i += 1
                # adjust function/functors
                col_adj.append([self.getAdjFunc(col.adj), col_indexes])
        #print 'FORMA' , [x[1] for x in col_formatters]
        #print 'COL  ' , [x[1] for x in col_adj]
        #
        # get variant and their info
        cur = self.db.cursor()
        try:
            cur.execute(query)
        except Exception as e:
            raise ValueError('Failed to generate output: {}\nIf your project misses one of the following fields {}, you might want to add them to the project (vtools update TABLE INPUT_FILE --var_info FIELDS) or stop exporting them using format parameters (if allowed).'\
                .format(e, ', '.join(var_fields)))
        prog = ProgressBar(self.filename)
        with open(self.filename, 'w') as output:
            # write header
            if self.header:
                print >> output, self.header.rstrip()
            for idx, rec in enumerate(cur):
                try:
                    # step one: apply formatters
                    fields = [fmt([rec[x] for x in col]) if fmt else str(rec[col[0]]) for fmt, col in col_formatters]
                    # step two: apply adjusters
                    columns = [adj([fields[x] for x in col]) if adj else fields[col[0]] for adj, col in col_adj]
                    # step three: output columns
                    print >> output, sep.join(columns)
                except Exception as e:
                    self.logger.debug('Failed to process record {}: {}'.format(rec, e))
                if idx % 10000 == 0:
                    prog.update(idx)
        prog.done()



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
    parser.add_argument('--header', nargs='*', default='',
        help='''If a valid file is specified, the header (leading comment lines starting with
            #) of this file will become the header of the exported file. If the file does not
            exist, variant tools will retrieve saved header of this file (if filename exists
            in the sample table) or the file from which the first sample is imported (if
            filename is empty and --samples are specified). No header will be used for
            the exported file is this parameter is left unspecified.''')

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            exporter = Exporter(proj=proj, table=args.table, filename=args.filename,
                samples=' AND '.join(args.samples), format=args.format, build=args.build, header=args.header,
                fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        sys.exit(e)


