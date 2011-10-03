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
        # how to process each column
        sep = '\t' if self.format.delimiter is None else self.format.delimiter
        col_fields = []
        for col in self.format.columns:
            f = [x.strip() for x in col.field.split(',')] if col.field else []
            e = eval(col.export_adj) if col.export_adj else None
            if hasattr(e, '__iter__'):
                # if there are multiple functors, use a sequence functor
                e = SequentialCollector(e)
            if hasattr(e, '__call__'):
                e = e.__call__
            col_fields.append((tuple([output_fields.index(x) for x in f]), e))
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
        prog = ProgressBar(self.filename, cur.rowcount)
        with open(self.filename, 'w') as output:
            # write header
            if self.header:
                print >> output, self.header.rstrip()
            for idx, rec in enumerate(cur):
                try:
                    print >> output, sep.join([str(adj([str(rec[x]) for x in col])) if adj else ''.join([str(rec[x]) for x in col]) \
                        for col,adj in col_fields])
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


