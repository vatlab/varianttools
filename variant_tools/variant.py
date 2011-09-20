#!/usr/bin/env python
#
# $File: variant.py $
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

import sys
from .project import Project
from .utils import ProgressBar, consolidateFieldName, typeOfValues, lineCount
from .sample import Sample


def outputArguments(parser):
    parser.add_argument('table', help='''variants to output.''')
    parser.add_argument('fields', nargs='+',
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')

def generalOutputArguments(parser):
    grp = parser.add_argument_group('Output options')
    grp.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter, default to tab, a popular alternative is ',' for csv output''')
    grp.add_argument('--na', default='NA',
        help='Output string for missing value')
    grp.add_argument('-l', '--limit', default=-1, type=int,
        help='''Number of record to display. Default to all record.''')
    grp.add_argument('--build',
        help='''Output reference genome. If set to alternative build, chr and pos
            in the fields will be replaced by alt_chr and alt_pos''')
    grp.add_argument('-g', '--group_by', nargs='*',
        help='''Group output by fields. This option is useful for aggregation output
            where summary statistics are grouped by one or more fields.''')

def outputVariants(proj, table, output_fields, args, query=None, reverse=False):
    '''Output selected fields'''
    # output
    out = sys.stdout
    #
    # table
    if not proj.isVariantTable(table):
        raise ValueError('Variant table {} does not exist.'.format(table))
    #
    # fields
    select_clause, fields = consolidateFieldName(proj, table, ','.join(output_fields),
        args.build and args.build == proj.alt_build)
    #
    # FROM clause
    from_clause = 'FROM {} '.format(table)
    fields_info = sum([proj.linkFieldToTable(x, table) for x in fields], [])
    #
    processed = set()
    for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
        if (tbl.lower(), conn.lower()) not in processed:
            from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
            processed.add((tbl.lower(), conn.lower()))
    # WHERE clause
    where_clause = ''
    if query is not None:
        # FIXME: if the query has a simple where clause, we should use that directly.
        where_clause = ' WHERE {}.variant_id {} IN ({})'.format(table, 'NOT' if reverse else '', query)
    # GROUP BY clause
    group_clause = ''
    if args.group_by:
        group_fields, tmp = consolidateFieldName(proj, table, ','.join(args.group_by))
        group_clause = ' GROUP BY {}'.format(group_fields)
    # LIMIT clause
    limit_clause = '' if args.limit < 0 else ' LIMIT 0,{}'.format(args.limit)
    query = 'SELECT {} {} {} {} {};'.format(select_clause, from_clause, where_clause, group_clause, limit_clause)
    proj.logger.debug('Running query {}'.format(query))
    # if output to a file
    cur = proj.db.cursor()
    cur.execute(query)
    proj.logger.info('Writing output')
    for count, rec in enumerate(cur):
        out.write(args.delimiter.join([args.na if x is None else str(x) for x in rec]) + '\n')


def output(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            outputVariants(proj, args.table, args.fields, args)
    except Exception as e:
        sys.exit(e) 

def selectArguments(parser):
    parser.add_argument('from_table', help='''Source variant table.''')
    parser.add_argument('condition', nargs='*', default=[],
        help='''Conditions by which variants are selected. Multiple arguments are
            automatically joined by 'AND' so 'OR' conditions should be provided by
            a single argument with conditions joined by 'OR'. If unspecified, all
            variants (except those excluded by parameter --samples) will be selected.''')
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-t', '--to_table',
        help='''Destination variant table. ''')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant, which is a shortcut to '--output count(1)'.''')
    grp.add_argument('-o', '--output', nargs='*', default=[],
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')


def select(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # table?
            if not proj.isVariantTable(args.from_table):
                raise ValueError('Variant table {} does not exist.'.format(args.from_table))
            if args.to_table == 'variant':
                raise ValueError('Cannot overwrite master variant table. Please choose another name for the variant table')
            if not args.to_table and not args.output and not args.count:
                proj.logger.warning('Neither --to_table and --output/--count is specified. Nothing to do.')
                return
            if len(args.condition) > 0:    
                # fields? We need to replace things like sift_score to dbNSFP.sift_score
                condition, fields = consolidateFieldName(proj, args.from_table, ' AND '.join(['({})'.format(x) for x in args.condition]))
                # 
                fields_info = sum([proj.linkFieldToTable(x, args.from_table) for x in fields], [])
                # WHERE clause: () is important because OR in condition might go beyond condition
                where_clause = ' WHERE ({})'.format(condition)
                # 
                # FROM clause
                from_clause = 'FROM {} '.format(args.from_table)
                # avoid duplicate
                processed = set()
                for table, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                    if (table.lower(), conn) not in processed:
                        from_clause += ', {} '.format(table)
                        where_clause += ' AND ({}) '.format(conn)
                        processed.add((table.lower(), conn))
            else:
                # select all variants
                where_clause = ' WHERE 1 '
                from_clause = 'FROM {} '.format(args.from_table)
            # if limiting to specified samples
            if args.samples:
                # we save genotype in a separate database to keep the main project size tolerable.
                proj.db.attach(proj.name + '_genotype')
                IDs = proj.selectSampleByPhenotype(' AND '.join(args.samples))
                if len(IDs) == 0:
                    proj.logger.warning('No sample is selected by condition: {}'.format(' AND '.join(args.samples)))
                    # nothing will be selected
                    where_clause += ' AND 0'
                #
                # This special case does not hold because sometimes variants are imported without sample information.
                #
                #elif len(IDs) == proj.db.numOfRows('sample'):
                #    proj.logger.info('All {} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                #    # we do not have to add anything to where_clause
                else:
                    proj.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                    where_clause += ' AND ({}.variant_id IN ({}))'.format(
                        args.from_table, 
                        '\nUNION '.join(['SELECT variant_id FROM {}_genotype.sample_variant_{}'.format(proj.name, id) for id in IDs])) 
            #
            # we are treating different outcomes different, for better performance
            #
            # NOTE: count/output do not co-exist
            #
            # case:       count,  to_table, output
            # 1:           Y       N       N   <- select count(variant)
            # 2:           N       Y       Y   <- generate table, then count and output
            #              Y       Y       N
            #              N       Y       N
            # 3:           N       N       Y   <- direct output
            #
            # Others:      N       N       N   <- do nothing.
            #              Y       Y       Y   <- not allowed
            #              Y       N       Y
            #
            # case 1: simple count.
            if args.count and not args.to_table and not args.output:
                query = 'SELECT COUNT(DISTINCT {}.variant_id) {} {};'.format(args.from_table,
                    from_clause, where_clause)
                proj.logger.debug('Running query {}'.format(query))
                proj.db.startProgress('Counting variants')
                cur = proj.db.cursor()
                cur.execute(query)
                count = cur.fetchone()[0]
                # exclude ...
                if reverse:
                    count = proj.db.numOfRows(args.from_table) - int(count)
                proj.db.stopProgress()
                #
                print(count)
            # case 2: to table
            elif args.to_table:
                if proj.db.hasTable(args.to_table):
                    new_table = proj.db.backupTable(args.to_table)
                    proj.logger.warning('Existing table {} is renamed to {}.'.format(args.to_table, new_table))
                #
                proj.createVariantTable(args.to_table)
                if not reverse:
                    query = 'INSERT INTO {0} SELECT DISTINCT {1}.variant_id {2} {3};'.format(args.to_table, args.from_table,
                        from_clause, where_clause)
                else:
                    query = 'INSERT INTO {0} SELECT DISTINCT {1}.variant_id FROM {1} WHERE {1}.variant_id NOT IN (SELECT {1}.variant_id {2} {3});'.\
                        format(args.to_table, args.from_table, from_clause, where_clause)
                proj.logger.debug('Running query {}'.format(query))
                #
                cur = proj.db.cursor()
                proj.db.startProgress('Running')
                cur.execute(query)
                proj.db.stopProgress()
                proj.db.commit()
                #
                count = proj.db.numOfRows(args.to_table)
                proj.logger.info('{} variants selected.'.format(count))
                if args.output:
                    outputVariants(proj, args.to_table, args.output, args)
                if args.count:
                    print(count)
            # case 3: output, but do not write to table, and not count
            elif args.output: 
                query = 'SELECT DISTINCT {}.variant_id {} {}'.format(args.from_table,
                    from_clause, where_clause)
                outputVariants(proj, args.from_table, args.output, args, query, reverse)
    except Exception as e:
        sys.exit(e) 

def excludeArguments(parser):
    parser.add_argument('from_table', help='''Source variant table.''')
    parser.add_argument('condition', nargs='*', default=[],
        help='''Conditions by which variants are execluded. Multiple arguments are
            automatically joined by 'AND' so 'OR' conditions should be provided by
            a single argument with conditions joined by 'OR'. If unspecified, all
            variants (except those excluded by parameter --samples) will be excluded.''')
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-t', '--to_table', 
        help='''Destination variant table.''')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant, which is a shortcut to '--output count(1)'.''')
    grp.add_argument('-o', '--output', nargs='*', default=[],
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')

def exclude(args):
    select(args, reverse=True)

def compareArguments(parser):
    parser.add_argument('table_A', help='''variant table A.''')
    parser.add_argument('table_B', help='''variant table B.''')
    parser.add_argument('--A_diff_B',  
        help='''Save variants in A but not in B to a table''')
    parser.add_argument('--B_diff_A',  
        help='''Save variants in B but not in A to a table''')
    parser.add_argument('--A_and_B',  
        help='''Save variants in both tables A and B.''')
    parser.add_argument('--A_or_B',  
        help='''Save variants in both tables A or B.''')
    parser.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant that are only in A, only in B, in both A and B, and in A or B.''')

def compare(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # table?
            if not proj.isVariantTable(args.table_A):
                raise ValueError('Variant table {} does not exist.'.format(args.table_A))
            if not proj.isVariantTable(args.table_B):
                raise ValueError('Variant table {} does not exist.'.format(args.table_B))
            if not args.A_diff_B and not args.B_diff_A and not args.A_and_B and not args.A_or_B and not args.count:
                proj.logger.warning('No action parameter is specified. Nothing to do.')
                return
            #
            cur = proj.db.cursor()
            # read variants in table_A
            proj.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(args.table_A), args.table_A))
            cur.execute('SELECT variant_id from {};'.format(args.table_A))
            variant_A = set([x[0] for x in cur.fetchall()])
            # read variants in table_B
            proj.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(args.table_B), args.table_B))
            cur.execute('SELECT variant_id from {};'.format(args.table_B))
            variant_B = set([x[0] for x in cur.fetchall()])
            #
            # output?
            if args.count:
                proj.logger.info('Output number of variants in A but not B, B but not A, A and B, and A or B')
                print('{}\t{}\t{}\t{}'.format(len(variant_A - variant_B), 
                    len(variant_B - variant_A),
                    len(variant_A & variant_B),
                    len(variant_A | variant_B)
                    ))
            #
            for var, table in [
                    (set() if args.A_diff_B is None else variant_A - variant_B, args.A_diff_B), 
                    (set() if args.B_diff_A is None else variant_B - variant_A, args.B_diff_A), 
                    (set() if args.A_and_B is None else variant_A & variant_B, args.A_and_B), 
                    (set() if args.A_or_B is None else variant_A | variant_B, args.A_or_B)]:
                if table is None:
                    continue
                if table == 'variant':
                    raise ValueError('Cannot overwrite master variant table')
                if proj.db.hasTable(table):
                    new_table = proj.db.backupTable(table)
                    proj.logger.warning('Existing table {} is renamed to {}.'.format(table, new_table))
                proj.createVariantTable(table)
                prog = ProgressBar('Writing to ' + table, len(var))
                query = 'INSERT INTO {} VALUES ({});'.format(table, proj.db.PH)
                # sort var so that variant_id will be in order, which might
                # improve database performance
                for count,id in enumerate(sorted(var)):
                    cur.execute(query, (id,))
                    if count % proj.db.batch == 0:
                        proj.db.commit()
                        prog.update(count)
                proj.db.commit()
                prog.done()
    except Exception as e:
        sys.exit(e) 

