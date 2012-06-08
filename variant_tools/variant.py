#!/usr/bin/env python
#
# $File: variant.py $
# $LastChangedDate$
# $Rev$
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

import sys
import re
from .project import Project
from .utils import ProgressBar, consolidateFieldName, typeOfValues, lineCount, delayedAction
from .phenotype import Sample


def outputArguments(parser):
    parser.add_argument('table', help='''variants to output.''')
    parser.add_argument('fields', nargs='+',
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')

def generalOutputArguments(parser):
    grp = parser.add_argument_group('Output options')
    grp.add_argument('--header', nargs='*',
        help='''A complete header or a list of names that will be joined by a delimiter
            (parameter --delimiter). If a special name - is specified, the header will
            be read from the standard input, which is the preferred way to specify large
            multi-line headers (e.g. cat myheader | vtools export --header -). If this
            parameter is given without parameter, a default header will be derived from
            field names.'''),
    grp.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter, default to tab, a popular alternative is ',' for csv output''')
    grp.add_argument('--na', default='NA',
        help='Output string for missing value')
    grp.add_argument('-l', '--limit', default=-1, type=int,
        help='''Number of record to display. Default to all record.''')
    grp.add_argument('--build',
        help='''Output reference genome. If set to alternative build, chr and pos
            in the fields will be replaced by alt_chr and alt_pos''')
    grp.add_argument('-g', '--group_by', nargs='*', metavar='FIELD',
        help='''Group output by fields. This option is useful for aggregation output
            where summary statistics are grouped by one or more fields.''')
    grp.add_argument('--order_by', nargs='*', metavar='FIELD',
        help='''Order output by specified fields in ascending order.''')

def outputVariants(proj, table, output_fields, args, query=None, reverse=False):
    '''Output selected fields'''
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
    where_conditions = []
    fields_info = sum([proj.linkFieldToTable(x, table) for x in fields], [])
    #
    processed = set()
    # the normal annotation databases that are 'LEFT OUTER JOIN'
    for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
        if (tbl.lower(), conn.lower()) not in processed and '.__' not in tbl:
            from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
            processed.add((tbl.lower(), conn.lower()))
    # temporary connection tables are appended as WHERE clause.
    #for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
    #    if (tbl.lower(), conn.lower()) not in processed and '.__' in tbl:
    #        from_clause += ' , {}'.format(tbl)
    #        where_conditions.append(conn)
    #        processed.add((tbl.lower(), conn.lower()))
    # WHERE clause
    if query is not None:
        # FIXME: if the query has a simple where clause, we should use that directly.
        where_conditions.append('{}.variant_id {} IN ({})'.format(table, 'NOT' if reverse else '', query))
    where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
    # GROUP BY clause
    group_clause = ''
    if args.group_by:
        group_fields, tmp = consolidateFieldName(proj, table, ','.join(args.group_by))
        group_clause = ' GROUP BY {}'.format(group_fields)
    order_clause = ''
    if args.order_by:
        order_fields, tmp = consolidateFieldName(proj, table, ','.join(args.order_by))
        order_clause = ' ORDER BY {}'.format(order_fields)
    # LIMIT clause
    limit_clause = '' if args.limit < 0 else ' LIMIT 0,{}'.format(args.limit)
    query = 'SELECT {} {} {} {} {} {};'.format(select_clause, from_clause, where_clause, group_clause, order_clause, limit_clause)
    proj.logger.debug('Running query {}'.format(query))
    # if output to a file
    cur = proj.db.cursor()
    cur.execute(query)
    if args.header is not None:
        if len(args.header) == 0:
            # if no real header is given, use output_fields, but replace things like (, ), and , to space
            sys.stdout.write(args.delimiter.join([re.sub('[\W]+', ' ', x) for x in output_fields]) + '\n')
        elif args.header == ['-']:
            proj.logger.info('Reading header from standard input')
            sys.stdout.write(sys.stdin.read().rstrip() + '\n')
        else:
            # other wise, use the user-provided header
            if len(args.header) != len(output_fields):
                proj.logger.warning('User-provided header ({}) does not match number of fields ({})'.format(len(args.header), len(output_fields)))
            sys.stdout.write(args.delimiter.join(args.header) + '\n')
    # output with a warning to potentially duplicated lines
    has_duplicate = False
    last_line = None
    for count, rec in enumerate(cur):
        line = args.delimiter.join([args.na if x is None else str(x) for x in rec]) + '\n'
        if not has_duplicate:
            if line == last_line:
                has_duplicate = True
            else:
                last_line = line
        sys.stdout.write(line)
    if has_duplicate:
        proj.logger.warning('Duplicate records are outputted. This may due to '
            'multiple entries of variants in the annotation database. Piping standard output '
            'to a uniq command ("vtools output ... | uniq") will remove these entries.')


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
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-t', '--to_table', nargs='*', metavar=('TABLE', 'DESC'),
        help='''Destination variant table. ''')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant, which is a shortcut to '--output count(1)'.''')
    grp.add_argument('-o', '--output', nargs='*', metavar='FIELDS', default=[],
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')


def select(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # separate table and message
            if args.to_table:
                if len(args.to_table) > 2:
                    raise ValueError('Only a table name and an optional message is allowed for parameter to_table')
                args.table_desc = args.to_table[1] if len(args.to_table) == 2 else ''
                args.to_table = args.to_table[0]
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
                for field in fields:
                    # indexing fields in annotation databases?
                    try:
                        # if table name is specified
                        db, fld = field.split('.', 1)
                        annoDB = [x for x in proj.annoDB if x.name.lower() == db.lower()][0]
                    except:
                        continue
                    #
                    # STOP automatic indexing fields used in vtools select because creating indexes can be 
                    # very time-consuming for sqlite databases, and the performance benefit is uncertain.
                    #
                    # db is one of the annotation database but fld has already been indexed
                    #
                    #if fld.lower() in [x.lower() for x in annoDB.linked_by] or \
                    #    (annoDB.build is not None and fld.lower() in [x.lower() for x in annoDB.build]) or \
                    #    (annoDB.alt_build is not None and fld.lower() in [x.lower() for x in annoDB.alt_build]) or \
                    #    (fld.lower() not in [x.lower() for x in proj.db.getHeaders('{0}.{0}'.format(db))]) or \
                    #    proj.db.hasIndex('{0}.{0}_{1}'.format(db, fld)):
                    #    continue
                    #
                    #s = delayedAction(proj.logger.info, 'Indexing {}'.format(field))
                    #cur = proj.db.cursor()
                    #try:
                    #    query = 'CREATE INDEX IF NOT EXISTS {0}.{0}_{1} ON {0} ({1} ASC);'.format(db, fld)
                    #    proj.logger.debug(query)
                    #    cur.execute(query)
                    #except Exception as e:
                    #    proj.logger.debug('Failed to create index: {}'.format(e))
                    #del s
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
                IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.samples]))
                if len(IDs) == 0:
                    proj.logger.warning('No sample is selected by condition: {}'.format(' AND '.join(['({})'.format(x) for x in args.samples])))
                    # nothing will be selected
                    where_clause += ' AND 0'
                #
                # This special case does not hold because sometimes variants are imported without sample information.
                #
                #elif len(IDs) == proj.db.numOfRows('sample'):
                #    proj.logger.info('All {} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                #    # we do not have to add anything to where_clause
                elif len(IDs) < 50:  
                    # we allow 14 tables in other 'union' or from condition...
                    proj.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                    where_clause += ' AND ({}.variant_id IN ({}))'.format(
                        args.from_table, 
                        '\nUNION '.join(['SELECT variant_id FROM {}_genotype.genotype_{}'.format(proj.name, id) for id in IDs])) 
                else:
                    # we have to create a temporary table and select variants sample by sample
                    # this could be done in parallel if there are a large number of samples, but that needs a lot more
                    # code, and perhaps RAM
                    proj.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                    cur = proj.db.cursor()
                    BLOCK_SIZE = 64
                    NUM_BLOCKS = len(IDs) // BLOCK_SIZE + 1
                    myIDs = list(IDs)
                    myIDs.sort()
                    merged_table = None
                    prog = ProgressBar('Collecting sample variants', len(IDs)) if NUM_BLOCKS > 1 else None
                    count = 0
                    for i in range(NUM_BLOCKS):
                        # step 1: create a table that holds all
                        block_IDs = myIDs[(i*BLOCK_SIZE):((i+1)*BLOCK_SIZE)]
                        if len(block_IDs) == 0:
                            continue
                        merged_table = '__variants_from_samples_{}'.format(i)
                        query = 'CREATE TEMPORARY TABLE {} (variant_id INT);'.format(merged_table)
                        # proj.logger.debug(query)
                        cur.execute(query)
                        query = 'INSERT INTO {} {} {};'.format(merged_table,
                            # also merge last batch
                            '\nSELECT variant_id FROM __variants_from_samples_{} UNION '.format(i-1) if i > 1 else '',
                            '\nUNION '.join(['SELECT variant_id FROM {}_genotype.genotype_{}'.format(proj.name, id) for id in block_IDs]))
                        #proj.logger.debug(query)
                        cur.execute(query)
                        if i > 1:
                            # remove last batch
                            query = 'DROP TABLE __variants_from_samples_{}'.format(i-1)
                            #proj.logger.debug(query)
                            cur.execute(query)
                        count += len(block_IDs)
                        if prog:
                            prog.update(count)
                    if prog:
                        prog.done()
                    where_clause += ' AND ({}.variant_id IN (SELECT variant_id FROM {}))'.format(
                        args.from_table, merged_table)
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
                proj.describeTable(args.to_table, args.table_desc, True, True)
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
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-t', '--to_table', 
        help='''Destination variant table.''')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant, which is a shortcut to '--output count(1)'.''')
    grp.add_argument('-o', '--output', nargs='*', metavar='FIELDS', default=[],
        help='''A list of fields that will be outputed. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')

def exclude(args):
    select(args, reverse=True)

def compareArguments(parser):
    parser.add_argument('tables', nargs='+', help='''variant tables to compare.''')
    parser.add_argument('--A_diff_B', metavar= 'TABLE',
        help='''deprecated, use --difference instead.''')
    parser.add_argument('--B_diff_A', metavar= 'TABLE',
        help='''deprecated, use --difference instead.''')
    parser.add_argument('--A_and_B', metavar= 'TABLE',
        help='''deprecated, use --intersect instead.''')
    parser.add_argument('--A_or_B', metavar= 'TABLE',
        help='''deprecated, use --union instead.''')
    parser.add_argument('--union', metavar=('TABLE', 'DESC'), nargs='*', default=[],
        help='''Save variants in any of the tables (T1 | T2 | T3 ...) to TABLE if a name
             is specified. An optional message can be added to describe the table.''')
    parser.add_argument('--intersection', metavar=('TABLE', 'DESC'), nargs='*', default=[],
        help='''Save variants in all the tables (T1 & T2 & T3 ...) to TABLE if a name
             is specified. An optional message can be added to describe the table.''')
    parser.add_argument('--difference', metavar=('TABLE', 'DESC'), nargs='*', default=[],
        help='''Save variants in the first, but in the others (T1 - T2 - T3...) to TABLE
              if a name is specified. An optional message can be added to describe the table.''')
    parser.add_argument('-c', '--count', action='store_true',
        help='''Output number of variants for specified option (e.g. --union -c).''')


def compareTwoTables(proj, args):
    # this is the old version of vtools compare, kept for compatibility reasons
    #
    # We can use a direct query to get diff/union/intersection of tables but we cannot
    # display a progress bar during query. We therefore only use that faster method (3m38s
    # instead of 2m33s) in the case of -v0.
    direct_query = proj.verbosity is not None and proj.verbosity.startswith('0')
    cur = proj.db.cursor()
    variant_A = set()
    variant_B = set()
    if args.count or not direct_query:
        # read variants in tables[0]
        proj.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(args.tables[0], exact=False), args.tables[0]))
        cur.execute('SELECT variant_id from {};'.format(args.tables[0]))
        variant_A = set([x[0] for x in cur.fetchall()])
        # read variants in tables[1]
        proj.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(args.tables[1], exact=False), args.tables[1]))
        cur.execute('SELECT variant_id from {};'.format(args.tables[1]))
        variant_B = set([x[0] for x in cur.fetchall()])
    #
    if args.count:
        proj.logger.info('Output number of variants in A but not B, B but not A, A and B, and A or B')
        print('{}\t{}\t{}\t{}'.format(len(variant_A - variant_B), 
            len(variant_B - variant_A),
            len(variant_A & variant_B),
            len(variant_A | variant_B)
            ))
    #
    for var, opt, table, tables_A, table_B in [
            (set() if args.A_diff_B is None else variant_A - variant_B, 'EXCEPT', args.A_diff_B, args.tables[0], args.tables[1]), 
            (set() if args.B_diff_A is None else variant_B - variant_A, 'EXCEPT', args.B_diff_A, args.tables[1], args.tables[0]), 
            (set() if args.A_and_B is None else variant_A & variant_B, 'INTERSECT', args.A_and_B, args.tables[0], args.tables[1]), 
            (set() if args.A_or_B is None else variant_A | variant_B, 'UNION', args.A_or_B, args.tables[0], args.tables[1])]:
        if table is None:
            continue
        if table == 'variant':
            raise ValueError('Cannot overwrite master variant table')
        if proj.db.hasTable(table):
            new_table = proj.db.backupTable(table)
            proj.logger.warning('Existing table {} is renamed to {}.'.format(table, new_table))
        proj.createVariantTable(table)
        if direct_query:
            #proj.db.startProgress('Creating table {}'.format(table))
            cur = proj.db.cursor()
            query = 'INSERT INTO {table} SELECT variant_id FROM {table_A} {opt} SELECT variant_id FROM {table_B}'.format(opt=opt, table=table, table_A=table_A, table_B=table_B)
            proj.logger.debug(query)
            cur.execute(query)
            #proj.db.stopProgress()
        else:
            prog = ProgressBar('Writing to ' + table, len(var))
            query = 'INSERT INTO {} VALUES ({});'.format(table, proj.db.PH)
            # sort var so that variant_id will be in order, which might
            # improve database performance
            for count,id in enumerate(sorted(var)):
                cur.execute(query, (id,))
                if count % 10000 == 0:
                    prog.update(count)
            prog.done()       
        proj.db.commit()

def compareMultipleTables(proj, args):
    # We can use a direct query to get diff/union/intersection of tables but we cannot
    # display a progress bar during query. We therefore only use that faster method (3m38s
    # instead of 2m33s) in the case of -v0.
    if args.count and sum([args.difference != '', args.union != '', args.intersection != '']) > 1:
        raise ValueError('Argument --count can be used only with one operation.')
    # args.difference is
    #    None  for --difference
    #    value for --difference value
    #    ''    for not specified
    if not args.count and (args.difference is None or args.union is None or args.intersection is None):
        raise ValueError('Please specify either a table to output variants, or --count')
    #
    cur = proj.db.cursor()
    variants = []
    for table in args.tables:
        # read variants in tables[0]
        proj.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(table, exact=False), table))
        cur.execute('SELECT variant_id from {};'.format(table))
        variants.append(set([x[0] for x in cur.fetchall()]))
    #
    var_diff = set()
    var_union = set()
    var_intersect = set()
    if args.difference is not None:
        var_diff = variants[0]
        for var in variants[1:]:
           var_diff = var_diff - var 
    if args.union is not None:
        var_union = variants[0]
        for var in variants[1:]:
           var_union = var_union | var 
    if args.intersection is not None:
        var_intersect = variants[0]
        for var in variants[1:]:
           var_intersect = var_intersect & var 

    # count
    if args.count:
        if args.difference is not None:
            print(len(var_diff))
        elif args.union is not None:
            print(len(var_union))
        elif args.intersection is not None:
            print(len(var_intersect))
    #
    for var, table_with_desc in [
            (var_intersect, args.intersection),
            (var_union, args.union),
            (var_diff, args.difference)]:
        if not table_with_desc:
            continue
        table = table_with_desc[0]
        desc = table_with_desc[1] if len(table_with_desc) == 2 else ''
        if len(table_with_desc) > 2:
            raise ValueError('Only a table name and an optional table description is allowed: %s provided'.format(table_with_desc))
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
            if count % 10000 == 0:
                prog.update(count)
        prog.describeTable(table, desc, True, True)
        prog.done()       
        proj.db.commit()

def compare(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # table?
            for table in args.tables:
                if not proj.isVariantTable(table):
                    raise ValueError('Variant table {} does not exist.'.format(table))
            # this is the old behavior
            if len(args.tables) == 2 and ((args.B_diff_A is None and args.A_diff_B is None and \
                args.A_and_B is None and args.A_or_B is None and args.intersection == '' and \
                args.union == '' and args.difference == '' and args.count) or args.B_diff_A is not None
                or args.A_diff_B is not None or args.A_and_B is not None or args.A_or_B is not None):
                if args.intersection != '' or args.union != '' or args.difference != '':
                    raise ValueError('Parameters for the old and new interface cannot be mixed.')
                compareTwoTables(proj, args)
            else:
                # new interface, ignores all the A_XX_B parameters
                if args.intersection == '' and args.union == '' and args.difference == '':
                    proj.logger.warning('No action parameter is specified. Nothing to do.')
                    return
                compareMultipleTables(proj, args)
    except Exception as e:
        sys.exit(e) 

