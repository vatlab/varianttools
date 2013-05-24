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
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org)
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
from .utils import ProgressBar, consolidateFieldName, typeOfValues, lineCount,\
    delayedAction, encodeTableName, decodeTableName, env
from .phenotype import Sample


def outputArguments(parser):
    parser.add_argument('table', help='''variants to output.''')
    parser.add_argument('fields', nargs='+',
        help='''A list of fields that will be outputted. SQL-compatible expressions
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
    grp.add_argument('-l', '--limit', metavar='N', type=int,
        help='''Limit output to the first N records.''')
    grp.add_argument('--build',
        help='''Output reference genome. If set to alternative build, chr and pos
            in the fields will be replaced by alt_chr and alt_pos''')
    grp.add_argument('-g', '--group_by', nargs='*', metavar='FIELD',
        help='''Group output by fields. This option is useful for aggregation output
            where summary statistics are grouped by one or more fields.''')
    grp.add_argument('--order_by', nargs='*', metavar='FIELD',
        help='''Order output by specified fields in ascending order.''')
    grp.add_argument('--dedup', default=False, action='store_true',
        help='''Remove adjacent duplicated records resulting from multiple
            entries of variants in annotation databases.''')

def outputVariants(proj, table_name, output_fields, args, query=None, reverse=False):
    '''Output selected fields'''
    # table
    table = encodeTableName(table_name)
    if not proj.isVariantTable(table):
        raise ValueError('Variant table {} does not exist.'.format(table_name))
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
    limit_clause = '' if args.limit is None or args.limit < 0 else ' LIMIT 0,{}'.format(args.limit)
    query = 'SELECT {} {} {} {} {} {};'.format(select_clause, from_clause, where_clause, group_clause, order_clause, limit_clause)
    env.logger.debug('Running query {}'.format(query))
    # if output to a file
    cur = proj.db.cursor()
    cur.execute(query)
    if args.header is not None:
        if len(args.header) == 0:
            # if no real header is given, use output_fields, but replace things like (, ), and , to space
            sys.stdout.write(args.delimiter.join([re.sub('[\W]+', ' ', x) for x in output_fields]) + '\n')
        elif args.header == ['-']:
            env.logger.info('Reading header from standard input')
            sys.stdout.write(sys.stdin.read().rstrip() + '\n')
        else:
            # other wise, use the user-provided header
            if len(args.header) != len(output_fields):
                env.logger.warning('User-provided header ({}) does not match number of fields ({})'.format(len(args.header), len(output_fields)))
            sys.stdout.write(args.delimiter.join(args.header) + '\n')
    # output with a warning to potentially duplicated lines
    last_line = None
    for count, rec in enumerate(cur):
        line = args.delimiter.join([args.na if x is None else str(x) for x in rec]) + '\n'
        if args.dedup:
            if line == last_line:
                continue
            else:
                last_line = line
        sys.stdout.write(line)


def output(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            outputVariants(proj, args.table, args.fields, args)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1) 

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
        help='''A list of fields that will be outputted. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')


def select(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # separate table and message
            if args.to_table:
                if len(args.to_table) > 2:
                    raise ValueError('Only a table name and an optional message is allowed for parameter to_table')
                if args.to_table[0] == 'variant':
                    raise ValueError('Cannot overwrite the master variant table.')
                if '*' in args.to_table[0] or '?' in args.to_table[0]:
                    env.logger.warning('Use of wildcard character * or ? in table names is not recommended because such names can be expanded to include other tables in some commands.')
                args.table_desc = args.to_table[1] if len(args.to_table) == 2 else ''
                args.to_table = args.to_table[0]
            # table?
            if not proj.isVariantTable(encodeTableName(args.from_table)):
                raise ValueError('Variant table {} does not exist.'.format(args.from_table))
            if not args.to_table and not args.output and not args.count:
                env.logger.warning('Neither --to_table and --output/--count is specified. Nothing to do.')
                return
            if len(args.condition) > 0:    
                # fields? We need to replace things like sift_score to dbNSFP.sift_score
                condition, fields = consolidateFieldName(proj, encodeTableName(args.from_table), ' AND '.join(['({})'.format(x) for x in args.condition]))
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
                    #s = delayedAction(env.logger.info, 'Indexing {}'.format(field))
                    #cur = proj.db.cursor()
                    #try:
                    #    query = 'CREATE INDEX IF NOT EXISTS {0}.{0}_{1} ON {0} ({1} ASC);'.format(db, fld)
                    #    env.logger.debug(query)
                    #    cur.execute(query)
                    #except Exception as e:
                    #    env.logger.debug('Failed to create index: {}'.format(e))
                    #del s
                # 
                fields_info = sum([proj.linkFieldToTable(x, encodeTableName(args.from_table)) for x in fields], [])
                # WHERE clause: () is important because OR in condition might go beyond condition
                where_clause = ' WHERE ({})'.format(condition)
                # 
                # FROM clause
                from_clause = 'FROM {} '.format(encodeTableName(args.from_table))
                # avoid duplicate
                #
                processed = set()
                for table, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                    if (table.lower(), conn.lower()) not in processed:
                        from_clause += ' LEFT OUTER JOIN {} ON {}'.format(table, conn)
                        processed.add((table.lower(), conn.lower()))
            else:
                # select all variants
                where_clause = ' WHERE 1 '
                from_clause = 'FROM {} '.format(encodeTableName(args.from_table))
            # if limiting to specified samples
            if args.samples:
                # we save genotype in a separate database to keep the main project size tolerable.
                proj.db.attach(proj.name + '_genotype')
                IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.samples]))
                if len(IDs) == 0:
                    env.logger.warning('No sample is selected by condition: {}'.format(' AND '.join(['({})'.format(x) for x in args.samples])))
                    # nothing will be selected
                    where_clause += ' AND 0'
                #
                # This special case does not hold because sometimes variants are imported without sample information.
                #
                #elif len(IDs) == proj.db.numOfRows('sample'):
                #    env.logger.info('All {} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                #    # we do not have to add anything to where_clause
                elif len(IDs) < 50:  
                    # we allow 14 tables in other 'union' or from condition...
                    env.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                    where_clause += ' AND ({}.variant_id IN ({}))'.format(
                        encodeTableName(args.from_table), 
                        '\nUNION '.join(['SELECT variant_id FROM {}_genotype.genotype_{}'.format(proj.name, id) for id in IDs])) 
                else:
                    # we have to create a temporary table and select variants sample by sample
                    # this could be done in parallel if there are a large number of samples, but that needs a lot more
                    # code, and perhaps RAM
                    env.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
                    cur = proj.db.cursor()
                    BLOCK_SIZE = 64
                    NUM_BLOCKS = len(IDs) // BLOCK_SIZE + 1
                    myIDs = list(IDs)
                    myIDs.sort()
                    merged_table = '__variants_from_samples'
                    query = 'CREATE TEMPORARY TABLE {} (variant_id INT);'.format(merged_table)
                    # env.logger.debug(query)
                    cur.execute(query)
                    prog = ProgressBar('Collecting sample variants', len(IDs)) if NUM_BLOCKS > 1 else None
                    count = 0
                    for i in range(NUM_BLOCKS):
                        # step 1: create a table that holds all
                        block_IDs = myIDs[(i*BLOCK_SIZE):((i+1)*BLOCK_SIZE)]
                        if len(block_IDs) == 0:
                            continue
                        query = 'INSERT INTO {} {};'.format(merged_table,
                            '\nUNION '.join(['SELECT variant_id FROM {}_genotype.genotype_{}'.format(proj.name, id) for id in block_IDs]))
                        #env.logger.debug(query)
                        cur.execute(query)
                        count += len(block_IDs)
                        if prog:
                            prog.update(count)
                    if prog:
                        prog.done()
                    where_clause += ' AND ({}.variant_id IN (SELECT variant_id FROM {}))'.format(
                        encodeTableName(args.from_table), merged_table)
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
                query = 'SELECT COUNT(DISTINCT {}.variant_id) {} {};'.format(encodeTableName(args.from_table),
                    from_clause, where_clause)
                env.logger.debug('Running query {}'.format(query))
                proj.db.startProgress('Counting variants')
                cur = proj.db.cursor()
                cur.execute(query)
                count = cur.fetchone()[0]
                # exclude ...
                if reverse:
                    count = proj.db.numOfRows(encodeTableName(args.from_table)) - int(count)
                proj.db.stopProgress()
                #
                print(count)
            # case 2: to table
            elif args.to_table:
                if proj.db.hasTable(encodeTableName(args.to_table)):
                    new_table = proj.db.backupTable(encodeTableName(args.to_table))
                    env.logger.warning('Existing table {} is renamed to {}.'.format(args.to_table, decodeTableName(new_table)))
                #
                proj.createVariantTable(encodeTableName(args.to_table))
                proj.describeTable(encodeTableName(args.to_table), args.table_desc, True, True)
                if not reverse:
                    query = 'INSERT INTO {0} SELECT DISTINCT {1}.variant_id {2} {3};'.format(encodeTableName(args.to_table), 
                        encodeTableName(args.from_table),
                        from_clause, where_clause)
                else:
                    query = 'INSERT INTO {0} SELECT DISTINCT {1}.variant_id FROM {1} WHERE {1}.variant_id NOT IN (SELECT {1}.variant_id {2} {3});'.\
                        format(encodeTableName(args.to_table), encodeTableName(args.from_table), from_clause, where_clause)
                env.logger.debug('Running query {}'.format(query))
                #
                cur = proj.db.cursor()
                proj.db.startProgress('Running')
                cur.execute(query)
                proj.db.stopProgress()
                proj.db.commit()
                #
                count = proj.db.numOfRows(encodeTableName(args.to_table))
                env.logger.info('{} variants selected.'.format(count))
                if args.output:
                    outputVariants(proj, args.to_table, args.output, args)
                if args.count:
                    print(count)
            # case 3: output, but do not write to table, and not count
            elif args.output: 
                query = 'SELECT DISTINCT {}.variant_id {} {}'.format(encodeTableName(args.from_table),
                    from_clause, where_clause)
                outputVariants(proj, args.from_table, args.output, args, query, reverse)
            # 
            # clean up temporary table
            if args.samples and proj.db.hasTable('__variants_from_samples'): 
                cur.execute('DROP TABLE __variants_from_samples')
    except Exception as e:
        env.logger.error(e)
        sys.exit(1) 

def excludeArguments(parser):
    parser.add_argument('from_table', help='''Source variant table.''')
    parser.add_argument('condition', nargs='*', default=[],
        help='''Conditions by which variants are excluded. Multiple arguments are
            automatically joined by 'AND' so 'OR' conditions should be provided by
            a single argument with conditions joined by 'OR'. If unspecified, all
            variants (except those excluded by parameter --samples) will be excluded.''')
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-t', '--to_table', nargs='*', metavar=('TABLE', 'DESC'),
        help='''Destination variant table.''')
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument('-c', '--count', action='store_true',
        help='''Output number of variant, which is a shortcut to '--output count(1)'.''')
    grp.add_argument('-o', '--output', nargs='*', metavar='FIELDS', default=[],
        help='''A list of fields that will be outputted. SQL-compatible expressions
            or functions such as "pos-1", "count(1)" or "sum(num)" are also allowed. ''')

def exclude(args):
    select(args, reverse=True)

