#!/usr/bin/env python
#
# $File: variant.py $
# $LastChangedDate: 2012-06-12 20:49:24 -0500 (Tue, 12 Jun 2012) $
# $Rev: 1216 $
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
from .project import Project
from .utils import ProgressBar, runOptions, validateTableName

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
    direct_query = runOptions.verbosity is not None and runOptions.verbosity.startswith('0')
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
        validateTableName(table, exclude=['variant'])
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
    if args.count and sum([args.difference != [], args.union != [], args.intersection != []]) > 1:
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
        validateTableName(table, exclude=['variant'])
        desc = table_with_desc[1] if len(table_with_desc) == 2 else ''
        if len(table_with_desc) > 2:
            raise ValueError('Only a table name and an optional table description is allowed: %s provided'.format(table_with_desc))
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
        proj.describeTable(table, desc, True, True)
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
            if args.intersection or args.union or args.difference:
                if args.B_diff_A or args.A_diff_B or args.A_and_B or args.A_or_B:
                    raise ValueError('Cannot mix deprecated and new parameters.')
                compareMultipleTables(proj, args)
            elif args.B_diff_A or args.A_diff_B or args.A_and_B or args.A_or_B:
                compareTwoTables(proj, args)
            elif args.count:
                compareTwoTables(proj, args)
            else:
                proj.logger.warning('No action parameter is specified. Nothing to do.')
                return
    except Exception as e:
        sys.exit(e) 

