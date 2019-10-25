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

import copy
import re
import sys
from argparse import SUPPRESS
from collections import defaultdict

from .geno_store import GenoStore
from .project import Project
from .utils import (ProgressBar, decodeTableName, encodeTableName, env,
                    matchName)


def compareArguments(parser):
    parser.add_argument(
        'tables',
        nargs='*',
        help='''variant tables to compare. Wildcard characters * and ? can be
            used to specify multiple tables. A table name will be automatically
            repeated for the comparison of genotype of multiple samples if only
            one table is specified.''')
    parser.add_argument(
        '--union',
        metavar=('TABLE', 'DESC'),
        nargs='*',
        help='''Print the number (default) or save variants with TYPE in the
            TYPE of any of the tables (T1 | T2 | T3 ...) to TABLE if a name is
            specified. An optional message can be added to describe the table. '''
    )
    parser.add_argument(
        '--intersection',
        metavar=('TABLE', 'DESC'),
        nargs='*',
        help='''Print the number (default) or save variants with TYPE in the TYPE
            of all the tables (T1 & T2 & T3 ...) to TABLE if a name is specified.
            An optional message can be added to describe the table.''')
    parser.add_argument(
        '--difference',
        metavar=('TABLE', 'DESC'),
        nargs='*',
        help='''Print the number (default) or save variants with TYPE in the TYPE
            of the first, but not in the TYPE of others (T1 - T2 - T3...) to TABLE
            if a name is specified. An optional message can be added to describe
            the table.''')
    parser.add_argument(
        '--expression',
        metavar=('EXPR', 'DESC'),
        nargs='+',
        help='''Evaluate a set expression with table names representing variants
            in these tables. Operators | (or), & (and), - (difference) and ^
            (A or B but not both) are allowed. The results will be saved to table
            if the result is assigned to a name (e.g. --expression 'D=A-(B&C)').
            The table names in the expression can be written as _1, _2 etc if
            tables are listed before the option, and be used to populate the
            list of tables if it was left unspecified.''')
    parser.add_argument(
        '--mode',
        choices=['variant', 'site', 'genotype'],
        help='''Compare variants (chr, pos, ref, alt), site (chr, pos), or
            genotype (chr, pos, ref, alt, GT for a specified sample) of
            variants. The results are variants from all input tables that match
            specified condition. The default comparison TYPE compares variants
            in input variant tables. For the comparison of sites, the position
            of all variants are collected and compared, and variants (in all tables)
            with site in resulting set of sites are returned. For the comparison
            of genotypes, the genotypes of specified samples for all variants
            (see option --samples) are collected and compared, and variants
            with genotype in resulting set of genotypes are returned. The results
            of genotype comparisons are affected by runtime option
            treat_missing_as_wildtype because items with missing genotype (chr,
            pos, ref, alt, NULL) are excluded if treat_missing_as_wildtype is
            false (default), and are treated as (chr, pos, ref, alt, 0) otherwise.
            The default comparison type is variant, or genotype if option
            --samples is specified.''')
    parser.add_argument(
        '--samples',
        nargs='*',
        help='''A list of sample names corresponding to the variant tables to
            compare. An error will be raised if a sample name matches no or multiple
            samples or if a sample does not have any genotype.''')
    parser.add_argument('-c', '--count', action='store_true', help=SUPPRESS)
    parser.add_argument(
        '--A_diff_B', '--A-diff-B', nargs='+', metavar='TABLE', help=SUPPRESS)
    parser.add_argument(
        '--B_diff_A', '--B-diff-A', nargs='+', metavar='TABLE', help=SUPPRESS)
    parser.add_argument(
        '--A_and_B', '--A-and-B', nargs='+', metavar='TABLE', help=SUPPRESS)
    parser.add_argument(
        '--A_or_B', '--A-or-B', nargs='+', metavar='TABLE', help=SUPPRESS)


def countVariantDifference(proj, args):
    cur = proj.db.cursor()
    variants = []
    if len(args.tables) > 2:
        env.logger.warning(
            'Only the first two specified tables will be compared for option --count.'
        )
    for table in args.tables[:2]:
        # read variants in tables[0]
        env.logger.info('Reading approximately {:,} variants in {}...'.format(
            proj.db.numOfRows(encodeTableName(table), exact=False), table))
        cur.execute('SELECT variant_id from {};'.format(encodeTableName(table)))
        variants.append(set([x[0] for x in cur.fetchall()]))
    #
    env.logger.info(
        'Number of variants in A but not B, B but not A, A and B, and A or B')
    print(('{}\t{}\t{}\t{}'.format(
        len(variants[0] - variants[1]), len(variants[1] - variants[0]),
        len(variants[0] & variants[1]), len(variants[0] | variants[1]))))


def countSiteDifference(proj, args):
    cur = proj.db.cursor()
    if len(args.tables) > 2:
        env.logger.warning(
            'Only the first two specified tables will be compared for option --count.'
        )
    #
    all_variants = set()
    sites = []
    for table in args.tables:
        s = defaultdict(set)
        env.logger.info(
            'Reading locations of approximately {:,} variants in {}...'.format(
                proj.db.numOfRows(encodeTableName(table), exact=False), table))
        if table == 'variant':
            cur.execute('SELECT variant_id, chr, pos FROM {};'.format(
                encodeTableName(table)))
        else:
            cur.execute(
                'SELECT {0}.variant_id, variant.chr, variant.pos '
                'FROM {0} LEFT OUTER JOIN variant ON {0}.variant_id = variant.variant_id'
                .format(encodeTableName(table)))
        for id, chr, pos in cur:
            s[(chr, pos)].add(id)
            all_variants.add(id)
        env.logger.debug('There are {} unique sites in table {}'.format(
            len(s), table))
        sites.append(s)
    #
    env.logger.info(
        'Number of sites in A only, B only, in A and B, and in A or B')
    site0 = set(sites[0].keys())
    site1 = set(sites[1].keys())
    site_counts = '{}\t{}\t{}\t{}'.format(
        len(site0 - site1), len(site1 - site0), len(site0 & site1),
        len(site0 | site1))
    env.logger.info(site_counts)
    env.logger.info('Number of variants in both tables with locations in '
                    'A only, B only, in A and B, and in A or B')
    print((
        '{}\t{}\t{}\t{}'.format(
            # variants in A, with location not in B
            sum([len(sites[0][x]) for x in site0 - site1]),
            # variants in B, with location not in A
            sum([len(sites[1][x]) for x in site1 - site0]),
            # variants with location in both A & B
            sum([len(sites[0][x] | sites[1][x]) for x in site0 & site1]),
            # variants with location in either A or B
            len(all_variants))))


def countGenotypeDifference(proj, args):
    cur = proj.db.cursor()
    geno = []
    NULL_to_0 = env.treat_missing_as_wildtype
    if len(args.tables) > 2:
        env.logger.warning(
            'Only the first two specified tables will be compared for option --count.'
        )
    # read variants in tables[0]
    for table, sample_name, sample_ID in zip(args.tables, args.samples,
                                             args.sample_IDs):
        # g = set()
        env.logger.info(
            'Reading genotype of sample {} of approximately {:,} variants in {}...'
            .format(sample_name,
                    proj.db.numOfRows(encodeTableName(table), exact=False),
                    table))
        # cur.execute('SELECT {0}.variant_id, geno.GT FROM {0} LEFT OUTER JOIN '
        #     '{1}_genotype.genotype_{2} geno ON {0}.variant_id = geno.variant_id'
        #         .format(encodeTableName(table), proj.name, sample_ID))
        # for id, GT in cur:
        #     if GT is None:
        #         if NULL_to_0:
        #             g.add((id, 0))
        #     else:
        #         g.add((id, GT))
        store = GenoStore(proj)
        g = store.get_Genotype(cur, encodeTableName(table), proj, sample_ID)
        geno.append(g)
        env.logger.debug(
            'There are {} genotypes with variants in table {} in sample {}'
            .format(len(g), table, sample_name))
    #
    env.logger.info(
        'Number of genotypes in A only, B only, in A and B, and in A or B')
    geno_counts = '{}\t{}\t{}\t{}'.format(
        len(geno[0] - geno[1]), len(geno[1] - geno[0]), len(geno[0] & geno[1]),
        len(geno[0] | geno[1]))
    env.logger.info(geno_counts)
    env.logger.info(
        'Number of variants with genotypes in A only, B only, in A and B, and in A or B'
    )
    print(('{}\t{}\t{}\t{}'.format(
        len(set([x[0] for x in (geno[0] - geno[1])])),
        len(set([x[0] for x in (geno[1] - geno[0])])),
        len(set([x[0] for x in (geno[0] & geno[1])])),
        len(set([x[0] for x in (geno[0] | geno[1])])))))


def compareTables(proj, args):
    # We can use a direct query to get diff/union/intersection of tables but we cannot
    # display a progress bar during query. We therefore only use that faster method (3m38s
    # instead of 2m33s) in the case of -v0.
    # args.difference is
    #    None  for --difference
    #    value for --difference value
    #    ''    for not specified
    var_diff = set()
    var_union = set()
    var_intersect = set()
    var_expression = set()
    cur = proj.db.cursor()
    if args.mode in [None, 'variant']:
        for idx, table in enumerate(args.tables):
            # read variants in tables[0]
            env.logger.info('Reading {:,} variants in {}...'.format(
                proj.db.numOfRows(encodeTableName(table), exact=False), table))
            cur.execute('SELECT variant_id from {};'.format(
                encodeTableName(table)))
            var = set([x[0] for x in cur.fetchall()])
            if idx == 0:
                if args.difference is not None:
                    var_diff = copy.copy(var)
                if args.union is not None:
                    var_union = copy.copy(var)
                if args.intersection is not None:
                    var_intersect = copy.copy(var)
                if args.expression is not None:
                    variants = [var]
            else:
                if args.difference is not None:
                    var_diff -= var
                if args.union is not None:
                    var_union |= var
                if args.intersection is not None:
                    var_intersect &= var
                if args.expression is not None:
                    variants.append(var)
        # only the expressionc ase needs to retain all variables
        if args.expression is not None:
            var_dict = {}
            for idx, var in enumerate(variants):
                var_dict['_{}'.format(idx + 1)] = var
            # evaluate the expression
            try:
                var_expression = eval(args.eval_expr, var_dict, var_dict)
            except Exception as e:
                raise RuntimeError(
                    'Failed to evaluate set expression {}: {}'.format(
                        args.eval_expr, e))
    elif args.mode == 'site':
        for idx, table in enumerate(args.tables):
            # read sites in tables[0]
            env.logger.info(
                'Reading locations of approximately {:,} sites in {}...'.format(
                    proj.db.numOfRows(encodeTableName(table), exact=False),
                    table))
            if table == 'variant':
                cur.execute('SELECT variant_id, chr, pos FROM {};'.format(
                    encodeTableName(table)))
            else:
                cur.execute(
                    'SELECT {0}.variant_id, variant.chr, variant.pos '
                    'FROM {0} LEFT OUTER JOIN variant ON {0}.variant_id = variant.variant_id'
                    .format(encodeTableName(table)))
            #
            # we can save some more memory if we use a dictionary of chromosome
            # (e.g. separate chromosomes), but let us stop here for now.
            v = defaultdict(set)
            for id, chr, pos in cur:
                v[(chr, pos)].add(id)
            env.logger.info('Unique sites in table {}: {}'.format(
                table, len(v)))
            if idx == 0:
                if args.difference is not None:
                    site_diff = copy.copy(v)
                if args.union is not None:
                    site_union = copy.copy(v)
                if args.intersection is not None:
                    site_intersect = copy.copy(v)
                if args.expression is not None:
                    sites = [v]
            else:
                if args.difference is not None:
                    # in diff but not in v
                    site_diff = {
                        x: y for x, y in list(site_diff.items()) if x not in v
                    }
                if args.union is not None:
                    # in any of them
                    site_union = {
                        x: (y if x not in v else (y | v[x]))
                        for x, y in list(site_union.items())
                    }
                    site_union.update({
                        x: y for x, y in list(v.items()) if x not in site_union
                    })
                if args.intersection is not None:
                    # in both, need to combine ids
                    site_intersect = {
                        x: (y | v[x])
                        for x, y in list(site_intersect.items())
                        if x in v
                    }
                if args.expression is not None:
                    sites.append(v)
        #
        if args.difference is not None:
            env.logger.info('Unique sites in table {} only: {}'.format(
                args.tables[0], len(site_diff)))
            for k in list(site_diff.values()):
                var_diff |= k
        if args.union is not None:
            env.logger.info('Sites in any of the tables: {}'.format(
                len(site_union)))
            for k in list(site_union.values()):
                var_union |= k
        if args.intersection is not None:
            env.logger.info('Sites in all tables: {}'.format(
                len(site_intersect)))
            for k in list(site_intersect.values()):
                var_intersect |= k
        if args.expression is not None:
            site_dict = {}
            for idx, s in enumerate(sites):
                site_dict['_{}'.format(idx + 1)] = list(s.keys())
            # evaluate the expression
            try:
                site_expr = eval(args.eval_expr, site_dict, site_dict)
            except Exception as e:
                raise RuntimeError(
                    'Failed to evaluate set expression {}: {}'.format(
                        args.eval_expr, e))
            env.logger.info('Sites in resulting tables: {}'.format(
                len(site_expr)))
            for k in site_expr:
                for var in sites:
                    if k in var:
                        var_expression |= var[k]
    else:
        # genotype
        NULL_to_0 = env.treat_missing_as_wildtype
        for idx, (table, sample, sample_ID) in enumerate(
                zip(args.tables, args.samples, args.sample_IDs)):
            # read geno in tables[0]
            env.logger.info(
                'Reading genotypes of sample {} of approximately {:,} geno in {}...'
                .format(sample,
                        proj.db.numOfRows(encodeTableName(table), exact=False),
                        table))
            # cur.execute('SELECT {0}.variant_id, geno.GT FROM {0} LEFT OUTER JOIN '
            #     '{1}_genotype.genotype_{2} geno ON {0}.variant_id = geno.variant_id'
            #         .format(encodeTableName(table), proj.name, id))
            # #
            # v = set()
            # for id, GT in cur:
            #     if GT is None:
            #         if NULL_to_0:
            #             v.add((id, 0))
            #     else:
            #         v.add((id, GT))

            store = GenoStore(proj)
            v = store.get_Genotype(cur, encodeTableName(table), proj, sample_ID)
            #
            if idx == 0:
                if args.difference is not None:
                    geno_diff = copy.copy(v)
                if args.union is not None:
                    geno_union = copy.copy(v)
                if args.intersection is not None:
                    geno_intersect = copy.copy(v)
                if args.expression is not None:
                    geno = [v]
            else:
                if args.difference is not None:
                    geno_diff -= v
                if args.union is not None:
                    geno_union |= v
                if args.intersection is not None:
                    geno_intersect &= v
                if args.expression is not None:
                    geno.append(v)
        #
        if args.difference is not None:
            env.logger.info('Genotypes in sample {} only: {} '.format(
                args.samples[0], len(geno_diff)))
            var_diff = set([x[0] for x in geno_diff])
        if args.union is not None:
            env.logger.info('Genotypes in any of the samples: {}'.format(
                len(geno_union)))
            var_union = set([x[0] for x in geno_union])
        if args.intersection is not None:
            env.logger.info('Genotypes in all samples: {}'.format(
                len(geno_intersect)))
            var_intersect = set([x[0] for x in geno_intersect])
        if args.expression is not None:
            geno_dict = {}
            for idx, g in enumerate(geno):
                geno_dict['_{}'.format(idx + 1)] = g
            # evaluate the expression
            try:
                var_expression = set(
                    [x[0] for x in eval(args.eval_expr, geno_dict, geno_dict)])
            except Exception as e:
                raise RuntimeError(
                    'Failed to evaluate set expression {}: {}'.format(
                        args.eval_expr, e))
    #
    for var, table_with_desc in [(var_intersect, args.intersection),
                                 (var_union, args.union),
                                 (var_diff, args.difference),
                                 (var_expression, args.expression)]:
        if not table_with_desc:
            continue
        table = table_with_desc[0]
        if table == 'variant':
            raise ValueError('Cannot overwrite the master variant table')
        if '*' in table or '?' in table:
            env.logger.warning(
                'Use of wildcard character * or ? in table '
                'names is not recommended because such names can be expanded '
                'to include other tables in some commands.')
        desc = table_with_desc[1] if len(table_with_desc) == 2 else ''
        if len(table_with_desc) > 2:
            raise ValueError(
                'Only a table name and an optional table '
                'description is allowed: %s provided'.format(table_with_desc))
        proj.createVariantTable(encodeTableName(table))
        prog = ProgressBar('Writing to ' + table, len(var))
        query = 'INSERT INTO {} VALUES ({},{});'.format(
            encodeTableName(table), proj.db.PH, proj.db.PH)

        # sort var so that variant_id will be in order, which might
        # improve database performance
        for count, id in enumerate(sorted(var)):
            cur.execute(query, (
                id,
                None,
            ))
            if count % 10000 == 0:
                prog.update(count)
        proj.describeTable(encodeTableName(table), desc, True, True)
        prog.done()
        proj.db.commit()
    # count by default
    if args.difference is not None:
        print((len(var_diff)))
    if args.union is not None:
        print((len(var_union)))
    if args.intersection is not None:
        print((len(var_intersect)))
    if args.expression is not None:
        print((len(var_expression)))


def parseExpression(expr, passed_tables):
    pieces = re.split(r'''("[^"]+"|'[^']+'|[a-zA-Z0-9_]+|\(|\)|\^|\||\&|-|=)''',
                      expr)
    result_table = None
    tables = [x for x in passed_tables]
    new_expr = []
    for idx, i in enumerate(pieces):
        if i.startswith('_'):
            try:
                i_value = int(i[1:]) - 1
            except:
                raise ValueError('Unacceptable table name {}'.format(i))
            if i_value >= len(tables):
                raise ValueError(
                    'Unacceptable table name {}: index out of range'.format(i))
            new_expr.append(i)
        elif (i.startswith('"') and i.endswith('"')) or (i.startswith("'") and
                                                         i.endswith("'")):
            if '=' in pieces and result_table is None:
                result_table = i[1:-1]
            elif i[1:-1] in tables:
                new_expr.append('_{}'.format(tables.index(i[1:-1]) + 1))
            else:
                tables.append(i[1:-1])
                new_expr.append('_{}'.format(len(tables)))
        elif i.replace('_', '').isalnum():
            if '=' in pieces and result_table is None:
                result_table = i
            elif i in tables:
                new_expr.append('_{}'.format(tables.index(i) + 1))
            else:
                tables.append(i)
                new_expr.append('_{}'.format(len(tables)))
        elif i.strip() == '=':
            pass
        elif i.strip() in ['(', ')', '&', '|', '-', '^']:
            new_expr.append(i)
        elif i.strip():
            raise ValueError('Invalid set expression for item {}: {})'.format(
                i, expr))
    return result_table, tables, ' '.join(new_expr).strip()


def compare(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # expand wildcard characters in args.tables
            tables = []
            allTables = proj.getVariantTables()
            for table in args.tables:
                if '?' in table or '*' in table:
                    match = False
                    for tbl in [decodeTableName(x) for x in allTables]:
                        if matchName(table, tbl):
                            tables.append(tbl)
                            match = True
                    if not match:
                        # * should match a table with '*' in its name.
                        env.logger.warning(
                            'Name {} does not match any existing table.'.format(
                                table))
                else:
                    tables.append(table)
            # set args.tables to its expanded version
            args.tables = tables
            #
            # if an expression is given but no table is specified
            # they will be extracted
            if args.expression:
                result_table, tables, new_expr = parseExpression(
                    args.expression[0], args.tables)
                if not args.tables:
                    args.tables = tables
                    env.logger.info('Comparing tables {}'.format(', '.join(
                        args.tables)))
                else:
                    if any([x not in args.tables for x in tables]):
                        raise ValueError(
                            'All tables in --expression must be '
                            'listed if not using the default table list.')
                if result_table is None:
                    if len(args.expression) > 1:
                        env.logger.warning(
                            'No table is created for expression {}, but a description is given'
                            .format(args.expression[0]))
                    args.expression = []
                else:
                    args.expression = [
                        result_table,
                        args.expression[1] if len(args.expression) > 1 else ''
                    ]
                args.eval_expr = new_expr
            elif not args.tables:
                raise ValueError('Please specify tables to compare.')
            # table?
            for table in args.tables:
                if not proj.isVariantTable(encodeTableName(table)):
                    raise ValueError(
                        'Variant table {} does not exist.'.format(table))
            #
            if not args.tables:
                raise ValueError('No table to compare')
            #
            # type of comparison
            if args.samples:
                if args.mode is not None and args.mode != 'genotype':
                    raise ValueError(
                        'Cannot specify samples for non-genotype comparisons.')
                args.mode = 'genotype'
                if len(args.samples) == 1:
                    raise ValueError(
                        'Please specify more than one sample to be compared.')
                if len(args.tables) == 1:
                    # automatically expand tables to multiple tables
                    args.tables = [args.tables[0]] * len(args.samples)
                elif len(args.tables) != len(args.samples):
                    raise ValueError(
                        'Please specify a variant table for all or each of the samples.'
                    )
                # check sample names
                proj.db.attach(proj.name + '_genotype')
                args.sample_IDs = []
                for name in args.samples:
                    IDs = proj.selectSampleByPhenotype(
                        "sample_name = '{}'".format(name))
                    if len(IDs) == 0:
                        raise ValueError(
                            "No sample with name '{}' is identified.".format(
                                name))
                    elif len(IDs) > 1:
                        raise ValueError(
                            "More than one sample with name '{}' is identified."
                            .format(name))
                    if proj.store == "sqlite":
                        if not 'GT' in proj.db.getHeaders(
                                '{}_genotype.genotype_{}'.format(
                                    proj.name, IDs[0])):
                            raise ValueError(
                                'Sample {} does not have any genotype.'.format(
                                    name))
                    args.sample_IDs.append(IDs[0])
            elif args.mode == 'genotype':
                raise ValueError(
                    'Please specify samples to be compared for genotype comparison.'
                )
            if len(args.tables) == 1:
                raise ValueError(
                    'Please specify at least two tables to compare.')
            #
            if args.B_diff_A or args.A_diff_B or args.A_and_B or args.A_or_B:
                raise ValueError(
                    'Options B_diff_A, A_diff_B, A_and_B and A_or_B are deprecated.'
                )
            if args.intersection is not None or args.union is not None or \
                args.difference is not None or args.expression is not None:
                compareTables(proj, args)
            else:
                if args.mode in [None, 'variant']:
                    countVariantDifference(proj, args)
                elif args.mode == 'site':
                    countSiteDifference(proj, args)
                elif args.mode == 'genotype':
                    countGenotypeDifference(proj, args)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
