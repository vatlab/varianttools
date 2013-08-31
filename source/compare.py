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
from argparse import SUPPRESS
from .project import Project
from .utils import ProgressBar, env, encodeTableName, decodeTableName
from collections import defaultdict

def compareArguments(parser):
    parser.add_argument('tables', nargs='+',
        help='''variant tables to compare. Wildcard characters * and ? can be 
            used to specify multiple tables. A table name will be automatically
            repeated for the comparison of genotype of multiple samples if only
            one table is specified.''')
    parser.add_argument('--union', metavar=('TABLE', 'DESC'), nargs='*', 
        help='''Save variants with TYPE in the TYPE of any of the tables (T1 |
            T2 | T3 ...) to TABLE if a name is specified. An optional message
            can be added to describe the table. ''')
    parser.add_argument('--intersection', metavar=('TABLE', 'DESC'), nargs='*', 
        help='''Save variants with TYPE in the TYPE of all the tables (T1 & T2 
            & T3 ...) to TABLE if a name is specified. An optional message can 
            be added to describe the table.''')
    parser.add_argument('--difference', metavar=('TABLE', 'DESC'), nargs='*',
        help='''Save variants with TYPE in the TYPE of the first, but not in the
            TYPE of others (T1 - T2 - T3...) to TABLE if a name is specified.
            An optional message can be added to describe the table.''')
    parser.add_argument('-c', '--count', action='store_true',
        help='''Output the number of variants for specified option (e.g. --union
            -c) without writing the variants to a table.''')
    parser.add_argument('--type', choices=['variant', 'site', 'genotype'],
        help='''Compare variants (chr, pos, ref, alt), site (chr, pos), or
            or genotype (chr, pos, ref, alt, GT for a specified sample) of 
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
            Because of the way variants are generated, intersection operations 
            of site comparison might produce more variants than the first table
            (because distinct variants from other tables might returned if they
            share locations with variants in the first table), and union
            operation of genotype comparison might produce fewer variants (because
            of missing genotypes). The default comparison type is variant, or 
            genotype if option --samples is specified.''')
    parser.add_argument('--samples', nargs='*',
        help='''A list of sample names corresponding to the variant tables to
            compare. An error will be raised if a sample name matches no or multiple
            samples or if a sample does not have any genotype.''')
    parser.add_argument('--A_diff_B', nargs='+', metavar= 'TABLE', help=SUPPRESS)
    parser.add_argument('--B_diff_A', nargs='+', metavar= 'TABLE', help=SUPPRESS)
    parser.add_argument('--A_and_B', nargs='+', metavar= 'TABLE', help=SUPPRESS)
    parser.add_argument('--A_or_B', nargs='+', metavar= 'TABLE', help=SUPPRESS)


def countVariantDifference(proj, args):
    cur = proj.db.cursor()
    variants = []
    if len(args.tables) > 2:
        env.logger.warning('Only the first two specified tables will be compared for option --count.')
    for table in args.tables[:2]:
        # read variants in tables[0]
        env.logger.info('Reading approximately {:,} variants in {}...'
            .format(proj.db.numOfRows(encodeTableName(table), exact=False), table))
        cur.execute('SELECT variant_id from {};'.format(encodeTableName(table)))
        variants.append(set([x[0] for x in cur.fetchall()]))
    #
    env.logger.info('Number of variants in A but not B, B but not A, A and B, and A or B')
    print('{}\t{}\t{}\t{}'.format(
        len(variants[0] - variants[1]), 
        len(variants[1] - variants[0]),
        len(variants[0] & variants[1]),
        len(variants[0] | variants[1])
        ))

def countSiteDifference(proj, args):
    cur = proj.db.cursor()
    if len(args.tables) > 2:
        env.logger.warning('Only the first two specified tables will be compared for option --count.')
    #
    all_variants = set()
    sites = []
    for table in args.tables:
        s = defaultdict(set)
        env.logger.info('Reading locations of approximately {:,} variants in {}...'
            .format(proj.db.numOfRows(encodeTableName(table), exact=False), table))
        if table == 'variant':
            cur.execute('SELECT variant_id, chr, pos FROM {};'.format(encodeTableName(table)))
        else:
            cur.execute('SELECT {0}.variant_id, variant.chr, variant.pos '
                'FROM {0} LEFT OUTER JOIN variant ON {0}.variant_id = variant.variant_id'
                .format(encodeTableName(table)))
        for id, chr, pos in cur:
            s[(chr, pos)].add(id)
            all_variants.add(id)
        env.logger.debug('There are {} unique sites in table {}'
            .format(len(s), table))
        sites.append(s)
    #
    env.logger.info('Number of variants in both tables with locations in '
        'A only, B only, in A and B, and in A or B')
    print('{}\t{}\t{}\t{}'.format(
        # variants in A, with location not in B
        sum([len(y) for x,y in sites[0].items() if x not in sites[1]]),
        # variants in B, with location not in A
        sum([len(y) for x,y in sites[1].items() if x not in sites[0]]),
        # variants with location in both A & B
        sum([len(y | sites[1][x]) for x,y in sites[0].items() if x in sites[1]]),
        # variants with location in either A or B
        len(all_variants)
    ))


def countGenotypeDifference(proj, args):
    cur = proj.db.cursor()
    geno = []
    NULL_to_0 = env.treat_missing_as_wildtype
    if len(args.tables) > 2:
        env.logger.warning('Only the first two specified tables will be compared for option --count.')
    # read variants in tables[0]
    for table, sample_name, sample_ID in zip(args.tables, args.samples, args.sample_IDs):
        g = set()
        env.logger.info('Reading genotype of sample {} of approximately {:,} variants in {}...'
            .format(sample_name, proj.db.numOfRows(encodeTableName(table),
                exact=False), table))
        cur.execute('SELECT {0}.variant_id, geno.GT FROM {0} LEFT OUTER JOIN '
            '{1}_genotype.genotype_{2} geno ON {0}.variant_id = geno.variant_id'
                .format(encodeTableName(table), proj.name, sample_ID))
        for id, GT in cur:
            if GT is None:
                if NULL_to_0:
                    g.add((id, 0))
            else:
                g.add((id, GT))
        geno.append(g)
        env.logger.debug('There are {} genotypes with variants in table {} in sample {}'
            .format(len(g), table, sample_name))
    #
    env.logger.info('Number of genotypes in A only, B only, in A and B, and in A or B')
    print('{}\t{}\t{}\t{}'.format(
        len(set([x[0] for x in (geno[0] - geno[1])])),
        len(set([x[0] for x in (geno[1] - geno[0])])),
        len(set([x[0] for x in (geno[0] & geno[1])])),
        len(set([x[0] for x in (geno[0] | geno[1])]))
    ))


def compareTables(proj, args):
    # We can use a direct query to get diff/union/intersection of tables but we cannot
    # display a progress bar during query. We therefore only use that faster method (3m38s
    # instead of 2m33s) in the case of -v0.
    if args.count and sum([args.difference is not None, args.union is not None, args.intersection is not None]) > 1:
        raise ValueError('Argument --count can be used only with one operation.')
    # args.difference is
    #    None  for --difference
    #    value for --difference value
    #    ''    for not specified
    if not args.count and (args.difference == [] or args.union == [] or args.intersection == []):
        raise ValueError('Please specify either a table to output variants, or --count')
    #
    var_diff = set()
    var_union = set()
    var_intersect = set()
    cur = proj.db.cursor()
    if args.type in [None, 'variant']:
        variants = []
        for table in args.tables:
            # read variants in tables[0]
            env.logger.info('Reading {:,} variants in {}...'.format(proj.db.numOfRows(encodeTableName(table), exact=False), table))
            cur.execute('SELECT variant_id from {};'.format(encodeTableName(table)))
            variants.append(set([x[0] for x in cur.fetchall()]))
        #
        if args.difference is not None:
            var_diff = variants[0]
            for var in variants[1:]:
               var_diff -= var 
        if args.union is not None:
            var_union = variants[0]
            for var in variants[1:]:
               var_union |=  var 
        if args.intersection is not None:
            var_intersect = variants[0]
            for var in variants[1:]:
               var_intersect &= var 
    elif args.type == 'site':
        sites = []
        record_site = args.intersection is not None or args.difference is not None
        record_all = args.union is not None
        for table in args.tables:
            # read sites in tables[0]
            env.logger.info('Reading locations of approximately {:,} sites in {}...'
                .format(proj.db.numOfRows(encodeTableName(table), exact=False), table))
            if table == 'variant':
                cur.execute('SELECT variant_id, chr, pos FROM {};'.format(encodeTableName(table)))
            else:
                cur.execute('SELECT {0}.variant_id, variant.chr, variant.pos '
                    'FROM {0} LEFT OUTER JOIN variant ON {0}.variant_id = variant.variant_id'
                .format(encodeTableName(table)))
            #
            v = defaultdict(set)
            for id, chr, pos in cur:
                if record_site:
                    v[(chr, pos)].add(id)
                if record_all:
                    var_union.add(id)
            if record_site:
                sites.append(v)
                env.logger.debug('There are {} unique sites in table {}'
                    .format(len(v), table))
        #
        if args.difference is not None:
            site_diff = set(sites[0].keys())
            for var in sites[1:]:
               site_diff -= set(var.keys())
            env.logger.debug('There are {} unique sites in table {} only'
                .format(len(site_diff), args.tables[0]))
            for k in site_diff:
                var_diff |= sites[0][k]
        if args.intersection is not None:
            site_intersect = set(sites[0].keys())
            for var in sites[1:]:
               site_intersect &= set(var.keys())
            env.logger.debug('There are {} unique sites in all tables'
                .format(len(site_intersect)))
            for k in site_intersect:
                for var in sites:
                    var_intersect |= var[k]
    else:
        # genotype
        geno = []
        NULL_to_0 = env.treat_missing_as_wildtype
        for table, sample, id in zip(args.tables, args.samples, args.sample_IDs):
            # read geno in tables[0]
            env.logger.info('Reading genotypes of sample {} of approximately {:,} geno in {}...'
                .format(sample, proj.db.numOfRows(encodeTableName(table), exact=False), table))
            cur.execute('SELECT {0}.variant_id, geno.GT FROM {0} LEFT OUTER JOIN '
                '{1}_genotype.genotype_{2} geno ON {0}.variant_id = geno.variant_id'
                    .format(encodeTableName(table), proj.name, id))
            #
            v = set()
            for id, GT in cur:
                if GT is None:
                    if NULL_to_0:
                        v.add((id, 0))
                else:
                    v.add((id, GT))
            geno.append(v)
        #
        if args.difference is not None:
            geno_diff = geno[0]
            for g in geno[1:]:
               geno_diff -= g
            env.logger.debug('There are {} genotypes in table {} only'
                .format(len(geno_diff), args.tables[0]))
            var_diff = set([x[0] for x in geno_diff])
        if args.union is not None:
            geno_union = geno[0]
            for g in geno[1:]:
               geno_union |= g 
            env.logger.debug('There are a total of {} genotypes'
                .format(len(geno_union)))
            var_union = set([x[0] for x in geno_union])
        if args.intersection is not None:
            geno_intersect = geno[0]
            for g in geno[1:]:
               geno_intersect &= g
            env.logger.debug('There are {} unique genotypes in all tables'
                .format(len(geno_intersect)))
            var_intersect = set([x[0] for x in geno_intersect])
    #
    for var, table_with_desc in [
            (var_intersect, args.intersection),
            (var_union, args.union),
            (var_diff, args.difference)]:
        if not table_with_desc:
            continue
        table = table_with_desc[0]
        if table == 'variant':
            raise ValueError('Cannot overwrite the master variant table')
        if '*' in table or '?' in table:
            env.logger.warning('Use of wildcard character * or ? in table '
                'names is not recommended because such names can be expanded '
                'to include other tables in some commands.')
        desc = table_with_desc[1] if len(table_with_desc) == 2 else ''
        if len(table_with_desc) > 2:
            raise ValueError('Only a table name and an optional table '
                'description is allowed: %s provided'.format(table_with_desc))
        if proj.db.hasTable(encodeTableName(table)):
            new_table = proj.db.backupTable(encodeTableName(table))
            env.logger.warning('Existing table {} is renamed to {}.'
                .format(table, decodeTableName(new_table)))
        proj.createVariantTable(encodeTableName(table))
        prog = ProgressBar('Writing to ' + table, len(var))
        query = 'INSERT INTO {} VALUES ({});'.format(encodeTableName(table), proj.db.PH)
        # sort var so that variant_id will be in order, which might
        # improve database performance
        for count,id in enumerate(sorted(var)):
            cur.execute(query, (id,))
            if count % 10000 == 0:
                prog.update(count)
        proj.describeTable(encodeTableName(table), desc, True, True)
        prog.done()       
        proj.db.commit()
    # count
    if args.count:
        if args.difference is not None:
            print(len(var_diff))
        elif args.union is not None:
            print(len(var_union))
        elif args.intersection is not None:
            print(len(var_intersect))

    

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
                        if re.match(table.replace('?', '.{1}').replace('*', '.*'), tbl, re.I):
                            tables.append(tbl)
                            match = True
                    if not match:
                        # * should match a table with '*' in its name.
                        env.logger.warning('Name {} does not match any existing table.'
                            .format(table))
                else:
                    tables.append(table)
            # table?
            for table in tables:
                if not proj.isVariantTable(encodeTableName(table)):
                    raise ValueError('Variant table {} does not exist.'.format(table))
            # set args.tables to its expanded version
            args.tables = tables
            #
            # type of comparison
            if args.samples:
                if args.type is not None and args.type != 'genotype':
                    raise ValueError('Cannot specify samples for non-genotype comparisons.')
                args.type = 'genotype'
                if len(args.samples) == 1:
                    raise ValueError('Please specify more than one sample to be compared.')
                if len(args.tables) == 1:
                    # automatically expand tables to multiple tables
                    args.tables = [args.tables[0]] * len(args.samples)
                elif len(args.tables) != len(args.samples):
                    raise ValueError('Please specify a variant table for all or each of the samples.')
                # check sample names
                proj.db.attach(proj.name + '_genotype')
                args.sample_IDs = []
                for name in args.samples:
                    IDs = proj.selectSampleByPhenotype("sample_name = '{}'".format(name))
                    if len(IDs) == 0:
                        raise ValueError("No sample with name '{}' is identified.".format(name))
                    elif len(IDs) > 1:
                        raise ValueError("More than one sample with name '{}' is identified.".format(name))
                    if not 'GT' in proj.db.getHeaders('{}_genotype.genotype_{}'.format(proj.name, IDs[0])):
                        raise ValueError('Sample {} does not have any genotype.'.format(name))
                    args.sample_IDs.append(IDs[0])
            elif args.type == 'genotype':
                raise ValueError('Please specify samples to be compared for genotype comparison.')
            if len(args.tables) == 1:
                raise ValueError('Please specify at least two tables to compare.')
            # 
            if args.B_diff_A or args.A_diff_B or args.A_and_B or args.A_or_B:
                raise ValueError('Options B_diff_A, A_diff_B, A_and_B and A_or_B are deprecated.')
            if args.intersection is not None or args.union is not None or args.difference is not None:
                compareTables(proj, args)
            elif args.count:
                if args.type in [None, 'variant']:
                    countVariantDifference(proj, args)
                elif args.type == 'site':
                    countSiteDifference(proj, args)
                elif args.type == 'genotype':
                    countGenotypeDifference(proj, args)
            else:
                env.logger.warning('No action parameter is specified. Nothing to do.')
                return
    except Exception as e:
        env.logger.error(e)
        sys.exit(1) 

