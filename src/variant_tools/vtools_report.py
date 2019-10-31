#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
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
import shlex
import argparse
import subprocess
from variant_tools.utils import RefGenome, env, OS_ENV, PrettyPrinter, \
     convertDoubleQuote, withinPseudoAutoRegion, encodeTableName, expandRegions, \
     genesInRegions, dissectGene, getRNASequence, getProteinSequence
from variant_tools.project import Project
# if sys.version_info.major == 2:
#     from variant_tools.assoTests_py2 import calculateInbreedingCoef
# else:
#     from variant_tools.assoTests_py3 import calculateInbreedingCoef
from variant_tools.assoTests import calculateInbreedingCoef
from variant_tools.plot import CTHEME, plotAssociation, rhist, rdot, rbox, \
     resolvePlotFilename

# save the command line that has been processed by the shell. It might not work as shown.
env.command_line = 'vtools_report ' + ' '.join([x if (x.isalnum() or x.startswith('-')) else repr(x) for x in sys.argv[1:]])
#
# These functions are utility functions used by all reports
#
def addCommonArgs(parser):
    parser.add_argument('-v', '--verbosity', choices=['0', '1', '2', '3'],
        help='''Output error and warning (0), info (1), debug (2) and trace (3) information
            to standard output (default to 1).'''),

def getoutput(cmd):
    if type(cmd) is str:
        cmd = shlex.split(cmd)
    # add -v option if unspecified
    if True not in [x.startswith('-v') for x in cmd]:
        cmd.append('-v{}'.format(env.verbosity))
    # call vtools and return its output
    env.logger.debug('Running: ' + ' '.join([x if (x.isalnum() or x.startswith('-')) else '"' + x.replace('"', '\"') + '"' for x in cmd]))
    # '.' is added to $PATH so that command (vtools) that is in the current directory
    # can be executed.
    return subprocess.check_output(cmd, env=OS_ENV).strip().decode()

def outputToFile(cmd, filename):
    # execute a command and send its output to a file
    if True not in [x.startswith('-v') for x in cmd]:
        cmd.append('-v{}'.format(env.verbosity))
    # call vtools and return its output
    env.logger.debug('Running: ' + ' '.join([x if (x.isalnum() or x.startswith('-')) else '"' + x.replace('"', '\"') + '"' for x in cmd]))
    # '.' is added to $PATH so that command (vtools) that is in the current directory
    # can be executed.
    with open(filename, 'w') as output:
        subprocess.call(cmd, stdout=output, env=OS_ENV)

#
# These functions call vtools to extract information that are needed by other
# reports.
#
def genotypeOfSample(name=None, ID=None, cond='GT != 0'):
    '''return a list of genotypes for a sample with specified name '''
    with Project(verbosity=env.verbosity) as proj:
        proj.db.attach('{}_genotype'.format(proj.name))
        cur = proj.db.cursor()
        if name is not None:
            cur.execute("SELECT sample_id FROM sample WHERE sample_name = '{}' "
                .format(name))
            IDs = cur.fetchall()
            if len(IDs) == 0:
                env.logger.error('No sample with name {} is identified'.format(name))
                return None
            elif len(IDs) > 1:
                env.logger.error('{} does not identify a sample uniquely'.format(name))
                return None
            else:
                ID = IDs[0][0]
        elif ID is None:
            env.logger.error('Please specify either name or ID of sample')
            return None
        #
        cur.execute("SELECT variant_id, GT FROM {}_genotype.genotype_{} {}"
            .format(proj.name, ID, 'WHERE {}'.format(cond) if cond else ''))
        return {x[0]:x[1] for x in cur.fetchall()}

def getSamples(samples, group_by=[]):
    '''return a generator that returns samples group by group. To use this functions, call it like
    for group, IDs in getSamples(condi, group_by):
        # do something with group and IDs.
    '''
    output = getoutput(['vtools', 'execute', 'SELECT {} FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id {} {};'\
        .format('sample_id, sample_name, {}'.format(', '.join(group_by)) if group_by else 'sample_id, sample_name',
            'WHERE {}'.format(' AND '.join(samples)) if samples else '',
            'GROUP BY {}'.format(', '.join(group_by)) if group_by else '')])
    if not output.strip():
        sys.exit('There are no available samples to analyze {}.\n'.format(('given {}'.format(samples)) if samples else ''))
    # we return results batch by batch as a generator
    IDs = []
    names = []
    group = ()
    for line in output.split('\n'):
        fields = [x.strip() for x in line.split('\t')]
        ID = fields[0]
        name = fields[1]
        g = tuple(fields[2:])
        if group == ():
            group = g
        if g == group:
            IDs.append(ID)
            names.append(name)
        else:
            yield (group, IDs, names)
            # start a new group
            group = g
            IDs = [ID]
            names = [name]
    yield (group, IDs, names)

#
# Command trans_ratio
#
def transRatioArguments(parser):
    parser.add_argument('-n', '--num_field', '--num-field', required=True,
        help='''Name of the field that holds sample variant count, which is the field name for
            command 'vtools update table --from_stat "num=#(alt)"'.''')
    parser.add_argument('--group_by', '--group-by', nargs='*', default=[],
        help='''Output transition/transversion rate for groups of variants. e.g. --group_by
            num for each sample variant frequency group.''')
    parser.add_argument('table',
        help='''Variant table for which transversion/transversion mutants are counted.''')

def transRatio(args):
    #
    prt = PrettyPrinter()
    if args.group_by:
        prt.write(args.group_by + ['num_of_transition', 'num_of_transversion', 'ratio'])
        transition = getoutput(['vtools', 'select', args.table,
            "((ref='A' AND alt='G') OR (ref='G' AND alt='A') OR (ref='C' AND alt='T') OR (ref='T' AND alt='C'))",
            '--output'] + args.group_by + ['sum({})'.format(args.num_field), '-d', '\t', '--group_by'] + args.group_by)
        transversion = getoutput(['vtools', 'select', args.table,
            "((ref='A' AND alt='C') OR (ref='C' AND alt='A') OR (ref='G' AND alt='T') OR " +
            " (ref='T' AND alt='G') OR (ref='A' AND alt='T') OR (ref='T' AND alt='A') OR " +
            " (ref='C' AND alt='G') OR (ref='G' AND alt='C'))", '--output'] + args.group_by +
            ['sum({})'.format(args.num_field), '-d', '\t', '--group_by'] +  args.group_by)
        values = {}
        for item in transition.split('\n'):
            g, v = item.rsplit('\t', 1)
            values[g] = [v, '0']
        for item in transversion.split('\n'):
            g, v = item.rsplit('\t', 1)
            if g in values:
                values[g][1] = v
            else:
                values[g] = ['0', v]
        keys = list(values.keys())
        keys.sort()
        for k in sorted(keys):
            prt.write([k, '{:,}'.format(int(values[k][0])), '{:,}'.format(int(values[k][1])),
                '{:.5f}'.format(int(values[k][0])/float(values[k][1]) if int(values[k][1]) > 0 else 0)])
    else:
        prt.write(['num_of_transition', 'num_of_transversion', 'ratio'])
        transition = int(getoutput(['vtools', 'select', args.table,
            "(ref='A' AND alt='G') OR (ref='G' AND alt='A') OR (ref='C' AND alt='T') OR (ref='T' AND alt='C')",
            '--output', 'sum({})'.format(args.num_field)]))
        transversion = int(getoutput(['vtools', 'select', args.table,
            "(ref='A' AND alt='C') OR (ref='C' AND alt='A') OR (ref='G' AND alt='T') OR " +
            "(ref='T' AND alt='G') OR (ref='A' AND alt='T') OR (ref='T' AND alt='A') OR " +
            "(ref='C' AND alt='G') OR (ref='G' AND alt='C')",
            '--output', 'sum({})'.format(args.num_field)]))
        prt.write(['{:,}'.format(transition), '{:,}'.format(transversion),
            '{:.5f}'.format(transition / float(transversion) if transversion != 0 else 0)])
    prt.write_rest()

#
# Command avg_depth
#
def avgDepthArguments(parser):
    parser.add_argument('-n', '--num_field', '--num-field', required=True,
        help='''Name of the field that holds sample variant count, which is the field name for
            command 'vtools update table --from_stat "num=#(alt)"'.''')
    parser.add_argument('-d', '--depth_field', '--depth-field', required=True,
        help='''Name of the field that holds average depth of each variant, which is the field
            name for command 'vtools update table --from_stat "meanDP=avg(DP_geno)"'.''')
    parser.add_argument('--group_by', '--group-by', nargs='*', default=[],
        help='''Output average depth for each group, for example,
            '--group_by NUM_FIELD to output depth for each sample variant frequency (count).''')
    parser.add_argument('table',
        help='''Variant table for which average depth are calculated.''')

def avgDepth(args):
    print(('{}num_of_variant\taverage_depth'.format(''.join([x+'\t' for x in args.group_by]))))
    print((getoutput(['vtools', 'output', args.table] + args.group_by +
        ['COUNT(1)', 'SUM({0}*{1})/SUM({0})'.format(args.num_field, args.depth_field)] +
        (['--group_by'] + args.group_by if args.group_by else []))))

#
# Command variant_stat
#
def variantStatArguments(parser):
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"'). If this parameter is specified without
            a value, variants belonging to any of the samples will be counted.
            If this parameter is left unspecified, all variants, including those
            that do not belong to any samples will be counted.''')
    parser.add_argument('-g', '--group_by', '--group-by', nargs='*', default=[],
        help='''Group samples by certain conditions such as 'aff=1'. A common
            usage is to group variants by 'filename' and 'sample_name' so that
            variant statistics are outputted for each sample.''')
    parser.add_argument('table',
        help='''Variant table for which variant metrics are calculated.''')

def variantStat(args):
    # 1) Get samples based on the conditional parameter --samples
    #    Exit the report if there are no samples to analyze.
    print(('\t'.join(list(args.group_by) + ['num_sample', 'num_snps', 'num_insertions', 'num_deletions',
        'num_substitutions', 'min_insertion_size', 'avg_insertion_size', 'max_insertion_size',
        'min_deletion_size', 'avg_deletion_size', 'max_deletion_size'])))
    if not args.samples and not args.group_by:   # no sample is specified:
        # 2a) Get the counts for snps and substitutions:
        #     command: vtools select __tmp_vs "ref != '-' and alt != '-' and (length(ref) = 1 and length(alt) = 1)" --samples 'sample_id in $sample_ids' --count
        num_snps = getoutput(['vtools', 'select', args.table, "ref != '-'", "alt != '-'", "(length(ref) = 1 and length(alt) = 1)",
            '--count'])
        num_substitutions = getoutput(['vtools', 'select', args.table, "ref != '-'", "alt != '-'", "(length(ref) > 1 or length(alt) > 1)",
            '--count'])
        #
        # 2b) Get the metrics to characterize the insertions
        #     command: vtools select variant "ref = '-'" --samples 'sample_id IN $sample_ids' --output 'count(alt)', 'avg(length(alt))' 'min(length(alt))' 'max(length(alt))'
        num_insertions, avg_insertion_size, min_insertion_size, max_insertion_size = getoutput(
            ['vtools', 'select', args.table, "ref='-'",
            '--output', 'count(alt)', 'avg(length(alt))', 'min(length(alt))', 'max(length(alt))']).split('\t')
        #
        # 2c) Get the metrics to characterize the deletions
        #     command: vtools select variant "alt = '-'" --samples 'sample_id IN $sample_ids' --output 'count(ref)', 'avg(length(ref))' 'min(length(ref))' 'max(length(ref))'
        num_deletions, avg_deletion_size, min_deletion_size, max_deletion_size = getoutput(
            ['vtools', 'select', args.table, "alt='-'",
            '--output', 'count(ref)', 'avg(length(ref))', 'min(length(ref))', 'max(length(ref))']).split('\t')
        num_samples = getoutput(['vtools', 'execute', 'SELECT COUNT(sample_id) FROM sample'])
        #
        print(('\t'.join([num_samples, num_snps, num_insertions, num_deletions,
            num_substitutions, min_insertion_size, avg_insertion_size, max_insertion_size,
            min_deletion_size, avg_deletion_size, max_deletion_size])))
    else:
        for group, sample_ids, sample_names in getSamples(args.samples, args.group_by):
            #
            # 2a) Get the counts for snps and substitutions:
            #     command: vtools select __tmp_vs "ref != '-' and alt != '-' and (length(ref) = 1 and length(alt) = 1)" --samples 'sample_id in $sample_ids' --count
            num_snps = getoutput(['vtools', 'select', args.table, "ref != '-'", "alt != '-'", "(length(ref) = 1 and length(alt) = 1)",
                '--samples', 'sample_id IN ({})'.format(','.join(sample_ids)), '--count'])
            num_substitutions = getoutput(['vtools', 'select', args.table, "ref != '-'", "alt != '-'", "(length(ref) > 1 or length(alt) > 1)",
                '--samples', 'sample_id IN ({})'.format(','.join(sample_ids)), '--count'])
            #
            # 2b) Get the metrics to characterize the insertions
            #     command: vtools select variant "ref = '-'" --samples 'sample_id IN $sample_ids' --output 'count(alt)', 'avg(length(alt))' 'min(length(alt))' 'max(length(alt))'
            num_insertions, avg_insertion_size, min_insertion_size, max_insertion_size = getoutput(
                ['vtools', 'select', args.table, "ref='-'", '--samples', 'sample_id IN ({})'.format(','.join(sample_ids)),
                '--output', 'count(alt)', 'avg(length(alt))', 'min(length(alt))', 'max(length(alt))']).split('\t')
            #
            # 2c) Get the metrics to characterize the deletions
            #     command: vtools select variant "alt = '-'" --samples 'sample_id IN $sample_ids' --output 'count(ref)', 'avg(length(ref))' 'min(length(ref))' 'max(length(ref))'
            num_deletions, avg_deletion_size, min_deletion_size, max_deletion_size = getoutput(
                ['vtools', 'select', args.table, "alt='-'", '--samples', 'sample_id IN ({})'.format(','.join(sample_ids)),
                '--output', 'count(ref)', 'avg(length(ref))', 'min(length(ref))', 'max(length(ref))']).split('\t')
            #
            print(('\t'.join(list(group) + [str(len(sample_ids)), num_snps, num_insertions, num_deletions,
                num_substitutions, min_insertion_size, avg_insertion_size, max_insertion_size,
                min_deletion_size, avg_deletion_size, max_deletion_size])))


#
# command discordance_rate
#
def discordanceRateArguments(parser):
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('--variants', default='variant', metavar='TABLE',
        help='''Limit variants to specified variant table. Default to all variants.''')
    parser.add_argument('--genotypes', nargs='*', default=[],
        help='''Limiting genotypes from samples that match conditions that
            involves genotype fields (e.g. filter by quality score, with fields
            shown in command 'vtools show genotypes'). If a variant is filtered
            for one sample but not another, it will be included if runtime option
            $treat_missing_as_wildtype is set to True, and discarded otherwise.''')

def discordanceRate(args):
    #
    # It is possible to use SQL to speed up this funciton. E.g. instead of
    # getting all variants and their genotypes, use something like
    #   SELECT count(a.variant_id) FROM genotype_1 a, genotype_2 b ON a.variant_id = b.variant_id;
    #
    for grp, sample_ids, sample_names in getSamples(args.samples):
        print(('\t'.join(sample_names)))
        # get list of variant for each sample
        rates = []
        treat_missing_as_wildtype = getoutput(['vtools', 'show', 'runtime_option', 'treat_missing_as_wildtype'])
        # geno_condition
        if not args.genotypes and not args.variants:
            geno_condition = ''
        else:
            geno_condition = []
            if args.genotypes:
                geno_condition.append(','.join(['({})'.format(x) for x in args.genotypes]))
            if args.variants:
                geno_condition.append('variant_id in (select variant_id from {})'.format(args.variants))
            geno_condition = ' AND '.join(geno_condition)
        for i, (id_i, name_i) in enumerate(zip(sample_ids, sample_names)):
            rates.append([0] * len(sample_ids))
            output_items = []
            var_i = genotypeOfSample(ID=id_i, cond=geno_condition)
            for j, (id_j, name_j) in enumerate(zip(sample_ids, sample_names)):
                if i > j:
                    output_items.append('{:.5f}'.format(float(rates[j][i][0]) / (rates[j][i][0] + rates[j][i][1])))
                elif i == j:
                    output_items.append('0/{0}'.format(len(var_i)))
                else:
                    var_j = genotypeOfSample(ID=id_j, cond=geno_condition)
                    if treat_missing_as_wildtype == 'False':
                        # find variant that exist in both samples
                        exist_in_both = [item for item in list(var_i.keys()) if item in var_j]
                        # find variant that has different variant
                        differ_in_two = [item for item in exist_in_both if var_i[item] != var_j[item]]
                        rates[i][j] = (len(differ_in_two), len(exist_in_both))
                    else:
                        exist_in_any = set(var_i.keys()) | set(var_j.keys())
                        # remove 0 items to facilitate search
                        var_i = {x:y for x,y in list(var_i.items()) if y != 0}
                        var_j = {x:y for x,y in list(var_j.items()) if y != 0}
                        #
                        differ_in_two = [item for item in exist_in_any if (not item in var_i) \
                            or (not item in var_j) or (var_i[item] != var_j[item])]
                        rates[i][j] = (len(differ_in_two), len(exist_in_any))
                    output_items.append('{}/{}'.format(rates[i][j][0], rates[i][j][1]))
            print(('\t'.join(output_items)))

#
# command transcript
#

def transcriptArguments(parser):
    parser.add_argument('regions', nargs='+', help=''''One or more chromosome regions in
        the format of chr:start-end (e.g. chr21:33,031,597-33,041,570), Field:Value
        from a region-based annotation database (e.g. refGene.name2:TRIM2 or
        refGene_exon.name:NM_000947), or set options of several regions (&, |, -,
        and ^ for intersection, union, difference, and symmetric difference).''')
    parser.add_argument('--build', nargs='?',
        help='''Output sequence at specified build of reference genome. The primary
            reference genome of the project will be used if by default.''')
    parser.add_argument('--strand_only', '--strand-only', action='store_true',
        help='''Only output strand of genes that covers the region.''')
    parser.add_argument('--first_transcript', '--first-transcript', action='store_true',
        help='''If set, only display the first transcript of RNA or Protein sequence''')
    parser.add_argument('--zero_based',  '--zero-based', default=False, action='store_true',
        help='''If set, user input is zero based and will be translated to 1-based
            coordinates before query. The output is always 1-based''')


def transcript(args):
    #
    with Project(verbosity=0) as proj:
        if args.build is None:
            args.build = proj.build
        #
        regions = expandRegions(';'.join(args.regions), proj, mergeRegions=False,
            zeroBased=args.zero_based)
        # if a region is specified
        printRNAInfo(args, regions)


def printRNAInfo(args, regions):
    # get chromosomal region
    with Project(verbosity=0) as proj:
        genes = genesInRegions(regions, proj)
        for gene in genes:
            structure = dissectGene(gene, proj)
            strand = structure['strand']
            if args.strand_only:
                if args.verbosity == '0':
                    print(strand)
                else:
                    print(('{}\t{}'.format(gene, strand)))
            else:
                coding = structure['coding']
                if not coding:
                    if args.verbosity is None or int(args.verbosity) > 0:
                        print(('>mRNA|{}|Non-coding'.format(gene)))
                    continue
                else:
                    chr = coding[0][0]
                if args.verbosity is None or int(args.verbosity) > 0:
                    if strand == '+':
                        print(('>mRNA|{}|{}:{} ({} bp)'.format(gene,
                            chr if chr.startswith('chr') else 'chr' + chr,
                            ','.join(['{}-{}'.format(x,y) for ch,x,y in coding]),
                            sum([y-x+1 for ch,x,y in coding]))))
                    else:
                        print(('>mRNA|{}|{}:{} ({} bp on reverse strand)'.format(gene,
                            chr if chr.startswith('chr') else 'chr' + chr,
                            ','.join(['{}-{}'.format(y,x) for ch,x,y in coding[::-1]]),
                            sum([y-x+1 for ch,x,y in coding]))))
            if args.first_transcript:
                break


#
# command sequence
#
def sequenceArguments(parser):
    parser.add_argument('regions', nargs='+', help='''One or more chromosome regions in
        the format of chr:start-end (e.g. chr21:33,031,597-33,041,570), Field:Value
        from a region-based annotation database (e.g. refGene.name2:TRIM2 or
        refGene_exon.name:NM_000947), or set options of several regions (&, |, -,
        and ^ for intersection, union, difference, and symmetric difference).
        Several regions will be printed if the name matches more than one regions.
        Chromosome positions are 1-based and are inclusive at both ends so the chromosome
        region has a length of end-start+1 bp. A reversed complementary sequence will be
        outputted if starting position is after the ending position.''')
    parser.add_argument('--build', nargs='?',
        help='''Output sequence at specified build of reference genome. The primary
            reference genome of the project will be used if by default.''')
    parser.add_argument('--numbered', nargs='?', const='left', choices=['left', 'right'],
        help='''If specified, add position of the first or last basepair of each line to
            the left or right of the line, and insert a space at every 10 basepair''')
    parser.add_argument('--char_per_line', '--char-per-line', type=int,
        help='''Number of characters (excluding space and leading numbers) per line. Default to
            70 in regular and 60 in numbered format.''')
    parser.add_argument('--transcribe', nargs='*', metavar='GENE',
        help='''Transcribe DNA sequence into RNA sequence. variant tools will look
            in the refGene database, find all genes that overlap with the region, locate
            exons and coding regions, transcribe the DNA sequence to RNA sequence (basically
            discard introns, and complement if on reverse strand). The complete mRNA sequence
            will be printed regardless of the bounaries of specified regions. If one or more
            names (refGene.name) are specified, only specified genes will be translated.''')
    parser.add_argument('--TtoU', action='store_true',
        help='''Print U for T for RNA sequence. Otherwise use T.''')
    parser.add_argument('--translate', nargs='*', metavar='GENE',
        help='''Translate DNA sequence into protein sequence. variant tools will look
            in the refGene database, find all genes that overlap with the region, locate
            exons and coding regions, transcribe and translate the DNA sequence to protein
            sequence. The complete protein sequence will be printed regardless of the
            boundaries of specified regions. If one or more names (refGene.name) are
            specified, only specified genes will be translated.''')
    parser.add_argument('--mark', nargs='*',
        help='''Mark a location (--mark chr pos), a variant (--mark chr pos ref alt),
            a region (e.g. refGene_exon.name:NM_000947), or a particular sequence (e.g. TCGGA)
            in red in the output. If a variant is specified, the changed nucleotide or amino
            acid will be printed. Currently only single nucleotide polymorphisms are supported.''')
    parser.add_argument('--mark_complement', '--mark-complement',  default=False, action='store_true',
        help='''If set, also try to mark the complement of the specified sequence''')
    parser.add_argument('--mark_reverse', '--mark-reverse', default=False, action='store_true',
        help='''If set, also try to mark the reverse sequence of the specified sequence.
            If both mark_complemnt and mark_reverse are set, four different sequences
            will be searched.''')
    parser.add_argument('--first_transcript', '--first-transcript', action='store_true',
        help='''If set, only display the first transcript of RNA or Protein sequence''')
    parser.add_argument('--show_transcript', '--show-transcript', action='store_true',
        help='''Put transcript name before transcript''')
    parser.add_argument('--marked_region', '--marked-region', nargs='*', default=None, type=int,
        help='''If set to a number or pair of number, print only n bp to the left and
            m (m=n if only one number is specified) of the marked region. The sequence
            itself is no longer marked. This option is not yet supported in all
            combinations of options. ''')
    parser.add_argument('--hide_unmatched', '--hide-unmatched', default=False, action='store_true',
        help='''If set, only display regions with marked variants or sequences''')
    parser.add_argument('--zero_based', '--zero-based', default=False, action='store_true',
        help='''If set, user input is zero based and will be translated to 1-based
            coordinates before query. The output is always 1-based''')

complement_table = {
        'A': 'T',
        'T': 'A',
        'a': 't',
        't': 'a',
        'U': 'A',
        'u': 'a',
        'G': 'C',
        'C': 'G',
        'g': 'c',
        'c': 'g',
        'N': 'N',
        'n': 'n',
    }

def sequence(args):
    #
    with Project(verbosity=0) as proj:
        if args.build is None:
            args.build = proj.build
        #
        regions = expandRegions(';'.join(args.regions), proj, mergeRegions=False,
            zeroBased=args.zero_based)
        # if a region is specified
        args.mark_sequence = ''
        if args.mark:
            if len(args.mark) == 1 and ':' in args.mark[0]:
                regs = expandRegions(args.mark[0], proj, mergeRegions=True,
                    zeroBased=args.zero_based)
                args.mark = []
                for r in regs:
                    args.mark.extend([r[0],x] for x in range(r[1], r[2]+1))
            elif len(args.mark) == 1:
                args.mark_sequence = args.mark[0]
                args.mark = []
            else:
                # only one position or variant will be marked
                if len(args.mark) > 2:
                    if args.mark[2] == '-' or args.mark[3] == '-' or len(args.mark[2]) != 1 or len(args.mark[3]) != 1:
                        env.logger.warning('Get protein sequence does not support indels yet.')
                args.mark = [[args.mark[0][3:] if args.mark[0].startswith('chr') else args.mark[0], int(args.mark[1])] + args.mark[2:]]
    #
    if not args.build:
        raise RuntimeError('Failed to get build of reference from the current project.')
    if args.translate is not None:
        printProteinSequence(args, regions)
        if args.transcribe is not None:
            raise ValueError('Please specify only one of the options of transcribe and translate')
    elif args.transcribe is not None:
        printRNASequence(args, regions)
    else:
        printDNASequence(args, regions)


def markSequence(seq):
    '''Given a string, mark small characters or specified sequence with color.'''
    return ''.join([('\033[41m' + x.upper() + '\033[0m') if x.islower() else x for x in seq])


def printDNASequence(args, regions):
    ref = RefGenome(args.build)
    # get chromosomal region
    for chr, start, end, comment in regions:
        # get reference seq
        seq = ref.getSequence(chr, start, end)
        # right now, we only support marking one locus with chr pos, or chr pos ref alt
        if args.mark:
            seq = [x for x in seq]
            for m in args.mark:
                if m[0] == chr and m[1] >= start and m[1] <= end:
                    idx = m[1] - start
                    if len(m) > 2:
                        if m[2] != seq[idx]:
                            env.logger.warning('Reference alleles mismatch (chr {}, pos {}, ref {}, mutant ref {})'
                                .format(m[0], m[1], seq[idx], m[2]))
                        seq[idx] = m[3].lower()
                    else:
                        seq[idx] = seq[idx].lower()
                else:
                    env.logger.debug('Failed to mark mutant {}'.format(m))
        #
        if 'reverse' in comment:
            start, end = end, start
            reversed_seq = True
            seq = ''.join([complement_table[x] for x in seq[::-1]])
        else:
            reversed_seq = False
            seq = ''.join(seq)
        if args.mark_sequence:
            seq = seq.replace(args.mark_sequence.upper(), args.mark_sequence.lower())
            if args.mark_complement:
                comp = ''.join([complement_table[x] for x in args.mark_sequence])
                seq = seq.replace(comp.upper(), comp.lower())
            if args.mark_reverse:
                rev = args.mark_sequence[::-1]
                seq = seq.replace(rev.upper(), rev.lower())
            if args.mark_complement and args.mark_reverse:
                comprev = ''.join([complement_table[x] for x in args.mark_sequence][::-1])
                seq = seq.replace(comprev.upper(), comprev.lower())
        unmatched = not any([x.islower() for x in seq])
        if args.verbosity is None or int(args.verbosity) > 1:
            print(('>ref|{}|{}:{}-{} {} {}'.format(args.build, chr if chr.startswith('chr') else 'chr' + chr,
                start, end, '({})'.format(comment) if comment else '',
                'unmatched' if unmatched and args.hide_unmatched else '')))
        if unmatched and args.hide_unmatched:
            continue
        # break into pieces of 70 bp
        if args.numbered:
            char_per_line = 60 if args.char_per_line is None else args.char_per_line
            block_per_line = char_per_line // 10
            if block_per_line * 10 != char_per_line:
                raise ValueError('char_per_line should be a multiple of 10 in numbered format')
            if args.numbered == 'left':
                fmt = '{{:>{0}}} '.format(len(str(end)))
                print(('\n'.join([
                    fmt.format(end - i if reversed_seq else start + i) + ' '.join([
                        markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)])
                    for i in range(0, len(seq), char_per_line)])))
            else:
                fmt = ' {{:>{0}}}'.format(len(str(end)))
                last_line = len(seq) // char_per_line * char_per_line
                residue = last_line + char_per_line - len(seq)
                print(('\n'.join([
                    ' '.join([markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)]) +
                    ((' ' * residue + fmt.format(start if reversed_seq else end)) if i == last_line \
                        else fmt.format(end - i - char_per_line + 1 if reversed_seq else start + i + char_per_line - 1)) \
                    for i in range(0, len(seq), char_per_line)])))
        elif args.marked_region is not None:
            # find marked regions
            if len(args.marked_region) == 0:
                l = r = 0
            elif len(args.marked_region) == 1:
                l = r = int(args.marked_region[0])
            elif len(args.marked_region) == 2:
                l, r = list(map(int, args.marked_region))
            else:
                raise ValueError('Value of option marked_region should be an integer or a pair of integers')
            import re
            # find marked regions
            pieces = re.split('([a-z]+)', seq)
            if len(pieces) == 1:
                print()
            elif len(pieces) == 3:
                print((('' if l == 0 else pieces[0][-l:]) + pieces[1].upper() + pieces[2][:r]))
            else:
                raise ValueError('More than one marked region')
        else:
            char_per_line = 70 if args.char_per_line is None else args.char_per_line
            print(('\n'.join([markSequence(seq[i:i+char_per_line]) for i in range(0, len(seq), char_per_line)])))


def printRNASequence(args, regions):
    # get chromosomal region
    with Project(verbosity=0) as proj:
        genes = genesInRegions(regions, proj)
        if len(args.transcribe) > 0:
            genes = sorted(list(set(genes) & set(args.transcribe)))
        for gene in genes:
            structure = dissectGene(gene, proj)
            strand = structure['strand']
            coding = structure['coding']
            if not coding:
                if args.verbosity is None or int(args.verbosity) > 0:
                    print(('>mRNA|{}|Non-coding'.format(gene)))
                continue
            else:
                chr = coding[0][0]
            seq = getRNASequence(structure, mutants=args.mark, TtoU=args.TtoU)
            if args.mark_sequence:
                seq = seq.replace(args.mark_sequence.upper().replace('T', 'U'),
                    args.mark_sequence.lower().replace('t', 'u'))
                if args.mark_complement:
                    comp = ''.join([complement_table[x] for x in args.mark_sequence])
                    seq = seq.replace(comp.upper().replace('T', 'U'), comp.lower().replace('t', 'u'))
                if args.mark_reverse:
                    rev = args.mark_sequence[::-1]
                    seq = seq.replace(rev.upper().replace('T', 'U'), rev.lower().replace('t', 'u'))
                if args.mark_complement and args.mark_reverse:
                    comprev = ''.join([complement_table[x] for x in args.mark_sequence][::-1])
                    seq = seq.replace(comprev.upper().replace('T', 'U'), comprev.lower().replace('t', 'u'))
            unmatch = not any([x.islower() for x in seq])
            if args.verbosity is None or int(args.verbosity) > 0:
                if strand == '+':
                    print(('>mRNA|{}|{}:{} ({} bp) {}'.format(gene,
                        chr if chr.startswith('chr') else 'chr' + chr,
                        ','.join(['{}-{}'.format(x,y) for ch,x,y in coding]),
                        sum([y-x+1 for ch,x,y in coding]),
                        'unmatch' if unmatch and args.hide_unmatched else '')))
                else:
                    print(('>mRNA|{}|{}:{} ({} bp on reverse strand) {}'.format(gene,
                        chr if chr.startswith('chr') else 'chr' + chr,
                        ','.join(['{}-{}'.format(y,x) for ch,x,y in coding[::-1]]),
                        sum([y-x+1 for ch,x,y in coding]),
                        'unmatch' if unmatch and args.hide_unmatched else '')))
            if unmatch and args.hide_unmatched:
                continue
            # break into pieces of 70 bp
            if args.numbered:
                char_per_line = 60 if args.char_per_line is None else args.char_per_line
                block_per_line = char_per_line // 10
                if block_per_line * 10 != char_per_line:
                    raise ValueError('char_per_line should be a multiple of 10 in numbered format')
                if args.numbered == 'left':
                    fmt = '{{:>{0}}} '.format(len(str(len(seq))))
                    print(('\n'.join([
                        fmt.format(i+1) + ' '.join([markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)])
                        for i in range(0, len(seq), char_per_line)])))
                else:
                    fmt = ' {{:>{0}}}'.format(len(str(len(seq))))
                    last_line = len(seq) // char_per_line * char_per_line
                    residue = last_line + char_per_line - len(seq)
                    print(('\n'.join([
                        ' '.join([markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)]) +
                        ((' ' * residue + fmt.format(len(seq))) if i == last_line \
                            else fmt.format(i + char_per_line)) \
                        for i in range(0, len(seq), char_per_line)])))
            elif args.marked_region is not None:
                # find marked regions
                if len(args.marked_region) == 0:
                    l = r = 0
                elif len(args.marked_region) == 1:
                    l = r = int(args.marked_region[0])
                elif len(args.marked_region) == 2:
                    l, r = list(map(int, args.marked_region))
                else:
                    raise ValueError('Value of option marked_region should be an integer or a pair of integers')
                import re
                # find marked regions
                pieces = re.split('([a-z]+)', seq)
                if len(pieces) == 1:
                    if args.show_transcript:
                        print('No marked region')
                elif len(pieces) == 3:
                    print(((gene.ljust(15) if args.show_transcript else '') + ('' if l == 0 else pieces[0][-l:]) + pieces[1].upper() + pieces[2][:r]))
                else:
                    raise ValueError('More than one marked region')
            else:
                char_per_line = 70 if args.char_per_line is None else args.char_per_line
                print(('\n'.join([markSequence(seq[i:i+char_per_line]) for i in range(0, len(seq), char_per_line)])))
            if args.first_transcript:
                break


def printProteinSequence(args, regions):
    # get chromosomal region
    with Project(verbosity=0) as proj:
        genes = genesInRegions(regions, proj)
        if len(args.translate) > 0:
            genes = sorted(list(set(genes) & set(args.translate)))
        for gene in genes:
            structure = dissectGene(gene, proj)
            strand = structure['strand']
            coding = structure['coding']
            if not coding:
                print(('>protein|{}|Non-coding'.format(gene)))
                continue
            else:
                chr = coding[0][0]
            seq = getProteinSequence(structure, mutants=args.mark)
            if args.mark_sequence:
                seq = seq.upper().replace(args.mark_sequence.upper(), args.mark_sequence.lower())
            unmatch = not any([x.islower() for x in seq])
            if strand == '+':
                print(('>protein|{}|{}:{} ({} bp) {}'.format(gene,
                    chr if chr.startswith('chr') else 'chr' + chr,
                    ','.join(['{}-{}'.format(x,y) for ch,x,y in coding]),
                    sum([y-x+1 for ch,x,y in coding]),
                    'unmatch' if unmatch and args.hide_unmatched else '')))
            else:
                print(('>protein|{}|{}:{} ({} bp on reverse strand)'.format(gene,
                    chr if chr.startswith('chr') else 'chr' + chr,
                    ','.join(['{}-{}'.format(y,x) for ch,x,y in coding[::-1]]),
                    sum([y-x+1 for ch,x,y in coding]),
                    'unmatch' if unmatch and args.hide_unmatched else '')))
            if unmatch and args.hide_unmatched:
                continue
            # break into pieces of 70 bp
            if args.numbered:
                char_per_line = 60 if args.char_per_line is None else args.char_per_line
                block_per_line = char_per_line // 10
                if block_per_line * 10 != char_per_line:
                    raise ValueError('char_per_line should be a multiple of 10 in numbered format')
                if args.numbered == 'left':
                    fmt = '{{:>{0}}} '.format(len(str(len(seq))))
                    print(('\n'.join([
                        fmt.format(i+1) + ' '.join([markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)])
                        for i in range(0, len(seq), char_per_line)])))
                else:
                    fmt = ' {{:>{0}}}'.format(len(str(len(seq))))
                    last_line = len(seq) // char_per_line * char_per_line
                    residue = last_line + char_per_line - len(seq)
                    print(('\n'.join([
                        ' '.join([markSequence(seq[i+10*j:i+10*j+10]) for j in range(block_per_line)]) +
                        ((' ' * residue + fmt.format(len(seq))) if i == last_line \
                            else fmt.format(i + char_per_line)) \
                        for i in range(0, len(seq), char_per_line)])))
            else:
                char_per_line = 70 if args.char_per_line is None else args.char_per_line
                print(('\n'.join([markSequence(seq[i:i+char_per_line]) for i in range(0, len(seq), char_per_line)])))
            if args.first_transcript:
                break



#
# command inbreeding_coefficient
#
def inbreedingCoefArguments(parser):
    parser.add_argument('table',
        help='''Variants based on which individual inbreeding coefficients are evaluated.''')
    parser.add_argument('--samples', nargs='*',
        help='''Conditions based on which samples are selected to have inbreeding coefficients
              calculated. Default to all samples.''')
    parser.add_argument('--maf_field', '--maf-field', required=True,
        help='''Name of the field that holds minor allele frequency for sample variants,
            which is the field name for command
            'vtools update table --from_stat "maf_field=maf()" --samples ...'.''')
    parser.add_argument('--skip_autosome', '--skip-autosome', action='store_true',
                # help='''With this switch, variants on autosomes as well as on pseudo-autosomal
                    # regions (PAR1 & PAR2) on chrX and chrY will be ignored.''')
                help=argparse.SUPPRESS)


def inbreedingCoef(args):
    table = args.table if args.table else 'variant'
    sarg = (['--samples'] + [convertDoubleQuote(x) for x in args.samples]) if args.samples else []
    build = getoutput(['vtools', 'execute', 'SELECT value FROM project WHERE name="build"'])
    print("sample_name\tF_stat")
    for name in getoutput('vtools phenotype --output sample_name ' + ' '.join(sarg)).split('\n'):
        # geno is a list of sample genotype and MAF, with each element being [chr, pos, GT, MAF]
        geno = list(zip(*[x.split() for x in getoutput('''vtools output {0} chr pos "genotype('{1}')" {2} --na -9'''.\
                                          format(table, name, args.maf_field)).split('\n')]))
        if not args.skip_autosome:
            gt = list(map(int, geno[2]))
            maf = list(map(float, geno[3]))
        else:
            gt = []
            maf = []
            for i in range(len(geno[0])):
                if not withinPseudoAutoRegion(geno[0][i], int(geno[1][i]), build):
                    gt.append(int(geno[2][i]))
                    maf.append(float(geno[3][i]))
        print(("{}\t{}".format(name, calculateInbreedingCoef(gt, maf))))

#
# transmission
#
def transmissionArguments(parser):
    parser.add_argument('--parents', nargs=2,
        help='''Names of parents, which should uniquely identify two samples.''')
    parser.add_argument('--offspring', nargs='+',
        help='''Names of one or more offspring samples.''')
    parser.add_argument('--denovo', nargs='*',
        help='''A list of tables to store denovo variants for each offspring.
            DeNovo variants are defined as offspring variants that do not exist
            in any of the parents, including the cases when the offspring have
            different variants from what parents have at the same genomic
            locations.''')
    parser.add_argument('--recessive', nargs='*',
        help='''A list of tables to store recessive variants for each offspring.
            Recessive variants are defined as variants that are homozygous in
            offspring, and heterozygous in both parents.''')
    parser.add_argument('--inconsistent', nargs='*',
        help='''A list of tables to store variants for each offspring that
            demonstrate mendelian inconsistencies, namely variants that are not
            passed from parents to offspring in a Mendelian fashion. Examples
            of inconsistent variants include de novo variants, homozygous variants
            in offspring with only one parental carrier, wildtype offspring
            variants with a homozygous parent, heterozygous offspring variants
            with two homozygous parents, and more complicated cases when multiple
            variants appear at the same sites.''')

def transmission(args):
    # first check if the number of tables match the number of offspring
    numOff = len(args.offspring)
    if (args.recessive and len(args.recessive) != numOff) or \
       (args.denovo and len(args.denovo) != numOff) or \
       (args.inconsistent and len(args.inconsistent) != numOff):
        env.logger.error('Please specify name of a variant table for each of the {} offspring'.format(numOff))
    if not (args.recessive or args.denovo or args.inconsistent):
        return
    # get genotypes for parents
    par1 = genotypeOfSample(args.parents[0])
    if par1 is None:
        return
    env.logger.info('{} genotypes are located for parent {}'.format(len(par1), args.parents[0]))
    par2 = genotypeOfSample(args.parents[1])
    if par2 is None:
        return
    env.logger.info('{} genotypes are located for parent {}'.format(len(par2), args.parents[1]))
    for idx, off in enumerate(args.offspring):
        geno = genotypeOfSample(off)
        if geno is None:
            return
        env.logger.info('{} genotypes are located for offspring {}'.format(len(geno), off))
        with Project(verbosity=env.verbosity) as proj:
            if args.recessive:
                # recessive (-1 means heterozygous for two alternative alleles)
                par1_het = set([x for x,y in list(par1.items()) if y in (1, -1)])
                env.logger.info('{} heterozygous variants are found in parent {}'
                    .format(len(par1_het), args.parents[0]))
                par2_het = set([x for x,y in list(par2.items()) if y in (1, -1)])
                env.logger.info('{} heterozygous variants are found in parent {}'
                    .format(len(par2_het), args.parents[1]))
                off_hom = set([x for x,y in list(geno.items()) if y==2])
                env.logger.info('{} homozygous variants are found in offspring {}'
                    .format(len(off_hom), off))
                off_recessive = off_hom & par1_het & par2_het
                env.logger.info('Writing {} variants to table {}'.format(len(off_recessive), args.recessive[idx]))
                proj.createVariantTable(encodeTableName(args.recessive[idx]),
                    variants=sorted(off_recessive))
                proj.describeTable(encodeTableName(args.recessive[idx]),
                    message='Recessive variants of sample {} with parents {} and {}'.format(
                        off, args.parents[0], args.parents[1]), save_date=True)
            if args.denovo:
                off_denovo = set([x for x in list(geno.keys()) if x not in par1 and x not in par2])
                env.logger.info('Writing {} variants to table {}'.format(len(off_denovo), args.denovo[idx]))
                proj.createVariantTable(encodeTableName(args.denovo[idx]),
                    variants=sorted(off_denovo))
                proj.describeTable(encodeTableName(args.denovo[idx]),
                    message='de novo variants of sample {} with parents {} and {}'.format(
                        off, args.parents[0], args.parents[1]), save_date=True)
            if args.inconsistent:
#
# If there is only one variant at a variant site
#
# Mendelian mismatch cases are:
#
# 001
# 002 -> de novo variant
# 012
# 022 -> single parent yield homozygous offspring
# 020
# 120
# 220 -> homozygous parent yield wildtype offspring
# 221 (or 22 -1)-> double homozygous parents not yield homozygous offspring
#
# where xyz are ordered parental (x,y) and offspring (z) genotypes.
#
# The situation is much more complicated when there are multiple variants
# involved.
#
# there can be at most 6 variants at a site for 3 individuals.  Let us forget
# about these for now.
                case1 = set([x for x in list(geno.keys()) if x not in par1 and x not in par2])
                case2 = set([x for x,y in list(geno.items()) if y==2 and (x not in par1 or x not in par2)])
                case31 = set([x for x,y in list(par1.items()) if y==2 and x not in geno])
                case32 = set([x for x,y in list(par2.items()) if y==2 and x not in geno])
                case4 = set([x for x,y in list(par1.items()) if y==2 and x in par2 and par2[x]==2 and (x in geno and geno[x] != 2)])
                off_inconsistent = case1 | case2 | case31 | case32 | case4
                env.logger.info('Writing {} variants to table {}'.format(len(off_inconsistent), args.inconsistent[idx]))
                proj.createVariantTable(encodeTableName(args.inconsistent[idx]),
                    variants=sorted(off_inconsistent))
                proj.describeTable(encodeTableName(args.inconsistent[idx]),
                    message='Mendelian inconsistent variants of sample {} with parents {} and {}'.format(
                        off, args.parents[0], args.parents[1]), save_date=True)


#
# command plot_fields
#
def plotFieldsCommonArguments(parser):
    parser.add_argument('--save_data', '--save-data', metavar='FILENAME', help='''Save data to file.''')
    parser.add_argument('--save_script', '--save-script', metavar='FILENAME', help='''Save R script to file.''')
    parser.add_argument('--width', metavar='px', type=int, default=800,
                        help='''Width of plot. Default to 800.''')
    parser.add_argument('--height', metavar='px', type=int, default=600,
                        help='''Height of plot. Default to 600.''')
    hist = parser.add_argument_group('Draw histogram')
    hist.add_argument('--hist',metavar='name',
        help='''File name of the outputted figure, which can have type PDF,
            EPS, or JPG. Multiple files might be produced if more than one
            figure is produced (e.g. MyFig_$FIELD1.pdf, MyFig_$FILED2.pdf
            if MyFig.pdf is specified)''')
    hist.add_argument('--norm_curve', '--norm-curve', action='store_true',
                help='''Add a normal distribution N(mean, stdev) density curve to the histogram.''')
    dot = parser.add_argument_group('''Draw dot plot. Allow up to 3 input fields: for single input
    field, the values will be plotted on y-axis with index being x-axis; for two input fields, the first
    field will be plotted on x-axis and the second field on y-axis; for three input fields, values of
    the third input field is represented by color of the dots.''')
    dot.add_argument('--dot', metavar='name',
        help='''File name of the outputted figure, which can have type PDF, EPS, or JPG.''')
    dot.add_argument('--dot_size', '--dot-size', metavar='pt', type=float, default=2.0,
                help='''Size of dots. Default is 2.0''')
    dot.add_argument('--discrete_color', '--discrete-color', type=str, choices=CTHEME,
                     help='''If specified, the third field of input will be treated as "factor" data.''')
    box = parser.add_argument_group('''Draw box plot. Allow one or more input fields and produce
    box plot for all fields in one graph. With --stratify option, box plot will be generated for field
    in different strata, if there is only one input field, or for the first field in different strata of
    the second field, if there are two input fields.''')
    box.add_argument('--box', metavar='name',
        help='''File name of the outputted figure, which can have type PDF, EPS, or JPG.''')
    box.add_argument('--stratify', metavar='C', nargs='+', type=float,
                     help='''Cutoff values to stratify data in the input field for box plot.
                     When this option is on, only one input field is allowed.''')
    box.add_argument('--outlier_size', '--outlier-size', metavar='pt', type=float, default=2.0,
                help='''Size of outlier dots. Default is 2.0''')
    box.add_argument('--color', type=str, choices=CTHEME,
                     help='''Color theme for boxes.''')

    # hist.add_argument('--title',
    #     help='''Title of the histogram. '$FIELD' in the title will be replaced
    #         by name of the field.''')
    # hist.add_argument('--group_by',
    #    help='''A field that will be used to group others. The histogram of this
    #        field will not be plotted.''')

    # cust = parser.add_argument_group('Draw plot using user-specified script.')
    # cust.add_argument('--script', nargs='+', metavar=('SCRIPT', 'OPT'),
    #    help='''Path to a user-provided script, which
    #        will be called by 'Rscript $script $name' where $name is the data file
    #        generated by this command. Additional arguments of this script will be passed
    #        directly the script.''')



def plotFieldsArguments(parser):
    parser.add_argument('fields', nargs='+', help='A list of fields that will be outputted.')
    parser.add_argument('--variants', default='variant', metavar='TABLE',
        help='''Limit value of fields to variant in specified variant table. Default to all variants.''')
    plotFieldsCommonArguments(parser)

def plotFields(args):
    env.logger.info('Gathering data for plot')
    output = getoutput(['vtools', 'output', args.variants] + args.fields + ['--na', 'NA', '--header'])
    if args.hist is not None:
        fns = resolvePlotFilename(args.hist, args.fields)
        rhist(output, fns, args.width, args.height, normcurve = args.norm_curve,
              save_data = args.save_data, save_script = args.save_script)
    if args.dot is not None:
        rdot(output, args.dot, args.width, args.height, args.dot_size, args.discrete_color,
             args.save_data, args.save_script)
    if args.box is not None:
        rbox(output, args.fields, args.stratify, args.box,args.width, args.height,
             args.outlier_size, args.color, args.save_data, args.save_script)


#
# command plot_pheno_fields
#
def plotPhenoFieldsArguments(parser):
    parser.add_argument('fields', nargs='+', help='A list of fields that will be outputted.')
    parser.add_argument('--samples', nargs='*',
        help='''Conditions based on which samples are selected. Default to all samples.''')
    plotFieldsCommonArguments(parser)

def plotPhenoFields(args):
    env.logger.info('Gathering data for plot')
    sarg = (['--samples'] + [convertDoubleQuote(x) for x in args.samples]) if args.samples else []
    output = getoutput('vtools phenotype --output ' + \
                       ' '.join(args.fields + sarg) + ' --na NA --header')
    if args.hist is not None:
        fns = resolvePlotFilename(args.hist, args.fields)
        rhist(output, fns, args.width, args.height,normcurve = args.norm_curve,
              save_data = args.save_data, save_script = args.save_script)
    if args.dot is not None:
        rdot(output, args.dot, args.width, args.height,
             args.dot_size, args.discrete_color, args.save_data, args.save_script)
    if args.box is not None:
        rbox(output, args.fields, args.stratify, args.box,args.width, args.height,
             args.outlier_size, args.color, args.save_data, args.save_script)

#
# command plot_geno_fields
#
def plotGenoFieldsArguments(parser):
    parser.add_argument('fields', nargs='+', help='A list of genotype fields that will be outputted.')
    parser.add_argument('--variants', metavar='TABLE',
        help='''Limit value of fields to variant in specified variant table. Default to all variants.''')
    parser.add_argument('--samples', nargs='*',
        help='''Conditions based on which samples are selected. Default to all samples.''')
    parser.add_argument('--genotypes', nargs='*',
        help='''Conditions based on which genotypes are selected. Default to all variants.''')
    plotFieldsCommonArguments(parser)

def plotGenoFields(args):
    output = {x:[] for x in args.fields}
    if not args.genotypes and not args.variants:
        where_clause = ''
    else:
        where_clause = []
        if args.genotypes:
            where_clause.append(','.join(['({})'.format(x) for x in args.genotypes]))
        if args.variants:
            where_clause.append('variant_id in (select variant_id from {})'.format(args.variants))
        where_clause = 'WHERE ' + ' AND '.join(where_clause)
    # get data into dict
    for g, IDs, names in getSamples(args.samples):
        for id, name in zip(IDs, names):
            env.logger.info('Gathering data from sample {}'.format(name))
            for field in args.fields:
                values = getoutput(['vtools', 'execute',
                    'SELECT {} FROM genotype.genotype_{} {}'.format(field, id, where_clause)])
                for line in values.split('\n'):
                    output[field].append(line)
    if args.hist is not None:
        fns = resolvePlotFilename(args.hist, args.fields)
        dfiles = resolvePlotFilename(args.save_data, args.fields)
        sfiles = resolvePlotFilename(args.save_script, args.fields)
        for f, p, d, s in zip(args.fields, fns, dfiles, sfiles):
            rhist(f + '\n' + '\n'.join(map(str,output[f])), [p], args.width, args.height,
                  normcurve = args.norm_curve,
              save_data = d, save_script = s)
    if args.dot is not None:
        rdot(output, args.dot, args.width, args.height,
             args.dot_size, args.discrete_color, args.save_data, args.save_script)
    if args.box is not None:
        rbox(output, args.fields, args.stratify, args.box,args.width, args.height,
             args.outlier_size, args.color, args.save_data, args.save_script)

class PlotAssociationOpt:
    def __init__(self, master_parser):
        self.master_parser = master_parser
        subparsers = self.master_parser.add_subparsers()
        # subparser 1
        parserQQ = subparsers.add_parser('qq', help='QQ plot via ggplot2')
        self.qqArguments(parserQQ)
        self.commonArguments(parserQQ)
        # subparser 2
        parserMan = subparsers.add_parser('manhattan', help='Manhattan plot via ggplot2')
        self.manArguments(parserMan)
        self.commonArguments(parserMan)
        # subparser 3
        parserManPlain = subparsers.add_parser('manhattan_plain',
                                               help='Manhattan plot implementation not using ggplot2')
        self.manArguments(parserManPlain)
        self.commonArguments(parserManPlain)

    def get(self):
        return self.master_parser

    def qqArguments(self, parser):
        parser.add_argument('--shape',
                metavar='INTEGER',
                            type=int,
                default=1,
                help='''Choose a shape theme
                (integer 1 to 16) for dots on QQ plot.
                Default set to 1.''')
        parser.add_argument('--fixed_shape', '--fixed-shape',
                            action='store_true',
                help='''Use the same dot-shape theme for all plots''')
        parser.add_argument('--no_slope', '--no-slope',
                            action='store_true',
                help='''Do not plot the diagonal line''')

    def manArguments(self, parser):
        parser.add_argument('--chrom',
                metavar = 'CHROMOSOME',
                nargs = '+',
                default=list(map(str, list(range(1,23)))) + ['X','Y','Un'],
                help='''Specify the particular chromosome(s) to display. Can be
                one or multiple in this list: "{}". Slicing syntax "?:?" is
                supported. For example "1:22" is equivalent to displaying
                all autosomes; "1:Y" is equivalent to displaying
                all mapped chromosomes. Default set to all including unmapped
                chromosomes.'''.format(' '.join(list(map(str, list(range(1,23)))) + ['X','Y','Un', '?:?'])))
        parser.add_argument('--chrom_prefix', '--chrom-prefix',
                metavar = 'PREFIX',
                type = str,
                default = 'chr',
                help='''Prefix chromosome ID with a string.
                Default is set to "chr" (X-axis will be displayed
                as "chr1", "chr2", etc). Use "None" for no prefix.
                ''')
        parser.add_argument('--gene_map', '--gene-map',
                metavar = 'FILE',
                type = str,
                help='''If the plot units are genes and the program fails to map certain genes to
                chromosomes, you can fix it by providing a text file of genomic coordinate
                information of these genes. Each gene in the file is a line of 3 columns
                specifying "GENENAME CHROM MIDPOINT_POSITION", e.g., "IKBKB 8 42128820".
                ''')

    def commonArguments(self, parser):
        parser.add_argument("--method",
                            default = sys.argv[2] if len(sys.argv) > 2 else '',
                            help = argparse.SUPPRESS)
        settings = parser.add_argument_group('graph properties')
        settings.add_argument('-t', '--title',
                            type=str,
                default='',
                            help='''Title of plot.''')
        settings.add_argument('--color',
                            type=str,
                choices=CTHEME,
                            help='''Choose a color theme from the list above to apply
                to the plot. (via the 'RColorBrewer' package:
                cran.r-project.org/web/packages/RColorBrewer)''')
        settings.add_argument('--width_height', '--width-height',
                metavar = 'INCHES',
                nargs = 2,
                help='''The width and height of the graphics region in inches''')
        settings.add_argument('-s', '--same_page', '--same-page',
                            action='store_true',
                            help='''Plot multiple groups of p-values on the same graph''')
        settings.add_argument('-o', '--output',
                metavar = 'FILE',
                type = str,
                help='''Specify output graph filename.
                Output is in pdf format. It can be converted to jpg format
                via the 'convert' command in Linux (e.g., convert -density 180 p.pdf p.jpg)''')
        labelling = parser.add_argument_group('variants/genes highlighting')
        labelling.add_argument('-b', '--bonferroni',
                            action='store_true',
                            help='''Plot the horizontal line at 0.05/N on Y-axis
                (significance level after Bonferroni correction)''')
        labelling.add_argument('-l', '--hlines',
                metavar = 'POSITION',
                nargs = '+',
                type=float,
                help='''Additional horizontal line(s) to
                be drawn on the Y-axis.''')
        labelling.add_argument('--label_top', '--label-top',
                metavar='INTEGER',
                            type=int,
                default=1,
                help='''Specify how many top hits (smallest p-values by rank)
                you want to highlight with their identifiers in text.''')
        labelling.add_argument('--label_these', '--label-these',
                metavar='NAME',
                            type=str,
                nargs = '+',
                help='''Specify the names of variants (chr:pos, e.g., 1:87463)
                or genes (genename, e.g., IKBKB) you want to
                highlight with their identifiers in text.''')
        labelling.add_argument('-f', '--font_size', '--font-size',
                metavar='SIZE',
                            type=float,
                default=2.5,
                help='''Font size of text labels. Default set to '2.5'.''')

#
# meta analysis
#

from variant_tools.meta import MetaAnalysis
def metaAnalysisArguments(parser):
    parser.add_argument('--beta',
            metavar = 'col',
            type = int,
            default = 0,
            help='''column number of beta''')
    parser.add_argument('--pval',
            metavar = 'col',
            type = int,
            default = 0,
            help='''column number of p-value''')
    parser.add_argument('--se',
            metavar = 'col',
            type = int,
            default = 0,
            help='''column number of standard error''')
    parser.add_argument('-n', '--size',
            metavar = 'col',
            type = int,
            default = 0,
            help='''column number of sample size''')
    parser.add_argument('--link',
            metavar = 'col',
            type=int,
            nargs = '+',
            default = [0],
            help='''columns that links entries of two datasets''')
    parser.add_argument('-m',
                        '--method',
                        metavar = "method",
                        default = "ssb",
                        choices = ['ssb','ivb'],
                        help='''Method (choose from "ssb" for sample based method and "ivb" for inverse variance based method), default set to "ssb"''')
    parser.add_argument('--to_db', '--to-db',
            metavar = 'database',
            type = str,
            help='''will write the results also to a sqlite3 database compatible
            with vtools associate result format''')
    parser.add_argument('files',
            nargs = '+',
            metavar = 'file',
            help='''Input text files in the format of $vtools associate
            output (supports plain text, gz or bz2 compressed text files)''')

def metaAnalysis(args):
    fs, beta, pval, size, linker = \
            args.files, args.beta - 1, args.pval - 1, args.size - 1, [x - 1 for x in args.link]
    for x in [beta, pval, size] + linker:
        if x < 0:
            raise ValueError('Invalid column specification for "--beta/--pval/-n/--link": '
                                 'should be positive integers')
    se = None if args.se is 0 else args.se - 1
    ma = MetaAnalysis(fs, beta, pval, se, size, linker, args.method)
    # calculate p_meta and print results
    print(('\t'.join(ma.header)))
    if args.to_db:
        ma.createDB(args.to_db)
    res = []
    for grp in ma.groups:
        b, p = ma.calculate(grp)
        s = ma.sample_size[grp]
        if p > 0:
            print(('\t'.join(list(grp) + [str(b), '{:.3E}'.format(p).replace('E+00', ''), str(s)] + ['\t'.join([str(x) for x in d[1][grp] if x is not None]) for d in ma.data])))
        if args.to_db:
            res.append(list(grp))
            res[-1].extend([b, p, s])
            for d in ma.data:
                res[-1].extend([x for x in d[1][grp] if x is not None])
    if args.to_db:
        try:
            ma.writeDB(res)
        except KeyboardInterrupt:
            ma.done()
        # sys.stderr.write("Tuning database ...\n")
        ma.done()
        # sys.stderr.write("Done!\n")


def main():
    master_parser = argparse.ArgumentParser(description='''A collection of functions that
        analyze data using vtools and generate various reports''',
        prog='vtools_report',
        #formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars='@',
        epilog='''Use 'vtools_report cmd -h' for details about each command.
        Please contact Bo Peng (bpeng at mdanderson.org) if you have any question.''')
    master_parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    subparsers = master_parser.add_subparsers(title='Available reports')
    #
    # command trans_ratio
    parser = subparsers.add_parser('trans_ratio',
        help='Transition count, transversion count and transition/transversion ratio',
        description='''This command counts the number of transition (A<->G and C<->T) and
            transversion variants (others) and calculate its ratio. A ratio of 2 is expected
            from a normal sample. If option '--by_count' is specified, it will calculate
            this ratio for variants with different sample allele frequency (count). This
            commands requires a field that stores the sample count for each variant, which
            should be prepared using command 'vtools update table --from_stat "num=#(alt)"'.''')
    transRatioArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=transRatio)
    #
    # command avg_depth
    parser = subparsers.add_parser('avg_depth',
        help='Average depth for each variant, can be divided by sample variant count',
        description='''Command 'vtools update table --from_stat "meanDP=avg(DP_geno)"' calculates the average
            depth of variants across sample (e.g. average depth of three variants if the
            variant appears three times in the sample). This command report average depth
            of all variants, or variants divided by sample allele count (output count,
            number of variant, and average depth for count from 1 to 2*#sample). This
            command requires a field that stores the sample count for each variant and
            a field to store average depth of each variant, which should be prepared
            using command 'vtools update table --from_stat "num=#(alt)" "meanDP=avg(DP_geno)"'.''')
    avgDepthArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=avgDepth)
    #
    # command variant_stat
    parser = subparsers.add_parser('variant_stat',
        help='''Reports number of snps, insertions, deletions and substitutions for
            groups of samples with some size metrics to characterize the indels''',
        description='''Command 'vtools variant_stat' calculates the number of
            snps, insertions, deletions and substitutions for groups of samples
            with some size metrics to characterize the indels. The statistics can
            be calculated for all samples (effectively for the master variant table
            when parameters --samples and --group_by are ignored), a subset of samples
            (e.g. --samples aff=1), grouped by samples (e.g. --group_by aff), or for
            each sample separately (e.g. --group_by filename sample_name, because those
            two fields in the sample table uniquely identify each sample.''')
    variantStatArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=variantStat)
    #
    # command discordance_rate
    parser = subparsers.add_parser('discordance_rate',
        help='''Calculate discordance rate between pairs of samples''',
        description='''Report discordance rate, namely the number of genotype calls that differ
            between a pair of samples divided by the total number of SNPs for which both
            calls are non-missing, between pairs of samples. The statistics can be
            calculated for all samples or selected samples specified by parameter --samples.
            This command output a n by n matrix with sample names in the header. Items (i,j)
            in this matrix is numbers in the format of diff/all for i >= j, and the actual
            ratio for i < j. This rate is affected by runtime option treat_missing_as_wildtype
            which assumes that variants that do not appear in a sample (or filtered by
            quality score etc) are wildtype alleles.''')
    discordanceRateArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=discordanceRate)
    #
    # command inbreeding_coefficient
    parser = subparsers.add_parser('inbreeding_coefficient',
        help='''Calculate inbreeding coefficent F for samples using selected variants''',
        description='''Report F statistic which describe the expected degree of a reduction in
             heterozygosity when compared to Hardy-Weinberg expectation. In simple two allele
             system with inbreeding, P(AA) = p^2(1-F)+pF, P(aa) = q^2(1-F)+qF and P(HET) = 2pq(1-F).
             For an individual F is estimated by F = 1 - #observed(HET) / #expected(HET).
             Tri-allelic loci, if identified, are excluded from calculation.''')
    inbreedingCoefArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=inbreedingCoef)
    #
    # command transcript
    parser = subparsers.add_parser('transcript',
        help = '''Obtain tramscripts that covers specified region''')
    transcriptArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=transcript)
    #
    # command sequence
    parser = subparsers.add_parser('sequence',
        help = '''Obtain DNA or protein sequence in specified chromosomal region. This command by default
            outputs nucleotide sequence at the reference genome.''')
    sequenceArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=sequence)
    #
    # command transmission
    parser = subparsers.add_parser('transmission',
        help = '''Given names of parents and offspring, this command locates,
            for each offspring, de novo and recessive variants, and variants
            that demonstrate Mendelian inconsistencies. ''')
    transmissionArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=transmission)
    #
    # command plot_fields
    parser = subparsers.add_parser('plot_fields',
        help='''Draw various plots of specified variant info field(s) and/or annotation fields.
            It essentially calls 'vtools output $table $fields' and utilizes R (www.r-project.org)
            to draw the plots. The output data and R script can be saved to disk and customied.''')
    plotFieldsArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=plotFields)
    # command plot_pheno_fields
    parser = subparsers.add_parser('plot_pheno_fields',
        help='''Draw various plots of specified sample phenotype fields.
            It essentially calls 'vtools output $table $fields' and utilizes R (www.r-project.org)
            to draw the plots. The output data and R script can be saved to disk and customied.''')
    plotPhenoFieldsArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=plotPhenoFields)
    #
    # command plot_geno_fields
    parser = subparsers.add_parser('plot_geno_fields',
        help='''Draw various plots of specified genotype variant info field(s) of all or selected
            samples. It essentially calls 'vtools output $table $fields' and utilizes R (www.r-project.org)
            to draw the plots. The output data and R script can be saved to disk and customied.''')
    plotGenoFieldsArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=plotGenoFields)
    #
    # command plot_association
    parser = subparsers.add_parser('plot_association',
        help = '''Generate QQ / Manhattan plots of p-values from association analysis.
        Input data should be in the format of the output from 'vtools associate' command
        and be piped to the program as stdardard input stream (stdin)''')
    parser = PlotAssociationOpt(parser).get()
    addCommonArgs(parser)
    parser.set_defaults(func=plotAssociation)
    # command meta_analysis
    parser = subparsers.add_parser('meta_analysis',
        help = '''Meta analysis on multiple sets of association testing results.''')
    metaAnalysisArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=metaAnalysis)
    #
    # getting args
    args = master_parser.parse_args()
    env.verbosity = args.verbosity
    # calling the associated functions
    try:
        args.func(args)
    except Exception as e:
        env.logger.error('{}'.format(e))
        sys.exit(1)
