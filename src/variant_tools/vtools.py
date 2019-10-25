#!/usr/bin/env python2.7
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

import argparse
import subprocess
import sys

import variant_tools.annotation as annotation
import variant_tools.compare as compare
import variant_tools.exporter as exporter
import variant_tools.importer as importer
import variant_tools.phenotype as phenotype
import variant_tools.pipeline as pipeline
import variant_tools.project as project
import variant_tools.update as update
import variant_tools.variant as variant
from variant_tools._version import VTOOLS_FULL_VERSION
from variant_tools.utils import env

if sys.platform != 'win32':
    import variant_tools.liftOver as liftOver
    import variant_tools.association as association

# save the command line that has been processed by the shell.
env.command_line = subprocess.list2cmdline(sys.argv[1:])

def addCommonArgs(parser):
    parser.add_argument(
        '-v',
        '--verbosity',
        choices=['0', '1', '2', '3'],
        help='''Output error and warning (0), info (1), debug (2) and trace (3)
            information to standard output (default to 1).'''),


def main():
    #
    master_parser = argparse.ArgumentParser(
        description='''A variant calling,
        processing, annotation and analysis tool for next-generation sequencing
        studies.''',
        prog='vtools',
        # formatter_class=argparse.RawDescriptionHelpFormatter,
        fromfile_prefix_chars='@',
        epilog='''Use 'vtools cmd -h' for details about each command.
        Please contact Bo Peng (bpeng at mdanderson.org) if you have any question.'''
    )
    master_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(VTOOLS_FULL_VERSION))
    subparsers = master_parser.add_subparsers(title='subcommands')
    #
    # command init
    parser = subparsers.add_parser(
        'init',
        help='''Create a new project, or a
        subproject from an existing parent project, or merge several existing
        projects into one''',
        description='''Create a new project in the current directory. This
        command will fail if another project already exists in this directory,
        unless option '--force' is used to remove the existing project.''')
    project.initArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=project.init)
    #
    # command import
    parser = subparsers.add_parser(
        'import',
        help='Import variants and related sample genotype from files in specified formats',
        description='''Import variants and related sample genotype from one or more
            delimiter separated files (e.g. VCF and a number of indel formats).'''
    )
    importer.importVariantsArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=importer.importVariants)
    #
    # command phenotype
    parser = subparsers.add_parser(
        'phenotype',
        help='Manage sample phenotypes',
        description='''Import phenotypes from a file, or set phenotypes to constants,
            or to summary statistics of sample genotype fields.''')
    phenotype.phenotypeArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=phenotype.phenotype)
    #
    # command show
    parser = subparsers.add_parser(
        'show',
        help='Display content of a project',
        description='''Output information of all system and project related items
            such as variant tables, samples, phenotypes, annotation databases
            and fields.''')
    project.showArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=project.show)
    #
    # command liftover
    if sys.platform != 'win32':
        parser = subparsers.add_parser(
            'liftover',
            help='''Set alternative reference genome and update alternative
                coordinates of all variant tables''',
            description='''Convert coordinates of existing variants to alternative
                coordinates in an alternative reference genome. The UCSC liftover
                tool will be automatically downloaded if it is not available.'''
        )
        liftOver.liftOverArguments(parser)
        addCommonArgs(parser)
        parser.set_defaults(func=liftOver.liftOver)
    #
    # command use
    parser = subparsers.add_parser(
        'use',
        help='Prepare (download or import if necessary) and use an annotation database',
        description='''Link an annotation database to the project, download it
            from the variant tools website or build it from source if needed.'''
    )
    annotation.useArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=annotation.use)
    #
    # command update
    parser = subparsers.add_parser(
        'update',
        help='''Add or update fields of existing variants and genotype using
            information from specified existing fields, sample genotype,
            or external files''',
        description='''Add or update fields of existing variants and genotype
            from other fields, statistics of genotypes and genotype info, or
            files that annotate variants or their locations (e.g. Read
            annotation from ANNOVAR output files, import additional variant
            or genotype fields from .vcf files).''')
    update.updateArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=update.update)
    #
    # command select
    parser = subparsers.add_parser(
        'select',
        help='''Output or save select variants that match specified conditions''',
        description='''Select variants according to properties (variant and
            annotation fields) and membership (samples) of variant. The result
            can be counted, outputted, or saved to a variant table.''')
    variant.selectArguments(parser)
    variant.generalOutputArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=variant.select)
    #
    # command exclude
    parser = subparsers.add_parser(
        'exclude',
        help='''Output or save variants after excluding variants that match
            specified conditions''',
        description='''Exclude variants according to properties (variant and
            annotation fields) and membership (samples) of variant. The result
            can be counted, outputted, or saved to a variant table.''')
    variant.excludeArguments(parser)
    variant.generalOutputArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=variant.exclude)
    #
    # command compare
    parser = subparsers.add_parser(
        'compare',
        help='''Compare sites, variants, or genotypes of variants in two or more
            variant tables''',
        description='''Get the difference, intersection and union of two or more
            variant tables, according to sites, variants, or genotypes of
            associated samples of these variants. Resulting variants can be
            counted or write to other variant tables.''')
    compare.compareArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=compare.compare)
    #
    # command output
    parser = subparsers.add_parser(
        'output',
        help='Output variants in tab or comma separated format',
        description='''Output variants, variant info fields, annotation fields
            and expressions that involve these fields in a tab or comma separated
            format.''')
    variant.outputArguments(parser)
    variant.generalOutputArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=variant.output)
    #
    # command export
    parser = subparsers.add_parser(
        'export',
        help='Export samples (variants and genotypes) in specified format',
        description='''Export variants and genotypes in text, vcf and other
            formats.''')
    exporter.exportArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=exporter.export)
    #
    # command remove
    parser = subparsers.add_parser(
        'remove',
        help='''Remove project or its
        contents such as variant tables, fields, and annotation databases.''',
        description='''Remove from the current project various items such as
            variants genotypes, and annotation fields.''')
    project.removeArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=project.remove)
    #
    # command associate
    #
    if sys.platform != 'win32':
        parser = subparsers.add_parser(
            'associate',
            help='''Test association between variants and phenotypes''',
            description='''Call one or more statistical association tests and
                return test results as fields to variants tested.''')
        association.associateArguments(parser)
        addCommonArgs(parser)
        parser.set_defaults(func=association.associate)
    #
    # command admin
    #
    parser = subparsers.add_parser(
        'admin',
        help='''Perform various administrative tasks including merge and
            rename samples.''',
        description='''Optimize or modify projects. Currently supports merging
            and rename of samples''')
    project.adminArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=project.admin)
    #
    # command execute
    parser = subparsers.add_parser(
        'execute',
        help='Execute a SQL query',
        description='''Execute a SQL query or pipeline.''')
    pipeline.executeArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=pipeline.execute)
    #
    # command simulate
    #
    parser = subparsers.add_parser(
        'simulate',
        help='''Simulate sequencing data using specified simulation
            models.''',
        description='''Simulate case control or family-based samples using
            specified simulation models.''')
    pipeline.simulateArguments(parser)
    addCommonArgs(parser)
    parser.set_defaults(func=pipeline.simulate)
    #
    # getting args, some commands accept arbitrary arguments so we need to
    # separate them into argv
    #
    # http://bugs.python.org/issue12353
    # parse_known_args cannot handle arguments with empty string '', we need
    # to manually change it to ' ' to get around of this problem.
    # Hopefully this will not cause any trouble in later argument processing.
    if '' in sys.argv:
        for idx, arg in enumerate(sys.argv):
            if arg == '':
                sys.argv[idx] = ' '
    #
    args, argv = master_parser.parse_known_args()
    # python3 does not raise an error "too few arguments" when no argument is
    # provided so we have to explicitly handle this case. The output of command
    # vtools is therefore different in python 2 (error message) and 3 (help message).
    if not hasattr(args, 'func'):
        master_parser.print_help()
        sys.exit(0)
    # commands that accepts format and other subparsers can pass arbitrary
    # parameter to subparsers
    if args.func in [
            association.associate, importer.importVariants, update.update,
            exporter.export, pipeline.execute, pipeline.simulate
    ]:
        args.unknown_args = argv
    elif len(argv) > 0:
        master_parser.print_usage(sys.stderr)
        sys.exit('vtools: error: unrecognized arguments: ' + ' '.join(argv))
    # calling the associated functions
    args.func(args)
