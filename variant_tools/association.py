#!/usr/bin/env python
#
# $File: association.py $
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
import argparse


class BaseAssociationStat:
    '''A base class that defines a common interface for association
    statistics calculator. Instances of these calculators will be created
    by vtools during command 'vtools asscoaition'. The results will be
    automatically inserted back to related variant tables. A association
    statistics calculator will have several data structures that stores
        1. phenotype
        2. genotype
        3. annotation
    These structures should not be changed during the calculation.
    '''
    def __init__(self, logger=None, name=None, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to 
        parse it'''
        self.phenotype = None
        self.genotype = None
        self.annotation = None
        self.logger = logger
        self.name = name
        self.parseArgs(*method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''A base association test method
            test does nothing, but can be used to measure the performance of 
            retrieving data from vtools.''',
            prog='vtools associate --method ' + self.name)
        # no argumant is added
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
            
    def setData(self, phenotype=None, genotype=None, annotation=None):
        '''Set or update internal data structure. These variables are set
        separately because it is likely that only some of them needs update.'''
        if self.phenotype is not None:
            self.phenotype = phenotype
        if self.genotype is not None:
            self.genotype = genotype
        if annotation is not None:
            self.annotation = annotation

    def caclulate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant'''
        return None
        

def associateArguments(parser):
    parser.add_argument('variants', help='''Variant table.''')
    parser.add_argument('phenotype', nargs='+',
        help='''A list of phenotypes that will be passed to the association
            statistics calculator''')
    parser.add_argument('--covariants', nargs='*',
        help='''Optional phenotypes that will be passed to statistical
            tests as covariants. Values of these phenotypes should be integer
            or float.''')
    parser.add_argument('-m', '--methods', nargs='+',
        help='''Method of one or more association tests. Parameters for each
            method should be specified together as a quoted long argument (e.g.
            --method "m --field mp" "m1 --permute 10000 --field m1p"), although
            the common method parameters can be specified separately, as long as
            they do not conflict with command arguments. (e.g. --method m1 m2 -p 10 
            is equivalent to --method "m1 -p 10" "m2 -p 10".)''')
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-g', '--group_by', nargs='*',
        help='''Group variants by fields. If specified, variants will be separated
            into groups and are tested one by one.''')


def associate(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # table?
            if not proj.isVariantTable(args.variants):
                raise ValueError('Variant table {} does not exist.'.format(args.variants))
            #
            # step 0: find a subclass of the BaseAssociationStat
            for m in args.methods:
                m_name = m.split()[0]
                m_args = m.split()[1:] + args.unknown_args
                stat = BaseAssociationStat(proj.logger, m_name, m_args)
            # step 1: handle --samples
            #
            # step 2: collect phenotype
            stat.setData(phenotype=None)
            #
            # step 3: handle group_by --- ignored for now
            # for each group
            for i in range(1):
                # prog = ProgressBar(i, 100)
            #
            #     step 4: collect genotype
                stat.setData(genotype=None) 
            #
            #     step 5: collect annotation
                stat.setData(annotation=None)
            #
            #     step 6: call stat.calculate
            #     NOTE: need to define a callback function
                values = stat.calculate()
            #
            #     step 7: update variant table
            #
    except Exception as e:
        sys.exit(e) 

