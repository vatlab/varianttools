#!/usr/bin/env python
#
# $File: association.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
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
from array import array
from .project import Project
from .utils import ProgressBar, consolidateFieldName
from .sample import Sample
import argparse


def associateArguments(parser):
    parser.add_argument('table', help='''Variant table.''')
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
            is equivalent to --method "m1 -p 10" "m2 -p 10".). Available statistical
            tests are {}.'''.format(', '.join(getAllTests())))
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-g', '--group_by', nargs='*',
        help='''Group variants by fields. If specified, variants will be separated
            into groups and are tested one by one.''')


class AssociationTester(Sample):
    def __init__(self, proj, table):
        Sample.__init__(self, proj)
        # table?
        if not self.proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        self.table = table

    def getAssoTests(self, methods, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common 
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least a statistical tests. Available statistical tests are {}'.format(', '.join(getAllTests())))
        self.tests = []
        for m in methods:
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                method = eval(name)
                self.tests.append(method(self.logger, name, args))
            except NameError as e:
                self.logger.debug(e)
                raise ValueError('Could not identify a statistical method {}. Please specify one of {}.'.format(name,
                    ', '.join(getAllTests())))

    def getSamples(self, condition):
        '''Get a list of samples from specified condition. This function sets self.IDs'''
        if condition:
            self.IDs = self.proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in condition]))
            if len(self.IDs) == 0:
                raise ValueError('No sample is selected by condition: {}'.format(' AND '.join(['({})'.format(x) for x in condition])))
            else:
                self.logger.info('{} condition are selected by condition: {}'.format(len(self.IDs), ' AND '.join(['({})'.format(x) for x in condition])))
        else:
            # select all condition
            self.IDs = self.proj.selectSampleByPhenotype('1')

    def getPhenotype(self, condition, phenotype):
        '''Get phenotype for specified samples (specified by condition).'''
        try:
            query = 'SELECT sample_id, {} FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id'.format(', '.join(phenotype)) + \
                (' WHERE {}'.format(' AND '.join(['({})'.format(x) for x in condition])) if condition else '') + ';'
            self.logger.debug('Select phenotype using query {}'.format(query))
            cur = self.db.cursor()
            cur.execute(query)
            self.phenotype = array('d', zip(*cur.fetchall())[1])
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to retrieve phenotype '.format(', '.join(phenotype)))

    def identifyGroups(self, group_by):
        '''Get a list of groups according to group_by fields'''
        self.group_by = group_by
        if not group_by:
            self.groups = ['all']
        else:
            group_fields, fields = consolidateFieldName(self.proj, self.table, ','.join(group_by))
            self.from_clause = self.table
            where_clause = []
            fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
            #
            processed = set()
            for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                if (tbl.lower(), conn.lower()) not in processed:
                    self.from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                    where_clause.append('({})'.format(conn))
                    processed.add((tbl.lower(), conn.lower()))
            #
            self.where_clause = ('WHERE ' + ' AND '.join(where_clause)) if where_clause else ''
            # select disinct fields
            query = 'SELECT {} FROM {} {};'.format(', '.join(['DISTINCT ' + x for x in fields]),
                self.from_clause, self.where_clause)
            self.logger.debug('Running query {}'.format(query))
            # get group by
            cur = self.db.cursor()
            cur.execute(query)
            self.groups = cur.fetchall()
            self.logger.info('Find {} groups'.format(len(self.groups)))
            self.logger.debug('Group by: {}'.format(', '.join([str(x) for x in self.groups])))

    def getVariants(self, group):
        '''Get variants for a certain group'''
        if not self.group_by:
            # group must be 'all'
            return self.table
        group_fields, fields = consolidateFieldName(self.proj, self.table, ','.join(self.group_by))
        vtable = '__asso_tmp'
        if self.db.hasTable(vtable):
            self.db.truncateTable(vtable)
        else:
            self.proj.createVariantTable(vtable, temporary=True)
        #
        where_clause = (self.where_clause if self.where_clause else 'WHERE ') + ' AND '.join(['{}={}'.format(x, self.proj.db.PH) for x in fields])
        query = 'INSERT INTO __asso_tmp SELECT variant_id FROM {} {};'.format(self.from_clause, where_clause)
        self.logger.debug('Running query {}'.format(query))
        cur = self.db.cursor()
        cur.execute(query, group)
        return vtable

    def getGenotype(self, vtable):
        '''Get genotype for variants in specified table'''
        genotype = []
        cur = self.db.cursor()
        cur.execute('SELECT min(variant_id), max(variant_id) FROM {}'.format(vtable))
        start, end = cur.fetchall()[0]
        for ID in self.IDs:
            query = 'SELECT variant_id, GT FROM {}_genotype.genotype_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                .format(self.proj.name, ID, vtable)
            self.logger.debug('Running query {}'.format(query))
            cur.execute(query)
            gtmp = {x[0]:x[1] for x in cur.fetchall()}
            genotype.append(array('d', [gtmp.get(x, -9) for x in range(start,end+1)]))
        return genotype

def associate(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            asso = AssociationTester(proj, args.table)
            # step 0: get testers
            asso.getAssoTests(args.methods, args.unknown_args)
            # step 1: get samples
            asso.getSamples(args.samples)
            # step 2: get phenotype and set it to everyone
            asso.getPhenotype(args.samples, args.phenotype)
            for test in asso.tests:
                test.setPhenotype(asso.phenotype)
            # step 3: handle group_by
            asso.identifyGroups(args.group_by)
            proj.db.attach(proj.name + '_genotype')
            for grp in asso.groups: 
                # select variants from each group:
                vtable = asso.getVariants(grp)
                # passing everything to association test
                genotype = asso.getGenotype(vtable)
                for test in asso.tests:
                    test.setGenotype(genotype)
                    # step 6: call stat.calculate
                    values = test.calculate()
                    # step 7: update variant table
                    fields = test.getOption('fields')
                    # ...
    except Exception as e:
        sys.exit(e) 

#
#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
def getAllTests():
    '''List all tests (all classes that subclasses of NullTest) in this module'''
    return [name for name,obj in globals().iteritems() if type(obj) == type(NullTest) and issubclass(obj, NullTest)]

import variant_tools.assoTests as t

class NullTest:
    '''A base class that defines a common interface for association
    statistics calculator. Instances of these calculators will be created
    by vtools during command 'vtools asscoaition'. The results will be
    automatically inserted back to related variant tables. An association
    statistics calculator will have several data structures that stores
        1. phenotype
        2. genotype
        3. annotation
    These structures should not be changed during the calculation.
    '''
    def __init__(self, logger=None, name=None, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to 
        parse it'''
        self.data = t.AssoData()
        self.logger = logger
        self.name = name
        self.parseArgs(*method_args)
        self.result = {}

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''A base association test method
            test does nothing, but can be used to measure the performance of 
            retrieving data from vtools.''',
            prog='vtools associate --method ' + self.name)
        # no argumant is added
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
    
    def getOption(self, name):
        '''Get option for the association test.'''
        if name == 'fields':
            # the null test does not set any field
            return []

    def setPhenotype(self, data):
        '''Set phenotype data'''
        self.data.setPhenotype(data)

    def setGenotype(self, data):
        self.data.setGenotype(data)
        self.data.setMaf()

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant'''
        self.logger.info('#samples: {}'.format(len(self.data.phenotype())))
        self.logger.info('#variants: {}'.format(len(self.data.raw_genotype()[0])))
        return 1

class ExternTest(NullTest):
    '''A test that exports data in standard formats, call an external program
    and prase its output. This is the simplest, but also the slowest method
    to utilize a third-party association test program.
    '''
    def __init__(self, logger=None, name=None, *method_args):
        pass

class LogisticBurdenTest(NullTest):
    '''Simple Logistic regression score test on collapsed genotypes
    within an association testing group
    '''
    def __init__(self, logger=None, name=None, *method_args):
        NullTest.__init__(self, logger, name, *method_args)
        self.ptime = 1000
        self.pvalue = None;
   
    def calculate(self):
        data = self.data.clone()
        a = t.ActionExecuter([t.SumToX(), t.SimpleLogisticRegression(), t.GaussianPval(1)])
        a.apply(data)
        self.pvalue = data.pvalue()
        self.logger.info('p-value (asymptotic) = {}'.format(self.pvalue))
        # logistic regression, permutation 
        p = t.PhenoPermutator(self.ptime, [t.SimpleLogisticRegression()])
        self.logger.info('p-value (permutation) = {}'.format(p.apply(data) / float(self.ptime)))
        return 1