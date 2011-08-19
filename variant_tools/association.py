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


def associate(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # table?
            if not proj.isVariantTable(args.table):
                raise ValueError('Variant table {} does not exist.'.format(args.table))
            #
            # step 0: get objects that subclasses NullTest
            if not args.methods:
                raise ValueError('Please specify at least a statistical tests. Available statistical tests are {}'.format(', '.join(getAllTests())))
            tests = []
            for m in args.methods:
                m_name = m.split()[0]
                m_args = m.split()[1:] + args.unknown_args
                try:
                    method = eval(m_name)
                    if not issubclass(method, NullTest):
                        raise ValueError('Class {} is not a subclass of NullTest'.format(m_name))
                    tests.append(method(proj.logger, m_name, m_args))
                except NameError as e:
                    proj.logger.debug(e)
                    raise ValueError('Could not identify a statistical method {}. Please specify one of {}.'.format(m_name,
                        ', '.join(getAllTests())))
            #
            cur = proj.db.cursor()
            # step 1: handle --samples
            p = Sample(proj)
            if args.samples:
                IDs = p.selectSampleByPhenotype(' AND '.join(args.samples))
                if len(IDs) == 0:
                    raise ValueError('No sample is selected by condition: {}'.format(' AND '.join(args.samples)))
                else:
                    p.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
            else:
                # select all samples
                IDs = p.selectSampleByPhenotype('1')
            #
            # step 2: collect phenotype
            pheno = p.phenotypeOfSamples(args.samples, args.phenotype)
            for test in tests:
                test.setPhenotype(pheno)
            #
            # step 3: handle group_by
            groups = ['all']
            if args.group_by:
                group_fields, fields = consolidateFieldName(proj, args.table, ','.join(args.group_by))
                from_clause = args.table
                where_clause = []
                fields_info = sum([proj.linkFieldToTable(x, args.table) for x in fields], [])
                #
                processed = set()
                for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                    if (tbl.lower(), conn.lower()) not in processed:
                        from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                        where_clause.append('({})'.format(conn))
                        processed.add((tbl.lower(), conn.lower()))
                where_clause = ('WHERE ' + ' AND '.join(where_clause)) if where_clause else ''
                # select disinct fields
                query = 'SELECT {} FROM {} {};'.format(', '.join(['DISTINCT ' + x for x in fields]),
                    from_clause, where_clause)
                proj.logger.debug('Running query {}'.format(query))
                # get group by
                cur = proj.db.cursor()
                cur.execute(query)
                groups = cur.fetchall()
                proj.logger.info('Find {} groups'.format(len(groups)))
                proj.logger.debug('Group by: {}'.format(', '.join([str(x) for x in groups])))
            #
            # select variants:
            # for each group
            proj.db.attach(proj.name + '_genotype')
            for i in range(len(groups)):
                if args.group_by:
                    vtable = '__asso_tmp'
                    if proj.db.hasTable(vtable):
                        proj.db.truncateTable(vtable)
                    else:
                        proj.createVariantTable(vtable, temporary=True)
                    #
                    group_fields, fields = consolidateFieldName(proj, args.table, ','.join(args.group_by))
                    from_clause = args.table
                    where_clause = []
                    fields_info = sum([proj.linkFieldToTable(x, args.table) for x in fields], [])
                    #
                    processed = set()
                    for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                        if (tbl.lower(), conn.lower()) not in processed:
                            from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                            where_clause.append('({})'.format(conn))
                            processed.add((tbl.lower(), conn.lower()))
                    link_where = ('WHERE ' + ' AND '.join(where_clause)) if where_clause else ''
                    where_clause = (link_where + ' AND ' if link_where else 'WHERE ') +\
                        ' AND '.join(['{}={}'.format(x, proj.db.PH) for x in fields])
                    query = 'INSERT INTO __asso_tmp SELECT variant_id FROM {} {};'.format(from_clause, where_clause)
                    proj.logger.debug('Running query {}'.format(query))
                    cur.execute(query, groups[i])
                else:
                    vtable = args.table
                # 
                # collect annotation
                #  --- ignored for now
                #
                nVariant = proj.db.numOfRows(vtable)
                if nVariant == 0:
                    # this should not happen
                    proj.logger.info('No variant exists for group {}. Ignoring.'.format(', '.join(groups[i])))
                # selecting variants from each sample
                genotype = {}
                for ID in IDs:
                    query = 'SELECT variant_id, variant_type FROM {}_genotype.sample_variant_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                        .format(proj.name, ID, vtable)
                    proj.logger.debug('Running query {}'.format(query))
                    cur.execute(query)
                    genotype[ID] = cur.fetchall()
                    # proj.logger.info(str(ID) + str(genotype[ID]))
                #
                # passing everything to association test
                for test in tests:
                    #
                    test.setGenotype(genotype)
                    #
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
        self.phenotype = None
        self.genotype = None
        self.annotation = None
        self.logger = logger
        self.name = name
        self.parseArgs(*method_args)
        self.t = c_level_class(self.phenotype, self.genotype)

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
        self.t.setPhenotype(data)

    def setGenotype(self, data):
        self.genotype = data

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant'''
        self.t.permute()
        return t.calculate()
        self.logger.info('Phenotype: {}'.format(len(self.phenotype)))
        self.logger.info('Genotype: {}'.format(len(self.genotype)))
        return 1

class ExternTest(NullTest):
    '''A test that exports data in standard formats, call an external program
    and prase its output. This is the simplest, but also the slowest method
    to utilize a third-party association test program.
    '''
    def __init__(self):
        pass

class WssTest(NullTest):
    def __init__(self, phenotypes, genotypes, covariates, annotations, mode, weightingTheme):
        self.y = phenotypes
        self.x = genotypes
        self.z = covariates
        self.a = annotations
        self.wt = weightingTheme
        if mode == 'R':
            scipy.where(x==2.0, 1.0, 0.0)
        elif mode == 'D':
            scipy.where(x!=0.0, 1.0, 0.0)
        else:
            pass
        
    def weighting(self):
        nsamples = [self.y.shape[0]]*self.y.shape[1]
        ncases = self.y.sum(axis=0)
        nctrls = [s - cs for s, cs in zip(nsamples, ncases)]
        nvariants = self.x.sum(axis=1)
        #!NOTICE: Now focus only on the primary phenotype, y[:,0]
        countsinctrls = [sum([m[i] for i in range(nsamples[0]) if not self.y[i, 0]]) for m in self.x.transpose()]
        countsincases = [all - ct for all, ct in zip(nvariants, countsinctrls)]
        # apply weighting theme based on annotation
        if self.wt == 1:
            weights = -scipy.log([1-i for i in self.a])
        #
        # apply weighting theme based on observed sample
        # there can be many more other weighting themes ... or combinations of these themes
        else:
            # weigthing by variants exclusive to ctrls
            q = [(m+1.0)/(nctrls[0]+2.0) for m in countsinctrls]
            weights = [1.0/scipy.sqrt(freq*(1.0-freq)) for freq in q]
        return weights, countsincases
    
    def simplePerm(self, nperm):
        weights = self.weighting()
        statistic = sum([w*c for w,c in zip(weights[0], weights[1])])
        permcount1, permcount2 = 0, 0
        for i in range(nperm):
            numpy.random.shuffle(self.y[:,0])
            weights = self.weighting()
            perm_statistic = sum([w*c for w,c in zip(weights[0], weights[1])])
            if perm_statistic >= statistic: permcount1 += 1
            if perm_statistic <= statistic: permcount2 += 1
        return permcount1, permcount2
    
    def regressionR(self, covariate=[1,2]):
        # logistic regression using R via rpy2
        # dirty codes for proof of concept
        weights = scipy.array(self.weighting()[0])
        regressors = scipy.matrix(self.x) * scipy.matrix(weights).transpose()
        #   input = scipy.concatenate((y, regressors, z), axis=1)
        R.assign('Y', Rfloats(self.y[:,0].tolist()))
        R.assign('X', Rfloats(regressors.transpose().tolist()[0]))
        for i in covariate:
            R.assign('Z'+str(i), Rfloats(self.z[:, i-1].tolist()))
        summary = R('glm(Y~X+Z1+Z2, family = "binomial")')
        beta_x = summary[0][1]
        return beta_x
    
    def regressionPerm(self, nperm, covariate=[1,2]):
        statistic = self.regressionR(covariate)
        permcount1, permcount2 = 0, 0
        for i in range(nperm):
            # permutation on X, permute the rows
            numpy.random.shuffle(self.x)
            perm_statistic = self.regressionR(covariate)
            if perm_statistic >= statistic: permcount1 += 1
            if perm_statistic <= statistic: permcount2 += 1
        return permcount1, permcount2
    
    def calcP1(self, covariate, nperm):
        # one-sided test
        if covariate:
            counts = self.regressionPerm(nperm, covariate)[0]
            return float(counts)/float(nperm)
        else:
            counts = self.simplePerm(nperm)[0]
            return float(counts)/float(nperm)        
        
    def calcP2(self, covariate, nperm):
        # two-sided test
        if covariate:
            counts = min(self.regressionPerm(nperm, covariate))
            return 2.0 * float(counts)/float(nperm)
        else:
            counts = min(self.simplePerm(nperm))
            return 2.0 * float(counts)/float(nperm)

