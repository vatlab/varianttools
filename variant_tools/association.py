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

import sys
import os
import threading
from multiprocessing import Process, Queue, Pipe, Lock
import time
from array import array
from copy import copy, deepcopy
from .project import Project
from .utils import ProgressBar, consolidateFieldName, DatabaseEngine
from .phenotype import Sample
import argparse
import variant_tools.assoTests as t

def associateArguments(parser):
    parser.add_argument('table', help='''Variant table.''')
    parser.add_argument('phenotype', nargs='+',
        help='''A list of phenotypes that will be passed to the association
            statistics calculator''')
    parser.add_argument('--covariates', nargs='*',
        help='''Optional phenotypes that will be passed to statistical
            tests as covariates. Values of these phenotypes should be integer
            or float.''')
    parser.add_argument('-m', '--methods', nargs='+',
        help='''Method of one or more association tests. Parameters for each
            method should be specified together as a quoted long argument (e.g.
            --method "m --alternative 2" "m1 --permute 1000"), although
            the common method parameters can be specified separately, as long as
            they do not conflict with command arguments. (e.g. --method m1 m2 -p 1000 
            is equivalent to --method "m1 -p 1000" "m2 -p 1000".). Available statistical
            tests are {}.'''.format(', '.join(getAllTests())))
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-g', '--group_by', nargs='*',
        help='''Group variants by fields. If specified, variants will be separated
            into groups and are tested one by one.''')
    parser.add_argument('-j', '--jobs', metavar='N', default=1, type=int,
        help='''Number of processes to carry out association tests.''')
    parser.add_argument('--to_file', metavar='FILE', nargs=1,
        help='''File name to which results from association tests will be written''')        


class AssociationTester(Sample):
    '''Parse command line and get data for association testing'''
    
    def __init__(self, proj, table):
        Sample.__init__(self, proj)
        # table?
        if not self.proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        self.table = table
        self.phenotype = None
        self.covariates = None


    def getAssoTests(self, methods, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common 
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least one statistical test. Available statistical tests are {}'.format(', '.join(getAllTests())))
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
            self.phenotype = [array('d', map(float, x)) for x in list(zip(*cur.fetchall()))[1:]]
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('''Failed to retrieve phenotype {}. Please 
                             make sure the specified phenotype names are correct and there 
                             is no missing value'''.format(', '.join(phenotype)))
    
    def getCovariate(self, condition, covariates):
        '''Get covariates for specified samples (specified by condition).'''
        if not covariates:
            return
        try:
            query = 'SELECT sample_id, {} FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id'.format(', '.join(covariates)) + \
                (' WHERE {}'.format(' AND '.join(['({})'.format(x) for x in condition])) if condition else '') + ';'
            self.logger.debug('Select phenotype covariates using query {}'.format(query))
            cur = self.db.cursor()
            cur.execute(query)
            self.covariates = [array('d', map(float, x)) for x in list(zip(*cur.fetchall()))[1:]]
            self.covariates.insert(0, array('d', [1]*len(self.covariates[0])))
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('''Failed to retrieve phenotype covariates {}. Please 
                             make sure the specified phenotype covariates names are correct and 
                             there is no missing value'''.format(', '.join(covariates)))
    
    def identifyGroups(self, group_by):
        '''Get a list of groups according to group_by fields'''
        # set default group_by to positions
        if not group_by:
          self.group_by = ['chr','pos']
        else:
          self.group_by = group_by
        #
        table_of_fields = [self.proj.linkFieldToTable(field, self.table)[-1].table for field in self.group_by]
        table_of_fields = [x if x else self.table for x in table_of_fields] 
        fields_names = [x.replace('.', '_') for x in self.group_by]
        fields_types = [self.db.typeOfColumn(y, x.rsplit('.', 1)[-1]) for x,y in zip(self.group_by, table_of_fields)]
        cur = self.db.cursor()
        # create a table that holds variant ids and groups, indexed for groups
        cur.execute('DROP TABLE IF EXISTS __asso_tmp;')
        cur.execute('DROP INDEX IF EXISTS __asso_tmp_index;')
        cur.execute('''\
            CREATE TABLE __asso_tmp (
              variant_id INT NOT NULL,
              {});
              '''.format(','.join(['{} {} NULL'.format(x,y) for x,y in zip(fields_names, fields_types)])))
        # select variant_id and groups for association testing
        group_fields, fields = consolidateFieldName(self.proj, self.table, ','.join(self.group_by))        
        from_clause = []
        from_clause.append(self.table)
        where_clause = []
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause.append('{}'.format(tbl))
                where_clause.append('({})'.format(conn))
                processed.add((tbl.lower(), conn.lower()))
        #
        self.from_clause = ', '.join(from_clause) 
        self.where_clause = ('WHERE ' + ' AND '.join(where_clause)) if where_clause else ''
        # This will be the tmp table to extract variant_id by groups
        query = 'INSERT INTO __asso_tmp SELECT DISTINCT {}.variant_id, {} FROM {} {};'.format(self.table, group_fields,
            self.from_clause, self.where_clause) 
        self.logger.debug('Running query {}'.format(query))
        cur.execute(query)
        cur.execute('''\
            CREATE INDEX __asso_tmp_index ON __asso_tmp ({});
            '''.format(','.join(['{} ASC'.format(x) for x in fields_names])))
        # get group by
        cur.execute('''\
            SELECT DISTINCT {} FROM __asso_tmp;
            '''.format(', '.join(fields_names)))
        self.groups = cur.fetchall() 
        self.logger.info('Find {} groups'.format(len(self.groups)))
        self.logger.debug('Group by: {}'.format(', '.join(map(str, self.groups))))


class UpdateResult:
    def __init__(self, db, logger, fn):
        self.grps = []
        self.lock = Lock()
        self.db = db
        self.logger = logger
        self.fn = fn[0] if fn else None

    def set(self, table, res):
        self.grps.append(res[0])

        if not self.fn:
            self.fn = 'vtoolsasso_{}.result'.format(time.strftime('%b%d_%H%M%S', time.gmtime()))
        try:
            #output = '[{0}]\ntest = {1}\np-value = {2}\nstatistic = {3}\nsamples = {4}\n'.format('__'.join(map(str, res[0])),
            #            res[1][0][0][0][0].split('_')[0], res[1][0][1]['pvalue'], res[1][0][1]['statistic'], res[1][0][1]['samples'])
            col_grp = '__'.join(map(str, res[0]))
            col_test = res[1][0][0][0][0].split('_')[0]
            cols_stat = res[1][0][1]
            output = '\t'.join(map(str, [col_grp, col_test, cols_stat['pvalue'], cols_stat['statistic'], cols_stat['samples']])) + '\n'
            if len(self.grps) == 1:
                output = '\t'.join(['#test_unit', 'method', 'pvalue', 'statistic', 'sample_size']) + '\n' + output
            self.lock.acquire()
            with open(self.fn, 'a') as out:
                out.write(output)
            self.lock.release()
        except IndexError as e:
            self.logger.warning('No association tests done for {}'.format(self.grps[-1]))
            pass
        return
        
    def getgrp(self):
        return self.grps
    
    def count(self):
        return len(self.grps)
        
class GroupAssociationCalculator(Process):
    '''Association test calculator'''
    def __init__(self, proj, table, samples, phenotypes, covariates, tests, group_by, grpQueue, output):
        self.proj = proj
        self.table = table
        self.IDs = samples
        self.phenotypes = phenotypes
        self.covariates = covariates
        self.tests = tests
        self.group_by = group_by
        self.queue = grpQueue
        self.output = output
        self.logger = proj.logger
        self.db = None
        Process.__init__(self, name='Phenotype association analysis for a group of variants')
    

    def getGenotype(self, group):
        '''Get genotype for variants in specified group'''
        #print 'Getting', group, self.table, self.group_by
        fields_names = [x.replace('.', '_') for x in self.group_by]
        where_clause =  'WHERE ' + \
            ' AND '.join(['{0}={1}'.format(x, self.db.PH) for x in fields_names])
        #
        # get variant_id
        query = 'SELECT variant_id FROM __asso_tmp {}'.format(where_clause)
        self.logger.debug('Running on group {0} query {1}'.format(group, query))
        cur = self.db.cursor()
        cur.execute(query, group)
        variant_id = [x[0] for x in cur.fetchall()]
        numSites = len(variant_id)
        self.logger.debug('Getting {0} variants for {1}'.format(numSites, group))
        #
        # get genotypes
#        cur.execute('DROP TABLE IF EXISTS __id_of_group;')
#        cur.execute('CREATE TEMPORARY TABLE __id_of_group (variant_id INT NOT NULL);')
#        cur.execute('INSERT INTO __id_of_group SELECT variant_id FROM __asso_tmp {};'.format(where_clause), group)
        genotype = []
        self.db.attach(self.proj.name+'_genotype.DB', '__fromGeno')
        for ID in self.IDs:
#            query = 'SELECT variant_id, GT FROM __fromGeno.genotype_{} WHERE variant_id IN (SELECT variant_id FROM __id_of_group);'\
            query = 'SELECT variant_id, GT FROM __fromGeno.genotype_{0} WHERE variant_id IN (SELECT variant_id FROM __asso_tmp {1});'\
                .format(ID, where_clause)
            cur.execute(query, group)
#            cur.execute(query)
            gtmp = {x[0]:x[1] for x in cur.fetchall()}
            # handle missing values
            gtmp = [gtmp.get(x, -9.0) for x in variant_id]
            # handle -1 coding (double heterozygotes)
            gtmp = [2.0 if int(x) == -1 else x for x in gtmp]
            genotype.append(array('d', gtmp))
        #
        self.logger.debug('Retrieved genotypes for {} samples'.format(len(genotype)))
        self.db.detach('__fromGeno')
        missing_counts = [x.count(-9.0) for x in genotype]
        # remove individuals having many missing genotypes, or have all missing variants
        # FIXME will pass it as an input arguement later
        #toKeep = [(x<0.5*numSites) for x in missing_counts]
        toKeep = [(x<numSites) for x in missing_counts]
        self.logger.debug('{} samples will be removed due to missing genotypes'.format(len(self.IDs)-sum(toKeep)))
        return genotype, toKeep
    

    def run(self):
        self.db = DatabaseEngine()
        self.db.connect(self.proj.name+'.proj', readonly=True)

        #
        while True:
            self.logger.debug('Getting group ...') 
            grp = self.queue.get()
            self.logger.debug('Got group {}'.format(grp))
            if grp is None:
                self.logger.debug('All groups are processed, terminating the queue')
                self.output.send(None)
                break
            # select variants from each group:
            genotype, which = self.getGenotype(grp)
            values = []
            try:
                for test in self.tests:
                    test.setPhenotype(which, self.phenotypes, self.covariates)
                    test.setGenotype(which, genotype)
                    test.setAttributes(grp)
                    test.calculate()
                    self.logger.debug('Finish test')
                    values.append([test.getFields(), test.result])
            except Exception as e:
                self.logger.info('Error processing group {}, {}'.format(grp, e))
            self.logger.debug('Finished group {}'.format(grp))
            self.output.send((grp, values))
        

def associate(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            asso = AssociationTester(proj, args.table)
            # step 0: get testers
            asso.getAssoTests(args.methods, args.unknown_args)
            # step 1: get samples
            asso.getSamples(args.samples)
            #
            proj.db.attach('{}_genotype.DB'.format(proj.name), '__fromGeno')
            unindexed_IDs = [x for x in asso.IDs if \
                not proj.db.hasIndex('__fromGeno.genotype_{}_index'.format(x))]
            if unindexed_IDs:
                cur = proj.db.cursor()
                prog = ProgressBar('Indexing genotypes', len(unindexed_IDs))
                for idx, ID in enumerate(unindexed_IDs):
                    cur.execute('CREATE INDEX __fromGeno.genotype_{0}_index ON genotype_{0} (variant_id ASC)'.format(ID))
                    prog.update(idx + 1)
                prog.done()
            # step 2: get phenotype and set it to everyone
            asso.getPhenotype(args.samples, args.phenotype)           
            asso.getCovariate(args.samples, args.covariates)
            # step 3: handle group_by
            asso.identifyGroups(args.group_by)
            nJobs = max(min(args.jobs, len(asso.groups)), 1)
            # step 4: start all workers
            grpQueue = Queue()
            results = UpdateResult(proj.db, proj.logger, args.to_file)
            readers = []
            for j in range(nJobs):
                r, w = Pipe(False)
                GroupAssociationCalculator(proj, args.table, asso.IDs, asso.phenotype[0], asso.covariates,
                    asso.tests, asso.group_by, grpQueue, w).start()
                readers.append(r)
            # put all jobs to queue, the workers will work on them
            for grp in asso.groups:
                grpQueue.put(grp)
            for j in range(nJobs):
                grpQueue.put(None)
            count = 0
            prog = ProgressBar('Testing for association', len(asso.groups))
            #
            proc_status = [True] * len(readers)
            while True:
                for idx, (s,r) in enumerate(zip(proc_status, readers)):
                    if not s:
                        continue
                    res = r.recv()
                    if res is None:
                        proc_status[idx] = False
                    else:
                        results.set(args.table, res)
                #
                if results.count() > count:
                    count = results.count()
                    prog.update(count)
                #
                if not any(proc_status):
                    # if everything is done
                    assert results.count() == len(asso.groups)
                    break
            prog.done()
    except Exception as e:
        sys.exit(e) 

#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
def getAllTests():
    '''List all tests (all classes that subclasses of NullTest) in this module'''
    return [name for name, obj in globals().iteritems() if type(obj) == type(NullTest) and issubclass(obj, NullTest)]


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
        self.group = None
        

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''A base association test method
            test does nothing, but can be used to measure the performance of 
            retrieving data from vtools.''',
            prog='vtools associate --method ' + self.name)
        # no argumant is added
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        
    
    def getFields(self):
        '''Get updated fields for the association test.'''
        fields = []
        for key in self.result.keys():
            if self.result[key] is not None:
                fields.append(('_'.join([self.__class__.__name__, key]), 'FLOAT'))
        return fields
        
        
    def setGenotype(self, which, data):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
    
    
    def setPhenotype(self, which, data, covariates=None):
        '''Set phenotype data'''
        phen = [x for idx, x in enumerate(data) if which[idx]]
        if covariates:
          covt = [[x for idx, x in enumerate(y) if which[idx]] for y in covariates]
          self.data.setPhenotype(phen, covt)
        else:
          self.data.setPhenotype(phen)
        
    def setAttributes(self, grp):
        self.group = '__'.join(map(str, grp))

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant. Will print
        data if NullTest is called'''
        self.logger.info('Currently no action is defined for NullTest')
        return 0

class ExternTest(NullTest):
    '''A test that exports data in standard formats, call an external program
    and prase its output. This is the simplest, but also the slowest method
    to utilize a third-party association test program.
    '''
    def __init__(self, logger=None, name=None, *method_args):
        pass
    
def freq(input):
    try:
        value = float(input)
        if not (value >=0 and value <= 1):
            msg = "%r is not valid input. Valid input should fall in range [0, 1]" % input
            raise ValueError(msg)
    except Exception as e:
        raise argparse.ArgumentTypeError(e)
    return value


class LinearBurdenTest(NullTest):
    '''Simple Linear regression score test on collapsed genotypes
    within an association testing group
    '''
    def __init__(self, logger=None, name=None, *method_args):
        NullTest.__init__(self, logger, name, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo code.''',
            prog='vtools associate --method ' + self.name)
        # no argumant is added
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 1.0''')  
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2 
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if envoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations.''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')        
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 100,000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommanded to
            specify a "C" that is slightly larger than the significance level for the study.
            To not using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if envoked, will apply variable thresholds method to the permutation routine in GENE based analysis.''')
        parser.add_argument('--weight_by_maf', action='store_true',
            help='''This option, if envoked, will apply Madsen&Browning weighting (based on observed allele frequencies in all samples)
            to the permutation routine in GENE based analysis. Note this option will be masked if --use_indicator is envoked''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def calculate(self):
        data = self.data.clone()
        doRegression =  t.MultipleLinearRegression() if data.covarcounts() > 0 else t.SimpleLinearRegression()
        codeX = t.BinToX() if self.use_indicator else t.SumToX()
        if self.permutations == 0:
            actions = [t.SetMaf(), t.SetSites(self.mafupper, self.maflower), codeX, doRegression, t.StudentPval(self.alternative)]
            if self.weight_by_maf and not self.use_indicator:
                actions = [t.SetMaf(), t.SetSites(self.mafupper, self.maflower), t.WeightByAllMaf(), codeX, doRegression, t.StudentPval(self.alternative)]
            a = t.ActionExecutor(actions)
            a.apply(data)
        else:
            actions = [t.SetMaf(), t.SetSites(self.mafupper, self.maflower)]
            if self.weight_by_maf and not self.use_indicator:
                actions.append(t.WeightByAllMaf())
            permute_actions = [codeX, doRegression]
            permutator = t.VariablePermutator()
            if not self.variable_thresholds:
                permutator = t.FixedPermutator()
                permute_actions = [doRegression]
                actions.append(codeX)
            # pre-processing data
            a = t.ActionExecutor(actions)
            a.apply(data)
            # permutation
            p = permutator(self.permute_by.upper(), self.alternative, self.permutations, self.adaptive, permute_actions)
            p.apply(data)
        
        self.result['pvalue'] = data.pvalue()
        self.result['statistic'] = data.statistic()
        self.result['samples'] = data.samplecounts()
        #print "TEST", self.__class__.__name__, '\n'
        #print "PHENOTYPES", data.phenotype(), '\n'
        #print "COVARIATES", data.covariates(), '\n'
        #print "GENOTYPES", data.raw_genotype(), '\n'
        #print "GROUP", self.group, '\n'
        #print "RESULTS", self.result
        return 0
