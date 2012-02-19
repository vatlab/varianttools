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
from .project import Project, Field, AnnoDBWriter
from .utils import ProgressBar, consolidateFieldName, DatabaseEngine, delayedAction
from .phenotype import Sample
import argparse

if sys.version_info.major == 2:
    import assoTests_py2 as t
else:
    import assoTests_py3 as t

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
            is equivalent to --method "m1 -p 1000" "m2 -p 1000".). You can use
            command 'vtools show tests' for a list of association tests, and
            'vtools show test TST' for details about a test.''')
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('-g', '--group_by', nargs='*',
        help='''Group variants by fields. If specified, variants will be separated
            into groups and are tested one by one.''')
    parser.add_argument('-j', '--jobs', metavar='N', default=1, type=int,
        help='''Number of processes to carry out association tests.''')
    parser.add_argument('--to_db', metavar='annoDB',
        help='''Name of a database to which results from association tests will be written''')        


class AssociationTestManager:
    '''Parse command line and get data for association testing. This class will provide
    the following attributes for others to use:

    tests:       a list of test objects
    IDs:         sample IDs
    table:       variant table (genotype)
    phenotype:   phenotype 
    covariates:  covariates
    groups:      a list of groups
    group_names: names of the group
    '''
    def __init__(self, proj, table, phenotype, covariates, methods, unknown_args, samples, group_by):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        # table?
        if not self.proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        self.table = table
        #
        # step 0: get testers
        self.tests = self.getAssoTests(methods, unknown_args)
        # step 1: get samples and related phenotypes
        self.IDs, self.phenotype, self.covariates = self.getPhenotype(samples, phenotype, covariates)
        # step 2: indexes genotype tables if needed
        proj.db.attach('{}_genotype.DB'.format(proj.name), '__fromGeno')
        unindexed_IDs = [x for x in self.IDs if \
            not proj.db.hasIndex('__fromGeno.genotype_{}_index'.format(x))]
        if unindexed_IDs:
            cur = proj.db.cursor()
            prog = ProgressBar('Indexing genotypes', len(unindexed_IDs))
            for idx, ID in enumerate(unindexed_IDs):
                cur.execute('CREATE INDEX __fromGeno.genotype_{0}_index ON genotype_{0} (variant_id ASC)'.format(ID))
                prog.update(idx + 1)
            prog.done()
        #
        # step 3: get groups
        self.group_names, self.group_types, self.groups = self.identifyGroups(group_by)

    def getAssoTests(self, methods, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common 
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least one statistical test. Please use command "vtools show tests" for a list of tests')
        tests = []
        for m in methods:
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                method = eval(name)(self.logger, args)
                tests.append(method)
            except NameError as e:
                self.logger.debug(e)
                raise ValueError('Could not identify association test {}. Please use command "vtools show tests" for a list of tests"')
        return tests

    def getPhenotype(self, condition, pheno, covar):
        '''Get a list of samples from specified condition. This function sets self.IDs'''
        try:
            query = 'SELECT sample_id, {} FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id'.format(
                ', '.join(pheno + (covar if covar is not None else []))) + \
                (' WHERE {}'.format(' AND '.join(['({})'.format(x) for x in condition])) if condition else '') + ' ORDER BY sample_id;'
            self.logger.debug('Select phenotype and covariates using query {}'.format(query))
            cur = self.db.cursor()
            cur.execute(query)
            data = list(zip(*cur.fetchall()))
            IDs, phenotypes, covariates = data[0], data[1: (1 + len(pheno))], data[(1 + len(pheno)) : ]
            if len(IDs) == 0:
                raise ValueError('No sample is selected by condition: {}'.format(' AND '.join(['({})'.format(x) for x in condition])))
            else:
                self.logger.info('{} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(['({})'.format(x) for x in condition])))
            # add intersection
            covariates.insert(0, array('d', [1]*len(IDs)))
            return IDs, phenotypes, covariates
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('''Failed to retrieve phenotype {} and covariate {}. Please 
                             make sure the specified phenotype names are correct and there 
                             is no missing value'''.format(', '.join(pheno), ', '.join(covar)))
    
    def identifyGroups(self, group_by):
        '''Get a list of groups according to group_by fields'''
        if not group_by:
            group_by = ['chr', 'pos']
        # find the source of fields in group_by
        table_of_fields = [self.proj.linkFieldToTable(field, self.table)[-1].table for field in group_by]
        table_of_fields = [x if x else self.table for x in table_of_fields] 
        # name the fields
        field_names = [x.replace('.', '_') for x in group_by]
        # type of the fields
        field_types = [self.db.typeOfColumn(y, x.rsplit('.', 1)[-1]) for x,y in zip(group_by, table_of_fields)]
        cur = self.db.cursor()
        #
        # create a table that holds variant ids and groups, indexed for groups.
        # The structure of this table will look like
        # 
        #   variant_id  INT NOT NULL,
        #   chr VARCHAR(20),
        #   pos INT
        #
        cur.execute('DROP TABLE IF EXISTS __asso_tmp;')
        cur.execute('DROP INDEX IF EXISTS __asso_tmp_index;')
        cur.execute('''\
            CREATE TABLE __asso_tmp (
              variant_id INT NOT NULL,
              {});
              '''.format(','.join(['{} {}'.format(x,y) for x,y in zip(field_names, field_types)])))
        #
        # select variant_id and groups for association testing
        group_fields, fields = consolidateFieldName(self.proj, self.table, ','.join(group_by))        
        from_clause = [self.table]
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
        s = delayedAction(self.logger.info, "Creating association test table (please be patient ...)")
        self.logger.debug('Running query {}'.format(query))
        cur.execute(query)
        cur.execute('''\
            CREATE INDEX __asso_tmp_index ON __asso_tmp ({});
            '''.format(','.join(['{} ASC'.format(x) for x in field_names])))
        del s
        # get group by
        cur.execute('''\
            SELECT DISTINCT {} FROM __asso_tmp;
            '''.format(', '.join(field_names)))
        groups = cur.fetchall()
        self.logger.info('{} groups are found'.format(len(groups)))
        #
        # the output can be too long
        #
        #self.logger.debug('Group by: {}'.format(', '.join(map(str, groups))))
        return field_names, field_types, groups


class ResultRecorder:
    def __init__(self, params, db_name=None, logger=None):
        self.completed = 0
        #
        self.group_names = params.group_names
        self.fields = []
        for n,t in zip(params.group_names, params.group_types):
            self.fields.append(Field(name=n, index=None, type=t, adj=None, comment=n))
        for test in params.tests:
            if test.name:
                self.fields.extend([
                    Field(name='{}_{}'.format(x.name, test.name), index=None,
                        type=x.type, adj=None, comment=x.comment) for x in test.fields])
            else:
                self.fields.extend(test.fields)
        if len(self.fields) != len(set([x.name for x in self.fields])):
            raise ValueError('Duplicate field names. Please rename one of the tests using parameter --name')
        print('#' + '\t'.join(self.group_names + [x.name for x in self.fields]))
        #
        self.writer = None
        if db_name:
            self.writer = AnnoDBWriter(db_name, self.fields, 
                'field',                       # field annotation databases 
                'Annotation database used to record results of association tests. Created on {}'.format(
                    time.strftime('%a, %d %b %Y %H:%M:%S', time.gmtime())),
                '1.0',                         # version 1.0 
                {'*': self.group_names},       # link by group fields
                logger, 
                True                           # allow updating an existing database
            )
            #
            self.cur = self.writer.db.cursor()
            if self.writer.update_existing:
                #
                self.update_query = 'UPDATE {0} SET {1} WHERE {2};'.format(db_name,
                    ', '.join(['{}={}'.format(x.name, self.writer.db.PH) for x in self.fields[len(self.group_names):]]),
                    ' AND '.join(['{}={}'.format(x, self.writer.db.PH) for x in self.group_names]))
                self.insert_query = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(db_name,
                    ','.join([x.name for x in self.fields]),
                    ','.join([self.writer.db.PH] * len(self.fields)))
            else:
                self.insert_query = 'INSERT INTO {0} VALUES ({1});'.format(db_name,
                    ','.join([self.writer.db.PH] * len(self.fields)))

    def record(self, res):
        self.completed += 1
        output = '\t'.join(map(str, res))
        print(output)
        # also write to an annotation database?
        if self.writer:
            if self.writer.update_existing:
                self.cur.execute(self.update_query, res[len(self.group_names):] + res[:len(self.group_names)])
                if self.cur.rowcount == 0:
                    self.cur.execute(self.insert_query, res)
            else:
                self.cur.execute(self.insert_query, res)
        
    def count(self):
        return self.completed

    def done(self):
        if self.writer:
            self.writer.finalize()
        

class AssoTestsWorker(Process):
    '''Association test calculator'''
    def __init__(self, param, grpQueue, output):
        self.proj = param.proj
        self.table = param.table
        self.IDs = param.IDs
        self.phenotypes = param.phenotype
        self.covariates = param.covariates
        self.tests = param.tests
        self.group_names = param.group_names
        self.queue = grpQueue
        self.output = output
        self.logger = self.proj.logger
        self.db = None
        Process.__init__(self, name='Phenotype association analysis for a group of variants')

    def getGenotype(self, group):
        '''Get genotype for variants in specified group'''
        # get variant_id
        where_clause = ' AND '.join(['{0}={1}'.format(x, self.db.PH) for x in self.group_names])
        query = 'SELECT variant_id FROM __asso_tmp WHERE {}'.format(where_clause) 
        #
        #self.logger.debug('Running on group {0} query {1}'.format(group, query))
        cur = self.db.cursor()
        cur.execute(query, group)
        variant_id = [x[0] for x in cur.fetchall()]
        numSites = len(variant_id)
        #self.logger.debug('Getting {0} variants for {1}'.format(numSites, group))
        #
        # get genotypes
        genotype = []
        for ID in self.IDs:
            query = 'SELECT variant_id, GT FROM __fromGeno.genotype_{0} WHERE variant_id IN (SELECT variant_id FROM __asso_tmp WHERE {1});'\
                .format(ID, where_clause)
            cur.execute(query, group)
            gtmp = {x[0]:x[1] for x in cur.fetchall()}
            # handle missing values
            gtmp = [gtmp.get(x, -9.0) for x in variant_id]
            # handle -1 coding (double heterozygotes)
            gtmp = [2.0 if int(x) == -1 else x for x in gtmp]
            genotype.append(array('d', gtmp))
        #
        self.logger.debug('Retrieved genotypes for {} samples'.format(len(genotype)))
        missing_counts = [x.count(-9.0) for x in genotype]
        # remove individuals having many missing genotypes, or have all missing variants
        # FIXME will pass it as an input arguement later
        #toKeep = [(x<0.5*numSites) for x in missing_counts]
        toKeep = [(x<numSites) for x in missing_counts]
        self.logger.debug('{} samples will be removed due to missing genotypes'.format(len(self.IDs)-sum(toKeep)))
        return genotype, toKeep
    

    def run(self):
        self.db = DatabaseEngine()
        self.db.connect(self.proj.name + '.proj', readonly=True)
        self.db.attach(self.proj.name + '_genotype.DB', '__fromGeno')
        #
        while True:
            grp = self.queue.get()
            self.logger.debug('Got group {}'.format(grp))
            if grp is None:
                self.output.send(None)
                break
            # select variants from each group:
            genotype, which = self.getGenotype(grp)
            values = list(grp)
            try:
                for test in self.tests:
                    test.setPhenotype(which, self.phenotypes, self.covariates)
                    test.setGenotype(which, genotype)
                    test.setAttributes(grp)
                    result = test.calculate()
                    self.logger.debug('Finish test')
                    values.extend(result)
            except Exception as e:
                self.logger.info('Error processing group {}, {}'.format(grp, e))
            self.logger.debug('Finished group {}'.format(grp))
            self.output.send(values)
        self.db.detach('__fromGeno')
        

def associate(args, reverse=False):
    try:
        with Project(verbosity=args.verbosity) as proj:
            asso = AssociationTestManager(proj, args.table, args.phenotype, args.covariates, args.methods, args.unknown_args,
                args.samples, args.group_by)
            #
            nJobs = max(min(args.jobs, len(asso.groups)), 1)
            # step 4: start all workers
            grpQueue = Queue()
            results = ResultRecorder(asso, args.to_db, proj.logger)
            readers = []
            for j in range(nJobs):
                r, w = Pipe(False)
                AssoTestsWorker(asso, grpQueue, w).start()
                readers.append(r)
            # put all jobs to queue, the workers will work on them
            for grp in asso.groups:
                grpQueue.put(grp)
            # the worker will stop once all jobs are finished
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
                        results.record(res)
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
            results.done()
    except Exception as e:
        sys.exit(e) 

#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
def getAllTests():
    '''List all tests (all classes that subclasses of NullTest) in this module'''
    return [(name, obj) for name, obj in globals().iteritems() if type(obj) == type(NullTest) and issubclass(obj, NullTest) and name != 'NullTest']


class NullTest:
    '''A base class that defines a common interface for association tests'''
    def __init__(self, logger=None, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to 
        parse it'''
        self.data = t.AssoData()
        self.logger = logger
        self.parseArgs(*method_args)
        #
        self.result = {}
        self.group = None

    def parseArgs(self, method_args):
        # this function should never be called.
        raise SystemError('All association tests should define their own parseArgs function')
    
    def setGenotype(self, which, data):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
    
    
    def setPhenotype(self, which, data, covariates=None):
        '''Set phenotype data'''
        phen = [x for idx, x in enumerate(data[0]) if which[idx]]
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
        #self.logger.debug('Printing out group name, phenotypes, covariates and genotypes')
        #self.logger.debug(self.group)
        #self.logger.debug(self.data.phenotype())
        #self.logger.debug(self.data.covariates())
        #self.logger.debug(self.data.raw_genotype())
        return 0

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
    '''Simple Linear regression score test on collapsed genotypes within an association testing group '''
    def __init__(self, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = [
            Field(name='p_value', index=None, type='FLOAT', adj=None, comment='p-value'),
            Field(name='statistic', index=None, type='FLOAT', adj=None, comment='statistic'),
            Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size')]

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
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
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')        
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To not using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in GENE based analysis''')
        parser.add_argument('--weight_by_maf', action='store_true',
            help='''This option, if evoked, will apply Madsen&Browning weighting (based on observed allele frequencies in all samples)
            to GENE based analysis. Note this option will be masked if --use_indicator is evoked''')
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
            p = t.VariablePermutator(self.permute_by.upper(), self.alternative, self.permutations, self.adaptive, permute_actions)
            if not self.variable_thresholds:
                permute_actions = [doRegression]
                p = t.FixedPermutator(self.permute_by.upper(), self.alternative, self.permutations, self.adaptive, permute_actions)
                actions.append(codeX)
            # pre-processing data
            a = t.ActionExecutor(actions)
            a.apply(data)
            # permutation
            p.apply(data)
        
        return data.pvalue(), data.statistic(), data.samplecounts()

class AliasTest(LinearBurdenTest):
    '''An example of a specialized linear burden test '''
    def __init__(self, logger=None, *method_args):
        LinearBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
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
        #
        # For example, if this parameter is fixed...
        #
        #parser.add_argument('--use_indicator', action='store_true',
        #    help='''This option, if evoked, will apply binary coding to genotype groups
        #    (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')        
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To not using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in GENE based analysis''')
        parser.add_argument('--weight_by_maf', action='store_true',
            help='''This option, if evoked, will apply Madsen&Browning weighting (based on observed allele frequencies in all samples)
            to GENE based analysis. Note this option will be masked if --use_indicator is evoked''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = True


