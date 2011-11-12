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
            self.phenotype = [array('d', map(float, x)) for x in zip(*cur.fetchall())[1:]]
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to retrieve phenotype {}'.format(', '.join(phenotype)))
    
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
            self.covariates = [array('d', map(float, x)) for x in zip(*cur.fetchall())[1:]]
            self.covariates.insert(0, array('d', [1]*len(self.covariates[0])))
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to retrieve phenotype covariate {}'.format(', '.join(covariates)))
    
    def identifyGroups(self, group_by):
        '''Get a list of groups according to group_by fields'''
        self.group_by = group_by
        if not group_by:
            self.groups = ['all']
            self.where_clause = None
            self.from_clause = None
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
            # FIXME the code here is buggy ...
            query = 'SELECT DISTINCT {} FROM {} {};'.format(group_fields,
                self.from_clause, self.where_clause) 
            self.logger.debug('Running query {}'.format(query))
            # get group by
            cur = self.db.cursor()
            cur.execute(query)
            self.groups = [x[0] for x in cur.fetchall()]
            self.logger.info('Find {} groups'.format(len(self.groups)))
            self.logger.debug('Group by: {}'.format(', '.join([str(x) for x in self.groups])))

#class StatStatus:
#    def __init__(self):
#        self.tasks = {}
#        self.lock = Lock()
#
#    def set(self, task, value):
#        self.tasks[task] = value
#    
#    def get(self, task):
#        return self.tasks[task]
#        
#    def count(self):
#        return len(self.tasks)

class UpdateResult:
    def __init__(self, db, logger, fn):
        self.grps = []
        self.lock = Lock()
        self.db = db
        self.logger = logger
        self.fn = fn[0] if fn else None

    def set(self, table, res):
        self.grps.append(res[0])
#        self.lock.acquire()
#        startID, endID = res[1][0]
#        fields = []
#        results = []
#        for x,y in res[1][1:]:
#            fields.extend(x)
#            results.append(y)
#        headers = self.db.getHeaders(table)
#        for field, fldtype in fields:
#            if field not in headers:
#                self.logger.debug('Adding field {}'.format(field))
#                self.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(table, field, fldtype))
#                self.db.commit()
#        # if the field exists it will be re-written
#        update_query = 'UPDATE {0} SET {2} WHERE variant_id>={1} AND variant_id<={1};'.format(table, self.db.PH,
#        ', '.join(['{}={}'.format(field, self.db.PH) for field, fldtype in fields]))
#        # fill up the update query and execute the update command
#        names = [x.split('_')[1] for x,y in fields]
#        values = []
#        for result in results:
#            values.extend([result[x] for x in names])
#            self.logger.debug('Running query {}'.format(update_query))
#            self.db.execute(update_query, values+[startID, endID])
#        self.db.commit()
#        self.lock.release()
        if not self.fn:
            self.fn = 'vtoolsasso_{}.result'.format(time.strftime('%b%d_%H%M%S', time.gmtime()))
        try:
            output = '[{0}]\ntest = {1}\np-value = {2}\nstatistic = {3}\nsamples = {4}\n'.format(res[0],
                        res[1][0][0][0][0].split('_')[0], res[1][0][1]['pvalue'], res[1][0][1]['statistic'], res[1][0][1]['samples'])
            self.lock.acquire()
            with open(self.fn, 'a') as out:
                out.write(output)
            self.lock.release()
        except IndexError as e:
            self.logger.info('No association tests done for {}'.format(self.grps[-1]))
            pass
        return
        
    def getgrp(self):
        return self.grps
    
    def count(self):
        return len(self.grps)
        
class GroupAssociationCalculator(Process):
    '''Association test calculator'''
    def __init__(self, proj, table, samples, phenotypes, covariates, tests, grpfields, grpQueue, where, from_clause, output):
        self.proj = proj
        self.table = table
        self.IDs = samples
        self.phenotypes = phenotypes
        self.covariates = covariates
        self.tests = tests
        self.pjname = proj.name
        self.group_fields = grpfields
        self.queue = grpQueue
        self.output = output
        self.logger = proj.logger
        self.db = None
        self.where_clause = where
        self.from_clause = from_clause
        Process.__init__(self, name='Phenotype association analysis for a group of variants')
        
    def getVariants(self, group):
        '''Get variants for a certain group'''
        if not self.group_fields:
            # group must be 'all'
            return self.table
        vtable = '__asso_tmp_{}'.format(group)
        if self.db.hasTable(vtable):
            self.db.truncateTable(vtable)
        else:
            self.db.execute('''CREATE TEMPORARY TABLE {} (
                variant_id INTEGER PRIMARY KEY);'''.format(vtable))
            self.db.commit()
        #
        where_clause = (self.where_clause + ' AND ' if self.where_clause else 'WHERE ') + \
            ' AND '.join(['{}={}'.format(x, self.db.PH) for x in self.group_fields])
        query = 'INSERT INTO __asso_tmp_{} SELECT {}.variant_id FROM {} {};'.format(group, self.table, self.from_clause, where_clause)
        self.logger.debug('Running query {}'.format(query))
        cur = self.db.cursor()
        cur.execute(query, (group,))
        return vtable

    def getGenotype(self, vtable):
        '''Get genotype for variants in specified table'''
        genotype = []
        self.db.attach(self.pjname+'_genotype.DB', '__fromGeno')
        cur = self.db.cursor()
        cur.execute('SELECT variant_id FROM {}'.format(vtable))
        variant_id = cur.fetchall()[0]
        numSites = len(variant_id)
        for ID in self.IDs:
            query = 'SELECT variant_id, GT FROM __fromGeno.genotype_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                .format(ID, vtable)
            #self.logger.debug('Running query {}'.format(query))
            cur.execute(query)
            gtmp = {x[0]:x[1] for x in cur.fetchall()}
            genotype.append(array('d', [gtmp.get(x, -9.0) for x in variant_id]))
        self.db.detach('__fromGeno')
        missing_counts = [x.count(-9.0) for x in genotype]
        # remove individuals having many missing genotypes
        toKeep = [(x<0.5*numSites) for x in missing_counts]
        self.logger.debug('{} samples will be removed due to missing genotypes'.format(len(self.IDs)-sum(toKeep)))
        return genotype, toKeep
    
    def run(self):
        self.db = DatabaseEngine()
        self.db.connect(self.pjname+'.proj')
        for annoDB in self.proj.annoDB:
           self.db.attach(os.path.join(annoDB.dir, annoDB.filename))
        #
        while True:
            self.logger.debug('Getting group ...') 
            grp = self.queue.get()
            self.logger.debug('Got group {}'.format(grp))
            if grp is None:
                self.logger.debug('All groups are processed, sending None')
                self.output.send(None)
                break
            # select variants from each group:
            vtable = self.getVariants(grp)
            genotype, which = self.getGenotype(vtable)
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
                self.logger.info('Error processing {}, {}'.format(grp, e))
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
                    asso.tests, args.group_by, grpQueue, asso.where_clause, asso.from_clause, w).start()
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
        self.result = {'pvalue':None, 'statistic':None}
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
        self.data.mean_phenotype()
        self.data.count_cases()
        self.data.count_ctrls()
        
    def setAttributes(self, grp):
        self.group = str(grp)

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant. Will print
        data if NullTest is called'''
        self.logger.info('Currently no action is defined for NullTest')
#        print('Group name: {}\n'.format(self.group)+'Phenotype-genotype data:')
#        print(' '.join(map(str, self.data.phenotype())))
#        print('\n'.join([' '.join(map(str, map(int, x))) for x in self.data.raw_genotype()])+'\n')
        return 0

class ExternTest(NullTest):
    '''A test that exports data in standard formats, call an external program
    and prase its output. This is the simplest, but also the slowest method
    to utilize a third-party association test program.
    '''
    def __init__(self, logger=None, name=None, *method_args):
        pass

class LinearBurdenTest(NullTest):
    '''Simple Linear regression score test on collapsed genotypes
    within an association testing group
    '''
    def __init__(self, logger=None, name=None, *method_args):
        NullTest.__init__(self, logger, name, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test
            which will collapse the variants within one group by counts and evaluate
            the significance of effect size (regression coefficient) of the group''',
            prog='vtools associate --method ' + self.name)
        # no argumant is added
        parser.add_argument('-p', '--permutations', type=int, default=0,
            help='''Number of permutations.''')
        parser.add_argument('-q1', '--mafupper', type=float, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 1.0''')  
        parser.add_argument('-q2', '--maflower', type=float, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2 
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', type=int, default=1,
            help='''Alternative hypothesis is one-sided (1) or two-sided (2).
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def calculate(self):
        data = self.data.clone()
        doRegression = t.SimpleLinearRegression()
        task_dbg = "-simple regression"
        if data.covarcounts() > 0:
            task_dbg = "-multiple regression"
            doRegression = t.MultipleLinearRegression()
        actions = [t.SetMaf(), t.FilterX(self.mafupper, self.maflower), t.SumToX(), doRegression, t.GaussianPval(self.alternative)]
        a = t.ActionExecuter(actions)
        a.apply(data)
        #print('{} on group {}, p-value (asymptotic) = {}'\
        #                 .format(self.__class__.__name__, self.group, data.pvalue()))
        # permutation 
        if not self.permutations == 0:
            #self.logger.info('permutation routine no ready')
            p = t.PhenoPermutator(self.permutations, [t.SimpleLinearRegression()])
            p.apply(data)
        #  print('{} on group {}, p-value (permutation) = {}'\
        #                .format(self.__class__.__name__, self.group, (p.apply(data)+1.0) / (self.permutations+1.0)))
        self.result['pvalue'] = data.pvalue()
        self.result['statistic'] = data.statistic()
        self.result['samples'] = data.samplecounts()
        #print data.raw_genotype(), '\n'
        #print data.phenotype(), '\n'
        #print data.covariates()
        #self.logger.debug('Finished test {0} for group {1}, {2}'.format(self.__class__.__name__+task_dbg, self.group, self.result))
        return 0
