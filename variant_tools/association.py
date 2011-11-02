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
import threading
import Queue
import time
from array import array
from .project import Project
from .utils import ProgressBar, consolidateFieldName
from .sample import Sample
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

class AssociationTester(Sample):
    '''Parse command line and get data for association testing'''
    
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
        return genotype, start, end
    
    def updateTestResult(self, test, startID, endID):

        if test.__class__.__name__ == 'NullTest':

        
class StatStatus:
    def __init__(self):
        self.tasks = {}
        self.lock = threading.Lock()

    def set(self, task, value):
        self.lock.acquire()
        self.tasks[task] = value
        self.lock.release()
    
    def get(self, task):
        return self.tasks[task]
        
    def count(self):
        return len(self.tasks)
        
class GroupAssociationCalculator(threading.Thread):
    def __init__(self, asso, grpQueue, status, logger):
        self.asso = asso
        self.queue = grpQueue
        self.status = status
        self.logger = logger
        threading.Thread.__init__(self, name='Phenotype association analysis for a group of variants')
    
    def run(self):
        while True:
            grp = self.queue.get()
            if grp is None:
                self.queue.task_done()
                break
            # select variants from each group:
            vtable = asso.getVariants(grp)
            genotype, startID, endID = asso.getGenotype(vtable)
            values = [(startID, endID)]
            for test in asso.tests:
                test.setGenotype(genotype)
                test.setAttributes(grp)
                test.calculate()
                values.append((test.getFields(), test.result))
            self.status.set(grp, values)
            self.queue.task_done()
        

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
            nJobs = max(min(args.jobs, len(asso.groups)), 1)
            # step 4: start all workers
            grpQueue = Queue.Queue()
            status = StatStatus()
            proj.db.attach(proj.name + '_genotype')
            for j in range(nJobs):
                GroupAssociationCalculator(asso, grpQueue, status, proj.logger).start()
            # put all jobs to queue, the workers will work on them
            for grp in asso.groups:
                idQueue.put(grp)

            count = 0
            prog = ProgressBar('Testing for association', len(asso.groups))
            while True:
                if status.count() > count:
                    count = status.count()
                    prog.update(count)
                # if everything is done
                if status.count() == len(IDs):
                    # stop all threads
                    for j in range(nJobs):
                        idQueue.put(None)
                    break
                # wait 1 sec to check status again
                time.sleep(1)
            prog.done()
        
            # step 5: update variant table by groups
            prog = ProgressBar('Updating {}'.format(args.table), len(asso.groups))
            for igrp, grp in enumerate(asso.groups):
                res = status.get(grp)
                startID, endID = res[0]
                fields = [x[0] for x in res[1:]]
                results = [x[1] for x in res[1:]]
                headers = self.db.getHeaders(args.table)
                for field, fldtype in fields:
                    if field not in headers:
                        proj.logger.info('Adding field {}'.format(field))
                        proj.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(proj.table, field, fldtype))
                        proj.db.commit()
                # if the field exists it will be re-written
                update_query = 'UPDATE {0} SET {2} WHERE variant_id={1};'.format(proj.table, proj.db.PH,
                ', '.join(['{}={}'.format(field, proj.db.PH) for field, fldtype in fields]))
                # fill up the update query and execute the update command
                names = [x.split('_')[1] for x,y in fields]
                for idx, id in enumerate(range(startID, endID+1)):
                    values = []
                    for result in results:
                        values.extend([result[x][idx] for x in names])
                    self.db.execute(update_query, values+[id])
                
                if igrp % 100 == 0:
                    proj.db.commit()
                    prog.update(idx)
            proj.db.commit()
            prog.done()                
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

    def setPhenotype(self, data):
        '''Set phenotype data'''
        self.data.setPhenotype(data)

    def setGenotype(self, data):
        self.data.setGenotype(data)
    
    def setAttributes(self, grp):
        self.group = repr(grp[0].encode('ascii'))
        self.data.setMaf()
        self.data.filterByMaf(upper=1.0, lower=0.0)
        self.data.mean_phenotype()
        self.data.count_cases()
        self.data.count_ctrls()

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant. Will print
        data if NullTest is called'''
        print('Group name: {}\n'.format(self.group)+'Phenotype-genotype data:')
        print(' '.join(map(str, self.data.phenotype())))
        print('\n'.join([' '.join(map(str, map(int, x))) for x in self.data.raw_genotype()])+'\n')
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
   
    def calculate(self):
        data = self.data.clone()
        actions = [t.SumToX(), t.SimpleLogisticRegression(), t.GaussianPval(1)]
        a = t.ActionExecuter(actions)
        a.apply(data)
        self.logger.debug('{} on group {}, p-value (asymptotic) = {}'\
                         .format(self.__class__.__name__, self.group, data.pvalue()))
        # logistic regression, permutation 
        p = t.PhenoPermutator(self.ptime, [t.SimpleLogisticRegression()])
        self.logger.debug('{} on group {}, p-value (permutation) = {}'\
                        .format(self.__class__.__name__, self.group, p.apply(data) / float(self.ptime)))
        self.result['pvalue'] = [data.pvalue()]*len(self.data.raw_genotype()[0])
        self.result['statistic'] = [data.statistic()]*len(self.data.raw_genotype()[0])
        return 1
