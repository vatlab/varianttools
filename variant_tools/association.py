#!/usr/bin/env python
#
# $File: association.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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
import math
from collections import OrderedDict
from copy import copy, deepcopy
from .project import Project, Field, AnnoDB, AnnoDBWriter
from .utils import ProgressBar, consolidateFieldName, DatabaseEngine, delayedAction
from .phenotype import Sample
import argparse

if sys.version_info.major == 2:
    import assoTests_py2 as t
else:
    import assoTests_py3 as t

def associateArguments(parser):
    parser.add_argument('table', help='''Variant table.''')
    parser.add_argument('phenotypes', nargs=1,
        help='''A list of phenotypes that will be passed to the association
            statistics calculator. Currently only a single phenotype is allowed.''')
    parser.add_argument('--covariates', nargs='*', default=[],
        help='''Optional phenotypes that will be passed to statistical
            tests as covariates. Values of these phenotypes should be integer
            or float.''')
    parser.add_argument('--var_info', nargs='*', default=[],
        help='''Optional variant information fields (e.g. minor allele frequency
            from 1000 genomes project) that will be passed to statistical tests.
            The fields could be any annotation fields of with integer or float
            values, including those from used annotation databases (use "vtools
            show fields" to see a list of usable fields). ''')
    parser.add_argument('--geno_info', nargs='*', default=[],
        help='''Optional genotype fields (e.g. quality score of genotype calls,
            c.f. "vtools show genotypes") that will be passed to statistical
            tests. Note that the fields should exist for all samples that are
            tested.''')
    parser.add_argument('-m', '--methods', nargs='+',
        help='''Method of one or more association tests. Parameters for each
            method should be specified together as a quoted long argument (e.g.
            --method "m --alternative 2" "m1 --permute 1000"), although
            the common method parameters can be specified separately, as long as
            they do not conflict with command arguments. (e.g. --method m1 m2 -p 1000 
            is equivalent to --method "m1 -p 1000" "m2 -p 1000".). You can use
            command 'vtools show tests' for a list of association tests, and
            'vtools show test TST' for details about a test. Customized association
            tests can be specified as mod_name.test_name where mod_name should
            be a Python module (system wide or in the current directory), and
            test_name should be a subclass of NullTest.''')
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"'). Each line of the sample table (vtools show
            samples) is considered as samples. If genotype of a physical sample
            is scattered into multiple samples (e.g. imported chromosome by 
            chromosome), they should be merged using command vtools admin.''')
    parser.add_argument('-g', '--group_by', nargs='*',
        help='''Group variants by fields. If specified, variants will be separated
            into groups and are tested one by one.''')
    parser.add_argument('-j', '--jobs', metavar='N', default=1, type=int,
        help='''Number of processes to carry out association tests.''')
    parser.add_argument('--to_db', metavar='annoDB',
        help='''Name of a database to which results from association tests will be written''')        
    parser.add_argument('--update', action='store_true',
        help='''When --to_db is specified, allow updating existing fields in the result database''')


class AssociationTestManager:
    '''Parse command line and get data for association testing. This class will provide
    the following attributes for others to use:

    tests:       a list of test objects
    sample_IDs:  sample IDs
    table:       variant table (genotype)
    phenotypes:  phenotypes 
    covariates:  covariates
    var_info:    variant information fields
    geno_info:   genotype information fields
    groups:      a list of groups
    group_names: names of the group
    '''
    def __init__(self, proj, table, phenotypes, covariates, var_info, geno_info, methods,
        unknown_args, samples, group_by):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        self.var_info = var_info
        self.geno_info = geno_info
        # table?
        if not self.proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        self.table = table
        #
        # step 0: get testers
        self.tests = self.getAssoTests(methods, len(covariates), unknown_args)
        # step 1: get samples and related phenotypes
        self.sample_IDs, self.phenotypes, self.covariates = self.getPhenotype(samples, phenotypes, covariates)
        # step 2: check if tests are compatible with phenotypes
        for idx, item in enumerate([list(set(x)) for x in self.phenotypes]):
            if (list(map(float, item)) == [2.0, 1.0] or list(map(float, item)) == [1.0, 2.0]):
                self.phenotypes[idx] = [i - 1.0 for i in self.phenotypes[idx]]
                item = [i - 1.0 for i in item]
            if not (list(map(float, item)) == [0.0, 1.0] or list(map(float, item)) == [1.0, 0.0]):
                for test in self.tests:
                    if test.model == 1:
                        raise ValueError("Cannot perform logistic regression on non-binary coding phenotype")
        # step 3: indexes genotype tables if needed
        proj.db.attach('{}_genotype.DB'.format(proj.name), '__fromGeno')
        unindexed_IDs = []
        for id in self.sample_IDs:
            if not proj.db.hasIndex('__fromGeno.genotype_{}_index'.format(id)):
                unindexed_IDs.append(id)
        if unindexed_IDs:
            cur = proj.db.cursor()
            prog = ProgressBar('Indexing genotypes', len(unindexed_IDs))
            for idx, ID in enumerate(unindexed_IDs):
                cur.execute('CREATE INDEX __fromGeno.genotype_{0}_index ON genotype_{0} (variant_id ASC)'.format(ID))
                prog.update(idx + 1)
            prog.done()
        #
        # step 4: get groups
        self.group_names, self.group_types, self.groups = self.identifyGroups(group_by)

    def getAssoTests(self, methods, ncovariates, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common 
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError('Please specify at least one statistical test. Please use command "vtools show tests" for a list of tests')
        tests = []
        for m in methods:
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                if '.' in name:
                    # if the method is defined elsewhere
                    m_module, m_name = name.split('.', 1)
                    # also search current working directory
                    my_dir = os.getcwd()
                    if my_dir not in sys.path:
                        sys.path.append(my_dir)
                        _temp = __import__(m_module, globals(), locals(), [m_name], -1)
                        sys.path.pop()
                    else:
                        _temp = __import__(m_module, globals(), locals(), [m_name], -1)
                    method = getattr(_temp, m_name)(ncovariates, self.logger, args)
                else:
                    method = eval(name)(ncovariates, self.logger, args)
                # check if method is valid
                if not hasattr(method, 'fields'):
                    raise ValueError('Invalid association test method {}: missing attribute fields'.format(name))
                if not method.fields:
                    self.logger.warning('Association test {} has invalid  or empty fields. No result will be generated.'.format(name))
                tests.append(method)
            except NameError as e:
                self.logger.debug(e)
                raise ValueError('Could not identify association test {}. Please use command "vtools show tests" for a list of tests'.format(name))
        return tests

    def getPhenotype(self, condition, pheno, covar):
        '''Get a list of samples from specified condition. This function sets self.sample_IDs, self.phenotypes and self.covariates'''
        try:
            query = 'SELECT sample_id, {} FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id'.format(
                ', '.join(pheno + (covar if covar is not None else []))) + \
                (' WHERE {}'.format(' AND '.join(['({})'.format(x) for x in condition])) if condition else '')
            self.logger.debug('Select phenotype and covariates using query {}'.format(query))
            cur = self.db.cursor()
            cur.execute(query)
            data = []
            for rec in cur:
                # get id, phenotype, and covariants
                data.append([rec[0], rec[1 : (1 + len(pheno))], rec[ (1 + len(pheno)) : ]])
            sample_IDs = []
            phenotypes = [[] for x in pheno]
            covariates = [[] for x in covar]
            for i, p, c in data:
                sample_IDs.append(i)
                [x.append(y) for x,y in zip(phenotypes, p)]
                [x.append(y) for x,y in zip(covariates, c)]
            if len(sample_IDs) == 0:
                raise ValueError('No sample is selected by condition: {}'.format(' AND '.join(['({})'.format(x) for x in condition])))
            else:
                if len(condition) > 0:
                    self.logger.info('{} samples are selected by condition: {}'.format(len(sample_IDs), ' AND '.join(['({})'.format(x) for x in condition])))
                else:
                    self.logger.info('{} samples are found'.format(len(sample_IDs)))
            # add intercept
            covariates.insert(0, [1]*len(sample_IDs))
            try:    
                phenotypes = [map(float, x) for x in phenotypes]
                covariates = [map(float, x) for x in covariates]
            except ValueError:
                raise ValueError('Invalid coding in phenotype/covariates values: '
                                 'missing values should be removed from analysis or '
                                 'inferred with numeric values')
            return sample_IDs, phenotypes, covariates
        except Exception as e:
            self.logger.debug(e)
            if str(e).startswith('Invalid coding'):
                raise ValueError(e)
            else:
                raise ValueError('Failed to retrieve phenotype {}{}. Please use command '
                    '"vtools show samples" to see a list of phenotypes'.format(', '.join(pheno), 
                    '' if covar is None else (' and/or covariate' + ', '.join(covar))))
    
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
    def __init__(self, params, db_name=None, update=False, logger=None):
        self.succ_count = 0
        self.failed_count = 0
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
        for field in self.fields:
            if '-' in field.name:
                raise ValueError('"-" is not allowed in field name {}'.format(field.name))
        if len(self.fields) != len(set([x.name for x in self.fields])):
            raise ValueError('Duplicate field names. Please rename one of the tests using parameter --name')
        print('#' + '\t'.join([x.name for x in self.fields]))
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
                True,                          # allow updating an existing database
                update                    # allow updating an existing field
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
        self.succ_count += 1
        if not res:
            self.failed_count += 1
            return
        str_output = '\t'.join(['{0:G}'.format(x, precision=5) if isinstance(x, float) else str(x) for x in res])
        print(str_output)
        # also write to an annotation database?
        if self.writer:
            if self.writer.update_existing:
                self.cur.execute(self.update_query, res[len(self.group_names):] + res[:len(self.group_names)])
                if self.cur.rowcount == 0:
                    self.cur.execute(self.insert_query, res)
            else:
                self.cur.execute(self.insert_query, res)
        
    def completed(self):
        return self.succ_count

    def failed(self):
        return self.failed_count

    def done(self):
        if self.writer:
            self.writer.finalize()
        

class AssoTestsWorker(Process):
    '''Association test calculator'''
    def __init__(self, param, grpQueue, output):
        self.proj = param.proj
        self.table = param.table
        self.sample_IDs = param.sample_IDs
        self.phenotypes = param.phenotypes
        self.covariates = param.covariates
        self.var_info = param.var_info
        self.geno_info = param.geno_info
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
        # get annotation
        var_info = {x:[] for x in self.var_info}
        if self.var_info:
            select_clause, fields = consolidateFieldName(self.proj, 'variant', ','.join(['variant_id'] + self.var_info))
            #
            # FROM clause
            from_clause = 'FROM variant '
            fields_info = sum([self.proj.linkFieldToTable(x, 'variant') for x in fields], [])
            #
            processed = set()
            # the normal annotation databases that are 'LEFT OUTER JOIN'
            for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
                if (tbl.lower(), conn.lower()) not in processed and '.__' not in tbl:
                    from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                    processed.add((tbl.lower(), conn.lower()))
            #
            # query
            query = 'SELECT {} {} WHERE variant.variant_id IN (SELECT variant_id FROM __asso_tmp WHERE {})'.format(
                select_clause, from_clause, where_clause)
            #self.logger.debug('Running query: {}'.format(query))
            cur.execute(query, group)
            #
            data = {x[0]:x[1:] for x in cur.fetchall()}
            for idx, key in enumerate(self.var_info):
                var_info[key] = [data[x][idx] for x in variant_id]
        #
        # get genotypes
        genotype = []
        geno_info = {x:[] for x in self.geno_info}
        for ID in self.sample_IDs:
            # handle the first ID
            query = 'SELECT variant_id, GT {2} FROM __fromGeno.genotype_{0} WHERE variant_id IN (SELECT variant_id FROM __asso_tmp WHERE {1});'\
                .format(ID, where_clause, ' '.join([', ' + x for x in self.geno_info]))
            try:
                cur.execute(query, group)
            except Exception as e:
                raise ValueError('Failed to retrieve genotype and genotype info ({}) for sample with ID {}: {}'.format(self.var_info, ID, e))
            data = {x[0]:x[1:] for x in cur.fetchall()}
            #
            # genotype belonging to the same sample name are put together 
            # 
            # handle missing values
            gtmp = [data.get(x, -9.0)[0] for x in variant_id]
            # handle -1 coding (double heterozygotes)
            gtmp = [2.0 if int(x) == -1 else x for x in gtmp]
            genotype.append(array('d', gtmp))
            #
            # handle genotype_info
            for idx, key in enumerate(self.geno_info):
                # FIXME: use -9.0 as missing data for genotype info?
                val = [data.get(x, -9.0)[idx + 1] for x in variant_id]
                geno_info[key].append(array('d', val))
        #
        missing_counts = [x.count(-9.0) for x in genotype]
        # remove individuals having many missing genotypes, or have all missing variants
        # FIXME will pass it as an input arguement later
        #toKeep = [(x<0.5*numSites) for x in missing_counts]
        toKeep = [(x<numSites) for x in missing_counts]
        numtoRemove = len(self.sample_IDs)-sum(toKeep)
        if numtoRemove > 0:
            self.logger.debug('{} out of {} samples will be removed due to missing genotypes'.format(numtoRemove, len(genotype)))
        return genotype, toKeep, var_info, geno_info

    def setGenotype(self, which, data):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
    
    def setPhenotype(self, which, data, covariates=None):
        '''Set phenotype data'''
        if len(data) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        phen = [x for idx, x in enumerate(data[0]) if which[idx]]
        if covariates:
          covt = [[x for idx, x in enumerate(y) if which[idx]] for y in covariates]
        if covariates:
          self.data.setPhenotype(phen, covt)
        else:
          self.data.setPhenotype(phen)

    def run(self):
        self.db = DatabaseEngine()
        self.db.connect(self.proj.name + '.proj', readonly=True)
        self.db.attach(self.proj.name + '_genotype.DB', '__fromGeno')
        #
        while True:
            grp = self.queue.get()
            try:
                grpname = '"' + ", ".join(map(str, grp)) + '"'
            except TypeError:
                grpname = None
            if grp is None:
                self.output.send(None)
                break
            self.logger.debug('Retrieved association unit {}'.format(grpname))
            #
            self.data = t.AssoData()
            values = list(grp)
            try:
                # select variants from each group:
                genotype, which, var_info, geno_info = self.getGenotype(grp)
                self.setGenotype(which, genotype)
                self.setPhenotype(which, self.phenotypes, self.covariates)
                for field in var_info.keys():
                    self.data.setVar('__var_' + field, var_info[field])
                for field in geno_info.keys():
                    self.data.setVar('__geno_' + field, geno_info[field])
                for test in self.tests:
                    test.setData(self.data)
                    result = test.calculate()
                    self.logger.debug('Finished association test on {}'.format(grpname))
                    values.extend(result)
            except Exception as e:
                self.logger.debug('An ERROR has occured while processing {}: {}'.format(grpname, e))
                # self.data might have been messed up, create a new one
                self.data = t.AssoData()
                # return no result for any of the tests if a test fails.
                values = []
            self.output.send(values)
        self.db.detach('__fromGeno')
        

def associate(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            try:
                asso = AssociationTestManager(proj, args.table, args.phenotypes, args.covariates,
                    args.var_info, args.geno_info, args.methods, args.unknown_args,
                    args.samples, args.group_by)
            except ValueError as e:
                sys.exit(e)
            #
            nJobs = max(min(args.jobs, len(asso.groups)), 1)
            # step 4: start all workers
            grpQueue = Queue()
            results = ResultRecorder(asso, args.to_db, args.update, proj.logger)
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
                if results.completed() > count:
                    count = results.completed()
                    prog.update(count)
                    proj.logger.debug('Processed: {}/{}'.format(count, len(asso.groups)))
                #
                if not any(proc_status):
                    # if everything is done
                    assert results.completed() == len(asso.groups)
                    break
            prog.done()
            results.done()
            # summary
            proj.logger.info('Association tests on {} groups have completed. {} failed.'.format(results.completed(), results.failed()))
            # use the result database in the project
            if args.to_db:
                proj.useAnnoDB(AnnoDB(proj, args.to_db, ['chr', 'pos'] if not args.group_by else args.group_by))
    except Exception as e:
        sys.exit(e) 

#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
def getAllTests():
    '''List all tests (all classes that subclasses of NullTest/GLMBurdenTest) in this module'''
    return sorted([(name, obj) for name, obj in globals().iteritems() \
        if type(obj) == type(NullTest) and issubclass(obj, NullTest) \
            and name not in ('NullTest', 'GLMBurdenTest')], key=lambda x: x[0])


class NullTest:
    '''A base class that defines a common interface for association tests'''
    def __init__(self, logger=None, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to 
        parse it'''
        self.logger = logger
        self.parseArgs(*method_args)
        #
        self.fields = []

    def parseArgs(self, method_args):
        # this function should never be called.
        raise SystemError('All association tests should define their own parseArgs function')
    
    def setData(self, data):
        self.data = data.clone()

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant. Will print
        data if NullTest is called'''
        self.logger.info('No action is defined for NullTest')
        #self.logger.debug('Printing out group name, phenotypes, covariates and genotypes')
        #self.logger.debug(self.group)
        #self.logger.debug(self.data.phenotype())
        #self.logger.debug(self.data.covariates())
        #self.logger.debug(self.data.raw_genotype())
        return 0

def freq(input):
    try:
        value = float(input)
        if not (value >= 0 and value <= 1):
            msg = "%r is not valid input. Valid input should fall in range [0, 1]" % input
            raise ValueError(msg)
    except Exception as e:
        raise argparse.ArgumentTypeError(e)
    return value


class GroupStat(NullTest):
    '''Calculates basic statistics for each testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        #
        NullTest.__init__(self, logger, *method_args)
        # set fields according to parameter --stat
        if 'num_variants' in self.stat:
            self.fields.append(
                Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group')
            )
        if 'sample_size' in self.stat:
            self.fields.append(
                Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size')
            )

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Group statistics calculator. This 'test'
            is a statistics calculator based on the association-test framework. It is usually used
            with other method to produce statistics for each testing group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('--stat', choices=['num_variants', 'sample_size'], nargs='+',
            help='''Statistics to calculate, which can be number of variants in this group,
                total sample size.''')  
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def setData(self, data):
        # do not clone data because this test does not change data
        self.data = data

    def calculate(self):
        res = []
        for field in self.fields:
            if field.name == 'num_variants':
                # FIXME: this is inefficient because we do not have to create this object
                # to count the number of variants
                res.append(len(self.data.raw_genotype()[0]))
            elif field.name == 'sample_size':
                # FIXME: this is inefficient because we do not have to create this object
                # to count the number of variants
                res.append(self.data.samplecounts())
        return res

class GLMBurdenTest(NullTest):
    '''Generalized Linear regression test on collapsed genotypes within an association testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='beta_x', index=None, type='FLOAT', adj=None, comment='Test statistic. In the context of regression, this is estimate of effect size for x'), 
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        if self.permutations == 0:
            self.fields.append(Field(name='wald_x', index=None, type='FLOAT', adj=None, comment='Wald statistic for x (beta_x/SE(beta_x))'))
            for i in range(ncovariates):
                self.fields.extend([Field(name='beta_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='estimate of beta for covariate {}'.format(str(i+2))),
                                    Field(name='beta_{}_pvalue'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='p-value for covariate {}'.format(str(i+2))),
                                    Field(name='wald_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='Wald statistic for covariate {}'.format(str(i+2)))])
        else:
            self.fields.append(Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated'))
        #
        # NullTest.__init__ will call parseArgs to get the parameters we need
        self.algorithm = self._determine_algorithm(ncovariates)
                
    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Generalized linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='GLMTest',
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
        parser.add_argument('--model', type=int, choices = [0,1], default=0,
            help='''Fit linear regression model ("0") or logistic regression model ("1").
            Default set to 0''')
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

    def _determine_algorithm(self, ncovariates):
        if ncovariates > 0:
            a_regression = t.MultipleRegression(self.permutations == 0, self.model)
        else:
            a_regression = t.SimpleLinearRegression() if self.model == 0 else t.SimpleLogisticRegression()
        a_scoregene = t.BinToX() if self.use_indicator else t.SumToX()
        # data pre-processing
        # FIXME have to implement WeightByCtrlMaf()
        if self.weight_by_maf and not self.use_indicator:
            algorithm = t.AssoAlgorithm([
                # calculate sample MAF
                t.SetMaf(),
                # filter out variants having MAF > mafupper or MAF <= maflower
                t.SetSites(self.mafupper, self.maflower),
                # weight genotype codings by w(MAF)
                t.WeightByAllMaf()
                ])
        else:
            algorithm = t.AssoAlgorithm([
                # calculate sample MAF
                t.SetMaf(),
                # filter out variants having MAF > mafupper or MAF <= maflower
                t.SetSites(self.mafupper, self.maflower)
                ])
        # association testing using analytic p-value
        if self.permutations == 0:
            algorithm.extend([
                # calculate genotype score for a set of variants
                a_scoregene,
                # fit regression model
                a_regression,
                # evaluate p-value for the Wald's statistic
                t.StudentPval(self.alternative)
                ])
        # association testing using permutation-based p-value
        else:
            if not self.variable_thresholds:
                a_permutationtest = t.FixedPermutator(
                        self.permute_by.upper(), 
                        self.alternative, 
                        self.permutations, 
                        self.adaptive, 
                        [a_regression]
                        )
                algorithm.extend([
                        a_scoregene, 
                        a_permutationtest
                        ])
            else:
                a_permutationtest = t.VariablePermutator(
                        self.permute_by.upper(), 
                        self.alternative, 
                        self.permutations, 
                        self.adaptive, 
                        [a_scoregene, a_regression]
                        )
                algorithm.append(a_permutationtest)
        return algorithm

    def calculate(self):
        if self.model == 1:
            self.data.countCaseCtrl()
        self.algorithm.apply(self.data)
        # get results
        pvalues = self.data.pvalue()
        regstats = self.data.statistic()
        regse = self.data.se()
        res = [self.data.samplecounts()]
        for (x, y, z) in zip(regstats, pvalues, regse):
            res.append(x)
            res.append(y)
            if self.permutations == 0:
                # Wald statistic
                res.append(x/z)
            else:
                # actual number of permutations
                if math.isnan(z): res.append(z)
                else: res.append(int(z))                
        return res


#
# Derived statistical Association tests. 
#

# quantitative traits
class LinRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for quantitative traits'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)
    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LinRegBurden',
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
        #
        # We add the fixed parameter here ...
        # 
        self.model = 0

class CollapseQt(GLMBurdenTest):
    '''Collapsing method for quantitative traits, Li & Leal 2008'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for quantitative traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 0.01''')  
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = False
        self.model = 0

class BurdenQt(GLMBurdenTest):
    '''Burden test for quantitative traits, Morris & Zeggini 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for quantitative traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 0.01''')  
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = False
        self.model = 0
        
class WeightedSumQt(GLMBurdenTest):
    '''Weighted sum statistic for quantitative traits, in the spirit of Madsen & Browning 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum statistic for quantitative traits (in the spirit of Madsen & Browning 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants will be weighted by 1/sqrt(P*(1-P)) and the weighted codings will be summed
            up as one regressor''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 0.01''')  
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = True
        self.model = 0

class VariableThresholdsQt(GLMBurdenTest):
    '''Variable thresholds method for quantitative traits, in the spirit of Price et al 2010'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for quantitative traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
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
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.use_indicator=False
        self.weight_by_maf = False
        self.model = 0

# binary traits
class LogitRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for binary traits'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Logistic regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LogitRegBurden',
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
        #
        # We add the fixed parameter here ...
        #
        self.model = 1

class CollapseBt(GLMBurdenTest):
    '''Collapsing method for binary traits, Li & Leal 2008'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for binary traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = False
        self.model = 1

class BurdenBt(GLMBurdenTest):
    '''Burden test for binary traits, Morris & Zeggini 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for binary traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 0.01''')  
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = False
        self.model = 1
        
class WeightedSumBt(GLMBurdenTest):
    '''Weighted sum statistic for binary traits, in the spirit of Madsen & Browning 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum statistic for binary traits (in the spirit of Madsen & Browning 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants will be weighted by 1/sqrt(P*(1-P)) and the weighted codings will be summed
            up as one regressor''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1 
            will be included in analysis. Default set to 0.01''')  
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        # 
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight_by_maf = True
        self.model = 1

class VariableThresholdsBt(GLMBurdenTest):
    '''Variable thresholds method for binary traits, in the spirit of Price et al 2010'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for binary traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
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
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))        
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.use_indicator=False
        self.weight_by_maf = False
        self.model = 1
