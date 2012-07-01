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
from .utils import ProgressBar, consolidateFieldName, DatabaseEngine, delayedAction, runOptions
from .phenotype import Sample
from .tester import *

import argparse

def istr(d):
    try:
        x = int(d)
        x = str(x)
    except:
        x = str(d)
    return x

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
            cf. "vtools show genotypes") that will be passed to statistical
            tests. Note that the fields should exist for all samples that are
            tested.''')
    parser.add_argument('--moi', type=str, choices = ['additive','dominant', 'recessive'],
            default='additive',
            help='''Mode of inheritance. Will code genotypes as 0/1/2/NA for additive mode, 
            0/1/NA for dominant or recessive model.
            Default set to quantitative''')
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
    parser.add_argument('--missing_filter', metavar='proportion', type=freq, default=1.0,
           help='''When -g is specified, will remove samples having missing genotypes
           exceeding the given proportion within a group''' )
    parser.add_argument('--to_db', metavar='annoDB',
        help='''Name of a database to which results from association tests will be written''')
    parser.add_argument('--update', action='store_true',
        help='''When --to_db is specified, allow updating existing fields in the result database''')
    parser.add_argument('-j', '--jobs', metavar='N', default=1, type=int,
        help='''Number of processes to carry out association tests.''')

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
    def __init__(self, proj, table, phenotypes, covariates, var_info, geno_info, moi, methods,
        unknown_args, samples, group_by, missing_filter):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        self.var_info = var_info
        self.geno_info = geno_info
        self.missing_filter = missing_filter
        self.moi = moi
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
                    if test.trait_type == 'disease':
                        raise ValueError("{0} cannot handle non-binary phenotype".format(test.__class__.__name__))
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
                    self.logger.warning('Association test {} has invalid or empty fields. No result will be generated.'.format(name))
                tests.append(method)
            except NameError as e:
                self.logger.debug(e)
                raise ValueError('Failed to load association test {0}: {1}. Please use command "vtools show tests" to list usable tests'.format(name, e))
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
                raise ValueError('Invalid (non-numeric) coding in phenotype/covariates values: '
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
            CREATE {} TABLE __asso_tmp (
              variant_id INT NOT NULL,
              {} {});
              '''.format('TEMPORARY' if runOptions.associate_genotype_cache_size > 0 else '',
              ','.join(['{} {}'.format(x,y) for x,y in zip(field_names, field_types)]),
              ''.join([', {} FLOAT'.format(x) for x in self.var_info])))
        #
        # select variant_id and groups for association testing
        group_fields, fields = consolidateFieldName(self.proj, self.table, ','.join(group_by + self.var_info))
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
        query = 'INSERT INTO __asso_tmp SELECT DISTINCT {}.variant_id, {} FROM {} {};'.format(
            self.table, group_fields,
            self.from_clause, self.where_clause)
        s = delayedAction(self.logger.info, "Grouping variants by {}, please be patient ...".format(', '.join(group_by)))
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

    def getVarInfo(self, grp):
        var_info = {x:[] for x in self.var_info}
        if not self.var_info:
            return var_info
        #
        where_clause = ' AND '.join(['{0}={1}'.format(x, self.db.PH) for x in self.group_names])
        query = 'SELECT {0} FROM __asso_tmp WHERE {2})'.format(
            ', '.join(self.var_info), from_clause, where_clause)
        #self.logger.debug('Running query: {}'.format(query))
        cur.execute(query, group)
        #
        data = {x[0]:x[1:] for x in cur.fetchall()}
        for idx, key in enumerate(self.var_info):
            var_info[key] = [data[x][idx] for x in variant_id]
        return var_info


class GenotypeGrabber:
    def __init__(self, param):
        #
        self.proj = param.proj
        self.table = param.table
        self.sample_IDs = param.sample_IDs
        self.phenotypes = param.phenotypes
        self.covariates = param.covariates
        self.var_info = param.var_info
        self.geno_info = param.geno_info
        self.moi = param.moi
        self.tests = param.tests
        self.group_names = param.group_names
        self.missing_filter = param.missing_filter
        #
        self.logger = param.proj.logger
        if runOptions.associate_genotype_cache_size == 0:
            self.db = DatabaseEngine()
            self.db.connect(param.proj.name + '.proj', readonly=True)
            self.db.attach(param.proj.name + '_genotype.DB', '__fromGeno')
        else:
            self.db = self.proj.db
            self.proj.logger.debug('Setting PRAGMA __fromGeno.cache_size=-{}'.format(runOptions.associate_genotype_cache_size))
            self.db.execute('PRAGMA __fromGeno.cache_size=-{}'.format(runOptions.associate_genotype_cache_size))
        
    def __del__(self):
        self.db.detach('__fromGeno')
        self.db.close()
        
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
        geno_info = {x:[] for x in self.geno_info}
        for ID in self.sample_IDs:
            query = 'SELECT variant_id, GT {2} FROM __fromGeno.genotype_{0} WHERE variant_id IN (SELECT variant_id FROM __asso_tmp WHERE {1});'\
                .format(ID,  where_clause, ' '.join([', ' + x for x in self.geno_info]))
            try:
                cur.execute(query, group) # , variant_id)
            except Exception as e:
                raise ValueError('Failed to retrieve genotype and genotype info ({0}) for sample with ID {1}: {2}'.format(self.var_info, ID, e))
            data = {x[0]:x[1:] for x in cur.fetchall()}
            #
            # genotype belonging to the same sample name are put together
            #
            # handle missing values
            gtmp = [data.get(x, [float('NaN')]*(len(self.geno_info)+1)) for x in variant_id]
            # handle -1 coding (double heterozygotes)
            genotype.append(array('d', [2.0 if x[0] == -1.0 else x[0] for x in gtmp]))
            #
            # handle genotype_info
            for idx, key in enumerate(self.geno_info):
                geno_info[key].append(array('d', [x[idx+1] for x in gtmp]))
        #
        missing_counts = [sum(list(map(math.isnan, x))) for x in genotype]
        # remove individuals having many missing genotypes, or have all missing variants
        toKeep = [(x < (self.missing_filter * numSites)) for x in missing_counts]
        numtoRemove = len(self.sample_IDs) - sum(toKeep)
        if numtoRemove > 0:
            self.logger.debug('{} out of {} samples will be removed due to missing genotypes'.format(numtoRemove, len(genotype)))
        return genotype, toKeep, geno_info


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
    def __init__(self, param, grpQueue, resQueue):
        Process.__init__(self, name='Phenotype association analysis for a group of variants')
        self.param = param
#        self.proj = param.proj
#        self.table = param.table
#        self.sample_IDs = param.sample_IDs
        self.phenotypes = param.phenotypes
        self.covariates = param.covariates
#        self.var_info = param.var_info
#        self.geno_info = param.geno_info
        self.moi = param.moi
        self.tests = param.tests
        self.num_extern_tests = sum([isinstance(x, ExternTest) for x in self.tests])
#        self.group_names = param.group_names
#        self.missing_filter = param.missing_filter
        self.queue = grpQueue
        self.resQueue = resQueue
        self.logger = param.proj.logger
        self.db = None

    def setGenotype(self, which, data, info):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
        self.data.setMOI(self.moi)
        for field in info.keys():
            self.data.setVar('__geno_' + field, [x for idx, x in enumerate(info[field]) if which[idx]])

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

    def setVarInfo(self, data):
        for field in data.keys():
            self.data.setVar('__var_' + field, data[field])

    def setPyData(self, which, geno, pheno, covar, var_info, geno_Info, missing_code=None, grpname=None):
        '''set all data to a python dictionary in str format'''
        self.pydata['name'] = grpname
        if len(pheno) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        if missing_code:
            self.pydata['genotype'] = [map(istr, [missing_code if math.isnan(e) else e for e in x]) for idx, x in enumerate(geno) if which[idx]]
        else:
            self.pydata['genotype'] = [map(istr, x) for idx, x in enumerate(geno) if which[idx]]
        self.pydata['phenotype'] = map(str, [x for idx, x in enumerate(pheno[0]) if which[idx]])
        if len(covar) > 0:
            # skip the first covariate, a vector of '1''s
            self.pydata['covariates'] = [map(str, [x for idx, x in enumerate(y) if which[idx]]) for y in covar[1:]]
        #FIXME: will not use var_info and geno_info for now
        if len(self.pydata['genotype']) == 0 or len(self.pydata['phenotype']) == 0 or len(self.pydata['genotype'][0]) == 0:
            raise ValueError("No input data")

    def run(self):
        # if genotypes are not cached, each worker will grab genotype from the database by itself
        if runOptions.associate_genotype_cache_size == 0:
            gg = GenotypeGrabber(self.param)
        while True:
            # if cached, get genotype from the main process
            if runOptions.associate_genotype_cache_size > 0:
                grp, var_info, (genotype, which, geno_info) = self.queue.get()
            else:
                # otherwise, only the group
                grp, var_info = self.queue.get()
            #
            try:
                grpname = ", ".join(map(str, grp))
            except TypeError:
                grpname = None
            if grp is None:
                break
            self.logger.debug('Retrieved association unit {}'.format(repr(grpname)))
            #
            self.data = t.AssoData()
            self.pydata = {}
            values = list(grp)
            try:
                if runOptions.associate_genotype_cache_size == 0:
                    # select variants from each group:
                    genotype, which, geno_info = gg.getGenotype(grp)
                # set C++ data object
                if (len(self.tests) - self.num_extern_tests) > 0:
                    self.setGenotype(which, genotype, geno_info)
                    self.setPhenotype(which, self.phenotypes, self.covariates)
                    self.setVarInfo(var_info)
                # set Python data object, for external tests
                if self.num_extern_tests:
                    self.setPyData(which, genotype, self.phenotypes, self.covariates, var_info, geno_info, 'NA', grpname)

                # association tests
                for test in self.tests:
                    test.setData(self.data, self.pydata)
                    result = test.calculate()
                    self.logger.debug('Finished association test on {}'.format(repr(grpname)))
                    values.extend(result)
            except Exception as e:
                self.logger.debug('An ERROR has occurred while processing {}: {}'.format(repr(grpname), e))
                # self.data might have been messed up, create a new one
                self.data = t.AssoData()
                # return no result for any of the tests if a test fails.
                values = []
            self.resQueue.put(values)


def associate(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            try:
                asso = AssociationTestManager(proj, args.table, args.phenotypes, args.covariates,
                    args.var_info, args.geno_info, args.moi, args.methods, args.unknown_args,
                    args.samples, args.group_by, args.missing_filter)
            except ValueError as e:
                sys.exit(e)
            #
            nJobs = max(min(args.jobs, len(asso.groups)), 1)
            # step 4: start all workers
            grpQueue = Queue()
            # the result queue is used by workers to return results
            resQueue = Queue()
            results = ResultRecorder(asso, args.to_db, args.update, proj.logger)
            for j in range(nJobs):
                AssoTestsWorker(asso, grpQueue, resQueue).start()
            #
            # put all jobs to queue, the workers will work on them
            if runOptions.associate_genotype_cache_size > 0:
                geno = GenotypeGrabber(asso)
                prog = ProgressBar('Getting genotypes', len(asso.groups))
                for count,grp in enumerate(asso.groups):
                    var_info = asso.getVarInfo(grp)
                    grpQueue.put((grp, var_info, geno.getGenotype(grp)))
                    prog.update(count)
                for j in range(nJobs):
                    grpQueue.put(None)
                prog.done()
            else:
                for grp in asso.groups:
                    var_info = asso.getVarInfo(grp)
                    grpQueue.put((grp, var_info))
                # the worker will stop once all jobs are finished
                for j in range(nJobs):
                    grpQueue.put(None)
            #
            count = 0
            # get initial completed and failed
            prog = ProgressBar('Testing for association', len(asso.groups), results.completed(), results.failed())
            while True:
                # if everything is done
                if count >= len(asso.groups):
                    break
                # not done? wait from the queue and write to the result recorder
                res = resQueue.get()
                results.record(res)
                # update progress bar
                count = results.completed()
                prog.update(count, results.failed())
                proj.logger.debug('Processed: {}/{}'.format(count, len(asso.groups)))
            # finished
            prog.done()
            results.done()
            # summary
            proj.logger.info('Association tests on {} groups have completed. {} failed.'.format(results.completed(), results.failed()))
            # use the result database in the project
            if args.to_db:
                proj.useAnnoDB(AnnoDB(proj, args.to_db, ['chr', 'pos'] if not args.group_by else args.group_by))
    except Exception as e:
        sys.exit(e)
