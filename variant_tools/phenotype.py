#!/usr/bin/env python
#
# $File: phenotype.py $
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
import threading
import Queue
import time
import re
from collections import defaultdict
from .project import Project
from .utils import DatabaseEngine, ProgressBar, typeOfValues, SQL_KEYWORDS

class GenotypeStatStatus:
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

class GenotypeStatCalculator(threading.Thread):
    def __init__(self, dbName, stat, idQueue, status, genotypes, logger):
        '''Use sql to process sample passed from queue, set results in status'''
        self.dbName = dbName
        # query, where, idx
        self.stat = []
        for field, item in stat:
            if item == '#(GT)':
                self.stat.append(('count(*)', None))
            elif item == '#(alt)':
                self.stat.append(('sum(abs(GT))', None))
            elif item == '#(hom)':
                self.stat.append(('count(*)', 'GT=2'))
            elif item == '#(het)':
                self.stat.append(('count(*)', 'GT=1'))
            elif item == '#(other)':
                self.stat.append(('count(*)', 'GT=-1'))
            else:
                self.stat.append((item, None))
        self.queue = idQueue
        self.status = status
        self.genotypes = genotypes
        self.logger = logger
        threading.Thread.__init__(self, name='Calculate genotype statistics')

    def run(self):
        db = DatabaseEngine()
        db.connect(self.dbName, readonly=True)
        projName = self.dbName.split('_')[0]
        db.attach(projName +'.proj', 'proj')
        cur = db.cursor()
        while True:
            ID = self.queue.get()
            # finish everything
            if ID is None:
                self.queue.task_done()
                break
            # if everything can be executed in a single query
            if all([x[1] is None for x in self.stat]):
                # run query
                res = [None] * len(self.stat)
                # regular stat
                query = 'SELECT {} FROM genotype_{} {};'\
                    .format(', '.join([x[0] for x in self.stat]), ID,
                        'WHERE {}'.format(self.genotypes) if self.genotypes.strip() else '')
                self.logger.debug(query)
                try:
                    cur.execute(query)
                    res = cur.fetchone()
                except Exception as e:
                    # some field might not exist, so we will have to run one by one
                    for idx, (expr, where) in enumerate(self.stat):
                        query = 'SELECT {} FROM genotype_{} {};'\
                            .format(expr, ID,
                                'WHERE {}'.format(self.genotypes) if self.genotypes.strip() else '')
                        self.logger.debug(query)
                        try:
                            cur.execute(query)
                            v = cur.fetchone()
                            if v is not None:
                                res[idx] = v[0]
                        except Exception as e:
                            self.logger.debug('Failed to evalulate {}: {}. Setting field to NULL.'.format(expr, e))
            else:
                res = [None] * len(self.stat)
                for idx, (expr, where) in enumerate(self.stat):
                    where_clause = 'WHERE {}'.format(self.genotypes) if self.genotypes.strip() else ''
                    if where:
                        if where_clause:
                            where_clause += ' AND ({})'.format(where)
                        else:
                            where_clause = 'WHERE {}'.format(where)
                    query = 'SELECT {} FROM genotype_{} {};'\
                        .format(expr, ID, where_clause)
                    self.logger.debug(query)
                    try:
                        cur.execute(query)
                        v = cur.fetchone()
                        if v is not None:
                            res[idx] = v[0]
                    except Exception as e:
                        self.logger.debug('Failed to evalulate {}: {}. Setting field to NULL.'.format(expr, e))
            #
            # set result
            self.status.set(ID, res)
            self.queue.task_done()
        db.close()  

class Sample:
    def __init__(self, proj, jobs=4):
        self.proj = proj
        self.jobs = jobs
        self.logger = proj.logger
        self.db = proj.db

    def load(self, filename, allowed_fields, samples):
        '''Load phenotype informaiton from a file'''
        if not self.db.hasTable('sample'):
            self.logger.warning('Project does not have a sample table.')
            return
        # num sample, num new field, num update field
        count = [0, 0, 0]
        with open(filename) as input:
            headers = input.readline().rstrip().split('\t')
            if len(headers) == 0:
                raise ValueError('Empty header line. No phenotype will be imported.')
            if headers[0] == 'sample_name':
                by_sample = True
                if len(headers) == 1:
                    raise ValueError('No phenotype to be imported')
            elif len(headers) >= 2 and headers[0] == 'filename' and headers[1] == 'sample_name':
                by_sample = False
                if len(headers) == 2:
                    raise ValueError('No phenotype to be imported')
            else:
                raise ValueError('The phenotype file must start with a header line with the first '
                    'column sample_name, or first two fields being filename and sample_name.')
            #
            records = {}
            nCol = len(headers)
            for idx, line in enumerate(input.readlines()):
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = [x.strip() for x in line.split('\t')]
                if len(fields) != nCol or ('' in fields):
                    raise ValueError('Invalid phenotype file: number of fields mismatch at line {}. Please use \'None\' for missing values.'.format(idx+2))
                #
                if by_sample:
                    key = (None if fields[0] == 'None' else fields[0],)
                    if key in records:
                        raise ValueError('Duplicate sample name ({}). Only the last record will be used'.format(key))
                    records[key] = fields[1:]
                else:
                    key = (fields[0], None if fields[1] == 'None' else fields[1])
                    if key in records:
                        raise ValueError('Duplicate filename and sample name ({},{}). Only the last record will be used'.format(key[0], key[1]))
                    records[key] = fields[2:]
            #
            new_fields = headers[(1 if by_sample else 2):]
            if allowed_fields:
                for f in allowed_fields:
                    if f.lower() not in [x.lower() for x in new_fields]:
                        raise ValueError('Field {} is not in specified input file {}'.format(f, filename))
        # 
        # get allowed samples
        cur = self.db.cursor()
        allowed_samples = self.proj.selectSampleByPhenotype(samples)
        if not allowed_samples:
            raise ValueError('No sample is selected using condition "{}"'.format(samples))
        #
        # get existing fields
        cur_fields = self.db.getHeaders('sample')[3:]
        # handle each field one by one
        for idx, field in enumerate(new_fields):
            if allowed_fields and field.lower() not in [x.lower() for x in allowed_fields]:
                self.logger.debug('Ignoring field {}'.format(field))
                continue
            # if adding a new field
            if field.lower() not in [x.lower() for x in cur_fields]:
                self.proj.checkFieldName(field, exclude='sample')
                fldtype = typeOfValues([x[idx] for x in records.values()])
                self.logger.info('Adding field {}'.format(field))
                self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(field, fldtype))
                count[1] += 1  # new
            else:
                count[2] += 1  # updated
            for key, rec in records.iteritems():
                # get matching sample
                if by_sample:
                    cur.execute('SELECT sample.sample_id FROM sample WHERE sample_name = {}'.format(self.db.PH), key)
                    ids = [x[0] for x in cur.fetchall()]
                    if len(ids) == 0:
                        self.logger.warning('Sample name {} does not match any sample'.format(key[0]))
                        continue
                    for id in [x for x in ids if x in allowed_samples]:
                        count[0] += 1
                        cur.execute('UPDATE sample SET {0}={1} WHERE sample_id={1};'.format(field, self.db.PH), [rec[idx], id])
                else:
                    cur.execute('SELECT sample.sample_id FROM sample LEFT JOIN filename ON sample.file_id = filename.file_id WHERE filename.filename = {0} AND sample.sample_name = {0}'.format(self.db.PH), key)
                    ids = [x[0] for x in cur.fetchall()]
                    if len(ids) == 0:
                        self.logger.warning('Filename {} and sample name {} does not match any sample'.format(key[0], key[1]))
                        continue
                    if len(ids) != 1:
                        raise ValueError('Filename and sample should unqiuely determine a sample')
                    for id in [x for x in ids if x in allowed_samples]:
                        count[0] += 1
                        cur.execute('UPDATE sample SET {0}={1} WHERE sample_id={1};'.format(field, self.db.PH), [rec[idx], id])
        self.logger.info('{} field ({} new, {} existing) phenotypes of {} samples are updated.'.format(
            count[1]+count[2], count[1], count[2], count[0]/(count[1] + count[2])))
        self.db.commit()

    def setPhenotype(self, field, expression, samples):
        '''Add a field using expression calculated from sample variant table'''
        IDs = self.proj.selectSampleByPhenotype(samples)
        if not IDs:
            raise ValueError('No sample is selected using condition "{}"'.format(samples))
        #
        count = [0, 0, 0]
        cur = self.db.cursor()
        cur.execute('SELECT {} FROM sample;'.format(expression))
        fldType = None
        for rec in cur:
            if fldType is None:
                fldType = type(rec[0])
                continue
            elif rec[0] is None: # missing
                continue
            if type(rec[0]) != fldType:
                if type(rec[0]) is float and fldType is int:
                    fltType = float
                else:
                    raise ValueError('Inconsistent type returned from different samples')
        if expression != 'NULL' and fldType is None:
            raise ValueError('Cannot determine the type of the expression')
        # if adding a new field
        cur_fields = self.db.getHeaders('sample')[3:]
        if field.lower() not in [x.lower() for x in cur_fields]:
            if field.upper in SQL_KEYWORDS:
                raise ValueError("Phenotype name '{}' is not allowed because it is a reserved word.".format(x))
            self.logger.info('Adding field {}'.format(field))
            self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(field,
                {int: 'INT',
                 float: 'FLOAT',
                 str: 'VARCHAR(255)',
                 unicode: 'VARCHAR(255)',
                 None: 'FLOAT'}[fldType]))
            count[1] += 1  # new
        else:
            # FIXME: check the case for type mismatch
            count[2] += 1  # updated
        #
        cur = self.db.cursor()
        for ID in IDs:
            cur.execute('UPDATE sample SET {0}={1} WHERE sample_id = {2}'.format(field, 
                None if expression == 'NULL' else expression, self.db.PH), (ID,))
            count[0] += 1
        self.logger.info('{} values of {} phenotypes ({} new, {} existing) of {} samples are updated.'.format(
            count[0], count[1]+count[2], count[1], count[2], len(IDs)))
        self.db.commit()

    def fromSampleStat(self, stat, genotypes, samples):
        '''Add a field using expression calculated from sample variant table'''
        IDs = self.proj.selectSampleByPhenotype(samples)
        if not IDs:
            raise ValueError('No sample is selected using condition "{}"'.format(samples))
        #
        # at least one, at most number of IDs
        nJobs = max(min(self.jobs, len(IDs)), 1)
        # start all workers
        idQueue = Queue.Queue()
        status = GenotypeStatStatus()
        for j in range(nJobs):
            GenotypeStatCalculator('{}_genotype.DB'.format(self.proj.name),
                stat, idQueue, status, genotypes, self.logger).start()
        #
        # put all jobs to queue, the workers will work on them
        for ID in IDs:
            idQueue.put(ID)
        #
        count = 0
        prog = ProgressBar('Calculating phenotype', len(IDs))
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
        # submit all results, these should be quick so no progress bar is used
        count = [0, 0, 0]
        cur_fields = self.db.getHeaders('sample')[3:]
        new_field = {}
        for field in [x[0] for x in stat]:
            if field.lower() not in [x.lower() for x in cur_fields]:
                if field.upper in SQL_KEYWORDS:
                    raise ValueError("Phenotype name '{}' is not allowed because it is a reserved word.".format(field))
                new_field[field] = True
            else:
                new_field[field] = False
                count[2] += 1  # updated
        cur = self.db.cursor()
        for ID in IDs:
            res = status.get(ID)
            for idx, (field, expr) in enumerate(stat):
                if new_field[field]:
                    self.logger.debug('Adding field {}'.format(field))
                    # determine the type of value
                    self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(field,
                        typeOfValues([str(status.get(x)[idx]) for x in IDs])))
                    new_field[field] = False
                    count[1] += 1  # new
                cur.execute('UPDATE sample SET {0}={1} WHERE sample_id = {1}'.format(field, self.db.PH), [res[idx], ID])
                count[0] += 1
        # report result
        self.logger.info('{} values of {} phenotypes ({} new, {} existing) of {} samples are updated.'.format(
            count[0], count[1]+count[2], count[1], count[2], len(IDs)))
        self.db.commit()

    def output(self, fields, samples):
        # output
        query = 'SELECT {} FROM sample LEFT JOIN filename ON sample.file_id = filename.file_id {}'.format(
            ','.join(fields), '' if not samples else 'WHERE ' + samples)
        self.logger.debug(query)
        cur = self.db.cursor()
        cur.execute(query)
        for rec in cur:
            print '\t'.join([str(x) for x in rec])
                
def phenotypeArguments(parser):
    '''Action that can be performed by this script'''
    parser.add_argument('-f', '--from_file', metavar='INPUT_FILE', nargs='*',
        help='''Import phenotype from a tab delimited file. The file should have
            a header, with either 'sample_name' as the first column, or 'filename'
            and 'sample_name' as the first two columns. In the former case, samples
            with the same 'sample_name' will share the imported phenotypes. If 
            a list of phenotypes (columns of the file) is specified after filename,
            only the specified phenotypes will be imported. Parameter --samples
            could be used to limit the samples for which phenotypes are imported.'''),
    parser.add_argument('--set', nargs='*', metavar='EXPRESSION', default=[],
        help='''Set a phenotype to a constant (e.g. --set aff=1), or an expression
            using other existing phenotypes (e.g. --set ratio_qt=high_qt/all_qt (the ratio
            of the number of high quality variants to the number of all variants, where
            high_qt and all_qt are obtained from sample statistics using parameter
            --from_stat). Parameter --samples could be used to limit the samples for
            which genotypes will be set.'''),
    parser.add_argument('--from_stat', nargs='*', metavar='EXPRESSION', default=[],
        help='''Set a phenotype to a summary statistics of a genotype field. For 
            example, "num=count(*)" sets phenotype num to be the number of genotypes
            of a sample, "GD=avg(DP)" sets phenotype DP to be the average depth (if
            DP is one of the genotype fields) of the sample. Multiple fields (e.g.
            '--from_stat "num=count(*)" "GD=avg(DP)"') are also allowed. In addition to
            standard SQL aggregation functions, variant tools supports special functions
            #(GT), #(alt), #(hom), #(het) and #(other), which calculates the number of
            genotypes (the same as count(*)), alternative alleles, homozygotes, 
            heterozygotes, and genotypes with two different alternative alleles.
            Parameters --genotypes and --samples could be used to limit the genotypes
            to be considered and the samples for which genotypes will be set.'''),
    parser.add_argument('--output', nargs='*', metavar='EXPRESSION', default=[],
        help='''A list of phenotype to be outputted. SQL-compatible expressions or
            functions such as "DP/DP_all" and "avg(DP)" are also allowed'''),
    parser.add_argument('-j', '--jobs', metavar='N', default=4, type=int,
        help='''Allow at most N concurrent jobs to obtain sample statistics for
            parameter --from_stat.''')
    parser.add_argument('-g', '--genotypes', nargs='*', metavar='COND', default=[],
        help='''Limit the operation to genotypes that match specified conditions.
            Use 'vtools show genotypes' to list usable fields for each sample.'''),
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Update phenotype for samples that match specified conditions.
            Use 'vtools show samples' to list usable fields in the sample table.''')

def phenotype(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            p = Sample(proj, args.jobs)
            if args.from_file:
                filename = args.from_file[0]
                fields = args.from_file[1:]
                p.load(filename, fields, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.set:
                for item in args.set:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError('Invalid parameter {}, which should have format field=expr_of_phenotype: {}'.format(item, e))
                    p.setPhenotype(field, expr, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.from_stat:
                stat = []
                for item in args.from_stat:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError('Invalid parameter {}, which should have format field=expr_of_field: {}'.format(item, e))
                    stat.append((field, expr))
                p.fromSampleStat(stat,
                        ' AND '.join(['({})'.format(x) for x in args.genotypes]),
                        ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.output:
                p.output(args.output, ' AND '.join(['({})'.format(x) for x in args.samples]))
        proj.close()
    except Exception as e:
        sys.exit(e)
                

