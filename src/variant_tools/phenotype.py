#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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

import csv
import os
import queue
import sys
import tempfile
import threading
import time
from collections import defaultdict
from multiprocessing import Manager, Process
from multiprocessing import Queue as mpQueue


from .geno_store import GenoStore
from .project import Project
from .utils import (SQL_KEYWORDS, DatabaseEngine, ProgressBar, env,
                    typeOfValues, validFieldName)


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

    def __init__(self, dbName, stat, idQueue, status, genotypes):
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
            elif item == '#(wtGT)':
                self.stat.append(('count(*)', 'GT=0'))
            elif item == '#(mutGT)':
                self.stat.append(('count(*)', 'GT!=0'))
            elif item.startswith('#'):
                raise ValueError(
                    '{} is not a valid special function (only #(GT), #(wtGT), #(mutGT), #(alt), #(hom), #(het), and #(other) are allowed).'
                    .format(item))
            else:
                self.stat.append((item, None))
        self.queue = idQueue
        self.status = status
        self.genotypes = genotypes
        threading.Thread.__init__(self, name='Calculate genotype statistics')

    def run(self):
        db = DatabaseEngine()
        db.connect(self.dbName, readonly=True)
        # dbName can be a_1_genotype so we have to split from the end
        projName = self.dbName.rsplit('_', 1)[0]
        db.attach(projName + '.proj', 'proj')
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
                env.logger.debug(query)
                try:
                    cur.execute(query)
                    res = cur.fetchone()
                except Exception:
                    # some field might not exist, so we will have to run one by one
                    for idx, (expr, where) in enumerate(self.stat):
                        query = 'SELECT {} FROM genotype_{} {};'\
                            .format(expr, ID,
                                'WHERE {}'.format(self.genotypes) if self.genotypes.strip() else '')
                        env.logger.debug(query)
                        try:
                            cur.execute(query)
                            v = cur.fetchone()
                            if v is not None:
                                res[idx] = v[0]
                        except Exception as e:
                            env.logger.debug(
                                'Failed to evalulate {}: {}. Setting field to NULL.'
                                .format(expr, e))
            else:
                res = [None] * len(self.stat)
                for idx, (expr, where) in enumerate(self.stat):
                    where_clause = 'WHERE {}'.format(
                        self.genotypes) if self.genotypes.strip() else ''
                    if where:
                        if where_clause:
                            where_clause += ' AND ({})'.format(where)
                        else:
                            where_clause = 'WHERE {}'.format(where)
                    query = 'SELECT {} FROM genotype_{} {};'\
                        .format(expr, ID, where_clause)
                    env.logger.debug(query)
                    try:
                        cur.execute(query)
                        v = cur.fetchone()
                        if v is not None:
                            res[idx] = v[0]
                    except Exception as e:
                        env.logger.debug(
                            'Failed to evalulate {}: {}. Setting field to NULL.'
                            .format(expr, e))
            #
            # set result
            self.status.set(ID, res)
            self.queue.task_done()
        db.close()


class GenotypeStatCalculator_HDF5(Process):

    def __init__(self, proj, stat, idQueue, status, genotypes):
        '''Use sql to process sample passed from queue, set results in status'''
        self.proj = proj
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
            elif item == '#(wtGT)':
                self.stat.append(('count(*)', 'GT=0'))
            elif item == '#(mutGT)':
                self.stat.append(('count(*)', 'GT!=0'))
            elif item.startswith('#'):
                raise ValueError(
                    '{} is not a valid special function (only #(GT), #(wtGT), #(mutGT), #(alt), #(hom), #(het), and #(other) are allowed).'
                    .format(item))
            else:
                self.stat.append((item, None))
        self.queue = idQueue
        self.status = status
        self.genotypes = genotypes

        # threading.Thread.__init__(self, name='Calculate genotype statistics')
        Process.__init__(self, name='Calculate genotype statistics')

    def run(self):
        # db = DatabaseEngine()
        # db.connect(self.proj.proj_file, readonly=True)
        # cur = db.cursor()
        while True:
            ID = self.queue.get()
            # finish everything
            if ID is None:
                break
            res = [None] * len(self.stat)
            try:
                for idx, (expr, where) in enumerate(self.stat):
                    store = GenoStore(self.proj)
                    if expr == "count(*)":
                        res[idx] = store.num_genotypes(ID, where,
                                                       self.genotypes)
                    elif expr == "sum(abs(GT))":
                        res[idx] = store.sum_genotypes(ID, where,
                                                       self.genotypes)
                    elif expr.startswith("avg") or expr.startswith(
                            "min") or expr.startswith("max"):
                        res[idx] = store.num_genoinfo(ID, expr, where)
            except Exception as e:
                print(e)
                env.logger.debug(
                    'Failed to evalulate {}: {}. Setting field to NULL.'.format(
                        expr, e))
            self.status[ID] = res
        # db.close()


class Sample:

    def __init__(self, proj, jobs=4):
        self.proj = proj
        self.jobs = jobs
        self.db = proj.db
        self.proj.db = proj.db
        self.proj.db.connect(self.proj.proj_file)

    def load(self, filename, header, allowed_fields, samples, na_str):
        '''Load phenotype information from a file'''
        if not self.db.hasTable('sample'):
            env.logger.warning('Project does not have a sample table.')
            return
        #
        from_stdin = filename == '-'
        if from_stdin:
            with tempfile.NamedTemporaryFile(delete=False) as tfile:
                tfile.write(sys.stdin.read())
                filename = tfile.name
        #
        if header is not None:
            if len(header) == 0:
                raise ValueError(
                    'A list of headers is required if parameter --header is specified'
                )
            elif header == ['-']:
                if from_stdin:
                    raise ValueError(
                        'Cannot read header from standard input because '
                        'input file is also read from standard input')
                env.logger.info('Reading header from standard input')
                headers = sys.stdin.read().rstrip().split()
            else:
                headers = header
        # num sample, num new field, num update field
        count = [0, 0, 0]
        # grab 4 lines from text for csv to sniff the delimiter
        example_lines = []
        for l in open(filename, 'rU'):
            if not l.startswith('#'):
                example_lines.append(l)
            if len(example_lines) >= 4:
                break
        csv_dialect = csv.Sniffer().sniff(''.join(example_lines))
        #
        with open(filename, 'rU') as input:
            reader = csv.reader(input, dialect=csv_dialect)
            if header is None:
                headers = next(reader)
            if len(set(headers)) != len(headers):
                for x in set(headers):
                    headers.remove(x)
                raise ValueError(
                    'Duplicated header names in input file {}: {}'.format(
                        filename, ', '.join(headers)))
            # determine fields to identify sample
            sample_idx = []
            if 'filename' in headers:
                # by filename
                sample_idx.append(headers.index('filename'))
                # sample_name should be the first, or the second if
                # filename is the first
                if headers[0] == 'filename':
                    sample_idx.append(1)
                else:
                    sample_idx.append(0)
            else:
                sample_idx.append(0)
            #
            if headers[sample_idx[-1]] != 'sample_name':
                env.logger.warning(
                    'Sample name field {} does not have header "sample_name"'
                    .format(headers[sample_idx[-1]]))
            #
            # determine phenotype fields
            phenotype_idx = [
                idx for idx, fld in enumerate(headers)
                if idx not in sample_idx and
                (not allowed_fields or fld in allowed_fields)
            ]
            #
            if not phenotype_idx:
                env.logger.error('No phenotype field to be imported')
                return
            #
            new_fields = []
            for idx in phenotype_idx:
                new_header = validFieldName(
                    headers[idx],
                    reserved=[
                        'filename', 'sample_name', 'sample_id', 'file_id'
                    ])
                if new_header != headers[idx]:
                    env.logger.warning(
                        'Phenotype "{}" is renamed to "{}".'.format(
                            headers[idx], new_header))
                new_fields.append(new_header)
            #
            records = {}
            nCol = len(headers)
            for fields in reader:
                if len(fields) != nCol:
                    env.logger.warning(
                        'Number of fields mismatch (expecting {}). Ignoring line "{}"'
                        .format(nCol, fields))
                    continue
                # ignore empty line (line with empty sample name)
                if not any([fields[x] for x in sample_idx]):
                    continue
                key = tuple([fields[x] for x in sample_idx])
                if key in records:
                    raise ValueError('Duplicate sample name ({}).'.format(
                        ', '.join(key)))
                records[key] = [fields[x] for x in phenotype_idx]
        # remove temporary file
        if from_stdin:
            os.remove(filename)
        #
        # get allowed samples
        cur = self.db.cursor()
        allowed_samples = self.proj.selectSampleByPhenotype(samples)
        if not allowed_samples:
            raise ValueError(
                'No sample is selected using condition "{}"'.format(samples))
        #
        # get existing fields
        cur_fields = self.db.getHeaders('sample')[3:]
        # handle each field one by one
        invalid_sample_names = set()
        for idx, field in enumerate(new_fields):
            # if adding a new field
            if field.lower() not in [x.lower() for x in cur_fields]:
                if field.upper() in SQL_KEYWORDS:
                    raise ValueError(
                        "Phenotype name '{}' is not allowed because it is a reserved word."
                        .format(field))
                fldtype = typeOfValues([x[idx] for x in list(records.values())])
                env.logger.info('Adding phenotype {} of type {}'.format(
                    field, fldtype))
                env.logger.debug(
                    'Executing ALTER TABLE sample ADD {} {} NULL;'.format(
                        field, fldtype))
                self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(
                    field, fldtype))
                count[1] += 1  # new
            else:
                fldtype = self.db.typeOfColumn('sample', field)
                count[2] += 1  # updated
            null_count = defaultdict(int)
            for key, rec in records.items():
                # by sample_name only
                if len(key) == 1:
                    # get matching sample
                    cur.execute(
                        'SELECT sample.sample_id FROM sample WHERE sample_name = {}'
                        .format(self.db.PH), key)
                else:  # by filename and sample_name
                    cur.execute(
                        'SELECT sample.sample_id FROM sample '
                        'LEFT JOIN filename ON sample.file_id = filename.file_id '
                        'WHERE filename.filename = {0} AND sample.sample_name = {0}'
                        .format(self.db.PH), key)
                ids = [x[0] for x in cur.fetchall()]
                if len(ids) == 0:
                    invalid_sample_names.add(key[0])
                    continue
                #
                for id in [x for x in ids if x in allowed_samples]:
                    count[0] += 1
                    if rec[idx] == na_str:
                        null_count[field] += 1
                        rec[idx] = None
                    elif fldtype.upper().startswith('INT'):
                        try:
                            int(rec[idx])
                        except:
                            env.logger.warning(
                                'Value "{}" is treated as missing in phenotype {}'
                                .format(rec[idx], field))
                            null_count[field] += 1
                            rec[idx] = None
                    elif fldtype.upper().startswith('FLOAT'):
                        try:
                            float(rec[idx])
                        except:
                            env.logger.warning(
                                'Value "{}" is treated as missing in phenotype {}'
                                .format(rec[idx], field))
                            null_count[field] += 1
                            rec[idx] = None
                    cur.execute(
                        'UPDATE sample SET {0}={1} WHERE sample_id={1};'.format(
                            field, self.db.PH), [rec[idx], id])
        for f, c in list(null_count.items()):
            env.logger.warning(
                '{} missing values are identified for phenotype {}'.format(
                    c, f))
        if invalid_sample_names:
            env.logger.warning(
                'Samples {} in input file does not match any sample.'.format(
                    ', '.join(sorted(invalid_sample_names))))
        env.logger.info(
            '{} field ({} new, {} existing) phenotypes of {} samples are updated.'
            .format(
                count[1] + count[2], count[1], count[2],
                int(count[0] / (count[1] + count[2])) if
                (count[1] + count[2]) else 0))
        self.db.commit()

    def setPhenotype(self, field, expression, samples):
        '''Add a field using expression calculated from sample variant table'''
        IDs = self.proj.selectSampleByPhenotype(samples)
        if not IDs:
            raise ValueError(
                'No sample is selected using condition "{}"'.format(samples))
        #
        count = [0, 0, 0]
        cur = self.db.cursor()
        cur.execute('SELECT {} FROM sample;'.format(expression))
        fldType = None
        for rec in cur:
            if rec[0] is None:  # missing
                continue
            if fldType is None:
                fldType = type(rec[0])
            if type(rec[0]) != fldType:
                if type(rec[0]) is float and fldType is int:
                    fldType = float
                else:
                    raise ValueError(
                        'Inconsistent type returned from different samples')
        if expression != 'NULL' and fldType is None:
            raise ValueError('Cannot determine the type of the expression')
        # if adding a new field
        cur_fields = self.db.getHeaders('sample')[3:]
        if field.lower() not in [x.lower() for x in cur_fields]:
            if field.upper() in SQL_KEYWORDS:
                raise ValueError(
                    "Phenotype name '{}' is not allowed because it is a reserved word."
                    .format(field))
            env.logger.info('Adding phenotype {}'.format(field))
            self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(
                field, {
                    int: 'INT',
                    float: 'FLOAT',
                    str: 'VARCHAR(255)',
                    str: 'VARCHAR(255)',
                    None: 'FLOAT'
                }[fldType]))
            count[1] += 1  # new
        else:
            # FIXME: check the case for type mismatch
            count[2] += 1  # updated
        #
        cur = self.db.cursor()
        for ID in IDs:
            cur.execute(
                'UPDATE sample SET {0}={1} WHERE sample_id = {2}'.format(
                    field, None if expression == 'NULL' else expression,
                    self.db.PH), (ID,))
            count[0] += 1
        env.logger.info(
            '{} values of {} phenotypes ({} new, {} existing) of {} samples are updated.'
            .format(count[0], count[1] + count[2], count[1], count[2],
                    len(IDs)))
        self.db.commit()

    def fromSampleStat(self, stat, genotypes, samples):
        '''Add a field using expression calculated from sample variant table'''
        IDs = self.proj.selectSampleByPhenotype(samples)
        if not IDs:
            raise ValueError(
                'No sample is selected using condition "{}"'.format(samples))
        #
        # at least one, at most number of IDs
        nJobs = max(min(self.jobs, len(IDs)), 1)
        # start all workers
        idQueue = queue.Queue()
        idmpQueue = mpQueue()
        status = GenotypeStatStatus()
        mpStatus = Manager().dict()

        if self.proj.store == "sqlite":
            for j in range(nJobs):
                GenotypeStatCalculator('{}_genotype.DB'.format(self.proj.name),
                                       stat, idQueue, status,
                                       genotypes).start()
            for ID in IDs:
                idQueue.put(ID)
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
                    if field.upper() in SQL_KEYWORDS:
                        raise ValueError(
                            "Phenotype name '{}' is not allowed because it is a reserved word."
                            .format(field))
                    new_field[field] = True
                else:
                    new_field[field] = False
                    count[2] += 1  # updated
            cur = self.db.cursor()
            for ID in IDs:
                res = status.get(ID)
                for idx, (field, expr) in enumerate(stat):
                    if new_field[field]:
                        fldtype = typeOfValues(
                            [str(status.get(x)[idx]) for x in IDs])
                        # determine the type of value
                        self.db.execute(
                            'ALTER TABLE sample ADD {} {} NULL;'.format(
                                field, fldtype))
                        env.logger.debug(
                            'Adding phenotype {} of type {}'.format(
                                field, fldtype))
                        new_field[field] = False
                        count[1] += 1  # new
                    cur.execute(
                        'UPDATE sample SET {0}={1} WHERE sample_id = {1}'
                        .format(field, self.db.PH), [res[idx], ID])
                    count[0] += 1
        elif self.proj.store == "hdf5":
            for j in range(nJobs):
                GenotypeStatCalculator_HDF5(self.proj, stat, idmpQueue,
                                            mpStatus, genotypes).start()
            for ID in IDs:
                idmpQueue.put(ID)
            count = 0
            prog = ProgressBar('Calculating phenotype', len(IDs))
            while True:
                if len(mpStatus) > count:
                    count = len(mpStatus)
                    prog.update(count)
                # if everything is done
                if len(mpStatus) == len(IDs):
                    # stop all threads
                    for j in range(nJobs):
                        idmpQueue.put(None)
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
                    if field.upper() in SQL_KEYWORDS:
                        raise ValueError(
                            "Phenotype name '{}' is not allowed because it is a reserved word."
                            .format(field))
                    new_field[field] = True
                else:
                    new_field[field] = False
                    count[2] += 1  # updated
            cur = self.db.cursor()
            for ID in IDs:
                res = mpStatus[ID]
                for idx, (field, expr) in enumerate(stat):
                    if new_field[field]:
                        fldtype = typeOfValues(
                            [str(mpStatus[x][idx]) for x in IDs])
                        # determine the type of value
                        self.db.execute(
                            'ALTER TABLE sample ADD {} {} NULL;'.format(
                                field, fldtype))
                        env.logger.debug(
                            'Adding phenotype {} of type {}'.format(
                                field, fldtype))
                        new_field[field] = False
                        count[1] += 1  # new
                    cur.execute(
                        'UPDATE sample SET {0}={1} WHERE sample_id = {1}'
                        .format(field, self.db.PH), [res[idx], ID])
                    count[0] += 1
        # report result
        env.logger.info(
            '{} values of {} phenotypes ({} new, {} existing) of {} samples are updated.'
            .format(count[0], count[1] + count[2], count[1], count[2],
                    len(IDs)))
        self.db.commit()

    def output(self, fields, header, delimiter, na, limit, samples):
        # output
        limit_clause = '' if limit < 0 else ' LIMIT 0,{}'.format(limit)
        query = 'SELECT {} FROM sample LEFT JOIN filename ON sample.file_id = filename.file_id {} {}'.format(
            ','.join(fields), '' if not samples else 'WHERE ' + samples,
            limit_clause)
        env.logger.debug(query)
        cur = self.db.cursor()
        cur.execute(query)
        if header is not None:
            if len(header) == 0:
                print((delimiter.join(fields)))
            elif header == ['-']:
                env.logger.info('Reading header from standard input')
                print((sys.stdin.read().rstrip()))
            else:
                if len(header) != len(fields):
                    env.logger.warning(
                        'User-provided header ({}) does not match number of fields ({})'
                        .format(len(header), len(fields)))
                print((delimiter.join(header)))
        for rec in cur:
            print((delimiter.join([na if x is None else str(x) for x in rec])))


def phenotypeArguments(parser):
    '''Action that can be performed by this script'''
    parser.add_argument(
        '-f',
        '--from_file',
        '--from-file',
        metavar='INPUT_FILE',
        nargs='*',
        help='''Import phenotype from a tab or space delimited file, which can be
            standard input if a name - is specified. Samples are
            determined by sample names in the first column, or jointly by sample
            name and filename if there is another column with header 'filename'.
            Names of phenotype fields are determined by header of the input file,
            or by names provided from option --header. Non-alphanumeric
            characters in input filed names will be replaced by '_'. If
            multiple samples in a project share the same names, they will shared
            the imported phenotypes. Optionally, a list of phenotypes (columns
            of the file) can be specified after filename, in which case only the
            specified phenotypes will be imported. Parameter --samples could be
            used to limit the samples for which phenotypes are imported. Values
            that match value of parameter --na and cannot be converted to the
            probed type of phenotype (e.g. '' in a column of numbers) are recorded
            as missing values.'''),
    parser.add_argument(
        '--set',
        nargs='*',
        metavar='EXPRESSION',
        default=[],
        help='''Set a phenotype to a constant (e.g. --set aff=1), or an expression
            using other existing phenotypes (e.g. --set ratio_qt=high_qt/all_qt (the ratio
            of the number of high quality variants to the number of all variants, where
            high_qt and all_qt are obtained from sample statistics using parameter
            --from_stat). Parameter --samples could be used to limit the samples for
            which genotypes will be set.'''),
    parser.add_argument(
        '--from_stat',
        '--from-stat',
        nargs='*',
        metavar='EXPRESSION',
        default=[],
        help='''Set a phenotype to a summary statistics of a genotype field. For
            example, "num=count(*)" sets phenotype num to be the number of genotypes
            of a sample, "GD=avg(DP)" sets phenotype DP to be the average depth (if
            DP is one of the genotype fields) of the sample. Multiple fields (e.g.
            '--from-stat "num=count(*)" "GD=avg(DP)"') are also allowed. In addition to
            standard SQL aggregation functions, variant tools supports special functions
            #(GT), #(wtGT), #(mutGT), #(alt), #(hom), #(het) and #(other), which
            counts the number of genotypes (the same as count(*)), wildtype genotypes,
            mutant genotypes alternative alleles, homozygotes, heterozygotes, and
            genotypes with two different alternative alleles. Parameters --genotypes
            and --samples could be used to limit the genotypes to be considered and
            the samples for which genotypes will be set.'''),
    parser.add_argument(
        '--output',
        nargs='*',
        metavar='EXPRESSION',
        default=[],
        help='''A list of phenotype to be outputted. SQL-compatible expressions or
            functions such as "DP/DP_all" and "avg(DP)" are also allowed'''),
    parser.add_argument(
        '-j',
        '--jobs',
        metavar='N',
        default=4,
        type=int,
        help='''Allow at most N concurrent jobs to obtain sample statistics for
            parameter --from-stat.''')
    parser.add_argument(
        '-g',
        '--genotypes',
        nargs='*',
        metavar='COND',
        default=[],
        help='''Limit the operation to genotypes that match specified conditions.
            Use 'vtools show genotypes' to list usable fields for each sample.'''
    ),
    parser.add_argument(
        '-s',
        '--samples',
        nargs='*',
        metavar='COND',
        default=[],
        help='''Update phenotype for samples that match specified conditions.
            Use 'vtools show samples' to list usable fields in the sample table.'''
    )
    grp = parser.add_argument_group('Input/Output options')
    grp.add_argument(
        '--header',
        nargs='*',
        help='''A list of header names for input file if used with option
            --from-file. Otherwise a complete header or a list of names that
            will be joined by a delimiter (parameter --delimiter), for option
            --output. If a special name - is specified, the header will
            be read from the standard input, which is the preferred way to specify large
            multi-line headers (e.g. cat myheader | vtools export --header -). If this
            parameter is given without parameter, a default header will be derived from
            field names.'''),
    grp.add_argument(
        '-d',
        '--delimiter',
        default='\t',
        help='''Delimiter, default to tab, a popular alternative is ',' for csv output'''
    )
    grp.add_argument(
        '--na', default='NA', help='Input or output string for missing value..')
    grp.add_argument(
        '-l',
        '--limit',
        default=-1,
        type=int,
        help='''Number of record to display. Default to all record.''')


def phenotype(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            p = Sample(proj, args.jobs)
            if args.from_file:
                filename = args.from_file[0]
                fields = args.from_file[1:]
                p.load(filename, args.header, fields,
                       ' AND '.join(['({})'.format(x) for x in args.samples]),
                       args.na)
            if args.set:
                for item in args.set:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError(
                            'Invalid parameter {}, which should have format field=expr_of_phenotype: {}'
                            .format(item, e))
                    if field == "sample_name" or field == "filename":
                        raise ValueError(
                            'Cannot alter phenotype field: {}'.format(field))
                    p.setPhenotype(
                        field, expr,
                        ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.from_stat:
                stat = []
                for item in args.from_stat:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError(
                            'Invalid parameter {}, which should have format field=expr_of_field: {}'
                            .format(item, e))
                    stat.append((field, expr))
                p.fromSampleStat(
                    stat,
                    ' AND '.join(['({})'.format(x) for x in args.genotypes]),
                    ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.output:
                p.output(args.output, args.header, args.delimiter, args.na,
                         args.limit,
                         ' AND '.join(['({})'.format(x) for x in args.samples]))
            if not (args.from_file or args.set or args.from_stat or
                    args.output):
                raise ValueError(
                    'Please add "-h" after phenotype to get more help for this command'
                )
        proj.close()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
