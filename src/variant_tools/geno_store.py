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

import glob
import os
import queue
import re
import sys
import time
from multiprocessing import Lock, Manager, Process, Queue, Value

import numpy as np

from .accessor import Engine_Access, Engine_Storage
from .exporter_reader import MultiVariantReader
from .importer_allele_hdf5 import (HDF5GenotypeImportWorker, manageHDF5,
                                   updateSample)
from .text_reader import TextReader
from .utils import (DatabaseEngine, ProgressBar, ProgressFileObj, chunks,
                    decodeTableName, delayedAction, encodeTableName, env)


class Base_Store(object):

    def __init__(self, proj):
        self.proj = proj

    def importGenotypes(self, importer):
        raise RuntimeError('A derived object should be used')


#
#  A status object to control the import process. I cannot use a simple
#  queue solution because different types of importers are needed to handle
#  different inputs.
#
class ImportStatus:

    def __init__(self):
        self.manager = Manager()
        self.tasks = self.manager.dict()
        self.pending_dedup = self.manager.list()
        self.pending_copied = self.manager.dict()
        self.lock = Lock()
        self.total_sample_count = 0
        self.total_genotype_count = 0
        self.all_done = Value('L', 0)

    def add(self, item, num_lines):
        '''Add a job, each item has
        input_filename:     source filename
        tmp_file:           temporary genotype file
        start_sample, end_sample: index of samples
        sample_ids:         actual sample ids

        Status of the items can be

        0:  initial input, not imported
        1:  being imported
        2:  imported and being dedupped
        3:  dedupped
        4:  being copied
        5:  copied (done)
        '''
        if item in self.tasks:
            raise RuntimeError('Item already added to ImportStatus')
        self.lock.acquire()
        self.tasks[item] = 0
        self.total_sample_count += len(item[-1])
        self.total_genotype_count += num_lines * len(item[-1])
        self.lock.release()

    def itemToImport(self, filelist):
        '''Try to find some item that can be imported. filelist is a list of files
        that the importer can work on
        '''
        ret = None
        self.lock.acquire()
        for item in sorted(self.tasks.keys()):
            # only check for status 0 items
            if self.tasks[item] == 0 and item[0] in filelist:
                self.tasks[item] = 1
                ret = item
                break
        self.lock.release()
        return ret

    def pendingImport(self):
        # return waiting (queued, to be imported) and ongoing import
        return sum([x == 0 for x in list(self.tasks.values())
                   ]), sum([x == 1 for x in list(self.tasks.values())])

    def addDedupItem(self, geno_file, geno_status, IDs):
        self.lock.acquire()
        self.pending_dedup.append((geno_file, geno_status, IDs))
        self.lock.release()

    def itemToDedup(self):
        ret = None
        self.lock.acquire()
        if self.pending_dedup:
            ret = self.pending_dedup.pop(0)
        self.lock.release()
        return ret

    def addCopyingItem(self, geno_file, geno_status, ID, dup_rows):
        #
        self.lock.acquire()
        if (geno_file, geno_status) in self.pending_copied:
            self.pending_copied[(geno_file, geno_status)] = \
                self.pending_copied[(geno_file, geno_status)] + [(ID, dup_rows)]
        else:
            self.pending_copied[(geno_file, geno_status)] = [(ID, dup_rows)]
        self.lock.release()

    def itemToCopy(self):
        ret = None
        self.lock.acquire()
        if not self.pending_copied:
            ret = None
        else:
            # return any item
            key = list(self.pending_copied.keys())[0]
            ret = key[0], key[1], self.pending_copied.pop(key)
        self.lock.release()
        return ret

    def set(self, item, status):
        self.lock.acquire()
        self.tasks[item] = status
        self.lock.release()


#
#
# Write genotype to disk
#
#


class BaseGenotypeWriter:

    def __init__(self, geno_info, genotype_status, sample_ids):
        self.sample_ids = sample_ids
        self.geno_info = geno_info
        self.geno_status = genotype_status

    def _createNewSampleVariantTable(self, cur, table, genotype=True):
        '''Create a table ``genotype_??`` to store genotype data'''
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS {0} (
                variant_id INT NOT NULL
            '''.format(table) + (
            ', GT INT' +
            ''.join([', {} {}'.format(f.name, f.type)
                     for f in self.geno_info]) if genotype else '') + ');')


class MultiGenotypeWriter(BaseGenotypeWriter):
    '''This class write genotypes to a genotype database, which does
    not have to be the project genotype database.'''

    def __init__(self, geno_db, geno_info, genotype_status, sample_ids):
        '''geno_db:   genotype database
        geno_info:    genotype information fields
        sample_ids:   ID of samples that will be written to this database
        '''
        BaseGenotypeWriter.__init__(self, geno_info, genotype_status,
                                    sample_ids)
        # we save genotypes to many small files, with a minimal of 5 samples
        # and a maximum of 10
        nDBs = max(min(10, len(sample_ids) // 5), 1)
        #
        self.dispatcher = {}
        self.geno_db = []
        self.db = []
        self.cur = []
        for x in range(nDBs):
            self.geno_db.append(geno_db.replace('.DB', '_{}.DB'.format(x)))
            self.db.append(DatabaseEngine())
            self.db[-1].connect(self.geno_db[-1])
            self.cur.append(self.db[-1].cursor())
        #
        if self.geno_status == 2:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'\
                .format(','.join([self.db[0].PH] * (2 + len(geno_info))))
        else:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'.format(
                self.db[0].PH)
        #
        with delayedAction(
                env.logger.debug,
                'Creating {} genotype tables'.format(len(self.sample_ids))):
            for idx, sid in enumerate(self.sample_ids):
                self.dispatcher[sid] = idx % nDBs
                # create table
                self._createNewSampleVariantTable(self.cur[idx % nDBs],
                                                  'genotype_{0}'.format(sid),
                                                  self.geno_status == 2)
                self.db[idx % nDBs].commit()
        #
        # number of genotype batches that has been written
        self.count = 0
        # cache of genotypes. This class will accumulate 1000 genotypes before
        # it writes to the disk using 'executemany', which will be faster than
        # calling 1000 'execute'.
        self.cache = {}

    def write(self, id, rec):
        try:
            if len(self.cache[id]) == 1000:
                # this will fail if id does not exist, so we do not need
                # extra test if id is valid
                self.cur[self.dispatcher[id]].executemany(
                    self.query.format(id), self.cache[id])
                self.cache[id] = [rec]
                self.count += 1
            else:
                self.cache[id].append(rec)
        except KeyError:
            # if this is a new id
            self.cache[id] = [rec]

    def commit_remaining(self):
        for id, val in self.cache.items():
            if len(val) > 0:
                self.cur[self.dispatcher[id]].executemany(
                    self.query.format(id), val)
        for db in self.db:
            db.commit()
            db.close()

    def dedup(self, status=None):
        # this part is done by a dedicated processor
        for idx, db in enumerate(self.geno_db):
            status.addDedupItem(
                db, self.geno_status,
                [id for id in self.sample_ids if self.dispatcher[id] == idx])


class InPlaceGenotypeWriter(BaseGenotypeWriter):
    '''This class write genotypes to a genotype database, which does
    not have to be the project genotype database.'''

    def __init__(self, geno_db, geno_info, genotype_status, sample_ids):
        '''geno_db:   genotype database
        geno_info:    genotype information fields
        sample_ids:   ID of samples that will be written to this database
        '''
        #
        BaseGenotypeWriter.__init__(self, geno_info, genotype_status,
                                    sample_ids)
        #
        self.geno_db = geno_db
        self.db = DatabaseEngine()
        self.db.connect(self.geno_db)
        if self.geno_status == 2:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'\
                .format(','.join([self.db.PH] * (2 + len(geno_info))))
        else:
            self.query = 'INSERT INTO genotype_{{}} VALUES ({0});'.format(
                self.db.PH)
        self.cur = self.db.cursor()
        with delayedAction(
                env.logger.info,
                'Creating {} genotype tables'.format(len(self.sample_ids))):
            for idx, sid in enumerate(self.sample_ids):
                # create table
                self._createNewSampleVariantTable(self.cur,
                                                  'genotype_{0}'.format(sid),
                                                  self.geno_status == 2)
            self.db.commit()
        #
        # number of genotype batches that has been written
        self.count = 0
        # cache of genotypes. This class will accumulate 1000 genotypes before
        # it writes to the disk using 'executemany', which will be faster than
        # calling 1000 'execute'.
        self.cache = {}

    def write(self, id, rec):
        try:
            if len(self.cache[id]) == 1000:
                # this will fail if id does not exist, so we do not need
                # extra test if id is valid
                self.cur.executemany(self.query.format(id), self.cache[id])
                self.cache[id] = [rec]
                self.count += 1
            else:
                self.cache[id].append(rec)
        except KeyError:
            # if this is a new id
            self.cache[id] = [rec]
        if self.count % 1000 == 0:
            self.db.commit()

    def commit_remaining(self):
        for id, val in self.cache.items():
            if len(val) > 0:
                self.cur.executemany(self.query.format(id), val)
        self.db.commit()
        self.db.close()

    def dedup(self, status=None):
        # checking if there are duplicated variant_ids in genotype tables
        # we do not do it during data insertion because (potentially) many tables
        # are handled simultenously, and keeping track of ids in each sample can
        # take a lot of ram.
        #
        db = DatabaseEngine()
        db.connect(self.geno_db)
        cur = db.cursor()
        for id in self.sample_ids:
            cur.execute('SELECT COUNT(variant_id), COUNT(DISTINCT variant_id) '
                        'FROM genotype_{}'.format(id))
            nRec, nVar = cur.fetchone()
            if nRec != nVar:
                cur.execute(
                    'DELETE from genotype_{0} WHERE rowid NOT IN '
                    '(SELECT MAX(rowid) FROM genotype_{0} GROUP BY variant_id)'
                    .format(id))
                if cur.rowcount != nRec - nVar:
                    raise SystemError(
                        'Failed to identify duplicated variants from '
                        'genotype table genotype_{}'.format(id))
                # cannot get variant id easily
                env.logger.debug(
                    '{} duplicated records have been removed from sample {}'
                    .format(cur.rowcount, id))
        db.commit()
        db.close()


def GenotypeWriter(geno_db, geno_info, genotype_status, sample_ids):
    if '/' in geno_db:
        return MultiGenotypeWriter(geno_db, geno_info, genotype_status,
                                   sample_ids)
    else:
        return InPlaceGenotypeWriter(geno_db, geno_info, genotype_status,
                                     sample_ids)


class GenotypeCopier(Process):

    def __init__(self, main_genotype_file, genotype_info, copied_samples,
                 status):
        '''copied_samples is a shared variable that should be increased with
        each sample copy.
        '''
        Process.__init__(self)
        self.daemon = True
        self.main_genotype_file = main_genotype_file
        self.genotype_info = genotype_info
        self.copied_samples = copied_samples
        self.status = status

    def run(self):
        db = None
        while True:
            try:
                item = self.status.itemToCopy()
            except Exception as e:
                env.logger.error('GenotypeCopier failed: {}'.format(e))
                sys.exit(1)
            if item is None:
                if self.status.all_done.value == 1:
                    if db is not None:
                        db.close()
                    env.logger.debug('Genotype of {} samples are copied'.format(
                        self.copied_samples.value))
                    break
                time.sleep(1)
                continue
            # only connect to the database engine when the main process is done
            if db is None:
                db = DatabaseEngine()
                db.connect(self.main_genotype_file)
            #
            cur = db.cursor()
            genotype_file, genotype_status, ID_and_DUPs = item
            #
            db.attach(genotype_file, '__from')
            # start copying genotype
            # copy genotype table
            start_copy_time = time.time()
            cur = db.cursor()
            for ID, rowids in ID_and_DUPs:
                if ID is None:
                    continue
                query = 'CREATE TABLE IF NOT EXISTS genotype_{0} (variant_id INT NOT NULL'.format(
                    ID)
                if genotype_status == 2:
                    query += ', GT INT' + ''.join([
                        ', {} {}'.format(f.name, f.type)
                        for f in self.genotype_info
                    ])
                query += ');'
                cur.execute(query)
                if rowids:
                    query = ('SELECT variant_id FROM  __from.genotype_{0} '
                             'WHERE rowid IN ({1});').format(
                                 ID, ','.join([str(x) for x in rowids]))
                    cur.execute(query)
                    var_ids = [x[0] for x in cur.fetchall()]
                    env.logger.debug(
                        'Removing {} records for variants {} from sample {}'
                        .format(
                            len(rowids), ', '.join([str(x) for x in var_ids]),
                            ID))
                    query = (
                        'INSERT INTO genotype_{0} SELECT * FROM __from.genotype_{0} '
                        'WHERE rowid NOT IN ({1});').format(
                            ID, ','.join([str(x) for x in rowids]))
                else:
                    query = 'INSERT INTO genotype_{0} SELECT * FROM __from.genotype_{0};'.format(
                        ID)
                cur.execute(query)
                db.commit()
                # update progress
                self.copied_samples.value += 1
            db.detach('__from')
            end_copy_time = time.time()
            env.logger.debug(
                'Copying {} samples from {} took {:.1f} seconds.'.format(
                    len([x for x in ID_and_DUPs if x[0] is not None]),
                    os.path.basename(genotype_file),
                    end_copy_time - start_copy_time))
            # if no IDs, all samples have been copied.
            if None in ID_and_DUPs[-1]:
                try:
                    # remove the temporary file
                    os.remove(self.genotype_file)
                except:
                    pass


class DedupWorker(Process):

    def __init__(self, status):
        Process.__init__(self)
        self.daemon = True
        self.status = status

    def run(self):
        while True:
            item = self.status.itemToDedup()
            if item is None:
                if self.status.all_done.value == 1:
                    break
                time.sleep(1)
                continue
            genotype_file, genotype_status, IDs = item
            db = DatabaseEngine()
            db.connect(genotype_file, readonly=True)
            cur = db.cursor()
            #
            for id in IDs:
                cur.execute(
                    'SELECT COUNT(variant_id), COUNT(DISTINCT variant_id) '
                    'FROM genotype_{}'.format(id))
                nRec, nVar = cur.fetchone()
                if nRec != nVar:
                    cur.execute(
                        'SELECT rowid from genotype_{0} WHERE rowid NOT IN '
                        '(SELECT MAX(rowid) FROM genotype_{0} GROUP BY variant_id)'
                        .format(id))
                    deleted_rows = [x[0] for x in cur.fetchall()]
                    if len(deleted_rows) != nRec - nVar:
                        raise SystemError(
                            'Failed to identify duplicated variants from '
                            'genotype table genotype_{}'.format(id))
                else:
                    deleted_rows = []
                #
                self.status.addCopyingItem(genotype_file, genotype_status, id,
                                           deleted_rows)
            self.status.addCopyingItem(genotype_file, genotype_status, None,
                                       None)
            db.close()


#
#   A working process to import genotype from a file, or part of a file
#   and write to a temporary genotype database.
#
#
class GenotypeImportWorker(Process):
    '''This class starts a process, import genotype to a temporary genotype database.'''

    def __init__(self, variantIndex, filelist, processor, encoding, header,
                 genotype_field, genotype_info, ranges, geno_count, proc_index,
                 status):
        '''
        variantIndex: a dictionary that returns ID for each variant.
        filelist: files from which variantIndex is created. If the passed filename
            is not in this list, this worker will suicide so that it can be replaced
            by a worker with more variants.
        encoding, genotypefield, genotype_info, ranges:  parameters to import data
        geno_count:  a shared variable to report number of genotype imported
        status:      an ImportStatus object to monitor the progress
        '''
        Process.__init__(self, name='GenotypeImporter')
        self.daemon = True
        self.variantIndex = variantIndex
        self.filelist = filelist
        self.encoding = encoding
        self.header = header
        self.processor = processor
        #
        self.genotype_field = genotype_field
        self.genotype_info = genotype_info
        self.ranges = ranges
        #
        self.geno_count = geno_count
        self.proc_index = proc_index
        #
        self.status = status

    def run(self):
        env.logger.debug(
            'Importer {} starts with variants from {} files'.format(
                self.proc_index, len(self.filelist)))
        while True:
            item = self.status.itemToImport(self.filelist)
            if item is None:
                # wait a second, make sure there is no job
                time.sleep(1)
                item = self.status.itemToImport(self.filelist)
                if item is None:
                    env.logger.debug('Importer {} exits normally'.format(
                        self.proc_index))
                    break
            # get parameters
            self.input_filename, self.genotype_file, self.genotype_status, start_sample, end_sample, self.sample_ids = item
            self.processor.reset(
                import_var_info=False,
                import_sample_range=[0, 0]
                if self.genotype_status == 1 else [start_sample, end_sample])
            self.count = [0, 0]
            # use the last value as start
            self.start_count = self.geno_count.value
            start_import_time = time.time()
            self._importData()
            # set the status to be imported (2) (and being dedupped)
            self.status.set(item, 2)
            env.logger.debug(
                'Importer {} starts deduplicating {} samples after importing genotypes in {:.1f} seconds'
                .format(self.proc_index, len(self.sample_ids),
                        time.time() - start_import_time))
            self._dedupData()
            self.status.set(item, 3)
            end_import_time = time.time()
            todo, going = self.status.pendingImport()
            env.logger.debug(
                'Importing {} samples ({} - {}) to {} took importer {} {:.1f} seconds, {} onging, {} to go.'
                .format(
                    len(self.sample_ids), min(self.sample_ids),
                    max(self.sample_ids), os.path.basename(self.genotype_file),
                    self.proc_index, end_import_time - start_import_time, going,
                    todo))

    def _importData(self):
        env.logger.debug(
            'Importer {} starts importing genotypes for {} samples ({} - {})'
            .format(self.proc_index, len(self.sample_ids), min(self.sample_ids),
                    max(self.sample_ids)))
        reader = TextReader(
            self.processor,
            self.input_filename,
            None,
            False,
            0,
            self.encoding,
            self.header,
            quiet=True)
        self.writer = GenotypeWriter(self.genotype_file, self.genotype_info,
                                     self.genotype_status, self.sample_ids)
        fld_cols = None
        last_count = 0
        for self.count[0], bins, rec in reader.records():
            try:
                variant_id = self.variantIndex[tuple(
                    (rec[0], rec[2], rec[3]))][rec[1]][0]
            except KeyError:
                env.logger.debug('Variant {} {} {} {} not found'.format(
                    rec[0], rec[1], rec[2], rec[3]))
                continue
            # if there is genotype
            if self.genotype_status == 2:
                if fld_cols is None:
                    col_rngs = [
                        reader.columnRange[x]
                        for x in range(self.ranges[2], self.ranges[4])
                    ]
                    fld_cols = []
                    for idx in range(len(self.sample_ids)):
                        fld_cols.append([
                            sc + (0 if sc + 1 == ec else idx)
                            for sc, ec in col_rngs
                        ])
                    if col_rngs[0][1] - col_rngs[0][0] != len(self.sample_ids):
                        env.logger.error(
                            'Number of genotypes ({}) does not match number of samples ({})'
                            .format(col_rngs[0][1] - col_rngs[0][0],
                                    len(self.sample_ids)))
                for idx, id in enumerate(self.sample_ids):
                    try:
                        # variant info is not read, ranges[1] should be used because ranges[2] is the index after
                        # variant index
                        if rec[self.ranges[1] + idx] is not None:
                            self.count[1] += 1
                            self.writer.write(id, [variant_id] +
                                              [rec[c] for c in fld_cols[idx]])
                    except IndexError:
                        env.logger.warning(
                            'Incorrect number of genotype fields: {} fields found, {} expected for record {}'
                            .format(len(rec), fld_cols[-1][-1] + 1, rec))
            else:
                # should have only one sample
                for id in self.sample_ids:
                    self.writer.write(id, [variant_id])
            if self.count[0] - last_count > 100:
                self.geno_count.value = self.start_count + self.count[0] * len(
                    self.sample_ids)
                last_count = self.count[0]
        self.writer.commit_remaining()

    def _dedupData(self):
        self.writer.dedup(self.status)


class Sqlite_Store(Base_Store):

    def __init__(self, proj):
        super(Sqlite_Store, self).__init__(proj)
        self.db = DatabaseEngine()
        self.db.connect('{}_genotype.DB'.format(self.proj.name))
        self.cur = self.db.cursor()

    def getGenotypeFileSize(self):
        if os.path.isfile('{}_genotype.DB'.format(self.proj.name)):
            filesize = os.path.getsize('{}_genotype.DB'.format(self.proj.name))
        return filesize

    def load_Genotype_From_SQLite(self, all_files, proj, importer=None):
        if os.path.isfile('{}_genotype.DB'.format(self.proj.name)):
            os.remove('{}_genotype.DB'.format(self.proj.name))
        if 'snapshot_genotype.DB' in all_files:
            os.rename(
                os.path.join(env.cache_dir, 'snapshot_genotype.DB'),
                '{}_genotype.DB'.format(self.proj.name))
            all_files.remove('snapshot_genotype.DB')
        elif '{}_genotype.DB'.format(self.proj.name) in all_files:
            # an old version of snapshot saves $name.proj
            os.rename(
                os.path.join(env.cache_dir,
                             '{}_genotype.DB'.format(self.proj.name)),
                '{}_genotype.DB'.format(self.proj.name))
            all_files.remove('{}_genotype.DB'.format(self.proj.name))
        else:
            # this is ok because a project might not have any genotype
            pass
        return all_files

    def addGenotypeToTar(self, tarfile, prog):
        if os.path.isfile('{}_genotype.DB'.format(self.proj.name)):
            tarinfo = tarfile.gettarinfo(
                '{}_genotype.DB'.format(self.proj.name),
                arcname='snapshot_genotype.DB')
            tarfile.addfile(
                tarinfo,
                ProgressFileObj(prog, '{}_genotype.DB'.format(self.proj.name),
                                'rb'))
        return tarfile

    def num_variants(self, sample_id):
        self.cur.execute('SELECT count(*) FROM genotype_{};'.format(sample_id))
        return self.cur.fetchone()[0]

    def geno_fields(self, sample_id):
        # sampleGenotypeHeader = [x.lower() for x in self.db.getHeaders('genotype_{}'.format(sample_id))]
        sampleGenotypeHeader = [
            x.lower()
            for x in self.db.getHeaders('genotype_{}'.format(sample_id))
        ]
        return sampleGenotypeHeader[
            1:]  # the first field is variant id, second is GT

    def geno_fields_nolower(self, sample_id):
        # sampleGenotypeHeader = [x.lower() for x in self.db.getHeaders('genotype_{}'.format(sample_id))]
        sampleGenotypeHeader = [
            x for x in self.db.getHeaders('genotype_{}'.format(sample_id))
        ]
        return sampleGenotypeHeader[
            1:]  # the first field is variant id, second is GT

    def get_typeOfColumn(self, sample_id, field):
        return self.db.typeOfColumn('genotype_{}'.format(sample_id), field)

    def remove_genofields(self, IDs, items):
        cnt = 0
        for table in ['genotype_{}'.format(id) for id in IDs]:
            header = [x.lower() for x in self.db.getHeaders(table)]
            items = [
                x for x in items
                if x.lower() in header and x.lower not in ['variant_id', 'gt']
            ]
            if items:
                cnt += 1
                env.logger.info(
                    'Removing fields {} from genotype table {}'.format(
                        ', '.join(items),
                        table.split('_')[-1]))
                self.db.removeFields(table, items)
        if cnt:
            env.logger.info(
                'Genotype info fields {} are removed from {} samples.'.format(
                    ', '.join(items), cnt))
        else:
            env.logger.warning(
                'Genotype info fields {} not found in any of the samples.'
                .format(', '.join(items)))

    def remove_sample(self, IDs):
        for sample_id in IDs:
            self.db.removeTable('genotype_{}'.format(sample_id))
            self.db.commit()

    def remove_variants(self, variantIDs, table):
        # get sample_ids
        self.db.attach(self.proj.name + ".proj", self.proj.name)
        self.cur.execute('SELECT sample_id, sample_name FROM {}.sample;'.format(
            self.proj.name))
        samples = self.cur.fetchall()

        for ID, name in samples:
            if not self.db.hasIndex('genotype_{0}_index'.format(ID)):
                self.cur.execute(
                    'CREATE INDEX genotype_{0}_index ON genotype_{0} (variant_id ASC)'
                    .format(ID))
            self.cur.execute('DELETE FROM genotype_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                .format(ID, self.proj.name+"."+table))

            env.logger.info('{} genotypes are removed from sample {}'.format(
                self.cur.rowcount, name))
        # remove the table itself
        env.logger.info('Removing table {} itself'.format(
            decodeTableName(table)))
        self.db.removeTable(self.proj.name + "." + table)
        self.db.commit()

    def remove_genotype(self, cond):
        # get sample_ids
        self.db.attach(self.proj.name + ".proj", self.proj.name)
        self.cur.execute('SELECT sample_id, sample_name FROM {}.sample;'.format(
            self.proj.name))
        samples = self.cur.fetchall()
        env.logger.info(
            'Removing genotypes from {} samples using criteria "{}"'.format(
                len(samples), cond))
        for ID, name in samples:
            try:
                # self.cur.execute('DELETE FROM {}_genotype.genotype_{} WHERE {};'\
                #     .format(self.proj.name, ID, cond))
                self.cur.execute('DELETE FROM genotype_{} WHERE {};'\
                    .format(ID, cond))
            except Exception as e:
                env.logger.warning(
                    'Failed to remove genotypes from sample {}: {}'.format(
                        name, e))
                continue
            env.logger.info('{} genotypes are removed from sample {}'.format(
                self.cur.rowcount, name))
        self.db.commit()

    def get_noWT_variants(self, IDs, proj, where_clause, args):
        hasGT = {
            id: 'GT' in [
                x[0]
                for x in proj.db.fieldsOfTable('{}_genotype.genotype_{}'.format(
                    proj.name, id))
            ] for id in IDs
        }
        noWT_clause = {
            id: 'WHERE {0}_genotype.genotype_{1}.GT != 0'.format(proj.name, id)
            if hasGT[id] else '' for id in IDs
        }
        #
        if len(IDs) == 0:
            env.logger.warning('No sample is selected by condition: {}'.format(
                ' AND '.join(['({})'.format(x) for x in args.samples])))
            # nothing will be selected
            where_clause += ' AND 0'
        #
        # This special case does not hold because sometimes variants are imported without sample information.
        #
        #elif len(IDs) == proj.db.numOfRows('sample'):
        #    env.logger.info('All {} samples are selected by condition: {}'.format(len(IDs), ' AND '.join(args.samples)))
        #    # we do not have to add anything to where_clause
        elif len(IDs) < 50:
            # we allow 14 tables in other 'union' or from condition...
            env.logger.info('{} samples are selected by condition: {}'.format(
                len(IDs), ' AND '.join(args.samples)))
            where_clause += ' AND ({}.variant_id IN ({}))'.format(
                encodeTableName(args.from_table), '\nUNION '.join([
                    'SELECT variant_id FROM {}_genotype.genotype_{} {}'.format(
                        proj.name, id, noWT_clause[id]) for id in IDs
                ]))
        else:
            # we have to create a temporary table and select variants sample by sample
            # this could be done in parallel if there are a large number of samples, but that needs a lot more
            # code, and perhaps RAM
            env.logger.info('{} samples are selected by condition: {}'.format(
                len(IDs), ' AND '.join(args.samples)))
            cur = proj.db.cursor()
            BLOCK_SIZE = 64
            NUM_BLOCKS = len(IDs) // BLOCK_SIZE + 1
            myIDs = list(IDs)
            myIDs.sort()
            merged_table = '__variants_from_samples'
            query = 'CREATE TEMPORARY TABLE {} (variant_id INT);'.format(
                merged_table)
            cur.execute(query)
            prog = ProgressBar('Collecting sample variants',
                               len(IDs)) if NUM_BLOCKS > 1 else None
            count = 0
            for i in range(NUM_BLOCKS):
                # step 1: create a table that holds all
                block_IDs = myIDs[(i * BLOCK_SIZE):((i + 1) * BLOCK_SIZE)]
                if len(block_IDs) == 0:
                    continue
                query = 'INSERT INTO {} {};'.format(
                    merged_table, '\nUNION '.join([
                        'SELECT variant_id FROM {}_genotype.genotype_{} {}'
                        .format(proj.name, id, noWT_clause[id])
                        for id in block_IDs
                    ]))
                cur.execute(query)
                count += len(block_IDs)
                if prog:
                    prog.update(count)
            if prog:
                prog.done()
            where_clause += ' AND ({}.variant_id IN (SELECT variant_id FROM {}))'.format(
                encodeTableName(args.from_table), merged_table)
        return where_clause

    def get_genoType_genoInfo(self, sampleDict, genotypes, variant_table,
                              genotypeFields, validGenotypeIndices,
                              validGenotypeFields, operations, fieldCalcs, prog,
                              prog_step):
        MEAN = 0
        SUM = 1
        MIN = 2
        MAX = 3
        variants = dict()
        id_idx = 0
        for id, sampleInfo in sampleDict.items():
            record_male_gt = sampleInfo[0]
            fieldSelect = sampleInfo[1]
            if not fieldSelect or all([x == 'NULL' for x in fieldSelect]):
                continue

            where_cond = []
            if genotypes is not None and len(genotypes) != 0:
                where_cond.extend(genotypes)
            if variant_table != 'variant':
                where_cond.append(
                    'variant_id in (SELECT variant_id FROM {})'.format(
                        self.proj.name + "." + variant_table))
            whereClause = 'WHERE ' + ' AND '.join(
                ['({})'.format(x) for x in where_cond]) if where_cond else ''

            query = 'SELECT variant_id {} FROM {}_genotype.genotype_{} {};'.format(
                ' '.join([',' + x for x in fieldSelect]), self.proj.name, id,
                whereClause)
            self.db.attach(self.proj.name + "_genotype.DB",
                           self.proj.name + "_genotype")
            self.db.attach(self.proj.name + ".proj", self.proj.name)
            env.logger.trace(query)
            self.cur.execute(query)
            result = self.cur.fetchall()
            self.db.commit()
            self.db.detach(self.proj.name + "_genotype")
            self.db.detach(self.proj.name)
            id_idx += 1

            for rec in result:
                if rec[0] not in variants:
                    # the last item is for number of genotype for male individual
                    variants[rec[0]] = [0, 0, 0, 0, 0]
                    variants[rec[0]].extend(list(fieldCalcs))
                # total valid GT
                if rec[1] is not None:
                    variants[rec[0]][3] += 1
                # if tracking genotype of males (for maf), and the sex is male
                if record_male_gt:
                    variants[rec[0]][4] += 1
                # type heterozygote
                if rec[1] == 1:
                    variants[rec[0]][0] += 1
                # type homozygote
                elif rec[1] == 2:
                    # if rec[0]==7:
                    #     print(id)
                    variants[rec[0]][1] += 1
                # type double heterozygote with two different alternative alleles
                elif rec[1] == -1:
                    variants[rec[0]][2] += 1
                elif rec[1] not in [0, None]:
                    env.logger.warning('Invalid genotype type {}'.format(
                        rec[1]))
                #
                # this collects genotype_field information
                if len(validGenotypeFields) > 0:
                    for index in validGenotypeIndices:
                        queryIndex = index + 2  # to move beyond the variant_id and GT fields in the select statement
                        recIndex = index + 5  # first 5 attributes of variants are het, hom, double_het, GT, GT in males
                        # ignore missing (NULL) values, and empty string that, if so inserted, could be returned
                        # by sqlite even when the field type is INT.
                        if rec[queryIndex] in [None, '', '.']:
                            continue
                        operation = operations[index]
                        #field = genotypeFields[index]
                        if operation == MEAN:
                            if variants[rec[0]][recIndex] is None:
                                # we need to track the number of valid records
                                # variants[rec[0]][recIndex] = [rec[queryIndex], 1]
                                variants[rec[0]][recIndex] = rec[queryIndex]
                            else:
                                # variants[rec[0]][recIndex][0] += rec[queryIndex]
                                # variants[rec[0]][recIndex][1] += 1
                                variants[rec[0]][recIndex] += rec[queryIndex]
                        elif operation == SUM:
                            if variants[rec[0]][recIndex] is None:
                                variants[rec[0]][recIndex] = rec[queryIndex]
                            else:
                                variants[rec[0]][recIndex] += rec[queryIndex]
                        elif operation == MIN:
                            if variants[rec[0]][recIndex] is None or rec[
                                    queryIndex] < variants[rec[0]][recIndex]:
                                variants[rec[0]][recIndex] = rec[queryIndex]
                        elif operation == MAX:
                            if variants[rec[0]][recIndex] is None or rec[
                                    queryIndex] > variants[rec[0]][recIndex]:
                                variants[rec[0]][recIndex] = rec[queryIndex]
            if id_idx % prog_step == 0:
                prog.update(id_idx + 1)
        prog.done()
        return variants

    def importGenotypes(self, importer):
        '''import files in parallel, by importing variants and genotypes separately, and in their own processes.
        More specifically, suppose that there are three files

        file1: variant1, sample1.1, sample1.2, sample1.3
        file2: variant2, sample2.1, sample2.2
        file3: variant3, sample3.1, sample3.2

        where variant1, 2, and 3 are three potentially overlapping sets of variants.
        sample1, sample2, sample3 are three groups of samples that are divided by the number of samples
        and number of processes (importer.jobs) (say, 3000 samples divided into three groups of 1000 samples).

        Then, there are

        A: 2 TextReader processes <--> main process read variant1, 2, 3 --> master variant table
        B: importer.jobs GenotypeImportWorker reads sample 1.1, 1.2, 1.3, 2.1 etc ... --> temporary genotype tables
           except for the first process, which writes to the main genotype table directly.
        C: 1 GenotypeCopier -> copy temporary genotype tables to the master genotype table

        The processes are organized so that
        1. sample A.x if read after variant A is read
        2. genotypes are copied after temporary tables are completed

        Because of the processes work together, although there are multiple progress bars, they
        might start from the middle because some work might have already been done in the previous step.
        '''
        importers = [None] * importer.jobs
        # number of genotypes each process have imported
        genotype_import_count = [Value('L', 0) for x in range(importer.jobs)]
        sample_copy_count = Value('L', 0)
        # import queue that accept jobs sample 1.1, 1.2, etc
        status = ImportStatus()
        # start copier
        copier = GenotypeCopier('{}_genotype.DB'.format(self.proj.name),
                                importer.genotype_info, sample_copy_count,
                                status)
        copier.start()
        #
        dedupier = []
        for i in range(max(2, min(importer.jobs, 4))):
            d = DedupWorker(status)
            d.start()
            dedupier.append(d)
        #
        # The logic of importer is complex here. Because an importer needs to know variantIndex to
        # write genotype tables, and because importers starts after each file is read, an importer
        # created earlier (eg. from file 1) cannot be used to import genotype for file 2. This is
        # why we
        #
        # 1. create importer only after a file is processed.
        # 2. If all workers are busy, no more new worker will be created.
        # 3. When an old importer finds that it cannot process new genotype file, it will commit succide.
        # 4. The master process will create new importers with updated variantIndex when
        #    there are empty slots.
        #
        # We do not pass variantIndex to importers because variantIndex is huge and it is slow
        # to pass it with other parameters.
        #
        # process each file
        for count, input_filename in enumerate(importer.files):
            env.logger.info('{} variants from {} ({}/{})'.format(
                'Importing', input_filename, count + 1, len(importer.files)))
            importer.importVariant(input_filename)
            env.logger.info('{:,} new variants {}{}{} from {:,} lines are imported.'\
                .format(importer.count[2], "(" if importer.count[2] else '',
                    ', '.join(['{:,} {}'.format(x, y) for x, y in \
                        zip(importer.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'unsupported']) if x > 0]),
                        ")" if importer.count[2] else '', importer.count[0]))
            # genotypes?
            if importer.genotype_field:
                importer.prober.reset()
            # if there are samples?
            sample_ids, genotype_status, names = importer.getSampleIDs(
                input_filename)
            #
            # we should have file line count from importVariant
            num_of_lines = importer.count[0]
            #
            for i in range(len(importer.count)):
                importer.total_count[i] += importer.count[i]
                importer.count[i] = 0
            #
            if len(sample_ids) == 0:
                continue
            #
            # determine workload:
            # from our benchmark, if there are a large number of jobs and if
            # we split jobs evenly, the last trunk will
            # take about double time to finish because the extra time to reach the end
            # if the lines are long. Therefore, we are using an algorithm that the last
            # piece will handle 2/3 of the samples of the first one.
            #
            # n -- len(sample_ids)
            # m -- number of processes
            # k -- workload of the first one
            # k - k/3(m-1) -- workload of the second one
            # ...
            # 2k/3 -- workload of the last one
            #
            # k = 6n/(5m)
            # d=k/3(m-1)
            #
            # k-xd (x=0, ..., m-1)
            #
            # each process handle at least 10 samples
            if len(sample_ids) > 2000:
                workload = [
                    max(
                        10,
                        int((1.2 * len(sample_ids) / importer.jobs) *
                            (1 - x / (3. * (importer.jobs - 1)))))
                    for x in range(importer.jobs)
                ]
            else:
                workload = [
                    max(10, int(float(len(sample_ids)) / importer.jobs))
                ] * importer.jobs
            # if there are missing ones, spread it across workers ...
            # less than 0 is possible because of the at least 10 policy
            unallocated = max(0, len(sample_ids) - sum(workload))
            for i in range(unallocated):
                workload[i % importer.jobs] += 1
            #
            env.logger.debug('Workload of processes: {}'.format(workload))

            start_sample = 0
            for job in range(importer.jobs):
                if workload[job] == 0:
                    continue
                end_sample = min(start_sample + workload[job], len(sample_ids))
                if end_sample <= start_sample:
                    continue
                # tell the processor do not import variant info, import part of the sample
                tmp_file = os.path.join(
                    env.temp_dir, 'tmp_{}_{}_genotype.DB'.format(count, job))
                if os.path.isfile(tmp_file):
                    os.remove(tmp_file)
                if os.path.isfile(tmp_file):
                    raise RuntimeError(
                        'Failed to remove existing temporary '
                        'database {}. Remove clean your cache directory'.format(
                            tmp_file))
                # send a import job to the importer workers
                status.add(
                    (input_filename, tmp_file, genotype_status, start_sample,
                     end_sample, tuple(sample_ids[start_sample:end_sample])),
                    num_of_lines)
                # start an importer if needed
                for i in range(importer.jobs):
                    if importers[i] is None or not importers[i].is_alive():
                        importers[i] = GenotypeImportWorker(
                            importer.variantIndex, importer.files[:count + 1],
                            importer.processor, importer.encoding,
                            importer.header, importer.genotype_field,
                            importer.genotype_info, importer.ranges,
                            genotype_import_count[i], i, status)
                        importers[i].start()
                        break
                start_sample = end_sample
        #
        # monitor the import of genotypes
        prog = ProgressBar(
            'Importing genotypes',
            status.total_genotype_count,
            initCount=sum([x.value for x in genotype_import_count]))
        while True:
            # each process update their passed shared value
            # the master process add them and get the total number of genotypes imported
            prog.update(sum([x.value for x in genotype_import_count]))
            #
            queued, importing = status.pendingImport()
            if queued + importing == 0:
                # the importer will kill themselves after there is no pending job
                prog.done()
                break
            # create importers if any of the importer is gone. This might be a waste of resource
            # but an importer that is pending does not cost much
            if queued > 0:
                new_count = 0
                for i in range(importer.jobs):
                    if importers[i] is None or not importers[i].is_alive():
                        worker = GenotypeImportWorker(
                            importer.variantIndex, importer.files,
                            importer.processor, importer.encoding,
                            importer.header, importer.genotype_field,
                            importer.genotype_info, importer.ranges,
                            genotype_import_count[i], i, status)
                        importers[i] = worker
                        worker.start()
                        new_count += 1
                    if new_count >= queued:
                        break
            time.sleep(2)
        # monitor the dedup of genotypes
        prog = ProgressBar(
            'Copying samples',
            status.total_sample_count,
            initCount=sample_copy_count.value)
        while True:
            prog.update(sample_copy_count.value)
            if sample_copy_count.value == status.total_sample_count:
                prog.done()
                status.all_done.value = 1
                copier.join()
                [d.join() for d in dedupier]
                break
            time.sleep(1)
        # final status line
        if len(importer.files) > 1:
            env.logger.info('{:,} new variants ({}) from {:,} lines ({:,} samples) are imported.'\
                .format(importer.total_count[2],
                    ', '.join(['{:,} {}'.format(x, y) for x, y in \
                        zip(importer.total_count[3:8], ['SNVs', 'insertions',
                        'deletions', 'complex variants', 'unsupported']) if x > 0]),
                    importer.total_count[0], status.total_sample_count))

    def get_Genotype(self, cur, table, proj, sample_ID):
        g = set()
        NULL_to_0 = env.treat_missing_as_wildtype
        cur.execute(
            'SELECT {0}.variant_id, geno.GT FROM {0} LEFT OUTER JOIN '
            '{1}_genotype.genotype_{2} geno ON {0}.variant_id = geno.variant_id'
            .format(encodeTableName(table), proj.name, sample_ID))
        for id, GT in cur:
            if GT is None:
                if NULL_to_0:
                    g.add((id, 0))
            else:
                g.add((id, GT))
        return g


class HDF5_Store(Base_Store):

    def __init__(self, proj):
        self.proj = proj
        self.cache = {}
        super(HDF5_Store, self).__init__(proj)

    def getGenotypeFileSize(self):
        hdf5files = glob.glob("tmp*h5")
        filesizes = 0
        if len(hdf5files) > 0:
            for hdf5file in hdf5files:
                filesizes += os.path.getsize(hdf5file)
        return filesizes

    def load_Genotype_From_SQLite(self, all_files, proj, importer=None):
        hdf5files = glob.glob("tmp*h5")
        cur = self.proj.db.cursor()
        allNames = manageHDF5(cur)
        proj.saveProperty('store', proj.store)

        if (os.path.isfile('{}_genotype.DB'.format(self.proj.name)) and
                os.path.getsize('{}_genotype.DB'.format(self.proj.name)) > 0
           ) or 'snapshot_genotype.DB' in all_files:
            # if os.path.isfile('{}_genotype.DB'.format(self.proj.name)):
            #     os.remove('{}_genotype.DB'.format(self.proj.name))
            if 'snapshot_genotype.DB' in all_files:
                os.rename(
                    os.path.join(env.cache_dir, 'snapshot_genotype.DB'),
                    '{}_genotype.DB'.format(self.proj.name))
                all_files.remove('snapshot_genotype.DB')
            elif '{}_genotype.DB'.format(self.proj.name) in all_files:
                # an old version of snapshot saves $name.proj
                tempFile = os.path.join(env.cache_dir,
                                        '{}_genotype.DB'.format(self.proj.name))
                if os.path.isfile(tempFile) and os.path.isfile(
                        '{}_genotype.DB'.format(self.proj.name)):
                    os.remove('{}_genotype.DB'.format(self.proj.name))
                    os.rename(tempFile, '{}_genotype.DB'.format(self.proj.name))
                    os.remove(tempFile)
                all_files.remove('{}_genotype.DB'.format(self.proj.name))
            else:
                pass
            print(
                "Current storage mode is HDF5, transfrom genotype storage mode....."
            )
            jobs = 8

            validGenotypeFields = ['GT']
            if importer is not None:
                jobs = importer.jobs
                if len(importer.genotype_field) == 0:
                    validGenotypeFields = []
                for info in importer.genotype_info:
                    validGenotypeFields.append(info.name.replace("_geno", ""))
            # self.proj.db = DatabaseEngine()
            # self.proj.db.connect(self.proj.proj_file)

            if proj.build is not None:
                self.proj.build = proj.build
            else:
                self.proj.build = "hg19"
            IDs = self.proj.selectSampleByPhenotype("hdf5 is null")
            IDs = list(IDs)

            if len(IDs) == 0:
                os.remove('{}_genotype.DB'.format(self.proj.name))
            elif len(IDs) > 0:
                IDs.sort()
                if len(IDs) < 100:
                    jobs = 1

                reader = MultiVariantReader(
                    self.proj, "variant", "chr,pos,ref", "",
                    ['chr', 'pos', 'ref', 'vcf_variant(chr,pos,ref,alt,".")'],
                    validGenotypeFields, False, IDs, 4, True)
                reader.start()

                chunk = {
                    'variants/CHROM': [],
                    "variants/POS": [],
                    "variants/REF": [],
                    "variants/ALT": [],
                    "variants/ID": []
                }
                for geno_field in validGenotypeFields:
                    # geno_field=geno_field.replace("_geno","")
                    chunk["calldata/" + geno_field] = []

                for idx, raw_rec in enumerate(reader.records()):
                    # print(idx,len(raw_rec),raw_rec[3])
                    chunk["variants/REF"].append(raw_rec[0])
                    chunk["variants/ALT"].append([raw_rec[1]])
                    chunk["variants/CHROM"].append(raw_rec[2])
                    chunk["variants/POS"].append(raw_rec[3])
                    chunk["variants/ID"].append(raw_rec[5])
                    for genoID, geno_field in enumerate(validGenotypeFields):
                        colpos = [
                            6 + genoID + pos * len(validGenotypeFields)
                            for pos in range(len(IDs))
                        ]
                        # geno_field=geno_field.replace("_geno","")
                        chunk["calldata/" + geno_field].append([
                            raw_rec[pos] if raw_rec[pos] is not None else -10
                            for pos in colpos
                        ])
                for key, value in chunk.items():
                    dtype = np.dtype(object)
                    if key == "variants/POS":
                        dtype = np.dtype(np.int32)
                    elif key == "calldata/DP":
                        dtype = np.dtype(np.int16)
                    elif key == "calldata/GT":
                        dtype = np.dtype(np.int8)
                    elif key == "calldata/GQ":
                        dtype = np.dtype(np.float32)
                    elif "calldata" in key:
                        dtype = np.dtype(np.int16)
                        value = [
                            [0] if val[0] == '' else [val[0]] for val in value
                        ]
                    chunk[key] = np.array(value, dtype=dtype)

                importers = [None] * jobs
                task = None
                taskQueue = queue.Queue()
                workload = [int(float(len(IDs)) / jobs)] * jobs

                unallocated = max(0, len(IDs) - sum(workload))
                for i in range(unallocated):
                    workload[i % jobs] += 1
                numTasks = len(workload)
                variantIndex = self.proj.createVariantMap('variant', False)
                start_sample = 0
                # if IDs[0]==1:
                #     start_sample =0
                # else:
                #     start_sample=IDs[0]-1

                for job in range(numTasks):
                    # readQueue[job].put(chunk)
                    if workload[job] == 0:
                        continue

                    end_sample = min(start_sample + workload[job],
                                     start_sample + len(IDs))

                    if end_sample <= start_sample:
                        continue
                    HDFfile_Merge = "tmp_" + str(start_sample + 1) + "_" + str(
                        end_sample) + "_genotypes.h5"
                    names = [key for key in allNames.keys()]
                    updateSample(cur, start_sample, end_sample, IDs, names,
                                 allNames, HDFfile_Merge)

                    taskQueue.put(
                        HDF5GenotypeImportWorker(chunk, variantIndex,
                                                 start_sample, end_sample, IDs,
                                                 0, job, HDFfile_Merge,
                                                 validGenotypeFields,
                                                 self.proj.build))
                    start_sample = end_sample
                while taskQueue.qsize() > 0:
                    for i in range(jobs):
                        if importers[i] is None or not importers[i].is_alive():
                            task = taskQueue.get()
                            importers[i] = task
                            importers[i].start()
                            break
                for worker in importers:
                    if worker is not None:
                        worker.join()
        else:
            for hdf5file in hdf5files:
                os.remove(hdf5file)
            all_files = list(set(all_files))

            for file in all_files:
                if re.search(r'tmp(.*)h5', file):
                    os.rename(os.path.join(env.cache_dir, file), file)
                    all_files.remove(file)
        # self.proj.db.close()
        return all_files

    def addGenotypeToTar(self, tarfile, prog):
        hdf5files = glob.glob("tmp*h5")
        if len(hdf5files) > 0:
            for hdf5file in hdf5files:
                tarinfo = tarfile.gettarinfo(hdf5file, arcname=hdf5file)
                tarfile.addfile(tarinfo, ProgressFileObj(prog, hdf5file, 'rb'))
        return tarfile

    def remove_sample(self, IDs):
        for ID in IDs:
            hdf5file = self.get_sampleFileName(ID)
            storageEngine = Engine_Storage.choose_storage_engine(hdf5file)
            storageEngine.remove_sample(ID)

            storageEngine.close()

    def remove_variants_worker(self, storageEngine, variantIDs):
        storageEngine.remove_variants(variantIDs)

    def remove_variants(self, variantIDs, table):
        HDFfileNames = glob.glob("tmp*_genotypes.h5")
        storageEngines = []
        procs = []
        for HDFfileName in HDFfileNames:
            storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
            storageEngines.append(storageEngine)
            p = Process(
                target=self.remove_variants_worker,
                args=(storageEngine, variantIDs))
            procs.append(p)
            p.start()
        for proc in procs:
            proc.join()
        for storageEngine in storageEngines:
            storageEngine.close()
        env.logger.info('Removing table {} itself'.format(
            decodeTableName(table)))
        self.proj.db.removeTable(table)
        self.proj.db.commit()

    def remove_genotype_workder(self, HDFfileName, cond):
        try:
            storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
            storageEngine.remove_genotype(cond)
        finally:
            storageEngine.close()

    def remove_genotype(self, cond):
        procs = []
        for HDFfileName in glob.glob("tmp*genotypes.h5"):
            p = Process(
                target=self.remove_genotype_workder, args=(HDFfileName, cond))
            procs.append(p)
            p.start()
        for proc in procs:
            proc.join()

    def remove_genofields(self, IDs, items):
        HDFfileNames = glob.glob("tmp*_genotypes.h5")
        for HDFfileName in HDFfileNames:
            storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
            storageEngine.remove_genofields(items)
            storageEngine.close()

    def num_genoinfo(self, sampleID, expr, cond):
        HDFfileName = self.get_sampleFileName(sampleID)
        num = 0
        if HDFfileName is not None:
            accessEngine = Engine_Access.choose_access_engine(HDFfileName)
            num = accessEngine.num_genoinfo(sampleID, expr, cond)
            accessEngine.close()
        return num

    def num_genotypes(self, sampleID, cond, genotypes):
        HDFfileName = self.get_sampleFileName(sampleID)
        num = 0
        if HDFfileName is not None:
            accessEngine = Engine_Access.choose_access_engine(HDFfileName, read_only=True)
            num = accessEngine.num_genotypes(sampleID, cond, genotypes)
            accessEngine.close()
        return num

    def sum_genotypes(self, sampleID, cond, genotypes):
        HDFfileName = self.get_sampleFileName(sampleID)
        num = 0
        if HDFfileName is not None:
            accessEngine = Engine_Access.choose_access_engine(HDFfileName, read_only=True)
            num = accessEngine.sum_genotypes(sampleID, cond, genotypes)
            accessEngine.close()
        return num

    def num_variants(self, sampleID):
        HDFfileName = self.get_sampleFileName(sampleID)
        num = 0
        if HDFfileName is not None:
            accessEngine = Engine_Access.choose_access_engine(HDFfileName, read_only=True)
            num, _ = accessEngine.num_variants(sampleID)
            accessEngine.close()
        return num

    def geno_fields(self, sampleID):
        HDFfileName = self.get_sampleFileName(sampleID)
        genoFields = ""
        if HDFfileName is not None and HDFfileName not in self.cache:
            storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
            genoFields = storageEngine.geno_fields(sampleID)
            genoFields.sort()
            # genoFields=[x.lower() for x in genoFields]
            storageEngine.close()
            self.cache[HDFfileName] = genoFields
        return self.cache[HDFfileName]

    def geno_fields_nolower(self, sampleID):
        HDFfileName = self.get_sampleFileName(sampleID)
        genoFields = ""
        if HDFfileName is not None:
            storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
            genoFields = storageEngine.geno_fields(sampleID)
            genoFields.sort()
            genoFields = [x.upper() for x in genoFields]
            storageEngine.close()
        return genoFields

    def get_sampleFileName(self, id):

        self.cur = self.proj.db.cursor()
        self.cur.execute('SELECT HDF5 FROM sample WHERE sample_id = ?;', (id,))
        res = self.cur.fetchone()
        hdf5file = res[0]
        self.proj.db.commit()
        return hdf5file

    def get_typeOfColumn(self, sampleID, field):
        return "INTEGER"

    # def get_genoType_genoInfo(self,sampleID,genotypes,variant_table,fieldSelect):
    #     HDFfileName=self.get_sampleFileName(sampleID)
    #     accessEngine=Engine_Access.choose_access_engine(HDFfileName)
    #     # print(fieldSelect)
    #     result=accessEngine.get_geno_field_from_table(sampleID,genotypes,fieldSelect)
    #     return result

    def get_HDF5_sampleMap(self):
        self.cur = self.proj.db.cursor()
        self.cur.execute('SELECT sample_id, HDF5 FROM sample')
        result = self.cur.fetchall()
        sampleFileMap = {}
        for res in result:
            if res[1] not in sampleFileMap:
                sampleFileMap[res[1]] = []
            sampleFileMap[res[1]].append(res[0])
        return sampleFileMap

    def get_selected_variants(self, variant_table):
        self.cur = self.proj.db.cursor()
        self.cur.execute('SELECT variant_id FROM {}'.format(variant_table))
        result = self.cur.fetchall()
        variants = []
        for res in result:
            variants.append(res[0])
        return variants

    def get_noWT_variants_worker(self, queue, accessEngine, samples):
        # queue.put(accessEngine.get_noWT_variants(samples))
        # queue.put(accessEngine.get_genotype([],samples))
        queue.put(accessEngine.get_all_genotype(samples))

    def get_noWT_variants(self, samples, proj, where_clause, args):
        sampleFileMap = self.get_HDF5_sampleMap()
        queue = Queue()
        procs = []
        noWT = set()
        accessEngines = []
        for HDFfileName in glob.glob("tmp*genotypes.h5"):
            filename = HDFfileName.split("/")[-1]
            if filename in sampleFileMap:
                samplesInfile = sampleFileMap[filename]
                accessEngine = Engine_Access.choose_access_engine(HDFfileName)
                accessEngines.append(accessEngine)
                if len(list(set(samples).intersection(samplesInfile))) > 0:
                    p = Process(
                        target=self.get_noWT_variants_worker,
                        args=(queue, accessEngine,
                              list(set(samples).intersection(samplesInfile))))
                    procs.append(p)
                    p.start()
        for _ in procs:
            # result=queue.get()
            # noWT=noWT.union(set(result))
            for rownames, colnames, genoinfo in queue.get():
                # deal with empty geno
                if ((genoinfo == -1).sum()) == len(rownames) * len(colnames):
                    noWTvariants = rownames.tolist()
                else:
                    genoinfo[genoinfo == -1] = 0
                    rowsum = np.nansum(genoinfo, axis=1)
                    noWTvariants = rownames[np.where(rowsum > 0)].tolist()
                noWT = noWT.union(set(noWTvariants))
        for p in procs:
            p.join()
        for accessEngine in accessEngines:
            accessEngine.close()
        cur = proj.db.cursor()
        merged_table = '__variants_from_samples'
        query = 'CREATE TEMPORARY TABLE {} (variant_id INT);'.format(
            merged_table)
        cur.execute(query)
        divData = chunks(list(noWT))
        for chunk in divData:
            cur.execute('BEGIN TRANSACTION')
            for id in chunk:
                # query = 'INSERT INTO {}({}) VALUES {};'.format(merged_table,
                #         'variant_id', ",".join([ "("+str(id)+")" for id in variantIDs[:1]]))
                query = 'INSERT OR IGNORE INTO {}({}) VALUES {};'.format(
                    merged_table, 'variant_id', "(" + str(id) + ")")
                cur.execute(query)
        where_clause += ' AND ({}.variant_id IN (SELECT variant_id FROM {}))'.format(
            encodeTableName(args.from_table), merged_table)
        return where_clause

    def get_genoType_genoInfo_worker(self, queue, accessEngine, sampleNames,
                                     variants, validGenotypeFields, genotypes):
        queue.put(
            accessEngine.get_all_genotype_genoinfo(sampleNames, variants,
                                                   validGenotypeFields,
                                                   genotypes))

    def get_genoType_genoInfo(self, sampleDict, genotypes, variant_table,
                              genotypeFields, validGenotypeIndices,
                              validGenotypeFields, operations, fieldCalcs, prog,
                              prog_step):
        # print(genotypes)
        sampleFileMap = self.get_HDF5_sampleMap()
        #fieldSelect = list(sampleDict.values())[0][1]
        variants = []
        if variant_table != 'variant':
            variants = self.get_selected_variants(variant_table)

        master = {}
        queue = Queue()
        accessEngines = []
        procs = []

        for HDFfileName in glob.glob("tmp*genotypes.h5"):

            filename = HDFfileName.split("/")[-1]
            if filename in sampleFileMap:
                samplesInfile = sampleFileMap[filename]
                accessEngine = Engine_Access.choose_access_engine(HDFfileName)
                accessEngines.append(accessEngine)
                p = Process(
                    target=self.get_genoType_genoInfo_worker,
                    args=(
                        queue, accessEngine,
                        list(
                            set(sampleDict.keys()).intersection(samplesInfile)),
                        variants, validGenotypeFields, genotypes))
                procs.append(p)
                p.start()

        minPos = [i + 5 for i, x in enumerate(operations) if x == 2]
        maxPos = [i + 5 for i, x in enumerate(operations) if x == 3]

        for _ in procs:
            vardict = {}
            try:
                for rownames, colnames, sub_all in queue.get():
                    numrow = len(rownames)
                    genotype = sub_all[0]
                    variants = np.zeros(
                        shape=(numrow, len(validGenotypeFields) + 5),
                        dtype=np.int64)

                    variants[:, 3] = np.nansum(~np.isnan(genotype), axis=1)

                    variants[:, 0] = np.nansum(genotype == 1, axis=1)
                    variants[:, 1] = np.nansum(genotype == 2, axis=1)
                    variants[:, 2] = np.nansum(genotype == -1, axis=1)

                    for fieldpos, field in enumerate(validGenotypeFields):
                        genoinfo = sub_all[fieldpos + 1]
                        operation = operations[fieldpos]
                        if operation == 0:
                            variants[:, 5 + fieldpos] = np.nansum(
                                genoinfo, axis=1)
                        if operation == 1:
                            variants[:, 5 + fieldpos] = np.nansum(
                                genoinfo, axis=1)
                        if operation == 2:
                            variants[:, 5 + fieldpos] = np.nanmin(
                                genoinfo, axis=1)
                        if operation == 3:
                            variants[:, 5 + fieldpos] = np.nanmax(
                                genoinfo, axis=1)

                    vardict.update(dict(zip(rownames, variants)))
            except Exception as e:
                print(e)
                pass

            result = vardict
            for key, value in result.items():
                if key not in master:
                    master[key] = value
                else:
                    master[key] = [sum(x) for x in zip(master[key], value)]

                    if len(minPos) > 0 or len(maxPos) > 0:
                        for pos in minPos:
                            preValue = master[key][pos] - value[pos]
                            master[key][pos] = value[
                                pos] if preValue >= value[pos] else preValue
                        for pos in maxPos:
                            preValue = master[key][pos] - value[pos]
                            master[key][pos] = value[
                                pos] if preValue <= value[pos] else preValue
        for p in procs:
            p.join()
        for accessEngine in accessEngines:
            accessEngine.close()
        # prog.done()
        return master

    # def get_genoType_genoInfo_worker(self,queue,accessEngine,samples,variants,genotypes,fieldSelect,validGenotypeFields,operations):
    #     queue.put(accessEngine.get_geno_field_from_HDF5(samples,variants,genotypes,fieldSelect,validGenotypeFields,operations))

    # def get_genoType_genoInfo(self,sampleDict,genotypes,variant_table,genotypeFields,validGenotypeIndices,validGenotypeFields,operations,fieldCalcs,prog,prog_step):
    #     # print(genotypes)
    #     # print(validGenotypeFields)
    #     sampleFileMap=self.get_HDF5_sampleMap()
    #     fieldSelect=list(sampleDict.values())[0][1]
    #     variants=[]
    #     if variant_table != 'variant':
    #         variants=self.get_selected_variants(variant_table)

    #     master={}
    #     queue=Queue()
    #     procs=[]
    #     for HDFfileName in glob.glob("tmp*genotypes.h5"):
    #         filename=HDFfileName.split("/")[-1]
    #         if filename in sampleFileMap:
    #             samplesInfile=sampleFileMap[filename]
    #             accessEngine=Engine_Access.choose_access_engine(HDFfileName)
    #             p=Process(target=self.get_genoType_genoInfo_worker,args=(queue,accessEngine,list(set(sampleDict.keys()).intersection(samplesInfile)),variants,genotypes,fieldSelect,validGenotypeFields,operations))
    #             procs.append(p)
    #             p.start()

    #     minPos=[i+5 for i, x in enumerate(operations) if x == 2]
    #     maxPos=[i+5 for i, x in enumerate(operations) if x == 3]

    #     for _ in procs:

    #         result=queue.get()
    #         for key,value in result.items():
    #             if key not in master:
    #                 master[key]=value
    #             else:
    #                 master[key]= [sum(x) for x in zip(master[key], value)]
    #                 if len(minPos)>0 or len(maxPos)>0:
    #                     for pos in minPos:
    #                         preValue=master[key][pos]-value[pos]
    #                         master[key][pos]=value[pos] if preValue>=value[pos] else preValue
    #                     for pos in maxPos:
    #                         preValue=master[key][pos]-value[pos]
    #                         master[key][pos]=value[pos] if preValue<=value[pos] else preValue
    #     for p in procs:
    #         p.join()
    #     return master

    def get_Genotype(self, cur, table, proj, sample_ID):
        HDFfileName = self.get_sampleFileName(sample_ID)
        ids = []
        chrs = []
        cur.execute('SELECT variant_id FROM {0}'.format(encodeTableName(table)))
        for id in cur:
            ids.append(id[0])
        cur.execute('SELECT distinct(chr) from {0}'.format(
            encodeTableName(table)))
        for chr in cur:
            chrs.append(chr[0])
        accessEngine = Engine_Access.choose_access_engine(HDFfileName)
        g = set()
        NULL_to_0 = env.treat_missing_as_wildtype
        rownames, colnames, genoinfo = accessEngine.get_genotype(
            ids, [sample_ID], chrs)
        accessEngine.close()
        for idx, id in enumerate(rownames.tolist()):
            GT = genoinfo[idx][0]
            if np.isnan(GT):
                if NULL_to_0:
                    g.add((id, 0))
            else:
                g.add((id, GT))
        return g
        # return accessEngine.get_geno_by_variant_IDs_sample(ids,sample_ID,chrs)

    def validate_sex(self, proj, sample_ID, sex):
        HDFfileName = self.get_sampleFileName(sample_ID)
        accessEngine = Engine_Access.choose_access_engine(HDFfileName)
        geno = []
        invalidate = []
        if sex == 1:
            # geno=accessEngine.get_geno_by_sample_ID(sample_ID,"GT_geno",["X"])
            for rownames, colnames, genoinfo in accessEngine.get_all_genotype(
                [sample_ID], ["X"]):
                for idx, rowname in enumerate(rownames):
                    genotype = genoinfo[idx]
                    if np.isnan(genotype):
                        genotype = -1
                    geno.append([rowname, genotype])
            geno = np.array(geno)
            invalidate = geno[geno == 2]
        if sex == 2:
            # invalidate=accessEngine.get_geno_by_sample_ID(sample_ID,"GT_geno",["Y"])
            for rownames, colnames, genoinfo in accessEngine.get_all_genotype(
                [sample_ID], ["Y"]):
                for idx, rowname in enumerate(rownames):
                    genotype = genoinfo[idx]
                    if np.isnan(genotype):
                        genotype = -1
                    geno.append([rowname, genotype])
            geno = np.array(geno)
            invalidate = geno
        accessEngine.close()
        return len(invalidate)

    def importGenotypes(self, importer):

        from .importer_allele_hdf5 import importGenotypesInParallel
        return importGenotypesInParallel(importer)


def GenoStore(proj, importer=None):
    if proj.store == 'sqlite':
        return Sqlite_Store(proj)
    elif proj.store == 'hdf5':
        if importer is None or importer.format == "vcf":
            return HDF5_Store(proj)
        else:
            return Sqlite_Store(proj)
    else:
        raise RuntimeError('Unsupported genotype storage model {}'.format(
            proj.store))
