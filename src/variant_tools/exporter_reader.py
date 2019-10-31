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

import glob as glob
import os
import re
import sys
import time
from multiprocessing import Lock, Pipe, Process

import numpy as np

# from .geno_store import GenoStore
from .accessor import Engine_Access
from .merge_sort_parallel import index_HDF5_rowIDs
from .preprocessor import *
from .utils import (DatabaseEngine, ProgressBar, consolidateFieldName,
                    delayedAction, env)

MAX_COLUMN = 62
def VariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
        export_alt_build, IDs, jobs):
    if jobs == 0 and len(IDs) < MAX_COLUMN and proj.store=="sqlite":
        # using a single thread
        return EmbeddedVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    elif jobs > 0 and len(IDs) < MAX_COLUMN and proj.store=="sqlite":
        # using a standalone process to read things and
        # pass information using a pipe
        return StandaloneVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    else:
        # using multiple process to handle more than 1500 samples
        if len(IDs) // MAX_COLUMN + 2 > jobs:
            env.logger.info('Using {} processes to handle {} samples'.format(len(IDs) // MAX_COLUMN + 2, len(IDs)))
        return MultiVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, jobs, False)

class BaseVariantReader:
    def __init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        self.proj = proj
        self.table = table
        self.export_by_fields = export_by_fields
        self.order_by_fields = order_by_fields
        self.var_fields = var_fields
        self.geno_fields = geno_fields
        self.export_alt_build = export_alt_build
        self.IDs = IDs

    def start(self):
        pass

    def getQuery(self):
        select_clause, select_fields = consolidateFieldName(self.proj, self.table,
            ','.join(['variant.ref', 'variant.alt'] + self.var_fields),
            self.export_alt_build)
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in select_fields], [])
        #
        processed = set()
        # the normal annotation databases that are 'LEFT OUTER JOIN'
        where_conditions = []
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        # WHERE clause
        where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
        # GROUP BY clause
        if self.order_by_fields:
            order_fields, order_field_names = consolidateFieldName(self.proj, self.table, self.order_by_fields)
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ''
        #
        #
        tmp_fields = list(set(select_fields + (order_field_names if self.order_by_fields else [])))
        if 'variant.variant_id' in tmp_fields:
            tmp_fields.remove('variant.variant_id')
        from_clause = ('FROM (SELECT min(variant.variant_id) AS variant_variant_ID, '
            '{} {} {} GROUP BY variant.variant_id) AS _TMP'.format(
                ', '.join(['{} AS {}'.format(x, x.replace('.', '_')) for x in tmp_fields]),
                from_clause, where_clause))
        for fld in select_fields:
            select_clause = re.sub(fld.replace('.', '\.'), '_TMP.' + fld.replace('.', '_'),
                select_clause, flags=re.IGNORECASE)
        if self.order_by_fields:
            for fld in order_field_names:
                order_clause = re.sub(fld.replace('.', '\.'), '_TMP.' + fld.replace('.', '_'),
                    order_clause, flags=re.IGNORECASE)
        if self.geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in self.geno_fields:
                    if fld.lower() in header:
                        select_clause += ', {}_genotype.genotype_{}.{}'.format(self.proj.name, id, fld)
                    else:
                        select_clause += ', NULL'
                from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} ON {0}_genotype.genotype_{1}.variant_id = variant_variant_id '\
                    .format(self.proj.name, id, self.table)

        return 'SELECT {} {} {};'.format(select_clause, from_clause, order_clause)

#    def getQuery(self):
#        select_clause, fields = consolidateFieldName(self.proj, self.table,
#            ','.join(self.var_fields), self.export_alt_build)
#        if self.geno_fields:
#            for id in self.IDs:
#                header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
#                for fld in self.geno_fields:
#                    if fld.lower() in header:
#                        select_clause += ', {}_genotype.genotype_{}.{}'.format(self.proj.name, id, fld)
#                    else:
#                        select_clause += ', NULL'
#        # FROM clause
#        from_clause = 'FROM {} '.format(self.table)
#        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
#        #
#        processed = set()
#        # the normal annotation databases that are 'LEFT OUTER JOIN'
#        where_conditions = []
#        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
#            if (tbl.lower(), conn.lower()) not in processed:
#                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
#                processed.add((tbl.lower(), conn.lower()))
#        if self.geno_fields:
#            for id in self.IDs:
#                from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} ON {0}_genotype.genotype_{1}.variant_id = {2}.variant_id '\
#                    .format(self.proj.name, id, self.table)
#        # WHERE clause
#        where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
#        # GROUP BY clause
#        if self.order_by_fields:
#            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.order_by_fields)
#            order_clause = ' ORDER BY {}'.format(order_fields)
#        else:
#            order_clause = ''
#        return 'SELECT variant.ref,variant.alt,{} {} {} {};'.format(select_clause,
#            from_clause, where_clause, order_clause)


    def getVariantQuery(self):

        select_clause, select_fields = consolidateFieldName(self.proj, self.table,
            ','.join(['variant_id', 'variant.ref', 'variant.alt'] + self.var_fields), self.export_alt_build)
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in select_fields], [])
        #
        processed = set()
        where_conditions = []
        # the normal annotation databases that are 'LEFT OUTER JOIN'
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        # WHERE clause
        where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
        # GROUP BY clause
        if self.order_by_fields:
            order_fields, order_field_names = consolidateFieldName(self.proj, self.table, self.order_by_fields + ',variant_id')
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            #self.order_by_fields = '{}.variant_id'.format(self.table)
            #order_field_names = ['{}.variant_id'.format(self.table)]
            #order_clause = ' ORDER BY {}.variant_id'.format(self.table)
            order_clause = ''
        #
        tmp_fields = list(set(select_fields + (order_field_names if self.order_by_fields else [])))
        if 'variant.variant_id' in tmp_fields:
            tmp_fields.remove('variant.variant_id')
        from_clause = ('FROM (SELECT min(variant.variant_id) AS variant_variant_ID, '
            '{} {} {} GROUP BY variant.variant_id) AS _TMP'.format(
                ', '.join(['{} AS {}'.format(x, x.replace('.', '_')) for x in tmp_fields]),
                from_clause, where_clause))
        for fld in select_fields:
            select_clause = re.sub(fld.replace('.', '\.'), '_TMP.' + fld.replace('.', '_'),
                select_clause, flags=re.IGNORECASE)
        if self.order_by_fields:
            for fld in order_field_names:
                order_clause = re.sub(fld.replace('.', '\.'), '_TMP.' + fld.replace('.', '_'),
                    order_clause, flags=re.IGNORECASE)
        #env.logger.error('SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause))
        return 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)

    def getSampleQuery(self, IDs):
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            'variant_id', False)
        for id in IDs:
            header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
            for fld in self.geno_fields:
                if fld.lower() in header:
                    select_clause += ', g{}.{}'.format(id, fld)
                else:
                    select_clause += ', NULL'
        # FROM clause
        from_clause = 'FROM {0}'.format(self.table)
        if self.table != 'variant':
            from_clause += ' LEFT OUTER JOIN variant on {0}.variant_id = variant.variant_id '.format(self.table)
        for id in IDs:
            from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} g{1} ON g{1}.variant_id = {2}.variant_id '\
                .format(self.proj.name, id, self.table)
        # WHERE clause
        where_clause = ''
        # GROUP BY clause
        if self.order_by_fields:
            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.order_by_fields + ', variant_id')
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ' ORDER BY {}.variant_id'.format(self.table)

        return 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)




class EmbeddedVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        self.proj = proj
        self.var_fields = var_fields
        BaseVariantReader.__init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)

    def records(self):
        env.logger.debug('Running query {}'.format(self.getQuery()))
        cur = self.proj.db.cursor()
        try:
            cur.execute(self.getQuery())
        except Exception as e:
            raise ValueError('Failed to generate output: {}\nIf your project misses one of the following fields {}, you might want to add them to the project (vtools update TABLE INPUT_FILE --var_info FIELDS) or stop exporting them using format parameters (if allowed).'\
                .format(e, ', '.join(self.var_fields)))
        for rec in cur:
            # the first two items are always ref and alt
            yield rec


class StandaloneVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs):
        BaseVariantReader.__init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)
        self.proj = proj
        ID_needed_idx = [id for id in IDs if not self.proj.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id))]
        if len(ID_needed_idx) > 0:
            prog = ProgressBar('Creating indexes', len(ID_needed_idx))
            cur = self.proj.db.cursor()
            for idx, id in enumerate(ID_needed_idx):
                if not self.proj.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id)):
                    cur.execute('CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'.format(self.proj.name, id))
                prog.update(idx)
            prog.done()
        self.var_fields = var_fields
        self.reader, w = Pipe(False)
        self.worker = VariantWorker(proj.name, proj.annoDB, self.getQuery(), w, None)
        self.worker.start()

    def start(self):
        # the first None, indicating ready to output
        with delayedAction(env.logger.info, 'Selecting genotypes...'):
            self.reader.recv()

    def records(self):
        while True:
            rec = self.reader.recv()
            if rec is None:
                break
            else:
                yield rec

class MultiVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, jobs,transformToHDF5):
        BaseVariantReader.__init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build,  IDs)
        self.proj = proj
        self.var_fields = var_fields
        # the first job for variants
        r, w = Pipe(False)
        lock = Lock()

        p = VariantWorker(proj.name, proj.annoDB, self.getVariantQuery(), w, lock)
        self.workers = [p]
        self.readers = [r]
        IDs = list(IDs)

        # we may need more jobs due to the limit of max columns
        # but we will only have self.jobs active jobs
        jobs = min(1, jobs)
        self.jobs = max(jobs, len(IDs) // MAX_COLUMN + 2)
        block = len(IDs) // (self.jobs-1) + 1


        if self.proj.store=="sqlite" or transformToHDF5:

            if transformToHDF5:
                self.proj.db.attach('{}_genotype.DB'.format(self.proj.name), '{}_genotype'.format(self.proj.name), lock=lock)

            with delayedAction(env.logger.info, 'Checking indexes'):
                ID_needed_idx = [id for id in IDs if not self.proj.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id))]
            if len(ID_needed_idx) > 0:
                prog = ProgressBar('Creating indexes', len(ID_needed_idx))
                cur = self.proj.db.cursor()
                for idx, id in enumerate(ID_needed_idx):
                    if not self.proj.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id)):
                        cur.execute('CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'.format(self.proj.name, id))
                    prog.update(idx)
                prog.done()
            for i in range(self.jobs - 1):
                r, w = Pipe(False)
                subIDs = IDs[(block*i):(block *(i + 1))]

                p = VariantWorker(proj.name, proj.annoDB, self.getSampleQuery(subIDs), w, lock)

                self.workers.append(p)
                self.readers.append(r)
        elif self.proj.store=="hdf5":
            # store = GenoStore(proj)
            # sampleFileMap=store.get_HDF5_sampleMap()
            cur = self.proj.db.cursor()
            cur.execute('SELECT sample_id, HDF5 FROM sample')
            result=cur.fetchall()
            sampleFileMap={}
            for res in result:
                if res[1] not in sampleFileMap:
                    sampleFileMap[res[1]]=[]
                sampleFileMap[res[1]].append(res[0])
            samplefiles=glob.glob("tmp*genotypes.h5")
            samplefiles.sort(key=lambda x:int(x.split("_")[1]))
            # self.jobs=len(samplefiles)+1

            self.jobs=1
            cur.execute('SELECT value FROM project WHERE name="multiVCF";')
            multiVCF=cur.fetchall()
            if len(multiVCF)==0:
                multiVCF=0
            else:
                multiVCF=int(multiVCF[0][0])
            for HDFfileName in samplefiles:
                filename=HDFfileName.split("/")[-1]
                if filename in sampleFileMap:
                    samplesInfile=sampleFileMap[filename]
                    overlapSamples=list(set(IDs).intersection(samplesInfile))
                    if len(overlapSamples)>0:
                        r, w = Pipe(False)
                        self.jobs+=1
                        if multiVCF==0:
                            p=VariantWorker_HDF5(HDFfileName,overlapSamples,self.geno_fields,w)
                        else:
                            p=VariantWorker_HDF5_multi(proj.name, proj.annoDB, self.getVariantQuery(),HDFfileName,overlapSamples,self.geno_fields,lock,w)
                        self.workers.append(p)
                        self.readers.append(r)




    def start(self):
        prog = ProgressBar('Selecting genotypes', len(self.readers))
        status = [0] * len(self.readers)
        while True:
            for idx, (w,r) in enumerate(zip(self.workers, self.readers)):

                if status[idx] == 2: # completed
                    continue
                elif status[idx] == 1: # started?
                    if r.poll():
                        r.recv()   # should have None
                        status[idx] = 2
                        prog.update(sum([x==2 for x in status]))
                else:    # not started
                    # start another process
                    if sum([x == 1 for x in status]) < self.jobs:
                        status[idx] = 1
                        w.start()
            if sum(status) == 2 * len(self.readers):
                break
            else:
                time.sleep(1)
        prog.done()

    def records(self):
        #
        rec = []
        id = None
        last = len(self.readers) - 1
        all_done = False
        if self.table=="variant":
            while True:
                try:
                    for idx, reader in enumerate(self.readers):
                        val = reader.recv()
                        if val is None:
                            all_done = True
                            break
                        sys.stdout.flush()
                        if idx == 0:
                            id = val[0]
                        elif id != val[0]:
                            raise ValueError('Read different IDs from multiple processes')
                        rec.extend(val[1:])
                        if idx == last:
                            yield rec
                            rec = []
                    if all_done:
                        break
                except Exception as e:
                    env.logger.debug('Failed to get record: {}'.format(e))
            for p in self.workers:
                p.terminate()
        else:

            while True:
                try:
                    val=self.readers[0].recv()
                    if val is None:
                        all_done=True
                        break
                    # sys.stdou.flush()
                    id=val[0]
                    rec.extend(val[1:])
                    notSameID=True
                    while notSameID:
                        for idx,reader in enumerate(self.readers[1:]):
                            val=reader.recv()
                            if id ==val[0]:
                                rec.extend(val[1:])
                                if idx+1==last:
                                    notSameID=False

                                    yield rec
                                    rec=[]
                                    break
                    if all_done:
                        break
                except Exception as e:
                    env.logger.debug('Failed to get record: {}'.format(e))
            for p in self.workers:
                p.terminate()


class VariantWorker_HDF5_multi(Process):
    # this class starts a process and used passed query to read variants
    def __init__(self,dbname, annoDB, query,HDFfileName,samples, geno_fields,lock,output):
        self.dbname = dbname
        self.annoDB = annoDB
        self.query = query
        self.lock = lock
        self.fileName=HDFfileName
        self.samples=samples
        self.geno_fields=geno_fields
        self.output = output
        Process.__init__(self)

    def run(self):
        accessEngine=Engine_Access.choose_access_engine(self.fileName)
        sortedID=index_HDF5_rowIDs(self.fileName)
        genoinfo_fields=[field.replace("_geno","") for field in self.geno_fields]
        if "GT" in genoinfo_fields:
            genoinfo_fields.remove("GT")
        db = DatabaseEngine()
        db.connect(self.dbname + '.proj', readonly=True, lock=self.lock)

        for anno in self.annoDB:
            db.attach(os.path.join(anno.dir, anno.filename), lock=self.lock)
        cur = db.cursor()
        cur.execute(self.query)
        self.output.send(None)
        last_id = None
        for rec in cur:
            if rec[0] != last_id:
                last_id = rec[0]
                chr=rec[3]
                # val=accessEngine.get_geno_by_row_pos(rec[0],chr,sortedID,self.samples,genoinfo_fields)
                sub_all=accessEngine.get_geno_by_row_pos(rec[0],chr,sortedID,self.samples,genoinfo_fields)
                genoType=sub_all[0][0]
                val=None
                # numrow,numcol=genoType.shape[0],genoType.shape[1]
                numcol=len(genoType)
                if len(self.geno_fields)==0:
                    val=genoType
                else:
                    val=np.zeros(shape=(numcol*(len(self.geno_fields))),dtype=float)
                    for col in range(numcol):
                        val[col*len(self.geno_fields)]=genoType[col]
                    if len(genoinfo_fields)>0:
                        for pos,field in enumerate(genoinfo_fields):
                            genoInfo=sub_all[pos+1][0]
                            for col in range(numcol):
                                val[col*len(self.geno_fields)+pos+1]=genoInfo[col]

                val=np.where(np.isnan(val), None, val)
                val=[ int(i) if i is not None else i for i in val]
                self.output.send([rec[0]]+val)
                # val=accessEngine.get_geno_by_row_pos(rec[0],chr,sortedID,self.samples,genoinfo_fields)
                # val=np.where(np.isnan(val), None, val)
                # self.output.send([rec[0]]+val.tolist())
        self.output.send(None)
        db.close()






class VariantWorker_HDF5(Process):
    # this class starts a process and used passed query to read variants
    def __init__(self,HDFfileName,samples, geno_fields,output):
        self.fileName=HDFfileName
        self.samples=samples
        self.geno_fields=geno_fields
        self.output = output
        Process.__init__(self)

    def run(self):
        accessEngine=Engine_Access.choose_access_engine(self.fileName)
        vardict={}
        genoinfo_fields=[field.replace("_geno","") for field in self.geno_fields]
        if "GT" in genoinfo_fields:
            genoinfo_fields.remove("GT")
        for rownames,colnames,sub_all in accessEngine.get_all_genotype_genoinfo(self.samples,[],genoinfo_fields):

            genoType=sub_all[0]
            numrow,numcol=genoType.shape[0],genoType.shape[1]
            if len(self.geno_fields)==0:
                info=np.zeros(shape=(numrow,numcol),dtype=int)
                for col in range(numcol):
                    info[:,col]=genoType[:,col]
            else:
                info=np.zeros(shape=(numrow,numcol*(len(self.geno_fields)+1)),dtype=float)
                for col in range(numcol):
                    info[:,col]=genoType[:,col]
                if len(genoinfo_fields)>0:
                    for pos,field in enumerate(genoinfo_fields):
                        genoInfo=sub_all[pos+1]
                        for col in range(numcol):
                            info[:,col*len(self.geno_fields)+pos+1]=genoInfo[:,col]
            vardict.update(dict(zip(rownames,info)))

        # result=accessEngine.get_genoType_forExport_from_HDF5(self.samples,self.geno_fields)
        self.output.send(None)
        last_id = None
        for key,val in vardict.items():
            if key!=last_id:
                last_id=key
                val=np.where(np.isnan(val), None, val)
                self.output.send([key]+val.tolist())
        self.output.send(None)




class VariantWorker(Process):
    # this class starts a process and used passed query to read variants
    def __init__(self, dbname, annoDB, query, output, lock):
        self.dbname = dbname
        self.annoDB = annoDB
        self.query = query
        self.output = output
        self.lock = lock
        Process.__init__(self)

    def run(self):
        db = DatabaseEngine()
        db.connect(self.dbname + '.proj', readonly=True, lock=self.lock)
        if '{}_genotype.'.format(self.dbname) in self.query:
            db.attach('{}_genotype.DB'.format(self.dbname), '{}_genotype'.format(self.dbname), lock=self.lock)
        for anno in self.annoDB:
            db.attach(os.path.join(anno.dir, anno.filename), lock=self.lock)
        cur = db.cursor()
        cur.execute(self.query)
        # reporting to the main process that SQL query is done
        self.output.send(None)
        last_id = None
        for rec in cur:
            if rec[0] != last_id:
                last_id = rec[0]

                self.output.send(rec)

        self.output.send(None)
        db.close()
