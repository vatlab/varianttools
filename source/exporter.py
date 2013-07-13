#!/usr/bin/env python
#
# $File: exporter.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org)
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

import os
import sys
import gzip
import re
import time
from multiprocessing import Process, Pipe, Lock
if sys.version_info.major == 2:
    from itertools import izip, repeat
else:
    izip = zip
    from itertools import repeat
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, normalizeVariant, \
    consolidateFieldName, DatabaseEngine, env, encodeTableName, decodeTableName

 
class JoinFields:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        try:
            if type(item) == str:
                return item
            else:
                return self.sep.join(item)
        except:
            return str(item)

class IfMulti:
    def __init__(self, ifFunc=None, elseFunc=None):
        if hasattr(ifFunc, '__call__'):
            self.ifFunc = ifFunc.__call__
        else:
            self.ifFunc = ifFunc
        if hasattr(elseFunc, '__call__'):
            self.elseFunc = elseFunc.__call__
        else:
            self.elseFunc = elseFunc

    def __call__(self, item):
        if type(item) == tuple:
            return item[0] if self.ifFunc is None else self.ifFunc(item)
        else:
            return item if self.elseFunc is None else self.elseFunc(item)
        
class JoinRecords:
    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep
    
    def __call__(self, item):
        try:
            if type(item) == tuple or type(item) == list:
                return self.sep.join([str(x) for x in item])
            else:
                return str(item)
        except:
            return str(item)

class ValueOfNull:
    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return self.val if item in ('', None) else item


class Formatter:
    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, item):
        try:
            return self.fmt.format(item)
        except:
            return str(item)

class CSVFormatter:
    def __init__(self):
        pass

    def __call__(self, item):
        if type(item) == str:
            if not item:
                return ''
            elif '"' in item:
                return '"' + item.replace('"', '""') + '"'
            # quote all strings, because sometimes excel will treat them differently.
            else:
                return '"' + item + '"'
        else:
            # not string...
            val = str(item)
            if '"' in val:
                return '"' + val.replace('"', '""') + '"'
            if ',' in val or '\n' in val:
                return '"' + val + '"'
            return val

rec_ref = '-'
rec_alt = '-'

class GenoFormatter:
    # representation for missing value is style dependent,
    # the default value None will cause each style to use its default value.
    def __init__(self, style='genotype', sep='\t', null='-', base=0):
        self.sep = sep
        self.null = null
        if style == 'numeric':
            self.missing = 'NA'
        elif style == 'genotype':
            self.missing = '.'
        elif style == 'plink':
            # PED format seems to use ACTG, and 0 for missing.
            # see http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html
            self.missing = '0'
        elif style :
            self.missing = ''
        self.base = base
        #
        self.vcf_map = {
            # When we output a single sample, it can have genotype
            #   0, 1, 2, and (-1, -1)
            # as we imported from .vcf file.
            #
            0: '0/0',
            1: '0/1',
            2: '1/1',
            (-1,-1): '1/2',
            #
            # if one of the two alternative variants is filtered out, we can have
            # a single -1 variant. Note that we cannot yet import partial missing
            # data.
            -1: './1',
            #
            # when two or more samples are outputted, one sample might not have
            # any genotype for a variant, or more than one variant. There can 
            # be more than one variant at a location though.
            #
            None: '.',
            (None, None): '.',
            (None, None, None): '.',
            # 
            # Having two valid and complete genotypes for both variants is
            # as far as I can imagine, not possible. However, because it is
            # true that the variant is homogenous wildtype at both variants,
            # (0,0) is listed here.
            (0,0): '0/0',
            #
            #(0,1): '0/1',
            #(0,2): '2/2',
            #
            # However, if there are two variants, one sample can have 0, 1, 2, and
            # the None for another variant (because that variant exists in other
            # samples). Note that the second item is for another variant
            (0, None): '0/0',
            (None,0):  '0/0',
            (1, None): '0/1',
            (None, 1): '0/2',
            (2,None):  '1/1', 
            (None, 2): '2/2',
            # the single -1 case 
            (-1, None): './1',
            (None, -1): './2',
            #
            # the same goes to the case with three alternative alleles
            (0, None, None): '0/0',
            (None, 0, None): '0/0',
            (None, None, 0): '0/0', 
            (1, None, None): '0/1',
            (None, 1, None): '0/2',
            (None, None, 1): '0/3', 
            (2, None, None): '1/1',
            (None, 2, None): '2/2',
            (None, None, 2): '3/3',
            # the -1 case is more complicated because there can be one or two -1
            (-1, None, None): './1',
            (None, -1, None): './2',
            (None, None, -1): './3',
            (-1, -1, None):   '1/2',
            (-1, None, -1):   '1/3',
            (None, -1, -1):   '2/3',
        }
        #
        self.numeric_map = {
            # number of non-wildtype alleles
            # FIXME: We are treating multiple alternative alleles as the same
            # non-wildtype allele, which might be wrong for some formats.
            0: 0,
            1: 1,
            2: 2,
            (-1, -1): 2,  # two DIFFERENT alternative alleles
            #
            -1: 1,
            #
            None: None,
            (None, None): None,
            (None, None, None): None,
            # 
            (0,0): 0,
            #
            (0, None): 0,
            (None,0):  0,
            (1, None): 1,
            (None, 1): 1,
            (2,None):  2, 
            (None, 2): 2,
            # 
            (-1, None): 1,
            (None, -1): 1,
            #
            (0, None, None): 0,
            (None, 0, None): 0,
            (None, None, 0): 0,
            (1, None, None): 1,
            (None, 1, None): 1,
            (None, None, 1): 1,
            (2, None, None): 2,
            (None, 2, None): 2,
            (None, None, 2): 2,
            #
            (-1, None, None): 1,
            (None, -1, None): 1,
            (None, None, -1): 1,
            (-1, -1, None):   2,
            (-1, None, -1):   2,
            (None, -1, -1):   2,
        }
        #
        if style == 'genotype':
            self.__call__ = self.fmt_genotype
        elif style == 'numeric':
            self.__call__ = self.fmt_numeric
        elif style == 'vcf':
            self.__call__ = self.fmt_vcf
        elif style == 'plink':
            self.__call__ = self.fmt_plink
        else:
            raise ValueError('Only genotype and numeric styles are allowed')

    def fmt_genotype(self, item):
        global rec_ref, rec_alt
        if type(item) == int:
            # single genotype case
            ref = self.null if rec_ref == '-' else rec_ref
            alt = self.null if rec_alt == '-' else rec_alt
            #
            # 0, 1, 2, -1
            if item == 0:
                return ref + self.sep + ref
            elif item == 1:
                return ref + self.sep + alt
            elif item == 2:
                return alt + self.sep + alt
            else: 
                return alt + self.sep + self.null
        elif type(item) == tuple:
            # Two aleternative alleles
            if item == (-1, -1):
                return (self.null if rec_alt[0] == '-' else rec_alt[0]) + self.sep + (self.null if rec_alt[1] == '-' else rec_alt[1])
            elif len(item) > 1 and item.count(item[0]) == len(item):
                # assume duplicate entry caused by annotation database
                ref = self.null if rec_ref[0] == '-' else rec_ref[0]
                alt = self.null if rec_alt[0] == '-' else rec_alt[0]
                if item[0] == 1:
                    return ref + self.sep + alt
                elif item[0] == 2:
                    return alt + self.sep + alt
                else:
                    return ref + self.sep + ref
            #
            # the cases for (None, None) and (None, None, None)
            elif all([x is None for x in item]):
                return self.missing + self.sep + self.missing
            else:
                raise ValueError('Failed to export genotype {} with ref {} and alt {}'.format(item, rec_ref, rec_alt))
        elif item is None:
            return self.missing + self.sep + self.missing
        else:
            raise ValueError('Failed to export genotype {}'.format(item))

    def fmt_plink(self, item):
        global rec_ref, rec_alt
        if type(item) == int:
            # single genotype case
            ref = self.null if rec_ref == '-' else rec_ref
            alt = self.null if rec_alt == '-' else rec_alt
            #
            # 0, 1, 2, -1
            if item == 0:
                return ref + self.sep + ref
            elif item == 1:
                return ref + self.sep + alt
            elif item == 2:
                return alt + self.sep + alt
            else:  # (./A), considered as missing
                return self.missing + self.sep + self.missing
        # multi-allele case, ignore
        elif type(item) == tuple:
            raise ValueError('plink format cannot handle multiple alleles {}'.format(item))
        elif item is None:
            return self.missing + self.sep + self.missing
        else:
            raise ValueError('Failed to export genotype {}'.format(item))


    def fmt_numeric(self, item):
        try:
            cnt = self.numeric_map[item]
            if cnt is None:
                return self.missing
            else:
                return cnt + self.base
        except:
            raise ValueError('Failed to export genotype {} in numeric style'.format(item))

    def fmt_vcf(self, item):
        try:
            return self.vcf_map[item]
        except:
            raise ValueError('Failed to export genotype {} in vcf style with ref {} and alt {}.'.format(item, rec_ref, rec_alt))

class Constant:
    def __init__(self, val=''):
        self.val = val

    def __call__(self, item):
        return self.val

class SequentialCollector:
    def __init__(self, extractors):
        '''Define an extractor that calls a list of extractors. The string extracted from
        the first extractor will be passed to the second, and so on.'''
        self.extractors = []
        for e in extractors:
            if hasattr(e, '__call__'):
                self.extractors.append(e.__call__)
            else:
                self.extractors.append(e)

    def __call__(self, item):
        for e in self.extractors:
            # if multiple records are returned, apply to each of them
            if type(item) is tuple:
                if type(item[0]) is tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item

MAX_COLUMN = 62
def VariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
        export_alt_build, IDs, jobs):
    if jobs == 0 and len(IDs) < MAX_COLUMN:
        # using a single thread
        return EmbeddedVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    elif jobs > 0 and len(IDs) < MAX_COLUMN:
        # using a standalone process to read things and
        # pass information using a pipe
        return StandaloneVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs)
    else:
        # using multiple process to handle more than 1500 samples
        if len(IDs) // MAX_COLUMN + 2 > jobs:
            env.logger.info('Using {} processes to handle {} samples'.format(len(IDs) // MAX_COLUMN + 2, len(IDs)))
        return MultiVariantReader(proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, jobs)

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
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(self.var_fields), self.export_alt_build)
        if self.geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in self.geno_fields:
                    if fld.lower() in header:
                        select_clause += ', {}_genotype.genotype_{}.{}'.format(self.proj.name, id, fld)
                    else:
                        select_clause += ', NULL'
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        # the normal annotation databases that are 'LEFT OUTER JOIN'
        where_conditions = []
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        if self.geno_fields:
            for id in self.IDs:
                from_clause += ' LEFT OUTER JOIN {0}_genotype.genotype_{1} ON {0}_genotype.genotype_{1}.variant_id = {2}.variant_id '\
                    .format(self.proj.name, id, self.table)
        # WHERE clause
        where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
        # GROUP BY clause
        if self.order_by_fields:
            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.order_by_fields)
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ''
        return 'SELECT variant.ref,variant.alt,{} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)

    def getVariantQuery(self):
        select_clause, fields = consolidateFieldName(self.proj, self.table,
            ','.join(['variant_id', 'variant.ref', 'variant.alt'] + self.var_fields), self.export_alt_build)
        # FROM clause
        from_clause = 'FROM {} '.format(self.table)
        fields_info = sum([self.proj.linkFieldToTable(x, self.table) for x in fields], [])
        #
        processed = set()
        where_conditions = []
        # the normal annotation databases that are 'LEFT OUTER JOIN'
        for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
            if (tbl.lower(), conn.lower()) not in processed and '.__' not in tbl:
                from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
                processed.add((tbl.lower(), conn.lower()))
        # temporary connection tables are appended as WHERE clause.
        #for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
        #    if (tbl.lower(), conn.lower()) not in processed and '.__' in tbl:
        #        from_clause += ' , {}'.format(tbl)
        #        where_conditions.append(conn)
        #        processed.add((tbl.lower(), conn.lower()))
        # WHERE clause
        where_clause = 'WHERE {}'.format(' AND '.join(['({})'.format(x) for x in where_conditions])) if where_conditions else ''
        # GROUP BY clause
        if self.order_by_fields:
            order_fields, tmp = consolidateFieldName(self.proj, self.table, self.order_by_fields + ',variant_id')
            order_clause = ' ORDER BY {}'.format(order_fields)
        else:
            order_clause = ' ORDER BY {}.variant_id'.format(self.table)
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
        processed = set()
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
        #print 'SELECT {} {} {} {};'.format(select_clause, from_clause, where_clause, order_clause)
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
        s = delayedAction(env.logger.info, 'Selecting genotypes...')
        self.reader.recv()
        del s
        
    def records(self):
        while True:
            rec = self.reader.recv()
            if rec is None:
                break
            else:
                yield rec

class MultiVariantReader(BaseVariantReader):
    def __init__(self, proj, table, export_by_fields, order_by_fields, var_fields, geno_fields,
            export_alt_build, IDs, jobs):
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
        #
        s = delayedAction(env.logger.info, 'Checking indexes')
        ID_needed_idx = [id for id in IDs if not self.proj.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id))]
        del s
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


class Exporter:
    '''A general class for importing variants'''
    def __init__(self, proj, table, filename, samples, format, build, header,
        jobs, fmt_args):
        self.proj = proj
        self.db = proj.db
        self.jobs = jobs
        #
        # table
        self.table = encodeTableName(table)
        if not self.proj.isVariantTable(self.table):
            raise ValueError('{} is not a valid variant table'.format(decodeTableName(self.table)))
        #
        # filename
        self.filename = filename
        #
        # samples
        self.IDs = self.proj.selectSampleByPhenotype(samples) if samples else []
        self.samples = []
        if samples:
            env.logger.info('Genotypes of {} samples are exported.'.format(len(self.IDs)))
            cur = self.db.cursor()
            for ID in self.IDs:
                cur.execute('SELECT filename, sample_name FROM sample, filename WHERE sample.file_id = filename.file_id AND sample.sample_id = {};'\
                    .format(self.db.PH), (ID,))
                for rec in cur:
                    self.samples.append('{}'.format(rec[1]))
                    env.logger.debug('\t'.join(['{}'.format(x) for x in rec]))
        # 
        # build
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                if self.proj.alt_build is not None:
                    env.logger.info('Using primary reference genome {} of the project.'.format(self.build))
            else:
                raise ValueError('This project does not have any data to export.')
            self.export_alt_build = False
        elif build in [self.proj.build, self.proj.alt_build]:
            self.build = build
            self.export_alt_build = build == self.proj.alt_build
        else:
            raise ValueError('Invalid parameter for --build, which should be either the primary ({}) or the alternative ({}) reference genome of the project.'\
                .format(self.proj.build, self.proj.alt_build))
        #
        # format
        if not format:
            filename = self.filename.lower() if self.filename is not None else ''
            if filename.lower().endswith('.vcf'):
                format = 'vcf'
            elif filename.lower().endswith('.csv'):
                format = 'csv'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            self.format = fileFMT(format, fmt_args)
        except Exception as e:
            env.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        if not self.format.columns:
            raise ValueError('Cannot output in format {} because no output column is defined for this format.'.format(format)) 
        #
        # header
        if header is None:
            self.header = ''
        elif header == ['-']:
            # read from standard input
            env.logger.info('Reading header from standard input')
            self.header = sys.stdin.read()
        else:
            self.header = self.format.delimiter.join(header)

    def getAdjFunc(self, code):
        if not code:
            return None
        e = eval(code)
        if hasattr(e, '__iter__'):
            # if there are multiple functors, use a sequence functor
            e = SequentialCollector(e)
        if hasattr(e, '__call__'):
            e = e.__call__
        return e

    def getFormatters(self, indexes, fields):
        #
        # indexes: indexes of data from the input data (returned by sql)
        # fields: field names
        #
        # formatter... : key: value on how to format field(s)
        #
        # Return:
        #    formatters to process all values
        #
        # get first keys for self.formatter
        fmt_first_keys = [x.lower() if ',' in x else x.split(',')[0].lower() for x in self.format.formatter.keys()]
        fmt_keys = [x.lower() for x in self.format.formatter.keys()]
        #
        formatters = []
        i = 0
        while i < len(indexes):
            if fields[i].lower() in fmt_keys:
                # use adj to handle value at index[i]
                formatters.append((self.getAdjFunc(self.format.formatter[fields[i].lower()]), indexes[i]))
                i += 1
            elif fields[i].lower() in fmt_first_keys:
                # start of a multi-field entry
                found = False
                for key,length in [(x,len(x)) for x in fmt_keys if type(x) != str]:
                    if ','.join(fields[i:i+length]).lower() == key.lower():
                        found = True
                        formatters.append(self.getAdjFunc(self.format.formatter[fields[i]]),
                            [indexes[j] for j in range(i, i+length)])
                        i += length
                if not found:
                    formatters.append((None, indexes[i]))
                    i += 1
            else:
                formatters.append((None, indexes[i]))
                i += 1
        return formatters

    def exportTfam(self, fname):
        fname = fname.replace('$table', self.table)
        fname += '.tfam' if not fname.endswith('.tfam') else ''
        if os.path.exists(fname):
            os.remove(fname)
        env.logger.info('Sample names are exported to {}'.format(fname))
        # a tfam file is the first 6 columns of a ped file
        # FID, ID, paternal ID, maternal ID, sex, phenotype
        # Will output ID only; other fields will be populated with placeholders
        with open(fname, 'a') as f:
            for i in range(len(self.samples)):
                f.write('\t'.join([str(i), self.samples[i], '0', '0', '0', '-9']) + '\n')
        return True
    
    def exportData(self):
        '''Export data in specified format'''
        if self.format.additional_exports:
            additionalFiles = [x.strip() for x in self.format.additional_exports.split(',')]
            for item in additionalFiles:
                base, ext = os.path.splitext(item)
                if len(ext) == 0 or len(base) == 0:
                    raise ValueError('Invalid filename specification "{}".'.format(item))
                try:
                    eval('{}'.format('self.export{}(base)'.format(ext[1:].capitalize())))
                except Exception as e:
                    raise ValueError('Additional export to file *{} is not supported'.format(ext))
        else:
            pass
        # get all fields
        var_fields = [x.strip() for x in self.format.export_by_fields.split(',')] if self.format.export_by_fields else []
        geno_fields = []
        for col in self.format.columns:
            col_fields = [x.strip() for x in col.field.split(',') if x] if col.field.strip() else []
            if 'GT' in col_fields:
                for fld in col_fields:
                    if fld not in geno_fields:
                        geno_fields.append(fld)
            else:
                for fld in col_fields:
                    if fld not in var_fields:
                        var_fields.append(fld)
        #
        # get indexes of fields
        #
        field_indexes = {}
        for idx,fld in enumerate(var_fields):
            # var_info, field = idx
            field_indexes[(-1, fld.lower())] = idx
        #
        idx = len(var_fields)
        if geno_fields:
            for id in self.IDs:
                header = [x.lower() for x in self.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                for fld in geno_fields:
                    field_indexes[(id, fld.lower())] = idx
                    idx += 1
        # 
        # how to process each column
        sep = '\t' if self.format.delimiter is None else self.format.delimiter
        formatters = [] # formatters that will be used to produce strings from values
        col_adj = []        # adjust functions to combine values to one column.

        default_formatter = str if '*' not in self.format.formatter.keys() else self.getAdjFunc(self.format.formatter['*'])
        #
        col_idx = 0  # index of things after formatter.
        for col in self.format.columns:
            # indexes to get values for each column
            fields = [x.strip() for x in col.field.split(',') if x] if col.field else []
            if 'GT' in fields:
                for id in self.IDs:
                    col_indexes = []
                    indexes = [field_indexes[(id, x.lower())] for x in fields]
                    fmt = self.getFormatters(indexes, fields)
                    formatters.extend(fmt if fmt else [(None, None)])
                    col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else range(col_idx, col_idx +  len(fmt))])
                    col_idx += max(1, len(fmt))
                    if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                        raise ValueError('Columns with multiple fields must have an adjust function to merge values')
            else:
                indexes = [field_indexes[(-1, x.lower())] for x in fields]
                fmt = self.getFormatters(indexes, fields)
                formatters.extend(fmt if fmt else [(None, None)])
                col_adj.append([self.getAdjFunc(col.adj), col_idx if len(fmt) <= 1 else range(col_idx, col_idx + len(fmt))])
                col_idx += max(1, len(fmt))
                if col_adj[-1][0] is None and type(col_adj[-1][1]) is not int:
                    raise ValueError('Columns with multiple fields must have an adjust function to merge values')

        # needs fmt and adj
        count = 0
        failed_count = 0
        nr = self.db.numOfRows(self.table)
        last_count = 0
        update_after = min(max(100, nr/100), 10000)
        rec_stack = []
        # if export_by_fields is empty
        nFieldBy = len([x for x in self.format.export_by_fields.split(',') if x])
        #
        reader = VariantReader(self.proj, self.table, self.format.export_by_fields, self.format.order_by_fields,
            var_fields, geno_fields, self.export_alt_build, self.IDs, max(self.jobs - 1, 0))
        reader.start()
        prog = ProgressBar(self.filename if self.filename else 'Writing', nr)
        #
        # we cannot use a with statement here because we cannot close sys.stdout
        # perhaps a better solution is available.
        output = open(self.filename, 'w') if self.filename else sys.stdout
        # write header
        if self.header:
            # we can support more variables later
            header_keywords = ['sample_names']
            # interpolate %(VAR)s with values
            for m in re.finditer('%\((\w+)\)s', self.header):
                if m.group(1) not in header_keywords:
                    env.logger.warning('variable {} is not supported in --header.'.format(m.group(1)))
            #
            # sample_names
            header = self.header.replace('%(sample_names)s',
                        self.format.delimiter.join(self.samples))
            print >> output, header.rstrip()
        global rec_ref, rec_alt
        for idx, raw_rec in enumerate(reader.records()):
            multi_records = False
            try:
                if nFieldBy != 0:
                    if not rec_stack:
                        rec_stack.append(raw_rec)
                        continue
                    # if the same, wait for the next record
                    elif rec_stack[-1][2:(nFieldBy+2)] == raw_rec[2:(2+nFieldBy)]:
                        # the record needs to be merged even with the same export_by_fields (e.g. chr pos ref with two alts)
                        if rec_stack[-1] != raw_rec:
                            rec_stack.append(raw_rec)
                        continue
                    elif len(rec_stack) == 1:
                        rec_ref = rec_stack[0][0]
                        rec_alt = rec_stack[0][1]
                        rec = rec_stack[0][2:]
                        rec_stack = [raw_rec]
                    else:
                        n = len(rec_stack)
                        rec_ref = [rec_stack[i][0] for i in range(n)]
                        rec_alt = [rec_stack[i][1] for i in range(n)]
                        rec = [tuple([rec_stack[i][x] for i in range(n)]) for x in range(2, len(raw_rec))]
                        multi_records = True
                        rec_stack = [raw_rec]
                else:
                    rec_ref = raw_rec[0]
                    rec_alt = raw_rec[1]
                    rec = raw_rec[2:]
                # step one: apply formatters 
                # if there is no fmt, the item must be either empty or a single item
                #
                # fmt: single or list
                # no fmt: single or None
                #
                # this is extremely ugly but are we getting any performance gain?
                if multi_records:
                    try:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col][0] is None) else default_formatter(rec[col][0])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'.format(
                                    rec[col] if type(col) is int else [rec[x] for x in col], col, e))
                else:
                    try:
                        fields = [fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col] is None) else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'.format(
                                    rec[col] if type(col) is int else [rec[x] for x in col], col, e))
                # step two: apply adjusters
                #
                # adj: single or list
                # no adj: must be single
                columns = [adj(fields[col] if type(col) is int else [fields[x] for x in col]) if adj else fields[col] for adj, col in col_adj]
                # step three: output columns
                line = sep.join(columns)
                output.write(line + '\n')
                count += 1
            except Exception as e:
                env.logger.debug('Failed to process record {}: {}'.format(rec, e))
                failed_count += 1
            if idx - last_count > update_after:
                last_count = idx
                prog.update(idx)
        # the last block
        if rec_stack:
            try:
                n = len(rec_stack)
                multi_records = n > 1
                if n == 1:
                    rec_ref = rec_stack[0][0]
                    rec_alt = rec_stack[0][1]
                    rec = rec_stack[0][2:]
                else:
                    rec_ref = [rec_stack[i][0] for i in range(n)]
                    rec_alt = [rec_stack[i][1] for i in range(n)]
                    rec = [tuple([rec_stack[i][x] for i in range(n)]) for x in range(2, len(raw_rec))]
                # step one: apply formatters
                if multi_records:
                    try:
                        fields = [fmt(None if col is None else (rec[col] \
                            if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col][0] is None)\
                            else default_formatter(rec[col][0])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] \
                                    if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'
                                    .format(rec[col] if type(col) is int else \
                                        [rec[x] for x in col], col, e))
                else:
                    try:
                        fields = [fmt(None if col is None else (rec[col] \
                            if type(col) is int else [rec[x] for x in col])) \
                            if fmt else ('' if (col is None or rec[col] is None)\
                            else default_formatter(rec[col])) for fmt, col in formatters]
                    except:
                        for fmt, col in formatters:
                            try:
                                if fmt:
                                    fmt(None if col is None else (rec[col] \
                                    if type(col) is int else [rec[x] for x in col]))
                            except Exception as e:
                                raise ValueError('Failed to format value {} at col {}: {}'
                                    .format(rec[col] if type(col) is int else \
                                        [rec[x] for x in col], col, e))
                # step two: apply adjusters
                columns = [adj(fields[col] if type(col) is int else \
                    [fields[x] for x in col]) if adj else fields[col] \
                    for adj, col in col_adj]
                # step three: output columns
                line = sep.join(columns)
                output.write(line + '\n')
                count += 1
            except Exception as e:
                env.logger.debug('Failed to process record {}: {}'.format(rec, e))
                failed_count += 1
        if self.filename is not None:
            output.close()
        prog.done()
        env.logger.info('{} lines are exported from variant table {} {}'
            .format(count, self.table, '' if failed_count == 0 else \
                'with {} failed records'.format(failed_count)))



# Functions provided by this script
#
#

def exportArguments(parser):
    parser.add_argument('table', help='''A variant table whose variants will be exported.
        If parameter --samples is specified, only variants belong to one or more of the
        samples will be exported.''')
    parser.add_argument('filename', nargs='?', help='''(deprecated) Name of output file.
        Output will be written to the standard output if this parameter is left unspecified.
        This parameter is deprecated and might be removed in a future version of
        variant tools.''')
    parser.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Samples that will be exported, specified by conditions such as 'aff=1'
            and 'filename like "MG%%"'. Multiple samples could be exported to a
            file if the output format allows. No sample will be outputted if this
            parameter is ignored.''')
    parser.add_argument('--format',
        help='''Format of the exported file. It can be one of the variant tools
            supported file types such as VCF (cf. 'vtools show formats') or a local
            format specification file (with extension .fmt). Some formats accept
            additional parameters (cf. 'vtools show format FMT') and allows you to
            export additional or alternative fields.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the exported data. It
            can only be one of the primary (default) of alternative (if exists) reference
            genome of the project.''')
    parser.add_argument('--header', nargs='*', 
        help='''A complete header or a list of names that will be joined by a
            delimiter specified by the file format to form a header. If a special
            name - is specified, the header will be read from the standard input,
            which is the preferred way to specify large multi-line headers (e.g.
            cat myheader | vtools export --header -). Strings in the form of
            %%(VAR)s will be interpolated to values of variable VAR, which can be
            "sample_names" for list of sample names.''')
    parser.add_argument('-j', '--jobs', type=int, default=1,
        help='''Number of processes to export data. Multiple threads will be automatically
            used if there are a large number of samples.''')

def export(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            exporter = Exporter(proj=proj, table=args.table, filename=args.filename,
                samples=' AND '.join(['({})'.format(x) for x in args.samples]),
                format=args.format, build=args.build, header=args.header,  
                jobs=args.jobs, fmt_args=args.unknown_args)
            exporter.exportData()
        proj.close()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
