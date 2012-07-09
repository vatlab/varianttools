#!/usr/bin/env python
#
# $File: update.py $
# $LastChangedDate: 2012-06-08 17:25:06 -0500 (Fri, 08 Jun 2012) $
# $Rev: 1189 $
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
import os
import re
import ConfigParser
import shutil
import urlparse
import gzip
import zipfile
from collections import defaultdict
from multiprocessing import Process, Pipe
import math
from .project import AnnoDB, Project, Field, AnnoDBWriter, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, DatabaseEngine, delayedAction, \
    consolidateFieldName, runOptions
from .importer import *

#
#
#  Command update
#
#
class Updater:
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, table, files, build, format, jobs, fmt_args=[]):
        # if update is None, recreate index
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
        #
        for filename in files:
            if not os.path.isfile(filename):
                raise ValueError('File {} does not exist'.format(filename))
        self.files = files
        # for #record, #genotype (new or updated), #new variant, SNV, insertion, deletion, complex variants, invalid record, updated record
        self.count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.total_count = [0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.import_alt_build = False
        if len(self.files) == 0:
            raise ValueError('No file to import')
        #
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                self.logger.info('Using primary reference genome {} of the project.'.format(self.build))
            else:
                raise ValueError('Please specify a reference genome using parameter --build')
        else:
            self.build = build
        #
        if self.proj.build is None:
            raise ValueError('Cannot update variants of a project without variants.')
        elif self.build == self.proj.build:
            # perfect case
            pass
        elif self.build == self.proj.alt_build:
            # troublesome
            self.import_alt_build = True
        elif self.proj.alt_build is None:
            # even more troublesome
            self.logger.warning('The new files uses a different reference genome ({}) from the primary reference genome ({}) of the project.'.format(self.build, self.proj.build))
            self.logger.info('Adding an alternative reference genome ({}) to the project.'.format(self.build))
            tool = LiftOverTool(self.proj)
            # in case of insert, the indexes will be dropped soon so do not build
            # index now
            tool.setAltRefGenome(self.build, build_index=True)
            self.import_alt_build = True

        #
        self.jobs = max(1, jobs)
        if not proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(table))
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        # try to guess file type
        if not format:
            filename = self.files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            fmt = fileFMT(format, fmt_args, logger=self.logger)
        except Exception as e:
            self.logger.debug(e)
            raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '.format(e, format))
        #
        # how to split processed records
        self.ranges = fmt.ranges
        self.variant_fields = [x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]]
        self.variant_info = [x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]]
        self.genotype_field = [x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]]
        self.genotype_info = [x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]]
        #
        if not self.variant_info and not self.genotype_info:
            raise ValueError('No variant or genotype info needs to be updated')
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for pos, ref, alt are 1, 2, 3 in fields.
            self.processor = LineProcessor(fmt.fields, [(1, 2, 3)], fmt.delimiter, self.ranges, self.logger)
        else:  # position or range type
            self.processor = LineProcessor(fmt.fields, [(1,)], fmt.delimiter, self.ranges, self.logger)
        # probe number of sample
        if self.genotype_field and self.genotype_info:
            self.prober = LineProcessor([fmt.fields[fmt.ranges[2]]], [], fmt.delimiter, None, self.logger)
        # there are variant_info
        if self.variant_info:
            cur = self.db.cursor()
            headers = self.db.getHeaders('variant')
            for f in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]:
                # either insert or update, the fields must be in the master variant table
                self.proj.checkFieldName(f.name, exclude='variant')
                if f.name not in headers:
                    s = delayedAction(self.logger.info, 'Adding column {}'.format(f.name))
                    cur.execute('ALTER TABLE variant ADD {} {};'.format(f.name, f.type))
                    del s
        #if len(self.variant_info) == 0 and len(self.genotype_info == 0:
        #    raise ValueError('No field could be updated using this input file')
        #
        self.input_type = fmt.input_type
        self.encoding = fmt.encoding
        fbin, fchr, fpos = ('alt_bin', 'alt_chr', 'alt_pos') if self.import_alt_build else ('bin', 'chr', 'pos')
        from_table = 'AND variant.variant_id IN (SELECT variant_id FROM {})'.format(table) if table != 'variant' else ''
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH)
        self.update_position_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} = {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        self.update_range_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} >= {0} AND variant.{4} <= {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        #
        self.variantIndex = self.proj.createVariantMap(table, self.import_alt_build)
        self.table = table

    def updateVariant(self, cur, bins, rec):
        if self.input_type == 'variant':
            var_key = (rec[0], rec[2], rec[3])
            if var_key in self.variantIndex and rec[1] in self.variantIndex[var_key]:
                variant_id = self.variantIndex[var_key][rec[1]][0]
                # update by variant_id, do not need bins
                if len(rec) > 4:
                    cur.execute(self.update_variant_query, rec[4:] + [variant_id])
                    self.count[8] += cur.rowcount
                return variant_id
        elif self.input_type == 'position':
            cur.execute(self.update_position_query, rec[2:] + bins + [rec[0], rec[1]])
            self.count[8] += cur.rowcount
        else:  # range based
            cur.execute(self.update_range_query, rec[3:] + bins + [rec[0], rec[1], rec[2]])
            self.count[8] += cur.rowcount
        return None

    def getSampleIDs(self, filename):
        if not self.genotype_field:
            # no genotype_field, good, do not have to worry about genotype
            return []
        # has the file been imported before?
        cur = self.db.cursor()
        cur.execute('SELECT filename from filename;')
        existing_files = [x[0] for x in cur.fetchall()]
        if filename not in existing_files:
            return []
        #
        # what are the samples related to this file?
        cur.execute('SELECT sample_id, sample_name FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id WHERE filename.filename = {}'\
            .format(self.db.PH), (filename,))
        sample_ids = []
        sample_names = []
        for rec in cur:
            sample_ids.append(rec[0])
            sample_names.append(rec[1])
        # what is the sample names get from this file?
        nSample, names = probeSampleName(filename, self.prober, self.encoding, self.logger)
        if nSample != len(sample_ids):
            self.logger.warning('Number of samples mismatch. Cannot update genotype')
            return []
        if nSample == 1:
            # if only one sample, update it regardless of sample name.
            return sample_ids
        if sample_names == names:
            # if sample name matches, get sample_ids
            return sample_ids
        else:
            self.logger.warning('Sample names mismatch. Cannot update genotype.')
            return []
        
    def updateFromFile(self, input_filename):
        self.processor.reset()
        if self.genotype_field and self.genotype_info:
            self.prober.reset()
        #
        # do we handle genotype info as well?
        sample_ids = self.getSampleIDs(input_filename) if self.genotype_info else []
        # one process is for the main program, the
        # other thread will handle input
        reader = TextReader(self.processor, input_filename,
            # in the case of variant, we filter from the reading stage to save some time
            None if (self.table == 'variant' or self.input_type != 'variant') else self.variantIndex,
            # getNew is False so we only get variants that are available in variantIndex
            False, self.jobs - 1, self.encoding, self.logger)
        #
        # do we need to add extra columns to the genotype tables
        if sample_ids:
            s = delayedAction(self.logger.info, 'Preparing genotype tables (adding needed fields and indexes)...')
            cur = self.db.cursor()
            for id in sample_ids:
                headers = [x.upper() for x in self.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
                if 'GT' not in headers:  # for genotype
                    self.logger.debug('Adding column GT to table genotype_{}'.format(id))
                    cur.execute('ALTER TABLE {}_genotype.genotype_{} ADD {} {};'.format(self.proj.name, id, 'GT', 'INT'))
                for field in self.genotype_info:
                    if field.name.upper() not in headers:
                        self.logger.debug('Adding column {} to table genotype_{}'.format(field.name, id))
                        cur.execute('ALTER TABLE {}_genotype.genotype_{} ADD {} {};'.format(self.proj.name, id, field.name, field.type))
            # if we are updating by variant_id, we will need to create an index for it
            for id in sample_ids:
                if not self.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.proj.name, id)):
                    cur.execute('CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'.format(self.proj.name, id))
            del s
            genotype_update_query = {id: 'UPDATE {0}_genotype.genotype_{1} SET {2} WHERE variant_id = {3};'\
                .format(self.proj.name, id,
                ', '.join(['{}={}'.format(x, self.db.PH) for x in [y.name for y in self.genotype_info]]),
                self.db.PH)
                for id in sample_ids}
        else:
            # do not import genotype even if the input file has them
            self.genotype_field = []
            self.genotype_info = []
            self.processor.reset(import_var_info=True, import_sample_range=[self.ranges[2], self.ranges[2]])
        #
        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc//200, 100), 100000)
        fld_cols = None
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0
        for self.count[0], bins, rec in reader.records():
            variant_id = self.updateVariant(cur, bins, rec[0:self.ranges[2]])
            # variant might not exist
            if variant_id is not None and sample_ids:
                if fld_cols is None:
                    col_rngs = [reader.columnRange[x] for x in range(self.ranges[3], self.ranges[4])]
                    fld_cols = []
                    for idx in range(len(sample_ids)):
                        fld_cols.append([sc + (0 if sc + 1 == ec else idx) for sc,ec in col_rngs])
                    if col_rngs[0][1] - col_rngs[0][0] != len(sample_ids):
                        self.logger.error('Number of genotypes ({}) does not match number of samples ({})'.format(
                            col_rngs[0][1] - col_rngs[0][0], len(sample_ids)))
                for idx, id in enumerate(sample_ids):
                    if rec[self.ranges[2] + idx] is not None:
                        cur.execute(genotype_update_query[id], [rec[c] for c in fld_cols[idx]] + [variant_id])
                        self.count[1] += 1
            if self.count[0] - last_count > update_after:
                last_count = self.count[0]
                self.db.commit()
                prog.update(self.count[0])
        self.count[7] = reader.skipped_lines
        self.db.commit()
        prog.done(self.count[0])

    def update(self):
        '''Start updating'''
        for count,f in enumerate(self.files):
            self.logger.info('{} variants from {} ({}/{})'.format('Updating', f, count + 1, len(self.files)))
            self.updateFromFile(f)
            self.logger.info('Field{} {} of {:,} variants{} are updated'.format('' if len(self.variant_info) == 1 else 's', ', '.join(self.variant_info), self.count[8],
                    '' if self.count[1] == 0 else ' and geno fields of {:,} genotypes'.format(self.count[1])))
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0
        if len(self.files) > 1:
            self.logger.info('Field{} {} of {:,} variants{} are updated'.format('' if len(self.variant_info) == 1 else 's', ', '.join(self.variant_info), self.total_count[8],
                    '' if self.total_count[1] == 0 else ' and geno fields of {:,} genotypes'.format(self.total_count[1])))


def setFieldValue(proj, table, items, build):
    # fields
    expr = ','.join([x.split('=',1)[1] for x in items])
    select_clause, fields = consolidateFieldName(proj, table, expr,
        build and build == proj.alt_build)
    #
    # FROM clause
    from_clause = 'FROM {} '.format(table)
    fields_info = sum([proj.linkFieldToTable(x, table) for x in fields], [])
    #
    processed = set()
    for tbl, conn in [(x.table, x.link) for x in fields_info if x.table != '']:
        if (tbl.lower(), conn.lower()) not in processed:
            from_clause += ' LEFT OUTER JOIN {} ON {}'.format(tbl, conn)
            processed.add((tbl.lower(), conn.lower()))
    # running query
    cur = proj.db.cursor()
    #
    if 'LEFT OUTER JOIN' not in from_clause:  # if everything is done in one table
        query = 'SELECT {} {};'.format(select_clause, from_clause)
        proj.logger.debug('Running {}'.format(query))
        cur.execute(query)
        fldTypes = [None] * len(items)
        for rec in cur:
            for i in range(len(items)):
                if rec[i] is None: # missing
                    continue
                elif fldTypes[i] is None:
                    fldTypes[i] = type(rec[i])
                    continue
                if type(rec[i]) != fldTypes[i]:
                    if type(rec[i]) is float and fldTypes[i] is int:
                        fltType[i] = float
                    else:
                        raise ValueError('Inconsistent type returned from different samples')
            if all([x is not None for x in fldTypes]):
                break
        #
        count = [0]*3
        # if adding a new field
        cur_fields = proj.db.getHeaders(table)[3:]
        for field, fldType in zip([x.split('=', 1)[0] for x in items], fldTypes):
            if field.lower() not in [x.lower() for x in cur_fields]:
                if fldType is None:
                    raise ValueError('Cannot determine the value of a new field {} because the values are all NULL'.format(field))
                proj.checkFieldName(field, exclude=table)
                proj.logger.info('Adding field {}'.format(field))
                query = 'ALTER TABLE {} ADD {} {} NULL;'.format(table, field,
                    {int: 'INT',
                     float: 'FLOAT',
                     str: 'VARCHAR(255)',
                     unicode: 'VARCHAR(255)'}[fldType])
                proj.logger.debug(query)
                cur.execute(query)
                count[1] += 1  # new
            else:
                # FIXME: check the case for type mismatch
                proj.logger.info('Updating field {}'.format(field))
                count[2] += 1  # updated
        proj.db.commit()
        # really update things
        query = 'UPDATE {} SET {};'.format(table, ','.join(items))
        proj.logger.debug('Running {}'.format(query))
        cur.execute(query)
    else:
        query = 'SELECT {}, {}.variant_id {};'.format(select_clause, table, from_clause)
        proj.logger.debug('Running {}'.format(query))
        cur.execute(query)
        fldTypes = [None] * len(items)
        s = delayedAction(proj.logger.info, 'Evaluating all expressions')
        results = cur.fetchall()
        del s
        #
        for res in results:
            for i in range(len(items)):
                if res[i] is None: # missing
                    continue
                elif fldTypes[i] is None:
                    fldTypes[i] = type(res[i])
                    continue
                if type(res[i]) != fldTypes[i]:
                    if type(res[i]) is float and fldTypes[i] is int:
                        fltType[i] = float
                    else:
                        raise ValueError('Inconsistent type returned from different samples')
            if all([x is not None for x in fldTypes]):
                break
        #
        count = [0]*3
        # if adding a new field
        cur_fields = proj.db.getHeaders(table)[3:]
        new_fields = [x.split('=', 1)[0] for x in items]
        for field, fldType in zip(new_fields, fldTypes):
            if field.lower() not in [x.lower() for x in cur_fields]:
                if fldType is None:
                    raise ValueError('Cannot determine the value of a new field {} because the values are all NULL'.format(field))
                proj.checkFieldName(field, exclude=table)
                proj.logger.info('Adding field {}'.format(field))
                query = 'ALTER TABLE {} ADD {} {} NULL;'.format(table, field,
                    {int: 'INT',
                     float: 'FLOAT',
                     str: 'VARCHAR(255)',
                     unicode: 'VARCHAR(255)'}[fldType])
                proj.logger.debug(query)
                cur.execute(query)
                count[1] += 1  # new
            else:
                # FIXME: check the case for type mismatch
                proj.logger.info('Updating field {}'.format(field))
                count[2] += 1  # updated
        proj.db.commit()
        # really update things
        query = 'UPDATE {} SET {} WHERE variant_id={};'.format(table,
            ','.join(['{}={}'.format(x, proj.db.PH) for x in new_fields]), proj.db.PH)
        proj.logger.debug('Using query {}'.format(query))
        prog = ProgressBar('Updating {}'.format(table), len(results))
        for count, res in enumerate(results):
            cur.execute(query, res)
            if count % 10000 == 0:
                prog.update(count)
        prog.done()


def calcSampleStat(proj, from_stat, IDs, variant_table, genotypes):
    '''Count sample allele count etc for specified sample and variant table'''
    if not proj.isVariantTable(variant_table):
        raise ValueError('"Variant table {} does not exist.'.format(variant_table))
    #
    #
    # NOTE: this function could be implemented using one or more query more
    # or less in the form of
    #
    # UPDATE variant SET something = something 
    # FROM 
    # (SELECT variant_id, avg(FIELD) FROM (
    #       SELECT FIELD FROM genotype_1 WHERE ...
    #       UNION SELECT FIELD FROM genotype_2 WHERE ...
    #       ...
    #       UNION SELECT FIELD FROM genotype_2 WHERE ) as total
    #   GROUP BY variant_id;
    #
    # This query can be faster because it is executed at a lower level, we cannot
    # really see the progress of the query though.
    #
    if not from_stat:
        proj.logger.warning('No statistics is specified')
        return

    # separate special functions...
    num = hom = het = other = cnt = None

    # keys to speed up some operations
    MEAN = 0
    SUM = 1
    MIN = 2
    MAX = 3
    operationKeys = {'avg': MEAN, 'sum': SUM, 'min': MIN, 'max': MAX}
    possibleOperations = operationKeys.keys()
    
    operations = []
    genotypeFields = []
    validGenotypeFields = []
    destinations = []
    fieldCalcs = []
    for stat in from_stat:
        f, e = [x.strip() for x in stat.split('=')]
        if e == '#(alt)':
            num = f
        elif e == '#(hom)':
            hom = f
        elif e == '#(het)':
            het = f
        elif e == '#(other)':
            other = f
        elif e == '#(GT)':
            cnt = f
        elif e.startswith('#('):
            raise ValueError('Unrecognized parameter {}: only parameters alt, hom, het, other and GT are accepted for special function #'.format(stat))
        else:
            m = re.match('(\w+)\s*=\s*(avg|sum|max|min)\s*\(\s*(\w+)\s*\)\s*', stat)
            if m is None:
                raise ValueError('Unrecognized parameter {}, which should have the form of FIELD=FUNC(GENO_INFO) where FUNC is one of #, avg, sum, max and min'.format(stat))
            dest, operation, field = m.groups()
            if operation not in possibleOperations:
                raise ValueError('Unsupported operation {}. Supported operations include {}.'.format(operation, ', '.join(possibleOperations)))
            operations.append(operationKeys[operation])
            genotypeFields.append(field)
            fieldCalcs.append(None)
            destinations.append(dest)
    #
    coreDestinations = [num, hom, het, other, cnt]
    cur = proj.db.cursor()
    if IDs is None:
        cur.execute('SELECT sample_id from sample;')
        IDs = [x[0] for x in cur.fetchall()]
    #
    numSample = len(IDs)
    if numSample == 0:
        proj.logger.warning('No sample is selected.')
        return
    
    # Error checking with the user specified genotype fields
    # 1) if a field does not exist within one of the sample genotype tables a warning is issued
    # 2) if a field does not exist in any sample, it is not included in validGenotypeFields
    # 3) if no fields are valid and no core stats were requested (i.e., num, het, hom, other), then sample_stat is exited
    genotypeFieldTypes = {}
    fieldInTable = defaultdict(list)
    for id in IDs:
        fields = [x.lower() for x in proj.db.getHeaders('{}_genotype.genotype_{}'.format(proj.name, id))]
        for field in fields:
            fieldInTable[field].append(id)
            if field not in genotypeFieldTypes:
                genotypeFieldTypes[field] = 'INT'
                fieldType = proj.db.typeOfColumn('{}_genotype.genotype_{}'.format(proj.name, id), field) 
                if fieldType.upper().startswith('FLOAT'):
                    genotypeFieldTypes[field] = 'FLOAT'
                elif fieldType.upper().startswith('VARCHAR'):
                    genotypeFieldTypes[field] = 'VARCHAR'
                    # We had been throwing an error here if a genotype field is a VARCHAR, but I think we should allow
                    # VARCHAR fields in the genotype tables.  We'll throw an error if someone wants to perform numeric operations on these fields
                    # further down in the code.
                    # raise ValueError('Genotype field {} is a VARCHAR which is not supported with sample_stat operations.'.format(field))

    validGenotypeIndices = []
    for index, field in enumerate(genotypeFields):
        if field.lower() not in [x.lower() for x in genotypeFieldTypes.keys()]:
            proj.logger.warning("Field {} is not an existing genotype field within your samples: {}".format(field, str(genotypeFieldTypes.keys())))
        else:
            if len(fieldInTable[field.lower()]) < len(IDs):
                proj.logger.warning('Field {} exists in {} of {} selected samples'.format(field, len(fieldInTable[field.lower()]), len(IDs))) 
            validGenotypeIndices.append(index)
            validGenotypeFields.append(field)
    # check GT field
    if not all([x is None for x in coreDestinations]):
        if 'gt' not in [x.lower() for x in genotypeFieldTypes.keys()]:
            proj.logger.warning('Genotype field does not exist in any of the selected samples')
        else:
            if len(fieldInTable['gt']) < len(IDs):
                proj.logger.warning('Genotype field GT exists in {} of {} selected samples'.format(len(fieldInTable[field.lower()]), len(IDs))) 

    if all([x is None for x in coreDestinations]) and len(validGenotypeFields) == 0:
        proj.logger.warning("No valid sample statistics operation has been specified.")
        return

    queryDestinations = coreDestinations
    for index in validGenotypeIndices:
        queryDestinations.append(destinations[index])
    for name in queryDestinations:
        if name is not None:
            proj.checkFieldName(name, exclude=variant_table)
    #
    variants = dict()
    prog = ProgressBar('Counting variants', len(IDs))
    prog_step = max(len(IDs) // 100, 1) 
    for id_idx, id in enumerate(IDs):
        where_cond = []
        if genotypes is not None and len(genotypes) != 0:
            where_cond.extend(genotypes)
        if variant_table != 'variant':
            where_cond.append('variant_id in (SELECT variant_id FROM {})'.format(variant_table))
        whereClause = 'WHERE ' + ' AND '.join(['({})'.format(x) for x in where_cond]) if where_cond else ''
        fieldSelect = ['GT' if ('gt' in fieldInTable and id in fieldInTable['gt']) else 'NULL']
        if validGenotypeFields is not None and len(validGenotypeFields) != 0:
            fieldSelect.extend([x if id in fieldInTable[x.lower()] else 'NULL' for x in validGenotypeFields])

        if not fieldSelect or all([x == 'NULL' for x in fieldSelect]):
            continue

        query = 'SELECT variant_id {} FROM {}_genotype.genotype_{} {};'.format(' '.join([',' + x for x in fieldSelect]),
                proj.name, id, whereClause)
        #proj.logger.debug(query)
        cur.execute(query)

        for rec in cur:
            if rec[0] not in variants:
                variants[rec[0]] = [0, 0, 0, 0]
                variants[rec[0]].extend(list(fieldCalcs))
            # total valid GT
            if rec[1] is not None:
                variants[rec[0]][3] += 1
            # type heterozygote
            if rec[1] == 1:
                variants[rec[0]][0] += 1
            # type homozygote
            elif rec[1] == 2:
                variants[rec[0]][1] += 1
            # type double heterozygote with two different alternative alleles
            elif rec[1] == -1:
                variants[rec[0]][2] += 1
            elif rec[1] not in [0, None]:
                proj.logger.warning('Invalid genotype type {}'.format(rec[1]))
            #
            # this collects genotype_field information
            if len(validGenotypeFields) > 0:
                for index in validGenotypeIndices:
                    queryIndex = index + 2     # to move beyond the variant_id and GT fields in the select statement
                    recIndex = index + 4       # first 4 attributes of variants are het, hom, double_het and GT
                    # ignore missing (NULL) values, and empty string that, if so inserted, could be returned
                    # by sqlite even when the field type is INT.
                    if rec[queryIndex] in [None, '']:
                        continue
                    operation = operations[index]
                    field = genotypeFields[index]
                    if operation == MEAN:
                        if variants[rec[0]][recIndex] is None:
                            # we need to track the number of valid records
                            variants[rec[0]][recIndex] = [rec[queryIndex], 1]
                        else:
                            variants[rec[0]][recIndex][0] += rec[queryIndex]
                            variants[rec[0]][recIndex][1] += 1
                    elif operation == SUM:
                        if variants[rec[0]][recIndex] is None:
                            variants[rec[0]][recIndex] = rec[queryIndex]
                        else:
                            variants[rec[0]][recIndex] += rec[queryIndex]
                    elif operation == MIN:
                        if variants[rec[0]][recIndex] is None or rec[queryIndex] < variants[rec[0]][recIndex]:
                            variants[rec[0]][recIndex] = rec[queryIndex]
                    elif operation == MAX:
                        if variants[rec[0]][recIndex] is None or rec[queryIndex] > variants[rec[0]][recIndex]:
                            variants[rec[0]][recIndex] = rec[queryIndex]  
        if id_idx % prog_step == 0:
            prog.update(id_idx + 1)
    prog.done()
    if len(variants) == 0:
        raise ValueError('No variant is updated')
    #
    headers = [x.lower() for x in proj.db.getHeaders(variant_table)]
    table_attributes = [(num, 'INT'), (hom, 'INT'),
            (het, 'INT'), (other, 'INT'), (cnt, 'INT')]
    fieldsDefaultZero = [num, hom, het, other, cnt]
    
    for index in validGenotypeIndices:
        field = genotypeFields[index]
        genotypeFieldType = genotypeFieldTypes.get(genotypeFields[index]) 
        
        if genotypeFieldType == 'VARCHAR':
            raise ValueError('Genotype field {} is a VARCHAR which is not supported with sample_stat operations.'.format(field))
        
        if operations[index] == MEAN:
            table_attributes.append((destinations[index], 'FLOAT'))
        else:                
            table_attributes.append((destinations[index], genotypeFieldType))
    for field, fldtype in table_attributes:
        if field is None:
            continue
        # We are setting default values on the count fields to 0. The genotype stat fields are set to NULL by default.
        defaultValue = 0 if field in fieldsDefaultZero else None
        if field.lower() in headers:
            if proj.db.typeOfColumn(variant_table, field) != (fldtype + ' NULL'):
                proj.logger.warning('Type mismatch (existing: {}, new: {}) for column {}. Please remove this column and recalculate statistics if needed.'\
                    .format(proj.db.typeOfColumn(variant_table, field), fldtype, field))
            proj.logger.info('Resetting values at existing field {}'.format(field))
            proj.db.execute('Update {} SET {} = {};'.format(variant_table, field, proj.db.PH), (defaultValue, ))
        else:
            proj.logger.info('Adding field {}'.format(field))
            proj.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(variant_table, field, fldtype))
            if defaultValue == 0:
                proj.db.execute ('UPDATE {} SET {} = 0'.format(variant_table, field))              
    #
    prog = ProgressBar('Updating {}'.format(variant_table), len(variants))
    update_query = 'UPDATE {0} SET {2} WHERE variant_id={1};'.format(variant_table, proj.db.PH,
        ' ,'.join(['{}={}'.format(x, proj.db.PH) for x in queryDestinations if x is not None]))
    warning = False
    for count,id in enumerate(variants):
        value = variants[id]
        res = []
        if num is not None:
            # het + hom * 2 + other
            res.append(value[0] + value[1] * 2 + value[2])
        if hom is not None:
            res.append(value[1])
        if het is not None:
            res.append(value[0])
        if other is not None:
            res.append(value[2])
        if cnt is not None:
            res.append(value[3])
        # for genotype_field operations, the value[operation_index] holds the result of the operation
        # except for the "mean" operation which needs to be divided by num_samples that have that variant
        try:
            for index in validGenotypeIndices:
                operationIndex = index + 4     # the first 4 indices hold the values for hom, het, double het and total genotype
                operationCalculation = value[operationIndex]
                if operations[index] == MEAN and operationCalculation is not None:
                    res.append(float(operationCalculation[0]) / operationCalculation[1])
                else:
                    res.append(operationCalculation)
            cur.execute(update_query, res + [id])
        except Exception as e:
            proj.logger.debug(e)
        if count % proj.db.batch == 0:
            proj.db.commit()
            prog.update(count)
    proj.db.commit()
    prog.done()
    proj.logger.info('{} records are updated'.format(count))


def updateArguments(parser):
    parser.add_argument('table', help='''variants to be updated.''')
    files = parser.add_argument_group('Update from files')
    files.add_argument('--from_file', nargs='+',
        help='''A list of files that will be used to add or update existing fields of
            variants. The file should be delimiter separated with format described by
            parameter --format. Gzipped files are acceptable. If input files contains
            genotype information, have been inputted before, and can be linked to the
            samples they created without ambiguity (e.g. single sample, or samples with
            detectable sample names), genotypes and their information will also be
            updated.''')
    files.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18) of the input files,
            which should be the primary (used by default) or alternative (if available)
            reference genome of the project. An alternative reference genome will be
            added to the project if needed.''')
    files.add_argument('--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as ANNOVAR_output (cf. 'vtools show formats'),
            or a local format specification file (with extension .fmt). Some formats 
            accept parameters (cf. 'vtools show format FMT') and allow you to update
            additional or alternative fields from the input file.''')
    files.add_argument('-j', '--jobs', metavar='N', default=1, type=int,
        help='''Number of processes to import input file. Variant tools by default
            uses a single process for reading and writing, and can use one or more
            dedicated reader processes (jobs=2 or more) to process input files. Due
            to the overhead of inter-process communication, more jobs do not
            automatically lead to better performance.''')
    field = parser.add_argument_group('Set value from existing fields')
    field.add_argument('--set', metavar='EXPR', nargs='*', default=[],
        help='''Add a new field or updating an existing field using a constant
            (e.g. mark=1) or an expression using other fields (e.g. freq=num/120,
            refgene=refGene.name). If multiple values are returned for a variant, only
            one of them will be used. Parameter samples could be used to limit the
            affected variants. In addition, special function are provided, including
            'HWE_exact' (exact test of Hardy-Weinberg Equilibrium) and
            'Fisher_exact' (Fisher's exact test for case/ctrl association).''')
    #field.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
    #    help='''Limiting variants from samples that match conditions that
    #        use columns shown in command 'vtools show sample' (e.g. 'aff=1',
    #        'filename like "MG%%"').''')
    #field.add_argument('--build',
    #    help='''Reference genome to use when chr and pos is used in expression.''')
    stat = parser.add_argument_group('Set fields from sample statistics')
    stat.add_argument('--from_stat', metavar='EXPR', nargs='*', default=[],
        help='''One or more expressions such as meanQT=avg(QT) that aggregate genotype info (e.g. QT)
            of variants in all or selected samples to specified fields (e.g. meanQT). Functions sum, avg,
            max, and min are currently supported. In addition, special functions #(GT), #(alt), #(hom),
            #(het) and #(other) are provided to count the number of valid genotypes (not NULL),
            alternative alleles, homozygotes, heterozygotes, and individuals with two different
            alternative alleles.''')
    stat.add_argument('-s', '--samples', nargs='*', metavar='COND', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    stat.add_argument('--genotypes', nargs='*', metavar='COND', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show genotypes' (e.g. 'GQ>15').''')


def update(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if not args.from_file and not args.from_stat and not args.set:
                proj.logger.warning('Please specify one of --from_file, --set and --from_stat for command vtools upate')
            if args.from_file:
                proj.db.attach(proj.name + '_genotype')
                updater = Updater(proj=proj, table=args.table, files=args.from_file,
                    build=args.build, format=args.format, jobs=args.jobs, fmt_args=args.unknown_args)
                updater.update()
            if args.set:
                for item in args.set:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError('Invalid parameter {}, which should have format field=expr_of_field: {}'.format(item, e))
                setFieldValue(proj, args.table, args.set, proj.build)
                #, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.from_stat:
                proj.db.attach(proj.name + '_genotype')
                variant_table = args.table if args.table else 'variant'
                if not proj.db.hasTable(variant_table):
                    raise ValueError('Variant table {} does not exist'.format(variant_table))
                IDs = None
                if args.samples:
                    IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.samples]))
                    if len(IDs) == 0:
                        proj.logger.info('No sample is selected (or available)')
                        return
                    else:
                        proj.logger.info('{} samples are selected'.format(len(IDs)))
                calcSampleStat(proj, args.from_stat, IDs, variant_table, args.genotypes)
        proj.close()
    except Exception as e:
        sys.exit(e)

