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
import os
import re
import sys
import time
from collections import defaultdict

from .geno_store import GenoStore
from .importer import LineProcessor, probeSampleName
from .importer_allele_hdf5 import UpdateGenotypeInParallel
from .liftOver import LiftOverTool
from .project import Project, fileFMT
from .text_reader import TextReader
from .utils import (ProgressBar, RefGenome, consolidateFieldName,
                    decodeTableName, delayedAction, determineSexOfSamples,
                    encodeTableName, env, getVariantsOnChromosomeX,
                    getVariantsOnChromosomeY, getVariantsOnManifolds,
                    lineCount)


#
#
#  Command update
#
#
class Updater:
    '''Import variants from one or more tab or comma separated files.'''

    def __init__(self,
                 proj,
                 table,
                 files,
                 build,
                 format,
                 jobs,
                 sample_name,
                 fmt_args=[]):
        # if update is None, recreate index
        self.proj = proj
        self.db = proj.db
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
        self.sample_name = sample_name
        if build is None:
            if self.proj.build is not None:
                self.build = self.proj.build
                env.logger.info(
                    'Using primary reference genome {} of the project.'.format(
                        self.build))
            else:
                raise ValueError(
                    'Please specify a reference genome using parameter --build')
        else:
            self.build = build
        #
        if self.proj.build is None:
            raise ValueError(
                'Cannot update variants of a project without variants.')
        elif self.build == self.proj.build:
            # perfect case
            pass
        elif self.build == self.proj.alt_build:
            # troublesome
            self.import_alt_build = True
        elif self.proj.alt_build is None:
            # even more troublesome
            env.logger.warning(
                'The new files uses a different reference genome ({}) from the primary reference genome ({}) of the project.'
                .format(self.build, self.proj.build))
            env.logger.info(
                'Adding an alternative reference genome ({}) to the project.'
                .format(self.build))
            tool = LiftOverTool(self.proj)
            # in case of insert, the indexes will be dropped soon so do not build
            # index now
            tool.setAltRefGenome(self.build, build_index=True)
            self.import_alt_build = True

        #
        self.jobs = max(1, jobs)
        if not proj.isVariantTable(table):
            raise ValueError('Variant table {} does not exist.'.format(
                decodeTableName(table)))
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError(
                'Please specify the reference genome of the input data.')
        #
        # try to guess file type
        if not format:
            filename = self.files[0].lower()
            if filename.endswith('.vcf') or filename.endswith('.vcf.gz'):
                format = 'vcf'
            else:
                raise ValueError('Cannot guess input file type from filename')
        try:
            fmt = fileFMT(format, fmt_args)
        except Exception as e:
            env.logger.debug(e)
            raise IndexError(
                'Unrecognized input format: {}\nPlease check your input parameters or configuration file *{}* '
                .format(e, format))
        #
        # how to split processed records
        self.ranges = fmt.ranges
        self.variant_fields = [
            x.name for x in fmt.fields[fmt.ranges[0]:fmt.ranges[1]]
        ]
        self.variant_info = [
            x.name for x in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]
        ]
        self.genotype_field = [
            x.name for x in fmt.fields[fmt.ranges[2]:fmt.ranges[3]]
        ]
        self.genotype_info = [
            x for x in fmt.fields[fmt.ranges[3]:fmt.ranges[4]]
        ]
        #
        if not self.variant_info and not self.genotype_info:
            raise ValueError('No variant or genotype info needs to be updated')
        #
        if fmt.input_type == 'variant':
            # process variants, the fields for chr, pos, ref, alt are 0, 1, 2, 3 in fields.
            self.processor = LineProcessor(
                fmt.fields, [(RefGenome(self.build), 0, 1, 2, 3)],
                fmt.delimiter, self.ranges)
        else:  # position or range type
            self.processor = LineProcessor(fmt.fields, [(1,)], fmt.delimiter,
                                           self.ranges)

        # probe number of sample
        if self.genotype_field or self.genotype_info:
            self.prober = LineProcessor([fmt.fields[fmt.ranges[2]]], [],
                                        fmt.delimiter, None)
        # there are variant_info
        if self.variant_info:
            cur = self.db.cursor()
            headers = self.db.getHeaders('variant')
            for f in fmt.fields[fmt.ranges[1]:fmt.ranges[2]]:
                # either insert or update, the fields must be in the master variant table
                self.proj.checkFieldName(f.name, exclude='variant')
                if f.name.upper() not in [x.upper() for x in headers]:
                    with delayedAction(env.logger.info,
                                       'Adding column {}'.format(f.name)):
                        cur.execute('ALTER TABLE variant ADD {} {};'.format(
                            f.name, f.type))
        #if len(self.variant_info) == 0 and len(self.genotype_info == 0:
        #    raise ValueError('No field could be updated using this input file')
        #
        self.input_type = fmt.input_type
        self.encoding = fmt.encoding
        self.header = fmt.header
        fbin, fchr, fpos = ('alt_bin', 'alt_chr',
                            'alt_pos') if self.import_alt_build else ('bin',
                                                                      'chr',
                                                                      'pos')
        from_table = 'AND variant.variant_id IN (SELECT variant_id FROM {})'.format(
            table) if table != 'variant' else ''
        self.update_variant_query = 'UPDATE variant SET {0} WHERE variant.variant_id = {1};'\
            .format(', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), self.db.PH)
        self.update_position_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} = {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        self.update_range_query = 'UPDATE variant SET {1} WHERE variant.{2} = {0} AND variant.{3} = {0} AND variant.{4} >= {0} AND variant.{4} <= {0} {5};'\
            .format(self.db.PH, ', '.join(['{}={}'.format(x, self.db.PH) for x in self.variant_info]), fbin, fchr, fpos, from_table)
        #
        self.variantIndex = self.proj.createVariantMap(table,
                                                       self.import_alt_build)
        self.table = table

    def updateVariant(self, cur, bins, rec):
        if self.input_type == 'variant':
            var_key = (rec[0], rec[2], rec[3])
            if var_key in self.variantIndex and rec[1] in self.variantIndex[
                    var_key]:
                variant_id = self.variantIndex[var_key][rec[1]][0]
                # update by variant_id, do not need bins
                if len(rec) > 4:
                    cur.execute(self.update_variant_query,
                                rec[4:] + [variant_id])
                    self.count[8] += cur.rowcount
                return variant_id
        elif self.input_type == 'position':
            cur.execute(self.update_position_query,
                        rec[2:] + bins + [rec[0], rec[1]])
            self.count[8] += cur.rowcount
        else:  # range based
            cur.execute(self.update_range_query,
                        rec[3:] + bins + [rec[0], rec[1], rec[2]])
            self.count[8] += cur.rowcount
        return None

    def getSampleIDs(self, input_filename):
        '''Return a list of sample ids, and a status field indicating
        0:  no sample. sample id has to be empty
        1:  no genotype but with sample (a list of variant ids, no genotype), there can be only one sample id.
        2:  genotype with samples.
        '''
        if not self.sample_name:
            try:
                numSample, names = probeSampleName(input_filename, self.prober,
                                                   self.encoding)
                self.sample_name = names
                if not names:
                    if numSample == 1:
                        env.logger.debug(
                            'Missing sample name (name None is used)')
                        self.sample_name = ['None']
                    elif numSample == 0:
                        return ([], 0)
                    else:
                        env.logger.warning(
                            'Failed to names for {} samples from input file'
                            .format(numSample))
            except ValueError as e:
                # cannot find any genotype column, perhaps no genotype is defined in the file (which is allowed)
                env.logger.warning(
                    'Failed to get sample names from input file: {}'.format(e))
                return ([], 0)
        else:
            # if sample names are provided, we try to check
            #
            # 1. if sample name is provided from input file, they should match
            # 2. if there is an sample that can be uniqly identified by
            #     sample_name and filename.
            try:
                numSample, names = probeSampleName(input_filename, self.prober,
                                                   self.encoding)
            except ValueError as e:
                env.logger.debug(e)
                numSample = 0
                #
            if numSample == 0:
                env.logger.warning(
                    'No genotype column could be found from the input file. Assuming no genotype.'
                )
                # remove genotype field from processor
                self.processor.reset(
                    import_var_info=True, import_sample_range=[0, 0])
                if len(self.sample_name) > 1:
                    raise ValueError(
                        "When there is no sample genotype, only one sample name is allowed."
                    )
            elif len(self.sample_name) != numSample:
                raise ValueError(
                    '{} sample detected but only {} sample names are specified'
                    .format(numSample, len(self.sample_name)))
            elif sorted(self.sample_name) != sorted(names):
                raise ValueError(
                    'Specified ({}) and detected sample names ({}) mismatch.'
                    .format(', '.join(self.sample_name), ', '.join(names)))
        #
        # now get sample ids:
        cur = self.db.cursor()
        sample_ids = []
        for name in self.sample_name:
            # get filename and sample_id
            cur.execute('SELECT sample_id, filename FROM sample '
                'LEFT OUTER JOIN filename ON sample.file_id = filename.file_id '
                'WHERE sample.sample_name = {}'\
                .format(self.db.PH), (name,))
            res = cur.fetchall()
            if len(res) == 0:
                env.logger.warning(
                    'Name {} does not match any existing sample'.format(name))
                sample_ids.append(-1)
            elif len(res) == 1:
                sample_ids.append(res[0][0])
            else:
                # if we can find one matching using filename
                res = [x for x in res if x[1] == input_filename]
                if len(res) == 1:
                    sample_ids.append(res[0][0])
                else:
                    raise ValueError(
                        'Cannot identify an unique sample using '
                        'sample name ({}) and filename ({})'.format(
                            name, input_filename))
        if all([x == -1 for x in sample_ids]):
            sample_ids = []
        return sample_ids

    def updateFromFile(self, input_filename):
        self.processor.reset()
        if self.genotype_field or self.genotype_info:
            self.prober.reset()

        #
        # do we handle genotype info as well?
        sample_ids = self.getSampleIDs(
            input_filename) if self.genotype_info else []
        # one process is for the main program, the
        # other thread will handle input
        reader = TextReader(
            self.processor,
            input_filename,
            # in the case of variant, we filter from the reading stage to save some time
            None if (self.table == 'variant' or self.input_type != 'variant')
            else self.variantIndex,
            # getNew is False so we only get variants that are available in variantIndex
            False,
            self.jobs - 1,
            self.encoding,
            self.header)
        #
        # do we need to add extra columns to the genotype tables
        if sample_ids:
            with delayedAction(
                    env.logger.info,
                    'Preparing genotype tables (adding needed fields and indexes)...'
            ):
                cur = self.db.cursor()
                for id in sample_ids:
                    if id == -1:
                        continue
                    headers = [
                        x.upper()
                        for x in self.db.getHeaders('{}_genotype.genotype_{}'
                                                    .format(self.proj.name, id))
                    ]
                    if 'GT' not in headers:  # for genotype
                        env.logger.trace(
                            'Adding column GT to table genotype_{}'.format(id))
                        cur.execute(
                            'ALTER TABLE {}_genotype.genotype_{} ADD {} {};'
                            .format(self.proj.name, id, 'GT', 'INT'))
                    for field in self.genotype_info:
                        if field.name.upper() not in headers:
                            env.logger.trace(
                                'Adding column {} to table genotype_{}'.format(
                                    field.name, id))
                            cur.execute(
                                'ALTER TABLE {}_genotype.genotype_{} ADD {} {};'
                                .format(self.proj.name, id, field.name,
                                        field.type))
                # if we are updating by variant_id, we will need to create an index for it
                for id in sample_ids:
                    if id == -1:
                        continue
                    if not self.db.hasIndex(
                            '{0}_genotype.genotype_{1}_index'.format(
                                self.proj.name, id)):
                        cur.execute(
                            'CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'
                            .format(self.proj.name, id))
            genotype_update_query = {id: 'UPDATE {0}_genotype.genotype_{1} SET {2} WHERE variant_id = {3};'\
                .format(self.proj.name, id,
                ', '.join(['{}={}'.format(x, self.db.PH) for x in [y.name for y in self.genotype_info]]),
                self.db.PH)
                for id in sample_ids if id != -1}
        else:
            # do not import genotype even if the input file has them
            self.genotype_field = []
            self.genotype_info = []
            self.processor.reset(
                import_var_info=True,
                import_sample_range=[self.ranges[2], self.ranges[2]])
        #

        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc // 200, 100), 100000)
        fld_cols = None
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0
        for self.count[0], bins, rec in reader.records():
            variant_id = self.updateVariant(cur, bins, rec[0:self.ranges[2]])
            # variant might not exist
            if variant_id is not None and sample_ids:
                if fld_cols is None:
                    col_rngs = [
                        reader.columnRange[x]
                        for x in range(self.ranges[3], self.ranges[4])
                    ]
                    fld_cols = []
                    for idx in range(len(sample_ids)):
                        fld_cols.append([
                            sc + (0 if sc + 1 == ec else idx)
                            for sc, ec in col_rngs
                        ])
                    if col_rngs[0][1] - col_rngs[0][0] != len(sample_ids):
                        env.logger.error(
                            'Number of genotypes ({}) does not match number of samples ({})'
                            .format(col_rngs[0][1] - col_rngs[0][0],
                                    len(sample_ids)))
                for idx, id in enumerate(sample_ids):
                    if id != -1 and rec[self.ranges[2] + idx] is not None:
                        cur.execute(genotype_update_query[id],
                                    [rec[c] for c in fld_cols[idx]] +
                                    [variant_id])
                        self.count[1] += 1
            if self.count[0] - last_count > update_after:
                last_count = self.count[0]
                self.db.commit()
                prog.update(self.count[0])
        self.count[7] = reader.unprocessable_lines
        self.db.commit()
        prog.done(self.count[0])

    def updateHDF5FromFile(self, input_filename):
        sample_ids = self.getSampleIDs(
            input_filename) if self.genotype_info else []
        reader = TextReader(
            self.processor,
            input_filename,
            # in the case of variant, we filter from the reading stage to save some time
            None if (self.table == 'variant' or self.input_type != 'variant')
            else self.variantIndex,
            # getNew is False so we only get variants that are available in variantIndex
            False,
            self.jobs - 1,
            self.encoding,
            self.header)
        cur = self.db.cursor()
        lc = lineCount(input_filename, self.encoding)
        update_after = min(max(lc // 200, 100), 100000)
        fld_cols = None
        prog = ProgressBar(os.path.split(input_filename)[-1], lc)
        last_count = 0

        for self.count[0], bins, rec in reader.records():
            variant_id = self.updateVariant(cur, bins, rec[0:self.ranges[2]])
            # variant might not exist
            if variant_id is not None and sample_ids:
                if fld_cols is None:
                    col_rngs = [
                        reader.columnRange[x]
                        for x in range(self.ranges[3], self.ranges[4])
                    ]
                    fld_cols = []
                    for idx in range(len(sample_ids)):
                        fld_cols.append([
                            sc + (0 if sc + 1 == ec else idx)
                            for sc, ec in col_rngs
                        ])

            if self.count[0] - last_count > update_after:
                last_count = self.count[0]
                self.db.commit()
                prog.update(self.count[0])
        self.count[7] = reader.unprocessable_lines
        hdf5Files = []
        for id in sample_ids:
            cur.execute('SELECT HDF5 FROM sample WHERE sample_id = ?;', (id,))
            res = cur.fetchone()
            hdf5file = res[0]
            hdf5Files.append(hdf5file)
        hdf5Files = list(set(hdf5Files))
        self.db.commit()
        prog.done(self.count[0])
        UpdateGenotypeInParallel(self, input_filename, sample_ids, hdf5Files)

    def update(self):
        '''Start updating'''
        for count, f in enumerate(self.files):

            env.logger.info('{} variants from {} ({}/{})'.format(
                'Updating', f, count + 1, len(self.files)))

            if self.proj.store == "sqlite" or (self.proj.store == "hdf5" and
                                               "vcf" not in f):
                self.updateFromFile(f)
                env.logger.info(
                    'Field{} {} of {:,} variants{} are updated'.format(
                        '' if len(self.variant_info) == 1 else 's',
                        ', '.join(self.variant_info), self.count[8],
                        '' if self.count[1] == 0 else
                        ' and geno fields of {:,} genotypes'.format(
                            self.count[1])))

            elif self.proj.store == "hdf5":
                self.updateHDF5FromFile(f)
            self.sample_name = ''
            for i in range(len(self.count)):
                self.total_count[i] += self.count[i]
                self.count[i] = 0

        # if len(self.files) > 1:
        #     env.logger.info('Field{} {} of {:,} variants{} are updated'.format('' if len(self.variant_info) == 1 else 's', ', '.join(self.variant_info), self.total_count[8],
        #             '' if self.total_count[1] == 0 else ' and geno fields of {:,} genotypes'.format(self.total_count[1])))


def setFieldValue(proj, table, items, build):
    # fields
    expr = ','.join([x.split('=', 1)[1] for x in items])
    select_clause, fields = consolidateFieldName(
        proj, table, expr, build and build == proj.alt_build)
    desc_of_field = {
        x.split('=', 1)[0].lower(): x.split('=', 1)[1] for x in items
    }
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
        env.logger.trace('Running {}'.format(query))
        cur.execute(query)
        fldTypes = [None] * len(items)
        for rec in cur:
            for i in range(len(items)):
                if rec[i] is None:  # missing
                    continue
                elif fldTypes[i] is None:
                    fldTypes[i] = type(rec[i])
                    continue
                if type(rec[i]) != fldTypes[i]:
                    if type(rec[i]) is float and fldTypes[i] is int:
                        fldTypes[i] = float
                    else:
                        raise ValueError(
                            'Inconsistent type returned from different samples')
            if all([x is not None for x in fldTypes]):
                break
        #
        count = [0] * 3
        # if adding a new field
        cur_fields = proj.db.getHeaders('variant')[3:]
        type_map = {
            int: 'INT',
            float: 'FLOAT',
            str: 'VARCHAR(255)',
            str: 'VARCHAR(255)'
        }
        for field, fldType in zip([x.split('=', 1)[0] for x in items],
                                  fldTypes):
            if field.lower() not in [x.lower() for x in cur_fields]:
                if fldType is None:
                    env.logger.warning(
                        'Use type VARCHAR for a new field {} because the values are all NULL'
                        .format(field))
                    fldType = str
                proj.checkFieldName(field, exclude=table)
                env.logger.info(
                    'Adding variant info field {} with type {}'.format(
                        field, type_map[fldType]))
                query = 'ALTER TABLE variant ADD {} {} NULL;'.format(
                    field, type_map[fldType])
                env.logger.trace(query)
                cur.execute(query)
                count[1] += 1  # new
            else:
                # FIXME: check the case for type mismatch
                env.logger.info('Updating field {}'.format(field))
                count[2] += 1  # updated
                # add an description
            try:
                proj.describeField(
                    field, 'Evaluated from "{}" with type {} on {}'.format(
                        desc_of_field[field.lower()], type_map[fldType],
                        time.strftime('%b%d', time.localtime())))
            except Exception as e:
                env.logger.warning(
                    'Failed to add a description to field {}: {}'.format(
                        field, e))
        proj.db.commit()
        # really update things
        # update clause is more difficult because it needs to consolidate the
        # right hand side, not the left hand side of expression
        update_clause = []
        for item in items:
            k, v = item.split('=', 1)
            v = consolidateFieldName(proj, table, v, build and
                                     build == proj.alt_build)[0]
            update_clause.append('{}={}'.format(k, v))
        if table == 'variant':
            query = 'UPDATE variant SET {};'.format(', '.join(update_clause))
        else:
            query = 'UPDATE variant SET {} WHERE variant_id IN (SELECT variant_id FROM {});'.format(
                ', '.join(update_clause), table)
        env.logger.trace('Running {}'.format(query))
        cur.execute(query)
    else:
        query = 'SELECT {}, {}.variant_id {};'.format(select_clause, table,
                                                      from_clause)
        env.logger.trace('Running {}'.format(query))
        cur.execute(query)
        fldTypes = [None] * len(items)
        with delayedAction(env.logger.info, 'Evaluating all expressions'):
            results = cur.fetchall()
        #
        for res in results:
            for i in range(len(items)):
                if res[i] is None:  # missing
                    continue
                elif fldTypes[i] is None:
                    fldTypes[i] = type(res[i])
                    continue
                if type(res[i]) != fldTypes[i]:
                    if type(res[i]) is float and fldTypes[i] is int:
                        fldTypes[i] = float
                    else:
                        raise ValueError(
                            'Inconsistent type returned from different samples')
            if all([x is not None for x in fldTypes]):
                break
        #
        count = [0] * 3
        # if adding a new field
        cur_fields = proj.db.getHeaders('variant')[1:]
        new_fields = [x.split('=', 1)[0] for x in items]
        for field, fldType in zip(new_fields, fldTypes):
            if field.lower() not in [x.lower() for x in cur_fields]:
                if fldType is None:
                    env.logger.warning(
                        'Use type VARCHAR for a new field {} '
                        'because the values are all NULL'.format(field))
                    fldType = str
                proj.checkFieldName(field, exclude=table)
                env.logger.info('Adding variant info field {}'.format(field))
                query = 'ALTER TABLE variant ADD {} {} NULL;'.format(
                    field, {
                        int: 'INT',
                        float: 'FLOAT',
                        str: 'VARCHAR(255)',
                        str: 'VARCHAR(255)'
                    }[fldType])
                env.logger.trace(query)
                cur.execute(query)
                count[1] += 1  # new
            else:
                # FIXME: check the case for type mismatch
                env.logger.info('Updating field {}'.format(field))
                count[2] += 1  # updated
        proj.db.commit()
        # really update things
        query = 'UPDATE variant SET {} WHERE variant_id={};'.format(
            ','.join(['{}={}'.format(x, proj.db.PH) for x in new_fields]),
            proj.db.PH)
        env.logger.trace('Using query {}'.format(query))
        prog = ProgressBar('Updating {}'.format(decodeTableName(table)),
                           len(results))
        # this particular query will return a bunch of bogus NULL values for range-based databases,
        # which needs to be filtered out.
        count = 0
        count_multi_value = 0
        last_id = None
        last_values = []
        for res in results:
            if last_id is None:
                last_id = res[-1]
                last_values = [res]
            elif last_id != res[-1]:
                # handle last record
                values = [x for x in last_values if x[0] is not None]
                if len(values) > 1:
                    # check duplicates
                    values = list(set(values))
                    if len(values) > 1:
                        env.logger.debug(
                            'Multiple field values {} exist for variant with id {}'
                            .format(
                                ', '.join([
                                    '\t'.join([str(y)
                                               for y in x[:-1]])
                                    for x in values
                                ]), last_id))
                        count_multi_value += 1
                    cur.execute(query, values[0])
                elif len(values) == 0:
                    # only none result exists
                    cur.execute(query, last_values[-1])
                else:
                    cur.execute(query, values[0])
                last_id = res[-1]
                last_values = [res]
            else:  # if id is the same
                last_values.append(res)
            count += 1
            if count % 10000 == 0:
                prog.update(count)
        # for the last one.
        if last_id is not None:
            # handle last record
            values = [x for x in last_values if x[0] is not None]
            if len(values) > 1:
                values = list(set(values))
                if len(values) > 1:
                    env.logger.debug(
                        'Multiple field values {} exist for variant with id {}'
                        .format(
                            ', '.join([
                                '\t'.join([str(y)
                                           for y in x[:-1]])
                                for x in values
                            ]), last_id))
                    count_multi_value += 1
                    cur.execute(query, values[0])
                elif len(values) == 0:
                    cur.execute(query, last_values[-1])
                else:
                    cur.execute(query, values[0])
            elif len(values) == 1:
                cur.execute(query, values[0])
        prog.done()
        if count_multi_value != 0:
            env.logger.warning(
                'Multiple field values are available for {} variants. Arbitrary valid values are chosen.'
                .format(count_multi_value))


def calcSampleStat(proj, from_stat, samples, variant_table, genotypes):
    '''Count sample allele count etc for specified sample and variant table'''
    if not proj.isVariantTable(variant_table):
        raise ValueError('"Variant table {} does not exist.'.format(
            decodeTableName(variant_table)))
    #
    IDs = None
    ID_sex = None
    if samples:
        IDs = proj.selectSampleByPhenotype(' AND '.join(
            ['({})'.format(x) for x in samples]))
        if len(IDs) == 0:
            env.logger.info('No sample is selected (or available)')
            return
        else:
            env.logger.info('{} samples are selected'.format(len(IDs)))

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
        env.logger.warning('No statistics is specified')
        return

    # separate special functions...
    alt = hom = het = other = GT = missing = wtGT = mutGT = maf = None

    # keys to speed up some operations
    MEAN = 0
    SUM = 1
    MIN = 2
    MAX = 3
    operationKeys = {'avg': MEAN, 'sum': SUM, 'min': MIN, 'max': MAX}
    possibleOperations = list(operationKeys.keys())

    operations = []
    genotypeFields = []
    validGenotypeFields = []
    destinations = []
    fieldCalcs = []
    desc_of_field = {}

    cur = proj.db.cursor()
    if IDs is None:
        cur.execute('SELECT sample_id from sample;')
        IDs = [x[0] for x in cur.fetchall()]
    #
    numSample = len(IDs)
    if numSample == 0:
        env.logger.warning('No sample is selected.')
        return
    prog = None
    if proj.store == "sqlite":
        prog = ProgressBar('Counting variants', len(IDs))
    else:
        env.logger.info("Reading genotype info for processing....")

    #
    # Error checking with the user specified genotype fields
    # 1) if a field does not exist within one of the sample genotype tables a warning is issued
    # 2) if a field does not exist in any sample, it is not included in validGenotypeFields
    # 3) if no fields are valid and no core stats were requested (i.e., alt, het, hom, other), then sample_stat is exited
    genotypeFieldTypes = {}
    fieldInTable = defaultdict(list)
    store = GenoStore(proj)

    for id in IDs:
        fields = ["variant_id"]
        #starttime = time.time()
        fields.extend(store.geno_fields(id))
        for field in fields:
            fieldInTable[field].append(id)
            if field not in genotypeFieldTypes:
                genotypeFieldTypes[field] = 'INT'
                # fieldType = proj.db.typeOfColumn('{}_genotype.genotype_{}'.format(proj.name, id), field)
                fieldType = store.get_typeOfColumn(id, field)
                if fieldType.upper().startswith('FLOAT'):
                    genotypeFieldTypes[field] = 'FLOAT'
                elif fieldType.upper().startswith('VARCHAR'):
                    genotypeFieldTypes[field] = 'VARCHAR'
                    # We had been throwing an error here if a genotype field is a VARCHAR, but I think we should allow
                    # VARCHAR fields in the genotype tables.  We'll throw an error if someone wants to perform numeric operations on these fields
                    # further down in the code.
                    # raise ValueError('Genotype field {} is a VARCHAR which is not supported with sample_stat operations.'.format(field))
    validGenotypeIndices = []

    for stat in from_stat:
        f, e = [x.strip() for x in stat.split('=', 1)]
        desc_of_field[f.lower()] = e
        if f in ['chr', 'pos', 'ref', 'alt']:
            raise ValueError('Cannot overwrite existing field "{}"!'.format(f))
        if e == '#(alt)':
            alt = f
        elif e == '#(hom)':
            hom = f
        elif e == '#(het)':
            het = f
        elif e == '#(other)':
            other = f
        elif e == '#(GT)':
            GT = f
        elif e == '#(missing)':
            missing = f
        elif e == '#(wtGT)':
            wtGT = f
        elif e == '#(mutGT)':
            mutGT = f
        elif e == 'maf()':
            maf = f
        elif e.startswith('#('):
            raise ValueError(
                'Unrecognized parameter {}: only parameters alt, '
                'wtGT, mutGT, missing, hom, het, other and GT are accepted for '
                'special functions'.format(stat))
        else:
            m = re.match('(\w+)\s*=\s*(avg|sum|max|min)\s*\(\s*(\w+)\s*\)\s*',
                         stat)
            if m is None:
                raise ValueError(
                    'Unrecognized parameter {}, which should have the form of FIELD=FUNC(GENO_INFO) where FUNC is one of #, avg, sum, max and min'
                    .format(stat))
            dest, operation, field = m.groups()
            if proj.store == "hdf5":
                field = field.replace("_geno", "")
            if operation not in possibleOperations:
                raise ValueError(
                    'Unsupported operation {}. Supported operations include {}.'
                    .format(operation, ', '.join(possibleOperations)))

            if field.lower() in [
                    x.lower() for x in list(genotypeFieldTypes.keys())
            ]:
                operations.append(operationKeys[operation])
                genotypeFields.append(field)
                fieldCalcs.append(None)
                destinations.append(dest)

    # for maf calculation, we need to know sex and chromosome name information
    if maf is not None:
        # find variants on sex chromosomes (exlcuding pseudi-autosomal regions)
        var_chrX = getVariantsOnChromosomeX(proj, variant_table)
        var_chrY = getVariantsOnChromosomeY(proj, variant_table)
        var_chrOther = getVariantsOnManifolds(proj, variant_table)
        #
        # if there are any variants on sex chromosome, we need to know the
        # sex of samples
        if len(var_chrX) > 0 or len(var_chrY) > 0:
            ID_sex = determineSexOfSamples(proj, IDs)
            numMales = len({x: y for x, y in list(ID_sex.items()) if y == 1})
            numFemales = len({x: y for x, y in list(ID_sex.items()) if y == 2})
            env.logger.info('{} males and {} females are identified'.format(
                numMales, numFemales))

    coreDestinations = [alt, hom, het, other, GT, missing, wtGT, mutGT, maf]
    for index, field in enumerate(genotypeFields):
        if field.lower() not in [
                x.lower() for x in list(genotypeFieldTypes.keys())
        ]:
            env.logger.warning(
                "Field {} does not exist in any of the samples.".format(field))
        else:
            if len(fieldInTable[field.lower()]) < len(IDs):
                env.logger.warning(
                    'Field {} exists in {} of {} selected samples'.format(
                        field, len(fieldInTable[field.lower()]), len(IDs)))
            validGenotypeIndices.append(index)
            validGenotypeFields.append(field)

    # check GT field
    if not all([x is None for x in coreDestinations]):
        if 'gt' not in [x.lower() for x in list(genotypeFieldTypes.keys())]:
            env.logger.warning(
                'Genotype field does not exist in any of the selected samples')
        else:
            if len(fieldInTable['gt']) < len(IDs):
                env.logger.warning(
                    'Genotype field GT exists in {} of {} selected samples'
                    .format(len(fieldInTable[field.lower()]), len(IDs)))

    if all([x is None for x in coreDestinations
           ]) and len(validGenotypeFields) == 0:
        env.logger.warning(
            "No valid sample statistics operation has been specified.")
        return

    queryDestinations = coreDestinations
    for index in validGenotypeIndices:
        queryDestinations.append(destinations[index])
    for name in queryDestinations:
        if name is not None:
            proj.checkFieldName(name, exclude='variant')

    prog_step = max(len(IDs) // 100, 1)

    sampleDict = dict()
    for id_idx, id in enumerate(IDs):
        record_male_gt = ID_sex is not None and ID_sex[id] == 1
        fieldSelect = [
            'GT' if ('gt' in fieldInTable and id in fieldInTable['gt']) or
            ('GT' in fieldInTable and id in fieldInTable['GT']) else 'NULL'
        ]
        if validGenotypeFields is not None and len(validGenotypeFields) != 0:
            fieldSelect.extend([
                x if id in fieldInTable[x.lower()] else 'NULL'
                for x in validGenotypeFields
            ])
        sampleDict[id] = [record_male_gt, fieldSelect]
        if not fieldSelect or all([x == 'NULL' for x in fieldSelect]):
            continue
    variants = store.get_genoType_genoInfo(sampleDict, genotypes, variant_table,
                                           genotypeFields, validGenotypeIndices,
                                           validGenotypeFields, operations,
                                           fieldCalcs, prog, prog_step)

    #
    # even if no variant is updated, we need set count 0 to fields.
    #if len(variants) == 0:
    #    raise ValueError('No variant is updated')
    #
    headers = [x.lower() for x in proj.db.getHeaders('variant')]
    table_attributes = [(alt, 'INT'), (hom, 'INT'), (het, 'INT'),
                        (other, 'INT'), (GT, 'INT'), (missing, 'INT'),
                        (wtGT, 'INT'), (mutGT, 'INT'), (maf, 'FLOAT')]
    fieldsDefaultZero = [alt, hom, het, other, GT, missing, wtGT, mutGT, maf]
    for index in validGenotypeIndices:
        field = genotypeFields[index]
        genotypeFieldType = genotypeFieldTypes.get(
            genotypeFields[index].lower())
        # genotypeFieldType = genotypeFieldTypes.get(genotypeFields[index])
        if genotypeFieldType == 'VARCHAR':
            raise ValueError(
                'Genotype field {} is a VARCHAR which is not supported with sample_stat operations.'
                .format(field))

        if operations[index] == MEAN:
            table_attributes.append((destinations[index], 'FLOAT'))
        else:
            table_attributes.append((destinations[index], genotypeFieldType))

    for field, fldtype in table_attributes:
        if field is None:
            continue
        # We are setting default values on the count fields to 0. The genotype stat fields are set to NULL by default.
        defaultValue = 0 if field in fieldsDefaultZero else None
        if variant_table != 'variant':
            where_clause = 'WHERE variant_id in (SELECT variant_id FROM {})'.format(
                variant_table)
        else:
            where_clause = ''

        if field.lower() in headers:
            if proj.db.typeOfColumn('variant', field) != (fldtype + ' NULL'):
                env.logger.warning('Type mismatch (existing: {}, new: {}) for column {}. Please remove this column and recalculate statistics if needed.'\
                    .format(proj.db.typeOfColumn('variant', field), fldtype, field))
            env.logger.info(
                'Resetting values at existing field {}'.format(field))

            proj.db.execute(
                'Update variant SET {} = {} {};'.format(field, proj.db.PH,
                                                        where_clause),
                (defaultValue,))
        else:
            env.logger.info('Adding variant info field {} with type {}'.format(
                field, fldtype))
            proj.db.execute('ALTER TABLE variant ADD {} {} NULL;'.format(
                field, fldtype))
            if defaultValue == 0:
                proj.db.execute('UPDATE variant SET {} = 0 {}'.format(
                    field, where_clause))
        # add an description

        try:
            proj.describeField(
                field, 'Created from stat "{}" {} with type {} on {}'.format(
                    desc_of_field[field.lower()],
                    'for samples {}'.format(samples) if samples else '',
                    fldtype, time.strftime('%b%d', time.localtime())))
        except Exception as e:
            env.logger.warning(
                'Failed to add a description to field {}: {}'.format(field, e))
    #
    prog = ProgressBar('Updating {}'.format(decodeTableName(variant_table)),
                       len(variants))
    update_query = 'UPDATE {0} SET {2} WHERE variant_id={1};'.format(
        'variant', proj.db.PH, ' ,'.join([
            '{}={}'.format(x, proj.db.PH)
            for x in queryDestinations
            if x is not None
        ]))

    count = 0
    for count, id in enumerate(variants):
        value = variants[id]
        res = []
        if alt is not None:
            # het + hom * 2 + other
            res.append(value[0] + value[1] * 2 + value[2])
        if hom is not None:
            res.append(value[1])
        if het is not None:
            res.append(value[0])
        if other is not None:
            res.append(value[2])
        if GT is not None:
            res.append(value[3])
        if missing is not None:
            # missing = # sample - #(GT)
            res.append(numSample - value[3])
        if wtGT is not None:
            # wtGT = #(GT) - hom - het - other
            res.append(value[3] - value[0] - value[1] - value[2])
        if mutGT is not None:
            # mutGT = hom + het + other
            res.append(value[0] + value[1] + value[2])
        if maf is not None:
            # numerator is the number of alternative alleles (het + hom * 2 + other)
            # denominator is total number of genotypes
            numerator = value[0] + value[1] * 2 + value[2]
            if id in var_chrX:
                if env.treat_missing_as_wildtype:
                    denominator = numFemales * 2 + numMales
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'chromosome X: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#male={}+2*#female={})'.format(
                                id, numerator, denominator, value[0], value[1],
                                value[2], numMales, numFemales))
                else:
                    denominator = (value[3] - value[4]) * 2 + value[4]
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'chromosome X: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#maleGT={}+2*#femaleGT={})'.format(
                                id, numerator, denominator, value[0], value[1],
                                value[2], value[4], value[3] - value[4]))
            elif id in var_chrY:
                if env.treat_missing_as_wildtype:
                    denominator = numMales
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'chromosome Y: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#male={})'.format(id, numerator,
                                                        denominator, value[0],
                                                        value[1], value[2],
                                                        numMales))
                else:
                    denominator = value[4]
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'chromosome X: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#maleGT={})'.format(
                                id, numerator, denominator, value[0], value[1],
                                value[2], value[4]))
            elif id in var_chrOther:
                if env.treat_missing_as_wildtype:
                    denominator = numSample
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'other chromosome: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#sample={})'.format(
                                id, numerator, denominator, value[0], value[1],
                                value[2], numSample))
                else:
                    denominator = value[3]
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'chromosome X: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=#GT={})'.format(id, numerator,
                                                      denominator, value[0],
                                                      value[1], value[2],
                                                      value[3]))
            else:
                # regular autosome
                if env.treat_missing_as_wildtype:
                    denominator = numSample * 2
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'autosome: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=2*#sample={})'.format(
                                id, numerator, denominator, value[0], value[1],
                                value[2], numSample))
                else:
                    denominator = value[3] * 2
                    if numerator > denominator:
                        env.logger.warning(
                            'Incorrect MAF for variant {} on '
                            'autosome: #alt={} > #alleles={} (#alt=#het={}+2*#hom={}+#other={} '
                            '#alleles=2*#GT={})'.format(id, numerator,
                                                        denominator, value[0],
                                                        value[1], value[2],
                                                        value[3]))
            ratio = 0 if denominator == 0 else numerator * 1.0 / denominator
            res.append(ratio if ratio < 0.5 else 1. - ratio)
        # for genotype_field operations, the value[operation_index] holds the result of the operation
        # except for the "mean" operation which needs to be divided by num_samples that have that variant
        try:
            for index in validGenotypeIndices:
                # the first 4 indices hold the values for hom, het, double het, total genotype and total genotype in males
                operationIndex = index + 5
                operationCalculation = value[operationIndex]
                # #temp fix
                # if proj.store=="hdf5":
                #     operationCalculation = value[operationIndex:]
                if operations[
                        index] == MEAN and operationCalculation is not None:
                    # res.append(float(operationCalculation[0]) / operationCalculation[1])
                    if value[3] != 0:
                        res.append(float(operationCalculation) / value[3])
                else:
                    res.append(operationCalculation)
            cur.execute(update_query, res + [id])
        except Exception as e:
            env.logger.debug(e)
        if count % 10000 == 0:
            proj.db.commit()
            prog.update(count)
    #
    # missing should be handled differently because variants only stores variants that
    # show up in at least one samples. If a variant does not exist in any of them, the
    # above queries will ignore them.
    if missing is not None:
        # first get all variants
        cur.execute('SELECT variant_id FROM {};'.format(variant_table))
        missing_variants = [
            x[0] for x in cur.fetchall() if x[0] not in variants
        ]
        if missing_variants:
            update_query = 'UPDATE {0} SET {2}={3} WHERE variant_id={1};'.format(
                'variant', proj.db.PH, missing, numSample)
            for id in missing_variants:
                cur.execute(update_query, [id])
                count += 1
                if count % 10000 == 0:
                    prog.update(count)
    prog.done()
    # done
    # result=cur.execute('SELECT * FROM {};'.format(variant_table))
    # for res in result:
    #     print(res)

    proj.db.commit()
    if proj.store == "sqlite":
        env.logger.info('{} records are updated'.format(count + 1))


def updateArguments(parser):
    parser.add_argument('table', help='''variants to be updated.''')
    files = parser.add_argument_group('Update from files')
    files.add_argument(
        '--from_file',
        '--from-file',
        nargs='+',
        help='''A list of files that will be used to add or update existing fields of
            variants. The file should be delimiter separated with format described by
            parameter --format. Gzipped files are acceptable. If input files contains
            genotype information, have been inputted before, and can be linked to the
            samples they created without ambiguity (e.g. single sample, or samples with
            detectable sample names), genotypes and their information will also be
            updated.''')
    files.add_argument(
        '--build',
        help='''Build version of the reference genome (e.g. hg18) of the input files,
            which should be the primary (used by default) or alternative (if available)
            reference genome of the project. An alternative reference genome will be
            added to the project if needed.''')
    files.add_argument(
        '--format',
        help='''Format of the input text file. It can be one of the variant tools
            supported file types such as ANNOVAR_output (cf. 'vtools show formats'),
            or a local format specification file (with extension .fmt). Some formats
            accept parameters (cf. 'vtools show format FMT') and allow you to update
            additional or alternative fields from the input file.''')
    files.add_argument(
        '-j',
        '--jobs',
        metavar='N',
        default=1,
        type=int,
        help='''Number of processes to import input file. Variant tools by default
            uses a single process for reading and writing, and can use one or more
            dedicated reader processes (jobs=2 or more) to process input files. Due
            to the overhead of inter-process communication, more jobs do not
            automatically lead to better performance.''')
    files.add_argument(
        '--sample_name',
        '--sample-name',
        nargs='*',
        default=[],
        help='''Name of the samples to be updated by the input files. If unspecified,
            headers of the genotype columns of the last comment line (line starts
            with #) of the input files will be used (and thus allow different sample
            names for input files). Sample names will be used to identify samples to
            be updated. Filename will be used to uniquely identify a sample if mutliple
            samples with the same name exist in the project. No genotype info will
            be updated if samples cannot be unquely determined.''')
    field = parser.add_argument_group('Set value from existing fields')
    field.add_argument(
        '--set',
        metavar='EXPR',
        nargs='*',
        default=[],
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
    stat.add_argument(
        '--from_stat',
        '--from-stat',
        metavar='EXPR',
        nargs='*',
        default=[],
        help='''One or more expressions such as meanQT=avg(QT) that aggregate
            genotype info (e.g. QT) of variants in all or selected samples to
            specified fields (e.g. meanQT). Functions sum, avg, max, and min
            are currently supported. In addition, special functions #(GT),
            #(hom), #(het), #(alt), #(other), #(missing), #(wtGT), #(mutGT),
            and maf(), are provided to count the number of valid genotypes (not
            missing), homozygote genotypes, heterozygote genotypes, alternative
            alleles (#(het) + 2*#(hom) + #(other)), genotypes with two different
            alternative alleles, missing genotypes (number of samples - #(GT)),
            number of non-missing wildtype genotypes (#(GT) - #(hom) - #(het)
            - #(other)), number of non-wildtype genotypes (#(hom) + #(het) +
            #(other)), and minor allele frequency. The maf() function treats
            chromosomes 1 to 22 as autosomes, X and Y as sex chromosomes, and
            other chromosomes as single-copy manifolds. It requires a phenotype
            named sex or gender that codes male/female by 1/2, M/F or Male/Female
            if maf of variants on sex chromosomes are calculated. This function
            by default calculates allele frequency among existing-alleles, but
            will treat all missing values as wild type alleles if runtime option
            treat_missing_as_wildtype is set to true.''')
    stat.add_argument(
        '-s',
        '--samples',
        nargs='*',
        metavar='COND',
        default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    stat.add_argument(
        '--genotypes',
        nargs='*',
        metavar='COND',
        default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show genotypes' (e.g. 'GQ>15').'''
    )


def update(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if not args.from_file and not args.from_stat and not args.set:
                raise ValueError('Please specify one of --from_file, --set '
                                 'and --from_stat for command vtools upate')
            if args.from_file:
                if (proj.store == "sqlite"):
                    proj.db.attach(proj.name + '_genotype')
                updater = Updater(
                    proj=proj,
                    table=encodeTableName(args.table),
                    files=args.from_file,
                    build=args.build,
                    format=args.format,
                    jobs=args.jobs,
                    sample_name=args.sample_name,
                    fmt_args=args.unknown_args)
                updater.update()
            if args.set:
                for item in args.set:
                    try:
                        field, expr = [x.strip() for x in item.split('=', 1)]
                    except Exception as e:
                        raise ValueError(
                            'Invalid parameter {}, which should have format field=expr_of_field: {}'
                            .format(item, e))
                setFieldValue(proj, encodeTableName(args.table), args.set,
                              proj.build)
                #, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.from_stat:
                if (proj.store == "sqlite"):
                    proj.db.attach(proj.name + '_genotype')
                variant_table = encodeTableName(
                    args.table) if args.table else 'variant'
                if not proj.db.hasTable(variant_table):
                    raise ValueError('Variant table {} does not exist'.format(
                        decodeTableName(variant_table)))
                calcSampleStat(proj, args.from_stat, args.samples,
                               variant_table, args.genotypes)
        proj.close()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
