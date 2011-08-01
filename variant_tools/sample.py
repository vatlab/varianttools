#!/usr/bin/env python
#
# $File: sample.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
from .project import Project
from .utils import ProgressBar, typeOfValues

class Sample:
    def __init__(self, proj):
        self.proj = proj
        self.logger = proj.logger
        self.db = proj.db

    def load(self, filename):
        '''Load phenotype informaiton from a file'''
        if not self.db.hasTable('sample'):
            self.logger.warning('Project does not have a sample table.')
            return
        with open(filename) as input:
            headers = input.readline().rstrip().split('\t')
            if len(headers) < 2 or headers[0] != 'filename' or headers[1] != 'sample_name':
                raise ValueError('The phenotype file must start with a header line with first two fields being filename and sample')
            # FIXME: check if headers are allowed because they will be names of tables.
            #
            nFields = len(headers)
            if nFields == 2:
                self.logger.warning('No phenotype is defined.')
                return
            record = {}
            for line in input.readlines():
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.rstrip().split('\t')
                if len(fields) != nFields:
                    raise ValueError('Invalid phenotype file: number of fields mismatch.')
                #
                record[(fields[0], fields[1])] = fields[2:]
        # now, try to guess the type
        types = []
        for i in range(nFields - 2):
            types.append(typeOfValues([x[i] for x in record.values()]))
        #
        # get existing table
        cur = self.db.cursor()
        cur.execute('SELECT sample.sample_id, filename.filename, sample.file_id,\
            sample.sample_name FROM sample, filename WHERE sample.file_id = filename.file_id;')
        table = []
        for rec in cur:
            if (rec[1], rec[3]) not in record:
                raise ValueError('No phenotype is defined for sample {} in file {}'.format(rec[3], rec[1]))
            table.append([rec[0], rec[2], rec[3]] + record[(rec[1], rec[3])])
        self.logger.info('Importing phenotypes into table sample')
        # check field names
        for name in headers[2:]:
            self.proj.checkFieldName(name, exclude='sample')
        # create a new table
        temp_table_name = '_tmp_sample'
        try:
            self.proj.createSampleTableIfNeeded(zip(headers[2:], types), table=temp_table_name)
            # insert rows
            query = 'INSERT INTO {} VALUES ({});'.format(temp_table_name, ','.join([self.db.PH]*(nFields+1)))
            for line in table:
                cur.execute(query, line)
        except Exception as e:
            self.proj.db.removeTable(temp_table_name)
            self.logger.error('Failed to create sample table. Perhaps an invalid field name is specified.')
            self.logger.debug(e)
            raise e
        # everything is OK.
        self.proj.db.removeTable('sample')
        self.proj.db.renameTable(temp_table_name, 'sample')
        self.db.commit()

    def selectSampleByPhenotype(self, cond):
        '''Select samples by conditions such as "aff=1"'''
        cur = self.db.cursor()
        try:
            query = 'SELECT sample_id FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id WHERE {};'.format(cond)
            self.logger.debug('Select samples using query')
            self.logger.debug(query)
            cur.execute(query)
            return set([x[0] for x in cur.fetchall()])
        except Exception as e:
            raise ValueError('Failed to retrieve samples by condition "{}"'.format(cond))

    def createVariantTableBySample(self, IDs, from_table, to_table):
        '''Create a variant table from genotype_table'''
        if not self.proj.isVariantTable(from_table):
            raise ValueError('Variant table {} does not exist'.format(from_table))
        if self.db.hasTable(to_table):
            self.logger.warning('Overwriting existing variant table {}, which can be slow for sqlite3 database.'.format(to_table))
            self.db.removeTable(to_table)
        if len(IDs) == 0:
            self.logger.warning('No ID is selected. Create an empty variant table {}'.format(to_table))
            self.proj.createVariantTable(to_table)
            return
        #
        cur = self.db.cursor()
        from_variants = set()
        if from_table != 'variant':
            self.logger.info('Getting variants from {}'.format(from_table))
            cur.execute('SELECT variant_id FROM {};'.format(from_table))
            from_variants = set([x[0] for x in cur.fetchall()])
        #
        prog = ProgressBar('Selecting variants',
            sum([self.db.numOfRows('{}_genotype.sample_variant_{}'.format(self.proj.name, id)) for id in IDs]))
        to_variants = set()
        count = 0
        for id in IDs:
            cur.execute('SELECT variant_id from {}_genotype.sample_variant_{};'.format(self.proj.name, id))
            for rec in cur:
                if len(from_variants) == 0 or rec[0] in from_variants:
                    to_variants.add(rec[0])
                count += 1
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        #
        self.proj.createVariantTable(to_table)
        query = 'INSERT INTO {} VALUES ({});'.format(to_table, self.db.PH)
        prog = ProgressBar('Generating table {}'.format(to_table), len(to_variants))
        # sort variant_id so that variant_id will be in order, which might
        # improve database performance
        for count,id in enumerate(sorted(to_variants)):
            cur.execute(query, (id,))
            if count % self.db.batch == 0:
                self.db.commit()
                prog.update(count)
        self.db.commit()
        prog.done()

    def calcSampleStat(self, IDs, variant_table, num, freq, hom, het, other, depth):
        '''Count sample allele frequency etc for specified sample and variant table'''
        if not self.proj.isVariantTable(variant_table):
            raise ValueError('"Variant_table {} does not exist.'.format(variant_table))
        #
        if num is None and freq is None and hom is None and het is None and other is None and depth is None:
            self.logger.warning('No statistics is specified')
            return
        #
        for name in (num, freq, hom, het, other, depth):
            if name is not None:
                self.proj.checkFieldName(name, exclude=variant_table)
        #
        cur = self.db.cursor()
        if IDs is None:
            cur.execute('SELECT sample_id from sample;')
            IDs = [x[0] for x in cur.fetchall()]
        #
        numSample = len(IDs)
        if numSample == 0:
            return
        #
        from_variants = set()
        if variant_table != 'variant':
            self.logger.info('Getting variants from table {}'.format(variant_table))
            cur.execute('SELECT variant_id FROM {};'.format(variant_table))
            from_variants = set([x[0] for x in cur.fetchall()])
        #
        # too bad that I can not use a default dict...
        variants = dict()
        prog = ProgressBar('Counting variants',
            sum([self.db.numOfRows('{}_genotype.sample_variant_{}'.format(self.proj.name, id)) for id in IDs]))
        count = 0
        for id in IDs:
            cur.execute('SELECT * FROM {}_genotype.sample_variant_{};'.format(self.proj.name, id))
            for rec in cur:
                if len(from_variants) == 0 or rec[0] in from_variants:
                    if rec[0] not in variants:
                        variants[rec[0]] = [0, 0, 0, 0] if depth is not None else [0, 0, 0]
                    if rec[1] == 1:
                        variants[rec[0]][0] += 1
                    elif rec[1] == 2:
                        variants[rec[0]][1] += 1
                    elif rec[1] == -1:
                        variants[rec[0]][2] += 1
                    else:
                        self.logger.warning('Invalid genotype type {}'.format(rec[1]))
                    if depth is not None:
                        try:
                            # FIXME: Can some depth information be missing? If so , we should test
                            # if rec[2] is None
                            variants[rec[0]][3] += rec[2]
                        except Exception as e:
                            self.logger.error('Cannot calculate mean depth. Is the vcf files imported with --import_depth option?')
                            raise e
                count += 1
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        #
        headers = self.db.getHeaders(variant_table)
        for field, fldtype in [(num, 'INT'), (freq, 'FLOAT'), (hom, 'INT'),
                (het, 'INT'), (other, 'INT'), (depth, 'FLOAT')]:
            if field is None:
                continue
            if field in headers:
                # NOTE: there is a possible problem of type mismatch 
                # e.g. saving frequency to an integer field
                self.logger.info('Updating existing field {}'.format(field))
                if fldtype == 'FLOAT':
                    self.logger.warning('Result will be wrong if field \'{}\' was created to hold integer values'.format(field))
            else:
                self.logger.info('Adding field {}'.format(field))
                self.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(variant_table, field, fldtype))
        #
        prog = ProgressBar('Updating table {}'.format(variant_table), len(variants))
        update_query = 'UPDATE {0} SET {2} WHERE variant_id={1};'.format(variant_table, self.db.PH,
            ' ,'.join(['{}={}'.format(x, self.db.PH) for x in [num, freq, hom, het, other, depth] if x is not None]))
        for count,id in enumerate(variants):
            value = variants[id]
            res = []
            if num is not None:
                res.append(value[0] + value[1] * 2 + value[2])
            if freq is not None:
                res.append((value[0] + value[1] * 2 + value[2])/(2. * numSample))
            if hom is not None:
                res.append(value[0])
            if het is not None:
                res.append(value[1])
            if other is not None:
                res.append(value[2])
            if depth is not None:
                res.append(value[3] / (value[0] + value[1] + value[2]))
            cur.execute(update_query, res + [id])
            if count % self.db.batch == 0:
                self.db.commit()
                prog.update(count)
        self.db.commit()
        prog.done()
                
def importPhenotypeArguments(parser):
    '''Action that can be performed by this script'''
    parser.add_argument('filename', 
        help='''Import phenotype from a tab delimited file. The file should have
            a header, and filename and sample_name as the first two columns.''')

def importPhenotype(args):
    try:
        with Project() as proj:
            p = Sample(proj)
            p.load(args.filename)
        proj.close()
    except Exception as e:
        sys.exit(e)
                
def subsampleArguments(parser):
    '''Arguments to select variants from certain sample'''
    parser.add_argument('condition',
        help='''Criteria by which samples are chosen. This parameter should be
            a SQL expression using one or more fields shown in 'vtools show sample'
            (e.g. 'aff=1' and 'filename like "MG%%"'. ''')
    parser.add_argument('-f', '--from_table', default='variant',
        help='''Source variant table. Default to the master variant table.''')
    parser.add_argument('-t', '--to_table', required=True,
        help='''Name of a new variant table.''')
    
def subsample(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            p = Sample(proj)
            # we save genotype in a separate database to keep the main project size tolerable.
            proj.db.attach(proj.name + '_genotype')
            IDs = p.selectSampleByPhenotype(args.condition)
            p.logger.info('{} samples are selected by condition {}'.format(len(IDs), args.condition))
            if len(IDs) == 0:
                return
            # create variant table
            p.createVariantTableBySample(IDs, args.from_table, args.to_table)
        # temporary tables will be removed
        proj.close()
    except Exception as e:
        sys.exit(e)


def sampleStatArguments(parser):
    '''Arguments to calculate sample statistics such as allele frequency'''
    parser.add_argument('-s', '--samples',
        help='''Criteria by which samples are chosen, the same as parameter 'condition'
            for command 'vtools subsample'. By default, all samples will be used.''')
    parser.add_argument('-t', '--table',
        help='''Variant table. The master variant table will be updated by default''')
    parser.add_argument('-n', '--num',
        help='''Name of the field to hold number of alternative alleles in the sample.''')
    parser.add_argument('-f', '--freq',
        help='''Name of the field to hold allele frequency, incorrect for variants on sex chromosomes.''')
    parser.add_argument('--hom',
        help='''Name of the field to hold number of samples with two identical alternative alleles.''')
    parser.add_argument('--het',
        help='''Name of the field to hold number of samples with one reference and one alternative alleles.''')
    parser.add_argument('--other',
        help='''Name of the field to hold number of samples with two different alternative alleles.''')
    parser.add_argument('--depth',
        help='''Name of the field to hold mean depth across samples. Different genotype
            types are given the same weight. To use this option, the vcf files must be
            imported with option --import_depth turned on.''')
    
def sampleStat(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # we save genotype in a separate database to keep the main project size tolerable.
            proj.db.attach(proj.name + '_genotype')
            variant_table = args.table if args.table else 'variant'
            if not proj.db.hasTable(variant_table):
                raise ValueError('Variant table {} does not exist'.format(variant_table))
            p = Sample(proj)
            IDs = None
            if args.samples:
                IDs = p.selectSampleByPhenotype(args.samples)
                if len(IDs) == 0:
                    p.logger.info('No sample is selected (or available)')
                    return
                else:
                    p.logger.info('{} samples are selected'.format(len(IDs)))
            p.calcSampleStat(IDs, variant_table, args.num, args.freq, args.hom,
                args.het, args.other, args.depth)
        # temporary tables will be removed
        proj.close()
    except Exception as e:
        sys.exit(e)

