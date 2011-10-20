#!/usr/bin/env python
#
# $File: sample.py $
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
from collections import defaultdict
from .project import Project
from .utils import ProgressBar, typeOfValues

class Sample:
    def __init__(self, proj):
        self.proj = proj
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
            self.proj.checkFieldName(field, exclude='sample')
            self.logger.info('Adding field {}'.format(field))
            self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(field,
                {int: 'INT',
                 float: 'FLOAT',
                 str: 'VARCHAR(255)',
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
        count = [0, 0, 0]
        # if adding a new field
        cur_fields = self.db.getHeaders('sample')[3:]
        new_field = {}
        for field in [x[0] for x in stat]:
            if field.lower() not in [x.lower() for x in cur_fields]:
                self.proj.checkFieldName(field, exclude='sample')
                new_field[field] = True
            else:
                new_field[field] = False
                count[2] += 1  # updated
        #
        cur = self.db.cursor()
        prog = ProgressBar('Calculating phenotype', len(IDs))
        for prog_count, ID in enumerate(IDs):
            try:
                res = [None] * len(stat)
                query = 'SELECT {} FROM {}_genotype.genotype_{} {};'\
                    .format(', '.join([x[1] for x in stat]), self.proj.name, ID,
                        'WHERE {}'.format(genotypes) if genotypes.strip() else '')
                self.logger.debug(query)
                try:
                    cur.execute(query)
                    res = cur.fetchone()
                except Exception as e:
                    # some field might not exist, so we will have to run one by one
                    for idx, (field, expr) in enumerate(stat):
                        query = 'SELECT {} FROM {}_genotype.genotype_{} {};'\
                            .format(expr, self.proj.name, ID,
                                'WHERE {}'.format(genotypes) if genotypes.strip() else '')
                        self.logger.debug(query)
                        try:
                            cur.execute(query)
                            v = cur.fetchone()
                            if v is not None:
                                res[idx] = v[0]
                        except Exception as e:
                            self.logger.debug('Failed: {}. Setting field {} to NULL.'.format(e, field))
                for val, (field, expr) in zip(res, stat):
                    if val is None:
                        if not new_field:
                            cur.execute('UPDATE sample SET {0}={1} WHERE sample_id = {1}'.format(field, self.db.PH), [None, ID])
                    else:
                        if new_field[field]:
                            self.logger.debug('Adding field {}'.format(field))
                            # determine the type of value from the first one.
                            self.db.execute('ALTER TABLE sample ADD {} {} NULL;'.format(field,
                                {int: 'INT',
                                 float: 'FLOAT',
                                 str: 'VARCHAR(255)',
                                 None: 'FLOAT'}[type(val)]))
                            new_field[field] = False
                            count[1] += 1  # new
                        cur.execute('UPDATE sample SET {0}={1} WHERE sample_id = {1}'.format(field, self.db.PH), [val, ID])
                        count[0] += 1
                prog.update(prog_count + 1)
            except Exception as e:
                self.logger.info('Updating phenotype {} failed: {}'.format(field, e))
                #cur.execute('UPDATE sample SET {0}={1} WHERE sample_id = {1}'.format(field, self.db.PH), [None, ID])
        prog.done()
        self.logger.info('{} values of {} phenotypes ({} new, {} existing) of {} samples are updated.'.format(
            count[0], count[1]+count[2], count[1], count[2], len(IDs)))
        self.db.commit()


    def calcSampleStat(self, IDs, variant_table, genotypes, num, hom, het, other, other_stats):
        '''Count sample allele count etc for specified sample and variant table'''
        if not self.proj.isVariantTable(variant_table):
            raise ValueError('"Variant_table {} does not exist.'.format(variant_table))
        
        if num is None and hom is None and het is None and other is None and not other_stats:
            self.logger.warning('No statistics is specified')
            return

        coreDestinations = [num, hom, het, other]  

        # keys to speed up some operations
        MEAN = 0
        SUM = 1
        MIN = 2
        MAX = 3
        operationKeys = {'mean': MEAN, 'sum': SUM, 'min': MIN, 'max': MAX}
        possibleOperations = operationKeys.keys()
        
        operations = []
        genotypeFields = []
        validGenotypeFields = []
        destinations = []
        fieldCalcs = []   
        for index in range(0, len(other_stats), 2):
            if other_stats[index].startswith('--'):
                if other_stats[index].find('_') == -1:
                    raise ValueError('Unsupported operation {}.  Supported operations include {}.'.format(other_stats[index][2:], ', '.join(possibleOperations)))
                operation, field = other_stats[index][2:].split('_',1)
                if operation not in possibleOperations:
                    raise ValueError('Unsupported operation {}.  Supported operations include {}.'.format(operation, ', '.join(possibleOperations)))
                operations.append(operationKeys[operation])
                genotypeFields.append(field)
                fieldCalcs.append(None)
                if index + 1 >= len(other_stats) or other_stats[index + 1].startswith('--'):
                    raise ValueError('Missing or invalid field name following parameter {}'.format(other_stats[index]))
                destinations.append(other_stats[index + 1])
            else:
                raise ValueError('Expected to see an argument (e.g., --mean_FIELD) here, but found {} instead.'.format(other_stats[index]))
          
        #
        cur = self.db.cursor()
        if IDs is None:
            cur.execute('SELECT sample_id from sample;')
            IDs = [x[0] for x in cur.fetchall()]
        #
        numSample = len(IDs)
        if numSample == 0:
            return
        
        # Error checking with the user specified genotype fields
        # 1) if a field does not exist within one of the sample genotype tables a warning is issued
        # 2) if a field does not exist in any sample, it is not included in validGenotypeFields
        # 3) if no fields are valid and no core stats were requested (i.e., num, het, hom, other), then sample_stat is exited
        genotypeFieldTypes = {}
        fieldInTable = defaultdict(list)
        for id in IDs:
            fields = [x.lower() for x in self.proj.db.getHeaders('{}_genotype.genotype_{}'.format(self.proj.name, id))]
            for field in fields:
                fieldInTable[field].append(id)
                if field not in genotypeFieldTypes:
                    genotypeFieldTypes[field] = 'INT'
                    fieldType = self.db.typeOfColumn('{}_genotype.genotype_{}'.format(self.proj.name, id), field) 
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
                self.logger.warning("Field {} is not an existing genotype field within your samples: {}".format(field, str(genotypeFieldTypes.keys())))
            else:
                if len(fieldInTable[field.lower()]) < len(IDs):
                    self.logger.warning('Field {} exists in {} of {} selected samples'.format(field, len(fieldInTable[field.lower()]), len(IDs))) 
                validGenotypeIndices.append(index)
                validGenotypeFields.append(field)

        if num is None and het is None and hom is None and other is None and len(validGenotypeFields) == 0:
            self.logger.warning("No valid sample statistics operation has been specified.")
            return
        
        queryDestinations = coreDestinations
        for index in validGenotypeIndices:
            queryDestinations.append(destinations[index])
        for name in queryDestinations:
            if name is not None:
                self.proj.checkFieldName(name, exclude=variant_table)
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
            sum([self.db.numOfRows('{}_genotype.genotype_{}'.format(self.proj.name, id)) for id in IDs]))
        count = 0
        for id in IDs:
            whereClause = ''
            if genotypes is not None and len(genotypes) != 0:
                whereClause = 'where ' + ' AND '.join(['({})'.format(x) for x in genotypes])
            
            fieldSelect = ''
            if validGenotypeFields is not None and len(validGenotypeFields) != 0:
                fieldSelect = ' '.join([', ' + x if id in fieldInTable[x.lower()] else ', NULL' for x in validGenotypeFields])
            
            query = 'SELECT variant_id, GT{} FROM {}_genotype.genotype_{} {};'.format(fieldSelect,
                self.proj.name, id, whereClause)
            cur.execute(query)

            for rec in cur:
                if len(from_variants) == 0 or rec[0] in from_variants:
                    if rec[0] not in variants:
                        variants[rec[0]] = [0, 0, 0, 0]
                        variants[rec[0]].extend(list(fieldCalcs))

                    # type heterozygote
                    if rec[1] == 1:
                        variants[rec[0]][0] += 1
                    # type homozygote
                    elif rec[1] == 2:
                        variants[rec[0]][1] += 1
                    # type double heterozygote with two different alternative alleles
                    elif rec[1] == -1:
                        variants[rec[0]][2] += 1
                    elif rec[1] in [0, None]:
                        pass
                    else:
                        self.logger.warning('Invalid genotype type {}'.format(rec[1]))
                
                    # this collects genotype_field information
                    if len(validGenotypeFields) > 0:
                        for index in validGenotypeIndices:
                            queryIndex = index + 2     # to move beyond the variant_id and GT fields in the select statement
                            recIndex = index + 4       # first 4 attributes of variants are het, hom, double_het and wildtype
                            # ignore missing (NULL) values
                            if rec[queryIndex] is None:
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
                count += 1
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        if len(variants) == 0:
            raise ValueError('No variant is updated')
        #
        headers = [x.lower() for x in self.db.getHeaders(variant_table)]
        table_attributes = [(num, 'INT'), (hom, 'INT'),
                (het, 'INT'), (other, 'INT')]
        fieldsDefaultZero = [num, hom, het, other]
        
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
            # We are setting default values on the count fields to 0.  The genotype stat fields are set to NULL by default.
            defaultValue = 0 if field in fieldsDefaultZero else None
            if field.lower() in headers:
                if self.proj.db.typeOfColumn(variant_table, field) != (fldtype + ' NULL'):
                    self.logger.warning('Type mismatch (existing: {}, new: {}) for column {}. Please remove this column and recalculate statistics if needed.'\
                        .format(self.proj.db.typeOfColumn(variant_table, field), fldtype, field))
                self.logger.info('Resetting values at existing field {}'.format(field))
                self.db.execute('Update {} SET {} = {};'.format(variant_table, field, self.db.PH), (defaultValue, ))
            else:
                self.logger.info('Adding field {}'.format(field))
                self.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(variant_table, field, fldtype))
                if defaultValue == 0:
                    self.db.execute ('UPDATE {} SET {} = 0'.format(variant_table, field))              
        #
        prog = ProgressBar('Updating {}'.format(variant_table), len(variants))
        update_query = 'UPDATE {0} SET {2} WHERE variant_id={1};'.format(variant_table, self.db.PH,
            ' ,'.join(['{}={}'.format(x, self.db.PH) for x in queryDestinations if x is not None]))
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
                
            # for genotype_field operations, the value[operation_index] holds the result of the operation
            # except for the "mean" operation which needs to be divided by num_samples that have that variant
            try:
                for index in validGenotypeIndices:
                    operationIndex = index + 4     # the first 3 indices hold the values for hom, het, double het and wildtype
                    operationCalculation = value[operationIndex]
                    if operations[index] == MEAN and operationCalculation is not None:
                        res.append(float(operationCalculation[0]) / operationCalculation[1])
                    else:
                        res.append(operationCalculation)
                cur.execute(update_query, res + [id])
            except Exception as e:
                self.logger.debug(e)
            if count % self.db.batch == 0:
                self.db.commit()
                prog.update(count)
        self.db.commit()
        prog.done()
        self.logger.info('{} records are updated'.format(count))
                
def phenotypeArguments(parser):
    '''Action that can be performed by this script'''
    parser.add_argument('-f', '--from_file', nargs='*',
        help='''Import phenotype from a tab delimited file. The file should have
            a header, with either 'sample_name' as the first column, or 'filename'
            and 'sample_name' as the first two columns. In the former case, samples
            with the same 'sample_name' will share the imported phenotypes. If 
            a list of phenotypes (columns of the file) is specified after filename,
            only the specified phenotypes will be imported. Parameter --samples
            could be used to limit the samples for which phenotypes are imported.'''),
    parser.add_argument('--set', nargs='*', default=[],
        help='''Set a phenotype to a constant (e.g. --set aff=1), or an expression
            using other existing phenotypes (e.g. --set ratio_qt=high_qt/all_qt (the ratio
            of the number of high quality variants to the number of all variants, where
            high_qt and all_qt are obtained from sample statistics using parameter
            --from_stat). Parameters --samples could be used to limit the samples for
            which genotypes will be set.'''),
    parser.add_argument('--from_stat', nargs='*', default=[],
        help='''Set a phenotype to a summary statistics of a genotype field. For 
            example, '--stat "num=count(*)"' sets phenotype num to be the number of
            genotypes of a sample, '--set "DP=avg(DP)"' sets phenotype DP to be the 
            average depth (if DP is one of the genotype fields) of the sample. Multiple
            fields (e.g. '--set "num=count(*)" "DP=avg(DP)"') are also allowed. 
            Parameters --genotypes and --samples could be used to limit the genotypes
            to be considered and the samples for which genotypes will be set.'''),
    parser.add_argument('-g', '--genotypes', nargs='*', default=[],
        help='''Limit the operation to genotypes that match specified conditions.
            Use 'vtools show genotypes' to list usable fields for each sample.'''),
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Update phenotype for samples that match specified conditions.
            Use 'vtools show samples' to list usable fields in the sample table.''')

def phenotype(args):
    try:
        with Project() as proj:
            proj.db.attach('{}_genotype'.format(proj.name))
            p = Sample(proj)
            if args.from_file:
                filename = args.from_file[0]
                fields = args.from_file[1:]
                p.load(filename, fields, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.set:
                for item in args.set:
                    field, expr = [x.strip() for x in item.split('=', 1)]
                    if not expr:
                        raise ValueError('Invalid parameter {}, which should have format field=expr_of_phenotype'.format(item))
                    p.setPhenotype(field, expr, ' AND '.join(['({})'.format(x) for x in args.samples]))
            if args.from_stat:
                stat = []
                for item in args.from_stat:
                    field, expr = [x.strip() for x in item.split('=', 1)]
                    if not expr:
                        raise ValueError('Invalid parameter {}, which should have format field=expr_of_field'.format(item))
                    stat.append((field, expr))
                p.fromSampleStat(stat,
                        ' AND '.join(['({})'.format(x) for x in args.genotypes]),
                        ' AND '.join(['({})'.format(x) for x in args.samples]))
        proj.close()
    except Exception as e:
        sys.exit(e)
                

def sampleStatArguments(parser):
    '''Arguments to calculate sample statistics such as allele count'''
    parser.add_argument('-s', '--samples', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show sample' (e.g. 'aff=1',
            'filename like "MG%%"').''')
    parser.add_argument('--genotypes', nargs='*', default=[],
        help='''Limiting variants from samples that match conditions that
            use columns shown in command 'vtools show genotypes' (e.g. 'GQ_INFO>15').''')
    parser.add_argument('table',
        help='''Variant table for which the statistics will be calculated and updated.''')
    parser.add_argument('-n', '--num',
        help='''Name of the field to hold number of alternative alleles in the sample.''')
    parser.add_argument('--hom',
        help='''Name of the field to hold number of samples with two identical alternative alleles.''')
    parser.add_argument('--het',
        help='''Name of the field to hold number of samples with one reference and one alternative alleles.''')
    parser.add_argument('--other',
        help='''Name of the field to hold number of samples with two different alternative alleles.''')
    
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
                IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.samples]))
                if len(IDs) == 0:
                    p.logger.info('No sample is selected (or available)')
                    return
                else:
                    p.logger.info('{} samples are selected'.format(len(IDs)))
            p.calcSampleStat(IDs, variant_table, args.genotypes, args.num, args.hom,
                args.het, args.other, args.unknown_args)
        # temporary tables will be removed
        proj.close()
    except Exception as e:
        sys.exit(e)

