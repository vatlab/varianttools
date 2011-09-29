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
            for idx, line in enumerate(input.readlines()):
                if line.startswith('#') or line.strip() == '':
                    continue
                fields = line.rstrip().split('\t')
                if len(fields) != nFields or ('' in fields):
                    raise ValueError('Invalid phenotype file: number of fields mismatch at line {}. Please use \'None\' for missing values.'.format(idx+2))
                #
                record[(fields[0], None if fields[1] == 'None' else fields[1])] = fields[2:]
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
                self.logger.warning('No phenotype is defined for sample: filename={}\tsample_name={}'.format(rec[1], rec[3]))
                table.append([rec[0], rec[2], rec[3]] + [None]*(nFields-2))
            else:
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
        operationKeys = {}
        operationKeys['mean'] = MEAN
        operationKeys['sum'] = SUM
        operationKeys['min'] = MIN
        operationKeys['max'] = MAX
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
        for id in IDs:
            fields = self.proj.db.getHeaders('{}_genotype.sample_variant_{}'.format(self.proj.name, id))
            for field in fields:
                if genotypeFieldTypes.get(field) is None:
                    genotypeFieldTypes[field] = 'INT'
                    fieldType = self.db.typeOfColumn('{}_genotype.sample_variant_{}'.format(self.proj.name, id), field) 
                    if fieldType.upper().startswith('FLOAT'):
                        genotypeFieldTypes[field] = 'FLOAT'
                    elif fieldType.upper().startswith('VARCHAR'):
                        raise ValueError('Genotype field {} is a VARCHAR which is not supported with sample_stat operations.'.format(field))
            
        validGenotypeIndices = []
        for index in range(0,len(genotypeFields)):
            field = genotypeFields[index]
            if field not in genotypeFieldTypes.keys():
                self.logger.warning("Field {} is not an existing genotype field within your samples: {}".format(field, str(genotypeFieldTypes.keys())))
            else:
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
            sum([self.db.numOfRows('{}_genotype.sample_variant_{}'.format(self.proj.name, id)) for id in IDs]))
        count = 0
        for id in IDs:
            whereClause = ''
            if genotypes is not None and len(genotypes) != 0:
                whereClause = 'where ' + ' AND '.join(genotypes)
            
            fieldSelect = ''
            if validGenotypeFields is not None and len(validGenotypeFields) != 0:
                fieldSelect = ', ' + ', '.join(validGenotypeFields)
            
            try:    
                query = 'SELECT variant_id, variant_type{} FROM {}_genotype.sample_variant_{} {};'.format(fieldSelect, self.proj.name, id, whereClause)
                cur.execute(query)
            except:
                self.logger.warning('Sample {} does not have all the requested genotype fields [{}].'.format(id, ', '.join(set(validGenotypeFields))))
                
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
                    # type wildtype
                    elif rec[1] == 0:
                        variants[rec[0]][3] += 1
                    elif rec[1] is None:
                        pass
                    else:
                        self.logger.warning('Invalid genotype type {}'.format(rec[1]))
                
                    # this collects genotype_field information
                    if len(validGenotypeFields) > 0:
                        for index in validGenotypeIndices:
                            queryIndex = index + 2     # to move beyond the variant_id and variant_type fields in the select statement
                            recIndex = index + 4       # first 4 attributes of variants are het, hom, double_het and wildtype
                            operation = operations[index]
                            field = genotypeFields[index]
                            if operation in [MEAN, SUM]:
                                if variants[rec[0]][recIndex] is None:
                                    variants[rec[0]][recIndex] = rec[queryIndex]
                                else:
                                    variants[rec[0]][recIndex] += rec[queryIndex]
                            if operation == MIN:
                                if variants[rec[0]][recIndex] is None or rec[queryIndex] < variants[rec[0]][recIndex]:
                                    variants[rec[0]][recIndex] = rec[queryIndex]
                            if operation == MAX:
                                if variants[rec[0]][recIndex] is None or rec[queryIndex] > variants[rec[0]][recIndex]:
                                    variants[rec[0]][recIndex] = rec[queryIndex]  
                    
                count += 1
                if count % self.db.batch == 0:
                    prog.update(count)
        prog.done()
        #
        headers = self.db.getHeaders(variant_table)
        [num, hom, het, other]
        table_attributes = [(num, 'INT'), (hom, 'INT'),
                (het, 'INT'), (other, 'INT')]
        for index in validGenotypeIndices:
            if operations[index] == MEAN:
                table_attributes.append((destinations[index], 'FLOAT'))
            else:
                table_attributes.append((destinations[index], genotypeFieldTypes.get(genotypeFields[index])))
        for field, fldtype in table_attributes:
            if field is None:
                continue
            if field in headers:
                self.logger.info('Updating existing field {}'.format(field))
                if fldtype == 'FLOAT' and operations[index] != MEAN:
                    self.logger.warning('Result will be wrong if field \'{}\' was created to hold integer values'.format(field))
            else:
                self.logger.info('Adding field {}'.format(field))
                self.db.execute('ALTER TABLE {} ADD {} {} NULL;'.format(variant_table, field, fldtype))               
        #

        prog = ProgressBar('Updating table {}'.format(variant_table), len(variants))
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
            for index in validGenotypeIndices:
                operationIndex = index + 4     # the first 3 indices hold the values for hom, het, double het and wildtype
                operationCalculation = value[operationIndex]
                if operations[index] == MEAN:
                    numSamples = value[0] + value[1] + value[2] + value[3]
                    operationCalculation /= float(numSamples)
                res.append(operationCalculation)
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
                IDs = proj.selectSampleByPhenotype(' AND '.join(args.samples))
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

