#!/usr/bin/env python
#
# $File: importer.py $
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

import os
import sys
import gzip
import re
from .project import Project
from .utils import ProgressBar, lineCount, getMaxUcscBin

class Importer:
    '''A general class for importing variants'''
    def __init__(self, proj, files, build):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
        #
        if build:
            self.proj.setRefGenome(build)
        #
        self.files = []
        cur = self.db.cursor()
        cur.execute('SELECT filename from filename;')
        existing_files = [x[0] for x in cur.fetchall()]
        for f in files:
            filename = os.path.split(f)[-1]
            if filename in existing_files:
                self.logger.info('Ignoring imported file {}'.format(filename))
            else:
                self.files.append(f)
        if len(self.files) == 0:
            return
        self.proj.dropIndexOnMasterVariantTable()
        #
        self.createLocalVariantIndex()

    def __del__(self):
        self.proj.createIndexOnMasterVariantTable()

    def openFile(self, filename):
        if filename.lower().endswith('.gz'):
            return gzip.open(filename, 'rb')
        else:
            # text file
            return open(filename, 'r')

    def createLocalVariantIndex(self):
        '''Create index on variant (chr, pos, alt) -> variant_id'''
        self.variantIndex = {}
        cur = self.db.cursor()
        numVariants = self.db.numOfRows('variant')
        if numVariants == 0:
            return
        self.logger.debug('Creating local indexes for {:,} variants'.format(numVariants));
        cur.execute('SELECT variant_id, chr, pos, alt FROM variant;')
        prog = ProgressBar('Getting existing variants', numVariants)
        for count, rec in enumerate(cur):
            self.variantIndex[(rec[1], rec[2], rec[3])] = rec[0]
            if count % self.db.batch == 0:
                prog.update(count)
        prog.done()

    def recordFileAndSample(self, filename, sampleNames, sampleFields = []):
        cur = self.db.cursor()
        cur.execute("INSERT INTO filename (filename) VALUES ({});".format(self.db.PH), (filename,))
        filenameID = cur.lastrowid
        self.proj.createSampleTableIfNeeded()
        sample_ids = []
        for samplename in sampleNames:
            cur.execute('INSERT INTO sample (file_id, sample_name) VALUES ({0}, {0});'.format(self.db.PH),
                (filenameID, samplename))
            sample_ids.append(cur.lastrowid)
            self.proj.createNewSampleVariantTable('{0}_genotype.sample_variant_{1}'.format(self.proj.name, cur.lastrowid),
                sampleFields)
        return sample_ids
        
    def importData(self):
        '''Start importing'''
        # 
        inserted = 0
        for count,f in enumerate(self.files):
            self.logger.info('Importing genotype from {} ({}/{})'.format(f, count + 1, len(self.files)))
            inserted += self.importFromFile(f)
        self.logger.info('All files imported. A total of {0:,} new records are inserted.'.format(inserted))
        return inserted


class vcfImporter(Importer):
    '''A vcf importer to import genotype from one or more vcf files.
    In case of vcf file, it records the type of variant for this sample, 
    which can be 1 for genotype 0/1, heterozygous, 2 for genotype 1/1,
    homozygous of alternative alleles, and -1 for gentoype 1/2, which 
    consists of two different alternative alleles. In the last case,
    an individual will have two variants with different alternative
    alleles, each with a type -1. That it to say, allele frequency should
    be calculated as sum (abs(type))/ (2*num_of_sample). '''
    def __init__(self, proj, files, build=None, variant_only=False, info=[]):
        '''see importVariant.py -h for details about parameters. Additional
        keyword paramters such as user, passwd and host are passed to
        MySQLdb.connect.
        '''
        Importer.__init__(self, proj, files, build)
        # vcf tools only support DP for now
        self.variant_only = variant_only
        self.import_depth = 'DP' in info
        # FIXME: self.infoFields and formatFields should initially read from
        # table variant_meta if this table already exists.        # 
        self.infoFields = None  # will be assigned when the first vcf file is read
        self.formatFields = None

    def getMetaInfo(self, filename):
        '''Probe vcf files for additional information to be put to the variant_meta table.
        '''
        infoFields = []
        formatFields = []
        samples = []
        build = None
        with self.openFile(filename) as input:
            line = input.readline()
            if not line.startswith('##fileformat=VCF'):
                self.logger.error('Invalid vcf file. file not started with line ##fileformat')
                raise ValueError('Invalid vcf file')
            if not line.strip().endswith('4.0') and not line.strip().endswith('4.1'):
                raise ValueError('This importer tool only supports VCF format v4.0 and v4.1. \
                    Please use vcftools to convert your vcf file to a supported format')
            #
            for line in input:
                if line.startswith('##reference'):
                    # guess reference genome from VCF header file
                    if 'NCBI36' in line.upper() or 'HG18' in line.upper() or 'HUMAN_B36' in line.upper():
                        build = 'hg18'
                    elif 'GRCH37' in line.upper() or 'HG19' in line.upper() or 'HUMAN_B37' in line.upper():
                        build = 'hg19'
                if line.startswith('INFO'):
                    # FIXME: properly handle INFO
                    pass
                if line.startswith('FORMAT'):
                    # FIXME: properly handle FORMAT
                    pass
                if line.startswith('#CHR'):
                    samples = line.split()[9:]
                if not line.startswith('#'):
                    break
        # set meta fields if this is the first file
        if self.infoFields is None:
            self.infoFields = infoFields
        elif self.infoFields != infoFields:
            # FIXME: give a warning?
            pass
        if self.formatFields is None:
            self.formatFields = formatFields
        elif self.formatFields != formatFields:
            # FIXME: give a warning?
            pass
        if self.proj.build is None:
            if build is None:
                raise ValueError('Cannot determine reference genome build from the meta information vcf file\n'
                        'Please use parameter --build to specify it.')
            else:
                self.logger.info('Reference genome build is determined to be {0}'.format(build))
                self.proj.setRefGenome(build)
        return samples

    def importFromFile(self, input_filename):
        '''Import a VCF file to sample_variant'''
        #
        # handle meta information and get sample names
        sampleNames = self.getMetaInfo(input_filename)
        #
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        no_sample = self.variant_only or len(sampleNames) == 0
        sample_ids = self.recordFileAndSample(os.path.split(input_filename)[-1], [None] if no_sample else sampleNames, 
            ['DP'] if self.import_depth else [])   # record individual depth, total depth is divided by number of sample in a file
        #
        all_records = 0
        skipped_records = 0
        inserted_variants = 0
        nSample = len(sample_ids)
        #
        DP_pattern = re.compile('.*DP=(\d+)')
        #
        cur = self.db.cursor()
        variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        # sample variants are inserted into different tables in a separate database.
        sample_variant_insert_query = {x: 'INSERT INTO {1}_genotype.sample_variant_{3} VALUES ({0}, {0} {2});'\
            .format(self.db.PH, self.proj.name, ',' + self.db.PH if self.import_depth else '', x) for x in sample_ids}
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                try:
                    # FIXME: # record sample meta information
                    if line.startswith('#'):
                        continue
                    # get data
                    tokens = [x.strip() for x in line.split('\t')]
                    chr = tokens[0][3:] if tokens[0].startswith('chr') else tokens[0]
                    pos = int(tokens[1])
                    ref = tokens[3]
                    if len(ref) > 1:
                        # this is a deletion
                        raise ValueError('vtools currently does not handle deletion')
                    # FIXME: handle INFO and FORMAT, here we only extract info
                    # get depth.
                    if self.import_depth:
                        m = DP_pattern.match(tokens[7])
                        DP = [None if m is None else float(m.group(1))/nSample]
                    else:
                        DP = []
                    if len(tokens[4]) == 1:
                        all_records += 1
                        # the easy case: there is only one alternative allele,
                        # all genotypes should be 0/0, 0/1, 1/0, or 1/1.
                        alt = tokens[4][0]
                        try:
                            variant_id = self.variantIndex[(chr, pos, alt)]
                        except:
                            bin = getMaxUcscBin(pos - 1, pos)
                            cur.execute(variant_insert_query, (bin, chr, pos, ref, alt))
                            variant_id = cur.lastrowid
                            self.variantIndex[(chr, pos, alt)] = variant_id
                            inserted_variants += 1
                        #
                        # FIXME: we should properly handle self.formatFields
                        if no_sample:
                            cur.execute(sample_variant_insert_query[sample_ids[0]], [variant_id, 1] + DP)
                        else:
                            variants = [x.split(':')[0].count('1') for x in tokens[-len(sample_ids):]]
                            for var_idx, var in enumerate(variants):
                                if var != 0:  # genotype 0|0 are ignored
                                    cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id, var] + DP)
                    else:
                        # there are two alternative alleles
                        variant_id = [0, 0]
                        for altidx, alt in enumerate(tokens[4].split(',')):
                            if len(alt) > 1:
                                raise ValueError('vtools currently does not handle insertion')
                            all_records += 1
                            try:
                                variant_id[altidx] = self.variantIndex[(chr, pos, alt)]
                            except:
                                bin = getMaxUcscBin(pos - 1, pos)
                                cur.execute(variant_insert_query, (bin, chr, pos, ref, alt))
                                variant_id[altidx] = cur.lastrowid
                                self.variantIndex[(chr, pos, alt)] = variant_id[altidx]
                                inserted_variants += 1
                        if no_sample:
                            cur.execute(sample_variant_insert_query[sample_ids[0]], [variant_id[0], 1] + DP)
                            cur.execute(sample_variant_insert_query[sample_ids[0]], [variant_id[1], 1] + DP)
                        else:
                            # process variants
                            for var_idx, var in enumerate([x.split(':')[0] for x in tokens[-len(sample_ids):]]):
                                if len(var) == 3:  # regular
                                    gt = var[0] + var[2]  # GT can be separated by / or |
                                    if gt in ['01', '10']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif gt in ['02', '20']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
                                    elif gt == '11':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 2] + DP)
                                    elif gt in ['12', '21']:
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], -1] + DP)
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], -1] + DP)
                                    elif gt == '22':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 2] + DP)
                                    elif gt == '00':
                                        pass
                                    else:
                                        raise ValueError('I do not know how to process genotype {}'.format(var))
                                else: # should have length 1
                                    if var == '1':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[0], 1] + DP)
                                    elif var == '2':
                                        cur.execute(sample_variant_insert_query[sample_ids[var_idx]], [variant_id[1], 1] + DP)
                                    elif var == '0':
                                        pass
                                    else:
                                        raise ValueError('I do not know how to process genotype {}'.format(var))
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                    skipped_records += 1
                if all_records % self.db.batch == 0:
                    self.db.commit()
                    prog.update(all_records)
            self.db.commit()
            prog.done()
        self.logger.info('{:,} new variants from {:,} records are imported, with {:,} invalid records.'\
            .format(inserted_variants, all_records, skipped_records))
        return inserted_variants


class txtImporter(Importer):
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, files, col, build, delimiter, zero):
        Importer.__init__(self, proj, files, build)
        # we cannot guess build information from txt files
        if build is None and self.proj.build is None:
            raise ValueError('Please specify the reference genome of the input data.')
        #
        self.col = [x - 1 for x in col]
        if len(self.col) != 4:
            raise ValueError('Four columns are required for each variant (chr, pos, ref, and alt)')
        self.delimiter = delimiter
        self.zero = zero

    def importFromFile(self, input_filename):
        '''Import a TSV file to sample_variant'''
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        filename = os.path.split(input_filename)[-1]
        # assuming one sample for each file
        sample_id = self.recordFileAndSample(filename, [None])[0]
        #
        all_records = 0
        skipped_records = 0
        inserted_variants = 0
        #
        cur = self.db.cursor()
        variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        sample_variant_insert_query = 'INSERT INTO {0}_genotype.sample_variant_{1} VALUES ({2}, {2});'\
            .format(self.proj.name, sample_id, self.db.PH)
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openFile(input_filename) as input_file:
            for line in input_file:
                all_records += 1
                try:
                    # get data
                    tokens = [x.strip() for x in line.split(self.delimiter)]
                    chr, pos, ref, alt = [tokens[x] for x in self.col]
                    if chr.startswith('chr'):
                        chr = chr[3:]
                    pos = int(pos) + 1 if self.zero else int(pos)
                    if len(ref) != 1:
                        raise ValueError('Incorrect reference allele: {}'.format(ref))
                    if len(alt) != 1:
                        raise ValueError('Incorrect alternative allele: {}'.format(alt))
                    # variant
                    try:
                        variant_id = self.variantIndex[(chr, pos, alt)]
                    except: # new variant
                        bin = getMaxUcscBin(pos - 1, pos)
                        cur.execute(variant_insert_query, (bin, chr, pos, ref, alt))
                        variant_id = cur.lastrowid
                        self.variantIndex[(chr, pos, alt)] = variant_id
                        inserted_variants += 1
                    # sample variant, the variant type is always hetero???
                    cur.execute(sample_variant_insert_query, (variant_id, 1))
                except Exception as e:
                    self.logger.debug('Failed to process line: ' + line.strip())
                    self.logger.debug(e)
                    skipped_records += 1
                if all_records % self.db.batch == 0:
                    self.db.commit()
                    prog.update(all_records)
            self.db.commit()
            prog.done()
        self.logger.info('{:,} new variants from {:,} records are imported, with {:,} invalid records.'\
            .format(inserted_variants, all_records, skipped_records))                
        if all_records == skipped_records:
          self.logger.warning('No valid record is imported')
          cur.execute("DELETE FROM filename WHERE (filename) = ({});".format(self.db.PH), (filename,))
        return inserted_variants

#
#
# Functions provided by this script
#
#

def importVCFArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''A list of files that will be imported. The file should be in 
            VCF 4.0 format and can be compressed in gzip format.''')
    parser.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18). This should be the reference
            genome of the input data and can only be the primary reference genome of the project.''')
    parser.add_argument('--variant_only', action='store_true',
        help='''Import only variants, and ignore sample and their genotypes.''')
    parser.add_argument('--info', nargs='*', default=['DP'],
        help='''Variant information fields to import. This command only support
            'DP' (total depth). When 'DP' is listed (default), vtools will look
            for total depth (DP=) in the INFO field of each variant and set average
            depth to each individual (DP/numSample in the vcf file) in field 'DP'.''')


def importVCF(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = vcfImporter(proj=proj, files=args.input_files, build=args.build,
                variant_only=args.variant_only, info=[] if args.variant_only else args.info)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)


def exportVCFArguments(parser):
    pass

def exportVCF(args):
    raise SystemError('This feature is currently not implemented')


def importTxtArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''A list of files that will be imported. The file should be in 
            tab or command separated value format. Gzipped files are acceptable.''')
    grp = parser.add_argument_group('Description of input files')
    grp.add_argument('--build',
        help='''Build version of the reference genome (e.g. hg18). This should be the
            reference genome of the input data and has to be the primary reference
            genome of the project. If unspecified, it assumes to be the primary
            reference genome of the exisiting project.''')
    grp.add_argument('-c', '--columns', default=[1,2,3,4], nargs='+', type=int,
        help='Columns for chromosome, position, reference and alternative alleles.')
    grp.add_argument('-d', '--delimiter', default='\t',
        help='''Delimiter, default to tab, a popular alternative is ',' for csv output''')
    grp.add_argument('-z', '--zero', action='store_true',
        help='''Whether or not specified file uses zero-based index. If unspecified, the
            position column is assumed to be 1-based.''')

def importTxt(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach(proj.name + '_genotype')
            importer = txtImporter(proj=proj, files=args.input_files,
                col=args.columns, build=args.build,
                delimiter=args.delimiter, zero=args.zero)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

