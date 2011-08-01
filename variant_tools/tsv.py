#!/usr/bin/env python
#
# $File: tsv.py $
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
from .project import Project
from .utils import ProgressBar, lineCount, getMaxUcscBin

class tsvImporter:
    '''Import variants from one or more tab or comma separated files.'''
    def __init__(self, proj, files, col, table, build, delimiter, zero):
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        #
        if len(files) == 0:
            raise IOError('Please specify the filename of the input data.')
            sys.exit(1)
        #
        if build is None:
            if self.proj.build is None:
                raise ValueError('Please specify the reference genome of the input data.')
        else:
            if self.proj.build is None:
                self.proj.setRefGenome(build)
            elif build != self.proj.build:
                raise ValueError('Specified build {} must be the primary build of the project'.format(build))
        #
        self.table = table
        if self.table == 'variant':
            raise ValueError('Value for parameter --table should not be the master variant table')
        #
        self.files = []
        cur = self.db.cursor()
        cur.execute('SELECT filename from filename;')
        existing_files = [x[0] for x in cur.fetchall()]
        for f in files:
            filename = os.path.split(f)[-1]
            if self.table is None and filename in existing_files:
                self.logger.info('Ignoring imported file {}'.format(filename))
            else:
                self.files.append(f)
        if len(self.files) == 0:
            return
        #
        self.col = [x - 1 for x in col]
        if len(self.col) != 4:
            raise ValueError('Four columns are required for each variant (chr, pos, ref, and alt)')
        self.delimiter = delimiter
        self.zero = zero
        #
        self.createLocalVariantIndex()
        self.proj.dropIndexOnMasterVariantTable()

    def __del__(self):
        self.proj.createIndexOnMasterVariantTable()

    def openTSV(self, filename):
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
            # 0 means NOT new variant, records with the 1-flag will be written to self.table
            # if applicable.
            self.variantIndex[(rec[1], rec[2], rec[3])] = [rec[0], 0]
            if count % self.db.batch == 0:
                prog.update(count)
        prog.done()

    def importTSV(self, input_filename):
        '''Import a TSV file to sample_variant'''
        #
        # record filename after getMeta because getMeta might fail (e.g. cannot recognize reference genome)
        cur = self.db.cursor()
        filename = os.path.split(input_filename)[-1]
        try:
            cur.execute("INSERT INTO filename (filename) VALUES ({});".format(self.db.PH), (filename,))
        except Exception as e:
            # filename might already exist, does not matter in this case.
            self.logger.debug(e)
        #
        all_records = 0
        skipped_records = 0
        inserted_variants = 0
        #
        variant_insert_query = 'INSERT INTO variant (bin, chr, pos, ref, alt) VALUES ({0}, {0}, {0}, {0}, {0});'.format(self.db.PH)
        prog = ProgressBar(os.path.split(input_filename)[-1], lineCount(input_filename))
        with self.openTSV(input_filename) as input_file:
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
                    try:
                        variant_info = self.variantIndex[(chr, pos, alt)]
                        # already exist
                        variant_info[1] = 1
                    except:
                        bin = getMaxUcscBin(pos - 1, pos)
                        cur.execute(variant_insert_query, (bin, chr, pos, ref, alt))
                        variant_id = cur.lastrowid
                        self.variantIndex[(chr, pos, alt)] = [variant_id, 1]
                        inserted_variants += 1
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

    def importData(self):
        '''Start importing'''
        # 
        imported = 0
        for count,f in enumerate(self.files):
            self.logger.info('Importing genotype from {} ({}/{})'.format(f, count + 1, len(self.files)))
            imported += self.importTSV(f)
        self.logger.info('All files imported. A total of {0:,} new records are inserted.'.format(imported))
        #
        if not self.table:
            return imported
        # write to table
        if self.db.hasTable(self.table):
            new_table = self.db.backupTable(self.table)
            self.logger.warning('Existing table {} is renamed to {}.'.format(self.table, new_table))
        self.proj.createVariantTable(self.table)
        prog = ProgressBar('Creating ' + self.table, imported)
        query = 'INSERT INTO {} VALUES ({});'.format(self.table, self.db.PH)
        # get variants with flag 1
        var = [x for x,y in self.variantIndex.values() if y == 1]
        cur = self.db.cursor()
        for count,id in enumerate(sorted(var)):
            cur.execute(query, (id,))
            if count % self.db.batch == 0:
                self.db.commit()
                prog.update(count)
        self.db.commit()
        prog.done()
        return imported

#
#
# Functions provided by this script
#
#

def importTSVArguments(parser):
    parser.add_argument('input_files', nargs='*',
        help='''A list of files that will be imported. The file should be in 
            tab or command separated value format. Gzipped files are acceptable.''')
    parser.add_argument('-t', '--table',
        help='''If specified, a variant table will be created with variants from
            the input files. This option is usually used to create a table of variants
            from a manually selected list of variants of a project.''')
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

def importTSV(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            importer = tsvImporter(proj=proj, files=args.input_files,
                col=args.columns, table=args.table, build=args.build,
                delimiter=args.delimiter, zero=args.zero)
            importer.importData()
        proj.close()
    except Exception as e:
        sys.exit(e)

