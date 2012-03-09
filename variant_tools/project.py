#!/usr/bin/env python
#
# $File: project.py $
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

import os
import sys
import glob
import logging
import getpass
import random
import textwrap
import tempfile
import shutil
import ConfigParser
import argparse
import threading
import Queue
import time
import re
from subprocess import Popen, PIPE
from collections import namedtuple, defaultdict
from .__init__ import VTOOLS_VERSION, VTOOLS_FULL_VERSION, VTOOLS_COPYRIGHT, VTOOLS_CITATION, VTOOLS_CONTACT
from .utils import DatabaseEngine, ProgressBar, setOptions, SQL_KEYWORDS, delayedAction, \
    filesInURL, downloadFile, makeTableName, getMaxUcscBin


# define a field type
Field = namedtuple('Field', ['name', 'index', 'adj', 'type', 'comment'])
Column = namedtuple('Column', ['index', 'field', 'adj', 'comment'])
#
# How field will be use in a query. For example, for field sift, it is
# connection clause will be:
#   field = dbNSFP.sift
#   table = dbNSFP.dbNSFP  # annotation
#   link = chr=dbNSFP.chr AND pos=dbNSFP.h18pos
#
FieldConnection = namedtuple('FieldConnection', ['field', 'table', 'link'])


class AnnoDB:
    '''A structure that is created from an existing annotation database.
    The annotation.py module is responsible of creating this structure from
    various sources.
    '''
    def __init__(self, proj, annoDB, linked_by=[], anno_type=None, linked_fields=None):
        proj.logger.debug('Loading annotation database {}'.format(annoDB))
        self.db = proj.db.newConnection()
        if self.db.hasDatabase(annoDB):
            self.db.connect(annoDB)
        else:
            raise ValueError("Cannot locate annotation database {}".format(annoDB))
        #
        self.dir = os.path.split(annoDB)[0]
        self.filename = os.path.splitext(os.path.split(annoDB)[-1])[0]
        if '-' in self.filename:
            self.version = self.filename.split('-', 1)[1]
            self.name = self.filename.split('-')[0]
        else:
            self.name = os.path.splitext(self.filename)[0]
            self.version = None
        for table in [self.name, self.name + '_field', self.name + '_info']:
            if not self.db.hasTable(table):
                raise ValueError('{} is not a valid annotation database. Missing table {}.'.format(annoDB, table))
        # read fields from DB
        self.fields = []
        cur = self.db.cursor()
        cur.execute('SELECT name, field, "", type, comment from {}_field;'.format(self.name))
        for rec in cur:
            self.fields.append(Field(*rec))
            # FIXME: We should enforce comment for all fields.
            #if not self.fields[-1].comment:
            #    proj.logger.warning('No comment for field {} in database {}'.format(self.fields[-1].name, annoDB))
        if len(self.fields) == 0:
            raise ValueError('Annotation database {} does not provide any field.'.format(annoDB))
        #
        self.anno_type = 'variant'
        self.linked_by = []
        for f in linked_by:
            self.linked_by.append(proj.linkFieldToTable(f, 'variant')[-1].field)
        self.description = ''
        self.refGenomes = None
        self.build = None
        self.alt_build = None
        self.version = None
        cur.execute('SELECT * from {}_info;'.format(self.name))
        for rec in cur:
            if rec[0] == 'description':
                self.description = rec[1]
            # stored name and version, if exist, will override name and version obtained from filename.
            elif rec[0] == 'name':
                self.name = rec[1]
            elif rec[0] == 'version':
                self.version = rec[1]
            elif rec[0] == 'anno_type':
                if anno_type is None:
                    self.anno_type = rec[1] if rec[1] != 'attribute' else 'field'
                else:
                    self.anno_type = anno_type
            elif rec[0] == 'build':
                self.refGenomes = eval(rec[1])
                if linked_fields is not None:
                    for key in self.refGenomes.keys():
                        self.refGenomes[key] = linked_fields
                for key in self.refGenomes.keys():
                    # no reference genome is needed
                    if key == '*':
                        self.build = self.refGenomes[key]
                    elif proj.build is None:
                        proj.logger.warning('Project does not have a primary build. Using {} from the annotation database'.format(key))
                        proj.setRefGenome(key)
                        self.build = self.refGenomes[key]
                        break
                    elif key == proj.build:
                        self.build = self.refGenomes[key]
                    elif key == proj.alt_build:
                        self.alt_build = self.refGenomes[key]
        #
        if self.anno_type == 'variant' and ((self.build is not None and len(self.build) != 4) or (self.alt_build is not None and len(self.alt_build) != 4)):
            raise ValueError('There should be four linking fields for variant annotation databases.')
        if self.anno_type == 'position' and ((self.build is not None and len(self.build) != 2) or (self.alt_build is not None and len(self.alt_build) != 2)):
            raise ValueError('There should be two linking fields for positional annotation databases.')
        if self.anno_type == 'range' and ((self.build is not None and len(self.build) != 3) or (self.alt_build is not None and len(self.alt_build) != 3)):
            raise ValueError('There should be three linking fields for range-based annotation databases.')
        if self.description == '':
            proj.logger.warning('No description for annotation database {}'.format(annoDB))
        if self.build is None and self.alt_build is None:
            raise ValueError('No reference genome information for annotation database {}'.format(annoDB))
        if self.anno_type == 'field' and len(self.linked_by) != len(self.build):
            raise RuntimeError('Please specify link fields for attributes {} using parameter --linked_by'.format(','.join(self.build)))
        if self.linked_by:
            s = delayedAction(proj.logger.info, 'Indexing linked field {}'.format(', '.join(self.linked_by)))
            self.indexLinkedField(proj, linked_by)
            del s
        if self.anno_type == 'range':
            s = delayedAction(proj.logger.info, 'Binning ranges')
            if self.build is not None:
                self.binningRanges(proj, proj.build, self.build)
            elif self.alt_build is not None:
                self.binningRanges(proj, proj.alt_build, self.alt_build)


    def binningRanges(self, proj, build, keys):
        cur = self.db.cursor()
        tbl = makeTableName([build] + keys)
        if self.db.hasTable(tbl):
            return
        cur.execute('SELECT rowid, {} FROM {}'.format(','.join(keys), self.name))
        ranges = cur.fetchall()
        cur.execute('CREATE TABLE {} (bin INT, chr VARCHAR(255), range_id INT)'.format(tbl))
        insert_query = 'INSERT INTO {0} VALUES ({1}, {1}, {1});'.format(tbl, self.db.PH)
        for rowid, chr, start, end in ranges:
            if start > end:
                raise ValueError('Start position {} greater than ending position {} in database {}'.format(start, end, self.name))
            sbin = getMaxUcscBin(start-1, start)
            ebin = getMaxUcscBin(end-1, end)
            if sbin > ebin:
                raise SystemError('Start bin greater than end bin...')
            cur.executemany(insert_query, [(bin, chr, rowid) for bin in range(sbin, ebin + 1)])
        self.db.commit()          
        cur.execute('CREATE INDEX {0}_idx ON {0} (bin ASC, chr ASC, range_id ASC);'.format(tbl))
        self.db.commit()

    def indexLinkedField(self, proj, linked_fields):
        '''Create index for fields that are linked by'''
        cur = proj.db.cursor()
        for linked_field in linked_fields:
            if '.' not in linked_field:
                linked_field = proj.linkFieldToTable(linked_field, 'variant')[-1].field
            table, field = linked_field.split('.')
            if proj.isVariantTable(table):
                try:
                    cur.execute('CREATE INDEX IF NOT EXISTS {0}_{1} ON {0} ({1} ASC);'.format(table, field))
                except Exception as e:
                    proj.logger.debug(e)
            else:
                # from an annotation database
                try:
                    # FIXME: this syntax is only valid for sqlite3, we will need to fix it for the mysql engine
                    cur.execute('CREATE INDEX IF NOT EXISTS {0}.{0}_{1} ON {0} ({1} ASC);'.format(table, field))
                except Exception as e:
                    proj.logger.debug(e)
        

    def describe(self, verbose=False):
        '''Describe this annotation database'''
        print('Annotation database {} {}'.format(self.name, '(version {})'.format(self.version) if self.version else ''))
        if self.description is not None:
            print('Description: {}'.format('\n'.join(textwrap.wrap(self.description,
                initial_indent='', subsequent_indent=' '*4))))
        print('Database type: {}'.format(self.anno_type))
        # get linking fields
        if verbose:
            # number of records
            cur = self.db.cursor()
            cur.execute('SELECT value FROM {}_info WHERE name="num_records";'.format(self.name))
            num_records = int(cur.fetchone()[0])
            print('Number of records: {:,}'.format(num_records))
            fields = self.refGenomes.values()[0]
            
            # Get number of unique keys
            cur.execute('SELECT value FROM {}_info WHERE name="distinct_keys";'.format(self.name))
            count = int(cur.fetchone()[0])
            
            #
            if self.anno_type == 'variant':
                print('Number of distinct variants: {:,}'.format(count))
            elif self.anno_type == 'position':
                print('Number of distinct positions: {:,}'.format(count))
            elif self.anno_type == 'range':
                print('Number of distinct ranges: {:,}'.format(count))
            elif self.anno_type == 'field':
                print('Number of distinct entries: {:,}'.format(count))
        #
        for key in self.refGenomes:
            print('Reference genome {}: {}'.format(key, self.refGenomes[key]))
        for field in self.fields:
            if not verbose:
                print('    {:<20}{}'.format(field.name, '\n'.join(textwrap.wrap(
                    field.comment, initial_indent=' ', subsequent_indent=' '*25))))
            else:
                print('\nField:           {}'.format(field.name))
                numeric = False
                if 'chromosome' in field.type.lower():
                    print('Type:            chromosome')
                elif 'position' in field.type.lower():
                    print('Type:            integer')
                    numeric = True
                elif 'int' in field.type.lower():
                    print('Type:            integer')
                    numeric = True
                elif 'float' in field.type.lower():
                    print('Type:            float')
                    numeric = True
                else:
                    print('Type:            string')
                if field.comment:
                    print('Comment: {}'.format('\n'.join(textwrap.wrap(
                        field.comment, initial_indent='        ', subsequent_indent=' '*17))))
                cur.execute('SELECT missing_entries FROM {0}_field WHERE name="{1}";'.format(self.name, field.name))
                missing = cur.fetchone()[0]
                #
                print('Missing entries: {:,} {}'.format(missing, '({:.1f}% of {:,} records)'.format(100. * missing/num_records, num_records) if missing else ''))
                if missing == num_records:
                    continue
                cur.execute('SELECT distinct_entries {2} FROM {0}_field WHERE name={1};'.format(
                    self.name, self.db.PH, ', min_value, max_value' if numeric else ''), (field.name,))
                res = cur.fetchone()
                print('Unique Entries:  {:,}'.format(res[0]))
                if numeric:
                    print('Range:           {} - {}'.format(res[1], res[2]))


class fileFMT:
    def __init__(self, name, fmt_args=[], logger=None):
        '''Input file format'''
        # locate a file format specification file
        self.description = None
        self.variant_fields = None
        self.position_fields = None
        self.range_fields = None
        self.variant_info = None
        self.genotype_fields = None
        self.genotype_info = None
        self.encoding = 'utf-8'
        self.logger = logger
        # for export only
        self.export_by_fields = ''
        self.order_by_fields = ''
        # for import only
        self.merge_by_cols = None
        #
        if os.path.isfile(name + '.fmt'):
            self.name = os.path.split(name)[-1]
            args = self.parseArgs(name + '.fmt', fmt_args)
            self.parseFMT(name + '.fmt', defaults=args) 
        elif name.endswith('.fmt') and os.path.isfile(name):
            self.name = os.path.split(name)[-1][:-4]
            args = self.parseArgs(name, fmt_args)
            self.parseFMT(name, defaults=args) 
        else:
            url = 'http://vtools.houstonbioinformatics.org/format/{}.fmt'.format(name)
            try:
                fmt = downloadFile(url, quiet=True)
            except Exception as e:
                raise ValueError('Failed to download format specification file {}.fmt'.format(name))
            self.name = name
            args = self.parseArgs(fmt, fmt_args)
            self.parseFMT(fmt, defaults=args)

    def parseArgs(self, filename, fmt_args):
        fmt_parser = ConfigParser.SafeConfigParser()
        fmt_parser.read(filename)
        parameters = fmt_parser.items('DEFAULT')
        parser = argparse.ArgumentParser(prog='vtools CMD --format {}'.format(os.path.split(filename)[-1]), description='''Parameters to override fields of
            existing format.''')
        self.parameters = []
        for par in parameters:
            # $NAME_comment is used for documentation only
            if par[0].endswith('_comment'):
                continue
            par_help = [x[1] for x in parameters if x[0] == par[0] + '_comment']
            self.parameters.append((par[0], par[1], par_help[0] if par_help else ''))
            parser.add_argument('--{}'.format(par[0]), help=self.parameters[-1][2],
                nargs='*', default=par[1])
        args = vars(parser.parse_args(fmt_args))
        for key in args:
            if type(args[key]) == list:
                args[key] = ','.join(args[key])
        return args

    def parseFMT(self, filename, defaults):
        parser = ConfigParser.SafeConfigParser()
        # this allows python3 to read .fmt file with non-ascii characters, but there is no
        # simple way to make it python2 compatible.
        #with open(filename, 'r', encoding='UTF-8') as inputfile:
        #    parser.readfp(inputfile)
        parser.read(filename)
        # sections?
        sections = parser.sections()
        if 'format description' not in sections:
            raise ValueError("Missing section 'format description'")
        #
        fields = []
        columns = []
        self.formatter = {}
        for section in sections:
            if section.lower() == 'format description':
                continue
            if section.lower() == 'field formatter':
                for item in parser.items(section, vars=defaults):
                    if item[0].startswith('fmt_'):
                        self.formatter[item[0][4:]] = item[1]
                continue
            if section.startswith('col_'):
                try:
                    items = [x[0] for x in parser.items(section, raw=True)]
                    for item in items:
                        if item.endswith('_comment'):
                            continue
                        if item not in ['field', 'adj', 'comment'] + defaults.keys():
                            raise ValueError('Incorrect key {} in section {}. Only field, adj and comment are allowed.'.format(item, section))
                    columns.append(
                        Column(index=int(section.split('_', 1)[1]),
                            field=parser.get(section, 'field', vars=defaults) if 'field' in items else '',
                            adj=parser.get(section, 'adj', vars=defaults) if 'adj' in items else None,
                            comment=parser.get(section, 'comment', raw=True) if 'comment' in items else '')
                        )
                except Exception as e:
                    raise ValueError('Invalid section {}: {}'.format(section, e))
            else:
                if not section.replace('_', '').isalnum():
                  raise ValueError('Illegal field name {}. Field names can only contain alphanumeric characters and underscores'.format(repr(section)))
                if section.upper() in SQL_KEYWORDS:
                  raise ValueError('Illegal field name. {} conflicts with SQL keywords'.format(repr(section)))
                try:
                    items = [x[0] for x in parser.items(section, raw=True)]
                    for item in items:
                        if item.endswith('_comment'):
                            continue
                        if item not in ['index', 'type', 'adj', 'comment'] + defaults.keys():
                            raise ValueError('Incorrect key {} in section {}. Only index, type, adj and comment are allowed.'.format(item, section))
                    fields.append(
                        Field(name=section,
                            index=parser.get(section, 'index', vars=defaults),
                            type=parser.get(section, 'type', vars=defaults),
                            adj=parser.get(section, 'adj', vars=defaults) if 'adj' in items else None,
                            comment=parser.get(section, 'comment', raw=True) if 'comment' in items else '')
                        )
                except Exception as e:
                    raise ValueError('Invalid section {}: {}'.format(section, e))
        #
        if len(fields) == 0:
            raise ValueError('No valid field is defined in format specification file {}'.format(self.name))
        #
        self.delimiter = '\t'
        #
        for item in parser.items('format description', vars=defaults):
            if item[0] == 'description':
                self.description = item[1]
            elif item[0] == 'delimiter':
                try:
                    # for None, '\t' etc
                    self.delimiter = eval(item[1])
                except:
                    # if failed, take things literally.
                    self.delimiter = item[1]
            elif item[0] == 'encoding':
                self.encoding = item[1]
            elif item[0] == 'merge_by':
                self.merge_by_cols = [x-1 for x in eval(item[1])]
            elif item[0] == 'export_by':
                self.export_by_fields = item[1]
            elif item[0] == 'sort_output_by':
                self.order_by_fields = item[1]
            elif item[0] in ['variant', 'position', 'range', 'genotype', 'variant_info', 'genotype_info']:
                setattr(self, item[0] if item[0].endswith('_info') else item[0]+'_fields', [x.strip() for x in item[1].split(',') if x.strip()])
        #
        # Post process all fields
        if (not not self.variant_fields) + (not not self.position_fields) + (not not self.range_fields) != 1:
            raise ValueError('Please specify one and only one of "variant=?", "position=?" or "range=?"')
        #
        if self.variant_fields:
            self.input_type = 'variant'
            self.ranges = [0, 4]
            self.fields = self.variant_fields
            if len(self.fields) != 4:
                raise ValueError('"variant" fields should have four fields for chr, pos, ref, and alt alleles')
        elif self.position_fields:
            self.input_type = 'position'
            self.ranges = [0, 2]
            self.fields = self.position_fields
            if len(self.fields) != 2:
                raise ValueError('"position" fields should have two fields for chr and pos')
        elif self.range_fields:
            self.input_type = 'range'
            self.ranges = [0, 3]
            self.fields = self.range_fields
            if len(self.fields) != 3:
                raise ValueError('"range" fields should have three fields for chr and starting and ending position')
        #
        if self.input_type != 'variant' and not self.variant_info:
            raise ValueError('Input file with type position or range must specify variant_info')
        if self.input_type != 'variant' and self.genotype_info:
            raise ValueError('Input file with type position or range can not have any genotype information.')
        if self.genotype_fields and len(self.genotype_fields) != 1:
            raise ValueError('Only one genotype field is allowed to input genotype for one or more samples.')
        #
        if self.variant_info:
            self.fields.extend(self.variant_info)
        self.ranges.append(self.ranges[-1] + (len(self.variant_info) if self.variant_info else 0))
        if self.genotype_fields:
            self.fields.extend(self.genotype_fields)
        self.ranges.append(self.ranges[-1] + (len(self.genotype_fields) if self.genotype_fields else 0))
        if self.genotype_info:
            self.fields.extend(self.genotype_info)
        self.ranges.append(self.ranges[-1] + (len(self.genotype_info) if self.genotype_info else 0))
        #
        # now, change to real fields
        for i in range(len(self.fields)):
            fld = [x for x in fields if x.name == self.fields[i]]
            if len(fld) != 1:
                #
                # This is a special case that allows users to use expressions as field....
                #
                if self.logger:
                    self.logger.warning('Field {} is not defined in format {}, some or all variants might fail to import.'.format(self.fields[i], filename))
                self.fields[i] = Field(name=self.fields[i], index=None, adj=None, type=None, comment='')
            else:
                self.fields[i] = fld[0]
        # other fields?
        self.other_fields = [x for x in fields if x not in self.fields]
        #
        # columns definition
        self.columns = []
        for idx in range(len(columns)):
            # find column
            try:
                col = [x for x in columns if x.index == idx + 1][0]
            except Exception as e:
                raise ValueError('Cannot find column {} from format specification: {}'.format(idx + 1, e))
            self.columns.append(col)

    def describe(self):
        print('Format:      {}'.format(self.name))
        if self.description is not None:
            print('Description: {}'.format('\n'.join(textwrap.wrap(self.description,
                initial_indent='', subsequent_indent=' '*2))))
        #
        print('\nColumns:')
        if self.columns:
            for col in self.columns:
                print('  {:12} {}'.format(str(col.index), '\n'.join(textwrap.wrap(col.comment,
                    subsequent_indent=' '*15))))
            if self.formatter:
                print('Formatters are provided for fields: {}'.format(', '.join(self.formatter.keys())))
        else:
            print('  None defined, cannot export to this format')
        #
        if self.input_type == 'variant':
            print('\n{0}:'.format(self.input_type))
        for fld in self.fields[self.ranges[0]:self.ranges[1]]:
            print('  {:12} {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                subsequent_indent=' '*15))))
        if self.ranges[1] != self.ranges[2]:
            print('\nVariant info:')
            for fld in self.fields[self.ranges[1]:self.ranges[2]]:
                print('  {:12} {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*15))))
        if self.ranges[2] != self.ranges[3]:
            print('\nGenotype:')
            for fld in self.fields[self.ranges[2]:self.ranges[3]]:
                print('  {:12} {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*15))))
        if self.ranges[3] != self.ranges[4]:
            print('\nGenotype info:')
            for fld in self.fields[self.ranges[3]:self.ranges[4]]:
                print('  {:12} {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*15))))
        if self.other_fields:
            print('\nOther fields (usable through parameters):')
            for fld in self.other_fields:
                print('  {:12} {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*15))))
        if self.parameters:
            print('\nFormat parameters:')
            for item in self.parameters:
                print('  {:12} {}'.format(item[0],  '\n'.join(textwrap.wrap(
                    '{} (default: {})'.format(item[2], item[1]),
                    subsequent_indent=' '*15))))
        else:
            print('\nNo configurable parameter is defined for this format.\n')



class AnnoDBWriter:
    '''
    A class to initiate and insert annotation database
    '''
    def __init__(self, name, fields, anno_type, description, version, build, logger, use_existing_db=False,
        overwrite_existing_fields=False):
        self.logger = logger
        self.name = name
        self.fields = fields
        self.anno_type = anno_type
        self.description = description
        self.version = version
        self.build = build
        # create database and import file
        self.db = DatabaseEngine()
        if not use_existing_db or not os.path.isfile(self.name + '.DB'):
            self.update_existing = False
            self.createAnnoDB()
        else:
            self.update_existing = True
            self.updateAnnoDB(overwrite_existing_fields)

    def createAnnoDB(self):
        # remove database if already exist
        self.db.removeDatabase(self.name)
        # create a new one
        self.db.connect(self.name)
        cur = self.db.cursor()
        #
        # creating the field table
        self.logger.debug('Creating {}_field table'.format(self.name))
        self.createFieldsTable()
        #
        for field in self.fields:
            cur.execute('INSERT INTO {0}_field (name, field, type, comment) VALUES ({1},{1},{1},{1});'.format(self.name, self.db.PH),
                (field.name, field.index, field.type, field.comment))
        self.db.commit()
        #
        # creating the info table
        self.logger.debug('Creating {}_info table'.format(self.name))
        query = 'INSERT INTO {0}_info VALUES ({1},{1});'.format(self.name, self.db.PH)
        self.createInfoTable()
        cur.execute(query, ('name', self.name))
        cur.execute(query, ('anno_type', self.anno_type))
        cur.execute(query, ('description', self.description))
        cur.execute(query, ('version', self.version))
        cur.execute(query, ('build', str(self.build)))
        self.db.commit()
        self.logger.debug('Creating table {}'.format(self.name))
        self.createAnnotationTable()
    
    def updateAnnoDB(self, overwrite_existing_fields):
        self.db.connect(self.name)
        for table in [self.name, self.name + '_field', self.name + '_info']:
            if not self.db.hasTable(table):
                raise ValueError('Existing file {}.DB is not a valid annotation database.'.format(self.name))
        # get linked fields
        cur = self.db.cursor()
        cur.execute('SELECT * from {}_info;'.format(self.name))
        for rec in cur:
            if rec[0] == 'anno_type':
                if rec[1] != 'field':
                    raise ValueError('Existing database is not field-based. Cannot add results to it.')
            elif rec[0] == 'build':
                if self.build != eval(rec[1]):
                    raise ValueError('Existing database has different linking fields (existing: {}, required: {}).'.format(self.build, rec[1]))
        # get existing fields
        cur.execute('SELECT name, field, "", type, comment from {}_field;'.format(self.name))
        cur_fields = []
        for rec in cur:
            cur_fields.append(Field(*rec))
        # add new fields
        s = delayedAction(self.logger.info, 'Adding fields to existing result database')
        for field in self.fields:
            # name already exist
            if field.name in [x.name for x in cur_fields]:
                cf = [x for x in cur_fields if x.name == field.name][0] 
                if field.type != cf.type:
                    raise ValueError('Type mismatch for new field {}: existing {}, new {}'.format(field.name, cf.type, field.type))
                if overwrite_existing_fields:
                    self.logger.warning('Results in field {} will be overwritten.'.format(field.name))
                else:
                    if field.name not in [x for tmp in self.build.values() for x in tmp]: 
                        raise ValueError('Cannot modify database {} because field {} already exists. '
                            'Please use test option --name to add a new suffix to this field, '
                            'write the results to a different database (option --to_db), or use '
                            'option --update to force updating the existing fields.'.format(self.name, field.name))
            else:
                # add new field
                cur.execute('INSERT INTO {0}_field (name, field, type, comment) VALUES ({1},{1},{1},{1});'.format(self.name, self.db.PH),
                    (field.name, field.index, field.type, field.comment))
                cur.execute('ALTER TABLE {} ADD {} {};'.format(self.name, field.name, field.type))
        #
        self.db.commit()
        del s

    def createFieldsTable(self):
        '''Create table name_fields'''
        cur = self.db.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS {}_field (
            name VARCHAR(40),
            field INT,
            type VARCHAR(80),
            comment VARCHAR(256),
            missing_entries INT,
            distinct_entries INT,
            min_value INT,
            max_value INT
        )'''.format(self.name))

    def createInfoTable(self):
        '''Create table name_fields'''
        cur = self.db.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS {}_info (
            name VARCHAR(40),
            value VARCHAR(1024)
        )'''.format(self.name)) 

    def createAnnotationTable(self):
        'Create an annotation table '
        items = []
        for build in self.build.keys():
            if build != '*':
                items.append('{0}_bin INTEGER'.format(build))
        for field in self.fields:
            items.append('{0} {1}'.format(field.name, field.type))
        query = '''CREATE TABLE IF NOT EXISTS {} ('''.format(self.name) + \
            ',\n'.join(items) + ');'
        self.logger.debug('Creating annotation table {} using query\n{}'.format(self.name, query))
        cur = self.db.cursor()
        try:
            cur.execute(query)
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to create table')
    
    def finalize(self):
        '''Create index and get statistics of the database'''
        cur = self.db.cursor()
        # creating indexes
        s = delayedAction(self.logger.info, 'Creating indexes (this can take quite a while)')
        # creates index for each link method
        for key in self.build.keys():
            if key != '*':
                if not self.db.hasIndex('{}_idx'.format(key)):
                    cur.execute('''CREATE INDEX {0}_idx ON {1} ({0}_bin ASC, {2});'''\
                      .format(key, self.name,  ', '.join(['{} ASC'.format(x) for x in self.build[key]])))
            else:
                if not self.db.hasIndex('{}_idx'.format(self.name)):
                    cur.execute('''CREATE INDEX {0}_idx ON {0} ({1});'''\
                        .format(self.name,  ', '.join(['{} ASC'.format(x) for x in self.build[key]])))
        del s
        s = delayedAction(self.logger.info, 'Analyzing and tuning database ...')
        # This is only useful for sqlite
        self.db.analyze()
        # calculating database statistics
        cur.execute('SELECT COUNT(*) FROM (SELECT DISTINCT {} FROM {});'.format(', '.join(self.build.values()[0]), self.name))
        count = cur.fetchone()[0]
        cur.execute('INSERT INTO {0}_info VALUES ({1}, {1});'.format(self.name, self.db.PH), ('distinct_keys', str(count)))
        cur.execute('INSERT INTO {0}_info VALUES ({1}, {1});'.format(self.name, self.db.PH), ('num_records', self.db.numOfRows(self.name)))
        del s
        for field in self.fields:
            s = delayedAction(self.logger.info, 'Calculating column statistics for field {}'.format(field.name))
            cur.execute('SELECT COUNT(*) FROM {1} WHERE {0} is NULL;'.format(field.name, self.name))
            missing = cur.fetchone()[0]
            cur.execute('UPDATE {0}_field SET missing_entries={1} WHERE name="{2}";'.format(self.name, self.db.PH, field.name),
                (missing,))
            if 'int' in field.type.lower() or 'float' in field.type.lower():
                cur.execute('SELECT COUNT(DISTINCT {0}), MIN({0}), MAX({0}) FROM {1} WHERE {0} IS NOT NULL;'.format(field.name, self.name))
                res = cur.fetchone()
                cur.execute('UPDATE {0}_field SET distinct_entries={1}, min_value={1}, max_value={1} WHERE name={1};'.format(
                    self.name, self.db.PH), (res[0], res[1], res[2], field.name))
            else:
                cur.execute('SELECT COUNT(DISTINCT {0}) FROM {1};'.format(field.name, self.name))
                res = cur.fetchone()
                cur.execute('UPDATE {0}_field SET distinct_entries={1} WHERE name={1};'.format(
                    self.name, self.db.PH), (res[0], field.name))
            del s
        self.db.commit()

#  Project management
#
class Project:
    '''
    A project maintains the following tables:

    1. Table "project", which has the following information:

            name:       name of property
            value:      value of property

    2. Table "filename", which stores names of all input files:

            file_id:    file ID
            filename:   filename

    3. Table "variant" has the following fields:

            variant_id:  varaint ID
            bin:        UCSC bin
            chr:        chromosome
            pos:        position
            ref:        reference allele
            alt:        alternative allele

            OPTIONAL (added by liftOver)
            alt_chr:    chromosome on alternative reference genome
            alt_pos:    position on alternative reference genome

            OPTIONAL (added by sampling tools)
            freq1:       Sample mutant frequency, if available
            freq2:       Sample mutant frequency, if available

       There can be mulitple variant tables in a project. The master
       variant table will be named "variant". The last three
       fields are maintained by liftOver and sample modules.

    If genetic samples are involved, the following tables will be
    created:

    1. Table "sample" is used to store the source of variant and genotype
       information. A file will not be re-imported if its name appears in this
       table. This table has the following fields:

            sample_id:   sample id
            file_id:     ID of filename
            sample_name: name of sample specified in .vcf file, will be empty
                         if no sample information is available (e.g. a bed file
                         is provided).
            PHENOTYPE    phenotype might be added by command 'import_phenotype'

    2. Table "genotype_$sample_id" stores variants of each sample
            variant_id:  ID of variant
            type:       A numeric (categorical) value

      These tables are stored in a separate database $name_genotype in order to
      keep the project database small.

    If there are meta information for sample, and or genotype, the following
    table will be created and used.


    1. Table "variant_meta" is used to store meta information for each variant.
       It has the same length as the "variant" table. The meta information includes
       INFO and FORMAT fields specified in the vcf files, rsname and other
       information stored in bed files.

            ID:         variant ID
            INFO1:      ..
            INFO2:      ..
            ...


    2. Table "sample_meta" stores additional information for each sample. For VCF
       files, a field vcf_meta is used to store all meta information. Future
       extension will allow the input of phenotype information (e.g. case
       control status).

            ID:          sample id
            vcf_meta:    basically header lines of the vcf file
            INFO1:       ..
            INFO2:       ..

    '''
    def __init__(self, name=None, build=None, new=False, verbosity=None, verify=True, **kwargs):
        '''Create a new project (--new=True) or connect to an existing one.'''
        files = glob.glob('*.proj')
        if new: # new project
            if len(files) > 0:
                if name + '.proj' in files:
                    raise ValueError('Project {0} already exists. Please use '.format(name) + \
                        'option --force to remove it if you would like to start a new project.')
                else:
                    raise ValueError('A project can only be created in a directory without another project.')
            if name is None:
                raise ValueError('A new project must have a name')
            # if a file is specified...
            elif '.' in name or os.path.split(name)[0]:
                raise ValueError('A project name cannot have extension or path')
            elif name[0].isdigit():
                raise ValueError('A project name cannot start with a number.')
        else: # exisitng project
            if len(files) == 0:
                raise ValueError('Cannot find any project in the current directory.')
            elif len(files) > 1:
                raise ValueError('More than one project exists in the current directory.')
            if name is None:
                name = files[0][:-5]
            elif name != files[0][:-5]:
                raise ValueError('Another project {} already exists in the current directory'.format(files[0]))
        #
        self.name = name
        self.proj_file = self.name + '.proj'
        # version of vtools, useful when opening a project created by a previous
        # version of vtools.
        self.version = VTOOLS_VERSION
        #
        # create a temporary directory
        try:
            if not os.path.isdir('cache'):
                os.mkdir('cache')
            self.temp_dir = 'cache'
        except:
            self.temp_dir = tempfile.mkdtemp()
        # set global verbosity level and temporary directory
        setOptions(verbosity=verbosity, temp_dir=self.temp_dir)
        #
        # create a logger
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        levels = {
            '0': logging.ERROR,
            '1': logging.INFO,
            '2': logging.DEBUG
            }
        #
        if verbosity is None and not new:
            # try to get saved verbosity level
            try:
                self.db = DatabaseEngine()
                self.db.connect(self.proj_file)
                self.verbosity = self.loadProperty('verbosity')
            except:
                self.verbosity = verbosity
        else:
            self.verbosity = verbosity
        #
        cout.setLevel(logging.INFO if self.verbosity is None else levels[self.verbosity[0]])
        cout.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        self.logger.addHandler(cout)
        # output to a log file
        ch = logging.FileHandler(self.name + '.log', mode='w' if new else 'a')
        # NOTE: debug informaiton is always written to the log file
        ch.setLevel(logging.DEBUG if self.verbosity is None or len(self.verbosity) == 1 else levels[self.verbosity[1]])
        ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
        self.logger.addHandler(ch)
        if new:
            self.create(build=build, **kwargs)
        else:
            self.open(verify)

    def create(self, build, **kwargs):
        '''Create a new project'''
        # open the project file
        self.logger.info(VTOOLS_COPYRIGHT)
        self.logger.info(VTOOLS_CITATION)
        self.logger.info(VTOOLS_CONTACT)
        self.logger.info('Creating a new project {}'.format(self.name))
        self.db = DatabaseEngine(engine='sqlite3', batch=kwargs.get('batch', 10000))
        self.db.connect(self.proj_file)
        #
        engine = kwargs.get('engine', 'sqlite3')
        if engine == 'mysql':
            # project information table
            self.createProjectTable()
            self.saveProperty('version', self.version)
            self.saveProperty('engine', engine)
            self.saveProperty('host', kwargs.get('host'))
            self.saveProperty('user', kwargs.get('user'))
            # FIXME: I am saving passwd as clear text here...
            self.saveProperty('passwd', kwargs.get('passwd'))
            self.saveProperty('batch', kwargs.get('batch', 10000))
            self.saveProperty('verbosity', self.verbosity)
            # turn to the real online engine
            self.logger.debug('Connecting to mysql server')
            self.db.commit()
            self.db = DatabaseEngine(**kwargs)
            if self.db.hasDatabase(self.name):
                raise ValueError('A mysql database named {} already exists. Please remove it from the mysql database before continuing.'.format(self.name))
            else:
                self.db.connect(self.name)
        #
        self.build = build
        self.alt_build = None
        # no meta information for now
        self.variant_meta = None
        self.sample_meta = None
        self.annoDB = []
        #
        # Initialize the core tables
        self.logger.debug('Creating core tables')
        self.createProjectTable()
        self.saveProperty('engine', engine)
        self.saveProperty('version', self.version)
        self.saveProperty('batch', kwargs.get('batch', 10000))
        self.saveProperty('verbosity', self.verbosity)
        self.saveProperty('name', self.name)
        self.saveProperty('build', self.build)
        self.saveProperty('alt_build', self.alt_build)
        self.saveProperty('annoDB', str(self.annoDB))
        #
        self.createFilenameTable()
        self.createMasterVariantTable()
        self.createSampleTableIfNeeded()

    def open(self, verify=True):
        '''Open an existing project'''
        # open the project file
        self.logger.debug('Opening project {}'.format(self.proj_file))
        self.db = DatabaseEngine()
        self.db.connect(self.proj_file)
        if not self.db.hasTable('project'):
            raise ValueError('Invalid project database (missing project table)')
        #
        self.batch = self.loadProperty('batch')
        self.batch = 10000 if self.batch is None else int(self.batch)
        #
        # get connection parameters
        if self.loadProperty('engine') == 'mysql':
            self.db = DatabaseEngine(engine='mysql',
                host=self.loadProperty('host'),
                user=self.loadProperty('user'),
                passwd=self.loadProperty('passwd'))
            self.db.connect(self.name)
        else: 
            # pass things like batch ... and re-connect
            self.db = DatabaseEngine(engine='sqlite3', batch=self.batch)
            self.db.connect(self.proj_file)
        #
        # existing project
        cur = self.db.cursor()
        self.version = self.loadProperty('version', '1.0')
        self.build = self.loadProperty('build')
        self.alt_build = self.loadProperty('alt_build')
        self.annoDB = []
        for db in eval(self.loadProperty('annoDB', '[]')):
            try:
                db_name = os.path.split(db)[-1].split('-')[0]
                linked_by = eval(self.loadProperty('{}_linked_by'.format(db_name), default='[]'))
                anno_type = self.loadProperty('{}_anno_type'.format(db_name), default='None')
                linked_fields = eval(self.loadProperty('{}_linked_fields'.format(db_name), default='None'))
                self.annoDB.append(AnnoDB(self, db, linked_by, anno_type, linked_fields))
                self.db.attach(db)
            except Exception as e:
                self.logger.warning(e)
                self.logger.warning('Cannot open annotation database {}'.format(db))
        #
        # get existing meta information
        # FIXME: these are not handled correctly for now
        self.variant_meta = self.db.getHeaders('variant_meta')
        self.sample_meta = self.db.getHeaders('sample_meta')
        #
        if verify:
            self.checkIntegrity()
        # 
        self.analyze()

    def checkIntegrity(self):
        '''Check if the project is ok...(and try to fix it if possible)'''
        for table in ['project', 'filename', 'sample', 'variant']:
            if not self.db.hasTable(table):
                raise RuntimeError('Corrupted project: missing table {}'.format(table))
        #
        headers = self.db.getHeaders('variant')
        if self.alt_build is not None:
            if not ('alt_bin' in headers and 'alt_chr' in headers and 'alt_pos' in headers):
                self.logger.warning('Disable alternative reference genome because of missing column {}.'.format(
                    ', '.join([x for x in ('alt_bin', 'alt_chr', 'alt_pos') if x not in headers])))
                self.alt_build = None
                self.saveProperty('alt_build', None)
        #
        # missing index on master variant table
        if not self.db.hasIndex('variant_index'):
            self.logger.warning('Missing index on master variant table. Trying to rebuild it.')
            self.createIndexOnMasterVariantTable()
            if not self.db.hasIndex('variant_index'):
                raise RuntimeError('Corrupted project: failed to create index on master variant table.')

    def analyze(self, force=False):
        '''Automatically analyze project to make sure queries are executed optimally.
        '''
        if self.db.engine != 'sqlite3':
            return
        cur = self.db.cursor()
        tables = self.db.tables()
        cur = self.db.cursor()
        s = delayedAction(self.logger.info, 'Analyzing project')
        for tbl in tables:
            analyzed = True
            if not force:
                # try to figure out if the table has been analyzed
                try:
                    cur.execute('SELECT tbl FROM sqlite_stat1 WHERE tbl=?', (tbl,))
                    analyzed = len(cur.fetchall()) == 1
                except:
                    analyzed = False
            if force or not analyzed:
                cur.execute('ANALYZE {}'.format(tbl))
        self.db.commit()
        del s
        

    def useAnnoDB(self, db):
        '''Add annotation database to current project.'''
        # DBs in different paths but with the same name are considered to be the same.
        self.logger.info('Using annotation DB {} in project {}.'.format(db.name, self.name))
        self.logger.info(db.description)
        if db.name not in [x.name for x in self.annoDB]:
            self.annoDB.append(db)
            self.saveProperty('annoDB', str([os.path.join(x.dir, x.filename) for x in self.annoDB]))
        # an annotation database might be re-used with a different linked_field
        self.saveProperty('{}_linked_by'.format(db.name), str(db.linked_by))
        self.saveProperty('{}_anno_type'.format(db.name), str(db.anno_type))
        self.saveProperty('{}_linked_fields'.format(db.name), str(db.build))

    def close(self):
        '''Write everything to disk...'''
        self.db.commit()
        if self.temp_dir != 'cache' and os.path.isdir(self.temp_dir):
            shutil.rmtree(self.temp_dir)
        
    def loadProperty(self, key, default=None):
        '''Retrieve property from the project table'''
        cur = self.db.cursor()
        cur.execute('SELECT value FROM project WHERE name={0};'.format(self.db.PH), (key,))
        try:
            return cur.fetchone()[0]
        except Exception as e:
            self.logger.debug(e)
            self.logger.warning('Failed to retrieve value for project property {}'.format(key))
            self.saveProperty(key, default)
            return default

    def saveProperty(self, key, value):
        '''Save property in the project table'''
        cur = self.db.cursor()
        try:
            cur.execute('SELECT value FROM project WHERE name={};'.format(self.db.PH), (key,))
            if cur.fetchall():
                cur.execute('UPDATE project SET value={0} WHERE name={0};'.format(self.db.PH), (value, key))
            else:
                cur.execute('INSERT INTO project VALUES ({0}, {0});'.format(self.db.PH), (key, value))
        except Exception as e:
            self.logger.debug(e)
            raise
        self.db.commit()

    def remove(self):
        '''Remove the current project'''
        # step 1: remove database
        if self.db.engine == 'mysql':
             self.logger.info('Removing database {}'.format(self.name))
             self.db.removeDatabase(self.name)
        # step 2: remove genotype
        if self.db.hasDatabase(self.name + '_genotype'):
            self.logger.info('Removing database {}_genotype'.format(self.name))
            self.db.removeDatabase(self.name + '_genotype')
        # step 3: remove temp database
        if self.db.hasDatabase(self.name + '_temp'):
            self.logger.info('Removing database {}_temp'.format(self.name))
            self.db.removeDatabase(self.name + '_temp')
        # step 4: remove project file
        self.logger.info('Removing project file {}'.format(self.proj_file))
        os.remove(self.proj_file)
        # step 5: remove log file
        self.logger.info('Removing log file {}'.format(self.proj_file[:-5] + '.log'))
        os.remove(self.proj_file[:-5] + '.log')
    

    #
    # Support for python with statement
    #
    #
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        '''Write everything to disk...'''
        self.close()
        return exc_val is None

    #
    # set property
    #
    def setRefGenome(self, build):
        if self.build is not None and self.build != build:
            self.logger.error('Cannot change reference genome of an existing project.')
            sys.exit(1)
        self.build = build
        self.saveProperty('build', build)

    #
    # Functions to create core and optional tables.
    #
    def createProjectTable(self):
        self.logger.debug('Creating table project')
        cur = self.db.cursor()
        cur.execute('''\
            CREATE TABLE project (
                name VARCHAR(40) NOT NULL PRIMARY KEY,
                value VARCHAR(256) NULL
            )''')

    def createFilenameTable(self):
        self.logger.debug('Creating table filename')
        cur = self.db.cursor()
        cur.execute('''\
            CREATE TABLE filename (
                file_id INTEGER PRIMARY KEY {0},
                filename VARCHAR(256) NOT NULL,
                header BLOB NULL
            )'''.format(self.db.AI))
        # create index
        try:
            cur.execute('''CREATE UNIQUE INDEX filename_index ON filename (filename ASC);''')
        except Exception as e:
            # the index might already exists
            return

    def createMasterVariantTable(self, fields=[]):
        '''Create a variant table with name. Fail if a table already exists.'''
        self.logger.debug('Creating table variant')
        #
        # ref and alt are 'VARCHAR' to support indels. sqlite database ignores VARCHAR length
        # so it can support really long indels. MySQL will have trouble handling indels that
        # are longer than 255 bp.
        #
        self.db.execute('''\
            CREATE TABLE variant (
                variant_id INTEGER PRIMARY KEY {0},
                bin INTEGER NULL,
                chr VARCHAR(20) NULL,
                pos INTEGER NULL,
                ref VARCHAR(255) NOT NULL,
                alt VARCHAR(255) NOT NULL {1});'''.format(self.db.AI,
                ''.join([', {} {}\n'.format(x,y) for x,y in fields])))
        self.createIndexOnMasterVariantTable()

    def createIndexOnMasterVariantTable(self, quiet=False):
        # create indexes
        #
        if not quiet:
            s = delayedAction(self.logger.info, 'Creating indexes on master variant table. This might take quite a while.')
        try:
            #
            # Index on the primary reference genome is UNIQUE when there is no alternative reference
            # genome. If there is, multiple variants from the alternative reference genome might
            # be mapped to the same coordinates on the primary reference genome. I have tried to 
            # set some of the coordinates to NULL, but the unqiueness problem becomes a problem
            # across projects. For example, 
            #
            # If c_19 and c_19a both map to c_18 in one project
            #
            #   c_18   c_19
            #   None   c_19a
            #
            # another project might have
            #
            #   c_18  c_19a
            #
            # although they are the same variant. A better solution is therefore to keep
            #
            #   c_18   c_19
            #   c_18   c_19a
            #
            # in both projects.
            #
            if self.alt_build is not None:
                self.db.execute('''CREATE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
            else:
                self.db.execute('''CREATE UNIQUE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists, this does not really matter
            self.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
        try:
            if self.alt_build:
                #
                # Index on alternative reference genome is not unique because several variants might
                # be mapped to the same coordinates in the alternative reference genome. 
                #
                self.db.execute('''CREATE INDEX variant_alt_index ON variant (alt_bin ASC, alt_chr ASC, alt_pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists, this does not really matter
            self.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
        if not quiet:
            del s

    def dropIndexOnMasterVariantTable(self):
        # before bulk inputting data, it is recommended to drop index.
        #
        # NOTE: for mysql, it might be better to use alt index disable/rebuild
        #
        s = delayedAction(self.logger.info, 'Dropping indexes of master variant table. This might take quite a while.')
        try:
            # drop index has different syntax for mysql/sqlite3.
            self.db.dropIndex('variant_index', 'variant')
        except Exception as e:
            # the index might not exist
            self.logger.debug(e)
        #
        try:
            if self.alt_build:
                self.db.dropIndex('variant_alt_index', 'variant')
        except Exception as e:
            # the index might not exist
            self.logger.debug(e)
        del s

    def createVariantTable(self, table, temporary=False):
        '''Create a variant table with name. Fail if a table already exists.
        '''
        if table == 'variant':
            raise ValueError('This function cannot be used to create a master variant table')
        self.logger.debug('Creating table {}'.format(table))
        self.db.execute('''CREATE {0} TABLE {1} (
                variant_id INTEGER PRIMARY KEY
            );'''.format('TEMPORARY' if temporary else '', table))
        self.db.commit()

    def createSampleTableIfNeeded(self, fields=[], table='sample'):
        if self.db.hasTable(table):
            return
        self.logger.debug('Creating table {}'.format(table))
        cur = self.db.cursor()
        query = '''\
            CREATE TABLE IF NOT EXISTS {0} (
                sample_id INTEGER PRIMARY KEY {1},
                file_id INTEGER NOT NULL,
                sample_name VARCHAR(256) NULL'''.format(table, self.db.AI)
        for (n, t) in fields:
            query += ',\n{} {} NULL'.format(n, t)
        query += ');'
        cur.execute(query)
        self.db.commit()

    def createNewSampleVariantTable(self, cur, table, genotype=True, fields=[]):
        '''Create a table ``sample_variant`` to store vcf data'''
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS {0} (
                variant_id INT NOT NULL
            '''.format(table) + 
            (', GT INT' if genotype else '') + 
            ''.join([', {} {}'.format(f.name, f.type) for f in fields]) + ');'
         )
        #
        # NOTE
        #  Adding an autoincrement field increase importing time from 1:30 to 3:00
        #
        #    ID INTEGER PRIMARY KEY {1},
        #
        # NOTE:
        #   Having such a key will significantly reduce the performance when the
        #   table gets larger. The speed will decrease from 25k/s to 1.7k/s after
        #   reading 44 files. Removing the index will lead to a stable performance.
        #
        # NOTE:
        #   Perhaps we should separate the files by sample ID? In this way,
        #   calculating sample allele frequency will be slower, but might be
        #   tolerable if we index by variant_id.
        #
        #try:
            #cur = self.db.cursor()
            #cur.execute('''CREATE INDEX {0}_index ON {0} (variant_id ASC, sample_id ASC);'''.format(table))
        #except Exception as e:
            # key might already exists
            #self.logger.debug(e)
            #pass

    def createVariantMetaTableIfNeeded(self):
        if self.variant_meta is None:
            self.logger.warning('Tried to create a variant meta table without valid meta fields')
            return
        if self.db.hasTable('variant_meta'):
            return
        self.logger.debug('Creating table variant_meta')
        cur = self.db.cursor()
        query = '''\
            CREATE TABLE IF NOT EXISTS variant_meta (
                variant_id INTEGER PRIMARY KEY {0},
                bin INTEGER NULL'''.format(self.db.AI)
        if self.alt_build is not None:
            query += ''',
            chr VARCHAR(20) NULL,
            pos INT NULL'''
        # FIXME: handle this
        if self.variant_meta is not None:
            pass
        #
        query += ');'
        cur.execute(query)

    def createSampleMetaTableIfNeeded(self):
        if self.sample_meta is None:
            self.logger.warning('Tried to create a sample meta table without valid meta fields')
            return
        if self.db.hasTable('sample_meta'):
            return
        self.logger.debug('Creating table sample_meta')
        cur = self.db.cursor()
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS sample_meta (
                ID INTEGER PRIMARY KEY {0},
                vcf_meta BLOB NULL
            );'''.format(self.db.AI))


    #
    # Project management
    #
    def isVariantTable(self, table):
        return (not table.startswith('__')) and self.db.hasTable(table) and self.db.getHeaders(table)[0] == 'variant_id'
        
    def getVariantTables(self):
        '''Return all variant tables'''
        return [table for table in self.db.tables() if self.isVariantTable(table)]

    def removeVariantTable(self, table):
        '''
        Remove specified variant table. Master table cannot be removed.
        '''
        if table == 'variant':
            raise ValueError('Master variant table cannot be removed.')
        if not self.isVariantTable(table):
            raise ValueError('{} is not found or is not a variant table.'.format(table))
        self.logger.info('Removing table {}'.format(table))
        self.db.removeTable(table)

    def selectSampleByPhenotype(self, cond):
        '''Select samples by conditions such as "aff=1", return IDs as a sorted list'''
        cur = self.db.cursor()
        try:
            query = 'SELECT sample_id FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id {};'.format(' WHERE {}'.format(cond) if cond.strip() else '')
            self.logger.debug('Select samples using query {}'.format(query))
            cur.execute(query)
            IDs = [x[0] for x in cur.fetchall()]
            IDs.sort()
            # return a tuple to avoid future change of order
            return tuple(IDs)
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to retrieve samples by condition "{}"'.format(cond))

    def removeSamples(self, IDs):
        '''Remove sample and their genotype, but not variants'''
        cur = self.db.cursor()
        samples = defaultdict(list)
        for ID in IDs:
            cur.execute('SELECT filename.filename, sample.sample_name FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id WHERE sample.sample_id = {};'.format(self.db.PH), (ID,))
            res = cur.fetchone()
            samples[res[0]].append(res[1])
        for f in samples:
            cur.execute('SELECT filename.filename, count(sample.sample_id) FROM filename LEFT OUTER JOIN sample on sample.file_id = filename.file_id WHERE filename.filename = {};'.format(self.db.PH),
                (f,))
            rec = cur.fetchone()
            if rec[1] == len(samples[f]):
                self.logger.info('Removing {1} and all its samples ({0})'.format('{} samples'.format(len(samples[f])) if len(samples[f]) > 1 else 'sample {}'.format(samples[f][0]), f)) 
                cur.execute('DELETE FROM filename WHERE filename = {}'.format(self.db.PH), (f,))
            else:
                self.logger.info('Removing {} from {}'.format('{} samples'.format(len(samples[f])) if len(samples[f]) > 1 else 'sample {}'.format(samples[f][0]), f)) 
        for ID in IDs:
            cur.execute('DELETE FROM sample WHERE sample_id = {};'.format(self.db.PH), (ID,))
            self.db.removeTable('{}_genotype.genotype_{}'.format(self.name, ID))        
        self.db.commit()
        
    def removeVariants(self, table):
        '''Remove variants from a project belong to table'''
        if table == 'variant':
            raise ValueError('Cannot remove variants from table variant because it will remove all variants')
        if not self.isVariantTable(table):
            raise ValueError('{} is not found or is not a variant table.'.format(table))
        cur = self.db.cursor()
        for t in self.getVariantTables():
            if t.lower() == table.lower():
                continue
            cur.execute('DELETE FROM {} WHERE variant_id IN (SELECT variant_id FROM {});'.format(t, table))
            self.logger.info('{} variants are removed from table {}'.format(cur.rowcount, t))
        # get sample_ids
        cur.execute('SELECT sample_id FROM sample;')
        IDs = [x[0] for x in cur.fetchall()]
        for ID in IDs:
            if not self.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.name, ID)):
                cur.execute('CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'.format(self.name, ID))
            cur.execute('DELETE FROM {}_genotype.genotype_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                .format(self.name, ID, table))
            self.logger.info('{} genotypes are removed from sample {}'.format(cur.rowcount, ID))
        # remove the table itself
        self.logger.info('Removing table {} itself'.format(table))
        self.db.removeTable(table)

    def removeGenotypes(self, cond):
        '''Remove genotype according to certain conditions'''
        # get sample_ids
        cur = self.db.cursor()
        cur.execute('SELECT sample_id FROM sample;')
        IDs = [x[0] for x in cur.fetchall()]
        self.logger.info('Removing genotypes from {} samples using criteria "{}"'.format(len(IDs), cond))
        for ID in IDs:
            cur.execute('DELETE FROM {}_genotype.genotype_{} WHERE {};'\
                .format(self.name, ID, cond))
            self.logger.info('Removing {} genotypes from sample {}'.format(cur.rowcount, ID))

    def summarize(self):
        '''Summarize key features of the project
        '''
        # FIXME: more summary
        info =  'Project name:                {}\n'.format(self.name)
        info += 'Primary reference genome:    {}\n'.format(self.build)
        info += 'Secondary reference genome:  {}\n'.format(self.alt_build)
        info += 'Database engine:             {}\n'.format(self.db.engine)
        info += 'Variant tables:              {}\n'.format('\n'.join([' '*29 + x for x in self.getVariantTables()]).lstrip())
        info += 'Annotation databases:        {}\n'.format('\n'.join([' '*29 + os.path.join(x.dir, x.name) \
            + (' ({})'.format(x.version) if x.version else '') for x in self.annoDB]).lstrip())
        return info

    #
    # temporary table which are created in separate database
    # because it is very slow for delite to remove temporary table
    # def getTempTable(self):
    #     '''Return a table name that does not exist in database'''
    #     self.db.attach('{}_temp'.format(self.name))
    #     while True:
    #         name = '{}_temp._tmp_{}'.format(self.name, random.randint(0, 10000))
    #         if not self.db.hasTable(name):
    #             return name   

    def checkFieldName(self, name, exclude=None):
        '''Check if a field name has been used, or is the SQL keyword'''
        if name.upper() in SQL_KEYWORDS:
            raise ValueError("Field name '{}' is not allowed because it is a reserved word.".format(name))
        for table in self.getVariantTables(): #  + ['sample', 'filename']:
            if table == exclude:
                continue
            if name.lower() in [x.lower() for x in self.db.getHeaders(table)]:
                raise ValueError("Field name '{}' is not allowed because it is already used in table {}".format(name, table))
        
    #
    # Handling field query
    #
    def linkFieldToTable(self, field, variant_table):
        '''Return one or more FieldConnections that link a field to a variant_table'''
        # if field is specified by table.field, good
        if '.' in field:
            # possibly two dots (db.table.field), but never seen them.
            table, field = [x.lower() for x in field.rsplit('.', 1)]
            if self.db.hasTable(table):
                if variant_table.lower() == table.lower():
                    # [] connection
                    # the field is in variant_table itself, no link is needed
                    return [FieldConnection(
                        field='{}.{}'.format(table, field),
                        table='',  # no need to join table
                        link='') ]
                elif variant_table.lower() == 'variant':
                    # [] <-> [master] connection
                    # outputting master variant table field from a non-master variant table
                    return [FieldConnection(
                        field='{}.{}'.format(table, field),
                        table='variant', # need to join table with the master variant table
                        link='{}.variant_id = variant.variant_id'.format(table))]
                else:
                    # [] <-> [] connection
                    # outputting non-master variant table field from another non-master variant table
                    # here we link table directly to another variant table by variant id.
                    return [FieldConnection(
                        field='{}.{}'.format(table, field),
                        table=table, # need to join table with another table
                        link='{}.variant_id = {}.variant_id'.format(variant_table, table))]
            # Annotation database
            if table.lower() in [x.name.lower() for x in self.annoDB]:
                db = [x for x in self.annoDB if x.name.lower() == table][0]
                if db.anno_type == 'field':
                    return sum([self.linkFieldToTable(x, variant_table) for x in db.linked_by], []) + [
                        FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= ' AND '.join(['{}.{}={}'.format(table, x, y) for x,y in zip(db.build, db.linked_by)]))]
                if db.build is not None:
                    if db.anno_type == 'position':  # chr and pos
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.bin = {0}.{1}_bin AND variant.chr = {0}.{2} AND variant.pos = {0}.{3}'\
                                .format(table, self.build, db.build[0], db.build[1]))]
                    elif db.anno_type == 'variant':  # chr, pos, alt and alt
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.bin = {0}.{1}_bin AND variant.chr = {0}.{2} AND variant.pos = {0}.{3} AND variant.ref = {0}.{4} AND variant.alt = {0}.{5}'\
                                    .format(table, self.build, db.build[0], db.build[1], db.build[2], db.build[3]))]
                    elif db.anno_type == 'range':  # chr, start, and end
                        binningTable = makeTableName([self.build] + db.build)
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.bin = {0}.bin AND variant.chr = {0}.chr AND {0}.range_id = {1}.rowid'
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.chr = {0}.{1} AND variant.pos >= {0}.{2} AND variant.pos <= {0}.{3}'
                                    .format(table, db.build[0], db.build[1], db.build[2]))
                            ]
                    else:
                        raise ValueError('Unsupported annotation type {}'.format(db.anno_type))
                if db.alt_build is not None:
                    if db.anno_type == 'position':  # chr and pos
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.alt_bin = {0}.{1}_bin AND variant.alt_chr = {0}.{2} AND variant.alt_pos = {0}.{3}'\
                                .format(table, self.alt_build, db.alt_build[0], db.alt_build[1]))]
                    elif db.anno_type == 'variant':  # chr, pos and alt
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.alt_bin = {0}.{1}_bin AND variant.alt_chr = {0}.{2} AND variant.alt_pos = {0}.{3} AND variant.ref = {0}.{4} AND variant.alt = {0}.{5}'\
                                    .format(table, self.alt_build, db.alt_build[0], db.alt_build[1], db.alt_build[2], db.alt_build[3]))]
                    elif db.anno_type == 'range':  # chr, start, and end
                        binningTable = makeTableName([self.alt_build] + db.alt_build)
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.alt_bin = {0}.bin AND variant.alt_chr = {0}.chr AND {0}.range_id = {1}.rowid'
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            # FIXME: how to use bin here?
                            link= 'variant.alt_chr = {0}.{1} AND variant.alt_pos >= {0}.{2} AND variant.alt_pos <= {0}.{3}'\
                                    .format(table, db.alt_build[0], db.alt_build[1], db.alt_build[2]))]
                    else:
                        raise ValueError('Unsupported annotation type {}'.format(db.anno_type))
            raise ValueError('Failed to locate field {}'.format(field))
        # get all fields
        if field.lower() not in ['chr', 'pos', 'ref', 'alt', 'variant_id']:
            matching_fields = []
            for table in self.getVariantTables():
                matching_fields.extend(['{}.{}'.format(table, f) for f in self.db.getHeaders(table) if f.lower() == field.lower()])
            for db in self.annoDB:
                matching_fields.extend(['{}.{}'.format(db.name, x.name) for x in db.fields if x.name.lower() == field.lower()])
            #
            # if no record
            if len(matching_fields) == 0:
                raise ValueError('Failed to locate field {}. Please use command "vtools show fields" to see a list of available fields.'.format(field))
            # if duplicate records
            elif len(matching_fields) > 1:
                raise RuntimeError('There are more than one matching fields {}. Please use table.field to avoid error.'.format(' ,'.join(matching_fields)))
        #
        #
        # in variant_table?
        if field.lower() in [x.lower() for x in self.db.getHeaders(variant_table)]:
            return [FieldConnection(
                field='{}.{}'.format(variant_table, field),
                table='',  # no link is necessary
                link='')]
        # in the master variant table
        if variant_table.lower() != 'variant' and field.lower() in [x.lower() for x in self.db.getHeaders('variant')]:
            return [FieldConnection(
                field='variant.{}'.format(field),
                table='variant',
                link='{}.variant_id = variant.variant_id'.format(variant_table))]
        # other variant tables
        for table in self.getVariantTables():
            # this case has been checked (in variant_table)
            if table.lower() == variant_table.lower():
                continue
            if field.lower() in [x.lower() for x in self.db.getHeaders(table)]:
                # direct link to another variant table
                return [FieldConnection(
                    field= '{}.{}'.format(table, field),
                    table= '{}.{}'.format(table, table),
                    link= '{}.variant_id = {}.variant_id'.format(table, variant_table))]
        # annotation database?
        for db in self.annoDB:
            if field.lower() in [x.name.lower() for x in db.fields]:
                table = db.name
                if db.anno_type == 'field':
                    return sum([self.linkFieldToTable(x, variant_table) for x in db.linked_by], []) + [
                        FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= ' AND '.join(['{}.{}={}'.format(table, x, y) for x,y in zip(db.build, db.linked_by)]))]
                if db.build is not None:
                    if db.anno_type == 'position':  # chr and pos
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.bin = {0}.{1}_bin AND variant.chr = {0}.{2} AND variant.pos = {0}.{3}'\
                                .format(table, self.build, db.build[0], db.build[1]))]
                    elif db.anno_type == 'variant':  # chr, pos and alt
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.bin = {0}.{1}_bin AND variant.chr = {0}.{2} AND variant.pos = {0}.{3} AND variant.ref = {0}.{4} AND variant.alt = {0}.{5}'\
                                    .format(table, self.build, db.build[0], db.build[1], db.build[2], db.build[3]))]
                    elif db.anno_type == 'range':  # chr, start, and end
                        binningTable = makeTableName([self.build] + db.build)
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.bin = {0}.bin AND variant.chr = {0}.chr AND {0}.range_id = {1}.rowid'
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.chr = {0}.{1} AND variant.pos >= {0}.{2} AND variant.pos <= {0}.{3}'\
                                    .format(table, db.build[0], db.build[1], db.build[2]))]
                    else:
                        raise ValueError('Unsupported annotation type {}'.format(db.anno_type))
                elif db.alt_build is not None:
                    if db.anno_type == 'position':  # chr and pos
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.alt_bin = {0}.{1}_bin AND variant.alt_chr = {0}.{2} AND variant.alt_pos = {0}.{3}'\
                                .format(table, self.alt_build, db.alt_build[0], db.alt_build[1]))]
                    elif db.anno_type == 'variant':  # chr, pos and alt
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.alt_bin = {0}.{1}_bin AND variant.chr = {0}.{2} AND variant.pos = {0}.{3} AND variant.ref = {0}.{4} AND alt = {0}.{5}'\
                                    .format(table, self.alt_build, db.alt_build[0], db.alt_build[1], db.alt_build[2], db.alt_build[3]))]
                    elif db.anno_type == 'range':  # chr, start, and end
                        binningTable = makeTableName([self.alt_build] + db.alt_build)
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.alt_bin = {0}.bin AND variant.alt_chr = {0}.chr AND {0}.range_id = {1}.rowid'
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= 'variant.chr = {0}.{1} AND variant.pos >= {0}.{2} AND variant.pos <= {0}.{3}'\
                                    .format(table, db.alt_build[0], db.alt_build[1], db.alt_build[2]))]
                    else:
                        raise ValueError('Unsupported annotation type {}'.format(db.anno_type))
                else:
                    raise ValueError('Database does not define any reference genome.')
        raise ValueError('Failed to locate field {}'.format(field))
        

class ProjCopier:
    def __init__(self, proj, dir, vtable, samples, genotypes):
        self.proj = proj
        self.logger = proj.logger
        self.db = proj.db
        #
        files = glob.glob('{}/*.proj'.format(dir))
        if len(files) != 1:
            raise ValueError('Directory {} does not contain a valid variant tools project'.format(dir))
        self.vtable = vtable
        self.samples = samples
        self.genotypes = genotypes
        self.proj_file = files[0]
        self.geno_file = self.proj_file.replace('.proj', '_genotype.DB')
        if not os.path.isfile(self.geno_file):
            self.geno_file = None

    def copyProject(self):
        self.db.attach(self.proj_file, '__fromDB')
        tables = self.db.tables('__fromDB')
        if self.vtable not in tables:
            raise ValueError('Table {} does not exist in project {}'.format(self.vtable, self.proj_file))
        headers = self.db.getHeaders('__fromDB.{}'.format(self.vtable))
        if headers[0] != 'variant_id':
            raise ValueError('Table {} is not a variant table'.format(args.table[1]))
        #
        prog = ProgressBar('Copying variant tables {}'.format(self.proj_file), len(tables))
        cur = self.db.cursor()
        for idx, table in enumerate(tables):
            # get schema
            cur.execute('SELECT sql FROM __fromDB.sqlite_master WHERE type="table" AND name={};'.format(self.db.PH),
                (table,))
            sql = cur.fetchone()[0]
            if self.db.hasTable(table):
                self.db.removeTable(table)
            try:
                #self.logger.debug(sql)
                cur.execute(sql)
            except Exception as e:
                self.logger.debug(e)
            # copying data over
            if self.proj.isVariantTable(table):
                if self.vtable == 'variant':
                    # copy all variants
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromDB.{0};'.format(table))
                else:
                    cur.execute('''INSERT INTO {0} SELECT * FROM __fromDB.{0} 
                        WHERE __fromDB.{0}.variant_id IN (SELECT variant_id FROM __fromDB.{1});'''\
                            .format(table, self.vtable))
            else:
                cur.execute('INSERT INTO {0} SELECT * FROM __fromDB.{0};'.format(table))
            #
            prog.update(idx)
        prog.done()
        # remove uncopied samples
        self.IDs = self.proj.selectSampleByPhenotype(self.samples)
        # if a condition is specified, remove unselected IDs
        if self.samples:
            allIDs = self.proj.selectSampleByPhenotype('')
            removed = [x for x in allIDs if x not in self.IDs]
            for ID in removed:
                cur.execute('DELETE FROM sample WHERE sample_id = {};'.format(self.db.PH),
                    (ID,))
            self.logger.debug('Removing {} unselected samples'.format(len(removed)))
        self.proj.saveProperty('annoDB', '[]')
        return self.db.numOfRows('variant', False)

    def copySamples(self):
        # copy genotype table
        db = DatabaseEngine()
        db.connect('{}_genotype'.format(self.proj.name))
        db.attach(self.proj_file, '__fromDB')
        db.attach(self.geno_file, '__fromGeno')
        cur = db.cursor()
        prog = ProgressBar('Copying samples', len(self.IDs))
        for idx, ID in enumerate(sorted(self.IDs)):
            table = 'genotype_{}'.format(ID)
            # get schema
            cur.execute('SELECT sql FROM __fromGeno.sqlite_master WHERE type="table" AND name={0};'.format(db.PH),
                (table,))
            sql = cur.fetchone()[0]
            if db.hasTable(table):
                db.removeTable(table)
            try:
                # self.logger.debug(sql)
                cur.execute(sql)
            except Exception as e:
                self.logger.debug(e)
            # copying data over
            if self.vtable == 'variant':
                if self.genotypes:
                    # copy selected genotypes
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0} WHERE {1};'.format(table, self.genotypes))
                else:
                    # copy all variants and genotypes
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0};'.format(table))
            else:
                if self.genotypes:
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0} WHERE ({2}) AND (__fromGeno.{0}.variant_id IN (SELECT variant_id FROM __fromDB.{1}));'.format(table, self.vtable, self.genotypes))
                else:
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0} WHERE __fromGeno.{0}.variant_id IN (SELECT variant_id FROM __fromDB.{1});'.format(table, self.vtable))
            prog.update(idx + 1)
        # remove all annotations
        prog.done()
        db.detach('__fromGeno')
        db.detach('__fromDB')
        db.close()
        return len(self.IDs)

    def createIndex(self, sqls):
        db = DatabaseEngine()
        db.connect(self.proj.proj_file)
        cur = db.cursor()
        for sql in sqls:
            # sql can be None ...
            if sql[0]:
                cur.execute(sql[0])
        db.close()

    def copy(self):
        copied_variants = self.copyProject()
        # create indexes in a separate thread
        cur = self.db.cursor()
        cur.execute('SELECT sql FROM __fromDB.sqlite_master WHERE type="index";')
        sqls = cur.fetchall()
        # close project because it will be opened again in a separate thread
        self.proj.db.detach('__fromDB')
        self.proj.db.close()
        #
        thread = threading.Thread(target=self.createIndex, args=(sqls,))
        thread.start()
        # start copying samples
        if self.geno_file is not None:
            copied_samples = self.copySamples()
        else:
            copied_samples = 0
        self.logger.info('{} variants and {} samples are copied'.format(copied_variants, copied_samples))
        # wait for the thread to close
        s = delayedAction(self.logger.info, 'Create indexes')
        if thread.is_alive():
            thread.join()
        del s
        # re-connect the main database for proer cleanup
        self.proj.db = DatabaseEngine()
        self.proj.db.connect(self.proj.proj_file)


class SortedVariantMapper(threading.Thread):
    '''The worker thread read variants from all projects, sort them, and
       figure out the distinct variants. It uses less RAM than another mapper
       but requires the Unix sort command. A huge disadvantage of this mapper
       if that the read and process threads cannot be run in parallel.'''
    def __init__(self, projects, alt_build, status):
        self.projects = projects
        self.alt_build = alt_build
        self.status = status
        threading.Thread.__init__(self, name='Read, sort and map variants')

    def run(self):
        if self.alt_build:
            psort = Popen(['sort', '-k3,3', '-k4,4n', '-k5,8'], stdin=PIPE, stdout=PIPE)
        else:
            psort = Popen(['sort', '-k3,3', '-k4,4n', '-k5,6'], stdin=PIPE, stdout=PIPE)
        for idx, proj in enumerate(self.projects):
            db = DatabaseEngine()
            db.connect(proj)
            cur = db.cursor()
            prog = ProgressBar('Loading {} ({}/{})'.format(proj[:-5], idx+1, len(self.projects)), db.numOfRows('variant',  False))
            if self.alt_build:
                cur.execute('SELECT variant_id, chr, pos, ref, alt, alt_chr, alt_pos FROM variant;')
                for count, (id, chr, pos, ref, alt, alt_chr, alt_pos) in enumerate(cur):
                    psort.stdin.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.\
                        format(idx, id, chr, pos, alt_chr, alt_pos, ref, alt).encode())
                    if count % 10000 == 0:
                        prog.update(count)
            else:
                cur.execute('SELECT variant_id, chr, pos, ref, alt FROM variant;')
                for count, (id, chr, pos, ref, alt) in enumerate(cur):
                    psort.stdin.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(idx, id, chr, pos, ref, alt).encode())
                    if count % 10000 == 0:
                        prog.update(count)
            prog.done()
            db.close()
        psort.stdin.close()
        #
        prog = ProgressBar('Sorting variants')
        idMaps = {x:{} for x in range(len(self.projects))}
        last_rec = None
        new_id = 1
        for count, line in enumerate(psort.stdout):
            rec = line.decode().rstrip().split('\t', 2)
            # rec[0] is source, rec[1] is source_id
            if last_rec is None or last_rec != rec[2]:
                # a new record
                idMaps[int(rec[0])][int(rec[1])] = (new_id, 0)
                last_rec = rec[2]
                new_id += 1
            else:
                # last_rec[2] == rec[2]
                # we use 1 to mark a duplicated entry
                idMaps[int(rec[0])][int(rec[1])] = (new_id, 1)
            if count % 10000 == 0:
                prog.update(count)
        prog.done()
        #
        prog = ProgressBar('Create mapping tables')
        for idx, proj in enumerate(self.projects):
            db = DatabaseEngine()
            db.connect(proj)
            cur = db.cursor()
            prog.reset('create mapping table {}'.format(idx + 1))
            cur.execute('DROP TABLE IF EXISTS __id_map;')
            cur.execute('CREATE TABLE __id_map (old_id INT, new_id INT, is_dup INT);')
            insert_query = 'INSERT INTO __id_map VALUES ({0}, {0}, {0});'.format(db.PH);
            the_same = True
            keep_all = True
            for old_id, (new_id, is_duplicate) in idMaps[idx].iteritems():
                if the_same and old_id != new_id:
                    the_same = False
                if keep_all and is_duplicate:
                    keep_all = False
                if not the_same and not keep_all:
                    break
            self.status.set(proj, 'the_same', the_same)
            self.status.set(proj, 'keep_all', keep_all)
            cur.executemany(insert_query, ([x, y[0], y[1]] for x, y in idMaps[idx].iteritems()))
            prog.update(count)
            db.commit()
            db.close()
        prog.done()
        # free some RAM
        idMaps.clear()
        # tells the dispatcher that we have created all maps
        for proj in self.projects:
            self.status.set(proj, 'stage', 1)
            self.status.set(proj, 'completed', 1)

class VariantMapper(threading.Thread):
    '''The worker thread read variants from all projects and create 
    id maps along the way. This mapper does not sort variants, uses 
    more RAM, but it is much faster because it can run in paralle..'''
    def __init__(self, projects, alt_build, status):
        self.projects = projects
        self.alt_build = alt_build
        self.status = status
        threading.Thread.__init__(self, name='Read and map variants')

    def run(self):
        existing = {}
        new_id = 1
        keep_all = True
        for idx, proj in enumerate(self.projects):
            db = DatabaseEngine()
            db.connect(proj)
            cur = db.cursor()
            idMaps = {}
            prog = ProgressBar('Loading {} ({}/{})'.format(proj[:-5], idx+1, len(self.projects)), db.numOfRows('variant',  False))
            if self.alt_build:
                count = 0
                cur.execute('SELECT variant_id, chr, pos, ref, alt, alt_chr, alt_pos FROM variant;')
                for id, chr, pos, ref, alt, alt_chr, alt_pos in cur:
                    if alt_chr is None:
                        sub_key = pos
                    elif alt_chr == chr:
                        sub_key = (pos, alt_pos)
                    else:
                        sub_key = (pos, alt_chr, alt_pos)
                    vv = existing.get((chr, ref, alt))
                    if vv is not None:
                        if sub_key in vv:
                            idMaps[id] = (vv[sub_key], 1)
                            keep_all = False
                        else:
                            vv[sub_key] = new_id
                            idMaps[id] = (new_id, 0)
                            new_id += 1
                    else:
                        existing[(chr, ref, alt)] = {sub_key: new_id}
                        idMaps[id] = (new_id, 0)
                        new_id += 1
                    count += 1
                    if count % 10000 == 0:
                        prog.update(count)
            else:
                count = 0
                cur.execute('SELECT variant_id, chr, pos, ref, alt FROM variant;')
                #
                for id, chr, pos, ref, alt in cur:
                    # a new record
                    vv = existing.get((chr, ref, alt))
                    if vv is not None:
                        if pos in vv:
                            idMaps[id] = (vv[pos], 1)
                            keep_all = False
                        else:
                            vv[pos] = new_id
                            idMaps[id] = (new_id, 0)
                            new_id += 1
                    else:
                        existing[(chr, ref, alt)] = {pos: new_id}
                        idMaps[id] = (new_id, 0)
                        new_id += 1
                    count += 1
                    if count % 10000 == 0:
                        prog.update(count)
            #
            prog.reset('mapping ids')
            cur.execute('DROP TABLE IF EXISTS __id_map;')
            cur.execute('CREATE TABLE __id_map (old_id INT, new_id INT, is_dup INT);')
            insert_query = 'INSERT INTO __id_map values ({0}, {0}, {0});'.format(db.PH);
            the_same = True
            for _old_id, (_new_id, _is_duplicate) in idMaps.iteritems():
                if _old_id != _new_id:
                    the_same = False
                    break
            #
            self.status.set(proj, 'the_same', the_same)
            self.status.set(proj, 'keep_all', keep_all)
            # insert into database
            cur.executemany(insert_query, ([x, y[0], y[1]] for x, y in idMaps.iteritems()))
            db.commit()
            #
            idMaps.clear()
            prog.done()
            #
            self.status.set(proj, 'completed', 1)
            self.status.set(proj, 'stage', 1)
            db.close()
        # free all the RAM
        existing.clear()

    
class VariantProcessor(threading.Thread):
    '''The worker thread to merge a project to the master project'''
    def __init__(self, queue, status, logger):
        self.queue = queue
        self.logger = logger
        self.status = status
        threading.Thread.__init__(self, name='cache variants')

    def run(self):
        while True:
            item = self.queue.get()
            if item is None:
                self.queue.task_done()
                break
            # get parameters
            self.src_proj = item
            self.all_the_same = self.status.get(self.src_proj, 'the_same')
            self.keep_all = self.status.get(self.src_proj, 'keep_all')
            if not (self.all_the_same and self.keep_all):
                d, p = os.path.split(self.src_proj)
                if not os.path.isdir(os.path.join(d, 'cache')):
                    os.mkdir(os.path.join(d, 'cache'))
                if not os.path.isdir(os.path.join(d, 'cache')):
                    raise RuntimeError('Cannot locate cache directory of project {}'.format(self.src_proj))
                self.cache_proj = os.path.join(d, 'cache', p)
                if os.path.isfile(self.cache_proj):
                    os.remove(self.cache_proj)
                #
                self.cacheVariants()
            self.status.set(self.src_proj, 'stage', 2)
            self.status.set(self.src_proj, 'completed', 2)
            # get be further processed
            self.status.set(self.src_proj, 'scheduled', False)
            self.queue.task_done()

    def cacheVariants(self):
        if self.all_the_same and self.keep_all:
            return
        db = DatabaseEngine()
        db.connect(self.cache_proj)
        db.attach(self.src_proj, '__fromDB')
        # create index on __fromDB
        cur = db.cursor()
        cur.execute('CREATE INDEX __fromDB.__id_map_idx ON __id_map (old_id ASC);')
        # create tables
        cur.execute('''SELECT sql FROM __fromDB.sqlite_master WHERE type="table" 
            AND name != "sqlite_sequence" AND name != "__id_map";''')
        for sql in cur.fetchall():
            # sql might be None
            if sql[0] is None:
                continue
            cur.execute(sql[0])
        # copy table
        for table in db.tables():
            # we should perhaps allow this so that the cache project become
            # a correct project.
            if table in ['sample', 'filename', 'project']:
                continue
            query = None
            headers = db.getHeaders(table)
            if self.all_the_same:
                # table variant:
                #     if keep_all: do not copy anything
                #     otherwise:   copy some of them (is_dup == 0)
                # 
                # Other variant tables:
                #     if keep_all: do not copy anything
                #     otherwise:   do not copy anything because it should keep all.
                # 
                if table == 'variant' and not self.keep_all:
                    query = '''INSERT INTO {0}  
                                SELECT {1}
                                FROM __fromDB.{0} v LEFT OUTER JOIN __fromDB.__id_map m
                                    ON v.variant_id = m.old_id 
                                WHERE m.is_dup = 0;'''.format(table, ', '.join(headers))
            else:
                # table variant:
                #     if keep_all: copy with id map
                #     otherwise:   copy with id map, and is_dup == 0
                #
                # Other variant table:
                #     if keep_all: copy with id_map
                #     otherwise:   copy with id_map 
                #
                if table != 'variant' or not self.keep_all:
                    query = '''INSERT INTO {1} (variant_id {0}) 
                                SELECT m.new_id {0} 
                                FROM __fromDB.{1} v LEFT OUTER JOIN __fromDB.__id_map m
                                    ON v.variant_id = m.old_id
                                WHERE m.is_dup = 0;'''.format(
                    ' '.join([', {}'.format(x) for x in headers[1:]]), table)
                else:
                    query = '''INSERT INTO {1} (variant_id {0}) 
                                SELECT m.new_id {0} 
                                FROM __fromDB.{1} v LEFT OUTER JOIN __fromDB.__id_map m
                                    ON v.variant_id = m.old_id;'''.format(
                    ' '.join([', {}'.format(x) for x in headers[1:]]), table)
            if query is not None:
                self.logger.debug('Caching table {} of project {} ({})'.format(table, self.src_proj,
                    query))
                cur.execute(query)
        db.detach('__fromDB')
        db.close()

class SampleProcessor(threading.Thread):
    def __init__(self, queue, status, logger):
        self.queue = queue
        self.logger = logger
        self.status = status
        threading.Thread.__init__(self, name='cache samples')
    
    def run(self):
        while True:
            item = self.queue.get()
            if item is None:
                self.queue.task_done()
                break
            # get parameters
            self.src_proj = item
            self.all_the_same = self.status.get(self.src_proj, 'the_same')
            self.samples = self.status.get(self.src_proj, 'old_ids')
            #
            if self.all_the_same:
                self.status.set(self.src_proj, 'completed', 3)
            else:
                self.status.set(self.src_proj, 'total_count', 2 + len(self.samples))
                self.src_geno = self.src_proj.replace('.proj', '_genotype.DB')
                if not os.path.isfile(self.src_geno):
                    self.src_geno = None
                d, p = os.path.split(self.src_proj)
                if not os.path.isdir(os.path.join(d, 'cache')):
                    os.mkdir(os.path.join(d, 'cache'))
                if not os.path.isdir(os.path.join(d, 'cache')):
                    raise RuntimeError('Cannot locate cache directory of project {}'.format(sec_proj))
                self.cache_geno = os.path.join(d, 'cache', p.replace('.proj', '_genotype.DB'))
                if os.path.isfile(self.cache_geno):
                    os.remove(self.cache_geno)
                #
                self.cacheSample()
            self.status.set(self.src_proj, 'stage', 3)
            self.status.set(self.src_proj, 'scheduled', False)
            self.queue.task_done()

    def cacheSample(self):
        # connect to cache geno
        db = DatabaseEngine()
        db.connect(self.cache_geno)
        db.attach(self.src_proj, '__proj')
        db.attach(self.src_geno, '__geno')
        cur = db.cursor()
        #
        for idx, _old_id in enumerate(self.samples):
            # create genotype table
            cur.execute('SELECT sql FROM __geno.sqlite_master WHERE type="table" AND name={0};'.format(db.PH),
                ('genotype_{}'.format(_old_id), ))
            sql = cur.fetchone()
            if sql is None:
                raise ValueError('Cannot recreate genotype table {} from {}, please check the integrity of the database.'.format(_old_id, self.src_geno))
            sql = sql[0]
            try:
                cur.execute(sql)
            except Exception as e:
                self.logger.debug(e)
            #
            # copy data over
            headers = db.getHeaders('genotype_{}'.format(_old_id))
            query = '''INSERT INTO genotype_{0} SELECT new_id {1} FROM __geno.genotype_{2} 
                LEFT JOIN __proj.__id_map ON __id_map.old_id = __geno.genotype_{2}.variant_id;'''\
                .format(_old_id, ''.join([', {}'.format(x) for x in headers[1:]]), _old_id)
            self.logger.debug('Caching sample {} of project {}'.format(_old_id, self.src_proj))
            cur.execute(query)
            db.commit()
            self.status.set(self.src_proj, 'completed', 3 + idx)
        # clean up
        #cur.execute('DROP TABLE IF EXISTS __proj.__id_map;')
        db.detach('__proj')
        db.detach('__geno')
        db.close()


class VariantCopier(threading.Thread):
    def __init__(self, proj, projects, status):
        self.proj = proj
        self.logger = self.proj.logger
        #
        self.projects = projects
        self.status = status
        threading.Thread.__init__(self, name='copy variants')
        
    def run(self):
        db = DatabaseEngine()
        db.connect(self.proj.proj_file)
        for idx, proj in enumerate(self.projects):
            all_the_same = self.status.get(proj, 'the_same')
            keep_all = self.status.get(proj, 'keep_all')
            # get cache_proj
            db.attach(proj, '__fromDB')
            if not (all_the_same and keep_all):
                d, p = os.path.split(proj)
                if not os.path.isdir(os.path.join(d, 'cache')):
                    raise RuntimeError('Cannot locate cache directory of project {}'.format(proj))
                cache_proj = os.path.join(d, 'cache', p)
                if not os.path.isfile(cache_proj):
                    raise RuntimeError('Cannot locate cache project {}'.format(proj))
                db.attach(cache_proj, '__cacheDB')
            #
            # create index on __fromDB
            cur = db.cursor()
            # copy table
            for table in db.tables():
                if table in ['sample', 'filename', 'project', '__id_map']:
                    continue
                # 
                # if ALL_THE_SAME:
                # table variant:
                #     if keep_all: copy from origin
                #     otherwise:   copy from cache
                # 
                # Other variant tables:
                #     if keep_all: copy from origin
                #     otherwise:   copy from origin
                # 
                # NOT ALL_THE_SAME:
                #     copy from cache
                #
                if all_the_same and (table != 'variant' or keep_all):
                    query = '''INSERT INTO {0} SELECT * FROM __fromDB.{0};'''.format(table)
                else:
                    query = '''INSERT INTO {0} SELECT * FROM __cacheDB.{0};'''.format(table)
                self.logger.debug('Copying table {} from project {} ({}, {})'.format(table, proj,
                    all_the_same, keep_all))
                cur.execute(query)
                db.commit()
            #
            db.detach('__fromDB')
            if not (all_the_same and keep_all):
                db.detach('__cacheDB')
            self.status.set('__copyVariants', 'completed', idx + 1)
        # create index
        #
        if self.proj.alt_build is not None:
            db.execute('''CREATE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
            db.execute('''CREATE INDEX variant_alt_index ON variant (alt_bin ASC, alt_chr ASC, alt_pos ASC, ref ASC, alt ASC);''')
        else:
            db.execute('''CREATE UNIQUE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
        db.close()
        self.status.set('__copyVariants', 'completed', len(self.projects) + 10)

class SampleCopier(threading.Thread):
    def __init__(self, proj, projects, status):
        self.proj = proj
        self.logger = proj.logger
        #
        self.projects = projects
        self.status = status
        #
        threading.Thread.__init__(self, name='copy samples')

    def run(self):
        # connect to cache geno
        db = DatabaseEngine()
        db.connect(self.proj.proj_file.replace('.proj', '_genotype.DB'))
        count = 0
        for proj in self.projects:
            all_the_same = self.status.get(proj, 'the_same')
            if all_the_same:
                src_geno = proj.replace('.proj', '_genotype.DB')
                db.attach(src_geno, '__geno')
            else:
                d, p = os.path.split(proj)
                if not os.path.isdir(os.path.join(d, 'cache')):
                    raise RuntimeError('Cannot locate cache directory of project {}'.format(proj))
                cache_geno = os.path.join(d, 'cache', p.replace('.proj', '_genotype.DB'))
                db.attach(cache_geno, '__geno')
            #
            cur = db.cursor()
            for _old_id, _new_id in zip(self.status.get(proj, 'old_ids'), self.status.get(proj, 'new_ids')):
                # create table
                cur.execute('SELECT sql FROM __geno.sqlite_master WHERE type="table" AND name={};'.format(db.PH),
                    ('genotype_{}'.format(_old_id),))
                sql = cur.fetchone()[0].replace('genotype_{}'.format(_old_id), 'genotype_{}'.format(_new_id))
                cur.execute(sql)
                query = 'INSERT INTO genotype_{} SELECT * FROM __geno.genotype_{};'\
                    .format(_new_id, _old_id)
                self.logger.debug('Copying sample {} from project {}: {}'.format(_old_id, proj, query))
                cur.execute(query)
                count += 1
                self.status.set('__copySamples', 'completed', count)
            db.detach('__geno')
        db.close()

class MergeStatus:
    def __init__(self):
        self.tasks = {}
        self.lock = threading.Lock()

    def add(self, name, items):
        self.tasks[name] = items

    def set(self, proj, item, value):
        self.lock.acquire()
        self.tasks[proj][item] = value
        self.lock.release()

    def get(self, proj, item):
        return self.tasks[proj][item]

    def canProcessVariant(self, proj):
        return (not self.tasks[proj]['scheduled']) and self.tasks[proj]['stage'] == 1

    def canProcessSample(self, proj):
        return (not self.tasks[proj]['scheduled']) and self.tasks[proj]['stage'] == 2

    def canCopyVariants(self):
        return (not self.tasks['__copyVariants']['scheduled']) and \
            all([y['stage'] >=2 for x,y in self.tasks.iteritems() if x not in ['__copyVariants', '__copySamples']])

    def canCopySamples(self):
        return (not self.tasks['__copySamples']['scheduled']) and \
            all([y['stage'] >=3 for x,y  in self.tasks.iteritems() if x not in ['__copyVariants', '__copySamples']])

    def count(self):
        return sum([x['completed'] for x in self.tasks.values()])
    
    def total_count(self):
        return sum([x['total_count'] for x in self.tasks.values()])

    def report(self):
        for key in self.tasks:
            print key, self.tasks[key]['completed'],
        print
        #print self.tasks

class ProjectsMerger:
    def __init__(self, proj, dirs, sort, jobs):
        '''Check if merge is possible, set primary and reference genome
            and set self.projects and self.structure '''
        self.proj = proj
        self.db = proj.db
        self.logger = proj.logger
        self.sort = sort
        self.jobs = jobs
        # valid projects
        self.projects = []
        # check if all subprojects have the same structure
        structure = {}
        # use the largest project as the first one
        proj_size = 0
        for idx, dir in enumerate(dirs):
            files = glob.glob('{}/*.proj'.format(dir))
            if len(files) != 1:
                raise ValueError('Directory {} does not contain a valid variant tools project'.format(dir))
            proj_file = files[0]
            if proj_file in [x[0] for x in self.projects]:
                self.logger.warning('Remove duplicate merge {}'.format(proj_file))
            #
            self.db.attach(files[0], '__fromDB')
            # 
            # get primary and alternative reference genome
            cur = self.db.cursor()
            # primary reference genome
            cur.execute('SELECT value FROM __fromDB.project WHERE name={}'.format(self.db.PH),
                ('build',))
            build = cur.fetchone()
            # alternative reference genome
            cur.execute('SELECT value FROM __fromDB.project WHERE name={}'.format(self.db.PH),
                ('alt_build',))
            alt_build = cur.fetchone()
            #
            if build is None or alt_build is None:
                self.logger.warning('Ignoring invalid or empty project.'.format(proj_file))
                continue
            #
            # if merging from an empty project, use the primary and alterantive reference genome of the first one
            if idx == 0:
                self.proj.build = build[0]
                self.proj.saveProperty('build', self.proj.build)
                self.proj.alt_build = alt_build[0]
                self.proj.saveProperty('alt_build', self.proj.alt_build)
                #
                if self.proj.alt_build is not None:
                    s = delayedAction(self.logger.info, 'Adding alternative reference genome {} to the project.'.format(self.proj.alt_build))
                    cur = self.db.cursor()
                    headers = self.db.getHeaders('variant')
                    for fldName, fldType in [('alt_bin', 'INT'), ('alt_chr', 'VARCHAR(20)'), ('alt_pos', 'INT')]:
                        if fldName in headers:
                            continue
                        self.db.execute('ALTER TABLE variant ADD {} {} NULL;'.format(fldName, fldType))
                    del s
            elif build[0] != self.proj.build:
                raise ValueError('Primary reference genome of project ({} of {}) does not match that of the current project ({}).'\
                    .format(build[0], proj_file, self.proj.build))
            elif alt_build[0] != self.proj.alt_build:
                raise ValueError('Alternative reference genome of project ({} of {}) does not match that of the current project ({})'\
                    .format(alt_build[0], proj_file, self.proj.alt_build))
            #
            # analyze and create project tables
            if idx == 0:
                tables = self.db.tables('__fromDB')
                cur.execute('DROP TABLE variant;')
                cur.execute('DROP TABLE sample;')
                for db_table in tables:
                    table = db_table.split('.')[-1]
                    if table.startswith('__'):
                        continue
                    structure[table] = self.db.fieldsOfTable('__fromDB.{}'.format(table))
                    if table in ['project', 'filename']:
                        continue
                    cur.execute('SELECT sql FROM __fromDB.sqlite_master WHERE type="table" AND name={0};'.format(self.db.PH),
                            (table, ))
                    sql = cur.fetchone()
                    cur.execute(sql[0])
            else:
                tables = [x for x in self.db.tables('__fromDB') if not x.startswith('__')]
                tables.sort()
                if tables != sorted(list(structure.keys())):
                    raise ValueError("Project {} does not have the same set of variant tables ({}) as other projects ({})"\
                        .format(proj_file,
                        ', '.join([x for x in tables if x not in ['sample', 'filename', 'project']]),
                        ', '.join([x for x in structure.keys() if x not in ['sample', 'filename', 'project']])))
                for table in tables:
                    if structure[table] != self.db.fieldsOfTable('__fromDB.{}'.format(table)):
                        raise ValueError('Columns of table {} ({}) in project {} does not match those are in others ({}).'\
                            .format(table, ', '.join(self.db.getHeaders('__fromDB.{}'.format(table))),
                            proj_file, ', '.join([x[0] for x in structure[table]])))
            # we put the largest project the first to improve efficiency, because the
            # first project is effectively copied instead of merged.
            if self.db.numOfRows('__fromDB.variant', False) * self.db.numOfRows('__fromDB.sample') > proj_size:
                self.projects.insert(0, proj_file)
                # we do not need an exact number
                proj_size = self.db.numOfRows('__fromDB.variant', False) * self.db.numOfRows('__fromDB.sample')
            else:
                self.projects.append(proj_file)
            #
            self.db.detach('__fromDB')
       
    def mapSamples(self, status):
        '''Population filename and sample table, return id maps
        '''
        cur = self.db.cursor()
        filenames = []
        #
        ignored_files = 0
        for proj in self.projects:
            self.db.attach(proj, '__proj')
            self.db.attach(proj.replace('.proj', '_genotype.DB'), '__geno')
            #
            # handling filename
            cur.execute('SELECT file_id, filename, header FROM __proj.filename;')
            filename_records = cur.fetchall()
            old_sample_id = []
            new_sample_id = []
            for old_file_id, filename, header in filename_records:
                if filename in filenames:
                    ignored_files += 1
                    continue
                else:
                    filenames.append(filename)
                #
                cur.execute('INSERT INTO filename (filename, header) VALUES ({0}, {0});'.format(self.db.PH),
                    (filename, header))
                new_file_id = cur.lastrowid
                # get samples
                cur.execute('SELECT sample_id FROM __proj.sample WHERE file_id={};'.format(self.db.PH),
                    (old_file_id, ))
                old_sid = [x[0] for x in cur.fetchall()]
                new_sid = []
                headers = self.db.getHeaders('sample')
                for sid in old_sid:
                    cur.execute('''INSERT INTO sample ('file_id', {0}) SELECT {1}, {0} 
                        FROM __proj.sample WHERE __proj.sample.sample_id = {1};'''\
                        .format(', '.join(headers[2:]), self.db.PH), (new_file_id, sid))
                    new_sid.append(cur.lastrowid)
                #
                old_sample_id.extend(old_sid)
                new_sample_id.extend(new_sid)
            status.set(proj, 'old_ids', old_sample_id)
            status.set(proj, 'new_ids', new_sample_id)
            self.db.detach('__proj')
            self.db.detach('__geno')
        return ignored_files
  
    def merge(self):
        if len(self.projects) == 0:
            return
        nJobs = max(min(self.jobs, len(self.projects) + 1), 1)
        #
        #
        # collect status
        status = MergeStatus()
        #
        # for projects, total_count means
        # 1: read
        # 2. variant processed
        # 3. sample processed
        for proj in self.projects:
            status.add(proj, {'completed': 0, 'stage': 0, 'scheduled': False, 'total_count': 3})
        # this will set for each project
        #  old_ids: sample id of the original projects
        #  new_ids: sample id in the new project
        ignored_files = self.mapSamples(status)
        # stop the database so that it can be opened in thread
        self.proj.db.close()
        #
        # two more tasks
        #
        # copy variants: performed after all variants are processed
        # the number of count is the number of projects + 10, which is the estimated
        # amount of work to build indexes
        status.add('__copyVariants', {'completed': 0, 'scheduled': False,
            'total_count': len(self.projects) + 10})
        # copy samples: performaned after all samples are processed,
        # the number of count is the number of samples
        status.add('__copySamples', {'completed': 0, 'scheduled': False,
            'total_count': sum([len(status.get(proj, 'old_ids')) for proj in self.projects])})
        # start all variant cachers
        self.vcQueue = Queue.Queue()
        for j in range(nJobs):
            VariantProcessor(self.vcQueue, status, self.logger).start()
        #
        # start all sample cachers
        self.scQueue = Queue.Queue()
        for j in range(nJobs):
            SampleProcessor(self.scQueue, status, self.logger).start()
        #
        # create a thread to read variants
        if self.sort:
            SortedVariantMapper(self.projects, self.proj.alt_build, status).start()
        else:
            VariantMapper(self.projects, self.proj.alt_build, status).start()
        #
        prog = None
        count = 0
        while True:
            for idx, proj in enumerate(self.projects):
                if status.canProcessVariant(proj) and self.vcQueue.qsize() < nJobs:
                    self.logger.debug('Process variant in {}'.format(proj))
                    status.set(proj, 'scheduled', True)
                    self.vcQueue.put(proj)
                if status.canProcessSample(proj) and self.vcQueue.qsize() + self.scQueue.qsize() < nJobs:
                    self.logger.debug('Process variant in {}'.format(proj))
                    status.set(proj, 'scheduled', True)
                    self.scQueue.put(proj)
            if status.canCopyVariants():
                status.set('__copyVariants', 'scheduled', True)
                VariantCopier(self.proj, self.projects, status).start()
                # stop all variant cachers
                for j in range(nJobs):
                    self.vcQueue.put(None)
                self.vcQueue.join()
                prog = ProgressBar('Merging all projects', status.total_count())
            if status.canCopySamples():
                status.set('__copySamples', 'scheduled', True)
                SampleCopier(self.proj, self.projects, status).start()
                # stop all sample cachers
                for j in range(nJobs):
                    self.scQueue.put(None)
                self.scQueue.join()
            #
            if status.count() > count and prog:
                count = status.count()
                prog.update(count)
            if status.count() == status.total_count():
                break
            time.sleep(1)
        #
        #for t in threading.enumerate():
        #    print 'Waiting for {}'.format(t.getName())
        prog.done()
        # re-hook up the project database for proper cleanup
        self.proj.db = DatabaseEngine()
        self.proj.db.connect(self.proj.proj_file)

        
#
#
# Functions provided by this script
#
#

def initArguments(parser):
    parser.add_argument('project',
        help='''Name of a new project. This will create a new file $project.proj under
            the current directory. Only one project is allowed in a directory.''')
    parent = parser.add_argument_group('Derive from a parent project')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Remove a project if it already exists.''')
    parent.add_argument('--parent', metavar='DIR',
        help='''Directory of a parent project (e.g. --parent ../main) from which all or part
            of variants (--variants), samples (--samples) and genotypes (--genotypes) will
            be copied to the newly created project.'''),
    parent.add_argument('--variants', nargs='?', metavar='TABLE', default='variant',
        help='''A variant table of the parental project whose variants will be copied to
            the new project. Default to variant (all variants).'''),
    parent.add_argument('--samples', nargs='*', metavar='COND', default=[],
        help='''Copy only samples of the parental project that match specified conditions.''')
    parent.add_argument('--genotypes', nargs='*', metavar='COND', default=[],
        help='''Copy only genotypes that match specified conditions.''')
    children = parser.add_argument_group('Merge from children projects')
    children.add_argument('--children', nargs='+', metavar='DIR',
        help='''A list of a subprojects (directories) that will be merged to create
            this new project. The subprojects must have the same structure (primary
            and alternative reference genome, variant info and phenotype). Samples
            imported from the same files will be ignored''')
    #sub.add_argument('--sort', action='store_true',
    #    help='''Sort variants read from subprojects, which takes less RAM but longer time. If
    #        unset, all variants will be read to RAM and perform a faster merge at a cost of
    #        higher memory usage'''),
    #sub.add_argument('-j', '--jobs', metavar='N', type=int, default=4,
    #    help='''Number of threads used to merge subprojects. Default to 4.'''),
    #parser.add_argument('--build',
    #    help='''Build of the primary reference genome of the project.'''),
    #engine = parser.add_argument_group('Database connection')
    #engine.add_argument('--engine', choices=['mysql', 'sqlite3'], default='sqlite3',
    #    help='''Database engine to use, can be mysql or sqlite3. Parameters --host, --user
    #        and --passwd will be needed for the creation of a new mysql project.''')
    #engine.add_argument('--host', default='localhost', 
    #    help='The MySQL server that hosts the project databases.')
    #engine.add_argument('--user', default=getpass.getuser(),
    #    help='User name to the MySQL server. Default to current username.')
    #engine.add_argument('--passwd',
    #    help='Password to the MySQL server.')
    #parser.add_argument('--batch', default=10000, 
    #    help='Number of query per transaction. Larger number leads to better performance but requires more ram.')


def init(args):
    try:
        if args.force:
            # silently remove all exiting project.
            try:
                proj = Project(verbosity='0', verify=False)
                proj.remove()
            except:
                pass
        # create a new project
        #
        # args.batch is temporarily removed to keep interface clean
        # args.build is temporarily removed to keep interface clean
        #
        # args.engine, host, user, passwd is disabled because MySQL support is
        # currently lagging behind.
        with Project(name=args.project, build=None, new=True, 
            verbosity='1' if args.verbosity is None else args.verbosity,
            engine='sqlite3', host=None, user=None, passwd=None,
            batch=10000) as proj:
            if args.parent:
                copier = ProjCopier(proj, args.parent, args.variants,
                    ' AND '.join(['({})'.format(x) for x in args.samples]),
                    ' AND '.join(['({})'.format(x) for x in args.genotypes]))
                copier.copy()
            elif args.children:
                # args.sort is temporarily removed to keep interface clean
                # because the only advantage (save RAM) does not really justify
                # its inclusion.
                #
                # A default value 4 is given for args.jobs because more threads usually
                # do not improve effiency
                if args.unknown_args:
                    raise ValueError('Unknown args for option --children: {}'.format(' '.join(args.unknown_args)))
                merger = ProjectsMerger(proj, args.children, False, 4)
                merger.merge()
    except Exception as e:
        sys.exit(e)


def mergeArguments(parser):
    parser.add_argument('projects', nargs='+', metavar='DIR',
        help='''Directories to a list of projects that will be merged''')

def merge(args):
    try:
        # create a new project
        with Project(verbosity=args.verbosity) as proj:
            proj.merge(args.projects)
    except Exception as e:
        sys.exit(e)

def removeArguments(parser):
    parser.add_argument('type', choices=['project', 'tables', 'samples', 'fields', 'geno_fields', 'annotations', 'variants', 'genotypes', 'phenotypes'],
        help='''Type of items to be removed.''')
    parser.add_argument('items', nargs='*',
        help='''Items to be removed, which should be, for 'project' the name of project to be
            removed (optional), for 'tables' names of one or more variant tables (wildcard characters
            ? and * are allowed), for 'samples' patterns using which matching samples are removed, for
            'fields' name of fields to be removed, for 'geno_fields' name of genotype fields to 
            be removed (cf. 'vtools show genotypes'), for 'annotations' names of annotation databases,
            for 'variants' variant tables whose variants will be removed from all variant
            tables and genotypes, for 'genotypes' conditions using which matching genotypes
            are removed, and for 'phenotypes' columns in the output of 'vtools show samples'.
            Note that removal of samples will only remove sample name, filename (if all related
            samples are removed), and related genotypes, but not variants themselves; removal
            of annotation databases will stop using these databases in the project, but will
            not delete them from disk.''')

def remove(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if args.type == 'project':
                if len(args.items) > 0 and args.items[0] != proj.name:
                    raise ValueError('Cannot remove project: Incorrect project name')
                proj.remove()
            elif args.type == 'tables':
                allTables = proj.getVariantTables()
                removed = []
                for table in args.items:
                    if '?' in table or '*' in table:
                        for tbl in allTables:
                            if re.match(table.replace('?', '.{1}').replace('*', '.*'), tbl, re.I) and tbl not in removed:
                                proj.removeVariantTable(tbl)
                                removed.append(table)
                    else:
                        proj.removeVariantTable(table)
                        removed.append(table)
            elif args.type == 'samples':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select samples to be removed')
                proj.db.attach(proj.name + '_genotype')
                IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.items]))
                if len(IDs) == 0:
                    proj.logger.warning('No sample is selected by condition {}'.format(' AND '.join(['({})'.format(x) for x in args.items])))
                proj.removeSamples(IDs)
            elif args.type == 'fields':
                from_table = defaultdict(list)
                for item in args.items:
                    if item.lower() in ['variant_id', 'chr', 'pos', 'alt']:
                        raise ValueError('Fields variant_id, chr, pos and alt cannot be removed')
                    found = False
                    for table in proj.getVariantTables():
                        if item.lower() in [x.lower() for x in proj.db.getHeaders(table)]:
                            from_table[table].append(item)
                            found = True
                            break
                    if not found:
                        raise ValueError('Field {} does not exist in any of the variant tables.'.format(item))
                # remove...
                for table, items in from_table.items():
                    proj.logger.info('Removing field {} from variant table {}'.format(', '.join(items), table))
                    proj.db.removeFields(table, items)
                    if 'alt_bin' in items or 'alt_chr' in items or 'alt_pos' in items:
                        proj.logger.info('Removing alternative reference genome because of removal of related fields')
                        proj.alt_build = None
                        proj.saveProperty('alt_build', None)
                # it is possible that new indexes are needed
                s = delayedAction(proj.logger.info, 'Rebuilding indexes')
                try:
                    proj.createIndexOnMasterVariantTable()
                except:
                    pass
                del s
            elif args.type == 'geno_fields':
                if len(args.items) == 0:
                    raise ValueError('Please specify name of genotype fields to be removed')
                if 'variant_id' in [x.lower() for x in args.items]:
                    proj.logger.warning('Genotype id variant_id cannot be removed')
                if 'gt' in [x.lower() for x in args.items]:
                    proj.logger.warning('Genotype field GT cannot be removed')
                proj.db.attach(proj.name + '_genotype')
                cur = proj.db.cursor()
                cur.execute('SELECT sample_id FROM sample;')
                IDs = [x[0] for x in cur.fetchall()]
                for table in ['{}_genotype.genotype_{}'.format(proj.name, id) for id in IDs]:
                    header = [x.lower() for x in proj.db.getHeaders(table)]
                    items = [x for x in args.items if x.lower() in header and x.lower not in ['variant_id', 'gt']]
                    if items:
                        proj.logger.info('Removing fields {} from genotype table {}'.format(', '.join(items), table.split('_')[-1]))
                        proj.db.removeFields(table, items)
            elif args.type == 'annotations':
                for item in args.items:
                    removed = False
                    for i in range(len(proj.annoDB)):
                        if proj.annoDB[i].name == item:
                            proj.logger.info('Removing annotation database {} from the project'.format(item))
                            proj.annoDB.pop(i)
                            removed = True
                            break
                    if not removed:
                        proj.logger.warning('Cannot remove annotation database {} from the project'.format(item))
                proj.saveProperty('annoDB', str([os.path.join(x.dir, x.filename) for x in proj.annoDB]))
            elif args.type == 'variants':
                if len(args.items) == 0:
                    raise ValueError('Please specify variant tables that contain variants to be removed')
                proj.db.attach(proj.name + '_genotype')
                for table in args.items:
                    proj.removeVariants(table)
            elif args.type == 'genotypes':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select genotypes to be removed')
                proj.db.attach(proj.name + '_genotype')
                proj.removeGenotypes(' AND '.join(['({})'.format(x) for x in args.items]))
            elif args.type == 'phenotypes':
                if len(args.items) == 0:
                    raise ValueError('Please specify one or more phenotypes (columns in the output of "vtools show samples") to be removed')
                phenos = [x.lower() for x in proj.db.getHeaders('sample')]
                toBeRemoved = []
                for item in args.items:
                    if item.lower() in ['filename', 'sample_name']:
                        raise ValueError('Cannot remove filename or sample_name from the sample table')
                    if item.lower() not in phenos:
                        raise ValueError('No phenotype {} exists in the sample table.'.format(item))
                    toBeRemoved.append(item)
                if len(toBeRemoved) == 0:
                    raise ValueError('No valid phenotype to be removed.')
                # remove
                proj.db.removeFields('sample', toBeRemoved)
    except Exception as e:
        sys.exit(e)


def showArguments(parser):
    parser.add_argument('type', choices=['project', 'tables', 'table',
        'samples', 'genotypes', 'fields', 'annotations', 'annotation', 'formats', 'format',
        'tests', 'test'], nargs='?', default='project',
        help='''Type of information to display, which can be 'project' for
            summary of a project, 'tables' for all variant tables (or all
            tables if --verbosity=2), 'table TBL' for details of a specific
            table TBL, 'samples' for a list of sample and associated phenotypes,
            'fields' for fields from variant tables and all used annotation
            databases, 'annotations' for a list of all available annotation
            databases, 'annotation ANN' for details about annotation database ANN,
            'formats' for all supported import and export formats, 'format FMT' for
            details of format FMT, 'tests' for a list of all association tests, and
            'test TST' for details of an associationt est TST. The default parameter
            of this command is 'project'.''')
    parser.add_argument('items', nargs='*',
        help='''Items to display, which can be names of tables for type 'table',
            name of an annotation database for type 'annotation', name of a format
            for type 'format', and name of an association test for type 'test'.''')
    parser.add_argument('-l', '--limit', default=10, type=int,
        help='''Number of record to display for option 'show table'. All
            records will be displayed if it is set to -1''')


def show(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            #
            limit_clause = ' LIMIT 0, {}'.format(args.limit) if args.limit > 0 else ''
            omitted = '({} records omitted, use parameter --limit to see more)'
            if args.type == 'project':
                print(proj.summarize())
            elif args.type == 'tables':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show tables"'.format(', '.join(args.items)))
                print('{:<20} {}'.format('table', '#variants'))
                for table in proj.db.tables() if args.verbosity=='2' else proj.getVariantTables():
                    print('{:<20} {:,}'.format(table, proj.db.numOfRows(table)))
            elif args.type == 'table':
                proj.db.attach('{}_genotype'.format(proj.name))
                if not args.items:
                    raise ValueError('Please specify a table to display')
                table = args.items[0]
                # showing an annotation database
                if table in [x.name for x in proj.annoDB]:
                    table = '{0}.{0}'.format(table)
                if not proj.db.hasTable(table):
                    if table.startswith('genotype_') and proj.db.hasTable('{}_genotype.{}'.format(proj.name, table)):
                        table = '{}_genotype.{}'.format(proj.name, table)
                    else:
                        raise ValueError('Table {} does not exist'.format(table))
                # print content of table
                headers = proj.db.getHeaders(table)
                print(', '.join(headers))
                cur = proj.db.cursor()
                cur.execute('SELECT * FROM {} {};'.format(table, limit_clause))
                for rec in cur:
                    print(', '.join([str(x) for x in rec]))
                nAll = proj.db.numOfRows(table)
                if args.limit > 0 and args.limit < nAll:
                    print omitted.format(nAll - args.limit)
            elif args.type == 'samples':
                if not proj.db.hasTable('sample'):
                    proj.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                fields = proj.db.getHeaders('sample')
                # headers are ID, file, sample, FIELDS
                print('sample_name\tfilename{}'.format(''.join(['\t'+x for x in fields[3:]])))
                cur.execute('SELECT sample_name, filename {} FROM sample, filename WHERE sample.file_id = filename.file_id ORDER BY sample_name {};'\
                    .format(' '.join([','+x for x in fields[3:]]), limit_clause))
                for rec in cur:
                    print('\t'.join(['{}'.format(x) for x in rec]))
                nAll = proj.db.numOfRows('sample')
                if args.limit > 0 and args.limit < nAll:
                    print omitted.format(nAll - args.limit)
            elif args.type == 'fields':
                if len(proj.annoDB) == 0:
                    proj.logger.info('No annotation database is attached.')
                for table in proj.getVariantTables():
                    print('\n'.join(['{}.{}'.format(table, field) for field in proj.db.getHeaders(table)]))
                for db in proj.annoDB:
                    if args.verbosity == '0':
                        print('\n'.join(['{}.{}'.format(db.name, x.name) for x in db.fields]))
                    else:
                        print('\n'.join(['{}.{} {}'.format(db.name, x.name,
                            '\n'.join(textwrap.wrap(x.comment, initial_indent=' '*(27-len(db.name)-len(x.name)),
                                subsequent_indent=' '*29))) for x in db.fields]))
            elif args.type == 'annotation':
                for item in args.items:
                    try:
                        annoDB = [x for x in proj.annoDB if x.name.lower() == item.lower()][0]
                    except Exception as e:
                        proj.logger.debug(e)
                        raise IndexError('Database {} is not currently used in the project'.format(item))
                    annoDB.describe(args.verbosity == '2')
            elif args.type == 'annotations':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show annotations"'.format(', '.join(args.items)))
                DBs = filesInURL('http://vtools.houstonbioinformatics.org/annoDB', ext='.ann')
                for db in DBs:
                    print(db)
            elif args.type == 'formats':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show formats"'.format(', '.join(args.items)))
                FMTs = filesInURL('http://vtools.houstonbioinformatics.org/format', ext='.fmt')
                for fmt in FMTs:
                    print(fmt)
            elif args.type == 'format':
                if not args.items:
                    raise ValueError('Please specify a format to display')
                for item in args.items:
                    try:
                        fmt = fileFMT(item)
                    except Exception as e:
                        proj.logger.debug(e)
                        raise IndexError('Unrecognized input format: {}\nPlease check your input parameters or configuration file "{}"'.format(e, item))
                    fmt.describe()
            elif args.type == 'genotypes':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show genotypes"'.format(', '.join(args.items)))
                # get sample ids and attach the genotypes database
                if not proj.db.hasTable('sample'):
                    proj.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                try:
                    proj.db.attach(proj.name + '_genotype') 
                except Exception as e:
                    # either the database doesn't exist or it's already been attached (should we add a method proj.db.isAttached(...)? and embedding this in the attach() method?
                    if not proj.db.hasDatabase(proj.name + '_genotype'):
                        proj.logger.debug('Trying to attach a database that doesn\'t exist' + e)
                # sample headers are ID, file, sample, FIELDS
                print('sample_name\tfilename\tnum_genotypes\tsample_genotype_fields')
                cur.execute('SELECT sample.sample_id, sample_name, filename FROM sample, filename WHERE sample.file_id = filename.file_id ORDER BY sample.sample_id {};'.format(limit_clause))
                records = cur.fetchall()
                for rec in records:
                    # sample fields
                    sampleFields = '\t'.join(['{}'.format(x) for x in rec[1:]])
                    # now get sample genotype counts and sample specific fields
                    sampleId = rec[0]
                    cur.execute('SELECT count(*) FROM {}_genotype.genotype_{};'.format(proj.name, sampleId))
                    numGenotypes = cur.fetchone()[0]
                    # get fields for each genotype table
                    sampleGenotypeHeader = proj.db.getHeaders('{}_genotype.genotype_{}'.format(proj.name, sampleId))
                    sampleGenotypeFields = ','.join(['{}'.format(x) for x in sampleGenotypeHeader[1:]])  # the first field is variant id, second is GT
                    print('{}\t{}\t{}'.format(sampleFields, numGenotypes, sampleGenotypeFields))
                nAll = proj.db.numOfRows('sample')
                if args.limit > 0 and args.limit < nAll:
                    print omitted.format(nAll - args.limit)
            elif args.type == 'tests':
                # it is very bad idea to use circular import, but I have no choice
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show tests"'.format(', '.join(args.items)))
                from .association import getAllTests
                print('\n'.join(['{} {}'.format(test,
                    '\n'.join(textwrap.wrap(obj.__doc__, initial_indent=' '*(27-len(test)),
                        subsequent_indent=' '*28))) for test, obj in getAllTests()]))
            elif args.type == 'test':
                from .association import getAllTests
                if len(args.items) == 0:
                    raise ValueError('Please specify the name of a test')
                if len(args.items) > 1:
                    raise ValueError('Please specify only one test')
                tests = getAllTests()
                if args.items[0].lower() not in [x[0].lower() for x in tests]:
                    raise ValueError('Unrecognized test name {}. A list of tests can be obtained from command "vtools show tests"'.format(args.items[0]))
                # test
                test = [y for x,y in tests if x.lower() == args.items[0].lower()][0]
                print 'Name:          {}'.format(args.items[0])
                print 'Description:   {}'.format('\n'.join(textwrap.wrap(test.__doc__, initial_indent='',
                        subsequent_indent=' '*15)))
                # create an instance of the test and pass -h to it
                test(None, ['-h']) 
    except Exception as e:
        sys.exit(e)


def executeArguments(parser):
    parser.add_argument('query', nargs='*',
        help='''A SQL query to be executed. The project genotype database is
        attached as $name_genotype. Annotation databases used in the project
        are attached and are available by their names.''')
    parser.add_argument('-d', '--delimiter', default='\t',
        help='Delimiter used to output results, default to tab.')

def execute(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            proj.db.attach('{}_genotype'.format(proj.name))
            cur = proj.db.cursor()
            query = ' '.join(args.query)
            if query.upper().startswith('SELECT'):
                proj.logger.debug('Analyze statement: "{}"'.format(query))
                cur.execute('EXPLAIN QUERY PLAN ' + query)
                for rec in cur:
                    proj.logger.debug('\t'.join([str(x) for x in rec]))
            # really execute the query
            cur.execute(query)
            proj.db.commit()
            sep = args.delimiter
            for rec in cur:
                print(sep.join(['{}'.format(x) for x in rec]))
    except Exception as e:
        sys.exit(e)


if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools.py
    pass

