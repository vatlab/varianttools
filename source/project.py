#!/usr/bin/env python
#
# $File: project.py $
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
import glob
import logging
import getpass
import textwrap
import tempfile
import shutil
import ConfigParser
import argparse
import threading
import Queue
import signal
import time
import re
import tarfile
import urllib
if sys.version_info.major == 2:
    from ucsctools_py2 import showTrack
else:
    from ucsctools_py3 import showTrack

from multiprocessing import Process
from subprocess import Popen, PIPE
from collections import namedtuple, defaultdict
from .__init__ import VTOOLS_VERSION, VTOOLS_FULL_VERSION, VTOOLS_COPYRIGHT, \
    VTOOLS_CITATION, VTOOLS_REVISION, VTOOLS_CONTACT
from .utils import DatabaseEngine, ProgressBar, SQL_KEYWORDS, delayedAction, \
    RefGenome, filesInURL, downloadFile, getMaxUcscBin, env, sizeExpr, \
    getSnapshotInfo, ResourceManager, decodeTableName, encodeTableName


# define a field type
Field = namedtuple('Field', ['name', 'index', 'adj', 'fmt', 'type', 'comment'])
Column = namedtuple('Column', ['index', 'field', 'adj', 'comment'])
#
# see http://varianttools.sourceforge.net/Calling/New for details
#
PipelineCommand = namedtuple('PipelineCommand', ['index', 'input',
    'input_emitter', 'action', 'comment'])
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
        env.logger.debug('Loading annotation database {}'.format(annoDB))
        self.db = proj.db.newConnection()
        if self.db.hasDatabase(annoDB):
            self.db.connect(annoDB)
        else:
            raise ValueError("Cannot locate annotation database {}".format(annoDB))
        # annoDB can be ~/.variant_tools/annoDB etc, need to expand ~
        self.filename = os.path.expanduser(annoDB)
        # when we save dir, we would better save ~ because $HOME might change ...
        self.dir = os.path.split(annoDB)[0]
        if self.dir.startswith(os.path.expanduser('~')):
            self.dir = self.dir.replace(os.path.expanduser('~'), '~', 1)
        self.filename = os.path.split(annoDB)[-1]
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
        cur.execute('SELECT name, field, "", "", type, comment from {}_field;'.format(self.name))
        for rec in cur:
            self.fields.append(Field(*rec))
            # FIXME: We should enforce comment for all fields.
            #if not self.fields[-1].comment:
            #    env.logger.warning('No comment for field {} in database {}'.format(self.fields[-1].name, annoDB))
        if len(self.fields) == 0:
            raise ValueError('Annotation database {} does not provide any field.'.format(annoDB))
        #
        self.anno_type = 'variant'
        self.linked_by = []
        for f in linked_by:
            try:
                self.linked_by.append(proj.linkFieldToTable(f, 'variant')[-1].field)
            except Exception as e:
                raise RuntimeError('Failed to locate linked_by field {}: {}'.format(f, e))
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
                # raw values used for description only
                self.raw_refGenomes = eval(rec[1])
                # refGenome that fits this project
                self.refGenomes = eval(rec[1])
                if linked_fields is not None:
                    for key in self.refGenomes.keys():
                        self.refGenomes[key] = linked_fields
                for key in self.refGenomes.keys():
                    # no reference genome is needed
                    if key == '*':
                        self.build = self.refGenomes[key]
                    elif proj.build is None:
                        env.logger.warning('Project does not have a primary build. Using {} from the annotation database'.format(key))
                        proj.setRefGenome(key)
                        self.build = self.refGenomes[key]
                        break
                    elif key == proj.build:
                        self.build = self.refGenomes[key]
                    elif key == proj.alt_build:
                        self.alt_build = self.refGenomes[key]
        #
        if self.anno_type == 'variant' and ((self.build is not None and
            len(self.build) != 4) or (self.alt_build is not None and
            len(self.alt_build) != 4)):
            raise ValueError('There should be four linking fields for variant annotation databases.')
        if self.anno_type == 'position' and ((self.build is not None and
            len(self.build) != 2) or (self.alt_build is not None and
            len(self.alt_build) != 2)):
            raise ValueError('There should be two linking fields for positional annotation databases.')
        if self.anno_type == 'range' and ((self.build is not None and
            len(self.build) != 3) or (self.alt_build is not None and
            len(self.alt_build) != 3)):
            raise ValueError('There should be three linking fields for range-based annotation databases.')
        if self.description == '':
            env.logger.warning('No description for annotation database {}'
                .format(annoDB))
        if self.build is None and self.alt_build is None:
            raise ValueError('No reference genome information for annotation database {}'
                .format(annoDB))
        if self.anno_type == 'field' and len(self.linked_by) != len(self.build):
            raise RuntimeError(
                'Please specify link fields for attributes {} using parameter --linked_by'
                .format(','.join(self.build)))
        if self.linked_by:
            s = delayedAction(env.logger.info, 'Indexing linked field {}'.format(', '.join(self.linked_by)))
            self.indexLinkedField(proj, linked_by)
            del s
        if self.anno_type == 'range':
            if self.build is not None:
                self.db.binningRanges(proj.build, self.build, self.name)
            elif self.alt_build is not None:
                self.db.binningRanges(proj.alt_build, self.alt_build, self.name)

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
                    env.logger.debug(e)
            else:
                # from an annotation database
                try:
                    # FIXME: this syntax is only valid for sqlite3, we will need to fix it for the mysql engine
                    cur.execute('CREATE INDEX IF NOT EXISTS {0}.{0}_{1} ON {0} ({1} ASC);'.format(table, field))
                except Exception as e:
                    env.logger.debug(e)

    def checkLinkedFields(self, proj):
        # 
        # check if all fields are correctly linked. If many of the values are not linked
        # it might be linked through wrong linked_by fields.
        #
        try:
            proj.db.attach(os.path.join(self.dir, self.filename))
        except:
            # the database might already been attached (repeated use ..)
            pass
        cur = proj.db.cursor()
        # how many values are there in the annotation database?
        cur.execute('SELECT DISTINCT {0} FROM {1}.{1}'.format(', '.join(self.build), self.name))
        val_annoDB = set(cur.fetchall())
        # how many values are there in the project?
        cur.execute('SELECT DISTINCT {0} FROM {1}'.format(', '.join(self.linked_by), 
            ', '.join(set(['{}'.format(x.rsplit('.', 1)[0] if '.' in x else 'variant') for x in self.linked_by if '.' in x]))))
        val_proj = set(cur.fetchall())
        #
        val_common = val_annoDB & val_proj 
        env.logger.info('{} out of {} {} are annotated through annotation database {}'
            .format(len(val_common), len(val_proj), ', '.join(self.linked_by), self.name))
        #
        # if not all values are used
        if len(val_common) < len(val_annoDB):
            env.logger.warning(
                '{} out of {} values in annotation database {} are not linked to the project.'
                .format(len(val_annoDB) - len(val_common), len(val_annoDB), self.name))
            val_unused = list(val_annoDB - val_common)[:100]
            val_unused.sort()
            env.logger.debug('The {} unlinked values are: {}'.format('first 100' if len(val_unused) == 100 else len(val_unused),
                ', '.join([':'.join([str(y) for y in x]) for x in val_unused])))
        

    def describe(self, verbose=False):
        '''Describe this annotation database'''
        print('Annotation database {} {}'.format(self.name, '(version {})'
            .format(self.version) if self.version else ''))
        if self.description is not None:
            print('\n'.join(textwrap.wrap(
                '{:<23} {}'.format('Description:', self.description),
                subsequent_indent=' '*2, width=78)))
        print('{:<23} {}'.format('Database type:', self.anno_type))
        # get linking fields
        if verbose:
            # number of records
            cur = self.db.cursor()
            cur.execute('SELECT value FROM {}_info WHERE name="num_records";'.format(self.name))
            num_records = int(cur.fetchone()[0])
            print('{:<23} {:,}'.format('Number of records:', num_records))
            fields = self.refGenomes.values()[0]
            
            # Get number of unique keys
            cur.execute('SELECT value FROM {}_info WHERE name="distinct_keys";'.format(self.name))
            count = int(cur.fetchone()[0])
            
            #
            if self.anno_type == 'variant':
                print('{:<23} {:,}'.format('Distinct variants:', count))
            elif self.anno_type == 'position':
                print('{:<23} {:,}'.format('Distinct positions:', count))
            elif self.anno_type == 'range':
                print('{:<23} {:,}'.format('Distinct ranges:', count))
            elif self.anno_type == 'field':
                print('{:<23} {:,}'.format('Distinct entries:', count))
        #
        for key in self.raw_refGenomes:
            print('{:<23} {}'.format(
                '{} {}:'.format('Reference genome', key), 
                ', '.join(self.raw_refGenomes[key])))
        for field in self.fields:
            if not verbose:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(field.name, field.comment),
                    subsequent_indent=' '*24, width=78)))
            else:
                print('\nField:                  {}'.format(field.name))
                numeric = False
                if 'chromosome' in field.type.lower():
                    print('Type:                   chromosome')
                elif 'position' in field.type.lower():
                    print('Type:                   integer')
                    numeric = True
                elif 'int' in field.type.lower():
                    print('Type:                   integer')
                    numeric = True
                elif 'float' in field.type.lower():
                    print('Type:                   float')
                    numeric = True
                else:
                    print('Type:                   string')
                if field.comment:
                    print('\n'.join(textwrap.wrap(
                        '{:<23} {}'.format('Comment:', field.comment),
                        width=78, subsequent_indent=' '*24)))
                cur.execute('SELECT missing_entries FROM {0}_field WHERE name="{1}";'
                    .format(self.name, field.name))
                missing = cur.fetchone()[0]
                #
                print('Missing entries:        {:,} {}'
                    .format(missing, '({:.1f}% of {:,} records)'
                        .format(100. * missing/num_records, num_records) if missing else ''))
                if missing == num_records:
                    continue
                cur.execute('SELECT distinct_entries {2} FROM {0}_field WHERE name={1};'.format(
                    self.name, self.db.PH, ', min_value, max_value' if numeric else ''), (field.name,))
                res = cur.fetchone()
                print('Unique Entries:         {:,}'.format(res[0]))
                if numeric:
                    print('Range:                  {} - {}'.format(res[1], res[2]))


class fileFMT:
    def __init__(self, name, fmt_args=[]):
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
        self.preprocessor = None
        self.header = None
        # for export only
        self.export_by_fields = ''
        self.order_by_fields = ''
        self.additional_exports = None
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
            url = 'format/{}.fmt'.format(name)
            try:
                fmt = downloadFile(url, quiet=True)
            except Exception as e:
                raise ValueError('Failed to download format specification file {}.fmt: {}'.format(name, e))
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
                            raise ValueError('Incorrect key {} in section {}. '
                                'Only field, adj, and comment are allowed.'.format(item, section))
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
                        if item not in ['index', 'type', 'adj', 'fmt', 'comment'] + defaults.keys():
                            raise ValueError('Incorrect key {} in section {}. Only index, type, adj, fmt, and comment are allowed.'.format(item, section))
                    fields.append(
                        Field(name=section,
                            index=parser.get(section, 'index', vars=defaults),
                            type=parser.get(section, 'type', vars=defaults),
                            adj=parser.get(section, 'adj', vars=defaults) if 'adj' in items else None,
                            fmt=parser.get(section, 'fmt', vars=defaults) if 'fmt' in items else None,
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
            elif item[0] == 'preprocessor':
                self.preprocessor = item[1]
            elif item[0] == 'merge_by':
                self.merge_by_cols = [x-1 for x in eval(item[1])]
            elif item[0] == 'export_by':
                self.export_by_fields = item[1]
            elif item[0] == 'additional_exports':
                # additional files to export from the variant/sample tables
                self.additional_exports = item[1]
            elif item[0] == 'sort_output_by':
                self.order_by_fields = item[1]
            elif item[0] == 'header':
                if item[1] in ('none', 'None'):
                    self.header = None
                else:
                    try:
                        self.header = int(item[1])
                    except:
                        # in this case header is a pattern
                        self.header = re.compile(item[1])
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
                env.logger.debug('Undefined field {} in format {}.'.format(self.fields[i], filename))
                self.fields[i] = Field(name=self.fields[i], index=None, adj=None, fmt=None, type=None, comment='')
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
        if self.description is not None:
            print('\n'.join(textwrap.wrap(self.description, width=78)))
        #
        if self.preprocessor is not None:
            print('Preprocessor: {}'.format(self.preprocessor))
        #
        print('\nColumns:')
        if self.columns:
            for col in self.columns:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(col.index, col.comment),
                    subsequent_indent=' '*24, width=78)))
            if self.formatter:
                print('Formatters are provided for fields: {}'.format(', '.join(self.formatter.keys())))
        else:
            print('  None defined, cannot export to this format')
        #
        if self.input_type == 'variant':
            print('\n{0}:'.format(self.input_type))
        for fld in self.fields[self.ranges[0]:self.ranges[1]]:
            print('\n'.join(textwrap.wrap(
                '  {:<21} {}'.format(fld.name, fld.comment),
                subsequent_indent=' '*24, width=78)))
        if self.ranges[1] != self.ranges[2]:
            print('\nVariant info:')
            for fld in self.fields[self.ranges[1]:self.ranges[2]]:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(fld.name, fld.comment),
                    subsequent_indent=' '*24, width=78)))
        if self.ranges[2] != self.ranges[3]:
            print('\nGenotype:')
            for fld in self.fields[self.ranges[2]:self.ranges[3]]:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(fld.name, fld.comment),
                    subsequent_indent=' '*24, width=78)))
        if self.ranges[3] != self.ranges[4]:
            print('\nGenotype info:')
            for fld in self.fields[self.ranges[3]:self.ranges[4]]:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(fld.name, fld.comment),
                    subsequent_indent=' '*24, width=78)))
        if self.other_fields:
            print('\nOther fields (usable through parameters):')
            for fld in self.other_fields:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {}'.format(fld.name, fld.comment),
                    subsequent_indent=' '*24, width=78)))
        if self.parameters:
            print('\nFormat parameters:')
            for item in self.parameters:
                print('\n'.join(textwrap.wrap(
                    '  {:<21} {} (default: {})'.format(item[0], item[2], item[1]),
                    subsequent_indent=' '*24, width=78)))
        else:
            print('\nNo configurable parameter is defined for this format.\n')



class PipelineDescription:
    def __init__(self, name, extra_args=[]):
        '''Pipeline configuration file'''
        self.description = None
        self.pipeline_vars = {}
        self.pipelines = {}
        self.pipeline_descriptions = {}
        #
        if name.endswith('.pipeline'):
            if os.path.isfile(name):
                self.name = os.path.split(name)[-1].rsplit('.', 1)[0]
                args = self.parseArgs(name, extra_args)
                self.parsePipeline(name, defaults=args) 
            else:
                raise ValueError('Pipeline description file not found: {}'
                    .format(name))
        elif os.path.isfile(name + '.pipeline'):
            self.name = os.path.split(name)[-1]
            args = self.parseArgs(name + '.pipeline', extra_args)
            self.parsePipeline(name + '.pipeline', defaults=args) 
        elif name.endswith('.pipeline') and os.path.isfile(name):
            self.name = os.path.split(name)[-1][:-9]
            args = self.parseArgs(name, extra_args)
            self.parsePipeline(name, defaults=args) 
        else:
            url = 'pipeline/{}.pipeline'.format(name)
            try:
                pipeline = downloadFile(url, quiet=True)
            except Exception as e:
                raise ValueError('Failed to download pipeline specification '
                    'file {}.pipeline: {}'.format(name, e))
            self.name = name
            args = self.parseArgs(pipeline, extra_args)
            self.parsePipeline(pipeline, defaults=args)

    def parseArgs(self, filename, fmt_args):
        fmt_parser = ConfigParser.SafeConfigParser()
        fmt_parser.read(filename)
        parameters = fmt_parser.items('DEFAULT')
        parser = argparse.ArgumentParser(prog='vtools CMD --pipeline {}'
            .format(os.path.split(filename)[-1]),
            description='Parameters to override parameters of existing steps.')
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

    def parsePipeline(self, filename, defaults):
        parser = ConfigParser.SafeConfigParser()
        # this allows python3 to read .pipeline file with non-ascii characters,
        # but there is no simple way to make it python2 compatible.
        #with open(filename, 'r', encoding='UTF-8') as inputfile:
        #    parser.readfp(inputfile)
        parser.read(filename)
        # sections?
        sections = parser.sections()
        if 'pipeline description' not in sections:
            raise ValueError("Missing section 'pipeline description' in "
                "configuration file {}".format(filename))
        #
        for section in sections:
            if section.lower() == 'pipeline description':
                for item in parser.items(section, vars=defaults):
                    if item[0] == 'description':
                        self.description = item[1]
                    elif item[0].endswith('_description'):
                        self.pipeline_descriptions[item[0].rsplit('_', 1)[0]] = item[1]
                    elif item[0] in defaults or item[0].endswith('_comment'):
                        pass
                    else:
                        self.pipeline_vars[item[0]] = item[1]
            else:
                try:
                    pname, pidx = section.rsplit('_', 1)
                except Exception as e:
                    raise ValueError('Invalid section name {} in pipeline description file {}'
                        .format(section, filename))
                if not pidx.isdigit():
                    raise ValueError('Index of a pipeline step should be an integer: {} provided'
                        .format(pidx))
                if pname not in self.pipelines:
                    self.pipelines[pname] = []
                try:
                    items = [x[0] for x in parser.items(section, raw=True)]
                    if 'action' not in items:
                        raise ValueError('Missing item "action" in section {}.'.format(section))
                    for item in items:
                        if item.endswith('_comment'):
                            continue
                        if item not in ['input', 'input_emitter', 'action', 'comment'] + defaults.keys():
                            raise ValueError('Incorrect key {} in section {}. '
                                'Only input, input_emitter, action, and comment are allowed.'
                                .format(item, section))
                    command = PipelineCommand(index=pidx,
                            input=parser.get(section, 'input', vars=defaults) if 'input' in items else '',
                            input_emitter=parser.get(section, 'input_emitter', vars=defaults) if 'input_emitter' in items else '',
                            action=parser.get(section, 'action', vars=defaults) if 'action' in items else '',
                            comment=parser.get(section, 'comment', raw=True) if 'comment' in items else '')
                    self.pipelines[pname].append(command)
                except Exception as e:
                    raise ValueError('Invalid section {}: {}'.format(section, e))
        # sort steps
        for pname in self.pipelines:
            self.pipelines[pname].sort(key=lambda x: int(x[0].rsplit('_')[-1]))
        # 
        # validate
        for pname in self.pipeline_descriptions:
            if pname not in [x.lower() for x in  self.pipelines.keys()]:
                raise ValueError('Invalid item {0}_description because pipeline '
                    '{0} is not defined in this file (available: {1}).'
                    .format(pname, ', '.join(self.pipelines.keys())))
        for pname, pipeline in self.pipelines.items():
            if pname.lower() not in self.pipeline_descriptions:
                self.pipeline_descriptions[pname.lower()] = ''
            for idx, cmd in enumerate(pipeline):
                if cmd is None:
                    raise ValueError('Invalid pipeline {}. Step {} is left unspecified.'
                        .format(pname, idx+1))
                if not cmd.action:
                    raise ValueError('Missing or empty action for step {} of pipeline {}'
                        .format(pname, idx + 1))
     
    def describe(self):
        if self.description is not None:
            print('\n'.join(textwrap.wrap(self.description, width=78)))
        #
        text = 'Available pipelines: {}'.format(', '.join(sorted(self.pipelines.keys())))
        print('\n' + '\n'.join(textwrap.wrap(text, width=78, subsequent_indent=' '*8)))
        for pname, pipeline in sorted(self.pipelines.items()):
            print('\n' + '\n'.join(textwrap.wrap('Pipeline "{}":  {}'
                .format(pname, self.pipeline_descriptions[pname.lower()]),
                width=78)))
            for idx, step in enumerate(pipeline):
                text = '{:<22}'.format('  {}_{}:'.format(pname, step.index)) + step.comment
                print('\n'.join(textwrap.wrap(text, width=78, subsequent_indent=' '*22)))
        #
        if self.parameters:
            print('\nPipeline parameters:')
            for item in self.parameters:
                #
                text = '  ' + item[0] + \
                    (' '*(22-len(item[0])-2) if len(item[0])<20 else ' ') + \
                    (item[2] + ' ' if item[2] else '') + \
                    ('(default: {})'.format(item[1]) if item[1] else '')
                print('\n'.join(textwrap.wrap(text, subsequent_indent=' '*22)))
        else:
            print('\nNo configurable parameter is defined for this format.\n')


class AnnoDBWriter:
    '''
    A class to initiate and insert annotation database
    '''
    def __init__(self, name, fields, anno_type, description, version, build, use_existing_db=False,
        overwrite_existing_fields=False):
        # remove extension .db
        self.name = name[:-3] if name.lower().endswith('.db') else name
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
        env.logger.debug('Creating {}_field table'.format(self.name))
        self.createFieldsTable()
        #
        for field in self.fields:
            cur.execute('INSERT INTO {0}_field (name, field, type, comment) VALUES ({1},{1},{1},{1});'.format(self.name, self.db.PH),
                (field.name, field.index, field.type, field.comment))
        self.db.commit()
        #
        # creating the info table
        env.logger.debug('Creating {}_info table'.format(self.name))
        query = 'INSERT INTO {0}_info VALUES ({1},{1});'.format(self.name, self.db.PH)
        self.createInfoTable()
        cur.execute(query, ('name', self.name))
        cur.execute(query, ('anno_type', self.anno_type))
        cur.execute(query, ('description', self.description))
        cur.execute(query, ('version', self.version))
        cur.execute(query, ('build', str(self.build)))
        self.db.commit()
        env.logger.debug('Creating table {}'.format(self.name))
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
        cur.execute('SELECT name, field, "", "", type, comment from {}_field;'.format(self.name))
        # cur_fields is made a class member to make others know what are available
        self.cur_fields = []
        for rec in cur:
            self.cur_fields.append(Field(*rec))
        # add new fields
        s = delayedAction(env.logger.info, 'Adding fields to existing result database')
        for field in self.fields:
            # name already exist
            if field.name in [x.name for x in self.cur_fields]:
                cf = [x for x in self.cur_fields if x.name == field.name][0] 
                if field.type != cf.type:
                    raise ValueError('Type mismatch for new field {}: existing {}, new {}'.format(field.name, cf.type, field.type))
                if overwrite_existing_fields:
                    env.logger.warning('Results in field {} will be overwritten.'.format(field.name))
                #else:
                #    if field.name not in [x for tmp in self.build.values() for x in tmp]: 
                #        raise ValueError('Cannot modify database {} because field {} already exists. '
                #            'Please use test option --name to add a new suffix to this field, '
                #            'write the results to a different database (option --to_db), or use '
                #            'option --update to force updating the existing fields.'.format(self.name, field.name))
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
        env.logger.debug('Creating annotation table {} using query\n{}'.format(self.name, query))
        cur = self.db.cursor()
        try:
            cur.execute(query)
        except Exception as e:
            raise ValueError('Failed to create table: {}'.format(e))

    def finalize(self):
        '''Create index and get statistics of the database'''
        cur = self.db.cursor()
        # creating indexes
        s = delayedAction(env.logger.info, 'Creating indexes (this can take quite a while)')
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
        # binning ranges
        if self.anno_type == 'range':
            for build, keys in self.build.items():
                self.db.binningRanges(build, keys, self.name)
        s = delayedAction(env.logger.info, 'Analyzing and tuning database ...')
        # This is only useful for sqlite
        self.db.analyze()
        # calculating database statistics
        cur.execute('SELECT COUNT(*) FROM (SELECT DISTINCT {} FROM {});'.format(', '.join(self.build.values()[0]), self.name))
        count = cur.fetchone()[0]
        cur.execute('INSERT INTO {0}_info VALUES ({1}, {1});'.format(self.name, self.db.PH), ('distinct_keys', str(count)))
        cur.execute('INSERT INTO {0}_info VALUES ({1}, {1});'.format(self.name, self.db.PH), ('num_records', self.db.numOfRows(self.name)))
        del s
        for field in self.fields:
            s = delayedAction(env.logger.info, 'Calculating column statistics for field {}'.format(field.name))
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
def unlock_proj(*args):
    env.unlock_all()
    env.logger.error('Killed by signal')
    sys.exit(1)

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
    # the following make Project a singleton class
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            # *args, **kwargs are not passed to avoid
            # DeprecationWarning: object.__new__() takes no parameters
            # cls._instance = super(Singleton, cls).__new__(cls, *args, **kwargs) 
            cls._instance = super(Project, cls).__new__(cls) #, *args, **kwargs)
        return cls._instance

    def __init__(self, name=None, build=None, new=False, verbosity=None, 
        verify=True, **kwargs):
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
            elif not name.replace('_', '').isalnum():
                raise ValueError('A project name can only contain alpha-numeric characters and underscores.')
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
        self.revision = VTOOLS_REVISION
        #
        # create a temporary directory
        self.db = DatabaseEngine()
        self.db.connect(self.proj_file)
        env.cache_dir = self.loadProperty('__option_cache_dir', None)
        #
        env.treat_missing_as_wildtype = self.loadProperty('__option_treat_missing_as_wildtype', None)
        env.association_timeout = self.loadProperty('__option_association_timeout', None)
        env.logfile_verbosity = self.loadProperty('__option_logfile_verbosity', None)
        #env.check_update = self.loadProperty('__option_check_update', True)
        env.associate_num_of_readers = self.loadProperty('__option_associate_num_of_readers', None)
        if verbosity is None and not new:
            # try to get saved verbosity level
            verbosity = self.loadProperty('__option_verbosity', None)
        # set global verbosity level and temporary directory
        env.verbosity = verbosity
        # env.verbosity will affect the creation of logger
        env.logger = '{}{}.log'.format('>' if new else '>>', self.name)
        env.logger.debug('')
        env.logger.debug(env.command_line)
        # if option temp_dir is set, the path will be used
        # if not, None will be passed, and a temporary directory will be used.
        try:
            env.temp_dir = self.loadProperty('__option_temp_dir', None)
        except Exception as e:
            env.logger.warning('Failed to use temporary directory as specified in runtime option temp_dir: {}'.format(e))
            env.temp_dir = None
        env.logger.debug('Using temporary directory {}'.format(env.temp_dir))
        #
        if new: 
            self.create(build=build, **kwargs)
            self.checkUpdate()
        else:
            self.open(verify)

    def create(self, build, **kwargs):
        '''Create a new project'''
        # open the project file
        env.logger.info(VTOOLS_COPYRIGHT)
        env.logger.info(VTOOLS_CITATION)
        env.logger.info(VTOOLS_CONTACT)
        env.logger.info('Creating a new project {}'.format(self.name))
        self.db = DatabaseEngine(engine='sqlite3', batch=kwargs.get('batch', 10000))
        self.db.connect(self.proj_file)
        #
        engine = kwargs.get('engine', 'sqlite3')
        if engine == 'mysql':
            # project information table
            self.createProjectTable()
            self.saveProperty('version', self.version)
            self.saveProperty('revision', self.revision)
            self.saveProperty('engine', engine)
            self.saveProperty('host', kwargs.get('host'))
            self.saveProperty('user', kwargs.get('user'))
            # FIXME: I am saving passwd as clear text here...
            self.saveProperty('passwd', kwargs.get('passwd'))
            self.saveProperty('batch', kwargs.get('batch', 10000))
            self.saveProperty('__option_verbosity', env.verbosity)
            # turn to the real online engine
            env.logger.debug('Connecting to mysql server')
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
        env.logger.debug('Creating core tables')
        self.createProjectTable()
        self.saveProperty('engine', engine)
        self.saveProperty('version', self.version)
        self.saveProperty('revision', self.revision)
        self.saveProperty('batch', kwargs.get('batch', 10000))
        self.saveProperty('__option_verbosity', env.verbosity)
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
        env.logger.debug('Opening project {}'.format(self.proj_file))
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
            # pragma, set to None if the key does not exist. In this case
            # the system default will be used.
            env.sqlite_pragma = self.loadProperty('__option_sqlite_pragma', None)
            # pass things like batch ... and re-connect
            # env['sqlite_pragma'] will be used 
            self.db = DatabaseEngine(engine='sqlite3', batch=self.batch)
            self.db.connect(self.proj_file)
        # loading other options if they have been set
        env.import_num_of_readers = self.loadProperty('__option_import_num_of_readers', None)
        env.local_resource = self.loadProperty('__option_local_resource', None)
        #
        # existing project
        cur = self.db.cursor()
        self.version = self.loadProperty('version', '1.0')
        self.revision = self.loadProperty('revision', '0')
        self.build = self.loadProperty('build')
        self.alt_build = self.loadProperty('alt_build')
        self.annoDB = []
        for db in eval(self.loadProperty('annoDB', '[]').replace('${local_resource}', env._local_resource)):
            try:
                # remove path, remove version string, and suffix
                db_name = os.path.split(db)[-1].split('-')[0]
                if db_name.endswith('.DB'):
                    db_name = db_name[:-3]
                linked_by = eval(self.loadProperty('{}_linked_by'.format(db_name), default='[]'))
                anno_type = self.loadProperty('{}_anno_type'.format(db_name), default='None')
                linked_fields = eval(self.loadProperty('{}_linked_fields'.format(db_name), default='None'))
                self.annoDB.append(AnnoDB(self, db, linked_by, anno_type, linked_fields))
                self.db.attach(db)
            except Exception as e:
                env.logger.warning('Cannot open annotation database {}: {}'.format(db, e))
        #
        # get existing meta information
        # FIXME: these are not handled correctly for now
        self.variant_meta = self.db.getHeaders('variant_meta')
        self.sample_meta = self.db.getHeaders('sample_meta')
        #
        if self.revision != VTOOLS_REVISION:
            proj_rev = int(self.revision)
            vtools_rev = int(VTOOLS_REVISION)
            if proj_rev > vtools_rev:
                env.logger.warning('Opening a project that is created by an '
                    'newer version of vtools ({}, rev {}) is dangerous.'
                    .format(self.version, self.revision))
            elif proj_rev < vtools_rev:
                # upgrade project
                try:
                    self.upgrade(proj_rev)
                except Exception as e:
                    env.logger.warning('Skip upgrading project: {}'.format(e))
        if verify:
            try:
                self.checkIntegrity()
            except Exception as e:
                env.logger.warning('Skip checking integrity of project: {}'.format(e))
        # 
        try:
            self.analyze()
        except Exception as e:
            env.logger.warning('Skip analyzing project: {}'.format(e))
        

    def checkIntegrity(self):
        '''Check if the project is ok...(and try to fix it if possible)'''
        for table in ['project', 'filename', 'sample', 'variant']:
            if not self.db.hasTable(table):
                raise RuntimeError('Corrupted project: missing table {}'.format(table))
        #
        headers = self.db.getHeaders('variant')
        if self.alt_build is not None:
            if not ('alt_bin' in headers and 'alt_chr' in headers and 'alt_pos' in headers):
                env.logger.warning('Disable alternative reference genome because of missing column {}.'.format(
                    ', '.join([x for x in ('alt_bin', 'alt_chr', 'alt_pos') if x not in headers])))
                self.alt_build = None
                self.saveProperty('alt_build', None)
        #
        # missing index on master variant table, this will happen after all data is imported
        if not self.db.hasIndex('variant_index'):
            self.createIndexOnMasterVariantTable()
            if not self.db.hasIndex('variant_index'):
                raise RuntimeError('Corrupted project: failed to create index on master variant table.')

    def checkUpdate(self):
        res = ResourceManager()
        try:
            res.getRemoteManifest()
            # check small files (.ann, .format and .pipelines) and
            # update at most 5 of them them silently
            changed = res.checkUpdate(5)
            if len(changed) == 1:
                env.logger.warning('Resource file {} has been updated. Please '
                    'update it using command "vtools admin --update_resource '
                    'existing".'.format(changed[0]))
            elif len(changed) > 1:
                env.logger.warning('Resource files {} have been updated. Please '
                    'update them using command "vtools admin --update_resource '
                    'existing".'.format(', '.join(changed)))
        except Exception as e:
            # if the machine is not connected to the internet,
            # do not get any update
            env.logger.debug('Failed to check update: {}'.format(e))
            return
        #
        # check current version of variant tools.
        try:
            (version_file, header) = urllib.urlretrieve('http://vtools.houstonbioinformatics.org/CURRENT_VERSION.txt')
            with open(version_file, 'r') as version:
                current_version = version.readline().strip()
            if [int(x) for x in re.sub('\D', ' ', current_version).split()] > \
                [int(x) for x in re.sub('\D', ' ', VTOOLS_VERSION).split()]:
                env.logger.warning('A new version of variant tools ({}) is available.'
                    .format(current_version))
        except Exception as e:
            env.logger.debug('Failed to check latest version: {}'.format(e))
        finally:
            # remove manifest_file
            urllib.urlcleanup()


    def analyze(self, force=False):
        '''Automatically analyze project to make sure queries are executed optimally.
        '''
        if self.db.engine != 'sqlite3':
            return
        cur = self.db.cursor()
        tables = self.db.tables()
        cur = self.db.cursor()
        s = delayedAction(env.logger.info, 'Analyzing project')
        for tbl in tables:
            # do not analyze these small tables
            if tbl in ['sample', 'filename', 'project']:
                continue
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
        
    def upgrade(self, proj_revision):
        for ver, rev, proc in project_format_history:
            # for example, if the project version if 1.0.6
            # it will call the upgrade version for 1.0.7 and higher
            if proj_revision < rev:
                proc(self)
        # mark the version of the project
        self.saveProperty('version', VTOOLS_VERSION)
        self.saveProperty('revision', VTOOLS_REVISION)

    def useAnnoDB(self, db):
        '''Add annotation database to current project.'''
        # DBs in different paths but with the same name are considered to be the same.
        env.logger.info('Using annotation DB {} in project {}.'.format(db.name, self.name))
        env.logger.info(db.description)
        if db.name not in [x.name for x in self.annoDB]:
            self.annoDB.append(db)
        # if db.name is in the list, put it to the end because it is the last
        # one to be used, and might reply on some other db
        else:
            i = [x.name for x in self.annoDB].index(db.name)
            self.annoDB.append(db)
            self.annoDB.pop(i)
        self.saveProperty('annoDB', str([os.path.join(x.dir, x.filename).replace(env._local_resource, '${local_resource}') for x in self.annoDB]))
        # an annotation database might be re-used with a different linked_field
        self.saveProperty('{}_linked_by'.format(db.name), str(db.linked_by))
        self.saveProperty('{}_anno_type'.format(db.name), str(db.anno_type))
        self.saveProperty('{}_linked_fields'.format(db.name), str(db.build))
        # 
        # if a field database, connect and check 
        if db.linked_by:
            db.checkLinkedFields(self)

    def close(self):
        '''Write everything to disk...'''
        self.db.commit()
        # temporary directories are cleared each time
        try:
            # __exit__() will also call close(), I do not know if I should remove
            # it from __exit__. Anyway, the second call to close() will try to 
            # remove temp_dir again...
            if os.path.isdir(env.temp_dir):
                shutil.rmtree(env.temp_dir)
        except Exception as e:
            env.logger.warning('Failed to remove temporary directory {0}: {1}'.format(env.temp_dir, e))
        
    def loadProperty(self, key, default=None):
        '''Retrieve property from the project table'''
        cur = self.db.cursor()
        try:
            cur.execute('SELECT value FROM project WHERE name={0};'.format(self.db.PH), (key,))
            return cur.fetchone()[0]
        except Exception as e:
            #env.logger.debug(e)
            #env.logger.warning('Failed to retrieve value for project property {}'.format(key))
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
            pass
        self.db.commit()

    def removeProperty(self, key):
        cur = self.db.cursor()
        try:
            cur.execute('DELETE FROM project WHERE name={};'.format(self.db.PH), (key,))
        except Exception as e:
            pass
        self.db.commit()

    def remove(self):
        '''Remove the current project'''
        # step 1: remove database
        if self.db.engine == 'mysql':
             env.logger.info('Removing database {}'.format(self.name))
             self.db.removeDatabase(self.name)
        # step 2: remove genotype
        if self.db.hasDatabase(self.name + '_genotype'):
            env.logger.info('Removing database {}_genotype'.format(self.name))
            self.db.removeDatabase(self.name + '_genotype')
        # step 3: remove temp database
        if self.db.hasDatabase(self.name + '_temp'):
            env.logger.info('Removing database {}_temp'.format(self.name))
            self.db.removeDatabase(self.name + '_temp')
        # step 4: remove project file
        env.logger.info('Removing project file {}'.format(self.proj_file))
        os.remove(self.proj_file)
        # step 5: remove lck and log file
        if os.path.isfile(self.proj_file[:-5] + '.lck'):
            os.remove(self.proj_file[:-5] + '.log')
        env.logger.info('Removing log file {}'.format(self.proj_file[:-5] + '.log'))
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
            env.logger.error('Cannot change reference genome of an existing project.')
            sys.exit(1)
        self.build = build
        self.saveProperty('build', build)

    #
    # Functions to create core and optional tables.
    #
    def createProjectTable(self):
        env.logger.debug('Creating table project')
        cur = self.db.cursor()
        cur.execute('''\
            CREATE TABLE project (
                name VARCHAR(40) NOT NULL PRIMARY KEY,
                value VARCHAR(256) NULL
            )''')

    def createFilenameTable(self):
        env.logger.debug('Creating table filename')
        cur = self.db.cursor()
        cur.execute('''\
            CREATE TABLE filename (
                file_id INTEGER PRIMARY KEY {0},
                filename VARCHAR(256) NOT NULL
            )'''.format(self.db.AI))
        # create index
        try:
            cur.execute('''CREATE UNIQUE INDEX filename_index ON filename (filename ASC);''')
        except Exception as e:
            # the index might already exists
            return

    def createMasterVariantTable(self, fields=[]):
        '''Create a variant table with name. Fail if a table already exists.'''
        env.logger.debug('Creating table variant')
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
        self.describeTable('variant', 'Master variant table', True, False)
        self.createIndexOnMasterVariantTable()

    def createIndexOnMasterVariantTable(self, quiet=False):
        # create indexes
        #
        if not quiet:
            s = delayedAction(env.logger.info, 'Creating indexes on master variant table. This might take quite a while.')
        try:
            #
            # Index on the primary reference genome is UNIQUE when there is no alternative reference
            # genome. If there is, multiple variants from the alternative reference genome might
            # be mapped to the same coordinates on the primary reference genome. I have tried to 
            # set some of the coordinates to NULL, but the uniqueness problem becomes a problem
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
            if not self.db.hasIndex('variant_index'):
                if self.alt_build is not None:
                    self.db.execute('''CREATE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
                else:
                    self.db.execute('''CREATE UNIQUE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists, this does not really matter
            env.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
        try:
            if self.alt_build and not self.db.hasIndex('variant_alt_index'):
                #
                # Index on alternative reference genome is not unique because several variants might
                # be mapped to the same coordinates in the alternative reference genome. 
                #
                self.db.execute('''CREATE INDEX variant_alt_index ON variant (alt_bin ASC, alt_chr ASC, alt_pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists, this does not really matter
            env.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
        if not quiet:
            del s

    def dropIndexOnMasterVariantTable(self):
        # before bulk inputting data, it is recommended to drop index.
        #
        # NOTE: for mysql, it might be better to use alt index disable/rebuild
        #
        s = delayedAction(env.logger.info, 'Dropping indexes of master variant table. This might take quite a while.')
        try:
            # drop index has different syntax for mysql/sqlite3.
            if self.db.hasIndex('variant_index'):
                self.db.dropIndex('variant_index', 'variant')
        except Exception as e:
            # the index might not exist
            env.logger.debug(e)
        #
        try:
            if self.alt_build and self.db.hasIndex('variant_alt_index'):
                self.db.dropIndex('variant_alt_index', 'variant')
        except Exception as e:
            # the index might not exist
            env.logger.debug(e)
        del s

    def createVariantTable(self, table, temporary=False):
        '''Create a variant table with name. Fail if a table already exists.
        '''
        if table == 'variant':
            raise ValueError('This function cannot be used to create a master variant table')
        env.logger.debug('Creating table {}'.format(table))
        self.db.execute('''CREATE {0} TABLE {1} (
                variant_id INTEGER PRIMARY KEY
            );'''.format('TEMPORARY' if temporary else '', table))
        self.db.commit()

    def describeTable(self, table, message, save_date=False, save_cmd=False):
        '''Attach a message to a table, optional date and command
        to create the table can also be saved.'''
        self.saveProperty('__desc_of_{}'.format(table), message)
        if save_date:
            self.saveProperty('__date_of_{}'.format(table), time.strftime('%b%d', time.gmtime()))
        if save_cmd:
            self.saveProperty('__cmd_of_{}'.format(table), env.command_line)

    def descriptionOfTable(self, table):
        '''Get description of table'''
        return (self.loadProperty('__desc_of_{}'.format(table), ''),
            self.loadProperty('__date_of_{}'.format(table), ''),
            self.loadProperty('__cmd_of_{}'.format(table), '') )
            
    def createSampleTableIfNeeded(self, fields=[], table='sample'):
        if self.db.hasTable(table):
            return
        env.logger.debug('Creating table {}'.format(table))
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
        '''Create a table ``genotype_??`` to store genotype data'''
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
            #env.logger.debug(e)
            #pass

    def createVariantMetaTableIfNeeded(self):
        if self.variant_meta is None:
            env.logger.warning('Tried to create a variant meta table without valid meta fields')
            return
        if self.db.hasTable('variant_meta'):
            return
        env.logger.debug('Creating table variant_meta')
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
            env.logger.warning('Tried to create a sample meta table without valid meta fields')
            return
        if self.db.hasTable('sample_meta'):
            return
        env.logger.debug('Creating table sample_meta')
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
            raise ValueError('{} is not found or is not a variant table.'
                .format(decodeTableName(table)))
        env.logger.info('Removing table {}'.format(decodeTableName(table)))
        self.db.removeTable(table)

    def selectSampleByPhenotype(self, cond):
        '''Select samples by conditions such as "aff=1", return IDs as a sorted list'''
        cur = self.db.cursor()
        try:
            cur.execute('SELECT sample_id FROM sample LEFT OUTER JOIN filename ON '
                'sample.file_id = filename.file_id {} ORDER BY sample.sample_name;'
                .format(' WHERE {}'.format(cond) if cond.strip() else ''))
            IDs = [x[0] for x in cur.fetchall()]
            # return a tuple to avoid future change of order
            return tuple(IDs)
        except Exception as e:
            raise ValueError('Failed to retrieve samples by condition "{}": {}'
                .format(cond, e))

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
                env.logger.info('Removing {1} and all its samples ({0})'.format('{} samples'.format(len(samples[f])) if len(samples[f]) > 1 else 'sample {}'.format(samples[f][0]), f)) 
                cur.execute('DELETE FROM filename WHERE filename = {}'.format(self.db.PH), (f,))
            else:
                env.logger.info('Removing {} imported from {}'.format('{} samples'.format(len(samples[f])) if len(samples[f]) > 1 else 'sample {}'.format(samples[f][0]), f)) 
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
        # first remove from the master variant table
        cur.execute('DELETE FROM variant WHERE variant_id IN (SELECT variant_id FROM {});'.format(table))
        if not cur.rowcount:
            env.logger.warning('Table {} is empty. No variant has been removed from the project.'.format(table))
            return
        else:
            env.logger.info('{} variants are removed from the master variant table.'.format(cur.rowcount))
        # remove from other variant tables
        for t in self.getVariantTables():
            if t.lower() in ['variant', table.lower()]:
                continue
            cur.execute('DELETE FROM {} WHERE variant_id IN (SELECT variant_id FROM {});'.format(t, table))
            env.logger.info('{} variants are removed from table {}'.format(cur.rowcount, t))
        # get sample_ids
        cur.execute('SELECT sample_id, sample_name FROM sample;')
        samples = cur.fetchall()
        for ID, name in samples:
            if not self.db.hasIndex('{0}_genotype.genotype_{1}_index'.format(self.name, ID)):
                cur.execute('CREATE INDEX {0}_genotype.genotype_{1}_index ON genotype_{1} (variant_id ASC)'.format(self.name, ID))
            cur.execute('DELETE FROM {}_genotype.genotype_{} WHERE variant_id IN (SELECT variant_id FROM {});'\
                .format(self.name, ID, table))
            env.logger.info('{} genotypes are removed from sample {}'.format(cur.rowcount, name))
        # remove the table itself
        env.logger.info('Removing table {} itself'.format(table))
        self.db.removeTable(table)

    def removeGenotypes(self, cond):
        '''Remove genotype according to certain conditions'''
        # get sample_ids
        cur = self.db.cursor()
        cur.execute('SELECT sample_id, sample_name FROM sample;')
        samples = cur.fetchall()
        env.logger.info('Removing genotypes from {} samples using criteria "{}"'.format(len(samples), cond))
        for ID, name in samples:
            cur.execute('DELETE FROM {}_genotype.genotype_{} WHERE {};'\
                .format(self.name, ID, cond))
            env.logger.info('{} genotypes are removed from sample {}'.format(cur.rowcount, name))

    def renameSamples(self, cond, name1, name2=None):
        '''If name2 is none, rename selected samples to specified name1. 
        Otherwise, replace the first occurance of name1 to name2'''
        cur = self.db.cursor()
        try:
            where_clause = ' WHERE {}'.format(cond) if cond.strip() else ''
            cur.execute('SELECT sample_name FROM sample WHERE sample_id IN '
                '(SELECT sample_id FROM sample LEFT OUTER JOIN filename ON '
                '  sample.file_id = filename.file_id {});'
                .format(where_clause))
            names = [x[0] for x in cur.fetchall()]
            if not names:
                env.logger.warning('No sample is selected using condition "{}"'
                    .format(cond))
            if name2 is None:
                # rename all names to name1
                cur.execute('UPDATE sample SET sample_name={} WHERE sample_id IN '
                    '(SELECT sample_id FROM sample LEFT OUTER JOIN filename ON '
                    'sample.file_id = filename.file_id {});'
                    .format(self.db.PH, where_clause), (name1, ))
                env.logger.info('{} samples with names {} are renamed to {}'
                    .format(cur.rowcount, ', '.join(sorted(set(names))), name1))
            else:
                for oldname in sorted(set(names)):
                    newname = oldname.replace(name1, name2, 1)
                    if cond.strip():
                        where_clause = ' WHERE {} AND sample.sample_name={}'.format(cond, self.db.PH)
                    else:
                        where_clause = ' WHERE sample.sample_name={}'.format(self.db.PH)
                    if newname == oldname:
                        continue
                    # rename all names to name1
                    cur.execute('UPDATE sample SET sample_name={} WHERE sample_id IN '
                        '(SELECT sample_id FROM sample LEFT OUTER JOIN filename ON '
                        'sample.file_id = filename.file_id {});'
                        .format(self.db.PH, where_clause),
                        (newname, oldname))
                    env.logger.info('Rename {} sample{} with name {} to {}'
                        .format(cur.rowcount, 's' if cur.rowcount > 1 else '',
                            oldname, newname))
        except Exception as e:
            raise ValueError('Failed to retrieve samples by condition "{}": {}'
                .format(cond, e))

    def mergeSamples(self):
        '''Merge samples with the same name to the same samples'''
        cur = self.db.cursor()
        query = 'SELECT sample_name, sample_id FROM sample ORDER BY sample_name'
        env.logger.debug('Executing {}'.format(query))
        cur.execute(query)
        # a map of sample_name to multiple sample ids
        samples = {}
        for rec in cur:
            if rec[0] in samples:
                samples[rec[0]].append(rec[1])
            else:
                samples[rec[0]] = [rec[1]]
        # who is going to be merged?
        merged = {}
        for name, ids in samples.iteritems():
            if len(ids) > 1:
                merged[name] = sorted(ids)
        if len(merged) == 0:
            env.logger.info('No sample is merged because all samples have unique names')
            return
        count = sum([len(x) for x in merged.values()])
        env.logger.info('{} samples that share identical names will be merged to {} samples'
            .format(count, len(merged)))
        # if all samples will be involved in merging, using
        # a separate database to speed up merging.
        inPlace=(count != sum([len(x) for x in samples.values()]))
        #
        self.db.attach(self.name + '_genotype')
        if not inPlace:
            if os.path.isfile('{}/{}_genotype_merged.DB'.format(env.cache_dir, self.name)):
                os.remove('{}/{}_genotype_merged.DB'.format(env.cache_dir, self.name))
            self.db.attach('{}/{}_genotype_merged.DB'.format(env.cache_dir, self.name),
                '{}_genotype_merged'.format(self.name))
        # merge genotypes
        prog = ProgressBar('Merging samples', count)
        copied = 0
        for name, ids in merged.iteritems():
            #
            # ranges of variant_id is used to check overlap of genotypes. If the ranges
            # overlap, a more thorough check will be used.
            if inPlace:
                new_table = '{}_genotype._tmp_{}'.format(self.name, ids[0])
            else:
                new_table = '{}_genotype_merged.genotype_{}'.format(self.name, ids[0])
            #
            try:
                # remove existing temporary table if exists
                self.db.removeTable(new_table)
            except:
                pass
            #
            # get schema
            new_fields = []
            # get a list of fields for old ids
            old_fields = {}
            cur.execute('SELECT name, sql FROM {}_genotype.sqlite_master '
                    'WHERE type="table" AND name in ({})'
                    .format(self.name, 
                    ', '.join(["'genotype_{}'".format(x) for x in ids])))
            for name, schema in cur:
                try:
                    fields = [x.strip() for x in schema.split(',')]
                    fields[0] = fields[0].split('(')[1].strip()
                    fields[-1] = fields[-1].rsplit(')', 1)[0].strip()
                    old_fields[int(name.rsplit('_', 1)[-1])] = ', '.join([fld.split(None, 1)[0] for fld in fields])
                    # if two fields have the same name and different types, they will be merged silently.
                    new_fields.extend([fld for fld in fields if fld.split(None, 1)[0] not in [x.split(None, 1)[0] for x in new_fields]])
                except Exception as e:
                    raise RuntimeError('Corrupted genotype database: Failed to '
                        'get structure of genotype table {}'.format(id))
            #
            query = 'CREATE TABLE {} ({})'.format(new_table, ', '.join(new_fields))
            env.logger.debug('Executing {}'.format(query))
            try:
                cur.execute(query)
            except Exception as e:
                raise RuntimeError('Failed to create new genotype table: {}'.format(e))
            # copying genotypes to the temp table
            ranges = []
            for id in ids:
                query = ('SELECT MIN(variant_id), MAX(variant_id) FROM {0}_genotype.genotype_{1}'
                    .format(self.name, id))
                try:
                    cur.execute(query)
                except Exception as e:
                    raise RuntimeError('Failed to get ID range of table {}: {}'
                        .format(id, e))
                ranges.append(cur.fetchone())
                try:
                    cur.execute('INSERT INTO {1} ({3}) SELECT * FROM {0}_genotype.genotype_{2};'.format(
                        self.name, new_table, id, old_fields[id]))
                except Exception as e:
                    try:
                        # remove existing temporary table if exists
                        self.db.removeTable(new_table)
                    except:
                        pass
                    prog.done()
                    raise RuntimeError('Failed to merge genotype tables: {}'.format(e))
                copied += 1
                prog.update(copied)
            #
            # verify table if the ranges overlap
            check_overlap = False
            for i in range(len(ranges) - 1):
                for j in range(i + 1, len(ranges)):
                    # skip when either range (a,b) or (c,d) are empty
                    # i.e., no variants for a sample from a particular file
                    if not any(ranges[i]) or not any(ranges[j]):
                        continue
                    # range overlap (a,b) with (c,d) <===> a <= d and b >= c 
                    if ranges[i][0] <= ranges[j][1] and ranges[i][1] >= ranges[j][0]:
                        check_overlap = True
                        break
                if check_overlap:
                    break
            #
            if check_overlap:
                query = ('SELECT COUNT(*), COUNT(DISTINCT variant_id) FROM {};'
                    .format(new_table))
                try:
                    cur.execute(query)
                except Exception as e:
                    raise RuntimeError('Failed to check overlap of table {}: {}'
                        .format(new_table, e))
                counts = cur.fetchone()
                if counts[0] != counts[1]:
                    try:
                        # remove existing temporary table if exists
                        self.db.removeTable(new_table)
                    except:
                        pass
                    prog.done()
                    raise ValueError('Failed to merge samples with name {} because '
                        'there are {} genotypes for {} unique variants.'
                        .format(name, counts[0], counts[1]))
        #
        prog.done()
        # the above steps are slow but will not affect project (process can be terminated)
        # the following will be fast
        for name, ids in merged.iteritems():
            # change filenames
            filenames = []
            for id in ids:
                query = ('SELECT filename FROM sample JOIN filename ON '
                    'sample.file_id = filename.file_id WHERE sample_id = {}'
                    .format(self.db.PH))
                try:
                    cur.execute(query, (id, ))
                except Exception as e:
                    raise RuntimeError('Failed to get filename for sample {}: {}'.
                        format(id, e))
                rec = cur.fetchone()
                filenames.extend(rec[0].split(','))
            #
            filenames = ','.join(sorted(list(set(filenames))))
            try:
                cur.execute('INSERT INTO filename (filename) VALUES ({0});'.format(
                    self.db.PH), (filenames,))
                file_id = cur.lastrowid
            except:
                # existing file id??
                cur.execute('SELECT file_id FROM filename WHERE filename = {};'.format(
                    self.db.PH), (filenames,))
                file_id = cur.fetchone()[0]
            #
            # step 4: if things are doing all right, remove/update existing tables
            for idx, id in enumerate(ids):
                if idx > 0:
                    cur.execute('DELETE FROM sample WHERE sample_id = {}'.format(self.db.PH), (id,))
                else:
                    cur.execute('UPDATE sample SET file_id={0} WHERE sample_id = {0}'.format(self.db.PH),
                        (file_id, ids[0]))
            # step 5: prepare tables to be removed
            if inPlace:
                for idx, id in enumerate(ids):
                    self.db.renameTable('{}_genotype.genotype_{}'.format(self.name, id),
                        '__genotype_{}'.format(id))
                self.db.renameTable('{}_genotype._tmp_{}'.format(self.name, ids[0]),
                    'genotype_{}'.format(ids[0]))
        # finally, remove filenames that are associated with no sample. The individual files can then
        # be imported again, which I am not sure is good or bad
        self.db.execute('DELETE FROM filename WHERE filename.file_id NOT IN (SELECT file_id FROM sample)')
        self.db.commit()
        if not inPlace:
            self.db.detach('{}_genotype_merged'.format(self.name))
            self.db.detach('{}_genotype'.format(self.name))
            os.remove('{}_genotype.DB'.format(self.name))
            os.rename('{}/{}_genotype_merged.DB'.format(env.cache_dir, self.name),
                '{}_genotype.DB'.format(self.name))
        # 
        # actually remove obsolete tables, this can be slow but interruption of
        # command will not break the database
        if inPlace:
            prog = ProgressBar('Removing obsolete tables', count)
            for idx, (name, ids) in enumerate(merged.iteritems()):
                for id in ids:
                    self.db.removeTable('{}_genotype.__genotype_{}'.format(self.name, id))
                    prog.update(idx + 1)
            self.db.commit()
            prog.done()


    def createVariantMap(self, table='variant', alt_build=False):
        '''Create a map of all variants from specified table.
        The dictionary looks like dict[(chr, ref, alt)][pos] = (id, 0)
        to avoid repeating ref, alt all the time.
        '''
        variantIndex = {}
        cur = self.db.cursor()
        numVariants = self.db.numOfRows(table)
        if numVariants == 0:
            return variantIndex
        env.logger.debug('Creating local indexes for {:,} variants'.format(numVariants));
        where_clause = 'WHERE variant_id IN (SELECT variant_id FROM {})'.format(table) if table != 'variant' else ''
        if alt_build:
            cur.execute('SELECT variant_id, alt_chr, alt_pos, ref, alt FROM variant {};'.format(where_clause))
        else:
            cur.execute('SELECT variant_id, chr, pos, ref, alt FROM variant {};'.format(where_clause))
        prog = ProgressBar('Getting existing variants', numVariants)
        for count, rec in enumerate(cur):
            # zero for existing loci
            key = (rec[1], rec[3], rec[4])
            if key in variantIndex:
                variantIndex[key][rec[2]] = (rec[0], 0)
            else:
                variantIndex[key] = {rec[2]: (rec[0], 0)}
            #variantIndex[(rec[1], rec[3], rec[4])][rec[2]] = (rec[0], 0)
            if count % self.db.batch == 0:
                prog.update(count)
        prog.done()
        return variantIndex

    def summarize(self):
        '''Summarize key features of the project
        '''
        # FIXME: more summary
        info =  'Project name:                {}\n'.format(self.name)
        info += 'Primary reference genome:    {}\n'.format(self.build)
        info += 'Secondary reference genome:  {}\n'.format(self.alt_build)
        #info += 'Database engine:             {}\n'.format(self.db.engine)
        #
        # list all runtime options as (name, val) pairs
        opts = [(x, self.loadProperty('__option_{}'.format(x), None)) for x in env.persistent_options]
        info += 'Runtime options:             {}\n'.format(
            ', '.join(['{}={}'.format(name, val) for name,val in opts if val is not None]))
        tables = [decodeTableName(x) for x in self.getVariantTables()]
        info += 'Variant tables:              {}\n'.format(
            '\n'.join(sorted([' '*29 + x for x in tables])).lstrip())
        info += 'Annotation databases:        {}\n'.format(
            '\n'.join([' '*29 + os.path.join(x.dir, x.name)
            + (' ({})'.format(x.version) if x.version else '') for x in self.annoDB]).lstrip())
        return info

    def saveSnapshot(self, name, message, files):
        '''Save snapshot'''
        if name.endswith('.tar') or name.endswith('.tar.gz') or name.endswith('.tgz'):
            filename = name
            mode = 'w' if name.endswith('.tar') else 'w:gz'
        elif name.replace('_', '').isalnum():
            filename = os.path.join(env.cache_dir, 'snapshot_{}.tar'.format(name))
            mode = 'w'
        else:
            raise ValueError('Snapshot name should be a filename with extension .tar, .tgz, or .tar.gz, or a name without any special character.')
        #
        s = delayedAction(env.logger.info, 'Creating snapshot')
        with tarfile.open(filename, mode) as snapshot:
            readme_file = os.path.join(env.cache_dir, '.snapshot.info')
            with open(readme_file, 'w') as readme:
                readme.write('Snapshot of variant tools project {}.\n'.format(self.name))
                readme.write('Name: {}\n'.format(name))
                readme.write('Date: {}\n'.format(time.strftime('%b%d %H:%M:%S', time.gmtime())))
                readme.write('Info: {}\n'.format(message))
            # add .snapshot.info file
            snapshot.add(readme_file, '.snapshot.info')
            env.logger.debug('Adding .snapshot.info'.format(self.name))
            s = delayedAction(env.logger.info, 'Copying project')
            snapshot.add('{}.proj'.format(self.name), arcname='snapshot.proj')
            env.logger.debug('Adding {}.proj as snapshot.proj'.format(self.name))
            del s
            if os.path.isfile('{}_genotype.DB'.format(self.name)):
                s = delayedAction(env.logger.info, 'Copying genotypes')
                snapshot.add('{}_genotype.DB'.format(self.name), arcname='snapshot_genotype.DB')
                env.logger.debug('Adding {}_genotype.DBj as snapshot_genotype.DB'.format(self.name))
                del s
            os.remove(readme_file)
            if files is not None:
                for f in files:
                    env.logger.debug('Adding {}'.format(f))
                    s = delayedAction(env.logger.info, 'Adding {}'.format(f))
                    snapshot.add(f)
                    del s
        # add a warning message if the snapshot starts with 'vt_'
        if name.startswith('vt_'):
            env.logger.warning('Snapshots with name starting with vt_ is reserved for public snapshots for documentation and training purposes.') 

    def loadSnapshot(self, name):
        '''Load snapshot'''
        #
        if name.endswith('.tar') or name.endswith('.tar.gz') or name.endswith('.tgz'):
            snapshot_file = name
            mode = 'r' if name.endswith('.tar') else 'r:gz'
        elif name.replace('_', '').isalnum():
            snapshot_file = os.path.join(env.cache_dir, 'snapshot_{}.tar'.format(name))
            mode = 'r'
        else:
            raise ValueError('Snapshot name should be a filename with extension .tar, .tgz, or .tar.gz, or a name without any special character.')
        #
        if not os.path.isfile(snapshot_file):
            # donload it from online?
            if name.startswith('vt_'):
                # only snapshots with name starting with vt_
                try:
                    print('Downloading snapshot {}.tar.gz from online'.format(name))
                    snapshot_file = downloadFile('snapshot/' + name + '.tar.gz', quiet=False)
                except:
                    raise ValueError('Failed to download snapshot {}.'.format(name))
            else:
                raise ValueError('Snapshot {} does not exist locally.'.format(name))
        #
        # close project
        self.db.close()
        try:
            with tarfile.open(snapshot_file, mode) as snapshot:
                all_files = snapshot.getnames()
                # old snapshot uses file README. The new format has .snapshot.info and will
                # treat README as user-provided file.
                info_file = '.snapshot.info' if '.snapshot.info' in all_files else 'README'
                all_files.remove(info_file)
                # project
                s = delayedAction(env.logger.info, 'Load project')
                if 'snapshot.proj' in all_files:
                    env.logger.debug('Extracting snapshot.proj as {}.proj'.format(self.name))
                    snapshot.extract('snapshot.proj')
                    os.rename('snapshot.proj', '{}.proj'.format(self.name))
                    all_files.remove('snapshot.proj')
                elif '{}.proj'.format(self.name) in all_files:
                    env.logger.debug('Extracting {}.proj'.format(self.name))
                    # an old version of snapshot saves $name.proj
                    snapshot.extract('{}.proj'.format(self.name))
                    all_files.remove('{}.proj'.format(self.name))
                else:
                    raise ValueError('Invalid snapshot. Missing project database')
                del s
                # genotype
                s = delayedAction(env.logger.info, 'Load genotypes')
                if 'snapshot_genotype.DB' in all_files:
                    env.logger.debug('Extracting snapshot_genotype.DB as {}_genotype.DB'.format(self.name))
                    snapshot.extract('snapshot_genotype.DB'.format(self.name))
                    os.rename('snapshot_genotype.DB', '{}_genotype.DB'.format(self.name))
                    all_files.remove('snapshot_genotype.DB')
                elif '{}_genotype.DB'.format(self.name) in all_files:
                    env.logger.debug('Extracting {}_genotype.DB'.format(self.name))
                    # an old version of snapshot saves $name.proj
                    snapshot.extract('{}_genotype.DB'.format(self.name))
                    all_files.remove('{}_genotype.DB'.format(self.name))
                else:
                    # this is ok because a project might not have any genotype
                    pass
                del s
                # other files
                for f in all_files:
                    s = delayedAction(env.logger.info, 'Extracting {}'.format(f))
                    env.logger.debug('Extracting {}'.format(f))
                    if os.path.isfile(f):
                        env.logger.warning('Ignore existing file {}.'.format(f))
                        continue
                    snapshot.extract(f)
                    del s
        except Exception as e:
            raise ValueError('Failed to load snapshot: {}'.format(e))
        finally:
            # re-connect the main database for proper cleanup
            self.db = DatabaseEngine()
            self.db.connect(self.proj_file)
        
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
                        binningTable = '__rng_' + table + '_' + encodeTableName('_'.join([self.build] + db.build))
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.bin = {0}.bin AND variant.chr = {0}.chr '
                                'AND variant.pos >= {0}.start AND variant.pos <= {0}.end '
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= '{0}.rowid = {1}.range_id'.format(table, binningTable))
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
                        binningTable = '__rng_' + table + '_' + encodeTableName('_'.join([self.alt_build] + db.alt_build))
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.alt_bin = {0}.bin AND variant.alt_chr = {0}.chr '
                                'AND variant.alt_pos >= {0}.start AND variant.alt_pos <= {0}.end '
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            # FIXME: how to use bin here?
                            link= '{0}.rowid = {1}.range_id'.format(table, binningTable))]
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
                        binningTable = '__rng_' + table + '_' +  encodeTableName('_'.join([self.build] + db.build))
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.bin = {0}.bin AND variant.chr = {0}.chr '
                                'AND variant.pos >= {0}.start AND variant.pos <= {0}.end '
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= '{0}.rowid = {1}.range_id'.format(table, binningTable))]
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
                        binningTable = '__rng_' + table + '_' + encodeTableName('_'.join([self.alt_build] + db.alt_build))
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, binningTable),
                            link= 'variant.alt_bin = {0}.bin AND variant.alt_chr = {0}.chr '
                                'AND variant.alt_pos >= {0}.start AND variant.alt_pos <= {0}.end '
                                    .format(binningTable, table)),
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= '{0}.rowid = {1}.range_id'.format(table, binningTable))]
                    else:
                        raise ValueError('Unsupported annotation type {}'.format(db.anno_type))
                else:
                    raise ValueError('Database does not define any reference genome.')
        raise ValueError('Failed to locate field {}'.format(field))
        

class MaintenanceProcess(Process):
    '''This class starts a separate process to tune the database, e.g.
    create indexes for genotypes. When active_flag is false, it will
    exit itself.'''
    def __init__(self, proj, jobs, active_flag):
        Process.__init__(self)
        self.name = proj.name
        if not jobs:  # if no job is specified
            self.jobs = {'genotype_index': []}
        else:
            self.jobs = jobs
        self.active_flag = active_flag

    def createIndexes(self, sample_IDs):
        # specified IDs will have high priority to be handled
        try:
            db = DatabaseEngine()
            db.connect('{}_genotype.DB'.format(self.name))
            cur = db.cursor()
            # get all tables, this is sqlite specific, and genotype table specific
            cur.execute('SELECT name FROM sqlite_master WHERE type="table" AND name LIKE "genotype_%"')
            all_indexes = set([x[0] + '_index' for x in cur.fetchall()])
            cur.execute('SELECT name FROM sqlite_master WHERE type="index" AND name LIKE "genotype_%"')
            cur_indexes = set([x[0] for x in cur.fetchall()])
            missing_indexes = all_indexes - cur_indexes
            if len(missing_indexes) == 0:
                return
        except KeyboardInterrupt as e:
            # interrupted just return, nothing harmful is done.
            db.close()
            return
        #
        env.logger.debug('Creating indexes for {} genotype tables'.format(len(missing_indexes)))
        try:
            # we process IDs in sample_IDs first ...
            for idx in [x for x in missing_indexes if int(x[9:-6]) in sample_IDs] + \
                [x for x in missing_indexes if int(x[9:-6]) not in sample_IDs]:
                if not self.active_flag.value:
                    break
                cur.execute('CREATE INDEX {0} ON {1} (variant_id)'.format(idx, idx[:-6]))
                db.commit()
        except KeyboardInterrupt as e:
            # if keyboard interrupt, stop, but not immediately
            self.active_flag.value = 0
        finally:
            # make sure the database is properly closed...
            db.commit()
            db.close()

    def run(self):
        if 'genotype_index' in self.jobs:
            self.createIndexes(set(self.jobs['genotype_index']))

class ProjCopier:
    def __init__(self, proj, dir, vtable, samples, genotypes):
        self.proj = proj
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
                #env.logger.debug(sql)
                cur.execute(sql)
            except Exception as e:
                env.logger.debug(e)
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
            env.logger.debug('Removing {} unselected samples'.format(len(removed)))
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
                # env.logger.debug(sql)
                cur.execute(sql)
            except Exception as e:
                env.logger.debug(e)
            # copying data over
            if self.vtable == 'variant':
                if self.genotypes:
                    # copy selected genotypes
                    try:
                        cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0} WHERE {1};'.format(table, self.genotypes))
                    except:
                        raise RuntimeError('Failed to copy {} with condition {}'.format(table, self.genotypes))
                else:
                    # copy all variants and genotypes
                    cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0};'.format(table))
            else:
                if self.genotypes:
                    try:
                        cur.execute('INSERT INTO {0} SELECT * FROM __fromGeno.{0} WHERE ({2}) AND (__fromGeno.{0}.variant_id IN (SELECT variant_id FROM __fromDB.{1}));'.format(table, self.vtable, self.genotypes))
                    except:
                        raise RuntimeError('Failed to copy {} with condition {}. The field might not exist.'.format(table, self.genotypes))
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
        try:
            if self.geno_file is not None:
                copied_samples = self.copySamples()
            else:
                copied_samples = 0
            env.logger.info('{} variants and {} samples are copied'.format(copied_variants, copied_samples))
            # wait for the thread to close
            s = delayedAction(env.logger.info, 'Create indexes')
            if thread.is_alive():
                thread.join()
            del s
        finally:
            # re-connect the main database for proer cleanup
            self.proj.db = DatabaseEngine()
            self.proj.db.connect(self.proj.proj_file)


class VariantMapper(threading.Thread):
    '''The worker thread read variants from all projects and create 
    id maps along the way. This mapper does not sort variants, uses 
    more RAM, but it is much faster because it can run in paralle..'''
    def __init__(self, projects, alt_build, status):
        self.projects = projects
        self.alt_build = alt_build
        self.status = status
        threading.Thread.__init__(self, name='Read and map variants')
        # set it to daemon so that it will stop after the master thread is killed
        self.daemon = True

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
                            # there exists a duplicate
                            idMaps[id] = (vv[pos], 1)
                            keep_all = False
                        else:
                            # a new variant is found
                            vv[pos] = new_id
                            idMaps[id] = (new_id, 0)
                            new_id += 1
                    else:
                        # a new variant is found 
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
            identical_ids = True
            for _old_id, (_new_id, _is_duplicate) in idMaps.iteritems():
                if _old_id != _new_id:
                    identical_ids = False
                    break
            # identical_ids means new_id and old_id are the same, no need to map
            self.status.set(proj, 'identical_ids', identical_ids)
            # keep_all means keeping all variants because there is no overlap (e.g. merge SNPs and INDELs)
            # it is set to False when there is at least one variant that overlap with a previous
            # project
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
    def __init__(self, queue, status):
        self.queue = queue
        self.status = status
        threading.Thread.__init__(self, name='cache variants')
        # set it to daemon so that it will stop after the master thread is killed
        self.daemon = True

    def run(self):
        while True:
            item = self.queue.get()
            if item is None:
                self.queue.task_done()
                break
            # get parameters
            self.src_proj = item
            self.identical_ids = self.status.get(self.src_proj, 'identical_ids')
            self.keep_all = self.status.get(self.src_proj, 'keep_all')
            if not (self.identical_ids and self.keep_all):
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
        if self.identical_ids and self.keep_all:
            return
        db = DatabaseEngine()
        db.connect(self.cache_proj)
        db.attach(self.src_proj, '__fromDB')
        # create index on __fromDB
        cur = db.cursor()
        cur.execute('CREATE INDEX __fromDB.__id_map_idx ON __id_map (old_id ASC);')
        # create tables
        cur.execute('''SELECT sql FROM __fromDB.sqlite_master WHERE type="table" 
            AND (NOT name LIKE "sqlite_%") AND name != "__id_map";''')
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
            # if identical_ids, no need to map ids,
            if self.identical_ids:
                # table variant:
                #     if keep_all: do not copy anything ...
                #     otherwise:   copy some of them (is_dup == 0)
                # 
                # Other variant tables:
                #     we should always copy all variants
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
                if table == 'variant' and not self.keep_all:
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
                env.logger.debug('Caching table {} of project {} ({})'.format(table, self.src_proj,
                    query))
                cur.execute(query)
        db.detach('__fromDB')
        db.close()

class SampleProcessor(threading.Thread):
    def __init__(self, queue, status):
        self.queue = queue
        self.status = status
        threading.Thread.__init__(self, name='cache samples')
        # set it to daemon so that it will stop after the master thread is killed
        self.daemon = True
    
    def run(self):
        while True:
            item = self.queue.get()
            if item is None:
                self.queue.task_done()
                break
            # get parameters
            self.src_proj = item
            self.identical_ids = self.status.get(self.src_proj, 'identical_ids')
            self.samples = self.status.get(self.src_proj, 'old_ids')
            #
            if self.identical_ids:
                self.status.set(self.src_proj, 'completed', 2 + len(self.samples))
            else:
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
                raise ValueError('Cannot recreate genotype table {} from {}, '
                    'please check the integrity of the database.'
                    .format(_old_id, self.src_geno))
            sql = sql[0]
            try:
                cur.execute(sql)
            except Exception as e:
                env.logger.debug(e)
            #
            # copy data over
            headers = db.getHeaders('genotype_{}'.format(_old_id))
            query = '''INSERT INTO genotype_{0} SELECT new_id {1} FROM __geno.genotype_{2} 
                LEFT JOIN __proj.__id_map ON __id_map.old_id = __geno.genotype_{2}.variant_id;'''\
                .format(_old_id, ''.join([', {}'.format(x) for x in headers[1:]]), _old_id)
            env.logger.debug('Caching sample {} of project {}'.format(_old_id, self.src_proj))
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
        #
        self.projects = projects
        self.status = status
        threading.Thread.__init__(self, name='copy variants')
        # set it to daemon so that it will stop after the master thread is killed
        self.daemon = True
        
    def run(self):
        db = DatabaseEngine()
        db.connect(self.proj.proj_file)
        for idx, proj in enumerate(self.projects):
            identical_ids = self.status.get(proj, 'identical_ids')
            keep_all = self.status.get(proj, 'keep_all')
            # get cache_proj
            db.attach(proj, '__fromDB')
            if not (identical_ids and keep_all):
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
                if not db.hasTable('__fromDB.{}'.format(table)):
                    continue
                # get fields of table because source and destination might have different fields
                fields = [x[0] for x in db.fieldsOfTable('__fromDB.{}'.format(table))]
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
                if identical_ids and (table != 'variant' or keep_all):
                    query = '''INSERT OR IGNORE INTO {0} ({1}) SELECT * FROM __fromDB.{0};'''.format(table, ', '.join(fields))
                else:
                    query = '''INSERT OR IGNORE INTO {0} ({1}) SELECT * FROM __cacheDB.{0};'''.format(table, ', '.join(fields))
                env.logger.debug('Copying table {} from project {} ({}, {})'.format(table, proj,
                    identical_ids, keep_all))
                cur.execute(query)
                db.commit()
            #
            db.detach('__fromDB')
            if not (identical_ids and keep_all):
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
        #
        self.projects = projects
        self.status = status
        #
        threading.Thread.__init__(self, name='copy samples')
        # set it to daemon so that it will stop after the master thread is killed
        self.daemon = True

    def run(self):
        # connect to cache geno
        db = DatabaseEngine()
        db.connect(self.proj.proj_file.replace('.proj', '_genotype.DB'))
        count = 0
        for proj in self.projects:
            identical_ids = self.status.get(proj, 'identical_ids')
            if identical_ids:
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
                env.logger.debug('Running {}'.format(sql))
                cur.execute(sql)
                query = 'INSERT INTO genotype_{} SELECT * FROM __geno.genotype_{};'\
                    .format(_new_id, _old_id)
                env.logger.debug('Copying sample {} from project {}: {}'.format(_old_id, proj, query))
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
            print (key, self.tasks[key]['completed'])
        #print self.tasks

class ProjectsMerger:
    def __init__(self, proj, dirs, jobs):
        '''Check if merge is possible, set primary and reference genome
            and set self.projects and self.structure '''
        self.proj = proj
        self.db = proj.db
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
                env.logger.warning('Remove duplicate merge {}'.format(proj_file))
            #
            self.db.attach(files[0], '__fromDB')
            # 
            # get primary and alternative reference genome
            cur = self.db.cursor()
            # primary reference genome
            cur.execute('SELECT value FROM __fromDB.project WHERE name={}'
                .format(self.db.PH), ('build',))
            build = cur.fetchone()
            # alternative reference genome
            cur.execute('SELECT value FROM __fromDB.project WHERE name={}'  
                .format(self.db.PH), ('alt_build',))
            alt_build = cur.fetchone()
            #
            if build is None or alt_build is None:
                env.logger.warning('Ignoring invalid or empty project.'.format(proj_file))
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
                    s = delayedAction(env.logger.info, 
                        'Adding alternative reference genome {} to the project.'
                        .format(self.proj.alt_build))
                    cur = self.db.cursor()
                    headers = self.db.getHeaders('variant')
                    for fldName, fldType in [('alt_bin', 'INT'), ('alt_chr', 'VARCHAR(20)'), ('alt_pos', 'INT')]:
                        if fldName in headers:
                            continue
                        self.db.execute('ALTER TABLE variant ADD {} {} NULL;'
                            .format(fldName, fldType))
                    del s
            elif build[0] != self.proj.build:
                raise ValueError('Primary reference genome of project ({} of {}) '
                    'does not match that of the current project ({}).'
                    .format(build[0], proj_file, self.proj.build))
            elif alt_build[0] != self.proj.alt_build:
                raise ValueError('Alternative reference genome of project '
                    '({} of {}) does not match that of the current project ({})'
                    .format(alt_build[0], proj_file, self.proj.alt_build))
            #
            # copy table message and runtime options
            cur.execute('SELECT name, value FROM __fromDB.project;')
            for name, value in cur:
                if name.startswith('__option_') or name.startswith('__desc_of_') or \
                    name.startswith('__date_of_') or name.startswith('__cmd_of_'):
                    self.proj.saveProperty(name, value)
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
            else:
                tables = [x for x in self.db.tables('__fromDB') if not x.startswith('__')]
                tables.sort()
                for table in tables:
                    # new table?
                    if table not in structure:
                        tbl = table.split('.')[-1]
                        if tbl.startswith('__'):
                            continue
                        structure[tbl] = self.db.fieldsOfTable('__fromDB.{}'.format(tbl))
                    else:
                        flds = self.db.fieldsOfTable('__fromDB.{}'.format(table))
                        for fld in flds:
                            if fld[0].lower() not in [x[0].lower() for x in structure[table]]:
                                env.logger.warning('Field {} in table {} does not exist in all projects.'
                                    .format(fld[0], table))
                                structure[table].append(fld)
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
        # create all variant tables
        for table in structure:
            if table in ['project', 'filename']:
                continue
            env.logger.debug('Creating table {} with columns {}'
                .format(table, ', '.join([x[0] for x in structure[table]])))
            cur.execute('CREATE TABLE {} ({});'.format(table, ', '.join([' '.join(x) for x in structure[table]])))
       
    def mapSamples(self, status):
        '''Population filename and sample table, return id maps
        '''
        cur = self.db.cursor()
        filenames = []
        #
        duplicated_samples = 0
        for proj in self.projects:
            self.db.attach(proj, '__proj')
            self.db.attach(proj.replace('.proj', '_genotype.DB'), '__geno')
            #
            # handling filename
            cur.execute('SELECT file_id, filename FROM __proj.filename;')
            filename_records = cur.fetchall()
            old_sample_id = []
            new_sample_id = []
            for old_file_id, filename in filename_records:
                if filename in filenames:
                    cur.execute('SELECT count(sample_id) FROM __proj.sample WHERE file_id={};'
                        .format(self.db.PH), (old_file_id, ))
                    cnt = cur.fetchone()[0]
                    duplicated_samples += int(cnt)
                    cur.execute('SELECT file_id FROM filename WHERE filename = {0};'
                        .format(self.db.PH), (filename,))
                    new_file_id = cur.fetchone()[0]
                    duplicate = True
                else:
                    filenames.append(filename)
                    #
                    cur.execute('INSERT INTO filename (filename) VALUES ({0});'
                        .format(self.db.PH), (filename, ))
                    new_file_id = cur.lastrowid
                    duplicate = False
                # get samples
                cur.execute('SELECT sample_id FROM __proj.sample WHERE file_id={};'
                    .format(self.db.PH), (old_file_id, ))
                old_sid = [x[0] for x in cur.fetchall()]
                new_sid = []
                # use __proj.sample because the project sample might have more
                # fields (phenotypes)
                headers = self.db.getHeaders('__proj.sample')
                for sid in old_sid:
                    cur.execute('''INSERT INTO sample ('file_id', {0}) SELECT {1}, {0} 
                        FROM __proj.sample WHERE __proj.sample.sample_id = {1};'''\
                        .format(', '.join(headers[2:]), self.db.PH), (new_file_id, sid))
                    new_sid.append(cur.lastrowid)
                #
                old_sample_id.extend(old_sid)
                new_sample_id.extend(new_sid)
            env.logger.debug('Mapping sample_ids of project {} from {} to {}'
                .format(proj, old_sample_id, new_sample_id))
            status.set(proj, 'old_ids', old_sample_id)
            status.set(proj, 'new_ids', new_sample_id)
            self.db.detach('__proj')
            self.db.detach('__geno')
        return duplicated_samples
  
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
        # 3. sample processed (need to be updated with samples to be cached)
        for proj in self.projects:
            status.add(proj, {'completed': 0, 'stage': 0, 'scheduled': False, 'total_count': 3})
        # this will set for each project
        #  old_ids: sample id of the original projects
        #  new_ids: sample id in the new project
        duplicated_samples = self.mapSamples(status)
        if duplicated_samples > 0:
            env.logger.warning('{} samples from the same source files have been '
                'copied, leading to potentially duplicated samples.'
                .format(duplicated_samples))
        #
        # update project status for total number of jobs
        for proj in self.projects:
            status.set(proj, 'total_count', 2 + len(status.get(proj, 'old_ids')))
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
            VariantProcessor(self.vcQueue, status).start()
        #
        # start all sample cachers
        self.scQueue = Queue.Queue()
        for j in range(nJobs):
            SampleProcessor(self.scQueue, status).start()
        #
        # create a thread to read variants
        VariantMapper(self.projects, self.proj.alt_build, status).start()
        #
        prog = None
        count = 0
        while True:
            for idx, proj in enumerate(self.projects):
                if status.canProcessVariant(proj) and self.vcQueue.qsize() < nJobs:
                    env.logger.debug('Mapping variants in {}'.format(proj))
                    status.set(proj, 'scheduled', True)
                    self.vcQueue.put(proj)
                if status.canProcessSample(proj) and self.vcQueue.qsize() + self.scQueue.qsize() < nJobs:
                    env.logger.debug('Mapping sample variants in {}'.format(proj))
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
# PROJECT UPGRADE TREE
#
def remove_duplicate_genotype(proj):
    #
    # remove all uplicate genotypes from sample tables
    # version 1.0.7 handles this automatically but projects created
    # by older version of variant tools might still have such genotypes
    db = DatabaseEngine()
    db.connect('{}_genotype'.format(proj.name))
    tables = [x for x in db.tables() if x.startswith('genotype_')]
    prog = ProgressBar('Upgrading project to 1.0.7', len(tables))
    duplicated_genotype = 0
    cur = db.cursor()
    for idx, table in enumerate(tables):
        cur.execute('SELECT COUNT(*), COUNT(DISTINCT variant_id) FROM {}'.format(table))
        nRec, nVar = cur.fetchone()
        if nRec != nVar:
            cur.execute('DELETE from {0} WHERE rowid NOT IN '
                '(SELECT MAX(rowid) FROM {0} GROUP BY variant_id)'
                .format(table))
            if cur.rowcount != nRec - nVar:
                raise SystemError('Failed to remove duplicated variants from '
                    'genotype table {}'.format(table))
        duplicated_genotype += nRec - nVar
        prog.update(idx + 1)
    db.commit()
    db.close()
    prog.done()
    if duplicated_genotype:
        env.logger.warning('{} duplicated genotypes are removed from {} sample{}.'
            .format(duplicated_genotype, len(tables),
                's' if len(tables) > 1 else ''))

project_format_history = [
    # version (for documentation purpose only), revision, upgrade function
    ['1.0.7', 1915, remove_duplicate_genotype],
]

#
#
# Functions provided by this script
#
#

def initArguments(parser):
    parser.add_argument('project',
        help='''Name of a new project. This will create a new .proj file under
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
            and alternative reference genome, variant info and phenotype). The same
            variant tables from multiple samples will be merged. Samples from the
            children projects will be copied even if they were identical samples
            imported from the same source files.''')
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
                # we might not be able to open a non-exist project if it is malformed.
                # however, this will not completely remove the project if the project
                # uses a MySQL engine.
                files = glob.glob('*.proj') + glob.glob('*.lck') + \
                    glob.glob('*.proj-journal') + glob.glob('*_genotype.DB')
                if len(files) > 0:
                    for f in files:
                        try:
                            os.remove(f)
                        except:
                            # we might not be able to remove files...
                            raise OSError('Failed to remove existing project {}'.format(f))
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
                # A default value 4 is given for args.jobs because more threads usually
                # do not improve effiency
                merger = ProjectsMerger(proj, args.children, 4)
                merger.merge()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


def removeArguments(parser):
    parser.add_argument('type', choices=['project', 'tables', 'samples', 'fields',
        'geno_fields', 'annotations', 'variants', 'genotypes', 'phenotypes'],
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
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select tables to be removed')
                allTables = proj.getVariantTables()
                removed = []
                for table in args.items:
                    if '?' in table or '*' in table:
                        matched = False
                        for tbl in [decodeTableName(x) for x in allTables]:
                            if re.match(table.replace('?', '.{1}').replace('*', '.*'), tbl, re.I) and tbl not in removed:
                                try:
                                    proj.removeVariantTable(encodeTableName(tbl))
                                except Exception as e:
                                    env.logger.warning('Failed to remove table "{}": {}'
                                        .format(tbl, e))
                                removed.append(tbl)
                                matched = True
                        if not matched:
                            env.logger.warning('Name {} does not match any existing variant table.'.format(table))
                    else:
                        proj.removeVariantTable(encodeTableName(table))
                        removed.append(table)
            elif args.type == 'samples':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select samples to be removed')
                proj.db.attach(proj.name + '_genotype')
                IDs = proj.selectSampleByPhenotype(' AND '.join(['({})'.format(x) for x in args.items]))
                if len(IDs) == 0:
                    env.logger.warning('No sample is selected by condition {}'
                        .format(' AND '.join(['({})'.format(x) for x in args.items])))
                proj.removeSamples(IDs)
            elif args.type == 'fields':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select fields to be removed')
                from_table = defaultdict(list)
                for item in args.items:
                    if item.lower() in ['variant_id', 'chr', 'pos', 'alt']:
                        raise ValueError('Fields variant_id, chr, pos and alt cannot be removed')
                    all_fields = [x.lower() for x in proj.db.getHeaders('variant')]
                    if item.lower() not in all_fields:
                        raise ValueError('{} is not a valid variant info field.'.format(item))
                # remove...
                env.logger.info('Removing variant info field {}'.format(', '.join(args.items)))
                proj.db.removeFields('variant', args.items)
                if 'alt_bin' in args.items or 'alt_chr' in args.items or 'alt_pos' in args.items:
                    env.logger.info('Removing alternative reference genome '
                        'because of removal of related fields')
                    proj.alt_build = None
                    proj.saveProperty('alt_build', None)
                # it is possible that new indexes are needed
                s = delayedAction(env.logger.info, 'Rebuilding indexes')
                try:
                    proj.createIndexOnMasterVariantTable()
                except:
                    pass
                del s
            elif args.type == 'geno_fields':
                if len(args.items) == 0:
                    raise ValueError('Please specify name of genotype fields to be removed')
                if 'variant_id' in [x.lower() for x in args.items]:
                    raise ValueError('Genotypes id variant_id cannot be removed')
                if 'gt' in [x.lower() for x in args.items]:
                    raise ValueError('Genotypes field GT cannot be removed')
                proj.db.attach(proj.name + '_genotype')
                cur = proj.db.cursor()
                cur.execute('SELECT sample_id FROM sample;')
                IDs = [x[0] for x in cur.fetchall()]
                cnt = 0
                for table in ['{}_genotype.genotype_{}'.format(proj.name, id) for id in IDs]:
                    header = [x.lower() for x in proj.db.getHeaders(table)]
                    items = [x for x in args.items if x.lower() in header and x.lower not in ['variant_id', 'gt']]
                    if items:
                        cnt += 1
                        env.logger.info('Removing fields {} from genotype table {}'
                            .format(', '.join(items), table.split('_')[-1]))
                        proj.db.removeFields(table, items)
                if cnt:
                    env.logger.info('Genotype info fields {} are removed from {} samples.'.format(', '.join(items), cnt))
                else:
                    env.logger.warning('Genotype info fields {} not found in any of the samples.'.format(', '.join(items)))
            elif args.type == 'annotations':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select annotation table to be removed')
                for item in args.items:
                    removed = False
                    for i in range(len(proj.annoDB)):
                        if proj.annoDB[i].name == item:
                            env.logger.info('Removing annotation database {} from the project'.format(item))
                            proj.annoDB.pop(i)
                            removed = True
                            break
                    if not removed:
                        env.logger.warning('Cannot remove annotation database {} from the project'.format(item))
                # use the un-expanded version of _local_resource to allow continued use of '~'
                proj.saveProperty('annoDB', str([os.path.join(x.dir, x.filename).replace(env._local_resource, '${local_resource}') for x in proj.annoDB]))
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
                    raise ValueError('Please specify one or more phenotypes '
                        '(columns in the output of "vtools show samples") to be removed')
                phenos = [x.lower() for x in proj.db.getHeaders('sample')]
                toBeRemoved = []
                for item in args.items:
                    if item.lower() in ['filename', 'sample_name']:
                        env.logger.warning('Cannot remove filename or sample_name from the sample table')
                        continue
                    if item.lower() not in phenos:
                        env.logger.warning('No phenotype {} exists in the sample table.'.format(item))
                        continue
                    toBeRemoved.append(item)
                if len(toBeRemoved) == 0:
                    env.logger.warning('No valid phenotype to be removed.')
                # remove
                proj.db.removeFields('sample', toBeRemoved)
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)

def showArguments(parser):
    parser.add_argument('type', choices=['project', 'tables', 'table',
        'samples', 'phenotypes', 'genotypes', 'fields', 'annotations',
        'annotation', 'track', 'formats', 'format', 'tests', 'test', 
        'runtime_options', 'runtime_option', 'snapshot', 'snapshots', 
        'pipeline', 'pipelines'],
        nargs='?', default='project',
        help='''Type of information to display, which can be 'project' for
            summary of a project, 'tables' for all variant tables (or all
            tables if --verbosity=2), 'table TBL' for details of a specific
            table TBL, 'samples [COND]' for sample name, files from which
            samples are imported, and associated phenotypes (can be supressed
            by option --verbosity 0) of all or selected samples, 
            'phenotypes [P1 P2...]' for all or specified phenotypes of samples, 
            'fields' for fields from variant tables and all used annotation
            databases, 'annotations' for a list of all available annotation
            databases, 'annotation ANN' for details about annotation database ANN,
            'track' for information of a track file in tabixed vcf, bigWig, or
            bigBed format, 'formats' for all supported import and export formats, 
            'format FMT' for details of format FMT, 'tests' for a list of all 
            association tests, and 'test TST' for details of an association test 
            TST, 'runtime_options' for a list of runtime options and their 
            descriptions, 'runtime_option OPT' for value of specified runtime 
            option OPT, 'snapshot' for a particular snapshot by name or
            filename, 'snapshots' for a list of publicly available snapshots,
            and snapshots of the current project saved by command 
            'vtools admin --save_snapshots', 'pipeline PIPELINE'
            for details of a particular align and variant calling pipeline, and
            'pipelines' for a list of available pipelines. The default parameter
            of this command is 'project'.''')
    parser.add_argument('items', nargs='*',
        help='''Items to display, which can be, for example, name of table for
            type 'table', conditions to select samples for type 'samples', 
            a list of phenotypes for type 'phenotypes', name of an annotation
            database for type 'annotation', a pattern to selected annotation
            databases for type 'annotations', name of a format for type 'format',
            and name of an association test for type 'test'.''')
    parser.add_argument('-l', '--limit', metavar='N', type=int,
        help='''Limit output to the first N records.''')


def show(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            #
            limit_clause = ' LIMIT 0, {}'.format(args.limit) if args.limit is not None and args.limit >= 0 else ''
            omitted = '({} records omitted)'
            if args.type == 'project':
                print(proj.summarize())
            elif args.type == 'tables':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show tables"'
                        .format(', '.join(args.items)))
                all_tables = sorted(proj.getVariantTables(), key=lambda x: decodeTableName(x))
                width = max([len(x) for x in all_tables])
                if args.verbosity != '0':
                    print(('{:<' + str(width+2) + '} {:>10} {:>8} {}')
                        .format('table', '#variants', 'date', 'message'))
                for idx, table in enumerate(all_tables):
                    if args.limit is not None and idx == args.limit:
                        break
                    if args.verbosity == '0':
                        print(decodeTableName(table))
                    else:
                        desc, date, cmd = proj.descriptionOfTable(table) 
                        print(('{:<' + str(width+2) + '} {: >10,} {:>8} {}')
                            .format(decodeTableName(table), proj.db.numOfRows(table), date, 
                            '\n'.join(textwrap.wrap(desc, initial_indent='',
                            subsequent_indent=' '*(width+23)))))
                nAll = len(all_tables)
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'table':
                proj.db.attach('{}_genotype'.format(proj.name))
                if not args.items:
                    raise ValueError('Please specify a variant table to display')
                if len(args.items) > 1:
                    raise ValueError('Only a single variant table is allowed.')
                table = args.items[0]
                if proj.isVariantTable(encodeTableName(table)):
                    table = encodeTableName(table)
                else:
                    raise ValueError('{} is not a valid variant table'.format(table))
                # print description of table
                desc, date, cmd = proj.descriptionOfTable(table)
                print('{:<23} {}'.format('Name:', decodeTableName(table)))
                print('\n'.join(textwrap.wrap(
                    '{:<23} {}'.format('Description:', desc),
                    width=78, subsequent_indent=' '*24)))
                print('{:<23} {}'.format('Creation date:', date))
                print('\n'.join(textwrap.wrap(
                    '{:<23} {}'.format('Command:', cmd),
                    width=78, subsequent_indent=' '*24)))
                # 
                headers = proj.db.getHeaders(table)
                print('\n'.join(textwrap.wrap(
                    '{:<23} {}'.format('Fields:', ', '.join(headers)),
                    width=78, subsequent_indent=' '*24)))
                print('{:<23} {}'.format('Number of variants:',
                    proj.db.numOfRows(table)))
            elif args.type == 'samples':
                if not proj.db.hasTable('sample'):
                    env.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                fields = proj.db.getHeaders('sample')
                # if -v0, do not show phenotypes
                if args.verbosity == '0':
                    fields = fields[:3]
                # headers are ID, file, sample, FIELDS
                print('sample_name\tfilename{}'.format(''.join(['\t'+x for x in fields[3:]])))
                cur.execute('SELECT sample_name, filename {} FROM sample, filename '
                    'WHERE sample.file_id = filename.file_id {} ORDER BY sample_name {};'
                    .format(' '.join([','+x for x in fields[3:]]),
                    ' '.join(['AND ({})'.format(x) for x in args.items]),
                    limit_clause))
                for rec in cur:
                    if args.verbosity != '2' and len(rec[1]) > 25:
                        rec = list(rec)
                        rec[1] = rec[1][:8] + '...' + rec[1][-14:]
                    print('\t'.join(['{}'.format(x) for x in rec]))
                nAll = proj.db.numOfRows('sample')
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'phenotypes':
                if not proj.db.hasTable('sample'):
                    env.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                fields = proj.db.getHeaders('sample')[3:]
                if args.items:
                    found = []
                    unfound = []
                    for item in args.items:
                        if item.lower() in [x.lower() for x in fields]:
                            idx = [x.lower() for x in fields].index(item.lower())
                            found.append(fields[idx])
                        else:
                            unfound.append(item)
                    #
                    if unfound:
                        env.logger.warning('Phenotypes {} not found.'
                            .format(', '.join(unfound)))
                    fields = found
                # headers are ID, file, sample, FIELDS
                print('sample_name{}'.format(''.join(['\t'+x for x in fields])))
                cur.execute('SELECT sample_name {} FROM sample, filename '
                    'WHERE sample.file_id = filename.file_id ORDER BY sample_name {};'
                    .format(' '.join([','+x for x in fields]),
                    limit_clause))
                for rec in cur:
                    print('\t'.join(['{}'.format(x) for x in rec]))
                nAll = proj.db.numOfRows('sample')
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'fields':
                if len(proj.annoDB) == 0:
                    env.logger.debug('No annotation database is attached.')
                for table in proj.getVariantTables():
                    tfields = [field for field in proj.db.getHeaders(table) if field not in ('variant_id', 'bin', 'alt_bin')]
                    if tfields:
                        print('\n'.join(['{}.{}'.format(table, field) for field in tfields]))
                for db in proj.annoDB:
                    if args.verbosity == '0':
                        print('\n'.join(['{}.{}'.format(db.name, x.name) for x in db.fields]))
                    else:
                        print('\n'.join(['{}.{} {}'.format(db.name, x.name,
                            '\n'.join(textwrap.wrap(x.comment,
                                initial_indent=' '*(27-len(db.name)-len(x.name)),
                                subsequent_indent=' '*29))) for x in db.fields]))
            elif args.type == 'annotation':
                if len(args.items) == 0:
                    raise ValueError('Please specify the annotation(s) to display')
                elif len(args.items) > 1:
                    raise ValueError('Only one annotation database is allowed.')
                try:
                    dbName = args.items[0].lower()
                    if dbName.endswith('.db'):
                        dbName = dbName[:-3]
                    annoDB = [x for x in proj.annoDB if x.name.lower() == dbName][0]
                except Exception as e:
                    raise IndexError('Database {} is not currently used in the project.'
                        .format(args.items[0]))
                annoDB.describe(args.verbosity == '2')
            elif args.type == 'annotations':
                res = ResourceManager()
                res.getRemoteManifest()
                res.selectFiles(resource_type='annotation')
                nAll = 0
                displayed = 0
                for idx, (annoDB, prop) in enumerate(sorted(res.manifest.iteritems())):
                    if not annoDB.endswith('.ann'):
                        continue
                    if args.items:
                        match = False
                        for item in args.items:
                            if item.lower() in annoDB[7:-4].lower():
                                match = True
                                break
                        if not match:
                            continue
                    nAll += 1
                    if args.limit is None or args.limit > displayed:
                        if args.verbosity == '0':
                            print(annoDB[7:-4])
                        else:
                            text = '{:<23} {}'.format(annoDB[7:-4], prop[3])
                            print('\n'.join(textwrap.wrap(text, width=78,
                                    subsequent_indent=' '*24)))
                        displayed += 1
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'track':
                if not args.items:
                    raise ValueError('Please provide a track file in tabixed vcf, '
                        'bigWig, or bigBed format')
                elif len(args.items) > 1:
                    raise ValueError('Only one track file is allowed.')
                showTrack(args.items[0])
            elif args.type == 'formats':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show formats"'
                        .format(', '.join(args.items)))
                res = ResourceManager()
                res.getRemoteManifest()
                res.selectFiles(resource_type='format')
                nAll = len(res.manifest)
                for idx, (fmt, prop) in enumerate(res.manifest.iteritems()):
                    if args.limit is not None and idx == args.limit:
                        break
                    if args.verbosity == '0':
                        print(fmt[7:-4])
                    else:
                        text = '{:<23} {}'.format(fmt[7:-4], prop[3])
                        print('\n'.join(textwrap.wrap(text, width=78,
                            subsequent_indent=' '*24)))
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'format':
                if not args.items:
                    raise ValueError('Please specify a format to display')
                elif len(args.items) > 1:
                    raise ValueError('Only one file format is allowed.')
                try:
                    fmt = fileFMT(args.items[0])
                except Exception as e:
                    raise IndexError('Unrecognized input format {}: {}'
                        .format(args.items[0], e))
                fmt.describe()
            elif args.type == 'genotypes':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show genotypes"'.format(', '.join(args.items)))
                # get sample ids and attach the genotypes database
                if not proj.db.hasTable('sample'):
                    env.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                try:
                    proj.db.attach(proj.name + '_genotype') 
                except Exception as e:
                    # either the database doesn't exist or it's already been attached (should we add a method proj.db.isAttached(...)? and embedding this in the attach() method?
                    if not proj.db.hasDatabase(proj.name + '_genotype'):
                        env.logger.debug('Trying to attach a database that doesn\'t exist' + e)
                # sample headers are ID, file, sample, FIELDS
                print('sample_name\tfilename\tnum_genotypes\tsample_genotype_fields')
                cur.execute('SELECT sample.sample_id, sample_name, filename FROM sample, filename WHERE sample.file_id = filename.file_id ORDER BY sample.sample_id {};'.format(limit_clause))
                records = cur.fetchall()
                for rec in records:
                    # sample fields
                    if args.verbosity != '2' and len(rec[2]) > 25:
                        sampleFields = '{}\t{}...{}'.format(rec[1], rec[2][:8], rec[2][-14:])
                    else:
                        sampleFields = '{}\t{}'.format(rec[1], rec[2])
                    # now get sample genotype counts and sample specific fields
                    sampleId = rec[0]
                    cur.execute('SELECT count(*) FROM {}_genotype.genotype_{};'.format(proj.name, sampleId))
                    numGenotypes = cur.fetchone()[0]
                    # get fields for each genotype table
                    sampleGenotypeHeader = proj.db.getHeaders('{}_genotype.genotype_{}'.format(proj.name, sampleId))
                    sampleGenotypeFields = ','.join(['{}'.format(x) for x in sampleGenotypeHeader[1:]])  # the first field is variant id, second is GT
                    print('{}\t{}\t{}'.format(sampleFields, numGenotypes, sampleGenotypeFields))
                nAll = proj.db.numOfRows('sample')
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'tests':
                # it is very bad idea to use circular import, but I have no choice
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show tests"'.format(', '.join(args.items)))
                from .association import getAllTests
                all_tests = getAllTests()
                nAll = len(all_tests)
                for idx, (test, obj) in enumerate(all_tests):
                    if args.limit is not None and idx == args.limit:
                        break
                    if args.verbosity == '0':
                        print(test)
                    else:
                        print('\n'.join(textwrap.wrap(
                            '{:<23} {}'.format(test, '' if obj.__doc__ is None else obj.__doc__),
                            subsequent_indent=' '*24, width=78)))
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
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
                print('Name:          {}'.format(args.items[0]))
                print('Description:   {}'.format('\n'.join(textwrap.wrap(test.__doc__, initial_indent='',
                        subsequent_indent=' '*15))))
                # create an instance of the test and pass -h to it
                test(1, ['-h']) 
            elif args.type == 'runtime_options':
                for idx, (opt, (def_value, description)) in enumerate(sorted(env.persistent_options.iteritems())):
                    if args.limit is not None and idx == args.limit:
                        break
                    # get the raw value of option (not the attribute, which might not be a string)
                    if args.verbosity == '0':
                        print(opt)
                    else:
                        val = str(getattr(env, '_' + opt))
                        print('{:<23} {} {}'.format(opt, val,
                            '(default)' if val == str(def_value) else '(default: {})'.format(def_value)))
                        print('\n'.join(textwrap.wrap(description, width=78,
                            initial_indent=' '*24, subsequent_indent=' '*24)))
                nAll = len(env.persistent_options)
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'runtime_option':
                if len(args.items) == 0:
                    raise ValueError('Please specify name of a runtime option')
                elif len(args.items) > 1:
                    raise ValueError('Please specify only one runtime option')
                print(getattr(env, '_' + args.items[0]))
            elif args.type == 'snapshot':
                if not args.items:
                    raise ValueError('Please provide a snapshot name or filename')
                elif len(args.items) > 1:
                    raise ValueError('Please specify only one snapshot')
                if args.items[0].startswith('vt_'): # only snapshot
                    source = 'online'
                    name = None
                    res = ResourceManager()
                    res.getRemoteManifest()
                    res.selectFiles(resource_type='snapshot')
                    for ss, prop in res.manifest.iteritems():
                        if ss[9:-7].lower() == args.items[0].lower():
                            name = ss[9:-7]
                            date = ''
                            desc = prop[3]
                            sz = prop[0]
                    if not name:
                        raise ValueError('Cannot locate an online snapshot with '
                            'name "{}"'.format(args.items[0]))
                else:
                    source = 'local'
                    name, date, desc = getSnapshotInfo(args.items[0])
                    sz = os.path.getsize(args.items[0])
                print('{:<23} {}'.format('Name:', name))
                print('{:<23} {}'.format('Source:', source))
                print('{:<23} {} ({})'.format('Size:', sz, sizeExpr(sz)))
                print('{:<23} {}'.format('Creation date:', date))
                print('\n'.join(textwrap.wrap('{:<23} {}'.format('Description:',
                    desc),  width=78, subsequent_indent=' '*24)))
            elif args.type == 'snapshots':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show snapshots"'.format(', '.join(args.items)))
                snapshots = []
                for snapshot_file in glob.glob(os.path.join(env.cache_dir, 'snapshot_*.tar')):
                    name, date, desc = getSnapshotInfo(snapshot_file)
                    sz = os.path.getsize(snapshot_file)
                    if name is not None:
                        snapshots.append((name, date, desc, sz))
                #
                for idx, (name, date, desc, sz) in enumerate(sorted(snapshots)):
                    if args.limit is not None and idx == args.limit:
                        break
                    if args.verbosity == '0':
                        print(name)
                    else:
                        text = '{:<23} {} ({}, created: {})'.format(name, desc,
                            sizeExpr(sz), date)
                        print('\n'.join(textwrap.wrap(text, width=78,
                            subsequent_indent=' '*24)))
                #
                nLocal = len(snapshots)
                res = ResourceManager()
                res.getRemoteManifest()
                res.selectFiles(resource_type='snapshot')
                for idx, (ss, prop) in enumerate(res.manifest.iteritems()):
                    if args.limit is not None and nLocal + idx >= args.limit:
                        break
                    if args.verbosity == '0':
                        print(ss[9:-7])
                    else:
                        text = '{:<23} {} ({}, online snapshot)'.format(ss[9:-7],
                            prop[3], sizeExpr(prop[0]))
                        print('\n'.join(textwrap.wrap(text, width=78,
                            subsequent_indent=' '*24)))
                nAll = nLocal + idx
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'pipelines':
                if args.items:
                    raise ValueError('Invalid parameter "{}" for command "vtools show pipelines"'
                        .format(', '.join(args.items)))
                res = ResourceManager()
                res.getRemoteManifest()
                res.selectFiles(resource_type='pipeline')
                nAll = len(res.manifest)
                for idx, (pipeline, prop) in enumerate(res.manifest.iteritems()):
                    if args.limit is not None and idx >= args.limit:
                        break
                    if args.verbosity == '0':
                        print(pipeline[9:-9])
                    else:
                        text = '{:<23} {}'.format(pipeline[9:-9], prop[3])
                        print('\n'.join(textwrap.wrap(text, width=78,
                            subsequent_indent=' '*24)))
                if args.limit is not None and args.limit >= 0 and args.limit < nAll:
                    print (omitted.format(nAll - args.limit))
            elif args.type == 'pipeline':
                if not args.items:
                    raise ValueError('Please specify a pipeline to display')
                elif len(args.items) > 1:
                    raise ValueError('Please specify only one pipeline to display')
                try:
                    pipeline = PipelineDescription(args.items[0])
                except Exception as e:
                    raise IndexError('Unrecognized pipeline {}: {}'
                        .format(args.items[0], e))
                pipeline.describe()
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


def adminArguments(parser):
    resource = parser.add_argument_group('Download or update resources')
    resource.add_argument('--update_resource', nargs='?', metavar='TYPE', 
        const='current', 
        choices=['current', 'all', 'existing', 'hg18', 'hg19', 'annotation', 
            'format', 'snapshot', 'pipeline'],
        help='''Download resources of specified type, which can be 'current' 
            (latest version of all resources excluding snapshots), 'all' (all 
            resources including obsolete databases), 'existing' (only update
            resources that exist locally), 'hg18' or 'hg19' (all resources for
            reference genome hg18 or hg19), 'annotation' (all current annotation
            databases), 'format' (all formats), and 'snapshot' (all online
            snapshots). Identical resources that are available locally (under
            ~/.variant_tools or runtime option $local_resource) are ignored.
            Note that option 'all' will download all versions of annotation 
            databases which can be slow and take a lot of disk spaces.''')
    merge = parser.add_argument_group('Merge samples')
    merge.add_argument('--merge_samples', action='store_true',
        help='''Merge samples with the same sample names by combining genotypes
        belonging to these samples. Phenotypes related to individual samples will
        be merged.''')
    rename_samples = parser.add_argument_group('Rename samples')
    rename_samples.add_argument('--rename_samples', nargs='+', metavar='COND',
        help='''This argument takes a condition by which samples are selected,
            followed by either a new sample name (assign a new name to selected
            samples) or an OLDNAME NEWNAME pair of patterns for which the first
            instance of OLDNAME in sample names will be replaced by NEWNAME.''')
    rename_table = parser.add_argument_group('Rename/describe tables')
    rename_table.add_argument('--rename_table', nargs=2, metavar=('NAME', 'NEW_NAME'),
        help='''Change table NAME to a NEW_NAME.''')
    rename_table.add_argument('--describe_table', nargs=2, metavar=('TABLE', 'NEW_DESCRIPTION'),
        help='''Update description for TABLE with a NEW_DESCRIPTION.''')
    validate = parser.add_argument_group('Validate reference genome')
    validate.add_argument('--validate_build', action='store_true',
        help='''Validate the build of project\'s reference genome by checking if
            the reference alleles of variants agree with the reference genome 
            sequence. A reference genome will be automatically downloaded if it
            does not exist in the local resource directory.''')
    snapshots = parser.add_argument_group('Save and load snapshots')
    snapshots.add_argument('--save_snapshot', nargs=2, metavar=('NAME', 'MESSAGE'),
        help='''Create a snapshot of the current project with NAME, which could be
            re-loaded using command 'vtools admin --load_snapshot'. A filename with
            extension .tar, .tgz or .tar.gz can be used to save the snapshot to a specific
            directory with compression but such snapshots are not listed by command
            'vtools show snapshots'. ''')
    snapshots.add_argument('--extra_files', nargs='*', metavar='FILE',
        help='''Additional files that will be saved along with the project and genotype
            databases. This could include customized format files, project-specific
            annotations, and results. Files outside of the current project directory
            are not allowed due to security considerations.''')
    snapshots.add_argument('--load_snapshot', metavar='NAME',
        help='''Revert the current project to specified snapshot. All changes since
        the that snapshot will be lost. The NAME should be one of the project snapshots
        or online snapshots listed by command 'vtools show snapshots', or name of a
        local snapshot file (with extension .tar, .tgz or .tar.gz).''')
    options = parser.add_argument_group('Set values for some various internal options.')
    options.add_argument('--set_runtime_option', nargs='+', metavar='OPTION',
        help='''Set value to internal options such as the batch size for database
        options. The default values of these options were chosen to fit most usage
        patterns but tweaking them might provide better performance under certain
        circumstances. Please use command "vtools show runtime_options" to list
        all currently supported runtime options.''')
    options.add_argument('--reset_runtime_option', metavar='OPT',
        help='''Reset value to a runtime option to its default value.''')


def admin(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if args.rename_samples:
                if len(args.rename_samples) not in [2, 3]:
                    raise ValueError('Option --rename_samples accept either '
                        'a new name or a pair of oldname newname patterns.')
                proj.renameSamples(args.rename_samples[0],
                    args.rename_samples[1].strip(),
                    None if len(args.rename_samples) == 2 else args.rename_samples[2].strip())
            elif args.merge_samples:
                proj.mergeSamples()
            elif args.rename_table:
                if args.rename_table[0] == 'variant':
                    raise ValueError('Cannot rename the master variant table')
                if args.rename_table[1] == 'variant':
                    raise ValueError('Cannot rename a table to the master variant table')
                if encodeTableName(args.rename_table[0]) not in proj.getVariantTables():
                    raise ValueError('Table {} does no exist or is not a variant table.'
                        .format(args.rename_table[0]))
                if encodeTableName(args.rename_table[1]) in proj.db.tables():
                    raise ValueError('Table {} already exists in the project'
                        .format(args.rename_table[1]))
                if args.rename_table[0] == args.rename_table[1]:
                    raise ValueError('Cannot rename a table to itself.')
                proj.db.renameTable(encodeTableName(args.rename_table[0]),
                    encodeTableName(args.rename_table[1]))
                env.logger.info('Table "{}" is renamed to "{}"'.format(args.rename_table[0],
                    args.rename_table[1]))
                # change the meta information of the table
                cur = proj.db.cursor()
                for key in ('desc', 'date', 'cmd'):
                    cur.execute('UPDATE project SET name="__{}_of_{}" WHERE name="__{}_of_{}"'
                        .format(key, encodeTableName(args.rename_table[1]),
                            key, encodeTableName(args.rename_table[0])))
            elif args.describe_table:
                if not proj.db.hasTable(encodeTableName(args.describe_table[0])):
                    raise ValueError('Table {} does not exist'.format(args.describe_table[0]))
                proj.describeTable(encodeTableName(args.describe_table[0]), args.describe_table[1])
                env.logger.info('Description of table {} is updated'
                    .format(args.describe_table[0]))
            elif args.validate_build:
                try:
                    refgenome = RefGenome(proj.build)
                except Exception as e:
                    raise RuntimeError('Failed to obtain reference genome for build {}: {}'
                        .format(proj.build, e))
                #
                cur = proj.db.cursor()
                prog = ProgressBar('Validate reference alleles',
                    proj.db.numOfRows('variant', exact=False))
                cur.execute('SELECT chr, pos, ref FROM variant WHERE ref != "-";')
                count = 0
                err_count = 0
                for chr, pos, ref in cur:
                    count += 1
                    if not refgenome.verify(chr, pos, ref):
                        err_count += 1
                        env.logger.debug('Ref allele mismatch: chr={}, pos={}, ref={}'
                            .format(chr, pos, ref))
                    prog.update(count + 1, err_count)
                prog.done()
                env.logger.info('{} non-insertion variants are checked. {} mismatch variants found.'
                    .format(count, err_count))
            elif args.set_runtime_option is not None:
                for option in args.set_runtime_option:
                    if '=' not in option:
                        raise ValueError('Runtime option should be specified as opt=value')
                    opt, value = option.split('=', 1)
                    if opt not in env.persistent_options:
                        raise ValueError('Only options {} are currently supported.'
                            .format(', '.join(env.persistent_options)))
                    try:
                        setattr(env, opt, value)
                    except Exception as e:
                        raise ValueError('Failed to set value {} to runtime options {}: {}'
                            .format(value, opt, e))
                    proj.saveProperty('__option_{}'.format(opt), getattr(env, opt))
                    env.logger.info('Option {} is set to {}'.format(opt, getattr(env, opt)))
            elif args.reset_runtime_option is not None:
                if args.reset_runtime_option not in env.persistent_options:
                    raise ValueError('Option {} is not a valid runtime option. '
                        'Use "vtools show runtime_options" to list currently '
                        'supported runtime options.'.format(args.reset_runtime_option))
                proj.removeProperty('__option_{}'.format(args.reset_runtime_option))
                env.logger.info('Option {} is set to its default value'.format(args.reset_runtime_option))
            elif args.save_snapshot is not None:
                if args.extra_files is not None:
                    cur_dir = os.path.realpath(os.getcwd())
                    for f in args.extra_files:
                        if f == '.snapshot.info':
                            raise ValueError('Cannot include ".snapshot.info" '
                                'in snapshot due to filename conflicts. '
                                'Please rename the file.')
                        if os.path.isdir(f):
                            raise ValueError('Cannot add directory "{0}" into '
                                'snapshot. You should use wildcard names (e.g., '
                                '{0}/*") if you want to save all files under '
                                'this directory.'.format(f))
                        if not os.path.isfile(f):
                            raise ValueError('Cannot include {} in snapshot. File does not exist.'.format(f))
                        # if the file is not under the current directory
                        if os.path.commonprefix([cur_dir, os.path.realpath(f)]) != cur_dir:
                            raise ValueError('Only files under the current project directory could be included in a snapshot.')
                        # if the file is automatically included
                        if f == proj.name + '.proj':
                            raise ValueError('Project database is already included.')
                        if f == proj.name + '_genotype.proj':
                            raise ValueError('Project genotype database is already included.')
                proj.saveSnapshot(args.save_snapshot[0], args.save_snapshot[1], args.extra_files)
                env.logger.info('Snapshot {} has been saved'.format(args.save_snapshot[0]))
            elif args.load_snapshot is not None:
                proj.loadSnapshot(args.load_snapshot)
                env.logger.info('Snapshot {} has been loaded'.format(args.load_snapshot))
            elif args.update_resource:
                res = ResourceManager()
                res.getRemoteManifest()
                res.selectFiles(resource_type=args.update_resource)
                res.excludeExistingLocalFiles()
                env.logger.info('{} files need to be downloaded or updated'.format(len(res.manifest)))
                res.downloadResources()
            else:
                env.logger.warning('Please specify an operation. Type `vtools admin -h\' for available options')
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)

if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools
    pass
