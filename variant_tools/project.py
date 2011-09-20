#!/usr/bin/env python
#
# $File: project.py $
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
import glob
import logging
import getpass
import random
import textwrap
import tempfile
import shutil
import ConfigParser
from collections import namedtuple, defaultdict
from .utils import DatabaseEngine, ProgressBar, setOptions, SQL_KEYWORDS, delayedAction, filesInURL, downloadFile

VTOOLS_VERSION = '1.0beta'
VTOOLS_COPYRIGHT = '''variant tools version {} : Copyright (c) 2011 Bo Peng.'''.format(VTOOLS_VERSION)
VTOOLS_CITE = '''Please cite Anthony et al ....''' # pending
VTOOLS_CONTACT = '''Please visit http://varianttools.sourceforge.net for more information.'''

# define a field type
Field = namedtuple('Field', ['name', 'index', 'adj', 'type', 'comment'])
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
    def __init__(self, proj, annoDB, linked_by=[]):
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
            self.version = '-'.join(self.filename.split('-')[1:])
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
                self.anno_type = rec[1]
            elif rec[0] == 'build':
                self.refGenomes = eval(rec[1])
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
        if self.anno_type == 'attribute' and len(self.linked_by) != len(self.build):
            raise RuntimeError('Please specify link fields for attributes {} using parameter --linked_by'.format(','.join(self.build)))
        if self.linked_by:
            s = delayedAction(proj.logger.info, 'Indexing linked field {}'.format(', '.join(self.linked_by)))
            self.indexLinkedField(proj, linked_by)
            del s

    def indexLinkedField(self, proj, linked_fields):
        '''Create index for fields that are linked by'''
        cur = proj.db.cursor()
        for linked_field in linked_fields:
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
    def __init__(self, name):
        '''Input file format'''
        # locate a file format specification file
        self.description = None
        #
        if os.path.isfile(name + '.fmt'):
            self.name = os.path.split(name)[-1]
            self.parseFMT(name + '.fmt') 
        elif name.endswith('.fmt') and os.path.isfile(name):
            self.name = os.path.split(name)[-1][:-4]
            self.parseFMT(name) 
        else:
            url = 'http://vtools.houstonbioinformatics.org/input_fmt/{}.fmt'.format(name)
            try:
                fmt = downloadFile(url, '.')
            except Exception as e:
                raise ValueError('Failed to download format specification file {}.fmt'.format(name))
            self.name = name
            self.parseFMT(fmt)

    def parseFMT(self, filename):
        parser = ConfigParser.SafeConfigParser()
        parser.read(filename)
        # sections?
        sections = parser.sections()
        if 'format description' not in sections:
            raise ValueError("Missing section 'format description'")
        #
        fields = []
        for section in sections:
            if section == 'format description':
                continue
            try:
                items = [x[0] for x in parser.items(section, raw=True)]
                fields.append(
                    Field(name=section,
                        index=parser.get(section, 'index', raw=True),
                        type=parser.get(section, 'type', raw=True),
                        adj=parser.get(section, 'adj', raw=True) if 'adj' in items else None,
                        comment=parser.get(section, 'comment', raw=True) if 'comment' in items else '')
                    )
            except Exception as e:
                raise ValueError('Invalid section {} in configuration file {}'.format(section, self.name))
        #
        if len(fields) == 0:
            raise ValueError('No valid field is defined in format specification file {}'.format(self.name))
        #
        self.delimiter = '\t'
        variant_fields = []
        position_fields = []
        range_fields = []
        variant_info = []
        genotype_fields = []
        genotype_info = []
        for item in parser.items('format description'):
            if item[0] == 'description':
                self.description = item[1]
            if item[0] == 'delimiter':
                self.delimiter = eval(item[1])
            if item[0] in ['variant_fields', 'position_fields', 'range_fields', 'genotype_fields', 'variant_info', 'genotype_info']:
                linked_fields = eval(item[0])
                for name in map(str.strip, item[1].split(',')):
                    fld = [x for x in fields if x.name == name]
                    if len(fld) != 1:
                        raise ValueError('Cannot find field "{}" in the format specification file, which defines fields {}.'\
                            .format(name, ', '.join([x.name for x in fields])))
                    linked_fields.append(fld[0])
        #
        # Post process all fields
        if (len(variant_fields) != 0) + (len(position_fields) != 0) + (len(range_fields) != 0) != 1:
            raise ValueError('Please specify one and only one of variant_fields, position_fields and range_fields')
        #
        if variant_fields:
            self.input_type = 'variant'
            self.ranges = [0, 4]
            self.fields = variant_fields
            if len(self.fields) != 4:
                raise ValueError('variant_fields should have four fields for chr, pos, ref, and alt alleles')
        elif position_fields:
            self.input_type = 'position'
            self.ranges = [0, 2]
            self.fields = position_fields
            if len(self.fields) != 2:
                raise ValueError('position_fields should have two fields for chr and pos')
        elif range_fields:
            self.input_type = 'range'
            self.ranges = [0, 3]
            self.fields = range_fields
            if len(self.fields) != 3:
                raise ValueError('range_fields should have three fields for chr and starting and ending position')
        #
        if self.input_type != 'variant' and not variant_info:
            raise ValueError('Input file with type position or range must specify variant_info')
        if self.input_type != 'variant' and genotype_info:
            raise ValueError('Input file with type position or range can not have any genotype information.')
        if genotype_fields and len(genotype_fields) != 1:
            raise ValueError('Variant tools currently only support input file with at most one sample')
        # extend self.fields to include variant_info
        self.fields.extend(variant_info)
        self.ranges.append(self.ranges[-1] + len(variant_info))
        self.fields.extend(genotype_fields)
        self.ranges.append(self.ranges[-1] + len(genotype_fields))
        self.fields.extend(genotype_info)
        self.ranges.append(self.ranges[-1] + len(genotype_info))
                

    def describe(self):
        print('Format:      {}'.format(self.name))
        if self.description is not None:
            print('Description: {}'.format('\n'.join(textwrap.wrap(self.description,
                initial_indent='', subsequent_indent=' '*4))))
        #
        print('\n{} fields:'.format(self.input_type))
        for fld in self.fields[self.ranges[0]:self.ranges[1]]:
            print('  {}:    {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                subsequent_indent=' '*8))))
        if self.ranges[1] != self.ranges[2]:
            print('\nVariant information fields:')
            for fld in self.fields[self.ranges[1]:self.ranges[2]]:
                print('   {}:    {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*8))))
        if self.ranges[2] != self.ranges[3]:
            print('\nGenotype field:')
            for fld in self.fields[self.ranges[2]:self.ranges[3]]:
                print('   {}:    {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*8))))
        if self.ranges[3] != self.ranges[4]:
            print('\nGenotype information fields:')
            for fld in self.fields[self.ranges[3]:self.ranges[4]]:
                print('   {}:    {}'.format(fld.name, '\n'.join(textwrap.wrap(fld.comment,
                    subsequent_indent=' '*8))))



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

    2. Table "sample_variant_$sample_id" stores variants of each sample
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
    def __init__(self, name=None, new=False, verbosity=None, **kwargs):
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
        # set global verbosity level
        setOptions(verbosity=verbosity)
        #
        # create a logger
        self.logger = logging.getLogger()
        self.logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        levels = {
            None: logging.INFO,
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
        cout.setLevel(levels[self.verbosity])
        cout.setFormatter(logging.Formatter('%(levelname)s: %(message)s'))
        self.logger.addHandler(cout)
        # output to a log file
        ch = logging.FileHandler(self.name + '.log', mode='w' if new else 'a')
        # NOTE: debug informaiton is always written to the log file
        ch.setLevel(logging.DEBUG)
        ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
        self.logger.addHandler(ch)
        if new:
            self.create(**kwargs)
        else:
            self.open()
        #
        # create a temporary directory
        self.temp_dir = tempfile.mkdtemp()

    def create(self, **kwargs):
        '''Create a new project'''
        # open the project file
        self.logger.info(VTOOLS_COPYRIGHT)
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
        self.build = None
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

    def open(self):
        '''Open an existing project'''
        # open the project file
        self.logger.info('Opening project {}'.format(self.proj_file))
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
                linked_by = eval(self.loadProperty('{}_linked_by'.format(os.path.split(db)[-1].split('-')[0]), default='[]'))
                self.annoDB.append(AnnoDB(self, db, linked_by))
                self.db.attach(db)
            except Exception as e:
                self.logger.warning(e)
                self.logger.warning('Cannot open annotation database {}'.format(db))
        #
        # get existing meta information
        # FIXME: these are not handled correctly for now
        self.variant_meta = self.db.getHeaders('variant_meta')
        self.sample_meta = self.db.getHeaders('sample_meta')

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

    def close(self):
        '''Write everything to disk...'''
        # This is no longer needed once we create temporary table outside
        # of the project database
        #self.removeTempTables()
        self.db.commit()
        if os.path.isdir(self.temp_dir):
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
                name VARCHAR(40) NOT NULL,
                value VARCHAR(256) NULL
            )''')
        # create index
        try:
            cur.execute('''CREATE UNIQUE INDEX project_index ON name (name ASC);''')
        except Exception as e:
            # the index might already exists
            return

    def createFilenameTable(self):
        self.logger.debug('Creating table filename')
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

    def createMasterVariantTable(self):
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
                alt VARCHAR(255) NOT NULL);'''.format(self.db.AI))
        self.createIndexOnMasterVariantTable()

    def createIndexOnMasterVariantTable(self):
        # create indexes
        #
        s = delayedAction(self.logger.info, 'Creating index on master variant table. This might take quite a while.')
        try:
            self.db.execute('''CREATE UNIQUE INDEX variant_index ON variant (bin ASC, chr ASC, pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists
            self.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
        try:
            if self.alt_build:
                self.db.execute('''CREATE UNIQUE INDEX variant_alt_index ON variant (alt_bin ASC, alt_chr ASC, alt_pos ASC, ref ASC, alt ASC);''')
        except Exception as e:
            # the index might already exists
            self.logger.debug(e)
        # the message will not be displayed if index is created within 5 seconds
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

    def createNewSampleVariantTable(self, table, genotype=True, fields=[]):
        '''Create a table ``sample_variant`` to store vcf data'''
        if self.db.hasTable(table):
            self.logger.info('Removing existing table {}, which can be slow for sqlite databases'.format(table))
            self.db.removeTable(table)
        cur = self.db.cursor()
        # self.logger.debug('Creating table {}'.format(table))
        #
        cur.execute('''\
            CREATE TABLE IF NOT EXISTS {0} (
                variant_id INT NOT NULL
            '''.format(table) + 
            (', variant_type INT' if genotype else '') + 
            ''.join([', {} FLOAT NULL'.format(f) for f in fields]) + ');'
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
        return self.db.hasTable(table) and self.db.getHeaders(table)[0] == 'variant_id'
        
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
        self.logger.info('Removing table {} (this might take a while)'.format(table))
        self.db.removeTable(table)

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
            self.logger.debug(e)
            raise ValueError('Failed to retrieve samples by condition "{}"'.format(cond))

    def removeSample(self, ID):
        '''Remove sample and their genotype, but not variants'''
        cur = self.db.cursor()
        cur.execute('SELECT filename.filename, sample.sample_name FROM sample LEFT OUTER JOIN filename ON sample.file_id = filename.file_id WHERE sample.sample_id = {};'.format(self.db.PH), (ID,))
        res = cur.fetchone()
        self.logger.info('Removing sample {} from file {}'.format(res[1], res[0]))
        cur.execute('DELETE FROM sample WHERE sample_id = {};'.format(self.db.PH), (ID,))
        self.db.commit()
        self.db.removeTable('{}_genotype.sample_variant_{}'.format(self.name, ID))        
        
    def summarize(self):
        '''Summarize key features of the project
        '''
        # FIXME: more summary
        info =  'Project name:                {}\n'.format(self.name)
        info += 'Primary reference genome:    {}\n'.format(self.build)
        info += 'Secondary reference genome:  {}\n'.format(self.alt_build)
        info += 'Database engine:             {}\n'.format(self.db.engine)
        info += 'Variant tables:              {}\n'.format(', '.join(self.getVariantTables()))
        info += 'Annotation databases:        {}\n'.format(', '.join([os.path.join(x.dir, x.name) \
            + (' ({})'.format(x.version) if x.version else '') for x in self.annoDB]))
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
        for table in self.getVariantTables() + ['sample', 'filename']:
            if table == exclude:
                continue
            if name.lower() in [x.lower() for x in self.db.getHeaders(table)]:
                raise ValueError("Field name '{}' is not allowed because it is already used in table {}".format(name, table))
        
    #def removeTempTables(self):
    #    '''Remove all temp tables'''
    #    for table in self.db.tables():
    #        if table.startswith('_tmp_'):
    #            self.db.removeTable(table)

    #
    # Handling field query
    #
    def linkFieldToTable(self, field, variant_table):
        '''Return one or more FieldConnections that link a field to a variant_table'''
        # if field is specified by table.field, good
        if '.' in field:
            # possibly two dots (db.table.field), but never seen them.
            table, field = '.'.join(field.split('.')[:-1]).lower(), field.split('.')[-1].lower()
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
                if db.anno_type == 'attribute':
                    return sum([self.linkFieldToTable(x, variant_table) for x in db.linked_by], []) + [
                        FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            link= ' AND '.join(['{}.{}={}'.format(table, x, y) for x,y in zip(db.build, db.linked_by)]))]
                if db.build is not None:
                    if db.anno_type == 'position':  # chr and pos
                        return self.linkFieldToTable('{}.variant_id'.format(varinat_table), 'variant') + [
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
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
                            FieldConnection(
                            field= '{}.{}'.format(table, field),
                            table= '{}.{}'.format(table, table),
                            # FIXME: how to use bin here?
                            link= 'variant.chr = {0}.{1} AND variant.pos >= {0}.{2} AND variant.pos <= {0}.{3}'
                                    .format(table, db.build[0], db.build[1], db.build[2]))]
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
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
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
                if db.anno_type == 'attribute':
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
                        # FIXME: how to use bin?
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
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
                        # FIXME: how to use bin here?
                        return self.linkFieldToTable('{}.variant_id'.format(variant_table), 'variant') + [
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
        
        
#
#
# Functions provided by this script
#
#

def initArguments(parser):
    parser.add_argument('project',
        help='''Name of a new project. This will create a new file $project.proj under
            the current directory. Only one project is allowed in a directory.''')
    parser.add_argument('-f', '--force', action='store_true',
        help='''Remove a project if it already exists.''')
    parser.add_argument('--engine', choices=['mysql', 'sqlite3'], default='sqlite3',
        help='''Database engine to use, can be mysql or sqlite3. Parameters --host, --user
            and --passwd will be needed for the creation of a new mysql project.''')
    parser.add_argument('--host', default='localhost', 
        help='The MySQL server that hosts the project databases.')
    parser.add_argument('--user', default=getpass.getuser(),
        help='User name to the MySQL server. Default to current username.')
    parser.add_argument('--passwd',
        help='Password to the MySQL server.')
    parser.add_argument('--batch', default=10000, 
        help='Number of query per transaction. Larger number leads to better performance but requires more ram.')


def init(args):
    try:
        if args.force:
            # silently remove all exiting project.
            try:
                proj = Project(verbosity='0')
                proj.remove()
            except:
                pass
        # create a new project
        proj = Project(name=args.project, new=True, verbosity=args.verbosity,
            engine=args.engine, host=args.host, user=args.user, passwd=args.passwd,
            batch=args.batch)
        proj.close()
    except Exception as e:
        sys.exit(e)


def removeArguments(parser):
    parser.add_argument('type', choices=['project', 'table', 'samples', 'field', 'annotation'],
        help='''Type of items to be removed.''')
    parser.add_argument('items', nargs='*',
        help='''Items to be removed. It can be the name of project for type
            project (optional), names of one or more variant tables for
            type table, a pattern for type 'samples', name of fields or annotation
            databases. Note that removal of samples will only remove sample information
            related to variants, not variants themselves; removal of annotation databases
            will stop using these databases in the project, but will not removing them
            from disk.''')
    

def remove(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            if args.type == 'project':
                if len(args.items) > 0 and args.items[0] != proj.name:
                    raise ValueError('Cannot remove project: Incorrect project name')
                proj.remove()
            elif args.type == 'table':
                for table in args.items:
                    proj.removeVariantTable(table)
            elif args.type == 'samples':
                if len(args.items) == 0:
                    raise ValueError('Please specify conditions to select samples to be removed')
                proj.db.attach(proj.name + '_genotype')
                IDs = proj.selectSampleByPhenotype(' AND '.join(args.items))
                if len(IDs) == 0:
                    proj.logger.warning('No sample is selected by condition {}'.format(' AND '.join(args.items)))
                for ID in IDs:
                    proj.removeSample(ID)
            elif args.type == 'field':
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
            elif args.type == 'annotation':
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
    except Exception as e:
        sys.exit(e)


def showArguments(parser):
    parser.add_argument('type', choices=['project', 'tables', 'table',
        'samples', 'fields', 'annotations', 'annotation', 'formats', 'format'],
        nargs='?', default='project',
        help='''Type of information to display, which can be project (summary
            of a project, tables (all variant tables, or all tables if
            verbosity=2), table (a specific table), samples (sample and
            phenotype information), fields (from variant tables and all used
            annotation databases), annotations (all available annotation
            databases for variant tools), specified annotation, all supported
            input formats, and details of specific format. Default to
            project.''')
    parser.add_argument('items', nargs='*',
        help='''Items to display, which can be name of a table for type 'table'.''')
    parser.add_argument('-l', '--limit', default=10, type=int,
        help='''Number of record to display for option 'show table'. All
            records will be displayed if it is set to -1''')


def show(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            #
            if args.type == 'project':
                print(proj.summarize())
            elif args.type == 'tables':
                print('{:<20} {}'.format('table', '#variants'))
                for table in proj.db.tables() if args.verbosity=='2' else proj.getVariantTables():
                    print('{:<20} {:,}'.format(table, proj.db.numOfRows(table)))
            elif args.type == 'table':
                if not args.items:
                    raise ValueError('Please specify a table to display')
                table = args.items[0]
                # showing an annotation database
                if table in [x.name for x in proj.annoDB]:
                    table = '{0}.{0}'.format(table)
                if not proj.db.hasTable(table):
                    raise ValueError('Table {} does not exist'.format(table))
                # print content of table
                headers = proj.db.getHeaders(table)
                print(', '.join(headers))
                cur = proj.db.cursor()
                if args.limit < 0:
                    cur.execute('SELECT * FROM {};'.format(table, 10))
                else:
                    cur.execute('SELECT * FROM {} LIMIT 0,{};'.format(table, args.limit))
                for rec in cur:
                    print(', '.join([str(x) for x in rec]))
            elif args.type == 'samples':
                if not proj.db.hasTable('sample'):
                    proj.logger.warning('Project does not have a sample table.')
                    return
                cur = proj.db.cursor()
                fields = proj.db.getHeaders('sample')
                # headers are ID, file, sample, FIELDS
                print('filename\tsample_name{}'.format(''.join(['\t'+x for x in fields[3:]])))
                cur.execute('SELECT filename, {} FROM sample, filename WHERE sample.file_id = filename.file_id;'\
                    .format(', '.join(fields[2:])))
                for rec in cur:
                    print('\t'.join(['{}'.format(x) for x in rec]))
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
                DBs = filesInURL('http://vtools.houstonbioinformatics.org/annoDB', ext='.ann')
                for db in DBs:
                    print(db)
            elif args.type == 'formats':
                FMTs = filesInURL('http://vtools.houstonbioinformatics.org/input_fmt', ext='.fmt')
                for fmt in FMTs:
                    print(fmt)
            elif args.type == 'format':
                for item in args.items:
                    try:
                        fmt = fileFMT(item)
                    except Exception as e:
                        proj.logger.debug(e)
                        raise IndexError('Input file format {} is not currently supported by variant tools'.format(item))
                    fmt.describe()
    except Exception as e:
        sys.exit(e)


def executeArguments(parser):
    parser.add_argument('query', nargs='*', help='A SQL query to execute')


def execute(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            cur = proj.db.cursor()
            query = ' '.join(args.query)
            proj.logger.debug('Executing SQL statement: "{}"'.format(query))
            cur.execute(query)
            proj.db.commit()
            for rec in cur:
                print(', '.join(['{}'.format(x) for x in rec]))
    except Exception as e:
        sys.exit(e)


if __name__ == '__main__':
    # for testing purposes only. The main interface is provided in vtools.py
    pass

