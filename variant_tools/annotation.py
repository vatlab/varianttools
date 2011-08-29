#!/usr/bin/env python
#
# $File: annotation.py $
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
import os
import ConfigParser
import shutil
import urlparse
import gzip
import zipfile
import tempfile
import re

from project import AnnoDB, Project, Field
from utils import ProgressBar, downloadFile, lineCount, DatabaseEngine, getMaxUcscBin, delayedAction, decompressIfNeeded, normalizeVariant


class AnnoDBConfiger:
    '''An annotation database can be created from either a configuration file
    or a database name. In the former case, the configuration file will be
    parsed to determine the reference genome and fields the database provides.
    The database will not be download or imported unless one of the fields
    is used. In the later case, the information will be directly read from
    the annotation table in the database. Please check annotation/dbNSFP.ann
    for a description of the format of the configuration file.'''
    def __init__(self, proj, annoDB):
        '''Create an annotation database.
        proj: current project
        annoDB: annotation database. It can be either a mysql or sqlite database
            name depend on the database engine used in the project, or a
            configuration file. In the latter case, a filename is allowed.
        '''
        self.proj = proj
        self.logger = proj.logger
        #
        # data that will be available after parsing
        self.path = ''
        self.name = None
        # point annotation table
        self.anno_type = 'variant'
        self.build = []
        self.descrption = ''
        self.direct_url = None
        self.source_type = None
        self.source_url = None
        self.source_files = None
        self.fields = []
        self.version = None
        # where is annoDB, in which form?
        self.parseConfigFile(annoDB)
        # some fields have to be determined.
        if self.name is None:
            raise ValueError('No valid name is set')
        if self.anno_type not in ['position', 'variant', 'range', 'attribute']:
            raise ValueError('vtools only support position, variant, range, and attribute types')
        if len(self.fields) == 0:
            raise ValueError('No valid field is located from database {}'.format(annoDB))
        if len(self.build) == 0:
            raise ValueError('No reference genome is specified for database {}'.format(annoDB))
        if not proj.build is None:
            if (not self.build.keys()[0] == '*') and (not proj.build in self.build.keys()) and (proj.alt_build is None or \
                proj.alt_build not in self.build.keys()):
                raise ValueError('Annotation database cannot be used for the existing project.')

    def remove(self):
        '''Remove an annotation database'''
        self.logger.info('Removing database {}'.format(self.name))
        self.proj.db.removeDatabase(self.name)

    def parseConfigFile(self, filename):
        '''Read from an ini style configuration file'''
        self.logger.debug('Checking configuration file {}'.format(filename))
        self.path = os.path.split(filename)[0]
        # get filename, remove extension, and keep version as things after -.
        self.name = os.path.splitext(os.path.split(filename)[-1])[0]
        self.version = None
        if '-' in self.name:
            self.version = '-'.join(self.name.split('-')[1:])
            self.name = self.name.split('-')[0]
        parser = ConfigParser.SafeConfigParser()
        parser.read(filename)
        # sections?
        sections = parser.sections()
        if 'linked fields' not in sections:
            raise ValueError('Ignore invalid annotation file {}. Missing build section.'.format(filename))
        if 'data sources' not in sections:
            raise ValueError('Ignore invalid annotation file {}. Missing source section.'.format(filename))
        # linked fields for each reference genome
        try:
            self.build = {}
            for item in parser.items('linked fields'):
                fields = [x.strip() for x in item[1].split(',')]
                for field in fields:
                    if field not in sections:
                        self.logger.error('Invalid reference genome: Unspecified field {}.'.format(field))
                self.build[item[0]] = fields
            self.logger.debug('Reference genomes: {}'.format(self.build))
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to obtain build information of annotation file {}'.format(filename))
        # source information
        self.anno_type = 'variant'
        self.description = None
        self.direct_url = None
        self.source_url = None
        self.source_pattern = None
        self.source_type = None
        for item in parser.items('data sources'):
            if item[0] == 'direct_url':
                self.direct_url = item[1]
            elif item[0] == 'source_url':
                self.source_url = item[1]
            elif item[0] == 'source_pattern':
                self.source_pattern = item[1]
            elif item[0] == 'source_type':
                self.source_type = item[1]
            elif item[0] == 'description':
                self.description = item[1]
            elif item[0] == 'anno_type':
                self.anno_type = item[1]
            elif item[0] == 'version':
                if self.version is not None and self.version != item[1]:
                    raise Warning('Version obtained from filename ({}) is inconsistent with version specified in the annotation file ({})'.format(self.version, item[1]))
                self.version = item[1]
            else:
                raise ValueError('Invalid keyword {} in section "data sources".'.format(item[0]) + 
                    'Only direct_url, source_url, source_type, version, source_pattern and description are allowed')
        # sections
        self.fields = []
        for section in sections:
            if section == 'linked fields' or section == 'data sources':
                continue
            try:
                items = [x[0] for x in parser.items(section, raw=True)]
                # 'index' - 1 because the input is in 1-based index
                self.fields.append(Field(name=section, index=parser.get(section, 'index', raw=True),
                    type=parser.get(section, 'type', raw=True),
                    null=parser.get(section, 'null', raw=True) if 'null' in items else None,
                    comment=parser.get(section, 'comment', raw=True) if 'comment' in items else ''))
            except Exception as e:
                self.logger.debug(e)
                self.logger.error('Invalid section {} in configuration file {}'.format(section, filename))
                return False
        self.fields.sort(key=lambda k: k.index)
        for field in self.fields:
            self.logger.debug("Field {0}: index {1},\t{2}".format(field.name, field.index, field.type))

    def createAnnotationTable(self, db):
        'Create an annotation table '
        items = []
        for build in self.build.keys():
            if build != '*':
                items.append('{0}_bin INTEGER'.format(build))
        for field in self.fields:
            if 'chromosome' in field.type:
                items.append('{0} VARCHAR(20)'.format(field.name))
            elif 'position' in field.type:
                items.append('{0} INTEGER'.format(field.name))
            else:
                items.append('{0} {1}'.format(field.name, field.type))
        query = '''CREATE TABLE IF NOT EXISTS {} ('''.format(self.name) + \
            ',\n'.join(items) + ');'
        self.logger.debug('Creating annotation table {} using query\n{}'.format(self.name, query))
        cur = db.cursor()
        try:
            cur.execute(query)
        except Exception as e:
            self.logger.debug(e)
            raise ValueError('Failed to create table')
    
    def createFieldsTable(self, db):
        '''Create table name_fields'''
        cur = db.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS {}_field (
            name VARCHAR(40),
            field INT,
            type VARCHAR(80),
            null_val VARCHAR(20),
            comment VARCHAR(256)
        )'''.format(self.name))

    def createInfoTable(self, db):
        '''Create table name_fields'''
        cur = db.cursor()
        cur.execute('''CREATE TABLE IF NOT EXISTS {}_info (
            name VARCHAR(40),
            value VARCHAR(1024)
        )'''.format(self.name))

    def getSourceFiles(self, tdir):
        '''Download a file from a URL and save to a directory. If it is a zip file,
        unzip it in the directory.'''
        # FIXME: if we detect a file in the local directory, we should not try to 
        # download from online? However, users might want to get the latest 
        # version from online, and if they want to use a local copy, they can always
        # use option --files.
        tempFile = None
        try:
            if os.path.isfile(self.source_url):
                if not self.source_url.lower().endswith('.zip'):
                    return self.source_url
                else:
                    tempFile = self.source_url
            elif os.path.isfile(os.path.join(self.path, self.source_url)):
                if not self.source_url.lower().endswith('.zip'):
                    return os.path.join(self.path, self.source_url)
                else:
                    tempFile = os.path.join(self.path, self.source_url)
            else:
                filename = os.path.split(self.source_url)[-1]
                self.logger.info('Downloading {}'.format(filename))
                tempFile = downloadFile(self.source_url, tdir)
        except Exception as e:
            self.logger.error(e)
            raise ValueError('Failed to download database source from {}'.format(self.source_url))
        #
        if not os.path.isfile(tempFile):
            raise ValueError('Could not find downloaded database source')
        # regular file?
        if not tempFile.endswith('.zip'):
            return [tempFile]
        # if zip file?
        bundle = zipfile.ZipFile(tempFile)
        self.logger.info('Extracting zip file to a temporary directory {}'.format(tdir))
        bundle.extractall(tdir)
        return [os.path.join(tdir, name) for name in bundle.namelist()]
    
    def openAnnoFile(self, filename):
        if filename.lower().endswith('.gz'):
            return gzip.open(filename, 'rb')
        else:
            # text file
            return open(filename, 'r')

    def importTxtRecords(self, db, source_files):
        # First: Ucsc bins calculated for position fields
        build_info = []
        for key,items in self.build.items():
            if key == '*':
                continue
            try:
                # items have chr/pos, chr/pos/ref/alt, chr/start/end for different annotation types
                pos_idx, pos_adj = [(i, 1) if x.type == '0-based position' else (i, 0) for i, x in enumerate(self.fields) if x.name == items[1]][0]
                if self.anno_type == 'variant':
                    # save indexes for pos, ref and alt
                    ref_idx = [i for i,x in enumerate(self.fields) if x.name == items[2]][0]
                    alt_idx = [i for i,x in enumerate(self.fields) if x.name == items[3]][0]
                    build_info.append((pos_idx, pos_adj, ref_idx, alt_idx))
                else:
                    build_info.append((pos_idx, pos_adj))
            except Exception as e:
                self.logger.error('No field {} for build {}'.format(items[1], key))
        #
        field_info = []
        variant_fields = [0] if self.anno_type == 'variant' else None
        # Other fields
        for field in self.fields:
            # adjust for zero-based position
            if field.type == '0-based position':
                field_info.append((int(field.index), 1, field.null))
            elif field.type == 'chromosome':
                field_info.append((int(field.index), 'c', field.null))
            # if index is not pure digital, it should be in the form of
            # 20:DEPTH where 20 is index, : is separator and DEPTH is field name
            # we need to get the fields...
            elif not field.index.isdigit():
                try:
                    i, s, n = re.match('(\d+)([\s\t:;,|])(\w+)', field.index).groups()
                except:
                     # wrong index...
                     self.logger.debug("Regex is not working for {}".format(field))
                pattern = re.compile('(^|{0}){1}((=[^{0}]+)*)($|{0})'.format(s, n)) # patter of [separator or begin]name=value[separator or end] )
                field_info.append((int(i), (n, pattern), field.null))
            else:
                field_info.append((int(field.index), None, field.null))
        # files?
        cur = db.cursor()
        insert_query = 'INSERT INTO {0} VALUES ('.format(self.name) + \
                            ','.join([db.PH] * (len(self.fields) + len(build_info))) + ');'
        for f in source_files:
            self.logger.info('Importing annotation data from {0}'.format(f))
            all_records = 0
            skipped_records = 0
            prog = ProgressBar(os.path.split(f)[-1], lineCount(f))
            with self.openAnnoFile(f) as input_file:
                for line in input_file:
                    all_records += 1
                    try:
                        if line.startswith('#'):
                            continue
                        # get data
                        tokens = [x.strip() for x in line.split('\t')]
                        records = []
                        #
                        for col, adj, null in field_info:
                            col -= 1
                            item = tokens[col]
                            if adj is not None:
                                if adj == 1:
                                    item = int(item) + adj
                                elif adj == 'c':
                                    if item.startswith('chr'):  # if the item doesn't start with chr, then item=chromsome_number
                                        item = item[3:]
                                else: # must be a name,pattern pair
                                    try:
                                        # a flag that isn't listed (therefore set to False or 0)
                                        if not adj[0] in item:
                                            item = 0
                                        else:
                                            item = adj[1].search(item).group(3)
                                            #self.logger.debug("matched item = " + item)
                                            if item is not None:
                                                item = item[1:]
                                            else: # a flag that exists (therefore set to True or 1)
                                                item = 1
                                    except:
                                        # this is the case missed by 'adj[0] in item'. For example
                                        # Pilot1 is in item with field Pilot123
                                        item = 0
                            if item == null:
                                item = None
                            #
                            records.append(item)
                        #
                        # handle records
                        if not build_info:
                            # there is no build information, this is 'field' annotation, nothing to worry about
                            cur.execute(insert_query, records)
                        elif len(build_info[0]) == 2:
                            # there is no ref and alt, easy
                            bins = []
                            for pos_idx, pos_adj in build_info:
                                try:
                                    # zero-based: v, v+1 (adj=1)
                                    # one-based:  v-1, v (adj=0)
                                    bin = getMaxUcscBin(int(records[pos_idx]) + pos_adj - 1, int(records[pos_idx]) + pos_adj)
                                except:
                                    # position might be None (e.g. dbNSFP has complete hg18 coordinates,
                                    # but incomplete hg19 coordinates)
                                    bin = None
                                bins.append(bin)
                            cur.execute(insert_query, bins + records)
                        else:
                            # variant... most troublesome
                            input_ref = records[build_info[0][2]]
                            all_alt = records[build_info[0][3]]
                            # support multiple alternative alleles
                            for input_alt in all_alt.split(','):
                                bins = []
                                for pos_idx, pos_adj, ref_idx, alt_idx in build_info:
                                    input_pos = int(records[pos_idx]) + pos_adj
                                    bin, pos, ref, alt = normalizeVariant(input_pos, input_ref, input_alt)
                                    # these differ build by build
                                    bins.append(bin)
                                    records[pos_idx] = pos
                                    # these should be kept unchanged
                                    records[ref_idx] = ref
                                    records[alt_idx] = alt
                                cur.execute(insert_query, bins + records)    
                    except Exception as e:
                        # if any problem happens, just ignore
                        self.logger.debug(e)
                        skipped_records += 1
                    if all_records % db.batch == 0:
                        prog.update(all_records)
                        db.commit()
            db.commit()
            prog.done()
            self.logger.info('{0} records handled, {1} ignored.'\
                .format(all_records, skipped_records))

    def importFromSource(self, source_files):
        '''Importing data from source files'''
        tdir = None
        if source_files == []:
            tdir = tempfile.mkdtemp()
            source_files = self.getSourceFiles(tdir)
        elif len(source_files) == 1 and source_files[0].endswith('.zip'):
            # trick the program to handle the file as one that has been downloaded.
            self.source_url = source_files[0]
            tdir = tempfile.mkdtemp()
            source_files = self.getSourceFiles(tdir)
        #
        self.logger.info('Importing database {} from source files {}'.format(self.name, source_files))
        # create database and import file
        db = self.proj.db.newConnection()
        # remove database if already exist
        db.removeDatabase(self.name)
        # create a new one
        db.connect(self.name)
        cur = db.cursor()
        
        #
        # creating the field table
        self.logger.debug('Creating {}_field table'.format(self.name))
        self.createFieldsTable(db)
        for field in self.fields:
            cur.execute('INSERT INTO {0}_field VALUES ({1},{1},{1},{1},{1});'.format(self.name, self.proj.db.PH),
                (field.name, field.index, field.type, field.null, field.comment))
        db.commit()
        #
        # creating the info table
        self.logger.debug('Creating {}_info table'.format(self.name))
        query = 'INSERT INTO {0}_info VALUES ({1},{1});'.format(self.name, self.proj.db.PH)
        self.createInfoTable(db)
        cur.execute(query, ('name', self.name))
        cur.execute(query, ('anno_type', self.anno_type))
        cur.execute(query, ('description', self.description))
        cur.execute(query, ('version', self.version))
        cur.execute(query, ('build', str(self.build)))
        db.commit()
        self.logger.debug('Creating table {}'.format(self.name))
        self.createAnnotationTable(db)
        #
        # read records from files
        if self.source_type == 'txt':
            self.importTxtRecords(db, source_files)
        else:
            raise ValueError('Unrecognizable source input type: {}'.format(self.source_type))
        #
        # creating indexes
        s = delayedAction(self.logger.info, 'Creating indexes (this can take quite a while)')
        #
        # creates index for each link method
        for key in self.build.keys():
            if key != '*':
                cur.execute('''CREATE INDEX {0}_idx ON {1} ({0}_bin ASC, {2});'''\
                    .format(key, self.name,  ', '.join(['{} ASC'.format(x) for x in self.build[key]])))
            else:
                cur.execute('''CREATE INDEX {0}_idx ON {0} ({1});'''\
                    .format(self.name,  ', '.join(['{} ASC'.format(x) for x in self.build[key]])))
        # This is only useful for sqlite
        db.analyze()
        del s
        if tdir is not None:
            shutil.rmtree(tdir) 

    def prepareDB(self, source_files=[], linked_by=[]):
        '''Importing data to database. If direct_url or source_url is specified,
        they will overwrite settings in configuraiton file. If successful, this
        function set self.db to a live connection.
        '''
        if not source_files:
            dbFile = self.name + ('-' + self.version if self.version else '') + '.DB'
            if os.path.isfile(dbFile):
                try:
                    return AnnoDB(self.proj, dbFile)
                except ValueError as e:
                    self.logger.debug(e)
                    self.logger.info('Existing database cannot be used.')
            # if there is a direct URL?
            if self.direct_url is not None:
                self.logger.info('Downloading annotation database from {}'.format(self.direct_url))
                try:
                    dbFile = downloadFile(self.direct_url, '.')
                    s = delayedAction(self.logger.info, 'Decompressing {}'.format(dbFile))
                    dbFile = decompressIfNeeded(dbFile, inplace=True)
                    del s
                    return AnnoDB(self.proj, dbFile, linked_by)
                except Exception as e:
                    self.logger.debug(e)
                    self.logger.info('Failed to download database or downloaded database unusable.')
        # have to build from source
        self.importFromSource(source_files)
        return AnnoDB(self.proj, self.name, linked_by)


def useArguments(parser):
    parser.add_argument('source',
        help='''Use an annotation database ($source.DB or $source.DB.gz) if it is available,
            download or build the database if a description file ($source.ann) is available.
            Otherwise, this command will download a description file and the corresponding database
            from web (http://vtools.houstonbioinformatics.org/annoDB/$source.ann and the latest
            version of the datavase). If all means fail, this command will try to download the
            source of the annotation database (or use source files provided by option --files).''')
    parser.add_argument('-f', '--files', nargs='*', default=[],
        help='''A list of source files. If specified, vtools will not try to
            download and select source files. This is used only when no local
            annotation database is located.''')
    parser.add_argument('-l', '--linked_by', nargs='*', default=[],
        help='''A list of fields that are used to link the annotation database to
            tables in the existing project. This parameter is required only for
            'attribute' type of annotation databases that link to fields of existing
            tables.''')


def use(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # try to get source.ann, source.DB, source.DB.gz or get source.ann from 
            # http://vtools.houstonbioinformatics.org/annoDB
            if os.path.isfile(args.source):
                # if a local file?
                s = delayedAction(proj.logger.info, 'Decompressing {}'.format(args.source))
                # do not remove local .gz file. Perhaps this is a script and we do not want to break that.
                annoDB = decompressIfNeeded(args.source, inplace=False)
                del s
            elif os.path.isfile(args.source + '.DB.gz'):
                s = delayedAction(proj.logger.info, 'Decompressing {}'.format(args.source + '.DB.gz'))
                annoDB = decompressIfNeeded(args.source + '.DB.gz', inplace=False)
                del s
            elif os.path.isfile(args.source + '.DB'):
                annoDB = args.source + '.DB'
            elif os.path.isfile(args.source + '.ann'):
                annoDB = args.source + '.ann'
            else:
                res = urlparse.urlsplit(args.source)
                if not res.scheme:
                    args.source = 'http://vtools.houstonbioinformatics.org/annoDB/{}.ann'.format(args.source)
                # download?
                if proj.db.engine == 'mysql':
                    raise RuntimeError('MySQL databases are not portable and cannot be downloaded.')
                #
                proj.logger.info('Downloading annotation database from {}'.format(args.source))
                try:
                    annoDB = downloadFile(args.source, '.')
                    s = delayedAction(proj.logger.info, 'Decompressing {}'.format(annoDB))
                    # for downloaded file, we decompress inplace
                    annoDB = decompressIfNeeded(annoDB, inplace=True)
                    del s
                except Exception as e:
                    proj.logger.debug(e)
                    proj.logger.info('Failed to download database.')
            #
            # annDB is now a local file
            if annoDB.endswith('.ann'):
                if os.path.isfile(annoDB):
                    cfg = AnnoDBConfiger(proj, annoDB)
                    return proj.useAnnoDB(cfg.prepareDB(args.files, args.linked_by))
                else:
                    raise ValueError('Failed to locate configuration file {}'.format(annoDB))
            elif annoDB.endswith('.DB'):
                if proj.db.engine != 'sqlite3':
                    raise ValueError('A sqlite3 annotation database cannot be used with a mysql project.')
                if os.path.isfile(annoDB):
                    return proj.useAnnoDB(AnnoDB(proj, annoDB, args.linked_by))
                else:
                    raise ValueError('Failed to locate annotation database {}'.format(annoDB))
            else: # missing file extension?
                # no extension? try mysql database, .ann and .DB
                if proj.db.engine == 'mysql' and proj.db.hasDatabase(annoDB):
                    return proj.useAnnoDB(AnnoDB(proj, annoDB, args.linked_by))
                if os.path.isfile(annoDB + '.DB'):
                    if proj.db.engine != 'sqlite3':
                        raise ValueError('A sqlite3 annotation database cannot be used with a mysql project.')
                    try:
                        return proj.useAnnoDB(AnnoDB(proj, annoDB + '.DB', args.linked_by))
                    except Exception as e:
                        proj.logger.debug(e)
                if os.path.isfile(annoDB + '.ann'):
                    cfg = AnnoDBConfiger(proj, annoDB + '.ann')
                    try:
                        return proj.useAnnoDB(cfg.prepareDB(args.files, args.linked_by))
                    except Exception as e:
                        proj.logger.debug(e)
                # do not know what to do
                #FIXME if line 470 fails due to "No reference genome information" then the error prompt below is misleading(wanggao) 
                raise ValueError('Cannot find annotation database {}'.format(annoDB))
    except Exception as e:
        sys.exit(e)



