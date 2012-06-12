#!/usr/bin/env python
#
# $File: annotation.py $
# $LastChangedDate$
# $Rev$
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
import os
import ConfigParser
import shutil
import urlparse
import gzip
import zipfile
from multiprocessing import Process, Pipe

from .project import AnnoDB, Project, Field, AnnoDBWriter
from .utils import ProgressBar, downloadFile, lineCount, \
    DatabaseEngine, getMaxUcscBin, delayedAction, decompressIfNeeded, \
    normalizeVariant, compressFile, SQL_KEYWORDS, extractField
from .importer import LineProcessor, TextReader
  
class AnnoDBConfiger:
    '''An annotation database can be created from either a configuration file
    or a database name. In the former case, the configuration file will be
    parsed to determine the reference genome and fields the database provides.
    The database will not be download or imported unless one of the fields
    is used. In the later case, the information will be directly read from
    the annotation table in the database. Please check annotation/dbNSFP.ann
    for a description of the format of the configuration file.'''
    def __init__(self, proj, annoDB, jobs=2):
        '''Create an annotation database.
        proj: current project
        annoDB: annotation database. It can be either a mysql or sqlite database
            name depend on the database engine used in the project, or a
            configuration file. In the latter case, a filename is allowed.
        '''
        self.proj = proj
        self.logger = proj.logger
        self.jobs = jobs
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
        self.delimiter = '\t'
        self.version = None
        self.encoding = 'utf-8'
        # where is annoDB, in which form?
        self.parseConfigFile(annoDB)
        # some fields have to be determined.
        if self.name is None:
            raise ValueError('No valid name is set')
        if self.anno_type not in ['position', 'variant', 'range', 'field']:
            raise ValueError('vtools only support position, variant, range, and field types')
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
            self.version = self.name.split('-', 1)[1]
            self.name = self.name.split('-')[0]
        parser = ConfigParser.SafeConfigParser()
        # this allows python3 to read .fmt file with non-ascii characters, but there is no
        # simple way to make it python2 compatible.
        #with open(filename, 'r', encoding='UTF-8') as inputfile:
        #    parser.readfp(inputfile)
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
                    # field can be expressions such as pos + 100
                    if extractField(field) not in sections:
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
            elif item[0] == 'encoding':
                self.encoding = item[1]
            elif item[0] == 'source_url':
                self.source_url = item[1]
            elif item[0] == 'source_pattern':
                self.source_pattern = item[1]
            elif item[0] == 'source_type':
                self.source_type = item[1]
            elif item[0] == 'description':
                self.description = item[1]
            elif item[0] == 'anno_type':
                self.anno_type = item[1] if item[1] != 'attribute' else 'field'
            elif item[0] == 'delimiter':
                self.delimiter = eval(item[1])
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
            if section.lower() == 'linked fields' or section.lower() == 'data sources':
                continue
            if not section.replace('_', '').isalnum():
                raise ValueError('Illegal field name {}. Field names can only contain alphanumeric characters and underscores'.format(repr(section)))
            if section.upper() in SQL_KEYWORDS:
                raise ValueError('Illegal field name. {} conflicts with SQL keywords'.format(repr(section)))
            try:
                items = [x[0] for x in parser.items(section, raw=True)]
                for item in items:
                    if item not in ['index', 'type', 'adj', 'comment']:
                        raise ValueError('Incorrect key {} in section {}. Only index, type, adj and comment are allowed.'.format(item, section))
                # 'index' - 1 because the input is in 1-based index
                self.fields.append(Field(name=section, index=parser.get(section, 'index', raw=True),
                    type=parser.get(section, 'type', raw=True),
                    adj=parser.get(section, 'adj', raw=True) if 'adj' in items else None,
                    comment=parser.get(section, 'comment', raw=True) if 'comment' in items else ''))
            except Exception as e:
                self.logger.debug(e)
                self.logger.error('Invalid section {} in configuration file {}'.format(section, filename))
                return False
        # Do not sort fields, use the order sepcified in the .ann file.
        #self.fields.sort(key=lambda k: k.index.split()[0])
        for field in self.fields:
            self.logger.debug("Field {0}: index {1},\t{2}".format(field.name, field.index, field.type))


    def getSourceFiles(self):
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
            elif self.source_url.upper().startswith('SQL:'):
                res = urlparse.urlparse(self.source_url)
                user = res.user
                password = res.password
                query = res.netloc
                db = DababaseEngine(engine='mysql', user=user, passwd=password)
                cur = db.cursor()
                res = db.execute(query)
                filename = os.path.join(self.proj.temp_dir, '{}_sql.txt'.format(self.name))
                with open(filename, 'w') as output:
                    for rec in res:
                        output.write('{}\n'.format(','.join([str(x) for x in rec])))
                return [filename]
            else:
                filename = os.path.split(self.source_url)[-1]
                self.logger.info('Downloading {}'.format(filename))
                tempFile = downloadFile(self.source_url)
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
        bundle.extractall(self.proj.temp_dir)
        return [os.path.join(self.proj.temp_dir, name) for name in bundle.namelist()]
    
    def importTxtRecords(self, db, source_files):
        #
        build_info = []
        for key,items in self.build.items():
            if key == '*':
                continue
            try:
                # items have chr/pos, chr/pos/ref/alt, chr/start/end for different annotation types
                pos_idx = [i for i,x in enumerate(self.fields) if x.name == extractField(items[1])][0]
                if self.anno_type == 'variant':
                    # save indexes for pos, ref and alt
                    ref_idx = [i for i,x in enumerate(self.fields) if x.name == extractField(items[2])][0]
                    alt_idx = [i for i,x in enumerate(self.fields) if x.name == extractField(items[3])][0]
                    build_info.append((pos_idx, ref_idx, alt_idx))
                else:
                    build_info.append((pos_idx, ))
            except Exception as e:
                self.logger.error('No field {} for build {}: {}'.format(items[1], key, e))
        #
        processor = LineProcessor(self.fields, build_info, self.delimiter, None, self.logger)
        # files?
        cur = db.cursor()
        insert_query = 'INSERT INTO {0} VALUES ('.format(self.name) + \
                            ','.join([db.PH] * (len(self.fields) + len(build_info))) + ');'
        for f in source_files:
            self.logger.info('Importing annotation data from {0}'.format(f))
            skipped_lines = 0
            lc = lineCount(f, self.encoding)
            update_after = min(max(lc//200, 100), 100000)
            p = TextReader(processor, f, None, self.jobs - 1, self.encoding, self.logger)
            prog = ProgressBar(os.path.split(f)[-1], lc)
            all_records = 0
            skipped_records = 0
            for all_records, bins, rec in p.records():
                try:
                    cur.execute(insert_query, bins + rec)
                except Exception as e:
                    skipped_records += 1
                    self.logger.debug('Failed to process record {}: {}'.format(rec, e))
                if all_records % update_after == 0:
                    prog.update(all_records)
                    db.commit()
            db.commit()
            prog.done()
            self.logger.info('{0} records are handled, {1} ignored.'\
                .format(all_records, processor.skipped_lines + skipped_records))

    def importFromSource(self, source_files):
        '''Importing data from source files'''
        if source_files == []:
            source_files = self.getSourceFiles()
        elif len(source_files) == 1 and source_files[0].endswith('.zip'):
            # trick the program to handle the file as one that has been downloaded.
            self.source_url = source_files[0]
            source_files = self.getSourceFiles()
        if self.source_pattern is not None:
            source_files = [x for x in source_files if self.source_pattern in x]
        #
        self.logger.info('Importing database {} from source files {}'.format(self.name, source_files))
        #
        writer = AnnoDBWriter(self.name, self.fields, self.anno_type, self.description,
            self.version, self.build, self.logger)
        # read records from files
        if self.source_type == 'txt':
            self.importTxtRecords(writer.db, source_files)
        else:
            raise ValueError('Unrecognizable source input type: {}'.format(self.source_type))
        #
        writer.finalize()

    def prepareDB(self, source_files=[], linked_by=[], rebuild=False, anno_type=None, linked_fields=None):
        '''Importing data to database. If direct_url or source_url is specified,
        they will overwrite settings in configuraiton file. If successful, this
        function set self.db to a live connection.
        '''
        if not source_files and not rebuild:
            dbFile = self.name + ('-' + self.version if self.version else '') + '.DB'
            if os.path.isfile(dbFile):
                try:
                    return AnnoDB(self.proj, dbFile, linked_by, anno_type, linked_fields)
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
                    return AnnoDB(self.proj, dbFile, linked_by, anno_type, linked_fields)
                except Exception as e:
                    self.logger.debug(e)
                    self.logger.info('Failed to download database or downloaded database unusable.')
        # have to build from source
        self.importFromSource(source_files)
        #
        if rebuild:
            dbFile = self.name + ('-' + self.version if self.version else '') + '.DB.gz'
            self.logger.info('Creating compressed database {}'.format(dbFile))
            compressFile(self.name + '.DB', dbFile)
        return AnnoDB(self.proj, self.name, linked_by, anno_type, linked_fields)


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
            download and select source files. These source files will be
            compiled into a local annotation database. This is used only
            when no local annotation database is specified.''')
    parser.add_argument('-l', '--linked_by', nargs='*', default=[],
        help='''A list of fields that are used to link the annotation database to
            tables in the existing project. This parameter is required only for
            'field' type of annotation databases that link to fields of existing
            tables.''')
    parser.add_argument('--anno_type', choices=['variant', 'position', 'range', 'field'],
        help='''This option overrides type of an existing annotation database when it
            is attached to a project. It corresponds to key anno_type of the data sources
            section of an annotation file (with suffix .ann) but does not affect the .ann file
            or the database built from it.''')
    parser.add_argument('--linked_fields', nargs='*',
        help='''An alternative set of fields that are used to link the annotation database to
            the master variant table. It should have four, two, and three values for database
            of type variant, position, and range. Similar to anno_type, this option does not
            affect the .ann file or the database built from it.''')
    parser.add_argument('--rebuild', action='store_true',
        help='''If set, variant tools will always rebuild the annotation database from source,
            ignoring existing local and online database. In addition to $name.DB, variant tools
            will also create $name-$version.DB.gz that can be readily distributed.'''),
    parser.add_argument('-j', '--jobs', metavar='N', type=int, default=2,
        help='''If need to build database from source, maximum number of processes to use.''')


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
                    annoDB = downloadFile(args.source, None if args.source.endswith('.ann') else '.', quiet=args.source.endswith('.ann'))
                    s = delayedAction(proj.logger.info, 'Decompressing {}'.format(annoDB))
                    # for downloaded file, we decompress inplace
                    annoDB = decompressIfNeeded(annoDB, inplace=True)
                    del s
                except Exception as e:
                    proj.logger.debug(e)
                    raise ValueError('Failed to download database.')
            #
            # annDB is now a local file
            if annoDB.endswith('.ann'):
                if os.path.isfile(annoDB):
                    cfg = AnnoDBConfiger(proj, annoDB, args.jobs)
                    return proj.useAnnoDB(cfg.prepareDB(args.files, args.linked_by, args.rebuild,
                        args.anno_type, args.linked_fields))
                else:
                    raise ValueError('Failed to locate configuration file {}'.format(annoDB))
            elif args.rebuild:
                raise ValueError('Only an .ann file can be specified when option --rebuild is set')
            elif annoDB.endswith('.DB'):
                if proj.db.engine != 'sqlite3':
                    raise ValueError('A sqlite3 annotation database cannot be used with a mysql project.')
                if os.path.isfile(annoDB):
                    return proj.useAnnoDB(AnnoDB(proj, annoDB, args.linked_by, args.anno_type, args.linked_fields))
                else:
                    raise ValueError('Failed to locate annotation database {}'.format(annoDB))
            else: # missing file extension?
                # no extension? try mysql database, .ann and .DB
                if proj.db.engine == 'mysql' and proj.db.hasDatabase(annoDB):
                    return proj.useAnnoDB(AnnoDB(proj, annoDB, args.linked_by, args.anno_type, args.linked_fields))
                if os.path.isfile(annoDB + '.DB'):
                    if proj.db.engine != 'sqlite3':
                        raise ValueError('A sqlite3 annotation database cannot be used with a mysql project.')
                    try:
                        return proj.useAnnoDB(AnnoDB(proj, annoDB + '.DB', args.linked_by, args.anno_type, args.linked_fields))
                    except Exception as e:
                        proj.logger.debug(e)
                if os.path.isfile(annoDB + '.ann'):
                    cfg = AnnoDBConfiger(proj, annoDB + '.ann', args.jobs)
                    try:
                        return proj.useAnnoDB(cfg.prepareDB(args.files, args.linked_by, args.rebuild, args.anno_type, args.linked_fields))
                    except Exception as e:
                        proj.logger.debug(e)
                # do not know what to do
                #FIXME if line 470 fails due to "No reference genome information" then the error prompt below is misleading(wanggao) 
                raise ValueError('Cannot find annotation database {}'.format(annoDB))
    except Exception as e:
        sys.exit(e)



