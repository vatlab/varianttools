#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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

import configparser
import os
import re
import sys
import tarfile
import urllib.parse
import zipfile

from .importer import LineProcessor
from .project import AnnoDB, AnnoDBWriter, Field, Project
from .text_reader import TextReader
from .utils import (SQL_KEYWORDS, DatabaseEngine, ProgressBar, RefGenome,
                    ResourceManager, calculateMD5, compressFile,
                    decompressGzFile, delayedAction, downloadFile, env,
                    extractField, isAnnoDB, lineCount)


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
        annoDB: annotation database. It can be either a sqlite database
            name depend on the database engine used in the project, or a
            configuration file. In the latter case, a filename is allowed.
        '''
        self.proj = proj
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
        self.db_md5 = None
        self.source_type = None
        self.source_url = None
        self.source_files = None
        self.fields = []
        self.delimiter = '\t'
        # the original version is 1
        # Databases created using normalized variant importer engine will be assigned
        # format 2 so that older version will be automatically upgraded.
        self.database_format = '2'
        self.version = None
        self.encoding = 'utf-8'
        self.preprocessor = None
        self.header = None
        # where is annoDB, in which form?
        self.parseConfigFile(annoDB)
        # some fields have to be determined.
        if self.name is None:
            raise ValueError('No valid name is set')
        if self.anno_type not in ['position', 'variant', 'range', 'field']:
            raise ValueError(
                'vtools only support position, variant, range, and field types')
        if len(self.fields) == 0:
            raise ValueError(
                'No valid field is located from database {}'.format(annoDB))
        if len(self.build) == 0:
            raise ValueError(
                'No reference genome is specified for database {}'.format(
                    annoDB))
        if proj.build is not None:
            if (list(self.build.keys())[0] != '*') and (proj.build not in list(self.build.keys())) \
                and (proj.alt_build is None or
                                                                                                            proj.alt_build not in list(self.build.keys())):
                raise ValueError(
                    'Annotation database cannot be used because it is based on a reference '
                    'genome that is different from the one used by the project. Please use a version of '
                    'annotation databse for the project (vtools show annotations), or liftover the existing '
                    'project (vtoos liftover) to make it compatible with the annotation database.'
                )

    def remove(self):
        '''Remove an annotation database'''
        env.logger.info('Removing database {}'.format(self.name))
        self.proj.db.removeDatabase(self.name)

    def parseConfigFile(self, filename):
        '''Read from an ini style configuration file'''
        env.logger.trace('Checking configuration file {}'.format(filename))
        self.path = os.path.split(filename)[0]
        # get filename, remove extension, and keep version as things after -.
        self.name = os.path.splitext(os.path.split(filename)[-1])[0]
        self.version = None
        if '-' in self.name:
            self.version = self.name.split('-', 1)[1]
            self.name = self.name.split('-')[0]
        parser = configparser.SafeConfigParser()
        # this allows python3 to read .fmt file with non-ascii characters, but there is no
        # simple way to make it python2 compatible.
        # with open(filename, 'r', encoding='UTF-8') as inputfile:
        #    parser.readfp(inputfile)
        parser.read(filename)
        # sections?
        sections = parser.sections()
        if 'linked fields' not in sections:
            raise ValueError(
                'Ignore invalid annotation file {}. Missing build section.'
                .format(filename))
        if 'data sources' not in sections:
            raise ValueError(
                'Ignore invalid annotation file {}. Missing source section.'
                .format(filename))
        # linked fields for each reference genome
        try:
            self.build = {}
            for item in parser.items('linked fields'):
                fields = [x.strip() for x in item[1].split(',')]
                for field in fields:
                    # field can be expressions such as pos + 100
                    if extractField(field) not in sections:
                        env.logger.error(
                            'Invalid reference genome: Unspecified field {}.'
                            .format(field))
                self.build[item[0]] = fields
            env.logger.trace('Reference genomes: {}'.format(self.build))
        except Exception as e:
            env.logger.debug(e)
            raise ValueError(
                'Failed to obtain build information of annotation file {}'
                .format(filename))
        # source information
        self.anno_type = 'variant'
        self.description = None
        self.direct_url = None
        self.source_url = None
        self.source_pattern = None
        self.source_type = None
        for item in parser.items('data sources'):
            if item[0] == 'direct_url':
                # for backward compatibility, change old houstonbioinformatics.org host
                # to new one
                if item[1].startswith(
                        'http://vtools.houstonbioinformatics.org/'):
                    self.direct_url = item[1][len(
                        'http://vtools.houstonbioinformatics.org/'):].strip()
                else:
                    self.direct_url = item[1].strip()
                #
                # if a md5 signature is put after the URL, in the format of
                # fileToGet md5=xxxxxxxxxx, assign self.db_md5
                if '\t' in self.direct_url:
                    url, md5 = self.direct_url.rsplit('\t', 1)
                    self.direct_url = url
                    self.db_md5 = md5
            elif item[0] == 'encoding':
                self.encoding = item[1]
            elif item[0] == 'preprocessor':
                self.preprocessor = item[1]
            elif item[0] == 'header':
                if item[1] in ('none', 'None'):
                    self.header = None
                else:
                    try:
                        self.header = int(item[1])
                    except:
                        # in this case header is a pattern
                        self.header = re.compile(item[1])
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
                    raise Warning(
                        'Version obtained from filename ({}) is inconsistent with version specified in the annotation file ({})'
                        .format(self.version, item[1]))
                self.version = item[1]
            else:
                raise ValueError(
                    'Invalid keyword {} in section "data sources".'.format(
                        item[0]) +
                    'Only direct_url, source_url, source_type, version, preprocessor, source_pattern and description are allowed'
                )
        # sections
        self.fields = []
        for section in sections:
            if section.lower() == 'linked fields' or section.lower(
            ) == 'data sources':
                continue
            if not section.replace('_', '').isalnum():
                raise ValueError(
                    'Illegal field name {}. Field names can only contain alphanumeric characters and underscores'
                    .format(repr(section)))
            if section.upper() in SQL_KEYWORDS:
                raise ValueError(
                    'Illegal field name. {} conflicts with SQL keywords'.format(
                        repr(section)))
            try:
                items = [x[0] for x in parser.items(section, raw=True)]
                for item in items:
                    if item not in ['index', 'type', 'adj', 'comment']:
                        raise ValueError(
                            'Incorrect key {} in section {}. Only index, type, adj and comment are allowed.'
                            .format(item, section))
                # 'index' - 1 because the input is in 1-based index
                self.fields.append(
                    Field(
                        name=section,
                        index=parser.get(section, 'index', raw=True),
                        type=parser.get(section, 'type', raw=True),
                        adj=parser.get(section, 'adj', raw=True)
                        if 'adj' in items else None,
                        fmt=None,
                        comment=parser.get(section, 'comment', raw=True)
                        if 'comment' in items else ''))
            except Exception as e:
                env.logger.debug(e)
                env.logger.error(
                    'Invalid section {} in configuration file {}'.format(
                        section, filename))
                return False
        # Do not sort fields, use the order sepcified in the .ann file.
        #self.fields.sort(key=lambda k: k.index.split()[0])
        # for field in self.fields:
        #    env.logger.trace("Field {0}: index {1},\t{2}".format(field.name, field.index, field.type))

    def getSourceFiles(self):
        '''Download a file from a URL and save to a directory. If it is a zip file,
        unzip it in the directory.'''
        # FIXME: if we detect a file in the local directory, we should not try to
        # download from online? However, users might want to get the latest
        # version from online, and if they want to use a local copy, they can always
        # use option --files.
        if not self.source_url:
            raise ValueError(
                'No source_url is specified in .ann file. Please modify your spec file or use option -f to specify source files.'
            )
        URLs = self.source_url.strip().split()
        ret_files = []
        for URL in URLs:
            tempFile = None
            try:
                if os.path.isfile(URL):
                    if not URL.lower().endswith('.zip'):
                        ret_files.append(URL)
                    else:
                        tempFile = URL
                elif os.path.isfile(os.path.join(self.path, URL)):
                    if not URL.lower().endswith('.zip'):
                        ret_files.append(os.path.join(self.path, URL))
                    else:
                        tempFile = os.path.join(self.path, URL)
                elif URL.upper().startswith('SQL:'):
                    res = urllib.parse.urlparse(URL)
                    query = res.netloc
                    db = DatabaseEngine()
                    res = db.execute(query)
                    filename = os.path.join(env.cache_dir,
                                            '{}_sql.txt'.format(self.name))
                    with open(filename, 'w') as output:
                        for rec in res:
                            output.write('{}\n'.format(','.join(
                                [str(x) for x in rec])))
                    ret_files.append(filename)
                else:
                    filename = os.path.split(URL)[-1]
                    if not filename.strip():
                        raise ValueError('No source_url is specified.')
                    env.logger.info('Downloading {}'.format(filename))
                    tempFile = downloadFile(URL, dest_dir=env.cache_dir)
            except Exception as e:
                raise ValueError(
                    'Failed to download database source from {}: {}'.format(
                        URL, e))
            #
            if tempFile is None:
                continue
            if not os.path.isfile(tempFile):
                raise ValueError('Could not find downloaded database source')
            # regular file?
            if tempFile.endswith('.zip'):
                # if zip file?
                bundle = zipfile.ZipFile(tempFile)
                bundle.extractall(env.cache_dir)
                ret_files.extend([
                    os.path.join(env.cache_dir, name)
                    for name in bundle.namelist()
                ])
            elif tempFile.endswith('.tar.gz') or tempFile.endswith('.tgz'):
                with tarfile.open(tempFile, 'r:gz') as tfile:
                    with delayedAction(env.logger.info,
                                       'Decompressing {}'.format(tempFile)):
                        names = tfile.getnames()
                        tfile.extractall(env.cache_dir)
                ret_files.extend(
                    [os.path.join(env.cache_dir, name) for name in names])
            else:
                ret_files.append(tempFile)
        return ret_files

    def importTxtRecords(self, db, source_files):
        #
        build_info = []
        for key, items in list(self.build.items()):
            if key == '*':
                continue
            try:
                # items have chr/pos, chr/pos/ref/alt, chr/start/end for
                # different annotation types
                chr_idx = [
                    i for i, x in enumerate(self.fields)
                    if x.name == extractField(items[0])
                ][0]
                pos_idx = [
                    i for i, x in enumerate(self.fields)
                    if x.name == extractField(items[1])
                ][0]
                if self.anno_type == 'variant':
                    # save indexes for pos, ref and alt
                    ref_idx = [
                        i for i, x in enumerate(self.fields)
                        if x.name == extractField(items[2])
                    ][0]
                    alt_idx = [
                        i for i, x in enumerate(self.fields)
                        if x.name == extractField(items[3])
                    ][0]
                    build_info.append(
                        (RefGenome(key), chr_idx, pos_idx, ref_idx, alt_idx))
                else:
                    build_info.append((pos_idx,))
            except Exception as e:
                env.logger.error('No field {} for build {}: {}'.format(
                    items[1], key, e))
        #
        processor = LineProcessor(self.fields, build_info, self.delimiter, None)
        # files?
        cur = db.cursor()
        insert_query = 'INSERT INTO {0} VALUES ('.format(self.name) + \
            ','.join([db.PH] * (len(self.fields) + len(build_info))) + ');'
        for idx, f in enumerate(source_files):
            env.logger.info('Importing annotation data from {0}'.format(f))
            lc = lineCount(f, self.encoding)
            update_after = min(max(lc // 200, 100), 100000)
            p = TextReader(processor, f, None, None, self.jobs - 1,
                           self.encoding, self.header)
            prog = ProgressBar(
                '{} ({}/{})'.format(
                    os.path.split(f)[-1], idx + 1, len(source_files)), lc)
            all_records = 0
            skipped_records = 0
            for all_records, bins, rec in p.records():
                try:
                    cur.execute(insert_query, bins + rec)
                except Exception as e:
                    skipped_records += 1
                    env.logger.debug('Failed to process record {}: {}'.format(
                        rec, e))
                if all_records % update_after == 0:
                    prog.update(all_records)
                    db.commit()
            db.commit()
            prog.done()
            env.logger.info('{0} records are handled, {1} ignored.'.format(
                all_records, p.unprocessable_lines + skipped_records))

    def importFromSource(self, source_files):
        '''Importing data from source files'''
        if source_files == []:
            source_files = self.getSourceFiles()
        elif len(source_files) == 1 and source_files[0].endswith('.zip'):
            # trick the program to handle the file as one that has been
            # downloaded.
            self.source_url = source_files[0]
            source_files = self.getSourceFiles()
        if self.source_pattern is not None:
            source_files = [x for x in source_files if self.source_pattern in x]
        #
        env.logger.info('Importing database {} from {} source files {}'.format(
            self.name, len(source_files), ', '.join(source_files)))
        #
        if self.preprocessor is not None:
            env.logger.info(
                'Preprocessing data [{}] to generate intermediate input files for import'
                .format(', '.join(source_files)))
            # if this is the case, only one input stream will be allowed.
            # process command line
            command = self.preprocessor
            # replace command with other stuff, if applicable
            command = command.replace('$build', "'{}'".format(self.build))
            #
            # create temp files
            temp_files = [
                os.path.join(env.cache_dir,
                             os.path.basename(x) + '.' + self.name)
                for x in source_files
            ]
            try:
                processor = eval(command)
                # intermediate files will be named as
                # "cache_dir/$inputfilename.$(self.name)"
                processor.convert(source_files, temp_files)
                for output in temp_files:
                    if not os.path.isfile(output):
                        raise ValueError(
                            "Preprocessed file {} does not exist.".format(
                                output))
            except Exception as e:
                raise ValueError(
                    "Failed to execute preprocessor '{}': {}".format(
                        re.sub(r'\((.*)\)', '', command), e))
            #
            # we record file as cache files
            source_files = temp_files
        writer = AnnoDBWriter(self.name, self.fields, self.anno_type,
                              self.description, self.version, self.build,
                              self.database_format)
        # read records from files
        if self.source_type == 'txt':
            self.importTxtRecords(writer.db, source_files)
        else:
            raise ValueError('Unrecognizable source input type: {}'.format(
                self.source_type))
        #
        writer.finalize()

    def prepareDB(self,
                  source_files=[],
                  linked_by=[],
                  rebuild=False,
                  anno_type=None,
                  linked_fields=None,
                  linked_name=None):
        '''Importing data to database. If direct_url or source_url is specified,
        they will overwrite settings in configuraiton file. If successful, this
        function set self.db to a live connection.
        '''
        if not source_files and not rebuild:
            dbFilename = self.name + \
                ('-' + self.version if self.version else '') + '.DB'
            for dbFile in (dbFilename,
                           os.path.join(env.local_resource, 'annoDB',
                                        dbFilename)):
                if os.path.isfile(dbFile):
                    if self.db_md5 and calculateMD5(
                            dbFile, partial=True) != self.db_md5:
                        env.logger.info(
                            '{}: MD5 signature mismatch, the database might have been upgraded locally.'
                            .format(dbFile))
                    try:
                        return AnnoDB(self.proj, dbFile, linked_by, anno_type,
                                      linked_fields, linked_name)
                    except ValueError as e:
                        env.logger.debug(e)
                        env.logger.info('Existing database cannot be used.')
            # if there is a direct URL?
            if self.direct_url is not None:
                env.logger.info(
                    'Downloading annotation database from {}'.format(
                        self.direct_url))
                try:
                    dbFile = downloadFile(self.direct_url)
                    with delayedAction(env.logger.info,
                                       'Decompressing {}'.format(dbFile)):
                        dbFile = decompressGzFile(
                            dbFile, inplace=False, md5=self.db_md5)
                    return AnnoDB(self.proj, dbFile, linked_by, anno_type,
                                  linked_fields, linked_name)
                except RuntimeError:
                    raise
                except Exception as e:
                    env.logger.warning(
                        'Failed to download database or downloaded database unusable: {}'
                        .format(e))
        # have to build from source
        self.importFromSource(source_files)
        #
        if rebuild:
            dbFile = self.name + \
                ('-' + self.version if self.version else '') + '.DB.gz'
            env.logger.info(
                'Creating compressed database {} for DB with md5 {}'.format(
                    dbFile, calculateMD5(self.name + '.DB', partial=True)))
            compressFile(self.name + '.DB', dbFile)
        return AnnoDB(self.proj, self.name, linked_by, anno_type, linked_fields,
                      linked_name)


def useArguments(parser):
    parser.add_argument(
        'source',
        help='''Use an annotation database ($source.DB or $source.DB.gz) if it is available,
            download or build the database if a description file ($source.ann) is available.
            Otherwise, this command will download a description file and the corresponding
            database from web (c.f. runtime variable $search_path) and the latest version
            of the datavase). If all means fail, this command will try to download the
            source of the annotation database (or use source files provided by option --files).'''
    )
    grp = parser.add_argument_group('Basic link options')
    grp.add_argument(
        '--as',
        metavar='NAME',
        help='''An alternative name for the linked database. This option allows
            the use of shorter field names (e.g. tg.chr instead of thousandGenomes.chr)
            and the use of multiple versions of the same database.''')
    grp.add_argument(
        '-l',
        '--linked_by',
        '--linked-by',
        nargs='*',
        default=[],
        metavar='FIELD',
        help='''A list of fields that are used to link the annotation database to
            tables in the existing project. This parameter is required only for
            'field' type of annotation databases that link to fields of existing
            tables.''')
    grp = parser.add_argument_group('Advanced link options')
    grp.add_argument(
        '--anno_type',
        '--anno-type',
        choices=['variant', 'position', 'range', 'field'],
        help='''This option overrides type of an existing annotation database when it
            is attached to a project. It corresponds to key anno_type of the data sources
            section of an annotation file (with suffix .ann) but does not affect the .ann file
            or the database built from it.''')
    grp.add_argument(
        '--linked_fields',
        '--linked-fields',
        nargs='*',
        help='''An alternative set of fields that are used to link the annotation database to
            the master variant table. It should have four, two, and three values for database
            of type variant, position, and range. Similar to anno_type, this option does not
            affect the .ann file or the database built from it.''')
    grp = parser.add_argument_group('Build database from source')
    grp.add_argument(
        '-f',
        '--files',
        nargs='*',
        default=[],
        help='''A list of source files. If specified, vtools will not try to
            download and select source files. These source files will be
            compiled into a local annotation database. This is used only
            when no local annotation database is specified.''')
    grp.add_argument(
        '--rebuild',
        action='store_true',
        help='''If set, variant tools will always rebuild the annotation database from source,
            ignoring existing local and online database. In addition to $name.DB, variant tools
            will also create $name-$version.DB.gz that can be readily distributed.'''
    ),
    grp.add_argument(
        '-j',
        '--jobs',
        metavar='N',
        type=int,
        default=2,
        help='''If need to build database from source, maximum number of processes to use.'''
    )


def use(args):
    try:
        # since as is a python keyword, we cannot use args.linked_name to access it and has to rename it
        # to linked_name
        setattr(args, 'linked_name', getattr(args, 'as'))
        with Project(verbosity=args.verbosity) as proj:
            if args.linked_name is not None and not args.linked_name.isalnum():
                raise ValueError(
                    'Invalid alias to annotation database: {}'.format(
                        args.linked_name))
            # try to get source.ann, source.DB, source.DB.gz or get source.ann
            # from
            annoDB = None
            if os.path.isfile(args.source):
                # if a local file?
                if args.source.endswith('.ann'):
                    annoDB = args.source
                elif args.source.endswith('.DB.gz'):
                    with delayedAction(env.logger.info,
                                       'Decompressing {}'.format(args.source)):
                        # do not remove local .gz file. Perhaps this is a script
                        # and we do not want to break that.
                        annoDB = decompressGzFile(args.source, inplace=False)
                    if not isAnnoDB(annoDB):
                        annoDB = None
                elif args.source.endswith('.DB'):
                    annoDB = args.source
                    if not isAnnoDB(annoDB):
                        annoDB = None
            if annoDB is None:
                if os.path.isfile(args.source + '.DB.gz'):
                    with delayedAction(
                            env.logger.info,
                            'Decompressing {}'.format(args.source + '.DB.gz')):
                        annoDB = decompressGzFile(
                            args.source + '.DB.gz', inplace=False)
                        if not isAnnoDB(annoDB):
                            annoDB = None
            if annoDB is None:
                if os.path.isfile(args.source +
                                  '.DB') and isAnnoDB(args.source + '.DB'):
                    annoDB = args.source + '.DB'
            if annoDB is None:
                if os.path.isfile(args.source + '.ann'):
                    annoDB = args.source + '.ann'
            if annoDB is None:
                res = urllib.parse.urlsplit(args.source)
                if not res.scheme:
                    # if this is a versioned annotation string, use the version
                    res = ResourceManager()
                    res.getLocalManifest()
                    res.selectFiles(resource_type='annotation')
                    avail_annoDBs = [
                        x for x in list(res.manifest.keys())
                        if x.endswith('.ann')
                    ]
                    if '-' in args.source:
                        if 'annoDB/' + args.source in avail_annoDBs:
                            args.source = 'annoDB/{}'.format(args.source)
                        elif 'annoDB/' + args.source + '.ann' in avail_annoDBs:
                            args.source = 'annoDB/{}.ann'.format(args.source)
                        else:
                            raise ValueError(
                                'Annotation database {} not found.'.format(
                                    args.source))
                    else:
                        # now, no version information is given, try to find the
                        # matching one.
                        all_versions = [
                            x for x in avail_annoDBs
                            if x.split('/')[1].split('-')[0] == args.source and
                            x.endswith('.ann')
                        ]
                        if not all_versions:
                            raise ValueError(
                                'Annotation database [] not found.'.format(
                                    args.source))
                        # find the matching reference genome
                        ref_filtered = [
                            x for x in all_versions
                            if res.manifest[x][2] == '*' or proj.build is None
                            or (proj.build is not None and
                                proj.build in res.manifest[x][2]) or
                            (proj.alt_build is not None and
                             proj.alt_build in res.manifest[x][2])
                        ]
                        #
                        if not ref_filtered:
                            raise ValueError(
                                'No database with matching reference genome is found: {}, available {}'
                                .format(args.source, ', '.join(all_versions)))
                        # now, we use the latest version
                        # search for date information ...
                        if all([
                                re.search('\d\d\d\d\d\d\d\d', x)
                                for x in ref_filtered
                        ]):
                            # date information information is used, use it.
                            args.source = sorted(
                                ref_filtered,
                                key=lambda x: re.search('\d\d\d\d\d\d\d\d', x).
                                group(0))[-1]
                        else:
                            args.source = sorted(ref_filtered)[-1]
                        env.logger.info(
                            'Choosing version {} from {} available databases.'
                            .format(
                                args.source.split('/')[-1].split('.')[0],
                                len(all_versions)))
                # download?
                env.logger.info('Downloading annotation database {}'.format(
                    args.source))
                try:
                    annoDB = downloadFile(
                        args.source,
                        None if args.source.endswith('.ann') else '.',
                        quiet=args.source.endswith('.ann'))
                    with delayedAction(env.logger.info,
                                       'Decompressing {}'.format(annoDB)):
                        # for downloaded file, we decompress inplace
                        annoDB = decompressGzFile(annoDB, inplace=True)
                except Exception as e:
                    raise ValueError('Database {} not found: {}'.format(
                        args.source, e))
            #
            # annDB is now a local file
            if annoDB.endswith('.ann'):
                if os.path.isfile(annoDB):
                    cfg = AnnoDBConfiger(proj, annoDB, args.jobs)
                    return proj.useAnnoDB(
                        cfg.prepareDB(args.files, args.linked_by, args.rebuild,
                                      args.anno_type, args.linked_fields,
                                      args.linked_name))
                else:
                    raise ValueError(
                        'Failed to locate configuration file {}'.format(annoDB))
            elif args.rebuild:
                raise ValueError(
                    'Only an .ann file can be specified when option --rebuild is set'
                )
            elif annoDB.endswith('.DB') and isAnnoDB(annoDB):
                if os.path.isfile(annoDB):
                    return proj.useAnnoDB(
                        AnnoDB(proj, annoDB, args.linked_by, args.anno_type,
                               args.linked_fields, args.linked_name))
                else:
                    raise ValueError(
                        'Failed to locate annotation database {}'.format(
                            annoDB))
            else:  # missing file extension?
                if os.path.isfile(annoDB + '.DB') and isAnnoDB(annoDB + '.DB'):
                    try:
                        return proj.useAnnoDB(
                            AnnoDB(proj, annoDB + '.DB', args.linked_by,
                                   args.anno_type, args.linked_fields,
                                   args.linked_name))
                    except Exception as e:
                        env.logger.debug(e)
                if os.path.isfile(annoDB + '.ann'):
                    cfg = AnnoDBConfiger(proj, annoDB + '.ann', args.jobs)
                    try:
                        return proj.useAnnoDB(
                            cfg.prepareDB(args.files, args.linked_by,
                                          args.rebuild, args.anno_type,
                                          args.linked_fields, args.linked_name))
                    except Exception as e:
                        env.logger.debug(e)
                # do not know what to do
                # FIXME if line 470 fails due to "No reference genome
                # information" then the error prompt below is
                # misleading(wanggao)
                raise ValueError(
                    'Cannot find annotation database {}'.format(annoDB))
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)
