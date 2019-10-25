#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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

import itertools as it
import os
import sys
from collections import defaultdict
from multiprocessing import Lock

from .plinkfile import PlinkFile
from .utils import (DatabaseEngine, ProgressBar, RefGenome, downloadFile, env,
                    openFile)


#
#
# Functors to process input
#
#
# Extractors to extract value from a field
class ExtractField:

    def __init__(self, index, sep=';', default=None):
        '''Define an extractor that returns the index-th (1-based) field of the fields
        separated by specified delimiter. Return default if unsuccessful.'''
        self.index = index - 1
        self.sep = sep
        self.default = default

    def __call__(self, item):
        try:
            return item.split(self.sep, self.index + 1)[self.index]
        except:
            return self.default


g_geneNameStandardizer = {}


class _GeneNameStandardizer:

    def __init__(self, convertTo='geneSymbol'):
        self.convertTo = convertTo
        # CREATE TABLE `kgAlias` (
        #  `kgID` varchar(40) NOT NULL default '',
        #  `alias` varchar(80) default NULL,
        #  KEY `kgID` (`kgID`),
        #  KEY `alias` (`alias`)
        #)
        kgAliasFile = downloadFile(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgAlias.txt.gz'
        )
        #CREATE TABLE `kgXref` (
        #  `kgID` varchar(255) NOT NULL,
        #  `mRNA` varchar(255) NOT NULL,
        #  `spID` varchar(255) NOT NULL,
        #  `spDisplayID` varchar(255) NOT NULL,
        #  `geneSymbol` varchar(255) NOT NULL,
        #  `refseq` varchar(255) NOT NULL,
        #  `protAcc` varchar(255) NOT NULL,
        #  `description` longblob NOT NULL,
        #  `rfamAcc` varchar(255) NOT NULL,
        #  `tRnaName` varchar(255) NOT NULL,
        #  KEY `kgID` (`kgID`),
        #  KEY `mRNA` (`mRNA`),
        #  KEY `spID` (`spID`),
        #  KEY `spDisplayID` (`spDisplayID`),
        #  KEY `geneSymbol` (`geneSymbol`),
        #  KEY `refseq` (`refseq`),
        #  KEY `protAcc` (`protAcc`),
        #  KEY `rfamAcc` (`rfamAcc`),
        #  KEY `tRnaName` (`tRnaName`)
        #) ENGINE=MyISAM DEFAULT CHARSET=latin1;
        #
        kgXRefFile = downloadFile(
            'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/kgXref.txt.gz'
        )
        self.nameMap = self.processAlias(kgAliasFile, kgXRefFile)

    def processAlias(self, kgAliasFile, kgXRefFile):
        aliasMap = {}
        XrefMap = {}
        with openFile(kgAliasFile) as kgAlias:
            for line in kgAlias:
                kgID, alias = line.decode('UTF8').strip().split('\t', 1)
                aliasMap[alias.upper()] = kgID
        try:
            resIndex = {
                'kgID': 0,
                'mRNA': 1,
                'spID': 2,
                'spDisplayID': 3,
                'geneSymbol': 4,
                'refseq': 5,
                'protAcc': 6,
                'description': 7,
                'rfamAcc': 8,
                'tRnaName': 9
            }[self.convertTo]
        except KeyError:
            raise ValueError(
                'Incorrect conversion type {}. Allowed types for GeneNameStandardizer '
                'are kgID, mRNA, spID, spDisplayID, geneSymbol, refseq, protAcc, description, rfamAcc, and tRnaName'
                .format(self.convertTo))
        with openFile(kgXRefFile) as kgXRef:
            for line in kgXRef:
                try:
                    fields = line.decode('UTF8').strip().split('\t', 9)
                except:
                    # ignore invalid lines
                    continue
                XrefMap[fields[0]] = fields[resIndex]
        # get the final map
        nameMap = {}
        for alias in aliasMap:
            # if not, just ignore the alias
            try:
                name = XrefMap[aliasMap[alias]]
            except:
                # if does not exist, just pass
                continue
            #
            if name == alias:
                continue
            # if there is alreay another name
            if alias in nameMap:
                #nameMap[alias].append(name)
                sys.stderr.write(
                    'Multiple {} names for alias {}: {} used\n'.format(
                        self.convertTo, alias, name))
            nameMap[alias] = name
        #
        return nameMap

    def __call__(self, item):
        try:
            #if item.upper() in self.nameMap:
            #    print item.upper(), '==> ', self.nameMap[item.upper()]
            return self.nameMap[item.upper()]
        except Exception:
            return item


def ConvertGeneName(convertTo='geneSymbol'):
    global g_geneNameStandardizer
    if convertTo not in g_geneNameStandardizer:
        g_geneNameStandardizer[convertTo] = _GeneNameStandardizer(convertTo)
    return g_geneNameStandardizer[convertTo]


class CheckSplit:

    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. Return default if unsuccessful. It differs from
        SplitField in that it will return the item itself (instead of a tuple
        of one element) when there is only one element. The item will then
        will copy if multiple items exist.'''
        self.sep = sep

    def __call__(self, item):
        return item if self.sep not in item else tuple(item.split(self.sep))


class SplitField:

    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep

    def __call__(self, item):
        return tuple(item.split(self.sep))


class ExtractFlag:

    def __init__(self, name, sep=';'):
        '''Define an extractor that returns 1 is item contains name as one of the fields,
        and 0 otherwise. No space is allowed between delimiter and flag.'''
        self.n = name
        self.s = name + sep
        self.e = sep + name
        self.m = sep + name + sep

    def __call__(self, item):
        # this should be faster than
        #
        #     if self.name in item.split(self.sep):
        #
        # because we do not have to split the whole string.
        #
        if self.n not in item:
            return '0'
        # put the most common case first
        if self.m in item or item.startswith(self.s) or item.endswith(
                self.e) or item == self.n:
            return '1'
        else:
            return '0'


class CommonLeading:
    '''Find the common leading piece of two input strings (ref and alt).
    '''

    def __init__(self):
        pass

    def _commonLeading(self, ref, alt):
        common_leading = 0
        for i in range(min(len(ref), len(alt))):
            if ref[i] == alt[i]:
                common_leading += 1
        return ref[:common_leading]

    def __call__(self, item):
        if ',' in item[1]:
            return tuple([
                self._commonLeading(item[0], alt) for alt in item[1].split(',')
            ])
        else:
            return self._commonLeading(item[0], item[1])


class CommonEnding:
    '''Find the common ending piece of two input strings (ref and alt).
    '''

    def __init__(self):
        pass

    def _commonEnding(self, ref, alt):
        common_leading = 0
        for i in range(min(len(ref), len(alt))):
            if ref[i] == alt[i]:
                common_leading += 1
        if common_leading > 0:
            ref = ref[common_leading:]
            alt = alt[common_leading:]
        common_ending = 0
        for i in range(-1, -min(len(ref), len(alt)) - 1, -1):
            if ref[i] == alt[i]:
                common_ending -= 1
            else:
                break
        if common_ending < 0:
            return ref[common_ending:]
        else:
            return ''

    def __call__(self, item):
        if ',' in item[1]:
            return tuple([
                self._commonEnding(item[0], alt) for alt in item[1].split(',')
            ])
        else:
            return self._commonEnding(item[0], item[1])


class __FieldFromFormat:

    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that return the value of a field according
        to a format string. This is used to extract stuff from the format
        string of vcf files.
        '''
        self.name = name
        self.sep = sep
        self.default = default
        self.factory = defaultdict(dict)
        self.index = {}

    def __call__(self, item):
        try:
            # first try to get from a global factory
            return self.factory[item[0]][item[1]]
        except:
            fmt, val = item
            try:
                # now split .... assuming the format has been handled before.
                # this should be the case most of the time
                res = val.split(self.sep)[self.index[fmt]]
                # we assume that the most common ones has been added...
                # and we do not want to add all sorts of rare values forever
                if len(self.factory[fmt]) < 10000:
                    self.factory[fmt][val] = res
                return res
            except:
                # if the format has not been handled before.
                if fmt not in self.index:
                    fields = fmt.split(self.sep)
                    if self.name in fields:
                        self.index[fmt] = fields.index(self.name)
                    else:
                        self.index[fmt] = None
                # try again
                try:
                    res = val.split(self.sep)[self.index[fmt]]
                    if len(self.factory[fmt]) < 10000:
                        self.factory[fmt][val] = res
                    return res
                # if still error
                except:
                    self.factory[fmt][val] = self.default
                    return self.default


__all_field_from_format = {}


def FieldFromFormat(name, sep=';', default=None):
    # this is a factory of __FieldFromFormat class
    #
    global __all_field_from_format
    if (name, sep, default) in __all_field_from_format:
        return __all_field_from_format[(name, sep, default)]
    else:
        obj = __FieldFromFormat(name, sep, default)
        __all_field_from_format[(name, sep, default)] = obj
        return obj


class VcfGenotype:

    def __init__(self, default=None):
        '''Define an extractor that extract genotype from a .vcf file'''
        #
        # FIXME: the current genotype handling mechanism cannot handle partial
        # missing data (e.g. ./1) and phase.
        #
        self.default = default
        self.map = {
            '0/0': default,
            '0|0': default,
            '0/1': ('1',),
            '1/0': ('1',),
            '0|1': ('1',),
            '1|0': ('1',),
            '1/1': ('2',),
            '1|1': ('2',),
            '0/2': (None, '1'),
            '2/0': (None, '1'),
            '0|2': (None, '1'),
            '2|0': (None, '1'),
            '1/2': ('-1', '-1'),
            '2/1': ('-1', '-1'),
            '1|2': ('-1', '-1'),
            '2|1': ('-1', '-1'),
            '2/2': (None, '2'),
            '2|2': (None, '2'),
            '0/3': (None, None, '1'),
            '3/0': (None, None, '1'),
            '0|3': (None, None, '1'),
            '3|0': (None, None, '1'),
            '1/3': ('-1', None, '-1'),
            '3/1': ('-1', None, '-1'),
            '1|3': ('-1', None, '-1'),
            '3|1': ('-1', None, '-1'),
            '2/3': (None, '-1', '-1'),
            '3/2': (None, '-1', '-1'),
            '2|3': (None, '-1', '-1'),
            '3|2': (None, '-1', '-1'),
            '3/3': (None, None, '2'),
            '3|3': (None, None, '2'),
            '0/4': (None, None, None, '1'),
            '4/0': (None, None, None, '1'),
            '0|4': (None, None, None, '1'),
            '4|0': (None, None, None, '1'),
            '1/4': ('-1', None, '-1'),
            '4/1': ('-1', None, '-1'),
            '1|4': ('-1', None, '-1'),
            '4|1': ('-1', None, '-1'),
            '2/4': (None, '-1', None, '-1'),
            '4/2': (None, '-1', None, '-1'),
            '2|4': (None, '-1', None, '-1'),
            '4|2': (None, '-1', None, '-1'),
            '3/4': (None, None, '-1', '-1'),
            '4/3': (None, None, '-1', '-1'),
            '3|4': (None, None, '-1', '-1'),
            '4|3': (None, None, '-1', '-1'),
            '4/4': (None, None, None, '2'),
            '4|4': (None, None, None, '2'),
            '0': default,
            '1': ('1',)
        }

    def __call__(self, item):
        # the most common and correct case...
        try:
            return self.map[item.partition(':')[0]]
        except KeyError:
            env.logger.debug('Genotype {} cannot be imported'.format(item))
            return None


class VcfGenoFromFormat:

    def __init__(self, default=None):
        '''Define an extractor that return genotype according to a format string.
        This is used to extract genotype from the format string of vcf files.
        '''
        self.fmt = '\t'
        self.idx = None
        self.default = default
        self.map = {
            '0/0': default,
            '0|0': default,
            '0/1': ('1',),
            '1/0': ('1',),
            '0|1': ('1',),
            '1|0': ('1',),
            '1/1': ('2',),
            '1|1': ('2',),
            '0/2': (None, '1'),
            '2/0': (None, '1'),
            '0|2': (None, '1'),
            '2|0': (None, '1'),
            '1/2': ('-1', '-1'),
            '2/1': ('-1', '-1'),
            '1|2': ('-1', '-1'),
            '2|1': ('-1', '-1'),
            '2/2': (None, '2'),
            '2|2': (None, '2'),
            '0/3': (None, None, '1'),
            '3/0': (None, None, '1'),
            '0|3': (None, None, '1'),
            '3|0': (None, None, '1'),
            '1/3': ('-1', None, '-1'),
            '3/1': ('-1', None, '-1'),
            '1|3': ('-1', None, '-1'),
            '3|1': ('-1', None, '-1'),
            '2/3': (None, '-1', '-1'),
            '3/2': (None, '-1', '-1'),
            '2|3': (None, '-1', '-1'),
            '3|2': (None, '-1', '-1'),
            '3/3': (None, None, '2'),
            '3|3': (None, None, '2'),
            '0/4': (None, None, None, '1'),
            '4/0': (None, None, None, '1'),
            '0|4': (None, None, None, '1'),
            '4|0': (None, None, None, '1'),
            '1/4': ('-1', None, '-1'),
            '4/1': ('-1', None, '-1'),
            '1|4': ('-1', None, '-1'),
            '4|1': ('-1', None, '-1'),
            '2/4': (None, '-1', None, '-1'),
            '4/2': (None, '-1', None, '-1'),
            '2|4': (None, '-1', None, '-1'),
            '4|2': (None, '-1', None, '-1'),
            '3/4': (None, None, '-1', '-1'),
            '4/3': (None, None, '-1', '-1'),
            '3|4': (None, None, '-1', '-1'),
            '4|3': (None, None, '-1', '-1'),
            '4/4': (None, None, None, '2'),
            '4|4': (None, None, None, '2'),
            '0': default,
            '1': ('1',)
        }

    def __call__(self, item):
        # the most common and correct case...
        try:
            if item[0][:2] == 'GT':
                return self.map[item[1].partition(':')[0]]
            elif item[0] != self.fmt:
                fmt, val = item
                self.fmt = fmt
                fields = fmt.split(':')
                if 'GT' in fields:
                    self.idx = fields.index('GT')
                    return self.map[val.split(':')[self.idx]]
                else:
                    self.idx = None
                    return self.default
            return self.map[item[1].split(
                ':', self.idx +
                1)[self.idx]] if self.idx is not None else self.default
        except KeyError:
            env.logger.debug('Genotype {} cannot be imported'.format(item))
            return None


class ExtractValue:

    def __init__(self, name, sep=';', default=None):
        '''Define an extractor that returns the value after name in one of the fields,
        and a default value if no such field is found. No space is allowed between
        delimiter and the name.'''
        self.name = name
        self.sep = sep
        self.pos = len(name)
        self.default = default

    def __call__(self, item):
        #
        # Using two partisions seems to be a tiny bit faster than
        # split and startswith
        #
        #for field in item.split(self.sep):
        #    if field.startswith(self.name):
        #        return field[self.pos:]
        if item.startswith(self.name):
            return item.split(self.sep, 1)[0][self.pos:]
        else:
            ss = item.partition(self.sep + self.name)
            return ss[2].partition(self.sep)[0] if ss[2] else self.default


class IncreaseBy:
    '''Increase passed value by a given number, will convert input to integer'''

    def __init__(self, inc=1):
        '''Adjust position'''
        self.inc = inc

    def __call__(self, item):
        return str(int(item) + self.inc) if item.isdigit() else None


class MapValue:
    '''Map value to another one, return the item itself if unmapped'''

    def __init__(self, map):
        self.map = map

    def __call__(self, item):
        try:
            return self.map[item]
        except:
            return item


class RemoveLeading:
    '''Remove specified leading string if the input string starts with it. Used
    for example to remove chr from inputs such as chr15'''

    def __init__(self, val):
        self.val = val
        self.vlen = len(val)

    def __call__(self, item):
        return item[self.vlen:] if item.startswith(self.val) else item


class EncodeGenotype:
    '''Encode 1/1, 1/2 etc to variant tools code'''

    def __init__(self, default=None):
        self.map = {
            '0/0': default,
            '0|0': default,
            '0/1': ('1',),
            '1/0': ('1',),
            '0|1': ('1',),
            '1|0': ('1',),
            '1/1': ('2',),
            '1|1': ('2',),
            '0/2': ('0', '1'),
            '2/0': ('0', '1'),
            '0|2': ('0', '1'),
            '2|0': ('0', '1'),
            '1/2': ('-1', '-1'),
            '2/1': ('-1', '-1'),
            '1|2': ('-1', '-1'),
            '2|1': ('-1', '-1'),
            '2/2': ('0', '2'),
            '2|2': ('0', '2'),
            '0': default,
            '1': ('1',)
        }

    def __call__(self, item):
        return self.map[item]


class Nullify:
    '''Change specified input value to NULL '''

    def __init__(self, val):
        self.val = val
        if type(self.val) == str:
            self.__call__ = self.nullify_single
        else:
            self.__call__ = self.nullify_multiple

    def nullify_single(self, item):
        return None if item == self.val else item

    def nullify_multiple(self, item):
        return None if item in self.val else item


class IgnoredRecord(Exception):

    def __init__(self, value=None):
        self.value = value

    def __str__(self):
        return repr(self.value)


class DiscardRecord:

    def __init__(self, val, keepMatched=False):
        self.val = val
        if hasattr(self.val, '__call__'):
            if keepMatched:
                self.__call__ = self.discard_unmatched_eval
            else:
                self.__call__ = self.discard_matched_eval
        else:
            if keepMatched:
                if type(self.val) == str:
                    self.__call__ = self.discard_unmatched_single
                else:
                    self.__call__ = self.discard_unmatched_multiple
            else:
                if type(self.val) == str:
                    self.__call__ = self.discard_matched_single
                else:
                    self.__call__ = self.discard_matched_multiple

    def discard_unmatched_eval(self, item):
        if self.val.__call__(item) is not True:
            raise IgnoredRecord()
        return item

    def discard_matched_eval(self, item):
        if self.val.__call__(item) is True:
            raise IgnoredRecord()
        return item

    def discard_unmatched_single(self, item):
        if item != self.val:
            raise IgnoredRecord()
        return item

    def discard_unmatched_multiple(self, item):
        if item not in self.val:
            raise IgnoredRecord()
        return item

    def discard_matched_single(self, item):
        if item == self.val:
            raise IgnoredRecord()
        return item

    def discard_matched_multiple(self, item):
        if item in self.val:
            raise IgnoredRecord()
        return item


__databases = {}
#
# lock the database until it is sure that indexes are created
__db_lock = Lock()


class _DatabaseQuerier:
    '''This query a field from an annotation database'''

    def __init__(self, cursor, name, res_field, cond_fields, default=None):
        '''Supose res_field is alt, cond_fields are chr,pos, this querier
        will get alt using query
          SELECT dbSNP.alt FROM dbSNP WHERE chr=VAL1 AND pos=VAL2
        '''
        self.default = default
        self.cur = cursor
        self.single_cond = len(cond_fields) == 1
        self.query = 'SELECT {} FROM {} WHERE {}'.format(
            res_field, name,
            ' AND '.join(['{}=?'.format(x) for x in cond_fields]))

    def __call__(self, item):
        if self.single_cond:
            self.cur.execute(self.query, (item,))
        else:
            self.cur.execute(self.query, item)
        res = self.cur.fetchall()
        #env.logger.error('{} {}, {}'.format(self.query, item, res))
        if len(res) == 1:
            return res[0][0]
        elif len(res) > 1:
            return tuple([x[0] for x in res])
        else:
            return self.default


def FieldFromDB(dbfile, res_field, cond_fields, default=None):
    global __databases
    global __db_lock
    __db_lock.acquire()
    if dbfile not in __databases:
        db = DatabaseEngine()
        if not os.path.isfile(os.path.expanduser(dbfile)):
            if os.path.isfile(
                    os.path.join(env.local_resource, 'annoDB', dbfile)):
                database_file = os.path.join(env.local_resource, 'annoDB',
                                             dbfile)
            else:
                raise ValueError('Database file {} does not exist locally or '
                                 'under resource directory'.format(dbfile))
        else:
            database_file = os.path.expanduser(dbfile)
        db.connect(database_file, readonly=True)
        cur = db.cursor()
        tables = db.tables()
        if not tables:
            raise ValueError(
                'Incorrect annotation database with tables {}'.format(
                    ', '.join(tables)))
        #
        try:
            name = [x for x in tables if x.endswith('_info')][0][:-5]
        except Exception as e:
            raise ValueError(
                'Incorrect database (missing info table): {}'.format(e))
        if not name in tables:
            raise ValueError(
                'Incorrect database (missing table {})'.format(name))
        if not name + '_field':
            raise ValueError('Incorrect database (missing field table)')
        #
        # we need to create indexes for the databases but we have to re-open the
        # database to make it writable. Because multiple processes might be used
        # we have to use a global lock.
        #
        for fld in cond_fields.split(','):
            if not db.hasIndex('{}_idx'.format(fld)):
                db.close()
                db.connect(database_file, readonly=False)
                env.logger.info(
                    'Creating index for field "{}" in database {}'.format(
                        fld, name))
                db.execute('CREATE INDEX {0}_idx ON {1} ({0} ASC);'.format(
                    fld, name))
                db.commit()
                db.close()
                db.connect(database_file, readonly=True)
                cur = db.cursor()
        __databases[dbfile] = (cur, name)
    __db_lock.release()
    return _DatabaseQuerier(__databases[dbfile][0], __databases[dbfile][1],
                            res_field, cond_fields.split(','), default)


class RefAtPos:
    '''This function returns the reference allele from a position'''

    def __init__(self, build):
        self.refGenome = RefGenome(build)

    def __call__(self, item):
        # enter chromosome and pos
        return self.refGenome.getBase(item[0], int(item[1]))


class AltAtPos:
    '''This function returns the alternative allele from two-given alleles, and reference allele from a position'''

    def __init__(self, build):
        self.refGenome = RefGenome(build)

    def __call__(self, item):
        # enter chromosome and pos
        ref = self.refGenome.getBase(item[0], int(item[1]))
        if item[2] == item[3]:
            raise ValueError('Identical alleles provided to AltAtPos')
        if ref == item[2]:
            return item[3]
        elif ref == item[3]:
            return item[2]
        else:
            raise ValueError(
                'Two non-ref alleles are provided: chr={}, pos={}, ref={}, observed {} {}'
                .format(item[0], item[1], ref, item[2], item[3]))


# this is a dictionary to save extractors for each file used
g_SeqExtractor = {}


def SeqAtLoc(filename):
    # return the same object for multiple instances of SeqAtLoc because
    # we do not want to read the fasta file multiple times
    if filename not in g_SeqExtractor:
        g_SeqExtractor[filename] = SequentialExtractor(filename)
    return g_SeqExtractor[filename]


class SequentialExtractor:

    def __init__(self, extractors):
        '''Define an extractor that calls a list of extractors. The string extracted from
        the first extractor will be passed to the second, and so on.'''
        self.extractors = []
        for e in extractors:
            if hasattr(e, '__call__'):
                self.extractors.append(e.__call__)
            else:
                self.extractors.append(e)

    def __call__(self, item):
        for e in self.extractors:
            # if multiple records are returned, apply to each of them
            if type(item) is tuple:
                if type(item[0]) is tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = tuple(e(x) for x in item)
            # if item is None or ''
            elif item is None or item == '':
                return item
            else:
                item = e(item)
        return item


class JoinFields:

    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep

    def __call__(self, item):
        try:
            if type(item) == str:
                return item
            else:
                return self.sep.join(item)
        except:
            return str(item)


class IfMulti:

    def __init__(self, ifFunc=None, elseFunc=None):
        if hasattr(ifFunc, '__call__'):
            self.ifFunc = ifFunc.__call__
        else:
            self.ifFunc = ifFunc
        if hasattr(elseFunc, '__call__'):
            self.elseFunc = elseFunc.__call__
        else:
            self.elseFunc = elseFunc

    def __call__(self, item):
        if type(item) == tuple:
            return item[0] if self.ifFunc is None else self.ifFunc(item)
        else:
            return item if self.elseFunc is None else self.elseFunc(item)


class JoinRecords:

    def __init__(self, sep=','):
        '''Define an extractor that returns all items in a field separated by
        specified delimiter. These items will lead to multiple records in
        the database.'''
        self.sep = sep

    def __call__(self, item):
        try:
            if type(item) == tuple or type(item) == list:
                return self.sep.join([str(x) for x in item])
            else:
                return str(item)
        except:
            return str(item)


class ValueOfNull:

    def __init__(self, val):
        self.val = val

    def __call__(self, item):
        return self.val if item in ('', None) else item


class Formatter:

    def __init__(self, fmt):
        self.fmt = fmt

    def __call__(self, item):
        try:
            return self.fmt.format(item)
        except:
            return str(item)


class PlainFormatter:

    def __init__(self, mode='vcf'):
        self.mode = mode

    def mergeVariants(self, item):
        #
        # This is UGLY code but it is not my fault, it is the vcf format's fault. :-(
        #
        # this is a special case for vcf. It is better to assign a special
        # adjuster but I put it here for backward compatibility of the .vcf
        # file.
        #
        # case 1:  X  10000  A  C
        #          X  10000  A  T
        #
        # merge to: X  10000  A  C,T
        #
        # case 2:  X  10000 TAC T  (from X  10001 AC -)
        #          X  10001 A  T
        #
        # merge to: X 10000 TAC  T,TT
        #
        # join each subitems
        uniq = [[x] for x in item[0].split('\t')]
        for x in item[1:]:
            for idx, y in enumerate(x.split('\t')):
                if idx >= len(uniq):
                    uniq.append([y])
                elif y and y not in uniq[idx]:
                    uniq[idx].append(y)
        # we can only handle case 2 with 2 items
        if len(uniq[1]) == 1 or len(uniq) != 5 or len(item) != 2 or len(
                uniq[0]) != 1 or len(uniq[3]) == 1:
            return '\t'.join([','.join(x) for x in uniq])
        else:
            # we have same chromosome, two position, two items, two ref, length 5
            try:
                min_pos = min([int(x) for x in uniq[1]])
                min_item = 0 if int(uniq[1][0]) == min_pos else 1
                ref = uniq[3][min_item]
                alt0 = item[0].split('\t')[4]
                alt1 = item[1].split('\t')[4]
                return '\t'.join([uniq[0][0], # chromosome shared
                        uniq[1][min_item],        # minimal pos
                        uniq[2][0],
                        ref,        # longer ref
                        #  X  10000 TAC T
                        #  X  10001 A  T
                        (alt0 + ',' + ref[:int(uniq[1][1]) - min_pos] + alt1) if min_item == 0 \
                        else (ref[:int(uniq[1][0]) - min_pos] + alt0 + ',' + alt1)])
            except:
                return '\t'.join([','.join(x) for x in uniq])

    def __call__(self, item):
        if type(item) == tuple:
            if '\t' in item[0]:
                return self.mergeVariants(item)
            else:
                # join unique ones
                uniq = [item[0]]
                for x in item[1:]:
                    if x and x not in uniq:
                        uniq.append(x)
                return ','.join([str(x) for x in uniq])
        else:
            return str(item)


class CSVFormatter:

    def __init__(self):
        pass

    def __call__(self, item):
        if type(item) == str:
            if not item:
                return ''
            elif '"' in item:
                return '"' + item.replace('"', '""') + '"'
            # quote all strings, because sometimes excel will treat them differently.
            else:
                return '"' + item + '"'
        elif type(item) == tuple:
            val = ','.join([str(x) for x in item])
            if '"' in val:
                return '"' + val.replace('"', '""') + '"'
            if ',' in val or '\n' in val:
                return '"' + val + '"'
            return val
        else:
            # not string...
            val = str(item)
            if '"' in val:
                return '"' + val.replace('"', '""') + '"'
            if ',' in val or '\n' in val:
                return '"' + val + '"'
            return val


rec_alleles = ['-', '-']


class InfoFormatter:

    def __init__(self, name, ignore=''):
        '''Output value as $name=val'''
        self.name = name
        self.ignore = ignore

    def __call__(self, item):
        return '' if item == self.ignore else '{}={}'.format(self.name, item)


class FlagFormatter:

    def __init__(self, name):
        '''Output value as $name=val'''
        self.name = name

    def __call__(self, item):
        return self.name if item else ''


class GenoFormatter:
    # representation for missing value is style dependent,
    # the default value None will cause each style to use its default value.
    def __init__(self, style='genotype', sep='\t', null='-', base=0):
        self.sep = sep
        self.null = null
        if style == 'numeric':
            self.missing = 'NA'
        elif style == 'genotype':
            self.missing = '.'
        elif style == 'plink':
            # PED format seems to use ACTG, and 0 for missing.
            # see http://www.sph.umich.edu/csg/abecasis/Merlin/tour/input_files.html
            self.missing = '0'
        elif style:
            self.missing = ''
        self.base = base
        #
        self.vcf_map = {
            # When we output a single sample, it can have genotype
            #   0, 1, 2, and (-1, -1)
            # as we imported from .vcf file.
            #
            0:
                '0/0',
            1:
                '0/1',
            2:
                '1/1',
            (-1, -1):
                '1/2',
            #
            # if one of the two alternative variants is filtered out, we can have
            # a single -1 variant. Note that we cannot yet import partial missing
            # data.
            -1:
                './1',
            #
            # when two or more samples are outputted, one sample might not have
            # any genotype for a variant, or more than one variant. There can
            # be more than one variant at a location though.
            #
            None:
                '.',
            (None, None):
                '.',
            (None, None, None):
                '.',
            #
            # Having two valid and complete genotypes for both variants is
            # as far as I can imagine, not possible. However, because it is
            # true that the variant is homogenous wildtype at both variants,
            # (0,0) is listed here.
            (0, 0):
                '0/0',
            #
            #(0,1): '0/1',
            #(0,2): '2/2',
            #
            # However, if there are two variants, one sample can have 0, 1, 2, and
            # the None for another variant (because that variant exists in other
            # samples). Note that the second item is for another variant
            (0, None):
                '0/0',
            (None, 0):
                '0/0',
            (1, None):
                '0/1',
            (None, 1):
                '0/2',
            (2, None):
                '1/1',
            (None, 2):
                '2/2',
            # the single -1 case
            (-1, None):
                './1',
            (None, -1):
                './2',
            #
            # the same goes to the case with three alternative alleles
            (0, None, None):
                '0/0',
            (None, 0, None):
                '0/0',
            (None, None, 0):
                '0/0',
            (1, None, None):
                '0/1',
            (None, 1, None):
                '0/2',
            (None, None, 1):
                '0/3',
            (2, None, None):
                '1/1',
            (None, 2, None):
                '2/2',
            (None, None, 2):
                '3/3',
            # the -1 case is more complicated because there can be one or two -1
            (-1, None, None):
                './1',
            (None, -1, None):
                './2',
            (None, None, -1):
                './3',
            (-1, -1, None):
                '1/2',
            (-1, None, -1):
                '1/3',
            (None, -1, -1):
                '2/3',
        }
        #
        self.numeric_map = {
            # number of non-wildtype alleles
            # FIXME: We are treating multiple alternative alleles as the same
            # non-wildtype allele, which might be wrong for some formats.
            0: 0,
            1: 1,
            2: 2,
            (-1, -1): 2,  # two DIFFERENT alternative alleles
            #
            -1: 1,
            #
            None: None,
            (None, None): None,
            (None, None, None): None,
            (None, None, None, None): None,
            #
            (0, 0): 0,
            #
            (0, None): 0,
            (None, 0): 0,
            (1, None): 1,
            (None, 1): 1,
            (2, None): 2,
            (None, 2): 2,
            #
            (-1, None): 1,
            (None, -1): 1,
            #
            (0, None, None): 0,
            (None, 0, None): 0,
            (None, None, 0): 0,
            (1, None, None): 1,
            (None, 1, None): 1,
            (None, None, 1): 1,
            (2, None, None): 2,
            (None, 2, None): 2,
            (None, None, 2): 2,
            #
            (-1, None, None): 1,
            (None, -1, None): 1,
            (None, None, -1): 1,
            (-1, -1, None): 2,
            (-1, None, -1): 2,
            (None, -1, -1): 2,
        }
        #
        if style == 'genotype':
            self.__call__ = self.fmt_genotype
        elif style == 'numeric':
            self.__call__ = self.fmt_numeric
        elif style == 'vcf':
            self.__call__ = self.fmt_vcf
        elif style == 'plink':
            self.__call__ = self.fmt_plink
        else:
            raise ValueError('Only genotype and numeric styles are allowed')

    def fmt_genotype(self, item):
        global rec_alleles
        if type(item) == float and item is not None:
            item = int(item)
        if type(item) == int:
            # single genotype case
            ref = self.null if rec_alleles[0] == '-' else rec_alleles[0]
            alt = self.null if rec_alleles[1] == '-' else rec_alleles[1]
            #
            # 0, 1, 2, -1
            if item == 0:
                return ref + self.sep + ref
            elif item == 1:
                return ref + self.sep + alt
            elif item == 2:
                return alt + self.sep + alt
            else:
                return alt + self.sep + self.null
        elif type(item) == tuple:
            # Two aleternative alleles
            if item == (-1, -1):
                return (self.null if rec_alleles[1][0] == '-' else
                        rec_alleles[1][0]) + self.sep + (
                            self.null
                            if rec_alleles[1][1] == '-' else rec_alleles[1][1])
            elif len(item) > 1 and item.count(item[0]) == len(item):
                # assume duplicate entry caused by annotation database
                ref = self.null if rec_alleles[0][0] == '-' else rec_alleles[0][
                    0]
                alt = self.null if rec_alleles[1][0] == '-' else rec_alleles[1][
                    0]
                if item[0] == 1:
                    return ref + self.sep + alt
                elif item[0] == 2:
                    return alt + self.sep + alt
                else:
                    return ref + self.sep + ref
            #
            # the cases for (None, None) and (None, None, None)
            elif all([x is None for x in item]):
                return self.missing + self.sep + self.missing
            else:
                raise ValueError(
                    'Failed to export genotype {} with ref {} and alt {}'
                    .format(item, rec_alleles[0], rec_alleles[1]))
        elif item is None:
            return self.missing + self.sep + self.missing
        else:
            raise ValueError('Failed to export genotype {}'.format(item))

    def fmt_plink(self, item):
        global rec_alleles
        if type(item) == float and item is not None:
            item = int(item)
        if type(item) == int:
            # single genotype case
            ref = self.null if rec_alleles[0] == '-' else rec_alleles[0]
            alt = self.null if rec_alleles[1] == '-' else rec_alleles[1]
            #
            # 0, 1, 2, -1
            if item == 0:
                return ref + self.sep + ref
            elif item == 1:
                return ref + self.sep + alt
            elif item == 2:
                return alt + self.sep + alt
            else:  # (./A), considered as missing
                return self.missing + self.sep + self.missing
        # multi-allele case, ignore
        elif type(item) == tuple:
            raise ValueError(
                'plink format cannot handle multiple alleles {}'.format(item))
        elif item is None:
            return self.missing + self.sep + self.missing
        else:
            raise ValueError('Failed to export genotype {}'.format(item))

    def fmt_numeric(self, item):
        try:
            cnt = self.numeric_map[item]
            if cnt is None:
                return self.missing
            else:
                return str(cnt + self.base)
        except Exception as e:
            raise ValueError(
                'Failed to export genotype {} in numeric style: {}'.format(
                    item, e))

    def fmt_vcf(self, item):
        try:
            return self.vcf_map[item]
        except:
            raise ValueError(
                'Failed to export genotype {} in vcf style with ref {} and alt {}.'
                .format(item, rec_alleles[0], rec_alleles[1]))


class Constant:

    def __init__(self, val=''):
        self.val = val

    def __call__(self, item):
        return self.val


class SequentialCollector:

    def __init__(self, extractors):
        '''Define an extractor that calls a list of extractors. The string extracted from
        the first extractor will be passed to the second, and so on.'''
        self.extractors = []
        for e in extractors:
            if hasattr(e, '__call__'):
                self.extractors.append(e.__call__)
            else:
                self.extractors.append(e)

    def __call__(self, item):
        for e in self.extractors:
            # if multiple records are returned, apply to each of them
            if type(item) is tuple:
                if type(item[0]) is tuple:
                    raise ValueError('Nested vector extracted is not allowed')
                item = [e(x) for x in item]
            else:
                item = e(item)
        return item


class BatchWriter:
    '''write text to file in batches (#lines)'''

    def __init__(self, fn, batch=1000):
        self.batch = batch
        if os.path.exists(fn):
            os.remove(fn)
        self.fs = open(fn, 'a')
        self.counter = 0
        self.swap = ''

    def write(self, line):
        if line is None:
            self.fs.write(self.swap)
            self.fs.close()
        else:
            self.swap += line
            self.counter += 1
            if self.counter == self.batch:
                # time to write
                self.fs.write(self.swap)
                self.counter = 0
                self.swap = ''


class PlinkBinaryToVariants:
    """
    Class to write a PLINK BED format genotype dataset (.bed, .fam and .bim files)
    to vtools compatible variant format

    c.f., http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

    @param    dataset: path + prefix for a .bed, .fam and .bim without extensions
    @param      build:
    @param     logger:

    Note: implementation of this class is based on libplinkio
    https://bitbucket.org/mattias_franberg/libplinkio

    The SNP will be encoded as follows:
    * 0 - Homozygous major
    * 1 - Heterozygous
    * 2 - Homozygous minor
    * 3 - Missing value

    Output of this program will retain such coding if the major allele is reference genotype
    in reference human genome; otherwise will re-code it by:
    * 0 --> 2
    * 2 --> 0

    WARNING
    -------
    The allele1 and allele2 in the *.bim file does not tell which is major and which is minor
    so treating allele1 as major and allele2 as minor will result in completely swapped
    genotype coding if it is not the case !!

    self.determineMajorAllele() Attempts to resolve this issue
    """

    def __init__(self, dataset, build, chrom_namemap={}):
        # check file path
        for ext in ['.fam', '.bed', '.bim']:
            if not os.path.exists(dataset + ext):
                raise RuntimeError('Cannot find file {0}'.format(dataset + ext))
        self.dataset = dataset
        self.build = build
        self.cur = PlinkFile(self.dataset)
        # a list of sample names (sample ID's in .fam file)'''
        self.samples = [x.iid for x in self.cur.get_samples()]
        if None in self.samples:
            raise ValueError("Cannot read sample ID from malformed '{0}.fam' file.".\
                             format(self.dataset))
        # iterator for variants info: chr, pos, allele1, allele2
        self.variants = it.chain(self.cur.get_loci())
        # reference genome object and ATCG dictionary
        self.hgref = RefGenome(build)
        self.CSTRANDS = {'A': 'T', 'G': 'C', 'T': 'A', 'C': 'G'}
        # status 0 for not flip, 1 for flip, -1 for bad match
        self.status = 0
        # chromosome naming convention map
        self.cmap = chrom_namemap
        # variant writer
        self.variant_writer = None
        self.data_writer = None

    def initWriter(self, ofile):
        self.variant_writer = BatchWriter(
            self.dataset + ".{0}adjusted".format(self.build), batch=1000)
        self.data_writer = BatchWriter(ofile, batch=5)

    def getValidatedLocusGenotype(self, chrom, pos, allele1, allele2, geno_cur):
        '''Use cgatools to obtain validated genotype for given locus.

        Assume allele1 is major allele and allele2 is minor allele
        - check which of the two input allele is reference allele
        - if none is reference, use the alternative strand and check again
        - if the locus is validated (self.status >= 0)
          - flip ref/alt coding if the minor allele is found to be reference allele
          - output its genotype
          otherwise return None
        @return a locus genotypes string
        '''
        if chrom in self.cmap:
            chrom = self.cmap[chrom]
        try:
            ref = self.hgref.getBase(chrom, pos)
        except Exception:
            env.logger.warning(
                'Cannot find genomic coordinate {0}:{1} in reference genome {2}. '
                'Input variant is ignored'.format(chrom, pos, self.build))
            self.status = -1
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(
                    chrom, pos, 0, 0))
            return None
        self.status, strand, allele1, allele2 = self._matchref(
            ref, allele1, allele2)
        if self.status < 0:
            if self.status == -2:
                env.logger.warning('All genotypes for variant "{0}:{1}" are missing'.\
                                    format(chrom, pos))
            else:
                env.logger.warning('Variant "{0}:{1} {2} {3}" failed '
                                   'to match reference genome {4}/(A,T,C,G)'.\
                                    format(chrom, pos, allele1, allele2, ref))
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(
                    chrom, pos, 0, 0))
            return None
        elif self.status == 0:
            if strand:
                env.logger.debug('Use alternative strand for {0}:{1}'.format(
                    chrom, pos))
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(
                    chrom, pos, allele1, allele2))
            return ','.join([chrom, str(pos), allele1, allele2
                            ]) + ',' + str(geno_cur)
        else:
            # have to flip the genotypes coding
            if strand:
                env.logger.debug('Use alternative strand for {0}:{1}'.format(
                    chrom, pos))
            # env.logger.debug('Allele coding flipped for {0}:{1}'.format(chrom, pos))
            # Very time consuming compare to not flipping the genotype codes
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(
                    chrom, pos, allele2, allele1))
            return ','.join([chrom, str(pos), allele2, allele1]) + ',' + \
                ','.join([str(x) if x == 3 or x == 'E' else str(2 - x) for x in geno_cur])

    def getLociCounts(self):
        # FIXME: not efficient
        return len(self.cur.get_loci())

    def getHeader(self):
        '''a line of headers for the output text file'''
        return ','.join(['#chr', 'pos', 'ref', 'alt'] + self.samples)

    def getLine(self, viter=None, giter=None, which_major=1):
        '''a line of validated genotypes for a locus'''
        if viter is None:
            viter = self.variants
        if giter is None:
            giter = self.cur
        try:
            locus = next(viter)
            genotypes = next(giter)
        except StopIteration:
            giter.close()
            return False, None
        except Exception as e:
            env.logger.error('Failed to retrieve locus {0}:{1} '
                             '(plinkio error "{2}")'.format(
                                 locus.chromosome, locus.bp_position, e))
            return True, None
        if which_major == 1:
            # allele 1 is the major allele
            return True, self.getValidatedLocusGenotype(
                str(locus.chromosome), int(locus.bp_position),
                locus.allele1.upper(), locus.allele2.upper(), genotypes)
        else:
            # allele 2 is the major allele
            return True, self.getValidatedLocusGenotype(
                str(locus.chromosome), int(locus.bp_position),
                locus.allele2.upper(), locus.allele1.upper(), genotypes)

    def determineMajorAllele(self, n=1000):
        '''The logic here is that for the first n loci we
        - assume the allele1 in *.bim is major allele
        - try to map to hg reference and see how many genotype coding have to be flipped
        - if too many (over half) have to be flipped we conclude that allele2 should be major allele
        @return: a guess of major allele based on the first n sample loci. -1 for bad matching
        '''
        # new temporary connection
        cur = PlinkFile(self.dataset)
        variants = it.chain(cur.get_loci())
        # counters
        m_ones = zeros = ones = 0
        i = 0
        while i < n:
            self.status = 0
            # self.status will be updated
            flag, line = self.getLine(viter=variants, giter=cur)
            if not flag:
                # end of line
                n = i
                break
            if self.status == 0:
                # not flipped
                zeros += 1
            elif self.status == 1:
                # flipped
                ones += 1
            elif self.status == -1:
                # bad match
                m_ones += 1
            else:
                # status = -2, monomorphic site
                # ignore this test
                if self.status == -2:
                    continue
            i += 1
        # check if so many are negative values
        if m_ones > float(n) / 2.0:
            return -9
        if ones > float(n - m_ones) / 2.0:
            # allele 2 seems to be major allele
            return 2
        else:
            return 1

    def _matchref(self, ref, major, minor):
        '''try best to match reference allele
        @param    ref: reference base coding
        @param  major: major allele coding
        @param  minor: minor allele coding
        @return self.status (0 for no need to flip, 1 for having to flip)
        @return strand (0 for original, 1 for alternative),
        @return major and minor alleles (might be from alternative strand)

        Notice
        ------
        In *.bim file allele1 or allele2 (major/minor) can be 0 when only one allele is found
        in data. Such loci is not considered a variant site and will be ignored
        '''
        # monomorphic, missing, or invalid coding
        if major not in ['A', 'T', 'C', 'G', '0'
                        ] or minor not in ['A', 'T', 'C', 'G', '0']:
            return -1, 0, major, minor
        if major == '0' and minor == '0':
            return -2, 0, major, minor
        if major == '0' and minor != '0':
            major = minor
        if major != '0' and minor == '0':
            minor = major
        status = strand = 0
        if ref in [major, minor]:
            # allele found, determine flip status
            status = 0 if ref == major else 1
        else:
            # allele not found, use alternative strand
            strand = 1
            major = self.CSTRANDS[major]
            minor = self.CSTRANDS[minor]
            # allele still not found
            if ref not in [major, minor]:
                status = -1
            # allele found, determine flip status
            else:
                status = 0 if ref == major else 1
        return status, strand, major, minor


#
#
# Preprocessors of input files
# They will convert input files to intermediate variant based text files for easy import
#
#


class Preprocessor:

    def __init__(self):
        '''Base preprocessor class that converts $files
        and write intermediate output files to $outdir/file.format'''
        pass

    def convert(self, files, output_files):
        for item, ofile in zip(files, output_files):
            env.logger.info('Convert {} to {}'.format(item, ofile))


class Dos2Unix(Preprocessor):

    def __init__(self):
        Preprocessor.__init__(self)

    def convert(self, files, output_files):
        for ifile, ofile in zip(files, output_files):
            env.logger.info(
                'Converting {} with \r newline charater to unix format.'.format(
                    ifile))
            with open(ifile, 'rU') as input, open(ofile, 'w') as output:
                for line in input:
                    output.write(line)


class PlinkConverter(Preprocessor):

    def __init__(self, build, chrom_namemap={}):
        Preprocessor.__init__(self)
        self.build = build
        self.cmap = chrom_namemap

    def convert(self, files, output_files):
        for item, ofile in zip(files, output_files):
            if os.path.exists(item + ".bed"):
                self.decode_plink(
                    PlinkBinaryToVariants(item, self.build, self.cmap), ofile)
            else:
                import glob
                files = '/'.join([x for x in glob.glob(item + '*')])
                if files:
                    supported = ['*.bed/*.bim/*.fam']
                    raise ValueError("Unsupported input file '{}' (supported file types are {})".\
                                         format(files, ';'.join(supported)))
                else:
                    raise ValueError(
                        "Cannot find input files '{}'".format(item + '*'))

    def decode_plink(self, p2v, ofile, n=1000):
        '''decode plink data from p2v object and output to ofile'''
        env.logger.info("Determining major/minor allele from data")
        # check major allele
        which_major = p2v.determineMajorAllele(n)
        # raise on bad match
        if which_major == -9:
            raise ValueError(
                'Invalid dataset {0}: too many unmatched loci to reference genome {1}. '
                'Perhaps you specified the wrong build, or have too many unsupported allele '
                'types (not A/T/C/G, e.g, indels I/D) in BED file which you have to remove before '
                'import'.format(p2v.dataset, p2v.build))
        env.logger.debug("allele{} is major allele".format(which_major))
        # output
        nloci = p2v.getLociCounts()
        batch = int(nloci / 100) + 1
        prog = ProgressBar('Decoding {0}'.format(p2v.dataset), nloci)
        p2v.initWriter(ofile)
        p2v.data_writer.write(p2v.getHeader() + '\n')
        p2v.variant_writer.write("#chr\tpos\tref\talt\n")
        count = 0
        while True:
            flag, line = p2v.getLine(which_major=which_major)
            count += 1
            if not flag:
                prog.done()
                p2v.data_writer.write(None)
                p2v.variant_writer.write(None)
                break
            else:
                if line is not None:
                    p2v.data_writer.write(line + '\n')
                if count % batch == 0 and count > batch:
                    prog.update(count)
