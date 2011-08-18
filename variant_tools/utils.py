#!/usr/bin/env python
#
# $File: utils.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2011 Bo Peng (bpeng@mdanderson.org)
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

'''
This module provides utility functions and a database engine that works
for both sqlite3 and mysql.
'''
import os
import sys
import glob
import logging
import subprocess
import urllib
import urlparse
import sqlite3
import getpass
import time
import tokenize
import cStringIO
import gzip
import sched


runOptions = {
    'verbosity': '1'
}

SQL_KEYWORDS = set([
    'ADD', 'ALL', 'ALTER', 'ANALYZE', 'AND', 'AS', 'ASC', 'ASENSITIVE', 'BEFORE',
    'BETWEEN', 'BIGINT', 'BINARY', 'BLOB', 'BOTH', 'BY', 'CALL', 'CASCADE', 'CASE',
    'CHANGE', 'CHAR', 'CHARACTER', 'CHECK', 'COLLATE', 'COLUMN', 'CONDITION',
    'CONSTRAINT', 'CONTINUE', 'CONVERT', 'CREATE', 'CROSS', 'CURRENT_DATE',
    'CURRENT_TIME', 'CURRENT_TIMESTAMP', 'CURRENT_USER', 'CURSOR', 'DATABASE',
    'DATABASES', 'DAY_HOUR', 'DAY_MICROSECOND', 'DAY_MINUTE', 'DAY_SECOND', 'DEC',
    'DECIMAL', 'DECLARE', 'DEFAULT', 'DELAYED', 'DELETE', 'DESC',
    'DESCRIBE', 'DETERMINISTIC', 'DISTINCT', 'DISTINCTROW', 'DIV', 'DOUBLE',
    'DROP', 'DUAL', 'EACH', 'ELSE', 'ELSEIF', 'ENCLOSED', 'ESCAPED', 'EXISTS',
    'EXIT', 'EXPLAIN', 'FALSE', 'FETCH', 'FLOAT', 'FLOAT4', 'FLOAT8', 'FOR',
    'FORCE', 'FOREIGN', 'FROM', 'FULLTEXT', 'GRANT', 'GROUP', 'HAVING', 'HIGH_PRIORITY',
    'HOUR_MICROSECOND', 'HOUR_MINUTE', 'HOUR_SECOND', 'IF', 'IGNORE', 'IN',
    'INDEX', 'INFILE', 'INNER', 'INOUT', 'INSENSITIVE', 'INSERT',
    'INT', 'INT1', 'INT2', 'INT3', 'INT4', 'INT8', 'INTEGER', 'INTERVAL', 'INTO',
    'IS', 'ITERATE', 'JOIN', 'KEY', 'KEYS', 'KILL', 'LEADING', 'LEAVE', 'LEFT',
    'LIKE', 'LIMIT', 'LINES', 'LOAD', 'LOCALTIME', 'LOCALTIMESTAMP',
    'LOCK', 'LONG', 'LONGBLOB', 'LONGTEXT', 'LOOP', 'LOW_PRIORITY', 'MATCH',
    'MEDIUMBLOB', 'MEDIUMINT', 'MEDIUMTEXT', 'MIDDLEINT', 'MINUTE_MICROSECOND',
    'MINUTE_SECOND', 'MOD', 'MODIFIES', 'NATURAL', 'NOT', 'NO_WRITE_TO_BINLOG',
    'NULL', 'NUMERIC', 'ON', 'OPTIMIZE', 'OPTION', 'OPTIONALLY', 'OR',
    'ORDER', 'OUT', 'OUTER', 'OUTFILE', 'PRECISION', 'PRIMARY', 'PROCEDURE',
    'PURGE', 'READ', 'READS', 'REAL', 'REFERENCES', 'REGEXP', 'RELEASE',
    'RENAME', 'REPEAT', 'REPLACE', 'REQUIRE', 'RESTRICT', 'RETURN',
    'REVOKE', 'RIGHT', 'RLIKE', 'SCHEMA', 'SCHEMAS', 'SECOND_MICROSECOND',
    'SELECT', 'SENSITIVE', 'SEPARATOR', 'SET', 'SHOW', 'SMALLINT',
    'SONAME', 'SPATIAL', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE',
    'SQLWARNING', 'SQL_BIG_RESULT', 'SQL_CALC_FOUND_ROWS', 'SQL_SMALL_RESULT',
    'SSL', 'STARTING', 'STRAIGHT_JOIN', 'TABLE', 'TERMINATED',
    'THEN', 'TINYBLOB', 'TINYINT', 'TINYTEXT', 'TO', 'TRAILING',
    'TRIGGER', 'TRUE', 'UNDO', 'UNION', 'UNIQUE', 'UNLOCK', 'UNSIGNED',
    'UPDATE', 'USAGE', 'USE', 'USING', 'UTC_DATE', 'UTC_TIME', 'UTC_TIMESTAMP', 'VALUES',
    'VARBINARY', 'VARCHAR', 'VARCHARACTER', 'VARYING', 'WHEN', 'WHERE', 'WHILE',
    'WITH', 'WRITE', 'XOR', 'YEAR_MONTH', 'ZEROFILL', 'ASENSITIVE', 'CALL', 'CONDITION',
    'CONNECTION', 'CONTINUE', 'CURSOR', 'DECLARE', 'DETERMINISTIC', 'EACH',
    'ELSEIF', 'EXIT', 'FETCH', 'GOTO', 'INOUT', 'INSENSITIVE', 'ITERATE', 'LABEL', 'LEAVE',
    'LOOP', 'MODIFIES', 'OUT', 'READS', 'RELEASE', 'REPEAT', 'RETURN', 'SCHEMA', 'SCHEMAS',
    'SENSITIVE', 'SPECIFIC', 'SQL', 'SQLEXCEPTION', 'SQLSTATE', 'SQLWARNING', 'TRIGGER',
    'UNDO', 'UPGRADE', 'WHILE', 'ABS', 'ACOS', 'ADDDATE', 'ADDTIME', 'ASCII', 'ASIN',
    'ATAN', 'AVG', 'BETWEEN', 'AND', 'BINARY', 'BIN', 'BIT_AND',
    'BIT_OR', 'CASE', 'CAST', 'CEIL', 'CHAR', 'CHARSET', 'CONCAT', 'CONV', 'COS', 'COT',
    'COUNT', 'DATE', 'DAY', 'DIV', 'EXP', 'IS', 'LIKE', 'MAX', 'MIN', 'MOD', 'MONTH',
    'LOG', 'POW', 'SIN', 'SLEEP', 'SORT', 'STD', 'VALUES', 'SUM'
])

def setOptions(verbosity=None):
    if verbosity is not None:
        runOptions['verbosity'] = verbosity

#
# Utility functions
#
def lineCount(filename):
    '''Estimate the number of lines using file size and line size. This
    function does not attemp to calculate line count exactly because files
    handled by variant tools can be huge. '''
    totalSize = os.path.getsize(filename)
    if totalSize < 100000:
        # small file, read the number of lines directly
        if filename.endswith('.gz'):
            return len(gzip.open(filename).readlines())
        else:
            return len(open(filename).readlines())
    elif filename.endswith('.gz'):
        input = gzip.open(filename)
        input.seek(1000, 0)
        content = input.read(99000)
        input.close()
        lineCount = len(content.split('\n'))
        input.close()
        # assuming an arbitrary compression ratio of 5. :-)
        return int(lineCount * (5 * totalSize / 99000.))
    else:
        input = open(filename)
        # count from the back because they tend to be records
        # with consistent size
        input.seek(-99000, 2)
        content = input.read()
        input.close()
        lineCount = len(content.split('\n'))
        input.close()
        return int(lineCount * (totalSize / 99000.))

def typeOfValues(vals):
    '''Figure out type of values and return INT, FLOAT or VARCHAR(maxLength)'''
    if len(vals) == 0:
        # a good default value?
        return 'VARCHAR(10)'
    try:
        map(int, vals)
        return 'INT'
    except:
        try:
            map(float, vals)
            return 'FLOAT'
        except:
            return 'VARCHAR({})'.format(max([len(x) for x in vals]))

class delayedAction:
    '''Call the passed function with param after a few seconds. It is most often 
    used to display certain message only if an action takes a long time.

        action = delayedAction(self.logger.info, 'This might take a while', 5)
        some_action_that_might_take_a_while
        del action

    if the action finishes very quick, the message will not be displayed.    
    '''
    def __init__(self, func, param, delay=5):
        self.scheduler = sched.scheduler(time.time, time.sleep)
        self.event = self.scheduler.enter(delay, 1, func, (param,))
        self.scheduler.run()

    def __del__(self):
        self.scheduler.cancel(self.event)

from array import array
try:
    from fcntl import ioctl
    import termios
except ImportError:
    pass
import signal

class ProgressBar:
    '''A text-based progress bar'''
    def __init__(self, message, totalCount = None):
        if runOptions['verbosity'] == '0':
            self.update = self.empty
            self.curlUpdate = self.empty
            self.urllibUpdate = self.empty
            self.sqliteUpdate = self.empty
            self.outputProgress = self.empty
            self.done = self.empty
            return
        self.message = message
        self.count = 0
        self.totalCount = totalCount
        self.start_time = None
        self.last_time = None
        try:
            self.handle_resize(None,None)
            signal.signal(signal.SIGWINCH, self.handle_resize)
            self.signal_set = True
        except:
            self.term_width = 79
        self.outputProgress()

    def empty(self, *args, **kwargs):
        return

    def handle_resize(self, signum, frame):
        # this is borrowed from python progressbar module
        h,w = array('h', ioctl(sys.stderr, termios.TIOCGWINSZ, '\0'*8))[:2]
        self.term_width = w

    def update(self, count):
        '''finished count jobs'''
        self.count = count
        self.outputProgress()

    def curlUpdate(self, total, existing, upload_t, upload_d):
        '''Update called from pycurl'''
        self.count = existing
        self.totalCount = total
        self.outputProgress()
        
    def urllibUpdate(self, count, blockSize, totalSize):
        '''Update called from urllib'''
        self.count = count * blockSize
        self.totalCount = totalSize
        self.outputProgress()

    def sqliteUpdate(self):
        self.count += 1
        if self.count % 1000 == 0:
            self.outputProgress()

    def outputProgress(self, done=False):
        '''Output progress'''
        if not self.start_time:
            self.start_time = time.time()
            self.last_time = self.start_time
        cur_time = time.time()
        # stop update progress bar more than once per second.
        if self.count > 0 and self.count != self.totalCount and cur_time - self.last_time < 1:
            return
        #
        msg = ['', '', '', '', '', '']
        # message
        msg[0] = self.message + ':'
        self.last_time = cur_time
        second_elapsed = cur_time - self.start_time
        cps = 0 if second_elapsed < 0.0001 else self.count / second_elapsed
        # speed
        if cps > 1000000:
            msg[4] = ' {:.1f}M/s'.format(cps/1000000)
        elif cps > 1000:
            msg[4] = ' {:.1f}K/s'.format(cps/1000)
        else:
            msg[4] = ' {:.1f}/s'.format(cps/1000)
        # estimated time left
        if done:
            msg[5] = time.strftime(' in %H:%M:%S', time.gmtime(second_elapsed))
        elif self.totalCount:
            perc = min(1, float(self.count) / self.totalCount)
            time_left = (second_elapsed / perc * (1 - perc)) if perc > 0 else 0
            msg[5] += time.strftime(' in %H:%M:%S', time.gmtime(time_left))
        # percentage / progress
        msg[3] = ' {:,}'.format(int(self.count))
        if self.totalCount:
            # percentage
            perc = min(1, float(self.count) / self.totalCount)
            msg[1] = ' {:5.1f}%'.format(perc * 100)
            width = self.term_width - len(msg[0]) - len(msg[1]) - len(msg[3]) - len(msg[4]) - len(msg[5])
            if width > 5:
                front = int(perc * (width - 5))
                back = width - 5 - front
                msg[2] = ' [{}>{}]'.format('=' * front,  ' ' * back)
        # use stderr to avoid messing up process output
        sys.stderr.write('\r' + ''.join(msg))
        # stderr is not buffer... so this is not needed
        #sys.stderr.flush()

    def done(self):
        '''Finish, output a new line'''
        if self.totalCount:
            self.count = self.totalCount
        self.outputProgress(done=True)
        sys.stderr.write('\n')
        sys.stderr.flush()

#
# Well, it is not easy to do reliable download
# 
def downloadFile(URL, dest_dir = None):
    '''Download file from URL to filename.'''
    filename = os.path.split(urlparse.urlsplit(URL).path)[-1]
    dest = filename if dest_dir is None else os.path.join(dest_dir, filename)
    # use libcurl? Recommended but not always available
    try:
        import pycurl
        prog = ProgressBar(filename)
        with open(dest, 'wb') as f:
            c = pycurl.Curl()
            c.setopt(pycurl.URL, URL)
            c.setopt(pycurl.WRITEFUNCTION, f.write)
            c.setopt(pycurl.NOPROGRESS, False)
            c.setopt(pycurl.PROGRESSFUNCTION, prog.curlUpdate)
            c.perform()
        prog.done()
        if os.path.isfile(dest):
            return dest
        else:
            raise RuntimeError('Failed to download {} using pycurl'.format(URL))
    except ImportError:
        # no pycurl module
        pass
    # use wget? Almost universally available under linux
    try:
        p = subprocess.Popen(['wget', '-O', dest, URL])
        p.wait()
        if os.path.isfile(dest):
            return dest
        else:
            raise RuntimeError('Failed to download {} using wget'.format(URL))
    except OSError:
        # no wget command
        pass
    #
    # use python urllib?
    prog = ProgressBar(filename)
    urllib.urlretrieve(URL, dest, reporthook=prog.urllibUpdate)
    prog.done()
    # all methods failed.
    if os.path.isfile(dest):
        return dest
    else:
        raise RuntimeError('Failed to download {}'.format(URL))
    
#
#
# Database engine
#

class DatabaseEngine:
    '''variant tools can make use of two database engines. One is mysql, the
    other is sqlite3. This class wraps around their DB API and provides an
    unified interface.'''
    def __init__(self, engine='sqlite3', batch=10000, **kwargs):
        '''
        engine
            Database engine, can be mysql or sqlite3 (default)
        batch
            Number of query per transaction. Larger number usually leads to better
            performance but requires more system resource.

        Additional keyword parameters such as 'host', 'user' and 'passwd' are
        passed directly to a MySQL database engine.
        '''
        self.engine = engine
        self.batch = batch
        # saved in case a new connection is needed
        self.connectionParams = kwargs
        # not connected to any database for now
        self.dbName = None
        if self.engine == 'mysql':
            import MySQLdb
            self.PH = '%s'
            self.AI = 'AUTO_INCREMENT'
            self.database = MySQLdb.connect(host=kwargs.get('host', 'localhost'),
                user=kwargs.get('user', getpass.getuser()),
                passwd=kwargs.get('passwd'))
        else:
            self.PH = '?'
            self.AI = 'AUTOINCREMENT'
            self.database = None

    #
    # Connection
    #
    def newConnection(self):
        '''Create a new connection from existing configuration'''
        return DatabaseEngine(engine=self.engine, batch=self.batch,
            **self.connectionParams)
        
    def connect(self, db):
        '''Connect to a database'''
        if self.engine == 'mysql':
            if '.' in db or os.sep in db:
                raise ValueError('Invalid database name: {}'.format(db))
            self.dbName = db
            cur = self.database.cursor()
            if not self.hasDatabase(self.dbName):
                cur.execute('CREATE DATABASE {};'.format(self.dbName))
            cur.execute('USE {};'.format(self.dbName))
        else:
            self.dbName = db if (db.endswith('.proj') or db.endswith('.DB')) else db + '.DB'
            self.database = sqlite3.connect(self.dbName)
            # set default cache size to a larger number to improve query performance
            cur = self.database.cursor()
            cur.execute('PRAGMA default_cache_size=2000;')
            self.database.commit()

    def attach(self, db):
        '''Attach another database to this one. Only needed by sqlite'''
        if self.engine == 'mysql':
            # create the database if needed
            if not self.hasDatabase(db):
                self.execute('CREATE DATABASE {};'.format(db))
            return
        if db.endswith('.DB'):
            self.execute('''ATTACH DATABASE '{0}' as {1};'''.format(
                db, os.path.split(db)[-1].split('.')[0]))
        else:
            self.execute('''ATTACH DATABASE '{0}' as {1};'''.format(
                db + '.DB', os.path.split(db)[-1]))

    def analyze(self):
        '''Analyze a database for better performance'''
        if self.engine == 'mysql':
            return
        else:
            self.execute('analyze;')

    #
    # query
    #
    def execute(self, *args, **kwargs):
        '''A shortcut to get cursor, execute and commit.'''
        cur = self.database.cursor()
        cur.execute(*args, **kwargs)
        return self.database.commit()

    def cursor(self):
        return self.database.cursor()

    def commit(self):
        return self.database.commit()

    #
    # Database
    #
    def hasDatabase(self, db):
        '''Test if a database exists in the current connection'''
        if self.engine == 'mysql':
            cur = self.database.cursor()
            cur.execute('SHOW DATABASES;')
            return db.lower() in [x[0].lower() for x in cur.fetchall()]
        else:
            return os.path.isfile(db if (db.endswith('.DB') or db.endswith('.proj')) else db + '.DB')

    def removeDatabase(self, db):
        if self.engine == 'mysql':
            cur = self.database.cursor()
            if self.hasDatabase(db):
                cur.execute('DROP DATABASE {};'.format(db))
            self.database.commit()
        else:
            # has to have file extension
            dbFile = db if (db.endswith('.proj') or db.endswith('.DB')) else db + '.DB'
            try:
                os.remove(dbFile)
            except:
                pass
            if os.path.isfile(dbFile):
                sys.exit('Failed to remove database {}'.format(db))

    #
    # Table
    #
    def tables(self, dbName=None):
        '''List all tables in a database'''
        cur = self.database.cursor()
        try:
            if self.engine == 'mysql':
                if dbName is None:
                    cur.execute("SHOW TABLES;")
                    return [x[0] for x in cur.fetchall()]
                else:
                    cur.execute("SHOW TABLES IN {};".format(dbName))
                    return [x[0] for x in cur.fetchall()]
            else:
                if dbName is None:
                    cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
                    return [x[0] for x in cur.fetchall() if not x[0].startswith('sqlite')]
                else:
                    cur.execute("SELECT name FROM {}.sqlite_master WHERE type='table';".format(dbName))
                    return [x[0] for x in cur.fetchall() if not x[0].startswith('sqlite')]
        except:
            return []

    def hasTable(self, table):
        '''Test if a table exists in the current database '''
        if '.' not in table:
            return table.lower() in [x.lower() for x in self.tables()]
        else:
            dbName, tableName = table.split('.')
            return tableName.lower() in [x.lower() for x in self.tables(dbName)]

    def dropIndex(self, index, table):
        if self.engine == 'mysql':
            self.execute('DROP INDEX {} ON {};'.format(index, table))
        else:
            self.execute('DROP INDEX {};'.format(index))

    def removeTable(self, table):
        '''Remove specified table'''
        cur = self.database.cursor()
        cur.execute('DROP TABLE {};'.format(table))
        # FIXME: should we automatically do VACUUM, this can be slow when the table is deletec
        # but can help performance for the creation of new tables.
        #if self.engine == 'sqlite3':
            # NOTE: It seems that re-generating a table can be VERY slow without vacuum.
        #    cur.execute('VACUUM;')
        self.database.commit()
    
    def renameTable(self, fromTable, toTable):
        '''Rename a table from fromTable to toTable'''
        cur = self.database.cursor()
        cur.execute('ALTER TABLE {} RENAME TO {};'.format(fromTable, toTable))
        self.database.commit()
        
    def backupTable(self, table):
        '''Backup a table to table_timestamp'''
        while True:
            new_table = '{}_{}'.format(table, time.strftime('%b%d_%H%M%S', time.gmtime()))
            if not self.hasTable(new_table):
                self.renameTable(table, new_table)
                return new_table
            time.sleep(1)
        

    def removeFields(self, table, cols):
        '''Remove fields from a table'''
        cur = self.database.cursor()
        if self.engine == 'mysql':
            cur.execute('ALTER TABLE {} {};'.format(table,
                ', '.join(['DROP COLUMN {}'.format(x) for x in cols])))
        else:
            # for my sqlite, we have to create a new table
            cur.execute('SELECT sql FROM sqlite_master WHERE name = "{}";'.format(table))
            schema = cur.fetchone()[0]
            fields = [x.strip() for x in schema.split(',')]
            fields[0] = fields[0].split('(')[1].strip()
            fields[-1] = (')'.join(fields[-1].split(')')[:-1])).strip()
            new_fields = [x for x in fields if x.split()[0].lower() not in [y.lower() for y in cols]]
            # rename existing table
            cur.execute('ALTER TABLE {0} RENAME TO _{0}_tmp_;'.format(table))
            # create a new table
            cur.execute('CREATE TABLE {} ('.format(table) + ',\n'.join(new_fields) + ');')
            # insert data back
            cur.execute('INSERT INTO {0} SELECT {1} FROM _{0}_tmp_;'.format(table, 
                ','.join([x.split()[0] for x in new_fields])))
            # remove old table
            cur.execute('DROP TABLE _{}_tmp_;'.format(table))

    def numOfRows(self, table):
        '''FIXME: We need a fast way to get number of rows...
        '''
        cur = self.database.cursor()
        # FIXME: How fast is this?
        cur.execute('SELECT count(*) FROM {};'.format(table))
        return cur.fetchone()[0]

    def startProgress(self, text):
        if self.engine == 'mysql':
            return
        self.prog = ProgressBar(text)
        self.database.set_progress_handler(self.prog.sqliteUpdate, self.batch)

    def stopProgress(self):
        if self.engine == 'mysql':
            return
        self.prog.done()
        self.database.set_progress_handler(None, self.batch)

    def getHeaders(self, table):
        '''Obtain field names of a table'''
        cur = self.database.cursor()
        try:
            if self.engine == 'mysql':
                cur.execute('SELECT column_name FROM information_schema.columns WHERE table_name={};'.format(self.PH),
                    table)
                return [x[0] for x in cur.fetchall()]
            else:
                cur.execute('SELECT * FROM {} LIMIT 0,0;'.format(table))
                return [x[0] for x in cur.description]
        except:
            return None

import token

def consolidateFieldName(proj, table, clause, alt_build=False):
    '''For input sift_score > 0.5, this function expand it to
    dbNSFP.sift_score > 0.5 and return a list of fields (dbNSFP.sift_score
    in this case). It also change pos to alt_pos if alt_build is true.
    We are using a Python tokenizer here so the result might be wrong.
    '''
    tokens = [x for x in tokenize.generate_tokens(cStringIO.StringIO(clause).readline)]
    res = []
    fields = []
    for i in range(len(tokens)):
        before_dot = (i + 1 != len(tokens)) and tokens[i+1][1] == '.'
        after_dot = i > 1 and tokens[i-1][1] == '.'
        #
        toktype, toval, _, _, _ = tokens[i]
        # replace chr by alt_chr if using an alternative reference genome.
        if alt_build and toval in ['chr', 'pos'] and not before_dot:
            toval = 'alt_' + toval
        #
        if toktype == token.NAME and toval.upper() not in SQL_KEYWORDS:
            if before_dot:
                # A.B, does not try to expand A
                res.append((toktype, toval))
            elif after_dot:
                # A.B, do not expand
                res.append((toktype, toval))
                # try to get fields:
                try:
                    for info in proj.linkFieldToTable('{}.{}'.format(tokens[i-2][1], toval), table):
                        fields.append(info.field)
                except ValueError as e:
                    proj.logger.debug(e)
            else:
                # A: try to expand A and identify fields
                try:
                    for info in proj.linkFieldToTable(toval, table):
                        fields.append(info.field)
                    # use expanded field, ONLY the last one should have the expanded fieldname
                    res.append((toktype, info.field))
                except ValueError as e:
                    proj.logger.debug(e)
                    res.append((toktype, toval))
        else:
            # fasttrack for symbols or function names
            res.append((toktype, toval))
    # a quick fix for a.b parsed to a .b. :-(
    return tokenize.untokenize(res).replace(' .', '.'), fields

#
# Utility function to calculate bins.
#
# This function implements a hashing scheme that UCSC uses (developed by Jim Kent) to 
# take in a genomic coordinate range and return a set of genomic "bins" that your range
# intersects.  I found a Java implementation on-line (I need to find the URL) and I
# simply manually converted the Java code into Python code.  
    
# IMPORTANT: Because this is UCSC code the start coordinates are 0-based and the end 
# coordinates are 1-based!!!!!!
        
# BINRANGE_MAXEND_512M = 512 * 1024 * 1024
# binOffsetOldToExtended = 4681; #  (4096 + 512 + 64 + 8 + 1 + 0)

_BINOFFSETS = (
    512+64+8+1,   # = 585, min val for level 0 bins (128kb binsize)    
    64+8+1,       # =  73, min val for level 1 bins (1Mb binsize) 
    8+1,          # =   9, min val for level 2 bins (8Mb binsize)  
    1,            # =   1, min val for level 3 bins (64Mb binsize)  
    0)            # =   0, only val for level 4 bin (512Mb binsize)
     
#    1:   0000 0000 0000 0001    1<<0       
#    8:   0000 0000 0000 1000    1<<3
#   64:   0000 0000 0100 0000    1<<6
#  512:   0000 0010 0000 0000    1<<9
 
_BINFIRSTSHIFT = 17;            # How much to shift to get to finest bin.
_BINNEXTSHIFT = 3;              # How much to shift to get to next larger bin.
_BINLEVELS = len(_BINOFFSETS)
  
#
# IMPORTANT: the start coordinate is 0-based and the end coordinate is 1-based.
#
def getUcscBins(start, end):
    bins = []
    startBin = start >> _BINFIRSTSHIFT
    endBin = (end-1) >> _BINFIRSTSHIFT
    for i in range(_BINLEVELS):
        offset = _BINOFFSETS[i];
        if startBin == endBin:
            bins.append(startBin + offset)
        else:
            for bin in range(startBin + offset, endBin + offset):
                bins.append(bin);
        startBin >>= _BINNEXTSHIFT
        endBin >>= _BINNEXTSHIFT
    return bins

def getMaxUcscBin(start, end):
    bin = 0
    startBin = start >> _BINFIRSTSHIFT
    endBin = (end-1) >> _BINFIRSTSHIFT
    for i in range(_BINLEVELS):
        offset = _BINOFFSETS[i];
        if startBin == endBin:
            if startBin + offset > bin:
                bin = startBin + offset
        else:
            for i in range(startBin + offset, endBin + offset):
                if i > bin:
                    bin = i 
        startBin >>= _BINNEXTSHIFT
        endBin >>= _BINNEXTSHIFT
    return bin
