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

import os
import time
from math import sqrt

from variant_tools.assoTests import gsl_cdf_ugaussian_P as pnorm
from variant_tools.assoTests import gsl_cdf_ugaussian_Pinv as qnorm
from variant_tools.project import AnnoDBWriter, Field
from variant_tools.utils import env, openFile


class MetaAnalysis:

    def __init__(self, files, beta, pval, se, size, linker, method):
        self.bcol = beta
        self.pcol = pval
        self.scol = size
        self.secol = se
        self.linker = linker
        self.method = method
        #
        self.data = [tuple(self.read(fn)) for fn in files]
        self.header = list(
            set([tuple([d[0][i] for i in linker]) for d in self.data]))
        if len(self.header) > 1:
            raise ValueError("Linker field mismatch in input data: {}".\
                                 format(' != '.join([','.join(x) for x in self.header])))
        self.link_by = list(self.header[0])
        self.header = list(self.header[0])
        self.header.extend(['beta_meta', 'pvalue_meta', 'sample_size_meta'])
        for idx, d in enumerate(self.data):
            self.header.extend([
                d[0][i] + '_{}'.format(idx + 1)
                for i in [self.bcol, self.pcol, self.scol]
            ])
        # find overlapping groups
        self.groups = list(
            set.intersection(
                *list(map(set, [list(d[1].keys()) for d in self.data]))))
        self.sample_size = {}

    def read(self, filename):
        try:
            f = openFile(filename)
            head = [
                x.strip()
                for x in f.readline().decode("utf-8").strip().split('\t')
            ]
            dat = {}
            for line in f.readlines():
                line = [
                    x.strip() for x in line.decode("utf-8").strip().split('\t')
                ]
                dat[tuple([str(line[i]) for i in self.linker])] = [
                    float(line[self.bcol]),
                    float(line[self.pcol]),
                    int(line[self.scol]),
                    float(line[self.secol]) if self.secol else None
                ]
            f.close()
        except Exception as e:
            raise ValueError('ERROR while loading data: {}'.format(e))
        return head, dat

    def calculate(self, grp):
        '''x=[beta_a, pval_a, size_a, se_a]; y=[beta_b, pval_b, size_b, se_b], etc'''
        beta = float('nan')
        z = 0
        bpne = [d[1][grp] for d in self.data]
        self.sample_size[grp] = sum([d[1][grp][2] for d in self.data])
        # skip nan p-value
        bpne = [x for x in bpne if x[1] == x[1]]
        if len(bpne) <= 1:
            # nothing to do; will not meta-analyze anything
            return beta, -9
        # meta analysis
        try:
            if self.method == "ssb":
                # sample size based method
                z = [
                    (abs(qnorm(x[1] / 2.0)) * (abs(x[0]) / x[0]) if x[0] else 1,
                     x[2]) for x in bpne
                ]
                z = sum([x[0] * sqrt(x[1]) for x in z]) / sqrt(
                    sum([x[1] for x in z]))
            else:
                # inverse variance based method
                # calculate std err from beta and p if se is not available
                #
                var = [
                    x[3]**2 if x[3] is not None else
                    (x[0] / qnorm(x[1] / 2.0))**2 for x in bpne
                ]
                w = [1 / x for x in var]
                sumw = sum(w)
                beta = sum([x[0] * y for x, y in zip(bpne, w)]) / sumw
                z = beta / sqrt(1 / sumw)
        except:
            # qnorm(..) is zero, i.e., not significant
            beta = float('nan')
            z = float('nan')
        return beta, min(pnorm(z), 1 - pnorm(z)) * 2.0

    def createDB(self, db_name):
        linker = [x for x in self.header if x in self.link_by]
        fields = [x for x in self.header if x not in self.link_by]
        db_name = db_name[:-3] if db_name.lower().endswith('.db') else db_name
        if os.path.exists(db_name + '.DB'):
            os.remove(db_name + '.DB')
        self.fields = []
        self.group_names = [x for x in linker]
        linker = [(x, 'INTEGER') if x.endswith('pos') else (x, 'VARCHAR(20)')
                  for x in linker]
        fields = [(x, 'FLOAT') if not x.startswith('sample_size') else
                  (x, 'INTEGER') for x in fields]
        for item in linker + fields:
            self.fields.append(
                Field(
                    name=item[0],
                    index=None,
                    type=item[1],
                    adj=None,
                    fmt=None,
                    comment=''))
        env.sqlite_pragma = 'synchronous=OFF,journal_mode=MEMORY'
        self.writer = AnnoDBWriter(
            db_name,
            self.fields,
            'field',  # field annotation databases
            'Combined association tests result database. Created on {}'.format(
                time.strftime('%a, %d %b %Y %H:%M:%S', time.gmtime())),
            '1.0',  # version 1.0
            {'*': self.group_names},  # link by group fields
            2,  # database format
            True,  # allow updating an existing database
            True  # allow updating an existing field
        )
        self.cur = self.writer.db.cursor()
        self.db_name = db_name
        self.insert_query = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(
            db_name, ','.join([x.name for x in self.fields]),
            ','.join([self.writer.db.PH] * len(self.fields)))

    def writeDB(self, res):
        self.cur.executemany(self.insert_query, res)
        return

    def done(self):
        self.writer.finalize()
