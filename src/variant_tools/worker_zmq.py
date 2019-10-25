#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
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

import glob
import json
import math
import os
import re
import sys
import threading
import time

import zmq

from variant_tools.accessor import Engine_Access
from variant_tools.assoTests import AssoData
from variant_tools.project import AnnoDBWriter, Field
from variant_tools.tester import *
from variant_tools.utils import (DatabaseEngine, PrettyPrinter, env,
                                 executeUntilSucceed)


class AssoTestsWorker:
    '''Association test calculator'''

    def __init__(self, param, grps, args, path, projName):
        param = json.loads(param)
        self.param = param
        self.table = param["table"]
        self.sample_IDs = param["sample_IDs"]
        self.phenotypes = param["phenotypes"]
        self.covariates = param["covariates"]
        self.phenotype_names = param["phenotype_names"]
        self.covariate_names = param["covariate_names"]
        self.var_info = param["var_info"]
        self.geno_info = param["geno_info"]
        self.tests = self.getAssoTests(args["methods"], len(args["covariates"]),
                                       args["unknown_args"])
        self.group_names = param["group_names"]
        self.missing_ind_ge = param["missing_ind_ge"]
        self.missing_var_ge = param["missing_var_ge"]
        self.sample_names = param["sample_names"]

        self.num_extern_tests = param["num_extern_tests"]
        self.grps = grps

        self.path = path
        self.param["tests"] = self.tests
        self.results = ResultRecorder(self.param, args["to_db"],
                                      args["delimiter"], args["force"])

        # self.shelves = {}
        #

        self.db = DatabaseEngine()
        self.db.connect(self.path + "/" + projName + '.proj', readonly=True)

        # self.db.attach(param.proj.name + '.proj', '__fromVariant', lock=self.shelf_lock)

        #
        self.g_na = float('NaN')
        if env.treat_missing_as_wildtype:
            self.g_na = 0.0
        # self.args=args

    # def __del__(self):
    #     # self.db.close()
    #     for val in list(self.shelves.values()):
    #         val.close()

    def filterGenotype(self, genotype, geno_info, var_info, gname):
        '''
        Filter genotypes for missing calls or lack of minor alleles. Not very efficient because
        it copies genotype, var_info and geno_info skipping the loci to be removed.
            - genotype is a Individual_list * Variants_list matrix of genotype values
            - var_info is a dictionary with each key being information corresponding Variant_list
            - geno_info is a dictionary with each key having a matrix of the same structure as genotype matrix
        '''
        # Step 1: filter individuals by genotype missingness at a locus

        missing_ratios = [
            sum(list(map(math.isnan, x))) / float(len(x)) for x in genotype
        ]
        which = [x < self.missing_ind_ge for x in missing_ratios]
        # check for non-triviality of phenotype data
        if sum(which) < 5:
            raise ValueError(
                "Sample size too small ({0}) to be analyzed for {1}.".format(
                    sum(which), repr(gname)))
        if len(which) - sum(which) > 0:
            env.logger.debug('In {}, {} out of {} samples will be removed due to '
                              'having more than {}% missing genotypes'.\
                              format(repr(gname), len(which) - sum(which), len(which),
                                     self.missing_ind_ge * 100))
        # Step 2: filter variants by genotype missingness at a locus
        keep_loci = []
        for i in range(len(genotype[0])):
            # tag individuals missing variant calls
            missingness_vi = list(
                map(math.isnan, [x[i] for x, y in zip(genotype, which) if y]))
            # unique genotype codings on the locus
            gt_codings = list(
                set([
                    x[i]
                    for x, y in zip(genotype, which)
                    if y and not math.isnan(x[i])
                ]))
            keep_loci.append(
                (float(sum(missingness_vi)) / float(len(missingness_vi))) <
                self.missing_var_ge and len(gt_codings) > 1)
        if len(keep_loci) - sum(keep_loci) > 0:
            for idx in range(len(genotype)):
                # filter genotype and geno_info
                genotype[idx] = [
                    i for i, j in zip(genotype[idx], keep_loci) if j
                ]
                for k in list(geno_info.keys()):
                    geno_info[k][idx] = [
                        i for i, j in zip(geno_info[k][idx], keep_loci) if j
                    ]
            # filter var_info
            for k in list(var_info.keys()):
                var_info[k] = [i for i, j in zip(var_info[k], keep_loci) if j]
            #
            env.logger.debug('In {}, {} out of {} loci will be removed due to '
                              'having no minor allele or having more than {}% missing genotypes'.\
                              format(repr(gname), len(keep_loci) - sum(keep_loci),
                                     len(keep_loci), self.missing_ind_ge * 100))
        # check for non-triviality of genotype matrix
        if len(genotype[0]) == 0:
            raise ValueError("No variant found in genotype data for {}.".format(
                repr(gname)))
            # raise ValueError("No variant found in genotype data for {}.".format(repr(gname)))
        return genotype, which, var_info, geno_info

    def getAssoTests(self, methods, ncovariates, common_args):
        '''Get a list of methods from parameter methods, passing method specific and common
        args to its constructor. This function sets self.tests as a list of statistical tests'''
        if not methods:
            raise ValueError(
                'Please specify at least one statistical test. '
                'Please use command "vtools show tests" for a list of tests')
        tests = []
        for m in methods:
            m = m.replace("\\", "")
            m = m.strip('"')
            common_args = [common_arg.strip('"') for common_arg in common_args]
            name = m.split()[0]
            args = m.split()[1:] + common_args
            try:
                if '.' in name:
                    # if the method is defined elsewhere
                    m_module, m_name = name.split('.', 1)
                    # also search current working directory
                    my_dir = os.getcwd()
                    env.logger.info('Loading {} from {}'.format(
                        m_module, my_dir))
                    if my_dir not in sys.path:
                        sys.path.append(my_dir)
                        # use the default level, which is -1 for python 2 and 0 for python 3
                        _temp = __import__(m_module, globals(), locals(),
                                           [m_name])
                        sys.path.pop()
                    else:
                        _temp = __import__(m_module, globals(), locals(),
                                           [m_name])
                    env.logger.info('Loading {}'.format(m_module))
                    method = getattr(_temp, m_name)(ncovariates, args)
                else:
                    method = eval(name)(ncovariates, args)
                # check if method is valid
                if not hasattr(method, 'fields'):
                    raise ValueError('Invalid association test method {}: '
                                     'missing attribute fields'.format(name))
                if not method.fields:
                    env.logger.warning(
                        'Association test {} has invalid or empty fields. '
                        'No result will be generated.'.format(name))
                tests.append(method)
            except Exception as e:
                raise ValueError(
                    'Failed to load association test {0}: {1}.'
                    'Please use command "vtools show tests" to list usable tests'
                    .format(name, e))
        return tests

    def setGenotype(self, which, data, info, grpname):
        geno = [x for idx, x in enumerate(data) if which[idx]]
        self.data.setGenotype(geno)
        self.data.setVar("gname", str(grpname))
        for field in list(info.keys()):
            self.data.setVar(
                '__geno_' + field,
                [x for idx, x in enumerate(info[field]) if which[idx]])

    def setPhenotype(self, which):
        '''Set phenotype data'''
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        # print(len(self.phenotypes[0]))
        # for idx, x in enumerate(self.phenotypes[0]):
        #     if which[idx]:
        #         print(idx,x)
        phen = [x for idx, x in enumerate(self.phenotypes[0]) if which[idx]]
        if self.covariates:
            covt = [[x
                     for idx, x in enumerate(y)
                     if which[idx]]
                    for y in self.covariates]
        if self.covariates:
            self.data.setPhenotype(phen, covt)
        else:
            self.data.setPhenotype(phen)

    def setVarInfo(self, data):
        for field in list(data.keys()):
            if field not in ['chr', 'pos']:
                self.data.setVar('__var_' + field, data[field])

    def getVarInfo(self, group, where_clause):
        var_info = {x: [] for x in self.var_info}
        query = 'SELECT variant_id {0} FROM __asso_tmp WHERE ({1})'.format(
            ',' + ','.join([x.replace('.', '_') for x in self.var_info])
            if self.var_info else '', where_clause)
        #env.logger.debug('Running query: {}'.format(query))
        cur = self.db.cursor()
        # SELECT can fail when the disk is slow which causes database lock problem.
        msg = 'Load variant info for group {} using association worker'.format(
            group)
        executeUntilSucceed(cur, query, 5, msg, group)
        #
        if not self.var_info:
            data = {x[0]: [] for x in cur.fetchall()}
        else:
            data = {x[0]: x[1:] for x in cur.fetchall()}
        variant_id = sorted(list(data.keys()), key=int)
        for idx, key in enumerate(self.var_info):
            if key not in ['variant.chr', 'variant.pos']:
                var_info[key] = [
                    data[x][idx] if
                    (type(data[x][idx]) in [int, float]) else float('NaN')
                    for x in variant_id
                ]
            else:
                var_info[key] = [data[x][idx] for x in variant_id]
        return var_info, variant_id

    def setPyData(self,
                  which,
                  geno,
                  var_info,
                  geno_info,
                  missing_code,
                  grpname,
                  recode_missing=True):
        '''set all data to a python dictionary'''

        def fstr(x):
            try:
                float(x)
            except:
                x = str(x)
            return x

        #
        if len(self.phenotypes) > 1:
            raise ValueError('Only a single phenotype is allowed at this point')
        #
        self.pydata['name'] = grpname
        #
        if recode_missing:
            self.pydata['genotype'] = [[
                missing_code if math.isnan(e) else e for e in x
            ] for idx, x in enumerate(geno) if which[idx]]
        else:
            self.pydata['genotype'] = [
                x for idx, x in enumerate(geno) if which[idx]
            ]
        #
        try:
            self.pydata['coordinate'] = [(str(x), str(y)) for x, y in zip(
                var_info['variant.chr'], var_info['variant.pos'])]
        except:
            self.pydata['coordinate'] = []
        # var_info
        self.pydata['var_info'] = []
        self.pydata['var_info_header'] = []
        for k, item in list(var_info.items()):
            if k != 'variant.chr' and k != 'variant.pos':
                self.pydata['var_info_header'].append(k)
                self.pydata['var_info'].append(list(map(fstr, item)))
        self.pydata['var_info'] = list(zip(*self.pydata['var_info']))
        # geno_info
        self.pydata['geno_info'] = []
        self.pydata['geno_info_header'] = []
        for k, item in list(geno_info.items()):
            self.pydata['geno_info_header'].append(k)
            if recode_missing:
                self.pydata['geno_info'].append(
                    [[missing_code if math.isnan(e) else e
                      for e in x]
                     for idx, x in enumerate(item)
                     if which[idx]])
            else:
                self.pydata['geno_info'].append(
                    [x for idx, x in enumerate(item) if which[idx]])
        # convert geno_info to 3 dimensions:
        # D1: samples
        # D2: variants
        # D3: geno_info
        self.pydata['geno_info'] = list(zip(*self.pydata['geno_info']))
        self.pydata['geno_info'] = [
            list(zip(*item)) for item in self.pydata['geno_info']
        ]
        unique_names = self.sample_names
        if len(self.sample_names) != len(set(self.sample_names)):
            env.logger.warning(
                "Duplicated sample names found. Using 'sample_ID.sample_name' as sample names"
            )
            unique_names = [
                "{0}.{1}".format(i, s)
                for i, s in zip(self.sample_IDs, self.sample_names)
            ]
        self.pydata['sample_name'] = [
            str(x) for idx, x in enumerate(unique_names) if which[idx]
        ]
        self.pydata['phenotype_name'] = self.phenotype_names
        self.pydata['phenotype'] = [
            x for idx, x in enumerate(self.phenotypes[0]) if which[idx]
        ]
        if self.covariates:
            self.pydata['covariate_name'] = self.covariate_names
            # skip the first covariate, a vector of '1''s
            self.pydata['covariates'] = [[
                x for idx, x in enumerate(y) if which[idx]
            ] for y in self.covariates[1:]]
        #
        if len(self.pydata['genotype']) == 0 or len(
                self.pydata['phenotype']) == 0 or len(
                    self.pydata['genotype'][0]) == 0:
            raise ValueError("No input data")
        if len(self.pydata['geno_info']) > 0 and len(
                self.pydata['genotype']) != len(self.pydata['geno_info']):
            raise ValueError("Genotype and genotype information do not match")

    def getChr(self, variantID, cur):
        find_chr = "SELECT chr from variant where variant_id={0}".format(
            variantID)
        chr = [rec[0] for rec in cur.execute(find_chr)]
        return chr[0]

    def getChrs(self, variantIDs, cur):
        idString = "(" + ",".join([str(variantID) for variantID in variantIDs
                                  ]) + ")"
        find_chr = "SELECT chr,variant_id from variant where variant_id in " + idString
        varDict = {}
        for rec in cur.execute(find_chr):
            if rec[0] not in varDict:
                varDict[rec[0]] = []
            varDict[rec[0]].append(rec[1])
        return varDict

    def transformGeneName(self, geneSymbol):
        if ("-" in geneSymbol):
            geneSymbol = geneSymbol.replace("-", "_")
        pattern = re.compile(r'\.')
        if pattern.findall(geneSymbol):
            geneSymbol = geneSymbol.replace(".", "_")
        return geneSymbol

    def getGenotype_HDF5(self, group, sample_IDs):
        """This function gets genotype of variants in specified group.
        """
        where_clause = ' AND '.join(
            ['{0}={1}'.format(x, self.db.PH) for x in self.group_names])
        cur = self.db.cursor()
        # variant info
        var_info, variant_ids = self.getVarInfo(group, where_clause)
        chr = self.getChr(variant_ids[0], cur)
        chrEnd = self.getChr(variant_ids[-1], cur)

        if chr != chrEnd:
            varDict = self.getChrs(variant_ids, cur)
            chrs = list(varDict.keys())
        else:
            varDict = {chr: variant_ids}
            chrs = [chr]

        # get genotypes / genotype info
        genotype = []
        geno_info = {x: [] for x in self.geno_info}
        # getting samples locally from my own connection

        geneSymbol = self.transformGeneName(group[0])
        files = glob.glob(self.path + "/tmp*_genotypes_multi_genes.h5")
        # files=glob.glob("tmp*_genotypes_multi_genes.h5")

        files = sorted(
            files, key=lambda name: int(name.split("/")[-1].split("_")[1]))

        for fileName in files:
            accessEngine = Engine_Access.choose_access_engine(fileName)

            if len(chrs) == 1:
                # colnames=accessEngine.get_colnames(chr,geneSymbol)
                # snpdict=accessEngine.get_geno_info_by_group(geneSymbol,chr)

                # for chr in chrs:

                colnames = accessEngine.get_colnames(chr)
                snpdict = accessEngine.get_geno_by_group(chr, geneSymbol)
                accessEngine.close()
                for ID in colnames:
                    data = snpdict[ID]

                    gtmp = [
                        data.get(x, [self.g_na] +
                                 [float('NaN')] * len(self.geno_info))
                        for x in varDict[chr]
                    ]
                    # handle -1 coding (double heterozygotes)
                    genotype.append(
                        [2.0 if x[0] == -1.0 else x[0] for x in gtmp])

                    #
                    # handle genotype_info
                    #
                    for idx, key in enumerate(self.geno_info):
                        geno_info[key].append([
                            x[idx + 1] if
                            (type(x[idx + 1]) in [int, float]) else float('NaN')
                            for x in gtmp
                        ])
            else:
                colnames = accessEngine.get_colnames(chrs[0])

                for ID in colnames:
                    gtmp = []
                    for chr in chrs:
                        snpdict = accessEngine.get_geno_by_group(
                            chr, geneSymbol)
                        data = snpdict[ID]

                        gtmp.extend([
                            data.get(x, [self.g_na] +
                                     [float('NaN')] * len(self.geno_info))
                            for x in varDict[chr]
                        ])
                        # handle -1 coding (double heterozygotes)
                    genotype.append(
                        [2.0 if x[0] == -1.0 else x[0] for x in gtmp])
                    #
                    # handle genotype_info
                    #
                    for idx, key in enumerate(self.geno_info):
                        geno_info[key].append([
                            x[idx + 1] if
                            (type(x[idx + 1]) in [int, float]) else float('NaN')
                            for x in gtmp
                        ])

                accessEngine.close()
        gname = ':'.join(list(map(str, group)))
        return self.filterGenotype(genotype, geno_info, var_info, gname)

    def run(self):

        grps = self.grps
        #
        valuePack = []
        for grp in grps:
            try:
                grpname = ":".join(list(map(str, grp)))
            except TypeError:
                grpname = "None"
            # if grp is None:
            #     break
            # env.logger.debug('Retrieved association unit {}'.format(repr(grpname)))
            #
            #
            self.data = AssoData()
            self.pydata = {}
            values = list(grp)

            try:
                genotype, which, var_info, geno_info = self.getGenotype_HDF5(
                    grp)
                # if I throw an exception here, the program completes in 5 minutes, indicating
                # the data collection part takes an insignificant part of the process.
                #
                # set C++ data object
                if (len(self.tests) - self.num_extern_tests) > 0:
                    self.setGenotype(which, genotype, geno_info, grpname)
                    self.setPhenotype(which)
                    self.setVarInfo(var_info)
                # set Python data object, for external tests
                if self.num_extern_tests:
                    self.setPyData(which, genotype, var_info, geno_info, None,
                                   grpname)
                # association tests
                for test in self.tests:
                    test.setData(self.data, self.pydata)
                    result = test.calculate(env.association_timeout)
                    # env.logger.debug('Finished association test on {}'.format(repr(grpname)))
                    values.extend(result)
                valuePack.append(values)
            # except KeyboardInterrupt as e:
            #     # die silently if stopped by Ctrl-C
            #     break
            except Exception as e:
                # env.logger.debug('An ERROR has occurred in process {} while processing {}: {}'.\
                #                   format(self.index, repr(grpname), e),exc_info=True)
                # self.data might have been messed up, create a new one
                print(e)
                self.data = AssoData()
                self.pydata = {}
                # return no result for any of the tests if an error message is captured.
                values.extend([
                    float('NaN')
                    for x in range(len(self.results.fields) - len(list(grp)))
                ])
                # print(values)
                valuePack.append(values)
        return valuePack


class ResultRecorder:

    def __init__(self,
                 params,
                 db_name=None,
                 delimiter=None,
                 update_existing=False):
        self.succ_count = 0
        self.failed_count = 0
        #
        self.group_names = params["group_names"]
        self.fields = []
        self.group_fields = []
        for n, t in zip(params["group_names"], params["group_types"]):
            self.group_fields.append(
                Field(
                    name=n, index=None, type=t, adj=None, fmt=None, comment=n))
        self.fields.extend(self.group_fields)

        for test in params["tests"]:
            if test.name:
                self.fields.extend([
                    Field(
                        name='{}_{}'.format(x.name, test.name),
                        index=None,
                        type=x.type,
                        adj=None,
                        fmt=None,
                        comment=x.comment) for x in test.fields
                ])
            else:
                self.fields.extend(test.fields)
        for field in self.fields:
            if '-' in field.name:
                raise ValueError('"-" is not allowed in field name {}'.format(
                    field.name))
        if len(self.fields) != len(set([x.name for x in self.fields])):
            raise ValueError(
                'Duplicate field names. Please rename one of the tests using parameter --name'
            )
        #
        self.printer = PrettyPrinter(delimiter=delimiter)
        self.printer.write([x.name for x in self.fields])
        self.writer = None
        if db_name:
            db_name = db_name if not db_name.lower().endswith(
                '.db') else db_name[:-3]
            old_pragma = env.sqlite_pragma
            # make sure each commit will write data to disk, the performance can be bad though.
            env.sqlite_pragma = 'synchronous=FULL,journal_mode=DELETE'
            self.writer = AnnoDBWriter(
                db_name,
                self.fields,
                'field',  # field annotation databases
                'Annotation database used to record results of association tests. Created on {}'
                .format(time.strftime('%a, %d %b %Y %H:%M:%S', time.gmtime())),
                '1.0',  # version 1.0
                {'*': self.group_names},  # link by group fields
                2,  # database format
                True,  # allow updating an existing database
                update_existing  # allow updating an existing field
            )
            # restore system sqlite_pragma
            env.sqlite_pragma = ','.join(old_pragma)
            #
            self.cur = self.writer.db.cursor()
            if self.writer.update_existing:
                #
                self.update_query = 'UPDATE {0} SET {1} WHERE {2};'.format(
                    db_name, ', '.join([
                        '{}={}'.format(x.name, self.writer.db.PH)
                        for x in self.fields[len(self.group_names):]
                    ]), ' AND '.join([
                        '{}={}'.format(x, self.writer.db.PH)
                        for x in self.group_names
                    ]))
                self.insert_query = 'INSERT INTO {0} ({1}) VALUES ({2});'.format(
                    db_name, ','.join([x.name for x in self.fields]),
                    ','.join([self.writer.db.PH] * len(self.fields)))
                self.select_query = 'SELECT {1} FROM {0};'.format(
                    db_name, ', '.join(self.group_names))
            else:
                self.insert_query = 'INSERT INTO {0} VALUES ({1});'.format(
                    db_name, ','.join([self.writer.db.PH] * len(self.fields)))
        #
        self.last_commit = time.time()

    def get_groups(self):
        '''Get groups that have been calculated'''
        self.cur.execute(self.select_query)
        return self.cur.fetchall()

    def record(self, valuePack):
        for res in valuePack:
            self.succ_count += 1
            if len([x for x in res if x != x
                   ]) == len(self.fields) - len(self.group_fields):
                # all fields are NaN: count this as a failure
                self.failed_count += 1
            else:
                self.printer.write([
                    '{0:G}'.format(x, precision=5)
                    if isinstance(x, float) else str(x) for x in res
                ])
            # also write to an annotation database?
            if self.writer:
                if self.writer.update_existing:
                    self.cur.execute(
                        self.update_query, res[len(self.group_names):] +
                        res[:len(self.group_names)])
                    # if no record to update, insert a new one
                    if self.cur.rowcount == 0:
                        self.cur.execute(self.insert_query, res)
                else:
                    # insert a new record
                    self.cur.execute(self.insert_query, res)
                # commit the records from time to time to write data to disk
        self.writer.db.commit()

    def completed(self):
        return self.succ_count

    def failed(self):
        return self.failed_count

    def done(self):
        self.printer.write_rest()
        if self.writer:
            self.writer.finalize()


def retryConnection(client, poll, context, retries_left, SERVER_ENDPOINT, msg):
    #print("W: No response from server, retrying…")
    # Socket is confused. Close and remove it.
    client.setsockopt(zmq.LINGER, 0)
    client.close()
    poll.unregister(client)
    retries_left -= 1
    if retries_left == 0:
        print("E: Server seems to be offline, abandoning")
        return [retries_left, client, poll]
    #print("I: Reconnecting and resending")
    # Create new connection
    client = context.socket(zmq.REQ)
    client.connect(SERVER_ENDPOINT)
    poll.register(client, zmq.POLLIN)
    client.send_json(msg)
    return [retries_left, client, poll]


def worker(REQUEST_TIMEOUT, REQUEST_RETRIES):
    # Setup ZMQ.
    context = zmq.Context()
    client = context.socket(zmq.REQ)
    pid = os.getpid()
    try:

        projectFolder = os.environ.get("PROJECTFOLDER")
        portFilePath = projectFolder + "/randomPort.txt"
        starttime = time.time()
        #wait for main program to start up
        while not os.path.exists(portFilePath) and time.time() - starttime < 20:
            time.sleep(1)
        selected_port = ""
        if os.path.isfile(portFilePath):
            with open(portFilePath, "r") as portFile:
                selected_port = portFile.read()
        else:
            raise ValueError("%s isn't a file!" % portFilePath)

        SERVER_ENDPOINT = "tcp://" + os.environ["ZEROMQIP"] + ":" + selected_port
        client.connect(SERVER_ENDPOINT)  # IP of master
        poll = zmq.Poller()
        poll.register(client, zmq.POLLIN)

        retries_left = REQUEST_RETRIES

        param = None
        grps = None
        args = None
        path = None
        projName = None
        result = ""
        work = {}

        while retries_left:

            # env.logger.info("send available")
            client.send_json({"msg": "available", "pid": pid})
            expect_reply = True
            while expect_reply:
                socks = dict(poll.poll(REQUEST_TIMEOUT))
                if socks.get(client) == zmq.POLLIN:
                    # Retrieve work and run the computation.
                    try:
                        work = client.recv_json(flags=zmq.NOBLOCK)
                        if work == {}:
                            continue
                        if "preprocessing" in work:
                            # print("preprocessing")
                            expect_reply = False
                            continue
                        if "noMoreWork" in work:
                            break
                        retries_left = REQUEST_RETRIES
                        expect_reply = False
                        param = work['param']
                        grps = work['grps']
                        grps_string = ','.join(
                            elems[0] for elems in work["grps"])
                        args = work['args']
                        path = work["path"]
                        projName = work["projName"]

                        # Running association
                        worker = AssoTestsWorker(param, grps, args, path,
                                                 projName)
                        result = worker.run()
                        result = json.dumps(result)

                    except zmq.error.Again:
                        pass

                    try:
                        if result != "":
                            client.send_json({
                                "msg": "result",
                                "result": result,
                                "pid": pid,
                                "grps": grps_string
                            })
                            expect_thanks = True
                            while expect_thanks:
                                socks = dict(poll.poll(REQUEST_TIMEOUT))
                                if socks.get(client) == zmq.POLLIN:
                                    client.recv()
                                    expect_thanks = False
                                else:
                                    msg = {
                                        "msg": "result",
                                        "result": result,
                                        "pid": pid,
                                        "grps": grps_string
                                    }
                                    retries_left, client, poll = retryConnection(
                                        client, poll, context, retries_left,
                                        SERVER_ENDPOINT, msg)
                                    if retries_left == 0:
                                        raise Exception(
                                            "Server connection lost")
                    except Exception as e:
                        print(e)
                        break

                else:
                    msg = {"msg": "available", "pid": pid}
                    retries_left, client, poll = retryConnection(
                        client, poll, context, retries_left, SERVER_ENDPOINT,
                        msg)
                    if retries_left == 0:
                        raise Exception("Server connection lost")
    except Exception as e:
        print(e)
    finally:
        client.close()
        context.term()


def worker_heartbeat():
    try:
        pid = os.getpid()
        context = zmq.Context()
        hb_socket = context.socket(zmq.REQ)
        poll = zmq.Poller()
        poll.register(hb_socket, zmq.POLLIN)

        projectFolder = os.environ.get("PROJECTFOLDER")
        portFilePath = projectFolder + "/randomPort_heartbeat.txt"
        starttime = time.time()
        while not os.path.exists(portFilePath) and time.time() - starttime < 20:
            time.sleep(1)
        selected_port = ""
        if os.path.isfile(portFilePath):
            with open(portFilePath, "r") as portFile:
                selected_port = portFile.read()
        else:
            raise ValueError("%s isn't a file!" % portFilePath)

        SERVER_ENDPOINT = "tcp://" + os.environ["ZEROMQIP"] + ":" + selected_port
        hb_socket.connect(SERVER_ENDPOINT)  # IP of master

        while True:
            hb_socket.send_json({"msg": "heartbeat", "pid": pid})
            sock = dict(poll.poll(2500))
            if sock.get(hb_socket) == zmq.POLLIN:
                reply = hb_socket.recv_json()
                if reply["msg"] == "heartbeat":
                    time.sleep(2)
                if reply["msg"] == "stop":
                    break
            else:
                # msg={ "msg": "heartbeat","pid":pid}
                # retries_left,client,poll=retryConnection(client,poll,context,retries_left,SERVER_ENDPOINT,msg)
                # if retries_left==0:
                break
    except ValueError as e:
        print(e)
    except Exception as e:
        print(e)
    finally:
        hb_socket.setsockopt(zmq.LINGER, 0)
        poll.unregister(hb_socket)
        hb_socket.close()
        context.term()


def main():
    try:
        REQUEST_TIMEOUT = 2500
        REQUEST_RETRIES = 3

        if os.environ.get("PROJECTFOLDER") is None:
            raise ValueError("Please set PROJECTFOLDER.")
        if os.environ.get("ZEROMQIP") is None:
            os.environ["ZEROMQIP"] = "127.0.0.1"
        time.sleep(5)

        thread = threading.Thread(target=worker_heartbeat)
        # thread.setDaemon(True)
        thread.start()
        worker(REQUEST_TIMEOUT, REQUEST_RETRIES)

        thread.join()

    except Exception as e:
        print(e)


if __name__ == "__main__":
    try:
        REQUEST_TIMEOUT = 2500
        REQUEST_RETRIES = 3

        if os.environ.get("PROJECTFOLDER") is None:
            raise ValueError("Please set PROJECTFOLDER.")
        if os.environ.get("ZEROMQIP") is None:
            os.environ["ZEROMQIP"] = "127.0.0.1"
        time.sleep(5)

        thread = threading.Thread(target=worker_heartbeat)
        # thread.setDaemon(True)
        thread.start()
        worker(REQUEST_TIMEOUT, REQUEST_RETRIES)

        thread.join()

    except Exception as e:
        print(e)
