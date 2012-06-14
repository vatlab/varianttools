#!/usr/bin/env python
#
# $File: tester.py $
# $LastChangedDate: 2012-06-05 12:31:19 -0500 (Tue, 05 Jun 2012) $
# $Rev: 1179 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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

import sys
import os
import math
import argparse
if sys.version_info.major == 2:
    import assoTests_py2 as t
else:
    import assoTests_py3 as t
from .project import Field

def freq(frequency):
    try:
        value = float(input)
        if not (value >= 0 and value <= 1):
            msg = "{0} is not valid input. Valid input should fall in range [0, 1]".format(frequency)
            raise ValueError(msg)
    except Exception as e:
        raise argparse.ArgumentTypeError(e)
    return value

#
# Statistical Association tests. The first one is a NullTest that provides
# some utility function and define an interface. All statistical tests should
# subclass from this class.
#
def getAllTests():
    '''List all tests (all classes that subclasses of NullTest/GLMBurdenTest) in this module'''
    return sorted([(name, obj) for name, obj in globals().iteritems() \
        if type(obj) == type(NullTest) and issubclass(obj, NullTest) \
            and name not in ('NullTest', 'GLMBurdenTest', 'CaseCtrlBurdenTest')], key=lambda x: x[0])


class NullTest:
    '''A base class that defines a common interface for association tests'''
    def __init__(self, logger=None, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to
        parse it'''
        self.logger = logger
        self.trait_type = None
        self.fields = []
        self.parseArgs(*method_args)
        #

    def parseArgs(self, method_args):
        # this function should never be called.
        raise SystemError('All association tests should define their own parseArgs function')

    def setData(self, data):
        self.data = data.clone()

    def calculate(self):
        '''Calculate and return p-values. It can be either a single value
        for all variants, or a list of p-values for each variant. Will print
        data if NullTest is called'''
        self.logger.info('No action is defined for NullTest')
        #self.logger.debug('Printing out group name, phenotypes, covariates and genotypes')
        #self.logger.debug(self.group)
        #self.logger.debug(self.data.phenotype())
        #self.logger.debug(self.data.covariates())
        #self.logger.debug(self.data.raw_genotype())
        return 0


class GroupStat(NullTest):
    '''Calculates basic statistics for each testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        #
        NullTest.__init__(self, logger, *method_args)
        # set fields according to parameter --stat
        if 'num_variants' in self.stat:
            self.fields.append(
                Field(name='num_variants', index=None, type='INT', adj=None, comment='number of variants in each group')
            )
        if 'sample_size' in self.stat:
            self.fields.append(
                Field(name='sample_size', index=None, type='INT', adj=None, comment='sample size')
            )

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Group statistics calculator. This 'test'
            is a statistics calculator based on the association-test framework. It is usually used
            with other method to produce statistics for each testing group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('--stat', choices=['num_variants', 'sample_size'], nargs='+',
            help='''Statistics to calculate, which can be number of variants in this group,
                total sample size.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def setData(self, data):
        # do not clone data because this test does not change data
        self.data = data

    def calculate(self):
        res = []
        for field in self.fields:
            if field.name == 'num_variants':
                # FIXME: this is inefficient because we do not have to create this object
                # to count the number of variants
                res.append(len(self.data.raw_genotype()[0]))
            elif field.name == 'sample_size':
                # FIXME: this is inefficient because we do not have to create this object
                # to count the number of variants
                res.append(self.data.samplecounts())
        return res

class GLMBurdenTest(NullTest):
    '''Generalized Linear regression test on collapsed genotypes within an association testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='beta_x', index=None, type='FLOAT', adj=None, comment='Test statistic. In the context of regression, this is estimate of effect size for x'), 
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        if self.permutations == 0:
            self.fields.append(Field(name='wald_x', index=None, type='FLOAT', adj=None, comment='Wald statistic for x (beta_x/SE(beta_x))'))
            for i in range(ncovariates):
                self.fields.extend([Field(name='beta_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='estimate of beta for covariate {}'.format(str(i+2))),
                                    Field(name='beta_{}_pvalue'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='p-value for covariate {}'.format(str(i+2))),
                                    Field(name='wald_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='Wald statistic for covariate {}'.format(str(i+2)))])
        else:
            self.fields.append(Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated'))
        #
        # NullTest.__init__ will call parseArgs to get the parameters we need
        self.regression_model = {'quantitative':0, 'disease':1}
        self.algorithm = self._determine_algorithm(ncovariates)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Generalized linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='GLMTest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--trait_type', type=str, choices = ['quantitative','disease'], default='quantitative',
            help='''Phenotype is quantitative trait or disease trait (0/1 coding).
            Default set to quantitative''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in GENE based analysis''')
        parser.add_argument('--weight', nargs='*', default=[],
        help='''Weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially. Additionally two special
            weights, i.e., 'MadsenBrowning' and 'MadsenBrownging_ctrl' are available if specified, which will apply a weighting theme
            based on observed allele frequencies from data. Note that all weights will be masked if --use_indicator is evoked
            ''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self, ncovariates):
        if ncovariates > 0:
            a_regression = t.MultipleRegression(self.permutations == 0, self.regression_model[self.trait_type])
        else:
            a_regression = t.SimpleLinearRegression() if self.trait_type == 'quantitative' else t.SimpleLogisticRegression()
        a_scoregene = t.BinToX() if self.use_indicator else t.SumToX()
        # data pre-processing
        algorithm = t.AssoAlgorithm([
            # code genotype matrix by MOI being 0 1 or 2
            t.CodeXByMOI(),
            # calculate sample MAF
            t.SetMaf(),
            # filter out variants having MAF > mafupper or MAF <= maflower
            t.SetSites(self.mafupper, self.maflower)
            ])
        # recode missing data
        if self.nan_adjust:
            algorithm.append(t.SetGMissingToMaf())
        # weight genotype codings by w(MAF)
        if 'MadsenBrowning' in self.weight and not self.use_indicator:
            algorithm.append(t.WeightByMaf("maf"))
        # weight with var_info/geno_info
        otherweights = [x for x in self.weight if not x.startswith('MadsenBrowning')]
        if len(otherweights) > 0:
            algorithm.append(t.WeightByInfo(otherweights))
        # association testing using analytic p-value
        if self.permutations == 0:
            algorithm.extend([
                # calculate genotype score for a set of variants
                a_scoregene,
                # fit regression model
                a_regression,
                # evaluate p-value for the Wald's statistic
                t.StudentPval(self.alternative)
                ])
        # association testing using permutation-based p-value
        else:
            if not self.variable_thresholds:
                a_permutationtest = t.FixedPermutator(
                        self.permute_by.upper(),
                        self.alternative,
                        self.permutations,
                        self.adaptive,
                        [a_regression]
                        )
                algorithm.extend([
                        a_scoregene,
                        a_permutationtest
                        ])
            else:
                a_permutationtest = t.VariablePermutator(
                        self.permute_by.upper(),
                        self.alternative,
                        self.permutations,
                        self.adaptive,
                        [a_scoregene, a_regression]
                        )
                algorithm.append(a_permutationtest)
        return algorithm

    def calculate(self):
        if self.trait_type == 'disease':
            self.data.countCaseCtrl()
        self.algorithm.apply(self.data)
        # get results
        pvalues = self.data.pvalue()
        regstats = self.data.statistic()
        regse = self.data.se()
        res = [self.data.samplecounts()]
        for (x, y, z) in zip(regstats, pvalues, regse):
            res.append(x)
            res.append(y)
            if self.permutations == 0:
                # Wald statistic
                res.append(x/z)
            else:
                # actual number of permutations
                if math.isnan(z): res.append(z)
                else: res.append(int(z))
        return res


#
# Derived statistical Association tests.
#

# quantitative traits
class LinRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for quantitative traits'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)
    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LinRegBurden',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in GENE based analysis''')
        parser.add_argument('--weight', nargs='*', default=[],
        help='''Weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially. Additionally two special
            weights, i.e., 'MadsenBrowning' and 'MadsenBrownging_ctrl' are available if specified, which will apply a weighting theme
            based on observed allele frequencies from data. Note that all weights will be masked if --use_indicator is evoked
            ''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.trait_type = 'quantitative'

class CollapseQt(GLMBurdenTest):
    '''Collapsing method for quantitative traits, Li & Leal 2008'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for quantitative traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = []
        self.trait_type = 'quantitative'

class BurdenQt(GLMBurdenTest):
    '''Burden test for quantitative traits, Morris & Zeggini 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for quantitative traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = []
        self.trait_type = 'quantitative'

class WeightedSumQt(GLMBurdenTest):
    '''Weighted sum statistic for quantitative traits, in the spirit of Madsen & Browning 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum statistic for quantitative traits (in the spirit of Madsen & Browning 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants will be weighted by 1/sqrt(P*(1-P)) and the weighted codings will be summed
            up as one regressor''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = ['MadsenBrowning']
        self.trait_type = 'quantitative'

class VariableThresholdsQt(GLMBurdenTest):
    '''Variable thresholds method for quantitative traits, in the spirit of Price et al 2010'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for quantitative traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.use_indicator=False
        self.weight = []
        self.trait_type = 'quantitative'

# binary traits
class LogitRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for binary traits'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Logistic regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='LogitRegBurden',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--use_indicator', action='store_true',
            help='''This option, if evoked, will apply binary coding to genotype groups
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in GENE based analysis''')
        parser.add_argument('--weight', nargs='*', default=[],
        help='''Weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially. Additionally two special
            weights, i.e., 'MadsenBrowning' and 'MadsenBrownging_ctrl' are available if specified, which will apply a weighting theme
            based on observed allele frequencies from data. Note that all weights will be masked if --use_indicator is evoked
            ''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.trait_type = 'disease'

class CollapseBt(GLMBurdenTest):
    '''Collapsing method for binary traits, Li & Leal 2008'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for binary traits (Li & Leal 2008).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants within a group will be collapsed into a single binary coding using an indicator function
            (coding will be "1" if ANY locus in the group has the alternative allele, "0" otherwise)''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='CBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = True
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = []
        self.trait_type = 'disease'

class BurdenBt(GLMBurdenTest):
    '''Burden test for binary traits, Morris & Zeggini 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for binary traits (Morris & Zeggini 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, the group of variants will be coded using the counts of variants within the group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='BBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = []
        self.trait_type = 'disease'

class WeightedSumBt(GLMBurdenTest):
    '''Weighted sum statistic for binary traits, in the spirit of Madsen & Browning 2009'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum statistic for binary traits (in the spirit of Madsen & Browning 2009).
            p-value is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, variants will be weighted by 1/sqrt(P*(1-P)) and the weighted codings will be summed
            up as one regressor''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.use_indicator = False
        self.maflower = 0.0
        self.permutations = 0
        self.variable_thresholds = False
        self.weight = 'MadsenBrowning_ctrl'
        self.trait_type = 'disease'

class VariableThresholdsBt(GLMBurdenTest):
    '''Variable thresholds method for binary traits, in the spirit of Price et al 2010'''
    def __init__(self, logger=None, *method_args):
        GLMBurdenTest.__init__(self, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for binary traits (in the spirit of Price et al 2010).
            The burden test statistic of a group of variants will be maximized over subsets of variants defined by applying different minor allele frequency
            thresholds. Significance of the statistic obtained is evaluated via permutation''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='VTBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--permute_by', metavar='XY', choices = ['X','Y','x','y'], default='Y',
            help='''Permute phenotypes ("Y") or genotypes ("X"). Default is "Y"''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.''')

        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.variable_thresholds = True
        self.weight = []
        self.use_indicator=False
        self.trait_type = 'disease'

###
# below are interface for some single covariate tests
# implemented as they were originally published
###

class CaseCtrlBurdenTest(NullTest):
    '''Single covariate case/ctrl burden test on collapsed genotypes within an association testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='statistic', index=None, type='FLOAT', adj=None, comment='Test statistic.'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        if self.permutations > 0:
            self.fields.append(Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated'))
        #
        # NullTest.__init__ will call parseArgs to get the parameters we need
        self.algorithm = self._determine_algorithm()

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Single covariate case/ctrl burden test including CMC, WSS, KBAC, RBT and aSum. p-value
            is calculated using exact/asymptotic distributions or permutation, depending on the input method. If --group_by
            option is specified, it will collapse the variants within a group into a single pseudo coding''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='SBurdenTest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--aggregation_theme', type=str, choices = ['CMC','WSS', 'KBAC', 'RBT', 'aSum'], default='CMC',
            help='''Choose from "CMC", "WSS", "KBAC", "RBT", "aSum".
            Default set to "CMC"''')
        parser.add_argument('--alternative', metavar='SIDED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To NOT using adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test. It is only applicatable to one sided test of CMC.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self):
        algorithm = t.AssoAlgorithm([
            # code genotype matrix by MOI being 0 1 or 2
            t.CodeXByMOI(),
            # calculate sample MAF
            t.SetMaf(),
            # filter out variants having MAF > mafupper or MAF <= maflower
            t.SetSites(self.mafupper, self.maflower)
            ])

        # Now write implementation of each association methods separately.
        # association testing using analytic p-value
        if self.permutations == 0:
            if self.aggregation_theme == 'CMC':
                algorithm.extend([t.BinToX(),
                    t.Fisher2X2(self.alternative, self.midp)])
            elif self.aggregation_theme == 'WSS':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        1000,
                        1,
                        [t.SetMafByCtrl("ctrlmaf", 0),
                            t.WeightByMaf("ctrlmaf"),
                            t.SumToX(),
                            t.MannWhitneyu(store=True)]
                        )
                algorithm.extend([
                    a_permutationtest,
                    t.MannWhitneyuPval(1)
                    ])
            else:
                raise ValueError('Please specify number of permutations for {0} test'.format(self.aggregation_theme))

        # association testing using permutation-based p-value
        else:
            if self.aggregation_theme == 'WSS':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        1000,
                        1,
                        [t.SetMafByCtrl("ctrlmaf", reverse=False),
                            t.WeightByMaf("ctrlmaf"),
                            t.SumToX(),
                            t.MannWhitneyu(store=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'KBAC':
                algorithm.append(t.FindGenotypePattern())
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.KBACtest(reverse=False, weightOnly=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'RBT':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.RBTtest(self.alternative, weightOnly=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'aSum':
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.AsumScore(),
                            t.SimpleLogisticRegression()]
                        )
                algorithm.append(a_permutationtest)
            else:
                raise ValueError('Invalid permutation test {0}'.format(self.aggregation_theme))
        return algorithm


    def calculate(self):
        self.data.countCaseCtrl()
        self.algorithm.apply(self.data)
        pvalues = self.data.pvalue()
        statistics = self.data.statistic()
        se = self.data.se()
        res = [self.data.samplecounts()]
        for (x, y, z) in zip(statistics, pvalues, se):
            res.append(x)
            res.append(y)
            if self.permutations > 0:
                # actual number of permutations
                if math.isnan(z): res.append(z)
                else: res.append(int(z))
        return res
