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
        value = float(frequency)
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
        '''run an association test'''
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
        parser = argparse.ArgumentParser(description='''Group statistics calculator, usually is
               used with other method to produce statistics for each association testing group.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('--stat', choices=['num_variants', 'sample_size'], nargs='+',
            help='''Statistics to calculate for the group (currently avaiable statistics are number of variants
                and total sample size).''')
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
                res.append(self.data.locicounts())
            elif field.name == 'sample_size':
                res.append(self.data.samplecounts())
        return res

#
# single covariate case/ctrl burden tests
# implemented as they were originally published
#

class CaseCtrlBurdenTest(NullTest):
    '''Single covariate case/ctrl burden test on aggregated genotypes within an association testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='statistic', index=None, type='FLOAT', adj=None, comment='Test statistic.'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        if ncovariates > 1:
            self.logger.warning("This association test cannot handle covariates. Input option '--covariates' will be ignored.")
        if self.permutations > 0:
            self.fields.append(Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated'))
        # specify the trait type for the AssociationManager to make sure the input phenotype is proper (binary coding)
        self.trait_type = 'disease'
        # no external weight in these tests (for now)
        self.extern_weight = []
        # NullTest.__init__ will call parseArgs to get the parameters we need
        self.algorithm = self._determine_algorithm()

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Single covariate case/ctrl burden test including CMC, WSS, VT, VT_Fisher, KBAC, RBT and aSum.
            p-value is calculated using exact/asymptotic distributions or permutation, depending on the input method. "--group_by"
            option has to be specified to define genetic loci the burden tests will be applied to''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='SingleGeneCaseCtrlBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--aggregation_theme', type=str, choices = ['CMC','WSS', 'KBAC', 'RBT', 'aSum', 'VT', 'VT_Fisher', 'RareCover', 'Calpha'], default='CMC',
            help='''Choose from "CMC", "WSS", "KBAC", "RBT", "aSum", "VT", "VT_Fisher".
            Default set to "CMC"''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test. It is only applicatable to one sided test of CMC and VT_Fisher.''')
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
                        [t.WeightedGenotypeTester(
                            self.alternative,
                            self.extern_weight,
                            [t.BrowningWeight(self.alternative),
                            t.MannWhitneyu(alternative=self.alternative, store=True)]
                            )]
                        )
                algorithm.extend([
                    a_permutationtest,
                    t.WSSPvalue(self.alternative)
                    ])
            elif self.aggregation_theme == 'Calpha':
                # this has to be a two-sided test
                # the analytical version of Calpha is not reliable
                # implemented here to reflect its original publication
                algorithm.extend([t.CalphaTest(),
                    t.GaussianPval(1)])
            else:
                raise ValueError('Please specify number of permutations for {0} test'.format(self.aggregation_theme))

        # association testing using permutation-based p-value
        else:
            if self.aggregation_theme == 'WSS':
                # the rank test version of WSS only supports one-sided test
                # this has to be a one-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.WeightedGenotypeTester(
                            1,
                            self.extern_weight,
                            [t.BrowningWeight(1),
                            t.MannWhitneyu()]
                            )]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'KBAC':
                algorithm.append(t.FillGMissing(method="mlg"))
                algorithm.append(t.FindGenotypePattern())
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.KBACtest(alternative=self.alternative, weightOnly=False)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'RBT':
                algorithm.append(t.FillGMissing(method="mlg"))
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
                        [t.AdaptiveRvSum(),
                            t.SimpleLogisticRegression()]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'VT_Fisher':
                algorithm.extend([
                    t.FindVariantPattern(),
                    t.VTFisher(self.adaptive, alternative=self.alternative, midp=self.midp)
                    ])
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.VTFisher(self.adaptive, alternative=self.alternative, midp=self.midp)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'VT':
                algorithm.append(t.FindVariantPattern())
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.VTTest(alternative=self.alternative)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'RareCover':
                # this has to be a two-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.RareCoverTest(difQ = 0.5)]
                        )
                algorithm.append(a_permutationtest)
            elif self.aggregation_theme == 'Calpha':
                algorithm.append(t.FillGMissing(method="mlg"))
                # this has to be a two-sided test
                a_permutationtest = t.FixedPermutator(
                        'Y',
                        1,
                        self.permutations,
                        self.adaptive,
                        [t.CalphaTest(permutation=True)]
                        )
                algorithm.append(a_permutationtest)
            else:
                raise ValueError('Invalid permutation test {0}'.format(self.aggregation_theme))
        return algorithm


    def calculate(self):
        if self.data.locicounts() <= 2:
            raise ValueError("Cannot apply burden test on input data (number of variant sites has to be at least 3).")
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


#
# regression framework for some burden tests
# modified extended tests on the basis of their original publications
# to be able to take advantage of regression analysis
# while maintaining the spirit of each association method
#

class GLMBurdenTest(NullTest):
    '''Generalized Linear regression test on aggregated genotypes within an association testing group'''
    def __init__(self, ncovariates, logger=None, *method_args):
        # NullTest.__init__ will call parseArgs to get the parameters we need
        NullTest.__init__(self, logger, *method_args)
        # set fields name for output database
        self.fields = [Field(name='sample_size', index=None, type='INT', adj=None, comment='Sample size'),
                        Field(name='beta_x', index=None, type='FLOAT', adj=None, comment='Test statistic. In the context of regression this is estimate of effect size for x'),
                        Field(name='pvalue', index=None, type='FLOAT', adj=None, comment='p-value')]
        self.ncovariates = ncovariates
        if self.permutations == 0:
            self.fields.append(Field(name='wald_x', index=None, type='FLOAT', adj=None, comment='Wald statistic for x (beta_x/SE(beta_x))'))
            for i in range(self.ncovariates):
                self.fields.extend([Field(name='beta_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='estimate of beta for covariate {}'.format(str(i+2))),
                                    Field(name='beta_{}_pvalue'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='p-value for covariate {}'.format(str(i+2))),
                                    Field(name='wald_{}'.format(str(i+2)), index=None, type='FLOAT', adj=None, comment='Wald statistic for covariate {}'.format(str(i+2)))])
        else:
            self.fields.append(Field(name='num_permutations', index=None, type='INTEGER', adj=None, comment='number of permutations at which p-value is evaluated'))
        # check for weighting theme
        if self.permutations == 0 and self.weight in ['Browning', 'KBAC', 'RBT']:
            raise ValueError("Weighting theme {0} requires the use of permutation tests. Please specify number of permutations".format(self.weight))
        self.regression_model = {'quantitative':0, 'disease':1}
        self.algorithm = self._determine_algorithm()

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Generalized linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. For quantitative traits weights will be based on
               pseudo case/ctrl status defined by comparison with the mean of the quantitative traits. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
            ''')

        parser.add_argument('--nan_adjust', action='store_true',
            help='''This option, if evoked, will replace missing genotype values with a score relative to sample allele frequencies. The association test will
            be adjusted to incorperate the information. This is an effective approach to control for type I error due to differential degrees of missing genotypes among samples.
            ''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

    def _determine_algorithm(self):
        # define aggregation method and regression method
        a_scoregene = t.BinToX() if self.use_indicator else t.SumToX()
        if self.ncovariates > 0:
            a_regression = t.MultipleRegression(self.permutations == 0, self.regression_model[self.trait_type])
        else:
            a_regression = t.SimpleLinearRegression() if self.trait_type == 'quantitative' else t.SimpleLogisticRegression()
        if self.use_indicator:
            self.logger.warning("Cannot use weights in loci indicator coding. Setting weights to None.")
            self.extern_weight = []
            self.weight = 'None'
        # weighting theme
        a_wtheme = None
        a_wtimes = self.alternative
        if self.weight == 'Browning_all':
            a_wtheme = t.BrowningWeight(0)
            a_wtimes = 1
        elif self.weight == 'Browning':
            a_wtheme = t.BrowningWeight(self.alternative)
        elif self.weight == 'KBAC':
            a_wtheme = t.KBACtest(alternative=self.alternative, weightOnly=True)
        elif self.weight == 'RBT':
            a_wtheme = t.RBTtest(alternative=self.alternative, weightOnly=True)
        elif self.weight == 'None':
            pass
        else:
            raise ValueError('Invalid weighting theme {0}'.format(self.weight))
        # data pre-processing
        algorithm = t.AssoAlgorithm([
            # code genotype matrix by MOI being 0 1 or 2
            t.CodeXByMOI(),
            # calculate sample MAF
            t.SetMaf(),
            # filter out variants having MAF > mafupper or MAF <= maflower
            t.SetSites(self.mafupper, self.maflower)
            ])
        # special actions for KBAC and RBT weighting themes
        if self.weight in ['KBAC', 'RBT']:
            if not self.nan_adjust:
                self.logger.warning("In order to use weighting theme {0}, missing genotypes will be \
                        replaced by the most likely genotype based on MAF".format(self.weight))
            algorithm.append(t.FillGMissing(method="mlg"))
            if self.weight == 'KBAC':
                algorithm.append(t.FindGenotypePattern())
        # recode missing data
        if self.nan_adjust:
            algorithm.append(t.FillGMissing(method="maf"))
        # prepare the weighted sum tester
        a_wtester = None
        if a_wtheme:
            a_wtester = t.WeightedGenotypeTester(
                        a_wtimes,
                        self.extern_weight,
                        [a_wtheme, a_regression]
                    )
        #
        # association testing using analytic p-value
        #
        if self.permutations == 0:
            if not a_wtester:
                # not using any weighting themes
                algorithm.extend([
                    # weight by var_info/geno_info
                    t.WeightByInfo(self.extern_weight),
                    # calculate genotype score for a set of variants
                    a_scoregene,
                    # fit regression model
                    a_regression,
                    # evaluate p-value for the Wald's statistic
                    t.StudentPval(self.alternative)
                    ])
            else:
                # using Browning_all as weight
                algorithm.extend([
                    a_wtester,
                    t.StudentPval(self.alternative)
                    ])
        #
        # association testing using permutation-based p-value
        #
        else:
            if not a_wtester:
                algorithm.append(t.WeightByInfo(self.extern_weight))
            if not self.variable_thresholds:
                if not a_wtester:
                    algorithm.append(a_scoregene)
                a_permutationtest = t.FixedPermutator(
                        self.permute_by.upper(),
                        # FIXME: logic is correct yet confusing here
                        # use self.alternative = 2 for Browning_all only
                        # because it resets a_wtimes to 1
                        1 if a_wtimes == 2 else self.alternative,
                        self.permutations,
                        self.adaptive,
                        [a_wtester if a_wtester else a_regression]
                        )
                algorithm.append(a_permutationtest)
            else:
                a_permutationtest = t.VariablePermutator(
                        self.permute_by.upper(),
                        1 if a_wtimes == 2 else self.alternative,
                        self.permutations,
                        self.adaptive,
                        [a_wtester] if a_wtester else [a_scoregene, a_regression]
                        )
                algorithm.append(a_permutationtest)
        return algorithm

    def calculate(self):
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
# Derived association tests
# Separating disease traits and quantitative traits
#

# case/ctrl single covariate tests
class CFisher(CaseCtrlBurdenTest):
    '''Fisher's exact test on collapsed variant loci, Li & Leal 2008'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.midp and self.alternative == 2:
            self.logger.warning("midp option will be ignored for two-tailed test")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Collapsing test for case-control data (CMC test, Li & Leal 2008).
            Different from the original publication which jointly test for common/rare variants using Hotelling's t^2 method,
            this version of CMC will binarize rare variants (default frequency set to 0.01) within a group defined by "--group_by" and calculate p-value via Fisher's exact test.
            A "mid-p" option is available for one-sided test to obtain a less conservative p-value estimate.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='CFisher',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.
                Default set to "CFisher"''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2").
            Default set to 1''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))
        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'CMC'
        self.permutations = 0


class WSSRankTest(CaseCtrlBurdenTest):
    '''Weighted sum method using rank test statistic, Madsen & Browning 2009'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)


    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted sum method using rank test statistic, Madsen & Browning 2009. p-value
            is based on the significance level of the Wilcoxon rank-sum test. Two methods are avaliable for evaluating p-value: a semi-asymptotic
            p-value based on normal distribution, or permutation based p-value. Variants will be weighted by 1/sqrt(nP*(1-P)) and the weighted codings
            will be summed up for rank test. Two-sided test is available for the asymptotic version, which will calculate two p-values based on weights
            from controls and cases respectively, and use the smaller of them with multiple testing adjustment. For two-sided permutation based p-value
            please refer to "vtools show test WeightedBurdenBt"''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WSSRankTest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        #
        # arguments that are used by this test
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2"). Note that two-sided test is only
            available for asymptotic version of the test. Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations. Set it to zero to use the asymptotic version. Default is zero''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'WSS'
        if self.permutations and self.alternative == 2:
            self.logger.warning("Two-sided test is only available for asymptotic version of the test. Setting permutations to zero ...")
            self.permutations = 0


class VTtest(CaseCtrlBurdenTest):
    '''VT statistic for disease traits, Price et al 2010'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for VTtest method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds test for disease traits, Price et al 2010. 
        The burden test statistic of a group of variants will be
        maximized over subsets of variants defined by applying different minor allele
        frequency thresholds. This implementation provides two different statistics: 
        the original VT statistics in Price et al 2010 (default) and an adaptive VT 
        statistic combining the CFisher method (via "--cfisher" option). p-value is estimated by permutation test. 
        The adaptive VT statistic will not generate uniformly distributed p-value. 
        For a more generalized version of VT test, type "vtools show test VariableThresholdsBt / 
        VariableThresholdsQt". ''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='VTtest',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=1.0,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 1.0''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
            help='''Alternative hypothesis is one-sided ("1") or two-sided ("2"). Two sided test is only
                valid with "--cfisher" option evoked (please use "VariableThresholdsBt" otherwise). 
                Default set to 1''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--cfisher', action='store_true',
            help='''This option, if evoked, will use an adaptive VT test via Fisher's exact statistic.
            For more details, please refer to the online documentation.''')
        parser.add_argument('--midp', action='store_true',
            help='''This option, if evoked, will use mid-p value correction for one-sided Fisher's exact test. 
            It is only applicatable to one sided test with "--cfisher" option.''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'VT'
        if self.cfisher:
            self.aggregation_theme = 'VT_Fisher'


class KBAC(CaseCtrlBurdenTest):
    '''Kernel Based Adaptive Clustering method, Liu & Leal 2010'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for KBAC method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Kernel Based Adaptive Clustering method, Liu & Leal 2010. 
            Genotype pattern frequencies, weighted by a hypergeometric density kernel function, is compared 
            for differences between cases and controls. p-value is calculated using permutation for consistent 
            estimate with different sample sizes (the approximation method of the original publication is not implemented).
            Two-sided KBAC test is implemented by calculating a second statistic with case/ctrl label swapped, and
            the larger of the two statistic is used as two-sided test statistic''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='KBAC',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'KBAC'


class RBT(CaseCtrlBurdenTest):
    '''Repication Based Test for protective and deleterious variants, Ionita-Laza et al 2011'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for RBT method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Repication Based Test for protective and deleterious variants, 
            Ionita-Laza et al 2011. Variant sites are scored based on -log transformation of probability of having 
            more than observed variants in cases/ctrls; the RBT statistic is defined as sum of the variant sites scores.
            One-sided RBT is implemented in addition to the two-sided statistic described in the RBT paper. 
            p-value is estimated via permutation test.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='RBT',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'RBT'


class aSum(CaseCtrlBurdenTest):
    '''Adaptive Sum score test for protective and deleterious variants, Han & Pan 2010'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for aSum method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Adaptive Sum score test for protective and deleterious variants, 
            Han & Pan 2010. In the first stage of the test, each variant site are evaluated for excess of minor alleles 
            in controls and genotype codings are flipped, and the second stage performs a burden test similar to BRV 
            (Morris & Zeggini 2009). This two-stage test is robust to a mixture of protective/risk variants 
            within one gene, yet is computationally intensive. aSum test is a two-tailed test.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='aSum',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'aSum'
        self.alternative = 2


class Calpha(CaseCtrlBurdenTest):
    '''c-alpha test for unusual distribution of variants between cases and controls, Neale et al 2011'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''c-alpha test for unusual distribution of variants between 
            cases and controls, Neale et al 2011. It tests for deviation of variance of minor allele counts in 
            cases/ctrls from its expection based on binomial distribution. The statistic is asymptotically normally 
            distributed. p-value can be evaluted using either permutation or asymptotic distribution as described 
            in Neale et al 2011, although it is recommanded to use permutation to estimate a reliable p-value. 
            Calpha test is a two-tailed test''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='Calpha',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'Calpha'
        self.alternative = 2


class RareCover(CaseCtrlBurdenTest):
    '''A "covering" method for detecting rare variants association, Bhatia et al 2010.'''
    def __init__(self, ncovariates, logger=None, *method_args):
        CaseCtrlBurdenTest.__init__(self, ncovariates, logger, *method_args)
        if self.permutations == 0:
            raise ValueError("Please specify number of permutations for RareCover method")

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''A "covering" method for detecting rare variants association,
            Bhatia et al 2010. The algorithm combines a disparate collection of rare variants and maximize the association
            signal over the collection using a heuristic adaptive approach, which can be computationally intensive. 
            Different from VT method, it does not require rare variants evaluated being adjacent in minor allele 
            frequency ranking. RareCover test is a two-tailed test.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        parser.add_argument('--name', default='RareCover',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        parser.add_argument('-q1', '--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('-q2', '--maflower', type=freq, default=0.0,
            help='''Minor allele frequency lower limit. All variants having sample MAF>m2
            will be included in analysis. Default set to 0.0''')
        # permutations arguments
        parser.add_argument('-p', '--permutations', metavar='N', type=int, default=0,
            help='''Number of permutations''')
        parser.add_argument('--adaptive', metavar='C', type=freq, default=0.1,
            help='''Adaptive permutation using Edwin Wilson 95 percent confidence interval for binomial distribution.
            The program will compute a p-value every 1000 permutations and compare the lower bound of the 95 percent CI
            of p-value against "C", and quit permutations with the p-value if it is larger than "C". It is recommended to
            specify a "C" that is slightly larger than the significance level for the study.
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        args = parser.parse_args(method_args)
        # incorporate args to this class
        self.__dict__.update(vars(args))

        #
        # We add the fixed parameter here ...
        #
        self.aggregation_theme = 'RareCover'
        self.alternative = 2


# quantitative traits in regression framework
class LinRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for quantitative traits'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)
    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Linear regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
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
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'

class BurdenQt(GLMBurdenTest):
    '''Burden test for quantitative traits, Morris & Zeggini 2009'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'


class WeightedBurdenQt(GLMBurdenTest):
    '''Weighted genotype burden tests for quantitative traits, using one or many arbitrary external weights as well as one of 4 internal weighting themes'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted genotype burden tests for quantitative traits,
        using one or many arbitrary external weights as well as one of 4 internal weighting themes.
        External weights (variant/genotype annotation field) are passed into the test by --var_info and --geno_info options.
        Internal weighting themes are one of "Browning_all", "Browning", "KBAC" or "RBT". p-value is based on linear regression analysis
        and permutation procedure has to be used for "Browning", "KBAC" or "RBT" weights.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WQt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
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
        self.use_indicator = False
        self.maflower = 0.0
        self.variable_thresholds = False
        self.trait_type = 'quantitative'

class VariableThresholdsQt(GLMBurdenTest):
    '''Variable thresholds method for quantitative traits, in the spirit of Price et al 2010'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
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
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'quantitative'

# disease traits
class LogitRegBurden(GLMBurdenTest):
    '''A versatile framework of association tests for disease traits'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Logistic regression test. p-value
            is based on the significance level of the regression coefficient for genotypes. If --group_by
            option is specified, it will collapse the variants within a group into a generic genotype score''',
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--variable_thresholds', action='store_true',
            help='''This option, if evoked, will apply variable thresholds method to the permutation routine in burden test on aggregated variant loci''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
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
    '''Collapsing method for disease traits, Li & Leal 2008'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold collapsing method for disease traits (Li & Leal 2008).
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'disease'


class BurdenBt(GLMBurdenTest):
    '''Burden test for disease traits, Morris & Zeggini 2009'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Fixed threshold burden test for disease traits (Morris & Zeggini 2009).
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
        self.extern_weight = []
        self.weight = 'None'
        self.trait_type = 'disease'

class WeightedBurdenBt(GLMBurdenTest):
    '''Weighted genotype burden tests for disease traits, using one or many arbitrary external weights as well as one of 4 internal weighting themes'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Weighted genotype burden tests for disease traits,
        using one or many arbitrary external weights as well as one of 4 internal weighting themes. 
        External weights (variant/genotype annotation field) are passed into the test by --var_info and --geno_info options. 
        Internal weighting themes are one of "Browning_all", "Browning", "KBAC" or "RBT". p-value is based on logistic regression analysis
        and permutation procedure has to be used for "Browning", "KBAC" or "RBT" weights.''',
            prog='vtools associate --method ' + self.__class__.__name__)
        # argument that is shared by all tests
        parser.add_argument('--name', default='WBt',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        # no argumant is added
        parser.add_argument('--mafupper', type=freq, default=0.01,
            help='''Minor allele frequency upper limit. All variants having sample MAF<=m1
            will be included in analysis. Default set to 0.01''')
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
        parser.add_argument('--extern_weight', nargs='*', default=[],
            help='''External weights that will be directly applied to genotype coding. Names of these weights should be in one of '--var_info'
            or '--geno_info'. If multiple weights are specified, they will be applied to genotypes sequencially.
            Note that all weights will be masked if --use_indicator is evoked.
            ''')
        parser.add_argument('--weight', type=str, choices = ['Browning_all', 'Browning', 'KBAC', 'RBT'], default = 'None',
            help='''Internal weighting themes inspired by various association methods. Valid choices are:
               'Browning_all', 'Browning', 'KBAC' and 'RBT'. Except for
               'Browning_all' weighting, tests using all other weighting themes has to calculate p-value via permutation.
               For details of the weighting themes, please refer to the online documentation.
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
        self.use_indicator = False
        self.maflower = 0.0
        self.variable_thresholds = False
        self.trait_type = 'disease'


class VariableThresholdsBt(GLMBurdenTest):
    '''Variable thresholds method for disease traits, in the spirit of Price et al 2010'''
    def __init__(self, ncovariates, logger=None, *method_args):
        GLMBurdenTest.__init__(self, ncovariates, logger, *method_args)

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''Variable thresholds in burden test for disease traits (in the spirit of Price et al 2010).
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
        parser.add_argument('--alternative', metavar='TAILED', type=int, choices = [1,2], default=1,
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
            To disable the adaptive procedure, set C=1. Default is C=0.1''')
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
        self.extern_weight = []
        self.weight = 'None'
        self.use_indicator=False
        self.trait_type = 'disease'
