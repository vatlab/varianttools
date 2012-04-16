#!/usr/bin/env python

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, outputOfCmd
from zipfile import ZipFile
import argparse

from variant_tools.association import NullTest, t

class ActionTester(NullTest):
    def __init__(self, ncovariates=0, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = []
        self.algorithm = t.AssoAlgorithm([])

    def parseArgs(self, method_args):
        parser = argparse.ArgumentParser(description='''test methods''')
        parser.add_argument('--name', default='',
            help='''Name of the test that will be appended to names of output fields, usually used to
                differentiate output of different tests, or the same test with different parameters.''')
        args = parser.parse_args(method_args)
        self.__dict__.update(vars(args))

    def calculate(self):
        self.algorithm.apply(self.data)
        return []

class PyActionTester(ActionTester):
    def __init__(self, ncovariates,logger=None, *method_args):
        ActionTester.__init__(self, ncovariates, logger, *method_args)
        self.algorithm.append(
            t.PyAction(func=self.func)
        )

    def func(self, data):
        print data
        return True

class SetMafTester(ActionTester):
    def __init__(self, ncovariates,logger=None, *method_args):
        ActionTester.__init__(self, ncovariates, logger, *method_args)
        self.algorithm.append(
            t.SetMaf(),
            t.PyAction(func=self.func)
        )

    def func(self, data):
        # here to test if maf is correctly set...
        print data
        return True


class TestAsso(ProcessTestCase):

    def removeProj(self):
        runCmd('vtools remove project')

    def testInterface(self):
        'Test association module interface'
        zip = ZipFile('proj/assoproj.zip')
        dir = os.getcwd()
        zip.extractall(dir)
        # basic commands
        self.assertSucc('vtools associate -h')
        self.assertSucc('vtools associate variant filename --method LNBT -h')
        self.assertSucc('vtools associate variant phen2 -m "LNBT"')
        self.assertSucc('vtools associate variant phen2 -m "LNBT" -g chr')
        # regression methods
        self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LNBT --alternative 2"')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -q1 0.05 -q2 0.001" -g chr')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -q1 5 -q2 0.001" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT --use_indicator" -g chr')
        self.assertFail('vtools associate variant phen2 -m "LNBT --alternative 8" -g chr')
        # permutation based tests
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 100" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 100 --permute_by x" -g chr')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 100 --permute_by M" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 10000000 --adaptive 0.000001"')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 10000000 --adaptive 24"')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 100 --variable_thresholds"')
        
    def testResult(self):
        'Test association results'
        zip = ZipFile('proj/assoproj2.zip')
        dir = os.getcwd()
        zip.extractall(dir)
        for i in range(3):
          runCmd('vtools update variant --from_file assodat/grp{}.txt --format assodat/addcol.fmt --var_info grpby'.format(str(i+1)))
          self.assertOutput('vtools associate variant_ex BMI --samples 1 --covariate aff sex -m "LNBT --alternative 2"', '', 0,'assodat/assores{}1.txt'.format(str(i+1)))
          self.assertOutput('vtools associate variant_ex BMI --samples 1 --covariate aff sex -m "LNBT --alternative 2" -g grpby', '', 0,'assodat/assores{}2.txt'.format(str(i+1)))

    def testPyAction(self):
        'Test action pyaction'
        # 
        zip = ZipFile('proj/assoproj.zip')
        dir = os.getcwd()
        zip.extractall(dir)
        self.assertSuc('vtools associate variant phen2 -m "test_associate.PyActionTester" -g chr')

if __name__ == '__main__':
    unittest.main()
