#!/usr/bin/env python

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, output2list
from zipfile import ZipFile
import argparse
from random import choice
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
        print(data)
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
        print(data)
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
        self.assertSucc('vtools associate variant filename --method LinRegBurden -h')
        self.assertSucc('vtools associate variant phen2 -m "LinRegBurden"')
        self.assertSucc('vtools associate variant phen2 -m "LinRegBurden" -g chr')
        # regression methods
        self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LinRegBurden --alternative 2"')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -q1 0.05 -q2 0.001" -g chr')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -q1 5 -q2 0.001" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden --use_indicator" -g chr')
        self.assertFail('vtools associate variant phen2 -m "LinRegBurden --alternative 8" -g chr')
        # permutation based tests
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --permute_by x" -g chr')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --permute_by M" -g chr')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 10000000 --adaptive 0.000001"')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 10000000 --adaptive 24"')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --variable_thresholds"')

    def testResult(self):
        'Test association results'
        zip = ZipFile('proj/assoproj.zip')
        dir = os.getcwd()
        zip.extractall(dir)
        for i in range(8):
            runCmd('vtools update variant --from_file output/assogrp{}.txt --format fmt/randcol.fmt --var_info grpby'.format(str(i+1)))
            vtoolsout = output2list('vtools associate variant phen2 --covariate phen1 phen3 phen4 -m "LinRegBurden --alternative 2" -g grpby')
            vtoolsout.sort()
            vtoolsout = ['\t'.join([j for jdx, j in enumerate(x.split()) if jdx in [0,2,3,4]]) for idx, x in enumerate(vtoolsout) if idx > 0 and 'NAN' not in x]
            self.assertOutput('cat output/assores{}.txt'.format(str(i+1)), '\n'.join(vtoolsout)+'\n')


    def testPyAction(self):
        'Test action pyaction'
        #
        zip = ZipFile('proj/assoproj.zip')
        dir = os.getcwd()
        zip.extractall(dir)
        self.assertSucc('vtools associate variant phen2 -m "test_associate.PyActionTester" -g chr')

    def testResultRand(self):
        'Test association results with R'
        totest = True
        try:
            runCmd('Rscript --help')
            runCmd('awk --help')
        except:
            totest = False
        if totest:
            nums = set([choice(range(100)) for x in range(4)])
            nums = [0] + list(nums)
            for item in nums:
                self.assertOutput('bash test_associate.sh variant_ex CEU proj/Rtest.zip test_associate.R {}'.format(str(item)), '')

if __name__ == '__main__':
    unittest.main()
