#!/usr/bin/env python

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase
from zipfile import ZipFile
import argparse
from random import choice
from variant_tools.association import NullTest
from variant_tools.assoTests import *



class ActionTester(NullTest):
    def __init__(self, ncovariates=0, logger=None, *method_args):
        NullTest.__init__(self, logger, *method_args)
        self.fields = []
        self.algorithm = AssoAlgorithm([])

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
            PyAction(func=self.func)
        )

    def func(self, data):
        print(data)
        return True

class SetMafTester(ActionTester):
    def __init__(self, ncovariates,logger=None, *method_args):
        ActionTester.__init__(self, ncovariates, logger, *method_args)
        self.algorithm.append(
            SetMaf(),
            PyAction(func=self.func)
        )

    def func(self, data):
        # here to test if maf is correctly set...
        print(data)
        return True

 
class TestAssociate(ProcessTestCase):
    def setUp(self):
        ProcessTestCase.setUp(self)
        self.runCmd('vtools admin --load_snapshot proj/assoproj.tar.gz')


    def testInterface(self):
        'Test association module interface'
        # basic commands
        self.assertSucc('vtools associate -h')
        self.assertSucc('vtools associate variant filename --method LinRegBurden -h')
        self.assertSucc('vtools associate variant phen2 -m "LinRegBurden"')
        self.assertSucc('vtools associate variant phen2 -m "LinRegBurden" -g chr --to_db utest --force')
        # regression methods
        self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LinRegBurden --alternative 2" --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden" -g chr --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -q1 0.05 -q2 0.001" -g chr --force')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -q1 5 -q2 0.001" -g chr --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden --use_indicator" -g chr --force')
        self.assertFail('vtools associate variant phen2 -m "LinRegBurden --alternative 8" -g chr --force')
        # permutation based tests
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100" -g chr --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --permute_by x" -g chr --force')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --permute_by M" -g chr --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 10000000 --adaptive 0.000001" --force')
        self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 10000000 --adaptive 24" --force')
        self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinRegBurden -p 100 --variable_thresholds" --force')

    def testBasic(self):
        'Test basic association results'
        for i in range(8):
            self.runCmd('vtools update variant --from_file output/assogrp{}.txt --format fmt/randcol.fmt --var_info grpby'.format(str(i+1)))
            vtoolsout = self.runCmd('vtools associate variant phen2 --covariate phen1 phen3 phen4 -m "LinRegBurden --alternative 2" -g grpby -f -j 8', ret='list')
            vtoolsout = ['\t'.join([j for jdx, j in enumerate(x.split()) if jdx in [0,4,5,6]]) for idx, x in enumerate(vtoolsout) if idx > 0 and 'NAN' not in x and "_LinRegBurden" not in x and "tmp" not in x]
            vtoolsout.sort()
            with open('output/assores'+str(i+1)+'.txt','r') as f:
                infile = f.readlines()
            self.assertEqual([x.rstrip() for x in infile], vtoolsout)

    # def testWeights(self):
    #     'Test for weighting theme'
    #     for i in range(8):
    #         self.runCmd('vtools update variant --from_file output/assogrp{}.txt --format fmt/randcol.fmt --var_info grpby'.format(str(i+1)))
    #         vtoolsout = self.runCmd('vtools associate variant phen2 --covariate phen1 phen3 phen4 -m "WeightedBurdenQt --alternative 2 --weight Browning_all" -g grpby', ret='list')
    #         vtoolsout = ['\t'.join([j for jdx, j in enumerate(x.split()) if jdx in [0,4,5,6]]) for idx, x in enumerate(vtoolsout) if idx > 0 and 'NAN' not in x]
    #         vtoolsout.sort()
    #         with open('output/assores_wss'+str(i+1)+'.txt','r') as f:
    #             infile = f.readlines()
    #         self.assertEqual([x.rstrip() for x in infile], vtoolsout)

    def testPyAction(self):
        'Test action pyaction'
        self.assertSucc('vtools associate variant phen2 -m "test_associate.PyActionTester" -g chr')

    # def testRext(self):
    #     'Test R extensions'
    #     self.assertSucc('vtools associate variant phen1 -m "SKAT --name SKAT disease -k IBS" -g chr')

if __name__ == '__main__':
    unittest.main()
