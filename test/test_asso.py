#!/usr/bin/env python

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, outputOfCmd
from zipfile import ZipFile

class TestAsso(ProcessTestCase):
  def setUp(self):
      'Create a project'
      #runCmd('vtools init test -f')
      #runCmd('vtools import --format txt/assoc.fmt txt/assoc.dat --build hg18')
      #runCmd('vtools phenotype --from_file txt/assoc.phen)
      #too slow to build the test database. Better unzip a pre-stored project
      zip = ZipFile('txt/assoproj.zip')
      dir = os.getcwd()
      zip.extractall(dir)

  def removeProj(self):
      runCmd('vtools remove project')

  def testBasicAsso(self):
      'Test basic regression routines'
      # basic commands
      self.assertSucc('vtools associate -h')
      self.assertSucc('vtools associate variant filename --method LinearBurdenTest -h')
      self.assertSucc('vtools associate variant phen2 -m "LinearBurdenTest"')
      self.assertSucc('vtools associate variant phen2 -m "LinearBurdenTest" -g chr')
      # regression methods
      self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LinearBurdenTest" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LinearBurdenTest"')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -q1 0.05 -q2 0.001" -g chr')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -q1 5 -q2 0.001" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest --use_indicator" -g chr')
      self.assertSucc('vtools associate variant phen2 -m "LinearBurdenTest --alternative 2" -g chr')
      self.assertFail('vtools associate variant phen2 -m "LinearBurdenTest --alternative 8" -g chr')
      # permutation based tests
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 1000" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 1000 --permute_by x" -g chr')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 1000 --permute_by M" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 10000000 --adaptive 0.000001"')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 10000000 --adaptive 24"')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LinearBurdenTest -p 1000 --variable_thresholds"')
      
      
      
      
if __name__ == '__main__':
    unittest.main()
