#!/usr/bin/env python

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, outputOfCmd
from zipfile import ZipFile

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
      self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LNBT" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 -m "LNBT"')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -q1 0.05 -q2 0.001" -g chr')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -q1 5 -q2 0.001" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT --use_indicator" -g chr')
      self.assertSucc('vtools associate variant phen2 -m "LNBT --alternative 2" -g chr')
      self.assertFail('vtools associate variant phen2 -m "LNBT --alternative 8" -g chr')
      # permutation based tests
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 1000" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 1000 --permute_by x" -g chr')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 1000 --permute_by M" -g chr')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 10000000 --adaptive 0.000001"')
      self.assertFail('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 10000000 --adaptive 24"')
      self.assertSucc('vtools associate variant phen2 --covariate phen1 phen3 -m "LNBT -p 1000 --variable_thresholds"')
      
  def testAssoSet2(self):
      'Test set2 for association module'
      zip = ZipFile('proj/assoproj2.zip')
      dir = os.getcwd()
      zip.extractall(dir)
      self.assertSucc('vtools associate variant BMI --samples "RACE=0" --covariate SEX -m "LNBT" -j8')

if __name__ == '__main__':
    unittest.main()
