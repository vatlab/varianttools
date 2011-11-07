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

  def testAsso(self):
      'Test command associate'
      self.assertSucc('vtools associate -h')
      
if __name__ == '__main__':
    unittest.main()