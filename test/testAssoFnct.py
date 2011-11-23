#!/usr/bin/env python

import os
import glob
import variant_tools.assoTests as t
import unittest
import subprocess
from array import array
from testUtils import ProcessTestCase

class TestAssoFnct(ProcessTestCase):

  def setUp(self):
    
    self.z = [array('d', [1,1,1,1,1,1,1,1,1,1,1,1,1,1]), 
    array('d', [56,104,114,22,124,179,25,15,107,170,176,83,261,-15]),
    array('d', [92,221,-3,143,170,-86,236,61,32,132,1,95,-206,83]),
    array('d', [98,204,145,73,223,93,-50,269,107,135,-4,-99,120,14])]
    
    self.y = array('d', [165,-22,228,166,277,-18,102,13,231,22,-54,66,94,36])
    self.x = [array('d', [i]) for i in [2,0,1,2,1,0,1,0,2,1,0,0,1,0]]
    self.x2 = [array('d', [1,0,2,1,0,0,0,0,1,1,0,0,0,0]),
               array('d', [1,0,0,1,2,0,1,0,1,0,0,0,1,0]),
               array('d', [0,0,0,0,0,0,0,0,0,0,0,0,0,0])]
    self.data = t.AssoData()
    
  def testInitData(self):
    'Tests for basic data initialization and actions'
    self.data.setGenotype(self.x)
    self.assertEqual([list(i) for i in self.x], [list(i) for i in self.data.raw_genotype()])
    data = self.data.clone()
    self.assertEqual([list(i) for i in self.x], [list(i) for i in data.raw_genotype()])
    
    self.data.setGenotype(self.x2)
    self.assertEqual([list(i) for i in self.x2], [list(i) for i in self.data.raw_genotype()])
    self.assertEqual([],list(self.data.genotype()))
    # test calculate MAF
    a = t.SetMaf()
    a.apply(self.data)
    maf = [sum(x)/6.0 for x in zip(*self.x2)]
    self.assertEqual(list(self.data.maf()), maf)
    # test SumToX
    a = t.SumToX()
    a.apply(self.data)
    self.assertEqual(list(self.data.genotype()), [sum(i) for i in self.x2])
    # test BinToX
    a = t.BinToX()
    a.apply(self.data)
    self.assertEqual(list(self.data.genotype()), [1 if sum(i)>0 else 0 for i in self.x2])
    # test filter by MAF
    self.assertEqual(list(self.data.sites()), [])
    a = t.FilterX(0.3, 0.1)
    a.apply(self.data)
    self.assertEqual([list(i) for i in self.data.raw_genotype()], \
      [list([j for idx, j in enumerate(i) if maf[idx]<0.3 and maf[idx]>0.1]) for i in self.x2])
    self.assertEqual(list(self.data.sites()), \
                     [idx for idx, x in enumerate(maf) if x <0.3 and x>0.1])
    
  def testSimpleReg(self):
    self.data.setGenotype(self.x)
    ybar = self.data.setPhenotype(self.y)
    # test for setPhenotype
    self.assertEqual(list(self.data.phenotype()), list(self.y))
    self.assertEqual(ybar, sum(self.y) / (1.0*len(self.y)))
    self.assertEqual(self.data.samplecounts(), len(self.y))
    
    data = self.data.clone()
    # simple linear regression
    actions = [t.SumToX(), t.SimpleLinearRegression(), t.StudentPval(2)]
    a = t.ActionExecuter(actions)
    a.apply(data)
    # compare with lm(y~x) in R
    self.assertEqual(round(data.statistic(), 3), round(3.898, 3))
    self.assertEqual(data.meanOfX(), sum([i[0] for i in self.x]) / (1.0*len(self.x)))
    self.assertEqual(round(data.pvalue(), 5), round(0.00212, 5))
    
  def testMultiReg(self):
    ybar = self.data.setPhenotype(self.y, self.z)
    self.assertEqual(ybar, sum(self.y) / (1.0*len(self.y)))
    self.assertEqual(self.data.covarcounts(), 3)
    self.assertEqual([list(x) for x in self.data.covariates()[:-1]], [list(x) for x in self.z])
    self.assertEqual(list(self.data.covariates()[-1]), list(self.z[0]))
    self.data.setGenotype(self.x)
    
    data = self.data.clone()
    
    # multiple linear regression
    actions = [t.SumToX(), t.MultipleLinearRegression(), t.StudentPval(2)]
    a = t.ActionExecuter(actions)
    a.apply(data)
    #y = c(165,-22,228,166,277,-18,102,13,231,22,-54,66,94,36)
    #x = c(2,0,1,2,1,0,1,0,2,1,0,0,1,0)
    #z2 = c(92,221,-3,143,170,-86,236,61,32,132,1,95,-206,83)
    #z3 = c(98,204,145,73,223,93,-50,269,107,135,-4,-99,120,14)
    #z1 = c(56,104,114,22,124,179,25,15,107,170,176,83,261,-15)
    #summary(lm(y~x+z1+z2+z3))
    self.assertEqual(round(data.statistic(), 3), round(3.371, 3)) 
    self.assertEqual(round(data.pvalue(), 5), round(0.00824, 5)) 
      
if __name__ == '__main__':
  unittest.main()
