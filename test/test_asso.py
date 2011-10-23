#!/usr/bin/env python

import sys

def linearBurden(test_data):
  data = test_data
  # preprocessing
  a = gw.SumToX()
  a.apply(data)
  #print(data.genotype())
  #
  # linear regression, no permutation
  a = gw.SimpleLinearRegression()
  a.apply(data)
  print data.statistic()
  a = gw.GaussianPval(1)
  a.apply(data)
  print data.pvalue()
 
  # linear regression, permutation based 
  p = gw.PhenoPermutator(1000, [gw.SimpleLinearRegression()])
  print p.permute(data) / 1000.0

def logisticBurden(test_data):
  data = test_data
  # preprocessing
  a = gw.SumToX()
  a.apply(data)
  # logistic regression score test, no permutation
  a = gw.SimpleLogisticRegression()
  a.apply(data)
  print data.statistic()
  a = gw.GaussianPval(1)
  a.apply(data)
  print data.pvalue()
 
  # logistic regression, permutation based 
  p = gw.PhenoPermutator(1000, [gw.SimpleLogisticRegression()])
  print p.permute(data) / 1000.0
  
  
if __name__ == '__main__':
    
  try:
    import variant_tools.assoTests as gw
  except:
    sys.exit(5)
  #
  # set data
  test_data = gw.AssoData()
  test_data.setGenotype([map(float, x.split()[1:]) for x in file('txt/assoc.dat')])
  test_data.setPhenotype([float(x.split()[0]) for x in file('txt/assoc.dat')])
  #print(test_data.phenotype())
  #print(test_data.raw_genotype())
  linearBurden(test_data)
  logisticBurden(test_data)
  

