#!/usr/bin/env python

import sys

def linearBurden(test_data):
  ptime = 1000
  data = test_data.clone()
  actions = [t.SumToX(), t.SimpleLinearRegression(), t.GaussianPval(1)]
  a = t.ActionExecuter(actions)
  a.apply(data)
  p = t.PhenoPermutator(ptime, [t.SimpleLinearRegression()])
  print(data.pvalue())
  print(data.statistic())
  print((p.apply(data)+1.0) / (ptime+1.0))
  return 1

def logisticBurden(test_data):
  ptime = 1000
  data = test_data.clone()
  actions = [t.SumToX(), t.SimpleLogisticRegression(), t.GaussianPval(1)]
  a = t.ActionExecuter(actions)
  a.apply(data)
  p = t.PhenoPermutator(ptime, [t.SimpleLogisticRegression()])
  print(data.pvalue())
  print(data.statistic())
  print((p.apply(data)+1.0) / (ptime+1.0))
  return 1  
  
if __name__ == '__main__':
    
  try:
    import variant_tools.assoTests as t
  except:
    sys.exit(5)
  #
  # set data
  test_data = t.AssoData()
  test_data.setGenotype([map(float, x.split()[1:]) for x in file('txt/assoc.dat')])
  test_data.setPhenotype([float(x.split()[0]) for x in file('txt/assoc.dat')])
  test_data.setMaf()
  test_data.filterByMaf(upper=0.05, lower=0.0)
  test_data.mean_phenotype()
  test_data.count_cases()
  test_data.count_ctrls()
  #print(test_data.maf())
  #print(test_data.phenotype())
  #print(test_data.raw_genotype())
  #logisticBurden(test_data)
  #linearBurden(test_data)
  

