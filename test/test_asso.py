#!/usr/bin/env python

import sys

def linearBurden():
  # set data
  data = gw.AssoData()
  data.setGenotype([map(float, x.split()[1:]) for x in file('txt/assoc.dat')])
  data.setPhenotype([float(x.split()[0]) for x in file('txt/assoc.dat')])
  #
  # preprocessing
  a = gw.SumToX()
  a.apply(data)
  #print(data.phenotype())
  #print(data.genotype())
  #print(data.raw_genotype())
  #
  # linear regression, no permutation
  a = gw.SimpleLinearRegression()
  a.apply(data)
  a = gw.GaussianPval(1)
  a.apply(data)
  print data.pvalue()
 
  # linear regression, permutation based 
  p = gw.PhenoPermutator(1000, [gw.SimpleLinearRegression()])
  print p.permute(data) / 1000.0


if __name__ == '__main__':
    
  try:
    import variant_tools.assoTests as gw
  except:
    sys.exit(5)
  #
  linearBurden()

