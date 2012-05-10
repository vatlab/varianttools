#!/usr/bin/env python

import os, sys, shlex, re
from variant_tools.association import t
import unittest
from subprocess import Popen, PIPE
from zipfile import ZipFile
from random import choice
from array import array
from testUtils import ProcessTestCase, outputOfCmd, output2list

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
        #
        self.data.setGenotype(self.x2)
        self.assertEqual([list(i) for i in self.x2], [list(i) for i in self.data.raw_genotype()])
        self.assertEqual([],list(self.data.genotype()))
        # test calculate MAF
        a = t.SetMaf()
        a.apply(self.data)
        maf = [sum(x)/6.0 for x in zip(*self.x2)]
        self.assertEqual(list(self.data.getArrayVar("maf")), maf)
        # test SumToX
        a = t.SumToX()
        a.apply(self.data)
        self.assertEqual(list(self.data.genotype()), [sum(i) for i in self.x2])
        # test BinToX
        a = t.BinToX()
        a.apply(self.data)
        self.assertEqual(list(self.data.genotype()), [1 if sum(i)>0 else 0 for i in self.x2])
        # test filter by MAF
        a = t.SetSites(0.3, 0.1)
        a.apply(self.data)
        self.assertEqual(list(self.data.getArrayVar("maf")), [x for x in maf if x < 0.3 and x > 0.1])
        #
        self.data.setGenotype(self.x2)
        actions = [t.SetMaf(), t.SetSites(0.3, 0.1), t.SumToX()]
        a = t.AssoAlgorithm(actions)
        a.apply(self.data)
        self.assertEqual(list(self.data.genotype()), [1.0,2.0,0.0])

    def testSimpleReg(self):
        self.data.setGenotype(self.x)
        self.data.setPhenotype(self.y)
        # test for setPhenotype
        self.assertEqual(list(self.data.phenotype()), list(self.y))
        self.assertEqual(self.data.getDoubleVar("ybar"), sum(self.y) / (1.0*len(self.y)))
        self.assertEqual(self.data.samplecounts(), len(self.y))
        #
        data = self.data.clone()
        # simple linear regression
        actions = [t.SumToX(), t.SimpleLinearRegression(), t.StudentPval(2)]
        a = t.AssoAlgorithm(actions)
        a.apply(data)
        # compare with lm(y~x) in R
        self.assertEqual(round(data.statistic()[0]/data.se()[0], 3), round(3.898, 3))
        self.assertEqual(data.getDoubleVar("xbar"), sum([i[0] for i in self.x]) / (1.0*len(self.x)))
        self.assertEqual(round(data.pvalue()[0], 5), round(0.00212, 5))

    def testMultiReg(self):
        self.data.setPhenotype(self.y, self.z)
        self.assertEqual(self.data.getDoubleVar("ybar"), sum(self.y) / (1.0*len(self.y)))
        self.assertEqual(self.data.getIntVar("ncovar"), 3)
        self.assertEqual([list(x) for x in self.data.covariates()[:-1]], [list(x) for x in self.z])
        self.assertEqual(list(self.data.covariates()[-1]), list(self.z[0]))
        data = self.data.clone()
        data.setGenotype(self.x)
        #
        # multiple linear regression
        actions = [t.SumToX(), t.MultipleRegression(True, 0), t.StudentPval(2)]
        a = t.AssoAlgorithm(actions)
        a.apply(data)
        #y = c(165,-22,228,166,277,-18,102,13,231,22,-54,66,94,36)
        #x = c(2,0,1,2,1,0,1,0,2,1,0,0,1,0)
        #z2 = c(92,221,-3,143,170,-86,236,61,32,132,1,95,-206,83)
        #z3 = c(98,204,145,73,223,93,-50,269,107,135,-4,-99,120,14)
        #z1 = c(56,104,114,22,124,179,25,15,107,170,176,83,261,-15)
        #summary(lm(y~x+z1+z2+z3))
        self.assertEqual(round(data.statistic()[0]/data.se()[0], 3), round(3.371, 3))
        self.assertEqual(round(data.pvalue()[0], 5), round(0.00824, 5))

    def testPermute(self):
        'Test permutator with fixed maf threshold'
        ybar = self.data.setPhenotype(self.y, self.z)
        self.data.setGenotype(self.x)
        a = t.SumToX()
        a.apply(self.data)
        # test if permutation works
        data = self.data.clone()
        p = t.FixedPermutator('Y', 1, 1, 0.5, [t.SimpleLinearRegression()])
        p.apply(data)
        self.assertEqual(round(data.statistic()[0], 5), 98.22222)
        self.assertNotEqual(list(self.data.phenotype()), list(data.phenotype()))
        self.assertEqual(list(self.data.genotype()), list(data.genotype()))
        data = self.data.clone()
        p = t.FixedPermutator('X', 1, 1, 0.5, [t.SimpleLinearRegression()])
        p.apply(data)
        #
        self.assertEqual(round(data.statistic()[0], 5), 98.22222)
        self.assertEqual(list(data.phenotype()), list(self.data.phenotype()))
        self.assertNotEqual(list(data.genotype()), list(self.data.genotype()))
        # test if adaptive permutation works
        #n = 100000
        #val = vector()
        #for(i in 1:n) {
        # val[i] = summary(lm(y~x))$coef[2,3]
        # y=sample(y)
        #}
        #length(which(pval>=val[1]))/n *2
        data = self.data.clone()
        p = t.FixedPermutator('Y', 2, 10000000, 0.002, [t.SimpleLinearRegression()])
        p.apply(data)
        self.assertEqual(round(data.statistic()[0], 5), 98.22222)
        #print data.pvalue()

    def testVariablePermute(self):
        'Test permutator with variable thresholds'
        ybar = self.data.setPhenotype(self.y, self.z)
        self.data.setGenotype(self.x2)
        # test if permutation with variable thresholds works
        data = self.data.clone()
        a = t.SetMaf()
        a.apply(data)
        p = t.VariablePermutator('Y', 1, 1, 0.5, [t.SumToX(), t.BinToX()])
        p.apply(data)
        # all maf: 0.333333 0 0.333333 0.333333 0.333333 0 0.166667 0 0.333333 0.166667 0 0 0.166667 0
        # sorted: 0 0 0 0 0 0 0.166667 0.166667 0.166667 0.333333 0.333333 0.333333 0.333333 0.333333
        # unique: 0 0.166667 0.333333
        # first vt: 0 0 0 0 0 0 1 0 0 1 0 0 1 0
        # 2nd vt: 1 0 1 1 1 0 1 0 1 1 0 0 1 0
        # test if it works on pre-trimmed data
        data = self.data.clone()
        a = t.AssoAlgorithm([t.SetMaf(), t.SetSites(1, 0.17)])
        a.apply(data)
        p = t.VariablePermutator('Y', 1, 1, 0.5, [t.SumToX(), t.BinToX()])
        p.apply(data)
        # all maf: 0.333333 0.333333 0.333333 0.333333 0.333333
        # sorted: 0.333333 0.333333 0.333333 0.333333 0.333333
        # unique: 0.333333
        # first (only) vt: 1 0 1 1 1 0 0 0 1 0 0 0 0 0
        #
        # test variable threshold permutation
        self.data.setGenotype(self.x)
        data = self.data.clone()
        a = t.AssoAlgorithm([t.SetMaf(), t.VariablePermutator('Y', 2,
            10000000, 0.002, [t.SumToX(), t.SimpleLinearRegression()])])
        a.apply(data)
        self.assertEqual(round(data.statistic()[0], 5), 98.22222)
###
###
###
class RandomAssocTest(ProcessTestCase):

    def setUp(self):
        z = ZipFile('proj/Rtest.zip')
        wdir = os.getcwd()
        z.extractall(wdir)
        z.close()
        self.sname_pattern = "CEU"
        self.proj_name = "unitest"
        self.table = "variant_ex"
        self.rfile = "cache/test_associate.R"

    def testResultRand(self):
        'Test association results with R'
        # check if required R software is available
        totest = True
        try:
            Popen(['Rscript','--help'], shell = False, stdout = PIPE, stderr = PIPE)
        except:
            sys.stderr.write("Nothing done for 'testResultRand': Rscript program is required for this test\n")
            return

        ###
        # prepare input
        ###
        # get genotype information from vtools tped
        geno_tmp = output2list('vtools export {0} --format tped --samples 1 --style\
        numeric'.format(self.table))
        #print(geno_tmp)
        geno_tmp = [[y for idx, y in enumerate(x.split('\t')) if idx != 2 and idx != 1] for x in geno_tmp]
        geno_tmp = [['_'.join(x[:2])] + x[2:] for x in geno_tmp]
        #print(geno_tmp)
        # add sample names && generate header line
        with open(self.proj_name + ".log", 'r') as f:
            logtxt = f.readlines()
        snames = [x.split()[-1] for x in logtxt if self.sname_pattern in x]
        #print(snames)
        geno_tmp.insert(0, ['sample_name'] + snames)
        #print(geno_tmp)
        tgeno_tmp = list(zip(*geno_tmp))
        #print(tgeno_tmp)
        ###
        # SNV analysis
        ###
        # By R
        with open('cache/phenotype.txt', 'r') as f:
            for i in range(len(tgeno_tmp)):
                tgeno_tmp[i] = f.readline().rstrip() + '\t' + '\t'.join(tgeno_tmp[i][1:])
        tgeno_tmp = '\n'.join(tgeno_tmp)
        #print(tgeno_tmp)
        p = Popen(shlex.split("Rscript {0} testSnv".format(self.rfile)), shell = False, stdout =
                PIPE, stderr = PIPE, stdin = PIPE)
        (cout, cerr) = p.communicate(tgeno_tmp.encode(sys.getdefaultencoding()))
        p = Popen(['python', 'cache/gw_round.py'], shell=False, stdout = PIPE,
                stderr = PIPE, stdin = PIPE)
        (cout, cerr) = p.communicate(cout)
        Rres = cout.decode(sys.getdefaultencoding()).rstrip().replace('\nX','\n')
        Rres = [x.replace('_', '\t') for x in Rres.split('\n')[1:] if not 'NA' in x]
        Rres.sort()
        Rres = '\n'.join(Rres)
        #print(Rres)
        # By vtools
        (cout, cerr) = Popen(shlex.split('vtools associate {0} BMI --covariate aff sex -m "LinRegBurden --alternative 2"'.format(self.table)),
                shell = False, stdout = PIPE, stderr = PIPE, stdin = PIPE).communicate()
        vtoolsres = ['\t'.join([y for idx, y in enumerate(x.split('\t')) if idx in [0,1,3,4,5]]) for x in cout.decode(sys.getdefaultencoding()).rstrip().split('\n')[1:] if not 'NA' in x]
        vtoolsres.sort()
        vtoolsres = '\n'.join(vtoolsres)
        #print(vtoolsres)
        self.assertEqual(vtoolsres==Rres, True)
	    ###
	    # association for groupby
	    ###
        # generate random grouping theme numbers
        nums = list(set([choice(range(100)) for x in range(4)]))
        # output allele information
        (cout, cerr) = Popen(shlex.split('vtools output {0} chr pos ref alt --header'.format(self.table)), shell = False, stdout = PIPE, stderr = PIPE, stdin =
                PIPE).communicate()
        for num in nums:
            # add random group
            (ncout, ncerr) = Popen(['python','cache/fakecols.py',
                '{0}'.format(num)], shell = False, stdout = PIPE, stderr = PIPE, stdin =
                PIPE).communicate(cout)
            grping = ncout.decode(sys.getdefaultencoding()).rstrip()
            with open('vtools.group.tmp', 'w') as f:
                f.write(grping)
            # add new group in vtools
            Popen(shlex.split('''vtools update {0} --from_file vtools.group.tmp --format cache/vtools_randcol.fmt --var_info grpby'''.format(self.table)),
                shell = False, stdout = PIPE, stderr = PIPE, stdin = PIPE)
            # vtools association by group
            (ncout, ncerr) = Popen(shlex.split('vtools associate {0} BMI --covariate aff sex -m "LinRegBurden --alternative 2" -g grpby'.format(self.table)),
                    shell = False, stdout = PIPE, stderr = PIPE, stdin = PIPE).communicate()
            vtoolsres = ['\t'.join([y for idx, y in enumerate(x.split('\t')) if idx in [0,2,3,4]]) for x in ncout.decode(sys.getdefaultencoding()).rstrip().split('\n')[1:] if not 'NA' in x]
            vtoolsres = [x for x in vtoolsres if not '\t0\t1\t0' in x]
            vtoolsres.sort()
            # permutation based approach for empirical p-values
            (ncout, ncerr) = Popen(shlex.split('vtools associate {0} BMI --covariate aff sex -m "LinRegBurden --alternative 2 -p 5000 --permute_by X --adaptive 0.00001" -g grpby'.format(self.table)),
                shell = False, stdout = PIPE, stderr = PIPE, stdin = PIPE).communicate()
            vtoolsres2 = ['\t'.join([y for idx, y in enumerate(x.split('\t')) if idx in [0,2]]) for x in ncout.decode(sys.getdefaultencoding()).rstrip().split('\n')[1:] if not ('NA' in x or x.endswith('\t0'))]
            vtoolsres2.sort()
            # vtools association by variable threshold method
            (ncout, ncerr) = Popen(shlex.split('vtools associate {0} BMI --covariate aff sex -m "VariableThresholdsQt --alternative 2 -p 5000 --permute_by X --adaptive 0.00001" -g grpby'.format(self.table)), shell = False, stdout = PIPE, stderr = PIPE, stdin = PIPE).communicate()
            vtoolsres3 = ['\t'.join([y for idx, y in enumerate(x.split('\t')) if idx in [0,2]]) for x in ncout.decode(sys.getdefaultencoding()).rstrip().split('\n')[1:] if not ('NA' in x or x.endswith('\t0'))]
            vtoolsres3.sort()
            #print(vtoolsres3)
            #compare t values
            vtoolsres1 = [x.split()[0]+'\t'+x.split()[1] for x in vtoolsres]
            self.assertEqual(vtoolsres1 == vtoolsres2, 1)
            res2 = [float(x.split()[1]) for x in vtoolsres2]
            res3 = [float(x.split()[1]) for x in vtoolsres3]
            #print([[x,y] for x, y in zip(res2, res3) if abs(x) > abs(y)])
            self.assertEqual(len([x for x, y in zip(res2, res3) if abs(x) > abs(y)]), 0)
            ###
            # group analysis by R
            ###
            grping = ['_'.join(x.split()[0:2]) + '\t' + x.split()[4] for x in grping.split('\n')]
            geno_grped = '\n'.join([x + '\t' + '\t'.join(y[1:]) for x, y in zip(grping, geno_tmp)])
            (ncout,ncerr) = Popen(shlex.split("Rscript {0} scoreRegion".format(self.rfile)), shell = False, stdout =
                PIPE, stderr = PIPE, stdin = PIPE).communicate(geno_grped.encode(sys.getdefaultencoding()))
            tgeno_score = ncout.decode(sys.getdefaultencoding()).rstrip().split('\n')
            #print(tgeno_score)
            with open('cache/phenotype.txt', 'r') as f:
                for i in range(len(tgeno_score)):
                    tgeno_score[i] = f.readline().rstrip() + '\t' + '\t'.join(tgeno_score[i].split('\t')[1:])
            (ncout,ncerr) = Popen(shlex.split("Rscript {0} testRegion".format(self.rfile)), shell = False, stdout =
                PIPE, stderr = PIPE, stdin = PIPE).communicate('\n'.join(tgeno_score).encode(sys.getdefaultencoding()))
            (ncout,ncerr) = Popen(['python', 'cache/gw_round.py'], shell=False, stdout = PIPE,
                stderr = PIPE, stdin = PIPE).communicate(ncout)
            Rres = [x for x in ncout.decode(sys.getdefaultencoding()).rstrip().split('\n')[1:] if not ('NA' in x or x.endswith("\t0\t1\t0"))]
            Rres = [x.replace('X','') for x in Rres]
            Rres.sort()
            print('\n'.join(vtoolsres))
            print('-------------')
            print('\n'.join(Rres))
            self.assertEqual('\n'.join(vtoolsres) == '\n'.join(Rres), 1)

if __name__ == '__main__':
    unittest.main()
