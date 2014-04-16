#!/usr/bin/env python
#
# Author: Bo Peng
#
# Purpose: Simulation of samples with rare variants.
#
'''
Simulating a population of sequences forward in time, subject to mutation,
natural selection and population expansion. Because most mutants are introduced
during the rapid population expansion, most of the alleles will be rare at the
end of the simulation. Samples simulated using this script can be used to study
genetic diseases caused by a large number of rare variants.
'''

#
# This script simulates mutants in mutational space. That is to say, mutants
# are stored as locations on chromosomes, instead of allele number. Because
# indiviudals have different number of mutants, zeros are used to fill the
# rest of the markers.
#

import simuOpt
simuOpt.setOptions(alleleType='mutant', optimized=False, quiet=False, version='1.0.5')

import simuPOP as sim
from simuPOP.utils import ProgressBar, migrIslandRates
#from simuPOP.sandbox import RevertFixedSites, revertFixedSites, MutSpaceSelector, MutSpaceMutator, MutSpaceRecombinator

from .utils import env

import os, sys, logging, math, time



class RandomFitness:
    def __init__(self):
        self.coefMap = {}

    def _newS(self, loc, alleles):
        raise ValueError("This function should be redefined in the subclass of this class")

    def __call__(self, loc, alleles):
        # because s is assigned for each locus, we need to make sure the
        # same s is used for fitness of genotypes 01 (1-s) and 11 (1-2s)
        # at each locus
        if loc in self.coefMap:
            s, h = self.coefMap[loc]
        else:
            res = self._newS(loc, alleles)
            if type(res) in [list, tuple]:
                if len(res) != 2:
                    raise ValueError("A list or tuple of s, h is expected")
                s, h = res
            else:
                s = res
                h = self.h
            #
            self.coefMap[loc] = s, h
            #env.logger.info('SEL: loc:{} allele: {} s:{} h:{}'.format(loc, alleles, s, h))
        if 0 in alleles:
            return 1. - s * h
        else:
            return 1. - s        


class ConstantFitness(RandomFitness):
    def __init__(self, s, h=0.5):
        RandomFitness.__init__(self)
        self.s = s
        self.h = h

    def _newS(self, loc, alleles):
        # each mutant gets a penalty of s
        return self.s, self.h

class CustomizedFitness(RandomFitness):
    def __init__(self, func):
        RandomFitness.__init__(self)
        self.func = func
        self.h = 0.5   # default value that can be overwritten

    def _newS(self, loc, alleles):
        # each mutant gets a penalty of s
        return self.func(loc, alleles)

class GammaDistributedFitness(RandomFitness):
    def __init__(self, k, theta, h=0.5):
        RandomFitness.__init__(self)
        self.k = k
        self.theta = theta
        self.h = h

    def _newS(self, loc, alleles):
        return sim.getRNG().randGamma(self.k, self.theta)

class BoundedMixedGammaDistributedFitness(RandomFitness):
    def __init__(self, k, theta, h=0.5, p=0, a=0, lower=-1000, upper=10000):
        RandomFitness.__init__(self)
        self.k = k
        self.theta = theta
        self.h = h
        self.p = p
        self.a = a
        self.lower = lower
        self.upper = upper
     
    def _newS(self, loc, alleles):
        if p > 0 and sim.getRNG().randUniform() < p:
            return a
        while True:
            s = sim.getRNG().randGamma(self.k, self.theta)
            if s >= self.lower and s <= self.upper:
                return s
        
class TrimmedMixedGammaDistributedFitness(RandomFitness):
    def __init__(self, k, theta, h=0.5, p=0, a=0, lower=-1000, upper=10000):
        RandomFitness.__init__(self)
        self.k = k
        self.theta = theta
        self.h = h
        self.p = p
        self.a = a
        self.lower = lower
        self.upper = upper
     
    def _newS(self, loc, alleles):
        if p > 0 and sim.getRNG().randUniform() < p:
            return a
        s = sim.getRNG().randGamma(self.k, self.theta)
        if s >= self.lower and s <= self.upper:
            return s
        elif s < self.lower:
            return 0
        else:
            return 1

def getSelector(selDist, selCoef, selModel='exponential'):
    # set default parameter
    if selCoef is None:
        # set default parameters
        if selDist == 'mixed_gamma':
            selCoef = (0.0186, 0.0001, 0.184, 0.160*2, 0.5, 0.0001, 0.1)
        elif selDist == 'mixed_gamma1':
            selCoef = (0, -1, 0.562341, 0.01, 0.5, 0.00001, 0.1)
        elif selDist == 'gamma1':
            selCoef = (0.23, 0.185*2, 0.5)
        elif selDist == 'gamma2':
            selCoef = (0.184, 0.160*2, 0.5)
        elif selDist == 'gamma3':
            selCoef = (0.206, 0.146*2, 0.5)
        elif selDist == 'constant':
            selCoef = [0.01, 0.5]
        elif not callable(selDist):
            raise ValueError("Unsupported random distribution")
    #
    # use a right selection operator.
    mode = {'multiplicative': sim.MULTIPLICATIVE,
        'additive': sim.ADDITIVE,
        'exponential': sim.EXPONENTIAL}[selModel]
    #
    #
    if callable(selDist):
        pyFunc = CustomizedFitness(selDist)
    elif selDist == 'mixed_gamma':
        pyFunc = BoundedMixedGammaDistributedFitness(k=selCoef[2], theta=selCoef[3],
            h=selCoef[4], p=selCoef[0], a=selfCoef[1], lower=selCoef[5], upper=selCoef[6])
    elif selDist == 'mixed_gamma1':
        pyFunc = TrimmedMixedGammaDistributedFitness(k=selCoef[2], theta=selCoef[3],
            h=selCoef[4], p=selCoef[0], a=selCoef[1], lower=selCoef[5], upper=selCoef[6])
    elif selDist.startswith('gamma'):
        pyFunc = GammaDistributedFitness(k=selCoef[0], theta=selCoef[1])
    elif selDist == 'constant':
        pyFunc = ConstantFitness(s=selCoef[0], h=selCoef[1])
    #
    return sim.PyMlSelector(func=pyFunc, mode=mode)



def simuRareVariants2(pop, refGenome, demoModel, mu, selector, recRate=0):
    #
    # Evolve
    env.logger.info('Add info fields fitness and migrate_to')
    pop.addInfoFields(['fitness', 'migrate_to'])
    startTime = time.clock()
    #
    base = {'A': [], 'C': [], 'G': [], 'T': [], 'N': []}
    for chr in range(pop.numChrom()):
        chr_name = pop.chromName(chr)
        for loc in range(pop.chromBegin(chr), pop.chromEnd(chr)):
            ref = refGenome.getBase(chr_name, int(pop.locusPos(loc)))
            base[ref].append(loc)
    env.logger.info('Simulated regions with {:,} basepair have {:,} A, '
        '{:,} C, {:,} G, {:,} T, and {:,} N on reference genome.'
        .format(pop.totNumLoci(), len(base['A']), len(base['C']), 
        len(base['G']), len(base['T']), len(base['N'])))
    #
    env.logger.info('Start evolving...')
    pop.dvars().last_time = time.time()
    exec('import time', pop.vars(), pop.vars())
    pop.evolve(
        initOps=sim.InitSex(),
        preOps=[
            sim.PyOutput('''Statistics outputted are
1. Generation number,
2. population size (a list),
3. number of segregation sites,
4. average number of mutants per individual
5. average allele frequency * 100
6. average fitness value
7. minimal fitness value of the parental population
''', at = 0),
            # revert alleles at fixed loci to wildtype
            sim.RevertFixedSites(),
            # 
            # 'A' is zero, no need to map in and out
            sim.AcgtMutator(model='K80', rate=[mu, 2], loci=base['A']),
            sim.AcgtMutator(model='K80', rate=[mu, 2], loci=base['C'],
                # base is C=0, 0,1,2,3 to C,G,T,A (1,2,3,0)
                mapIn=[1,2,3,0], mapOut=[3,0,1,2]),
            sim.AcgtMutator(model='K80', rate=[mu, 2], loci=base['G'],
                # base is G=0, 0,1,2,3 to G,T,A,C (2,3,0,1)
                mapIn=[2,3,0,1], mapOut=[2,3,0,1]),
            sim.AcgtMutator(model='K80', rate=[mu, 2], loci=base['T'],
                # base is T=0, 0,1,2,3 to T,A,C,G (3,0,1,2)
                mapIn=[3,0,1,2], mapOut=[1,2,3,0]),
            # selection on all loci
            selector,
            sim.IfElse('time.time() - last_time > 30', [
                sim.PyExec('last_time = time.time()'),
                # output statistics in verbose mode
                sim.Stat(popSize=True, meanOfInfo='fitness', minOfInfo='fitness',
                   numOfSegSites=sim.ALL_AVAIL, numOfMutants=sim.ALL_AVAIL),
                sim.PyEval(r'"%5d %s %5d %.6f %.6f %.6f %.6f" '
                   '% (gen, subPopSize, numOfSegSites, float(numOfMutants)/popSize, '
                    '(numOfMutants * 50. / numOfSegSites / popSize) if numOfSegSites else 0, '
                    ' meanOfInfo["fitness"], minOfInfo["fitness"])',
                    output=env.logger.info),
                ])
        ],
        matingScheme=sim.RandomMating(ops=
            sim.Recombinator(rates=recRate, loci=sim.ALL_AVAIL) if recRate != 0 else sim.MendelianGenoTransmitter(),
            subPopSize=demoModel),
        finalOps=[
            # revert fixed sites so that the final population does not have fixed sites
            sim.RevertFixedSites(),
            # statistics after evolution
            sim.Stat(popSize=True, meanOfInfo='fitness', minOfInfo='fitness',
                numOfSegSites=sim.ALL_AVAIL, numOfMutants=sim.ALL_AVAIL),
            sim.PyEval(r'"%5d %s %5d %.6f %.6f %.6f %.6f\n" '
                '% (gen, subPopSize, numOfSegSites, float(numOfMutants)/popSize, '
                '(numOfMutants * 50. / numOfSegSites / popSize) if numOfSegSites else 0, '
                ' meanOfInfo["fitness"], minOfInfo["fitness"])',
                output=env.logger.info),
            sim.PyEval(r'"Simulated population has %d individuals, %d segregation sites. '
                       r'There are on average %.1f mutants per individual. Mean allele frequency is %.4f%%.\n"'
                       r'% (popSize, numOfSegSites, numOfMutants / popSize, (numOfMutants * 50. / numOfSegSites/ popSize) if numOfSegSites else 0)',
                output=env.logger.info),
        ],
        gen = demoModel.num_gens
    )
    #
    env.logger.info('Population simulation takes %.2f seconds' % (time.clock() - startTime))
    return pop


