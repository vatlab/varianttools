#
# random fitness effect
#
from simuOpt import setOptions
setOptions(alleleType='mutant', quiet=True)
import simuPOP as sim

class RandomFitness:
    def __init__(self, scale=1):
        self.coefMap = {}
        self.scale = scale

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
                h = 0.5
            #
            self.coefMap[loc] = s, h
        if 0 in alleles:
            return 1. - self.scale * s * h
        else:
            return 1. - self.scale * s

class ConstantFitness(RandomFitness):
    def __init__(self, s=0.001, h=0.5, scale=1):
        RandomFitness.__init__(self, scale)
        self.s = s
        self.h = h

    def _newS(self, loc, alleles):
        # each mutant gets a penalty of s
        return self.s, self.h

class CustomizedFitness(RandomFitness):
    def __init__(self, func, scale=1):
        RandomFitness.__init__(self, scale)
        self.func = func
        self.h = 0.5   # default value that can be overwritten

    def _newS(self, loc, alleles):
        # each mutant gets a penalty of s
        return self.func(loc, alleles)

class GammaDistributedFitness(RandomFitness):
    def __init__(self, k=0.184, theta=0.160*2, h=0.5, scale=1):
        RandomFitness.__init__(self, scale)
        self.k = k
        self.theta = theta
        self.h = h

    def _newS(self, loc, alleles):
        return sim.getRNG().randGamma(self.k, self.theta)

class BoundedMixedGammaDistributedFitness(RandomFitness):
    #
    # k, theta: gamma distribution
    # h:        h=0.5 for additive model
    # p, a:     return constant a with given probability p
    #
    def __init__(self, k=0.184, theta=0.160*2,
        h=0.5, p=0.0186, a=0.0001, lower=0.0001, upper=0.1, scale=1):
        RandomFitness.__init__(self, scale)
        self.k = k
        self.theta = theta
        self.h = h
        self.p = p
        self.a = a
        self.lower = lower
        self.upper = upper
     
    def _newS(self, loc, alleles):
        if self.p > 0 and sim.getRNG().randUniform() < self.p:
            return self.a
        while True:
            s = sim.getRNG().randGamma(self.k, self.theta)
            if s >= self.lower and s <= self.upper:
                return s, self.h
        
class TrimmedMixedGammaDistributedFitness(RandomFitness):
    def __init__(self, k=0.562341, theta=0.01, h=0.5,
        p=0, a=0, lower=0.00001, upper=0.1, scale=1):
        RandomFitness.__init__(self, scale)
        self.k = k
        self.theta = theta
        self.h = h
        self.p = p
        self.a = a
        self.lower = lower
        self.upper = upper
     
    def _newS(self):
        if self.p > 0 and sim.getRNG().randUniform() < self.p:
            return self.a
        s = sim.getRNG().randGamma(self.k, self.theta)
        if s >= self.lower and s <= self.upper:
            return s
        elif s < self.lower:
            return self.lower, self.h
        else:
            return self.upper, self.h

class RandomFitnessSelector(sim.PyMlSelector):
    def __init__(self, model='gamma', coef=[], scale=1, mode=sim.MULTIPLICATIVE):
        #
        if callable(model):
            pyFunc = CustomizedFitness(model, **{'scale':scale})
        elif model == 'constant':
            pyFunc = ConstantFitness(*coef, **{'scale':scale})
        elif model.startswith('gamma'):
            pyFunc = GammaDistributedFitness(*coef, **{'scale':scale})
        elif model == 'bounded_gamma':
            pyFunc = BoundedMixedGammaDistributedFitness(*coef, **{'scale':scale})
        elif model == 'trimmed_gamma1':
            pyFunc = TrimmedMixedGammaDistributedFitness(*coef, **{'scale':scale})
        sim.PyMlSelector.__init__(self, func=pyFunc)

