#!/usr/bin/env python
#
# $File: simulation.py $
# $LastChangedDate: 2014-01-14 10:38:56 -0600 (Tue, 14 Jan 2014) $
# $Rev: 2505 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2014 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

import simuOpt
simuOpt.setOptions(alleleType='short', optimized=True, quiet=True, version='1.0.5')

import simuPOP as sim
from simuPOP.demography import *

from .utils import env, expandRegions, ProgressBar, RefGenome, existAndNewerThan, \
    calculateMD5, genesInRegions, codon_table, codon_table_reverse_complement, \
    dissectGene

import os, sys, math, time, random
from .pipeline import SkiptableAction
from .project import Project

if sys.version_info.major == 2:
    from ucsctools_py2 import tabixFetch
else:
    from ucsctools_py3 import tabixFetch


class ExtractFromVcf(SkiptableAction):
    '''Extract gentotypes at a specified region from a vcf file.'''
    def __init__(self, filenameOrUrl, regions, output):
        self.filenameOrUrl = filenameOrUrl
        self.regions = regions
        SkiptableAction.__init__(self, cmd='ExtractFromVcf {} {} {}'.format(filenameOrUrl, regions, output),
            output=output, ignoreInput=True)

    def _execute(self, ifiles, pipeline):
        tabixFetch(self.filenameOrUrl, [], self.output[0], True)
        for r in expandRegions(self.regions, pipeline.proj):
            region = '{}:{}-{}'.format(r[0], r[1], r[2])
            env.logger.info('Retriving genotype for region chr{}{}'.format(region,
                ' ({})'.format(r[3] if r[3] else '')))
            tabixFetch(self.filenameOrUrl, [region], self.output[0], False)
        
class VcfToPop(SkiptableAction):
    '''Check out of of an command, and check if it matches a particular
    pattern. The pipeline will exit if fail is set to True (default).'''
    def __init__(self, regions, output):
        self.regions = regions
        SkiptableAction.__init__(self, cmd='VcfToPop {} {}\n'.format(regions, output),
            output=output)

    def _execute(self, ifiles, pipeline):
        # translate regions to simuPOP ...
        lociPos = {}
        for r in expandRegions(self.regions, pipeline.proj):
            if r[0] in lociPos:
                lociPos[r[0]].extend(range(r[1], r[2] + 1))
            else:
                lociPos[r[0]] = range(r[1], r[2] + 1)
        chroms = lociPos.keys()
        chroms.sort()
        #
        # create a dictionary of lociPos->index on each chromosome
        lociIndex = {}
        for chIdx,ch in enumerate(chroms):
            lociIndex[chIdx] = {pos:idx for idx,pos in enumerate(lociPos[ch])}
        #
        # number of individuals? 629
        pop = sim.Population(size=629, loci=[len(lociPos[x]) for x in chroms],
            chromNames = chroms, lociPos=sum([lociPos[x] for x in chroms], []))
        #
        # we assume 0 for wildtype, 1 for genotype
        #
        # extract genotypes
        allele_map = {'0': 0, '1': 1, '2': 1, '.': 0}
        mutantCount = 0
        with open(ifiles[0], 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                chr = pop.chromNames().index(fields[0])
                pos = lociIndex[chr][int(fields[1])]
                #
                for ind, geno in enumerate(fields[10:]):
                    if geno[0] != 0:
                        pop.individual(ind).setAllele(1, pos, 0, chr)
                        mutantCount += 1
                    if geno[2] != 0:
                        pop.individual(ind).setAllele(1, pos, 1, chr)
                        mutantCount += 1
        env.logger.info('{} mutants imported'.format(mutantCount))
        pop.save(self.output[0])

class PopToVcf(SkiptableAction):
    def __init__(self, output, sample_names=[]):
        self.sample_names = sample_names
        SkiptableAction.__init__(self, cmd='PopToVcf {} {}\n'.format(sample_names, output),
            output=output)

    def _execute(self, ifiles, pipeline):
        # import variants to the current project
        # translate 0,1,2,3 to A,C,G,T
        alleleMap = {
            'A': {0: 'A', 1: 'C', 2: 'G', 3: 'T'},
            'C': {0: 'C', 1: 'G', 2: 'T', 3: 'A'},
            'G': {0: 'G', 1: 'T', 2: 'A', 3: 'C'},
            'T': {0: 'T', 1: 'A', 2: 'C', 3: 'G'},
        }
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        # output genotype
        with open(self.output[0], 'w') as vcf:
            vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
            if self.sample_names:
                if len(self.sample_names) != pop.popSize():
                    raise ValueError('Sample names, if specified, should be assigned to'
                        'all {} individuals.'.format(pop.popSize()))
                vcf.write('\t'.join(self.sample_names) + '\n')
            else:
                vcf.write('\t'.join(['S_{}'.format(x) for x in range(1, pop.popSize()+1)]) + '\n')
            #
            # get reference genome
            refGenome = RefGenome(pipeline.proj.build)
            sim.stat(pop, alleleFreq=sim.ALL_AVAIL, vars='alleleNum')
            nAlleles = 2*pop.popSize()
            segregated = [loc for loc,nums in pop.dvars().alleleNum.items() if nums[0] == nAlleles]
            env.logger.info('Genetic variants identified on {} loci'.format(len(segregated)))
            prog = ProgressBar('Exporting simulated population', pop.totNumLoci())
            for chr in range(pop.numChrom()):
                chr_name = pop.chromName(chr)
                for loc in range(pop.chromBegin(chr), pop.chromEnd(chr)):
                    if pop.dvars().alleleNum[loc][0] == nAlleles:
                        continue
                    pos = int(pop.locusPos(loc))
                    ref = refGenome.getBase(chr_name, pos)
                    # genotypes
                    geno1 = [ind.allele(loc, 0) for ind in pop.individuals()]
                    geno2 = [ind.allele(loc, 1) for ind in pop.individuals()]
                    #
                    alt = list((set(geno1) | set(geno2)) - set([0]))
                    #
                    if len(alt) == 0:
                        env.logger.error('No alternative allele at locus {}:{}'.format(chr_name, pos))
                        continue
                    elif len(alt) == 1:
                        # easier...
                        vcf.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t'.format(
                            chr_name, pos, ref, alleleMap[ref][alt[0]]))
                        vcf.write('\t'.join(['{}/{}'.format(x if x==0 else 1, 
                            y if y == 0 else 1) for x,y in zip(geno1, geno2)]) + '\n')
                    else:
                        # we need to figure out the index of geno in alt
                        vcf.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t'.format(
                            chr_name, pos, ref, ','.join([alleleMap[ref][a] for a in alt])))
                        vcf.write('\t'.join(['{}/{}'.format(x if x==0 else alt.index(x) + 1, 
                            y if y == 0 else alt.index(y) + 1) for x,y in zip(geno1, geno2)]) + '\n')
                    #
                    prog.update(loc)
            prog.done()
        #
        # output phenotype
        with open(self.output[1], 'w') as phe:
            fields = pop.infoFields()
            phe.write('sample_name\t' + '\t'.join(fields) + '\n')


class ProteinSelector(sim.PySelector):
    def __init__(self, regions, s_missense=0.001, s_stoploss=0.002, s_stopgain=0.01):
        '''A protein selection operator that, for specified regions
            1. find coding regions and pass them to PySelector
            2. find amino acid change of each individual
            3. return fitness caused by change of amino acid

        s_missense: selection coefficient for missense (nonsynonymous mutations)
        s_stoploss: selection coefficient for stoploss muation (elongate protein) 
        s_stopgain: selection coefficient for stopgan muation (premature coding of protein)

        Selection coefficient should be a single number (fixed s, with fitness 1-s).
        The fitness of multiple amino acid change will be Prod(1-si) even if two changes
        are at the same location (that is to say, a homozygote change will have fitness
        1-2*s-s*2, which is close to an additive model for small s.
        '''
        self.s_missense = s_missense
        self.s_stoploss = s_stoploss
        self.s_stopgain = s_stopgain
        #
        # 
        #                         |<- coding                                coding    ->|
        #  exon 1     intron 1   exon 2        intron 2     exon 3      intron2     exon 4
        # [---------]---------[---|--------]-------------[------------]----------[------|-----]
        #
        # We are only interested in coding sequence (mRNA). 
        #
        #                          xxxxxxxx              xxxxxxxxxxxxx            xxxxxx
        # 
        # We pass location of coding sequences to PySelector, the call back function will receive
        # genotypes at these loci, and we need to remap these genotypes to the genome.
        #
        # all_loci: pos->idx in region
        # ......................................................................................
        #
        # coding_loci: positions of xxxxx
        #
        # coding_idx2pos:    index within coding loci ==> pos
        # coding_pos2idx:    pos within coding loci ==> index
        # nCoding:           number of coding loci
        # 
        #
        # FIXME: Right now we only handle regions on the same chromosome. Will fix it later.
        #
        # self.codon_info stores information  about 1 or more codon at codng regions (excluding introns)
        self.codon_info = {}
        # self.coding_base stores reference sequence at coding regions
        self.coding_base = {}
        #
        with Project(verbosity='1') as proj:
            # expand user provided regions to one or more (chr,start,end,comment)
            self.regions = expandRegions(regions, proj)
            #
            # all loci contrains all pos in the region and their indexes
            all_loci = set()
            for reg in self.regions:
                all_loci = all_loci.union(range(reg[1], reg[2]+1))
            all_loci = {y:x for x,y in enumerate(sorted(all_loci))}
            #
            ref = RefGenome(proj.build)
            genes = genesInRegions(self.regions, proj)
            env.logger.info('{} genes are identified in the simulated region.'
                .format(len(genes)))
            # there can be multiple genes (isoforms) in the same regions
            for gene in genes:
                stru = dissectGene(gene, proj)
                pos = []
                seq = ''
                # find the coding regions (xxxx regions) of each gene
                for reg in stru['coding']:
                    # get reference sequence and positions
                    seq += ref.getSequence(reg[0], reg[1], reg[2])
                    pos.extend(range(reg[1], reg[2]+1))
                env.logger.info('Length of coding regions of {}: {}'.format(
                    gene, len(pos)))
                # now, try to divide coding regions by codon. Note that a codon
                # can consist of nucleotie across two exon.
                for idx, (p, s) in enumerate(zip(pos, seq)):
                    if idx % 3 == 0:
                        # the complete codon must be in the simulated region. We do not handle
                        # partial codon because we cannot control mutations outside of the specified
                        # region
                        if (p not in all_loci) or (pos[idx+1] not in all_loci) or (pos[idx+2] not in all_loci):
                            continue
                        # information about the codon: p0, p1, p2, aa, strand
                        codon = (p, pos[idx+1], pos[idx+2],
                            codon_table[s + seq[idx+1] + seq[idx+2]] if stru['strand'] == '+' else
                            codon_table_reverse_complement[s + seq[idx+1] + seq[idx+2]] ,
                            stru['strand'])
                    # record reference sequence
                    self.coding_base[p] = s
                    # other two positions share the same codon
                    # because of there can be multiple genes, one basepair can be in multiple codon
                    if p in self.codon_info:
                        if codon not in self.codon_info[p]:
                            self.codon_info[p].append(codon)
                    else:
                        self.codon_info[p] = [codon]
        #
        # mutated loci
        coding_loci = sorted(self.codon_info.keys())
        self.coding_idx2pos = {x:y for x,y in enumerate(coding_loci)}
        self.coding_pos2idx = {y:x for x,y in enumerate(coding_loci)}
        self.nCoding = len(coding_loci)
        # remove all positions that are not in regions
        env.logger.info('{} out of {} bp are in coding regions of genes {}'.format(
            self.nCoding, len(all_loci), ', '.join(genes)))
        #
        if not genes:
            env.logger.warning('Specified region does not contain any gene. A neutral model will be used.')
            sim.PySelector.__init__(self, func=self._neutral, loci=[])
        else:
            # we need to send simuPOP to indexes of coding loci within the specified region
            sim.PySelector.__init__(self, func=self._select, loci=[all_loci[x] for x in coding_loci])
        # 
        # a cache for all fitness values
        #self.fitness_cache = {}
        #
        # the meaning of mutation is different according to ref sequence
        self.mutant_map = {
            'A': {0: 'A', 1: 'C', 2: 'G', 3: 'T'},
            'C': {0: 'C', 1: 'G', 2: 'T', 3: 'A'},
            'G': {0: 'G', 1: 'T', 2: 'A', 3: 'C'},
            'T': {0: 'T', 1: 'A', 2: 'C', 3: 'G'},
        }

    def _neutral(self):
        return 1

    def _select(self, geno):
        # 
        #try:
        #    return self.fitness_cache[tuple([(x,y) for x,y in enumerate(geno) if y!=0])]
        #except KeyError:
        #    # if fitness for gene is not cached, calculate it
        #    pass
        # we can not divide the sequence into triplets because it is possible that a nucleotide
        # belong to multiple codon with different locations.
        mutated = [idx for idx,x in enumerate(geno) if x != 0]
        # 
        # the same aa change can be caused by two mutations, we need to keep only one
        # aa_change with the same (aa, naa, ploidy, p0)
        aa_change = set()
        for m in mutated:
            # first homologous copy
            aa_change_at_m = []
            if m < self.nCoding:
                # we need to map index in coding region to their original position
                for p0, p1, p2, aa, s in self.codon_info[self.coding_idx2pos[m]]:
                    # p0: location
                    # coding_pos2idx[p0]: index in coding region
                    # geno: mutation (0 for wildtype)
                    # mutant_map: map mutation to nucleotide
                    codon = self.mutant_map[self.coding_base[p0]][geno[self.coding_pos2idx[p0]]] + \
                            self.mutant_map[self.coding_base[p1]][geno[self.coding_pos2idx[p1]]] + \
                            self.mutant_map[self.coding_base[p2]][geno[self.coding_pos2idx[p2]]]
                    naa = codon_table[codon] if s == '+' else codon_table_reverse_complement[codon]
                    # if this is a real change
                    if naa != aa:
                        aa_change_at_m.append((aa, naa, 0, p0))
            else:
                for p0, p1, p2, aa, s in self.codon_info[self.coding_idx2pos[m - self.nCoding]]:
                    # p0: location
                    # coding_pos2idx[p0]: index in coding region
                    # geno: mutation (0 for wildtype)
                    # mutant_map: map mutation to nucleotide
                    codon = self.mutant_map[self.coding_base[p0]][geno[self.nCoding + self.coding_pos2idx[p0]]] + \
                            self.mutant_map[self.coding_base[p1]][geno[self.nCoding + self.coding_pos2idx[p1]]] + \
                            self.mutant_map[self.coding_base[p2]][geno[self.nCoding + self.coding_pos2idx[p2]]]
                    naa = codon_table[codon] if s == '+' else codon_table_reverse_complement[codon]
                    # if this is a real change
                    if naa != aa:
                        aa_change_at_m.append((aa, naa, 1, p0))
            #
            # The same mutation can cause multiple aa change for different genes, we need to
            # keep only the most damaging one.
            if len(aa_change_at_m) == 1:
                aa_change.add(aa_change_at_m[0])
            elif len(aa_change_at_m) > 1:
                # find the most damaging one
                stopgain = [x for x in aa_change_at_m if x[1] == '*']
                if stopgain:
                    aa_change.add(stopgain[0])
                stoploss = [x for x in aa_change_at_m if x[0] == '*']
                if stoploss:
                    aa_change.add(stoploss[0])
                else:
                    aa_change.add(aa_change_at_m[0])
        #
        # we assume a multiplicative model
        fitness = 1
        for aa, naa, ploidy, start in aa_change:
            # stoploss
            if aa == '*': 
                fitness *= 1 - self.s_stoploss
            elif naa == '*':
                fitness *= 1 - self.s_stopgain
            else:
                fitness *= 1 - self.s_missense
        #self.fitness_cache[tuple([(x,y) for x,y in enumerate(geno) if y!=0])] = fitness
        return fitness
    

class EvolvePop:
    def __init__(self, mu, recRate, selector, demoModel, output):
        self.mu = mu
        self.selector = selector
        self.recRate = recRate
        self.demoModel = demoModel
        self.output = [output]
        #SkiptableAction.__init__(self, cmd='EvolvePop {} {} {}\n'
        #    .format(mu, recRate, output),
        #    output=output)

    def __call__(self, ifiles, pipeline):
        pop = sim.loadPopulation(ifiles[0])
        refGenome = RefGenome(pipeline.proj.build)
        
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
                sim.AcgtMutator(model='K80', rate=[self.mu, 2], loci=base['A']),
                sim.AcgtMutator(model='K80', rate=[self.mu, 2], loci=base['C'],
                    # base is C=0, 0,1,2,3 to C,G,T,A (1,2,3,0)
                    mapIn=[1,2,3,0], mapOut=[3,0,1,2]),
                sim.AcgtMutator(model='K80', rate=[self.mu, 2], loci=base['G'],
                    # base is G=0, 0,1,2,3 to G,T,A,C (2,3,0,1)
                    mapIn=[2,3,0,1], mapOut=[2,3,0,1]),
                sim.AcgtMutator(model='K80', rate=[self.mu, 2], loci=base['T'],
                    # base is T=0, 0,1,2,3 to T,A,C,G (3,0,1,2)
                    mapIn=[3,0,1,2], mapOut=[1,2,3,0]),
                # selection on all loci
                self.selector,
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
                sim.Recombinator(rates=self.recRate, loci=sim.ALL_AVAIL) if self.recRate != 0 else sim.MendelianGenoTransmitter(),
                subPopSize=self.demoModel),
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
            gen = self.demoModel.num_gens
        )
        #
        env.logger.info('Population simulation takes %.2f seconds' % (time.clock() - startTime))
        pop.save(self.output[0])
        return self.output



electedCase = 0
selectedControl = 0
discardedInds = 0

alpha = -5
beta1 = 0.20
beta2 = 0.4
beta3 = 0.4
gamma1 = 0.2
gamma2 = 0.4


def _selectInds(off, param):
    'Determine if the offspring can be kept.'
    e = random.randint(0, 1)
    g1 = off.allele(param[0], 0) + off.allele(param[0], 1)
    g2 = off.allele(param[1], 0) + off.allele(param[1], 1)
    logit = alpha + beta1*g1 + beta2*g2 + beta3*g1*g2 + gamma1*e*g1 + gamma2*e*g2
    affected = random.random() < (1 / (1. + math.exp(-logit)))
    global selectedCase, selectedControl, discardedInds
    if affected:
        if selectedCase < 1000:
            off.setAffected(True)
            selectedCase += 1
            return True
    else:
        if selectedControl < 1000:
            selectedControl += 1
            off.setAffected(False)
            return True
    discardedInds += 1
    return False

def generateSample(numCase=1000, numCtrl=1000, logger=None):
    if logger:
        logger.info('Generating %d cases and %d controls...' % (numCase, numCtrl))
    pop = loadPopulation('ex2_expanded.pop')
    loci = pop.lociByNames(DPL)
    pop.evolve(
        matingScheme=RandomMating(
            ops=[
                MendelianGenoTransmitter(),
                # an individual will be discarded if _selectInds returns False
                PyOperator(func=_selectInds, param=loci)
            ], 
            subPopSize=2000
        ),
        gen = 1
    )
    pop.save('ex2_sample.pop')
    if logger:
        logger.info('Disease prevalence: %.4f' % (1000. / (2000 + discardedInds)))



 

# IDs of the parents of selected offspring
parentalIDs = set()
# number of discarded individuals, used to calculate disease prevalence.
discardedInds = 0


alpha = -0.5
beta = -1.0
def _myPenetrance(geno):
    g = geno[0] + geno[1]
    logit = alpha + beta*g
    return math.exp(logit) / (1. + math.exp(logit))
    #return [0.5, 0.4, 0.01][g]

def _selectTrio(off, param):
    'Determine if the offspring can be kept.'
    global discardedInds, parentalIDs
    if off.affected() and not (off.father_id in parentalIDs or off.mother_id in parentalIDs):
        off.setAffected(True)
        parentalIDs |= set([off.father_id, off.mother_id])
        return True
    else:
        discardedInds += 1
        return False

def generateTrioSamples(logger=None):
    '''
    Generate trio samples to be analyzed by HelixTree.
    '''
    for popFile in ['ex3_shortsweep.pop', 'ex3_longsweep.pop']:
        if logger:
            logger.info('Loading population ' + popFile)
        pop = loadPopulation(popFile)
        locus = pop.locusByName(DPL)
        # save a map file
        map = open('ex3.map', 'w')
        print >> map, 'CHROMOSOME MARKER POSITION'
        for loc in range(pop.totNumLoci()):
            print >> map, pop.chromName(0), pop.locusName(loc), pop.locusPos(loc)/1e6
        map.close()
        dat = open('ex3.dat', 'w')
        print >> dat, 'A disease'
        for loc in range(pop.totNumLoci()):
            print >> dat, 'M', pop.locusName(loc)
        dat.close()
        pop.addInfoFields(['ind_id', 'father_id', 'mother_id'])
        # keep parental generation
        pop.setAncestralDepth(1)
        # give everyone an unique ID
        tagID(pop, reset=True)
        stat(pop, genoFreq=locus)
        print 'Genotype frequency: ', pop.dvars().genoFreq[locus]
        pop.evolve(
            preOps = PyPenetrance(func=_myPenetrance, loci=[locus]),
            matingScheme=RandomMating(
                ops=[
                    # pass genotype
                    MendelianGenoTransmitter(),
                    # assign new ID to offspring
                    IdTagger(),
                    # record the parent of each offspring
                    PedigreeTagger(),
                    # determine offspring affection status
                    PyPenetrance(func=_myPenetrance, loci=[locus]),
                    # discard the offspring if it is not affected
                    # or if parents have been chosen
                    PyOperator(func=_selectTrio, param=locus)
                ], subPopSize=1000
            ),
            gen = 1
        )
        global selectedInds, discardedInds, parentalIDs
        sample = pop.extractIndividuals(IDs=list(parentalIDs) + list(pop.indInfo('ind_id')))
        if logger:
            logger.info('Saving sample to file ' + popFile.replace('.pop', '_sample.pop'))
        sample.save(popFile.replace('.pop', '_sample.pop'))
        if logger:
            logger.info('Disease prevalence: %.4f' % (1000. / (1000 + discardedInds)))
        # save file in PED format so I need to change family ID...
        if logger:
            logger.info('Saving sample to file ' + popFile.replace('.pop', '.ped'))
        # write to merlin format
        csv = open(popFile.replace('.pop', '.ped'), 'w')
        #print >> csv, ' '.join(pop.lociNames())
        def genoString(ind):
            alleles = []
            for i in range(ind.totNumLoci()):
                alleles.extend([str(ind.allele(i, 0)+1), str(ind.allele(i, 1)+1)])
            return ' '.join(alleles)
        famid = 1
        id = 1
        for ind in sample.individuals():
            father = sample.indByID(ind.father_id)
            mother = sample.indByID(ind.mother_id)
            print >> csv, famid, id, 0, 0, '1', '2' if father.affected() else '1', genoString(father)
            print >> csv, famid, id+1, 0, 0, '2', '2' if mother.affected() else '1', genoString(mother)
            print >> csv, famid, id+2, id, id+1, '1' if ind.sex() == MALE else '2', '2' if ind.affected() else '1', genoString(ind)
            famid += 1
            id += 3
        csv.close()
        # for next population
        parentalIDs = set()
        selectedInds = 0
        discardedInds = 0
