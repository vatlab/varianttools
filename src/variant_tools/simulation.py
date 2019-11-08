#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
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
'''This module defines functions and actions for Variant Simulation
Tools.
'''

import datetime
import os
import random
import time
from collections import defaultdict

import simuOpt
simuOpt.setOptions(
    alleleType='mutant', optimized=True, quiet=True, version='1.1.4')

import simuPOP as sim
from simuPOP.sampling import drawRandomSample

from .pipeline import PipelineAction
from .project import Project
from .ucsctools import tabixFetch
from .utils import (ProgressBar, RefGenome, codon_table,
                    codon_table_reverse_complement, decompressGzFile,
                    dissectGene, downloadFile, env, expandRegions,
                    genesInRegions)








def FineScaleRecombinator(regions=None, scale=1, defaultRate=1e-8, output=None):
    '''For specified regions of the chromosome, find genetic locations of
    all loci using a genetic map downloaded from HapMap. If no genetic
    map is used, a default recombination rate (per bp) is used. If a
    output file is specified, the physical/genetic map will be written
    to the file. If each element in regions has only length two, it is
    assumed to be a single-locus region. Finally, if a population object
    is specified, the regions will be obtained automatically from all loci
    of the population object.
    '''
    if scale == 0:
        return sim.MendelianGenoTransmitter()
    #
    lociPos = {}
    if isinstance(regions, str):
        regions = expandRegions(regions)
    if isinstance(regions, (tuple, list)):
        if not regions:
            return sim.MendelianGenoTransmitter()
        if len(regions[0]) > 2:
            for reg in regions:
                if reg[0] in lociPos:
                    lociPos[reg[0]].extend(list(range(reg[1], reg[2] + 1)))
                else:
                    lociPos[reg[0]] = list(range(reg[1], reg[2] + 1))
        elif len(regions[0]) == 2:
            # single-locus regions
            for (chr, pos) in regions:
                if chr in lociPos:
                    lociPos[chr].append(pos)
                else:
                    lociPos[chr] = [pos]
        else:
            raise ValueError('Incorrect parameter regions: {}'.format(regions))
    else:
        for ch in range(regions.numChrom()):
            lociPos[regions.chromName(ch)] = [
                int(regions.locusPos(x))
                for x in range(regions.chromBegin(ch), regions.chromEnd(ch))
            ]
    #
    chroms = list(lociPos.keys())
    chroms.sort()
    #
    # download the map
    Recom_URL = 'ftp://ftp.hapmap.org/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz'
    recom = downloadFile(Recom_URL)
    recom_files = decompressGzFile(recom)
    rec_rates = []
    for ch, ch_name in enumerate(chroms):
        try:
            recom_file = [
                x for x in recom_files if os.path.basename(x) ==
                'genetic_map_GRCh37_chr{}.txt'.format(ch_name)
            ][0]
        except Exception as e:
            raise RuntimeError(
                'Could not find recombination map for chromosome {}: {}'.format(
                    ch_name, e))
        #
        physicalPos = lociPos[ch_name]
        geneticPos = {}
        mapPoints = []
        with open(recom_file) as recom:
            recom.readline()
            for line in recom:
                try:
                    fields = line.split()
                    pos = int(fields[1])
                    if pos >= physicalPos[0] and pos <= physicalPos[-1]:
                        geneticPos[pos] = float(fields[3])
                        mapPoints.append((pos, float(fields[3])))
                except Exception as e:
                    env.logger.warning(e)
        env.logger.info(
            'Map distance of {} markers within region {}:{}-{} are found'
            .format(len(geneticPos), ch_name, physicalPos[0], physicalPos[-1]))
        if len(geneticPos) <= 1:
            # using generic distance
            env.logger.info(
                'Using a default recombination rate {}'.format(defaultRate))
            geneticPos = {
                x: (x - physicalPos[0]) * defaultRate for x in physicalPos
            }
            rec_rates.extend([
                geneticPos[physicalPos[i + 1]] - geneticPos[physicalPos[i]]
                for i in range(len(physicalPos) - 1)
            ] + [0.5])
            continue
        #
        # step 1: find the first two positions and calculate the leading portion
        idx = 0
        if physicalPos[0] not in geneticPos:
            #
            # ...........P0.....P1
            p0 = physicalPos[1]
            while p0 not in geneticPos:
                p0 += 1
            p1 = p0 + 1
            while p1 not in geneticPos:
                p1 += 1
            #
            # now, rate
            rate = (geneticPos[p1] - geneticPos[p0]) / (p1 - p0)
            #
            # fill up to P1
            while idx < len(physicalPos):
                p = physicalPos[idx]
                if p <= p1:
                    geneticPos[p] = geneticPos[p1] - (p1 - p) * rate
                    idx += 1
                else:
                    break
        else:
            p1 = physicalPos[0]
            idx = 1
        #
        # idx is the currently processed physicalPos
        #
        # step 2: find the middle ones, starting from p1
        #
        # p1 ... p2 .....
        while p1 <= physicalPos[-1]:
            # can we found P2?
            p2 = p1 + 1
            while p2 not in geneticPos and p2 < physicalPos[-1]:
                p2 += 1
            # two cases
            if p2 in geneticPos:
                # good
                rate = (geneticPos[p2] - geneticPos[p1]) / (p2 - p1)
                while idx < len(physicalPos):
                    p = physicalPos[idx]
                    if p <= p2:
                        geneticPos[p] = geneticPos[p1] + (p - p1) * rate
                        idx += 1
                    else:
                        break
                p1 = p2
            else:
                # not in use the last rate
                for p in physicalPos[idx:]:
                    geneticPos[p] = geneticPos[p1] + (p - p1) * rate
                break
        #
        start_pos = len(rec_rates)
        rec_rates.extend([
            geneticPos[physicalPos[i + 1]] - geneticPos[physicalPos[i]]
            for i in range(len(physicalPos) - 1)
        ] + [0.5])
        env.logger.info(
            'Total map distance from {}:{}-{} ({} bp) is {:.4e} cM (r={:.4e}/bp)'
            .format(ch_name, physicalPos[0], physicalPos[-1],
                    physicalPos[-1] - physicalPos[0] + 1,
                    geneticPos[physicalPos[-1]] - geneticPos[physicalPos[0]],
                    (geneticPos[physicalPos[-1]] - geneticPos[physicalPos[0]]) /
                    (physicalPos[-1] - physicalPos[0])))
        #
        assert (len(lociPos[ch_name]) == len(rec_rates[start_pos:]))
        if output is not None:
            with open(output, 'w' if ch == 0 else 'a') as gmap:
                if ch == 0:
                    gmap.write('chr\tpos\tgenetic_pos\trate\n')
                lp, lg = None, None
                for p, g in mapPoints:
                    if lp is None:
                        gmap.write('#{}\t{}\t{}\t0\n'.format(ch_name, p, g))
                    else:
                        gmap.write('#{}\t{}\t{}\t{}\n'.format(
                            ch_name, p, g, (g - lg) / (p - lp)))
                    lp, lg = p, g
                for p, r in zip(lociPos[ch_name],
                                [0] + rec_rates[start_pos:start_pos +
                                                len(lociPos[ch_name]) - 1]):
                    gmap.write('{}\t{}\t{}\t{}\n'.format(
                        ch_name, p, geneticPos[p], r))
    #
    rec_rates = [x * scale if x * scale < 0.5 else 0.5 for x in rec_rates]
    # use loci position to return list of loci
    return sim.Recombinator(
        rates=rec_rates,
        loci=sum([[(ch, pos) for pos in lociPos[ch]] for ch in chroms], []))


class OutputPopulationStatistics(PipelineAction):

    def __init__(self, mut_count=None):
        output = []
        self.mut_count = mut_count
        if self.mut_count:
            output.append(self.mut_count[0])
        self.output = output
        PipelineAction.__init__(
            self, cmd='PopStat output={}\n'.format(output), output=output)

    def _execute(self, ifiles, pipeline):
        #
        env.logger.info('Loading population from {}'.format(ifiles[0]))
        self.pop = sim.loadPopulation(ifiles[0])
        if self.mut_count is not None:
            self.count_mutants()
        #
        #result = self.LD_curve(pop)
        #with open(self.output[0], 'w') as ld_out:
        #    for ld_stat in result.keys:
        #        for sp in result[ld_stat].keys():
        #            ld_out.write('{}\t{}\t{}\n'.format(ld_stat, sp,
        #                '\t'.join(['{:.4f}'.format(x) for x in result[ld_stat][sp]])))
        return True

    def count_mutants(self):
        #
        sz = min(self.pop.popSize(), self.mut_count[1])
        sample = drawRandomSample(self.pop, sz)
        #
        sim.stat(sample, alleleFreq=sim.ALL_AVAIL, vars='alleleNum')
        #
        with open(self.mut_count[0], 'w') as out:
            out.write('#index\tchr\tposition\t#0\t#1\t#2\t#3\n')
            # output allele counts at each locus
            for ch in range(self.pop.numChrom()):
                ch_name = self.pop.chromName(ch)
                for i in range(self.pop.chromBegin(ch), self.pop.chromEnd(ch)):
                    if int(sample.dvars().alleleNum[i][0]) != sz * 2:
                        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                            i + 1, ch_name, int(self.pop.locusPos(i)),
                            int(sample.dvars().alleleNum[i][0]),
                            int(sample.dvars().alleleNum[i][1]),
                            int(sample.dvars().alleleNum[i][2]),
                            int(sample.dvars().alleleNum[i][3])))

    def LD_curve(self, pop, maxDist=500000, binDist=1000, maxPairs=10000):
        '''Calculating LD as a function of distance. If multiple chromosomes are
        present, they are averaged.

        pop
            population to observe

        maxDist
            maximal distance. If two markers are located greater than
            this distance, their LD is not calculated. If markers are
            located by basepair, this value is usually 500kbp ~ 0.5 cM.

        binDist
            maxDist is separated into bins of this width. maxDist/binDist
            determines the length of return value. This value is usually
            1k for a maxDist of 500k.

        return a dictionary of res['LD_prime' or 'R2'][spName].
        '''
        env.logger.info('Calculating LD values for a population, with maxDist %d and binDist %d' \
                % (maxDist, binDist))
        sim.stat(pop, alleleFreq=sim.ALL_AVAIL)
        seg_sites = set([
            x for x in range(pop.totNumLoci())
            if len(pop.dvars().alleleFreq[x]) == 1
        ])
        env.logger.info('LD from {} segregating sites'.format(len(seg_sites)))
        numBin = maxDist / binDist
        result = {'LD_prime': {}, 'R2': {}}
        for sp in pop.subPopNames():
            result['LD_prime'][sp] = []
            result['R2'][sp] = []
        pairs = [[] for dist in range(numBin)]
        for ch in range(pop.numChrom()):
            seg_sites_ch = [
                x for x in range(pop.chromBegin(ch), pop.chromEnd(ch))
                if x in seg_sites
            ]
            # chName = pop.chromName(ch)
            for idx, i in enumerate(seg_sites_ch):
                for j in seg_sites_ch[idx + 1:]:
                    dist = int(pop.locusPos(j) - pop.locusPos(i))
                    if dist < maxDist and len(pairs[dist / binDist]) < maxPairs:
                        pairs[dist / binDist].append([i, j])
        sim.Stat(pop, popSize=True)
        progress = ProgressBar('Calculating average LD as a function of distance (%d %dbp bins)' \
            % (numBin, binDist), numBin)
        for dist in range(numBin):
            # clear previous value
            pop.dvars().LD_prime = None
            pop.dvars().R2 = None
            for sp in range(pop.numSubPop()):
                pop.dvars(sp).LD_prime = None
                pop.dvars(sp).R2 = None
            if len(pairs[dist]) == 0:
                for sp in range(pop.numSubPop()):
                    result['LD_prime'][pop.subPopName(sp)].append(0)
                    result['R2'][pop.subPopName(sp)].append(0)
                continue
            sim.Stat(pop, LD=pairs[dist], LD_param={'stat': ['LD_prime', 'R2']})
            avgLD = [0] * pop.numSubPop()
            avgR2 = [0] * pop.numSubPop()
            for sp in range(pop.numSubPop()):
                LD_var = pop.dvars(sp).LD_prime
                R2_var = pop.dvars(sp).R2
                for pair in pairs[dist]:
                    avgLD[sp] += LD_var[pair[0]][pair[1]]
                    avgR2[sp] += R2_var[pair[0]][pair[1]]
                result['LD_prime'][pop.subPopName(sp)].append(avgLD[sp] /
                                                              len(pairs[dist]))
                result['R2'][pop.subPopName(sp)].append(avgR2[sp] /
                                                        len(pairs[dist]))
            progress.update(dist + 1)
        progress.done()
        return result


class RefGenomeMutator(sim.PyOperator):
    '''This operator uses allele 0 as the reference allele at different loci and
    applies different @@ActgMutator@@ according to the actual nuclotides on the
    reference genome. For example, if the reference genome is
    @@ACCCCTTAGG@@, it is represented by haplotype @@0000000000@@, and
    @@1000000000@@ and @@0200000000@@ represents @@CCCCCTTAGG@@ and
    @@ATCCCTTAGG@@ respectively (A->C->G->T). If you apply a Kimura's
    2-parameter (K80) model to the reference genome, it will act differently
    at different location of the genome.
    '''

    def __init__(self, regions, model, rate):
        #
        base = {'A': [], 'C': [], 'G': [], 'T': [], 'N': []}
        with Project() as proj:
            refGenome = RefGenome(proj.build)
            idx = 0
            for reg in expandRegions(regions, proj):
                for s in refGenome.getSequence(reg[0], reg[1], reg[2]):
                    base[s].append(idx)
                    idx += 1
        #
        env.logger.info(
            'Simulated regions with {:,} basepair have {:,} A, '
            '{:,} C, {:,} G, {:,} T, and {:,} N on reference genome.'.format(
                idx, len(base['A']), len(base['C']), len(base['G']),
                len(base['T']), len(base['N'])))
        self.mutators = [
            sim.AcgtMutator(model=model, rate=rate, loci=base['A']),
            sim.AcgtMutator(
                model=model,
                rate=rate,
                loci=base['C'],
                # base is C=0, 0,1,2,3 to C,G,T,A (1,2,3,0)
                mapIn=[1, 2, 3, 0],
                mapOut=[3, 0, 1, 2]),
            sim.AcgtMutator(
                model=model,
                rate=rate,
                loci=base['G'],
                # base is G=0, 0,1,2,3 to G,T,A,C (2,3,0,1)
                mapIn=[2, 3, 0, 1],
                mapOut=[2, 3, 0, 1]),
            sim.AcgtMutator(
                model=model,
                rate=rate,
                loci=base['T'],
                # base is T=0, 0,1,2,3 to T,A,C,G (3,0,1,2)
                mapIn=[3, 0, 1, 2],
                mapOut=[1, 2, 3, 0])
        ]
        sim.PyOperator.__init__(self, func=self._mutate)

    def _mutate(self, pop):
        for m in self.mutators:
            m.apply(pop)
        return True


class MutantInfo:

    def __init__(self, raw_regions):
        '''Identify codon of a region'''
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
        # ......................................................................................
        # coding_loci: positions of xxxxx
        #
        #
        self._codon_info = None
        self._coding_base = None
        # codon_info stores information  about 1 or more codon at codng regions (excluding introns)
        self.codon_info = {}
        # self.coding_base stores reference sequence at coding regions
        self.coding_base = {}
        #
        # the meaning of mutation is different according to ref sequence
        self.mutant_map = {
            'A': {
                0: 'A',
                1: 'C',
                2: 'G',
                3: 'T'
            },
            'C': {
                0: 'C',
                1: 'G',
                2: 'T',
                3: 'A'
            },
            'G': {
                0: 'G',
                1: 'T',
                2: 'A',
                3: 'C'
            },
            'T': {
                0: 'T',
                1: 'A',
                2: 'C',
                3: 'G'
            },
        }
        with Project() as proj:
            # expand user provided regions to one or more (chr,start,end,comment)
            regions = expandRegions(raw_regions, proj)
            #
            # all loci contrains all pos in the region and their indexes
            all_loci = defaultdict(set)
            nLoci = 0
            for reg in regions:
                all_loci[reg[0]] = all_loci[reg[0]].union(
                    list(range(reg[1], reg[2] + 1)))
                nLoci += len(all_loci[reg[0]])
            # find the number of loci on each chromosome
            # chroms = sorted(all_loci.keys())
            #
            ref = RefGenome(proj.build)
            genes = genesInRegions(regions, proj)
            env.logger.info(
                '{} genes are identified in the simulated region.'.format(
                    len(genes)))
            # there can be multiple genes (isoforms) in the same regions
            for gene in genes:
                stru = dissectGene(gene, proj)
                pos = []
                seq = ''
                # find the coding regions (xxxx regions) of each gene
                for reg in stru['coding']:
                    # get reference sequence and positions
                    seq += ref.getSequence(reg[0], reg[1], reg[2])
                    pos.extend(list(range(reg[1], reg[2] + 1)))
                ch = reg[0]
                # now, try to divide coding regions by codon. Note that a codon
                # can consist of nucleotie across two exon.
                within = 0
                skip_codon = True
                for idx, (p, s) in enumerate(zip(pos, seq)):
                    if idx % 3 == 0:
                        # the complete codon must be in the simulated region. We do not handle
                        # partial codon because we cannot control mutations outside of the specified
                        # region
                        try:
                            if (p not in all_loci[ch]) or (
                                    pos[idx + 1] not in all_loci[ch]) or (
                                        pos[idx + 2] not in all_loci[ch]):
                                skip_codon = True
                                continue
                        except:
                            skip_codon = True
                            continue
                        # information about the codon: p0, p1, p2, aa, strand
                        codon = ((ch, p), (ch, pos[idx + 1]), (ch,
                                                               pos[idx + 2]),
                                 codon_table[s + seq[idx + 1] + seq[idx + 2]]
                                 if stru['strand'] == '+' else
                                 codon_table_reverse_complement[s +
                                                                seq[idx + 1] +
                                                                seq[idx + 2]],
                                 stru['strand'])
                        skip_codon = False
                    if skip_codon:
                        continue
                    within += 1
                    # record reference sequence
                    self.coding_base[(ch, p)] = s
                    # other two positions share the same codon
                    # because of there can be multiple genes, one basepair can be in multiple codon
                    if (ch, p) in self.codon_info:
                        if codon not in self.codon_info[(ch, p)]:
                            self.codon_info[(ch, p)].append(codon)
                    else:
                        self.codon_info[(ch, p)] = [codon]
                env.logger.info(
                    'Length of all coding regions of {}: {} ({} within specified regions)'
                    .format(gene, len(pos), within))
        #
        # remove all positions that are not in regions
        env.logger.info(
            '{} out of {} bp ({:.2f}%) are in coding regions of genes {}'
            .format(
                len(self.codon_info), nLoci,
                100. * len(self.codon_info) / nLoci, ', '.join(genes)))
        #with open('coding.txt', 'w') as coding:
        #    coding.write(''.join([str(x)+'\n' for x in sorted(self.codon_info.keys())]))

    def _updatePosWithIndexes(self, pop):
        # The class saves loci information as (chr, pos), which is very slow to access
        # when a true population is sent, we are translating them to idx.
        #
        all_pos = sorted(self.coding_base.keys())
        indexes = pop.indexesOfLoci(all_pos)
        posMap = {x: y for x, y in zip(all_pos, indexes)}
        self._coding_base = {posMap[x]: y for x, y in self.coding_base.items()}
        self._codon_info = {}
        for key, codons in self.codon_info.items():
            new_codons = []
            for codon in codons:
                new_codons.append((posMap[codon[0]], posMap[codon[1]],
                                   posMap[codon[2]], codon[3], codon[4]))
            self._codon_info[posMap[key]] = new_codons

    def findAlteredAminoAcid(self, mut, pop):
        # geno is arranged locus by locus (A1,A2,B1,B2 etc)
        #
        # we can not divide the sequence into triplets because it is possible that a nucleotide
        # belong to multiple codon with different locations.
        #
        # the same aa change can be caused by two mutations, we need to keep only one
        # aa_change with the same (aa, naa, ploidy, p0)
        if self._codon_info is None:
            self._updatePosWithIndexes(pop)
        aa_change = set()
        N = pop.totNumLoci()
        for m in mut.keys():
            p = m / N  # ploidy
            loc = m % N
            aa_change_at_m = []
            # we need to map index in coding region to their original position
            for p0, p1, p2, aa, s in self._codon_info[loc]:
                # p0: location
                # mut: mutation (0 for wildtype)
                codon = self.mutant_map[self._coding_base[p0]][mut[p0+p*N]] + \
                        self.mutant_map[self._coding_base[p1]][mut[p1+p*N]] + \
                        self.mutant_map[self._coding_base[p2]][mut[p2+p*N]]
                naa = codon_table[
                    codon] if s == '+' else codon_table_reverse_complement[codon]
                # if this is a real change
                if naa != aa:
                    aa_change_at_m.append((aa, naa, p, p0))
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
                else:
                    stoploss = [x for x in aa_change_at_m if x[0] == '*']
                    if stoploss:
                        aa_change.add(stoploss[0])
                    else:
                        aa_change.add(aa_change_at_m[0])
        return aa_change

    def mutantType(self, mut, pop):
        if self._codon_info is None:
            self._updatePosWithIndexes(pop)
        for k in mut.keys():
            if k not in self._coding_base:
                return 'non-coding'
            else:
                break
        aa_change = self.findAlteredAminoAcid(mut, pop)
        if not aa_change:
            return 'synonymous'
        res = []
        for aa, naa, ploidy, start in aa_change:
            # stoploss
            if aa == '*':
                res.append('stoploss')
            elif naa == '*':
                res.append('stopgain')
            else:
                res.append('missense')
        return ','.join(res)


class ProteinSelector(sim.PySelector, MutantInfo):

    def __init__(self,
                 regions,
                 s_missense=0.001,
                 s_stoploss=0.002,
                 s_stopgain=0.01):
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
        MutantInfo.__init__(self, regions)
        if (not self.codon_info) or (s_missense == 0 and s_stoploss == 0 and
                                     s_stopgain == 0):
            env.logger.warning(
                'Specified region does not contain any gene or all selection coefficient is zero. A neutral model will be used.'
            )
            sim.PySelector.__init__(self, func=self._neutral, loci=[])
        else:
            # we need to send simuPOP to indexes of coding loci within the specified region
            sim.PySelector.__init__(
                self,
                func=self._select,
                loci=[(str(x[0]), x[1]) for x in sorted(self.codon_info.keys())
                     ])

    def _neutral(self):
        return 1

    def _select(self, mut, pop):
        #
        # we assume a multiplicative model
        aa_change = self.findAlteredAminoAcid(mut, pop)
        fitness = 1
        for aa, naa, ploidy, start in aa_change:
            # stoploss
            if aa == '*':
                fitness *= 1 - self.s_stoploss
            elif naa == '*':
                fitness *= 1 - self.s_stopgain
            else:
                fitness *= 1 - self.s_missense
        return fitness


class ProteinPenetrance(sim.PyPenetrance, MutantInfo):

    def __init__(self,
                 regions,
                 s_sporadic=0.0001,
                 s_missense=0.001,
                 s_stoploss=0.002,
                 s_stopgain=0.01):
        '''A protein penetrance model that is identical to ProteinSelector, but
            use 1 minus calculated fitness value as pentrance probability.
        '''
        self.s_sporadic = s_sporadic
        self.s_missense = s_missense
        self.s_stoploss = s_stoploss
        self.s_stopgain = s_stopgain
        #
        MutantInfo.__init__(self, regions)
        if not self.codon_info:
            env.logger.warning(
                'Specified region does not contain any gene. A neutral model will be used.'
            )
            sim.PyPenetrance.__init__(self, func=self._neutral, loci=[])
        else:
            # we need to send simuPOP to indexes of coding loci within the specified region
            sim.PyPenetrance.__init__(
                self,
                func=self._penetrance,
                loci=sorted(self.codon_info.keys()))

    def _neutral(self):
        return 1 - self.s_sporadic

    def _penetrance(self, mut, pop):
        #
        aa_change = self.findAlteredAminoAcid(mut, pop)
        # we assume a multiplicative model
        fitness = 1 - self.s_sporadic
        for aa, naa, ploidy, start in aa_change:
            # stoploss
            if aa == '*':
                fitness *= 1 - self.s_stoploss
            elif naa == '*':
                fitness *= 1 - self.s_stopgain
            else:
                fitness *= 1 - self.s_missense
        return 1 - fitness


class ExtractVCF(PipelineAction):
    '''Extract variants and genotypes in specified regions of a local or online VCF
    file and save to a local VCF file.'''

    def __init__(self, regions, sourceURL, output):
        '''
        Parameters:
            regions (string):
                One or more chromosome regions in the format of chr:start-end
                (e.g. chr21:33,031,597-33,041,570), Field:Value from a region-based
                annotation database (e.g. refGene.name2:TRIM2 or refGene_exon.name:NM_000947),
                or set options of several regions (&, |, -, and ^ for intersection,
                union, difference, and symmetric difference).

            sourceURL (string or list of strings):
                URL or filename of the VCF file from which variants and genotypes
                will be extracted. A list of URLs can be specified to extract
                genotypes from multiple source files.

            output (string):
                Output vcf file.
        '''
        if isinstance(sourceURL, (list, tuple)):
            self.sourceURL = sourceURL
        else:
            self.sourceURL = [sourceURL]
        self.regions = regions
        PipelineAction.__init__(
            self,
            cmd='ExtractVCF {} {}'.format(sourceURL, regions),
            output=output)

    def _execute(self, ifiles, pipeline):
        '''
        Parameters:
            ifiles: unused
            pipeline: unused

        Results:
            Return the resulting vcf file as output.
        '''
        with Project(
                mode=['ALLOW_NO_PROJ', 'READ_ONLY'],
                verbosity=pipeline.verbosity) as proj:
            proj.setRefGenome('hg19')
            tabixFetch(self.sourceURL[0], [], self.output[0], True)
            for r in expandRegions(self.regions):
                region = '{}:{}-{}'.format(r[0], r[1], r[2])
                env.logger.info('Retriving genotype for region chr{}{}'.format(
                    region, ' ({})'.format(r[3] if r[3] else '')))
                for URL in self.sourceURL:
                    tabixFetch(URL, [region], self.output[0], False)
        return True


class CreatePopulation(PipelineAction):
    '''Create a simuPOP population from specified regions and
    number of individuals.
    '''

    def __init__(self,
                 regions,
                 size=None,
                 importGenotypeFrom=None,
                 infoFields=[],
                 build=None,
                 output=[],
                 **kwargs):
        '''
        Parameters:
            regions: (string):
                One or more chromosome regions in the format of chr:start-end
                (e.g. chr21:33,031,597-33,041,570), Field:Value from a region-based
                annotation database (e.g. refGene.name2:TRIM2 or refGene_exon.name:NM_000947),
                or set options of several regions (&, |, -, and ^ for intersection,
                union, difference, and symmetric difference).

            size (None, integer, list of integers):
                Size of the population. This parameter can be ignored (``None``) if
                parameter ``importGenotypeFrom`` is specified to import genotypes from
                external files.

            importGenotypeFrom (None or string):
                A file from which genotypes are imported. Currently a file with extension
                ``.vcf`` or ``.vcf.gz`` (VCF format) or ``.ms`` (MS format) is supported.
                Because the ms format does not have explicit location of loci, the loci are
                spread over the specified regions according to their relative locations.

            infoFields (string or list of strings):
                information fields of the population, if needed by particular operators
                during evolution.

            build (string):
                build of the reference genome. Default to hg19.

            output (string):
                Name of the created population in simuPOP's binary format.

            kwargs (arbitrary keyword parameters):
                Additional parameters that will be passed to the constructor of
                population (e.g. ``ploidy=1`` for haploid population). Please refer
                to the ``Population()`` function of simuPOP for details.
        '''
        self.regions = regions
        self.size = size
        self.infoFields = infoFields
        self.sourceFile = importGenotypeFrom
        self.extra_args = kwargs
        # this parameter is currently not used. Need to get it working later
        self.build = build
        PipelineAction.__init__(
            self,
            cmd='CreatePopulation {} {}\n'.format(regions, size),
            output=output)

    def _execute(self, ifiles, pipeline):
        # translate regions to simuPOP ...
        lociPos = {}
        regions = expandRegions(self.regions)
        for r in regions:
            if r[0] in lociPos:
                lociPos[r[0]].extend(list(range(r[1], r[2] + 1)))
            else:
                lociPos[r[0]] = list(range(r[1], r[2] + 1))
        #
        if self.sourceFile is not None:
            if self.sourceFile.lower().endswith(
                    '.vcf') or self.sourceFile.lower().endswith('.vcf.gz'):
                pop = self._importFromVcf(lociPos)
            elif self.sourceFile.lower().endswith('.ms'):
                pop = self._importFromMS(lociPos)
            else:
                raise ValueError(
                    'CreatePopulaton can only import genotypes from files in '
                    'vcf (with extension .vcf and .vcf.gz) or ms (with extension .ms) formats'
                )
            if self.size is not None and pop.popSize() != self.size:
                env.logger.warning(
                    'Population imported from {} has {} individuals where a '
                    'population of size {} is requested'.format(
                        self.sourceFile, pop.popSize(), self.size))
        else:
            chroms = list(lociPos.keys())
            chroms.sort()
            pop = sim.Population(
                size=self.size,
                loci=[len(lociPos[x]) for x in chroms],
                chromNames=chroms,
                lociPos=sum([lociPos[x] for x in chroms], []),
                infoFields=self.infoFields,
                **self.extra_args)
            #chromTypes=[sim.CHROMOSOME_X if x=='X' else (sim.CHROMOSOME_Y if x=='Y' else sim.AUTOSOME) for x in chroms])
        # save regions for later use.
        pop.dvars().regions = self.regions
        env.logger.info('Saving created population to {}'.format(
            self.output[0]))
        pop.save(self.output[0])
        return True

    def _importFromVcf(self, lociPos):
        chroms = list(lociPos.keys())
        chroms.sort()
        #
        mutantCount = 0
        # create a dictionary of lociPos->index on each chromosome
        lociIndex = {}
        for chIdx, ch in enumerate(chroms):
            lociIndex[chIdx] = {pos: idx for idx, pos in enumerate(lociPos[ch])}
        #
        pop = None
        with open(self.sourceFile, 'r') as vcf:
            for line in vcf:
                if line.startswith('#'):
                    continue
                fields = line.split('\t')
                if pop is None:
                    if len(fields) <= 10:
                        raise ValueError(
                            'Input vcf file does not contain any genotype information'
                        )
                    pop = sim.Population(
                        size=len(fields) - 10,
                        loci=[len(lociPos[x]) for x in chroms],
                        chromNames=chroms,
                        lociPos=sum([lociPos[x] for x in chroms], []),
                        infoFields=self.infoFields)
                #
                chr = pop.chromNames().index(fields[0])
                try:
                    # due to a bug in tabix, loci position might be outside of the specified region
                    # so there is no "pos" for this line of genotype, but that is ok anyway. (issue #62)
                    pos = lociIndex[chr][int(fields[1])]
                except:
                    continue
                for ind, geno in enumerate(fields[10:]):
                    if geno[0] not in ('0', '.'):
                        pop.individual(ind).setAllele(1, pos, 0, chr)
                        mutantCount += 1
                    if geno[2] not in ('0', '.'):
                        pop.individual(ind).setAllele(1, pos, 1, chr)
                        mutantCount += 1
        return pop

    def _importFromMS(self, lociPos):
        #
        with open(self.sourceFile, 'r') as ms:
            _ = ms.readline()
            _ = ms.readline()
            ms.readline()
            # current we do not handle multiple populations etc
            chroms = list(lociPos.keys())
            chroms.sort()
            #
            all_indexes = []
            all_geno = []
            indexes = []
            geno = []
            for line in ms:
                if line.strip() in ['' or '//']:
                    continue
                elif line.startswith('segsites:'):
                    segsites = int(line[len('segsites:'):])
                    # new block?
                    if indexes:
                        all_indexes.append(indexes)
                        all_geno.append(geno)
                        indexes = []
                        geno = []
                elif line.startswith('positions:'):
                    positions = [
                        float(x) for x in line[len('positions:'):].split()
                    ]
                    if len(positions) != segsites:
                        raise ValueError(
                            'Number of segsites do not match number of positions'
                        )
                    #
                    # check number of loci
                    nLoci = len(lociPos[chroms[len(all_indexes)]])
                    if nLoci < segsites:
                        raise ValueError(
                            'Specified region cannot accomendate {} segregating sites'
                            .format(segsites))
                    # the last index is not allowed
                    indexes = [int(nLoci * x) for x in positions if x != 1.0]
                    if len(set(indexes)) != len(indexes):
                        env.logger.warning(
                            '{} loci at identical location needs to be re-located'
                            .format(len(indexes) - len(set(indexes))))
                        existing_indexes = set(indexes)
                        acceptable_indexes = list(
                            set(range(nLoci)) - existing_indexes)
                        random.shuffle(acceptable_indexes)
                        indexes = list(
                            existing_indexes
                            | set(acceptable_indexes[:segsites -
                                                     len(existing_indexes)]))
                        indexes.sort()
                    # read?
                    # population size?
                else:
                    geno.append(line.strip())
            # add everything to all_index etc
            all_indexes.append(indexes)
            del indexes
            all_geno.append(geno)
            del geno
            #
            pop = sim.Population(
                size=[len(x) / 2 for x in all_geno],
                loci=[len(lociPos[x]) for x in chroms],
                chromNames=chroms,
                lociPos=sum([lociPos[x] for x in chroms], []),
                infoFields=self.infoFields)
            # add genotype
            prog = ProgressBar('Importing from {}'.format(self.sourceFile),
                               pop.popSize() * len(chroms))
            processed = 0
            for ch in range(len(chroms)):
                index = all_indexes[ch]
                for ind, geno in enumerate(all_geno[ch]):
                    for idx, g in enumerate(geno):
                        if g != '0':
                            pop.individual(ind / 2).setAllele(
                                int(g), index[idx], ind % 2, ch)
                    processed += 1
                    prog.update(processed)
            prog.done()
            return pop


class EvolvePopulation(PipelineAction):
    '''Evolve a population, subject to passed operators.'''

    def __init__(self,
                 selector=None,
                 mutator=None,
                 transmitter=None,
                 initOps=[],
                 preOps=[],
                 taggers=[],
                 postOps=[],
                 finalOps=[],
                 demoModel=None,
                 output=[]):
        '''
        Parameters:
            selector (None, one or a list of simuPOP operators):
                Operator for natural selection that will be applied before mating
                at each generation. Default to no natural selection.

            mutator (None, or one or a list of simuPOP operators):
                Mutation operator that will be applied before mating at each
                generation. Default to no mutation.

            transmitter (None, or one or a list of simuPOP operators):
                Genotype transmision operators. Default to ``MedelianGenoTransmitter``.

            initOps, preOps, taggers, postOps, finalOps (one or a list of simuPOP operators):
                Additional operaors that will be applied during initiation, before, during
                and after mating, and after the complection of evolutionary process.

            demoModel (A demographic model):
                A demographic model as defined in the simuPOP.demography module.

            output (string):
                Evolved population in simuPOP binary format.
        '''
        self.mutator = [] if mutator is None else (
            mutator if isinstance(mutator, (list, tuple)) else [mutator])
        self.selector = [] if selector is None else (
            selector if isinstance(selector, (list, tuple)) else [selector])
        self.demoModel = demoModel
        self.transmitter = [
            sim.MendelianGenoTransmitter()
        ] if transmitter is None else (
            transmitter if isinstance(transmitter,
                                      (list, tuple)) else [transmitter])
        self.output = [output]
        self.taggers = taggers if isinstance(taggers,
                                             (list, tuple)) else [taggers]
        self.initOps = initOps if isinstance(initOps,
                                             (list, tuple)) else [initOps]
        self.preOps = preOps if isinstance(preOps, (list, tuple)) else [preOps]
        self.postOps = postOps if isinstance(postOps,
                                             (list, tuple)) else [postOps]
        self.finalOps = finalOps if isinstance(finalOps,
                                               (list, tuple)) else [finalOps]
        PipelineAction.__init__(
            self, cmd='EvolvePop output={}\n'.format(output), output=output)

    def _execute(self, ifiles, pipeline):
        '''
        Parameters:
            ifiles (list of strings):
                A population file in simuPOP's binary format.

            pipeline: unused.

        Results:
            Evolved population
        '''
        pop = sim.loadPopulation(ifiles[0])
        pop.addInfoFields(['fitness', 'migrate_to'])
        startTime = time.clock()
        # reset RNG to use pipeline seed
        sim.setRNG(seed=int(pipeline.VARS['seed']))
        #
        # Evolve
        env.logger.info('Start evolving...')
        pop.dvars().last_time = time.time()
        pop.evolve(
            initOps=[sim.InitSex(), sim.PyExec('import time')] + self.initOps,
            preOps=[
                sim.PyOutput(
                    '''Statistics outputted are
    1. Generation number,
    2. population size (a list),
    3. number of segregation sites,
    4. average number of mutants per individual
    5. average allele frequency * 100
    6. average fitness value
    7. minimal fitness value of the parental population
    ''',
                    at=0),
                # revert alleles at fixed loci to wildtype
                sim.RevertFixedSites()
            ] +
            #
            # 'A' is zero, no need to map in and out
            self.mutator +
            #sim.TicToc(),
            self.selector + [
                sim.IfElse(
                    'time.time() - last_time > 30',
                    [
                        sim.PyExec('last_time = time.time()'),
                        # output statistics in verbose mode
                        sim.Stat(
                            popSize=True,
                            meanOfInfo='fitness',
                            minOfInfo='fitness',
                            numOfSegSites=sim.ALL_AVAIL,
                            numOfMutants=sim.ALL_AVAIL),
                        sim.PyEval(
                            r'"%5d %s %5d %.6f %.6f %.6f %.6f" '
                            '% (gen, subPopSize, numOfSegSites, float(numOfMutants)/popSize, '
                            '(numOfMutants * 50. / numOfSegSites / popSize) if numOfSegSites else 0, '
                            ' meanOfInfo["fitness"], minOfInfo["fitness"])',
                            output=env.logger.info),
                    ])
            ] + self.preOps,
            matingScheme=sim.RandomMating(
                ops=self.transmitter + self.taggers, subPopSize=self.demoModel),
            postOps=self.postOps,
            finalOps=[
                # revert fixed sites so that the final population does not have fixed sites
                sim.RevertFixedSites()
            ] +
            #
            # apply fitness to get final statistics and facilitate sampling
            self.selector + [
                # statistics after evolution
                sim.Stat(
                    popSize=True,
                    meanOfInfo='fitness',
                    minOfInfo='fitness',
                    numOfSegSites=sim.ALL_AVAIL,
                    numOfMutants=sim.ALL_AVAIL),
                sim.PyEval(
                    r'"%5d %s %5d %.6f %.6f %.6f %.6f\n" '
                    '% (gen, subPopSize, numOfSegSites, float(numOfMutants)/popSize, '
                    '(numOfMutants * 50. / numOfSegSites / popSize) if numOfSegSites else 0, '
                    ' meanOfInfo["fitness"], minOfInfo["fitness"])',
                    output=env.logger.info),
                sim.PyEval(
                    r'"Simulated population has %d individuals, %d segregation sites. '
                    r'There are on average %.1f mutants per individual. Mean allele frequency is %.4f%%.\n"'
                    r'% (popSize, numOfSegSites, numOfMutants / popSize, (numOfMutants * 50. / numOfSegSites/ popSize) if numOfSegSites else 0)',
                    output=env.logger.info),
            ] + self.finalOps,
            gen=int(self.demoModel.num_gens))
        #
        env.logger.info('Population simulation takes %.2f seconds' %
                        (time.clock() - startTime))
        pop.save(self.output[0])
        return True


class DrawCaseControlSample(PipelineAction):
    '''Draw case control samples from simulated population. If there are subpopulations and
    cases and controls have the same demension, cases and controls are drawn from each
    subpopulation. If more than one output files are provided, multiple samples are
    draw from the population.
    '''

    def __init__(self, cases, controls, penetrance, output):
        '''
        Parameters:
            cases (integer):
                Number of affected individuals to sample

            controls (integer):
                Number of unaffected individuals to sample

            penetrance (a simuPOP operator):
                A penetrance operator that determines the affection status
                of individuals.

            output:
                Sample population saved in simuPOP binary format.
        '''
        if isinstance(cases, int):
            self.cases = [cases]
        else:
            self.cases = cases
        if isinstance(controls, int):
            self.controls = [controls]
        else:
            self.controls = controls
        #
        if len(self.cases) != len(self.controls):
            raise ValueError(
                'Please specify cases and controls with the same dimensions.')
        self.penetrance = penetrance
        #
        PipelineAction.__init__(
            self, cmd='DrawCaseCtrlSample {}\n'.format(output), output=output)

    def _execute(self, ifiles, pipeline):
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        for idx, output in enumerate(self.output):
            self.selectedCases = [0] * len(self.cases)
            self.selectedCtrls = [0] * len(self.controls)
            self.parentsOfCases = []
            self.parentsOfCtrls = []
            self.numOfOffspring = 0
            self.numOfAffected = 0
            # if there is only one sample or the last sample, use the loaded population.
            if idx == len(self.output) - 1:
                self._drawSample(pop, output)
            else:
                # else use a copy of the loaded population
                self._drawSample(pop.clone(), output)
        return True

    def _drawSample(self, pop, output):
        #
        env.logger.info(
            'Draw {} cases and {} controls from a population of sizes sizes {}'
            .format(self.cases, self.controls, pop.subPopSizes()))
        if len(self.cases) > 1 and pop.numSubPop() != len(self.cases):
            raise ValueError(
                'If an array of cases and controls are specified, they should match the number of subpopulations.'
            )
        if len(self.cases) == 1 and pop.numSubPop() > 1:
            pop.mergeSubPops()
        self.prog = ProgressBar(
            'Generating %d cases and %d controls' %
            (sum(self.cases), sum(self.controls)),
            sum(self.cases) + sum(self.controls))
        pop.addInfoFields(['father_idx', 'mother_idx'])
        self.parentalPopSize = pop.popSize()
        pop.evolve(
            matingScheme=sim.RandomMating(
                ops=[
                    sim.MendelianGenoTransmitter(),
                    sim.ParentsTagger(),
                    # apply a penetrance model
                    self.penetrance
                ] + [
                    # an individual will be discarded if _selectInds returns False
                    sim.PyOperator(func=self._selectInds, subPops=sp, param=sp)
                    for sp in range(pop.numSubPop())
                ],
                subPopSize=[x + y for x, y in zip(self.cases, self.controls)],
            ),
            gen=1)
        self.prog.done()
        #
        for field in ['migrate_to', 'fitness', 'father_idx', 'mother_idx']:
            if field in pop.infoFields():
                pop.removeInfoFields(field)
        uParents = len(set(self.parentsOfCases))
        env.logger.info(
            '{} cases are produced from {} ({:.2f}% of {}) unique parents'
            .format(
                sum(self.selectedCases), uParents,
                uParents * 100. / self.parentalPopSize, self.parentalPopSize))
        uParents = len(set(self.parentsOfCtrls))
        env.logger.info(
            '{} controls are produced from {} ({:.2f}% of {}) unique parents'
            .format(
                sum(self.selectedCtrls), uParents,
                uParents * 100. / self.parentalPopSize, self.parentalPopSize))
        env.logger.info(
            'Observed prevalence of the disease is {:.3f}% ({}/{})'.format(
                self.numOfAffected * 100. / self.numOfOffspring,
                self.numOfAffected, self.numOfOffspring))
        env.logger.info('Saving samples to population {}'.format(output))
        pop.save(output)
        del pop

    def _selectInds(self, off, param):
        'Determine if the offspring can be kept.'
        self.numOfOffspring += 1
        if off.affected():
            self.numOfAffected += 1
            if self.selectedCases[param] < self.cases[param]:
                self.selectedCases[param] += 1
                self.parentsOfCases.extend([off.father_idx, off.mother_idx])
                self.prog.update(
                    sum(self.selectedCases) + sum(self.selectedCtrls))
                return True
        else:
            if self.selectedCtrls[param] < self.controls[param]:
                self.selectedCtrls[param] += 1
                self.parentsOfCtrls.extend([off.father_idx, off.mother_idx])
                self.prog.update(
                    sum(self.selectedCases) + sum(self.selectedCtrls))
                return True
        return False


class DrawRandomSample(PipelineAction):
    '''Draw random sample from an input population and save
    samples (population) in simuPOP's binary format. '''

    def __init__(self, sizes, output):
        '''
        Parameters:
            sizes (integer):
                Number of individuals to sample

            output:
                Sample population saved in simuPOP binary format.
        '''
        self.sizes = sizes
        PipelineAction.__init__(
            self, cmd='DrawRandomSampler {}\n'.format(output), output=output)

    def _execute(self, ifiles, pipeline):
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        for idx, output in enumerate(self.output):
            # if there is only one sample or the last sample, use the loaded population.
            if idx == len(self.output) - 1:
                self._drawSample(pop, output)
            else:
                # else use a copy of the loaded population
                self._drawSample(pop.clone(), output)
        return True

    def _drawSample(self, pop, output):
        pop.evolve(matingScheme=sim.RandomMating(subPopSize=self.sizes), gen=1)
        #
        for field in ['migrate_to', 'fitness', 'father_idx', 'mother_idx']:
            if field in pop.infoFields():
                pop.removeInfoFields(field)
        env.logger.info('Saving samples to population {}'.format(output))
        pop.save(output)
        del pop


class DrawQuanTraitSample(PipelineAction):
    '''Draw samples according to individual trait values.
    '''

    def __init__(self, sizes, conditions, qtrait, output):
        '''
        Parameters:
            sizes (integer):
                Number of individuals to sample.

            conditions (string):
                Conditions for selecting sample, which is a Python
                expression that involves name of the quantitative
                trait.

            qtrait (a simuPOP operator)
                A simuPOP quantitative trait operator to assign
                quantitative traits to individuals.

            output (string):
                Samples (population) saved in simuPOP binary format.
        '''
        if isinstance(sizes, int):
            self.sizes = [sizes]
        else:
            self.sizes = sizes
        #
        if isinstance(conditions, str):
            self.conditions = [conditions]
        else:
            self.conditions = conditions
        #
        self.qtrait = qtrait
        PipelineAction.__init__(
            self, cmd='DrawQuanTraitSample {}\n'.format(output), output=output)

    def _execute(self, ifiles, pipeline):
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        for idx, output in enumerate(self.output):
            # if there is only one sample or the last sample, use the loaded population.
            if idx == len(self.output) - 1:
                self._drawSample(pop, output)
            else:
                # else use a copy of the loaded population
                self._drawSample(pop.clone(), output)
        return True

    def _drawSample(self, pop, output):
        if pop.numSubPop() > 1:
            pop.mergeSubPops()
        #
        if len(self.sizes) > 1:
            pop.resize(0, pop.popSize() * len(self.sizes), propagate=True)
            pop.splitSubPop(0, [pop.popSize()] * len(self.sizes()))
        pop.evolve(
            matingScheme=sim.RandomMating(
                ops=[
                    sim.MendelianGenoTransmitter(),
                    # apply a penetrance model
                    self.qtrait
                ] + [
                    # an individual will be discarded if _selectInds returns False
                    sim.DiscardIf('not (' + cond + ')', subPops=sp)
                    for sp, cond in self.condistions
                ],
                subPopSize=self.sizes),
            gen=1)
        #
        for field in ['migrate_to', 'fitness', 'father_idx', 'mother_idx']:
            if field in pop.infoFields():
                pop.removeInfoFields(field)
        env.logger.info('Saving samples to population {}'.format(output))
        pop.save(output)


class ExportPopulation(PipelineAction):
    '''Export simulated populatons in vcf, and fasta format.
    '''

    def __init__(self, output, sample_names=[], var_info=[]):
        '''
        Parameters:
            output (string):
                Name of the output file. The variants and genotypes will be saved in ``VCF``
                format is the filename ends with ``.vcf``. Otherwise the sequences will be
                saved in fasta format.

            sample_names (list of strings):
                Name of samples. Default to S_i for i=1, 2, ...


            var_info (string or list of strings):
                Variant information fields that will be outputted if exporting in vcf format.
        '''
        self.sample_names = sample_names
        self.var_info = var_info
        for vi in var_info:
            if vi != 'MT':
                env.logger.warning(
                    'Unrecognized variant info field. Only MT (mutant type) is supported now.'
                )
        PipelineAction.__init__(
            self,
            cmd='ExportPopulation {} {}\n'.format(sample_names, output),
            output=output)

    def _execute(self, ifiles, pipeline):
        '''
        Parameters:
            ifiles (string):
                Population in simuPOP's binary format.

            pipeline: unused
        '''
        if self.output[0].endswith('.vcf'):
            self._exportVCF(ifiles, pipeline)
        else:
            self._exportFasta(ifiles, pipeline)
        #
        if len(self.output) == 1:
            return
        # output phenotype
        pop = sim.loadPopulation(ifiles[0])

        with open(self.output[1], 'w') as phe:
            fields = pop.infoFields()
            phe.write('sample_name\taff' + '\t'.join(fields) + '\n')
            prog = ProgressBar(
                'Exporting phenotype {}'.format(','.join(fields)),
                pop.popSize())
            for i in range(pop.popSize()):
                if self.sample_names:
                    phe.write(self.sample_names[i])
                else:
                    phe.write('S_{}'.format(i + 1))
                #
                if pop.individual(i).affected():
                    phe.write('\t2')
                else:
                    phe.write('\t1')
                for info in fields:
                    phe.write('\t{}'.format(pop.individual(i).info(info)))
                phe.write('\n')
                prog.update(i + 1)
            prog.done()
        return True

    def _exportFasta(self, ifiles, pipeline):
        # import variants to the current project
        # translate 0,1,2,3 to A,C,G,T
        alleleMap = {
            'A': {
                0: 'A',
                1: 'C',
                2: 'G',
                3: 'T'
            },
            'C': {
                0: 'C',
                1: 'G',
                2: 'T',
                3: 'A'
            },
            'G': {
                0: 'G',
                1: 'T',
                2: 'A',
                3: 'C'
            },
            'T': {
                0: 'T',
                1: 'A',
                2: 'C',
                3: 'G'
            },
            'N': {
                0: 'A',
                1: 'C',
                2: 'G',
                3: 'T'
            },
        }
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        # if we need to output variant information
        #aaInfo = None
        if self.var_info:
            if 'regions' not in pop.vars():
                env.logger.warning(
                    'No mutation type can be obtained because population object does not have variable "regions".'
                )
            #else:
                # figure out regions of the population
                #aaInfo = MutantInfo(pop.dvars().regions)
                #pass
        # output genotype
        with Project(
                mode=['ALLOW_NO_PROJ', 'READ_ONLY'],
                verbosity=pipeline.verbosity) as proj:
            build = proj.build
        if build is None:
            build = 'hg19'
        with open(self.output[0], 'w') as fasta:
            #
            # get reference genome
            refGenome = RefGenome(build)
            # get reference sequence
            refSeq = []
            for chr in range(pop.numChrom()):
                chr_name = pop.chromName(chr)
                refSeq.extend([
                    refGenome.getBase(chr_name, int(pop.locusPos(loc)))
                    for loc in range(pop.chromBegin(chr), pop.chromEnd(chr))
                ])
            #
            prog = ProgressBar('Exporting simulated population', pop.popSize())
            for idx, ind in enumerate(pop.individuals()):
                for p in range(2):
                    geno = [
                        alleleMap[ref][allele]
                        for ref, allele in zip(refSeq, ind.genotype(p))
                    ]
                    fasta.write(''.join(geno) + '\n')
                prog.update(idx + 1)
            prog.done()

    def _exportVCF(self, ifiles, pipeline):
        # import variants to the current project
        # translate 0,1,2,3 to A,C,G,T
        alleleMap = {
            'A': {
                0: 'A',
                1: 'C',
                2: 'G',
                3: 'T'
            },
            'C': {
                0: 'C',
                1: 'G',
                2: 'T',
                3: 'A'
            },
            'G': {
                0: 'G',
                1: 'T',
                2: 'A',
                3: 'C'
            },
            'T': {
                0: 'T',
                1: 'A',
                2: 'C',
                3: 'G'
            },
            'N': {
                0: 'A',
                1: 'C',
                2: 'G',
                3: 'T'
            },
        }
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        # if we need to output variant information
        aaInfo = None
        if self.var_info:
            if 'regions' not in pop.vars():
                env.logger.warning(
                    'No mutation type can be obtained because population object does not have variable "regions".'
                )
            else:
                # figure out regions of the population
                aaInfo = MutantInfo(pop.dvars().regions)
        # output genotype
        with Project(
                mode=['ALLOW_NO_PROJ', 'READ_ONLY'],
                verbosity=pipeline.verbosity) as proj:
            build = proj.build
        if build is None:
            build = 'hg19'
        with open(self.output[0], 'w') as vcf:
            vcf.write('##fileformat=VCFv4.1\n')
            vcf.write('##fileData={}\n'.format(
                datetime.date.today().strftime("%Y%m%d")))
            vcf.write('##source=Variant Simulation Tools\n')
            vcf.write('##reference={}\n'.format(build))
            if self.var_info:
                for vi in self.var_info:
                    if vi == 'MT':
                        vcf.write(
                            '#INFO=<ID=MT,Number=1,Type=String,Description="Mutant type">\n'
                        )
            vcf.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            )
            vcf.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t')
            if self.sample_names:
                if len(self.sample_names) != pop.popSize():
                    raise ValueError(
                        'Sample names, if specified, should be assigned to'
                        'all {} individuals.'.format(pop.popSize()))
                vcf.write('\t'.join(self.sample_names) + '\n')
            else:
                vcf.write('\t'.join(
                    ['S_{}'.format(x)
                     for x in range(1,
                                    pop.popSize() + 1)]) + '\n')
            #
            # get reference genome
            refGenome = RefGenome(build)
            sim.stat(pop, alleleFreq=sim.ALL_AVAIL, vars='alleleNum')
            nAlleles = 2 * pop.popSize()
            segregated = [
                loc for loc, nums in pop.dvars().alleleNum.items()
                if nums[0] == nAlleles
            ]
            env.logger.info('Genetic variants identified on {} loci'.format(
                len(segregated)))
            prog = ProgressBar(self.output[0], pop.totNumLoci())
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
                    if alt and aaInfo is not None:
                        # output mutant info, always assume first homologus copy
                        varInfo = 'MT=' + aaInfo.mutantType(
                            sim.defdict({loc: x for x in alt}), pop)
                    else:
                        varInfo = '.'
                    if len(alt) == 0:
                        env.logger.error(
                            'No alternative allele at locus {}:{}'.format(
                                chr_name, pos))
                        continue
                    elif len(alt) == 1:
                        # easier...
                        vcf.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t'.format(
                            chr_name, pos, ref, alleleMap[ref][alt[0]],
                            varInfo))
                        vcf.write('\t'.join([
                            '{}/{}'.format(x if x == 0 else 1, y if y ==
                                           0 else 1)
                            for x, y in zip(geno1, geno2)
                        ]) + '\n')
                    else:
                        # we need to figure out the index of geno in alt
                        vcf.write('{}\t{}\t.\t{}\t{}\t.\tPASS\t{}\tGT\t'.format(
                            chr_name, pos, ref,
                            ','.join([alleleMap[ref][a] for a in alt]),
                            varInfo))
                        vcf.write('\t'.join([
                            '{}/{}'.format(x if x == 0 else alt.index(x) +
                                           1, y if y == 0 else alt.index(y) + 1)
                            for x, y in zip(geno1, geno2)
                        ]) + '\n')
                    #
                    prog.update(loc)
            prog.done()
