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
simuOpt.setOptions(alleleType='mutant', optimized=True, quiet=True, version='1.0.5')

import simuPOP as sim
from simuPOP.demography import *
from simuPOP.sampling import *

from .utils import env, expandRegions, ProgressBar, RefGenome, existAndNewerThan, \
    calculateMD5, genesInRegions, codon_table, codon_table_reverse_complement, \
    dissectGene, downloadFile, decompressGzFile

import os, sys, math, time, random
from collections import defaultdict
from .pipeline import SkiptableAction
from .project import Project

if sys.version_info.major == 2:
    from ucsctools_py2 import tabixFetch
else:
    from ucsctools_py3 import tabixFetch

Recom_URL = 'ftp://ftp.hapmap.org/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz'

def set_map_dist(pop, recRate=1e-8):
    '''Set map distance for each locus'''
    recom = downloadFile(Recom_URL)
    recom_files = decompressGzFile(recom)
    recom_map = {}
    for ch in range(pop.numChrom()):
        ch_name = pop.chromName(ch)
        l_pos = int(pop.locusPos(pop.chromBegin(ch)))
        h_pos = int(pop.locusPos(pop.chromEnd(ch)-1) + 1)
        try:
            #env.logger.error([os.path.basename(x) for x in recom_files])
            #env.logger.error('genetic_map_GRCh37_chr{}.txt'.format(ch_name))
            recom_file = [x for x in recom_files if os.path.basename(x) == 'genetic_map_GRCh37_chr{}.txt'.format(ch_name)][0]
        except Exception as e:
            raise RuntimeError('Could not find recombination map for chromosome {}: {}'.format(ch_name, e))
        dist = {}
        with open(recom_file) as recom:
            recom.readline()
            for line in recom:
                try:
                    fields = line.split()
                    pos = int(fields[1])
                    if pos >= l_pos and pos <= h_pos:
                        dist[pos] = float(fields[3])
                except Exception as e:
                    env.logger.warning(e)
        env.logger.info('Map distance of %d markers within region %s:%d-%d are found' % (len(dist), ch_name, l_pos, h_pos))
        if not dist:
            # using generic distance
            for idx,loc in enumerate(range(pop.chromBegin(ch), pop.chromEnd(ch))):
                recom_map[loc] = idx*recRate
            continue
        nLoci = pop.numLoci(ch)
        endLoc = pop.chromEnd(ch)
        # now, try to set genetic map
        map_dist = [-1]*endLoc;
        cnt = 0
        for loc in range(pop.chromBegin(ch), pop.chromEnd(ch)):
            pos = int(int(pop.locusPos(loc)))
            try:
                map_dist[loc] = dist[pos]
                cnt += 1
            except:
                pass
        if cnt != nLoci:
            prev = -1     # name of the previous marker with map distance
            next = -1     # name of the next marker with map distance
            loc = 0
            while (loc < endLoc):
                # already has value
                if map_dist[loc] != -1:
                    prev = next = loc
                    loc += 1
                    continue
                # find the ext one
                if prev == next:
                    next = -1
                    for n in range(loc+1, endLoc):
                         if map_dist[n] != -1:
                             next = n
                             next_pos = pop.locusPos(next)
                             next_dis = map_dist[next]
                             break
                    if prev == -1:
                        for n in range(next):
                            # rough estimation: distance (in cM) proportional to 0.01 recombination rate
                            map_dist[n] = map_dist[next] - (pop.locusPos(next) - pop.locusPos(n)) * 1e-8
                    # if not found, this is at the end of a chromosome
                    elif next == -1:
                        for n in range(loc, endLoc):
                            map_dist[n] = map_dist[prev] + (pop.locusPos(n) - pop.locusPos(prev)) * 1e-8
                        break
                    # if found, but no previous, this is the first one
                    else:
                        prev_pos = pop.locusPos(prev)
                        prev_dis = map_dist[prev]
                        for n in range(loc, next):
                            map_dist[n] = map_dist[prev] + (pop.locusPos(n) - prev_pos) / (next_pos - prev_pos) * (next_dis - prev_dis)
                    prev = next
                    loc = next + 1
        for loc in range(pop.chromBegin(ch), pop.chromEnd(ch)):
            recom_map[loc] = map_dist[loc]
        env.logger.info('Total map distance from {}:{}-{} ({} bp) is {:.4e} cM (r={:.4e}/bp)'.format(ch_name,
            l_pos, h_pos, h_pos - l_pos, recom_map[pop.chromEnd(ch)-1] - recom_map[pop.chromBegin(ch)],
               (recom_map[pop.chromEnd(ch)-1] - recom_map[pop.chromBegin(ch)])/pop.numLoci(ch)))
    pop.dvars().geneticMap = recom_map


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
        
class PopFromRegions(SkiptableAction):
    '''Create a simuPOP population from specified regions and number of individuals.
    '''
    def __init__(self, regions, size, output):
        self.regions = regions
        self.size = size
        SkiptableAction.__init__(self, cmd='PopFromRegions {} {}\n'.format(regions, output),
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
        pop = sim.Population(size=self.size, loci=[len(lociPos[x]) for x in chroms],
            chromNames = chroms, lociPos=sum([lociPos[x] for x in chroms], []))
        set_map_dist(pop)
        pop.save(self.output[0])


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
        #
        pop = None
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
                if pop is None:
                    pop = sim.Population(size=len(fields)-10, loci=[len(lociPos[x]) for x in chroms],
                        chromNames = chroms, lociPos=sum([lociPos[x] for x in chroms], []))
                #
                chr = pop.chromNames().index(fields[0])
                pos = lociIndex[chr][int(fields[1])]
                for ind, geno in enumerate(fields[10:]):
                    if geno[0] not in ('0', '.'):
                        pop.individual(ind).setAllele(1, pos, 0, chr)
                        mutantCount += 1
                    if geno[2] not in ('0', '.'):
                        pop.individual(ind).setAllele(1, pos, 1, chr)
                        mutantCount += 1
        env.logger.info('{} mutants imported'.format(mutantCount))
        set_map_dist(pop)
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
            phe.write('sample_name\taff' + '\t'.join(fields) + '\n')
            prog = ProgressBar('Exporting phenotype {}'.format(','.join(fields)), pop.popSize())
            for i in range(pop.popSize()):
                if self.sample_names:
                    phe.write(self.sample_names[i])
                else:
                    phe.write('S_{}'.format(i+1))
                #
                if pop.individual(i).affected():
                    phe.write('\t2')
                else:
                    phe.write('\t1')
                for info in fields:
                    phe.write('\t{}'.format(pop.individual(i).info(info)))
                phe.write('\n')
                prog.update(i+1)
            prog.done()


class OutputPopStat(SkiptableAction):
    def __init__(self, mut_count=None, output=[]):
        self.output = output
        self.mut_count = mut_count
        SkiptableAction.__init__(self, cmd='PopStat output={}\n'
            .format(output), output=output)

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


    def count_mutants(self):
        #
        sz = min(self.pop.popSize(), self.mut_count[1])
        sample = drawRandomSample(self.pop, sz)
        #
        sim.stat(sample, alleleFreq=sim.ALL_AVAIL, vars='alleleNum')
        #
        num = [int(sz*2-x[0]) for x in sample.dvars().alleleNum.values()]
        cnt = [0]*(sz*2 + 1)
        for n in num:
            cnt[n] += 1
        with open(self.mut_count[0], 'w') as out:
            # output number of loci in each frequency class
            out.write('## number of mutant\t#loci in this mutant class\n')
            for idx,c in enumerate(cnt):
                if c != 0:
                    out.write('#{}\t{}\n'.format(idx, c))
            out.write('## index, chr, position, #0, #1, #2, #3\n')
            # output allele counts at each locus
            for ch in range(self.pop.numChrom()):
                ch_name = self.pop.chromName(ch)
                for i in range(self.pop.chromBegin(ch), self.pop.chromEnd(ch)):
                    if int(sample.dvars().alleleNum[i][0]) != sz*2:
                        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(i+1, 
                            ch_name,
                            int(self.pop.locusPos(i)),
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
        seg_sites = set([x for x in range(pop.totNumLoci()) if len(pop.dvars().alleleFreq[x]) == 1])
        env.logger.info('LD from {} segregating sites'.format(len(seg_sites)))
        numBin = maxDist/binDist
        result = {'LD_prime': {}, 'R2': {}}
        for sp in pop.subPopNames():
            result['LD_prime'][sp] = []
            result['R2'][sp] = []
        pairs = [[] for dist in range(numBin)]
        for ch in range(pop.numChrom()):
            seg_sites_ch = [x for x in range(pop.chromBegin(ch), pop.chromEnd(ch)) if x in seg_sites]
            chName = pop.chromName(ch)
            for idx,i in enumerate(seg_sites_ch):
                for j in seg_sites_ch[idx+1:]:
                    dist = int(pop.locusPos(j) - pop.locusPos(i))
                    if dist < maxDist and len(pairs[dist/binDist]) < maxPairs:
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
            sim.Stat(pop, LD = pairs[dist], LD_param = {'stat': ['LD_prime', 'R2']})
            avgLD = [0]*pop.numSubPop()
            avgR2 = [0]*pop.numSubPop()
            for sp in range(pop.numSubPop()):
                LD_var = pop.dvars(sp).LD_prime
                R2_var = pop.dvars(sp).R2
                for pair in pairs[dist]:
                    avgLD[sp] += LD_var[pair[0]][pair[1]]
                    avgR2[sp] += R2_var[pair[0]][pair[1]]
                result['LD_prime'][pop.subPopName(sp)].append(avgLD[sp] / len(pairs[dist]))
                result['R2'][pop.subPopName(sp)].append(avgR2[sp] / len(pairs[dist]))
            progress.update(dist + 1)
        progress.done()
        return result



def _identifyCodonInRegions(raw_regions):
    '''Identify codon of a region'''
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
    # pos2idx: pos->idx in region
    # ......................................................................................
    # coding_loci: positions of xxxxx
    #
    #
    # codon_info stores information  about 1 or more codon at codng regions (excluding introns)
    codon_info = {}
    # coding_base stores reference sequence at coding regions
    coding_base = {}
    #
    with Project(verbosity='1') as proj:
        # expand user provided regions to one or more (chr,start,end,comment)
        regions = expandRegions(raw_regions, proj)
        #
        # all loci contrains all pos in the region and their indexes
        all_loci = defaultdict(set)
        for reg in regions:
            all_loci[reg[0]] = all_loci[reg[0]].union(range(reg[1], reg[2]+1))
        # find the number of loci on each chromosome
        chroms = sorted(all_loci.keys())
        start_idx = 0
        pos2idx = {}
        for ch in sorted(all_loci.keys()):
            pos2idx.update({(ch,y):(start_idx + x) for x,y in enumerate(sorted(all_loci[ch]))})
            start_idx += len(pos2idx)
        #
        ref = RefGenome(proj.build)
        genes = genesInRegions(regions, proj)
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
            ch = reg[0]
            env.logger.info('Length of coding regions of {}: {}'.format(
                gene, len(pos)))
            # now, try to divide coding regions by codon. Note that a codon
            # can consist of nucleotie across two exon.
            skip_codon = True
            for idx, (p, s) in enumerate(zip(pos, seq)):
                if idx % 3 == 0:
                    # the complete codon must be in the simulated region. We do not handle
                    # partial codon because we cannot control mutations outside of the specified
                    # region
                    if (p not in all_loci[ch]) or (pos[idx+1] not in all_loci[ch]) or (pos[idx+2] not in all_loci[ch]):
                        skip_codon = True
                        continue
                    # information about the codon: p0, p1, p2, aa, strand
                    codon = (pos2idx[(ch, p)], pos2idx[(ch, pos[idx+1])], pos2idx[(ch, pos[idx+2])],
                        codon_table[s + seq[idx+1] + seq[idx+2]] if stru['strand'] == '+' else
                        codon_table_reverse_complement[s + seq[idx+1] + seq[idx+2]] ,
                        stru['strand'])
                    skip_codon = False
                if skip_codon:
                    continue
                # record reference sequence
                coding_base[pos2idx[(ch,p)]] = s
                # other two positions share the same codon
                # because of there can be multiple genes, one basepair can be in multiple codon
                if pos2idx[(ch,p)] in codon_info:
                    if codon not in codon_info[pos2idx[(ch,p)]]:
                        codon_info[pos2idx[(ch,p)]].append(codon)
                else:
                    codon_info[pos2idx[(ch,p)]] = [codon]
    #
    # remove all positions that are not in regions
    env.logger.info('{} out of {} bp ({:.2f}%) are in coding regions of genes {}'.format(
        len(codon_info), len(pos2idx), 
        100. * len(codon_info) / len(pos2idx), 
        ', '.join(genes)))
    return len(pos2idx), coding_base, codon_info

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
        self.num_loci, self.coding_base, self.codon_info = _identifyCodonInRegions(regions)
        if (not self.codon_info) or (s_missense == 0 and s_stoploss == 0 and s_stopgain == 0):
            env.logger.warning('Specified region does not contain any gene or all selection coefficient is zero. A neutral model will be used.')
            sim.PySelector.__init__(self, func=self._neutral, loci=[])
        else:
            # we need to send simuPOP to indexes of coding loci within the specified region
            sim.PySelector.__init__(self, func=self._select, loci=sorted(self.codon_info.keys()))

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

    def _select(self, mut):
        # geno is arranged locus by locus (A1,A2,B1,B2 etc)
        #
        # we can not divide the sequence into triplets because it is possible that a nucleotide
        # belong to multiple codon with different locations.
        # 
        # the same aa change can be caused by two mutations, we need to keep only one
        # aa_change with the same (aa, naa, ploidy, p0)
        aa_change = set()
        N = self.num_loci
        for m in mut.keys():
            p = m / N  # ploidy
            loc = m % N
            aa_change_at_m = []
            # we need to map index in coding region to their original position
            for p0, p1, p2, aa, s in self.codon_info[loc]:
                # p0: location
                # mut: mutation (0 for wildtype)
                codon = self.mutant_map[self.coding_base[p0]][mut[p0+p*N]] + \
                        self.mutant_map[self.coding_base[p1]][mut[p1+p*N]] + \
                        self.mutant_map[self.coding_base[p2]][mut[p2+p*N]]
                naa = codon_table[codon] if s == '+' else codon_table_reverse_complement[codon]
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
        return fitness

class ProteinPenetrance(sim.PyPenetrance):
    def __init__(self, regions, s_sporadic=0.0001, s_missense=0.001, 
            s_stoploss=0.002, s_stopgain=0.01):
        '''A protein penetrance model that is identical to ProteinSelector, but 
            use 1 minus calculated fitness value as pentrance probability.
        '''
        self.s_sporadic = s_sporadic
        self.s_missense = s_missense
        self.s_stoploss = s_stoploss
        self.s_stopgain = s_stopgain
        #
        self.num_loci, self.coding_base, self.codon_info = _identifyCodonInRegions(regions)
        if not self.codon_info:
            env.logger.warning('Specified region does not contain any gene. A neutral model will be used.')
            sim.PyPenetrance.__init__(self, func=self._neutral, loci=[])
        else:
            # we need to send simuPOP to indexes of coding loci within the specified region
            sim.PyPenetrance.__init__(self, func=self._penetrance,
                loci=sorted(self.codon_info.keys()))

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
        return 1 - s_sporadic

    def _penetrance(self, mut):
        # geno is arranged locus by locus (A1,A2,B1,B2 etc)
        #
        # we can not divide the sequence into triplets because it is possible that a nucleotide
        # belong to multiple codon with different locations.
        # 
        # the same aa change can be caused by two mutations, we need to keep only one
        # aa_change with the same (aa, naa, ploidy, p0)
        aa_change = set()
        N = self.num_loci
        for m in mut.keys():
            p = m / N  # ploidy
            loc = m % N
            aa_change_at_m = []
            # we need to map index in coding region to their original position
            for p0, p1, p2, aa, s in self.codon_info[loc]:
                # p0: location
                # mut: mutation (0 for wildtype)
                codon = self.mutant_map[self.coding_base[p0]][mut[p0+p*N]] + \
                        self.mutant_map[self.coding_base[p1]][mut[p1+p*N]] + \
                        self.mutant_map[self.coding_base[p2]][mut[p2+p*N]]
                naa = codon_table[codon] if s == '+' else codon_table_reverse_complement[codon]
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
        #
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


class EvolvePop(SkiptableAction):
    def __init__(self,
        selector = None, demoModel=None, 
        mutModel='K80', mutRate=[1.8e-8, 2], 
        recScale=1, output=[]):
        self.mutModel = mutModel
        self.mutRate = mutRate
        self.selector = selector
        self.demoModel = demoModel
        self.recScale = recScale
        self.output = [output]
        SkiptableAction.__init__(self, cmd='EvolvePop mutModel={} mutRate={}, recScale={} output={}\n'
            .format(mutModel, mutRate, recScale, output), output=output)

    def _execute(self, ifiles, pipeline):
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
        # try to use a genetic map
        genetic_pos = [pop.dvars().geneticMap[x] for x in range(pop.totNumLoci())]
        rec_rate = [genetic_pos[i+1] - genetic_pos[i] for i in range(pop.totNumLoci()-1)] + [0.5]
        rec_rate = [x * self.recScale if x >=0 else 0.5 for x in rec_rate]
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
                sim.AcgtMutator(model=self.mutModel, rate=self.mutRate, loci=base['A']),
                sim.AcgtMutator(model=self.mutModel, rate=self.mutRate, loci=base['C'],
                    # base is C=0, 0,1,2,3 to C,G,T,A (1,2,3,0)
                    mapIn=[1,2,3,0], mapOut=[3,0,1,2]),
                sim.AcgtMutator(model=self.mutModel, rate=self.mutRate, loci=base['G'],
                    # base is G=0, 0,1,2,3 to G,T,A,C (2,3,0,1)
                    mapIn=[2,3,0,1], mapOut=[2,3,0,1]),
                sim.AcgtMutator(model=self.mutModel, rate=self.mutRate, loci=base['T'],
                    # base is T=0, 0,1,2,3 to T,A,C,G (3,0,1,2)
                    mapIn=[3,0,1,2], mapOut=[1,2,3,0]),
                # selection on all loci
                #sim.TicToc(),
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
                sim.MendelianGenoTransmitter() if self.recScale == 0 else sim.Recombinator(rates=rec_rate, loci=sim.ALL_AVAIL),
                subPopSize=self.demoModel),
            finalOps=[
                # revert fixed sites so that the final population does not have fixed sites
                sim.RevertFixedSites(),
                #
                # apply fitness to get final statistics and facilitate sampling
                self.selector,
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


class DrawCaseControlSample(SkiptableAction):
    def __init__(self, cases, controls, penetrance, output):
        self.cases = cases
        self.controls = controls
        self.penetrance = penetrance
        self.selectedCases = 0
        self.selectedCtrls = 0
        SkiptableAction.__init__(self, cmd='CaseCtrlSampler {}\n'.format(output),
            output=output)
    
    def _execute(self, ifiles, pipeline):
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        self.prog = ProgressBar('Generating %d cases and %d controls...' % (self.cases, self.controls), self.cases + self.controls)
        pop.evolve(
            matingScheme=sim.RandomMating(
                ops=[
                    sim.MendelianGenoTransmitter(),
                    # apply a penetrance model 
                    self.penetrance,
                    # an individual will be discarded if _selectInds returns False
                    sim.PyOperator(func=self._selectInds)
                ],
                subPopSize=self.cases + self.controls,
            ),
            gen = 1
        )
        self.prog.done()
        # 
        if 'migrate_to' in pop.infoFields():
            pop.removeInfoFields('migrate_to')
        if 'fitness' in pop.infoFields():
            pop.removeInfoFields('fitness')
        env.logger.info('Saving samples to population {}'.format(self.output[0]))
        pop.save(self.output[0])

    def _selectInds(self, off):
        'Determine if the offspring can be kept.'
        if off.affected():
            if self.selectedCases < self.cases:
                self.selectedCases += 1
                self.prog.update(self.selectedCases + self.selectedCtrls)
                return True
        else:
            # we keep number of ctrls less than cases to keep the progress even 
            # because it is generally much easier to find controls
            if self.selectedCtrls < self.controls and self.selectedCtrls <= self.selectedCases:
                self.prog.update(self.selectedCases + self.selectedCtrls)
                self.selectedCtrls += 1
                return True
        return False

class DrawRandomSample(SkiptableAction):
    def __init__(self, sizes, output):
        self.sizes = sizes
        SkiptableAction.__init__(self, cmd='CaseCtrlSampler {}\n'.format(output),
            output=output)
    
    def _execute(self, ifiles, pipeline):
        env.logger.info('Loading {}'.format(ifiles[0]))
        pop = sim.loadPopulation(ifiles[0])
        pop.evolve(
            matingScheme=sim.RandomMating(
                subPopSize=self.sizes
            ),
            gen = 1
        )
        # 
        if 'migrate_to' in pop.infoFields():
            pop.removeInfoFields('migrate_to')
        if 'fitness' in pop.infoFields():
            pop.removeInfoFields('fitness')
        env.logger.info('Saving samples to population {}'.format(self.output[0]))
        pop.save(self.output[0])

