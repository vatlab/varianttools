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

import sys, os, re
from .project import Project
from .utils import ProgressBar, DatabaseEngine, delayedAction, env,\
    consolidateFieldName

if sys.version_info.major == 2:
    from ucsctools_py2 import tabixFetch
else:
    from ucsctools_py3 import tabixFetch

import argparse

'''Draft design of the simulation feature

Goals:
1. Simulate 'real' data in the sense that the simulated datasets should have
  chromsome, locations, build, and preferrable realistic regions for 
  exome, etc.

2. Users provide description of samples, preferrably NOT the way to simulate
  it.

3. Results are written as variant tools projects, and should be easy to analyze
  using the subsequent tools. Results can be exported in vcf format.

4. Simulate SNPs, and Indels if model allows.

Design:
1. Use the existing 'show' feature to show available models.
  a. vtools show simulations/models/simu_models
  b. vtools show model SOME_MODEL

2. Use the existing 'snapshot' feature to distribute pre-simulated
  datasets.
  a. vtools show snapshots
  b. vtools admin --load_snapshot

3. Use the existing 'export' feature to export simulated data in vcf format
  a. vtools export --output simulated.vcf

4. Use the existing 'pipeline' feature to simulate complex samples if that
  evolves multiple steps.

Implementations:



Models:
1. Sequencing error model: manipulate existing (simulated) genotype
2. De Nova mutation model: create novel mutations for offspring
3. Phenotype model: create phenotype based on genotypes on selected
   variants.
4. Indel models.
5. Resampling model: easy to implement
6. Coalescane model: ... many stuff with python interface
7. forward-time model: simupop, srv

'''


class NullModel:
    '''A base class that defines a common interface for simulation models'''
    def __init__(self, *method_args):
        '''Args is arbitrary arguments, might need an additional parser to
        parse it'''
        # trait type
        self.trait_type = None
        # group name
        self.gname = None
        self.fields = []
        self.parseArgs(*method_args)
        #

    def parseArgs(self, method_args):
        # this function should never be called.
        raise SystemError('All association tests should define their own parseArgs function')


def getAllModels():
    '''List all simulation models (all classes that subclasses of NullModel) in this module'''
    return sorted([(name, obj) for name, obj in globals().iteritems() \
        if type(obj) == type(NullModel) and issubclass(obj, NullModel) \
            and name not in ('NullModel',)], key=lambda x: x[0])


def simuRareVariants(regions, N, G, mu, selDist, selCoef, selModel='exponential', recRate=0, 
        splitTo=[1], splitAt=0, migrRate=0, steps=[100], mutationModel='finite_sites',
        initPop='', extMutantFile='', addMutantsAt=0, postHook=None,
        statFile='', popFile='', markerFile='', mutantFile='', genotypeFile='',
        verbose=1, logger=None):
    '''
    Please refer to simuRareVariants.py -h for a detailed description of all parameters.
    Note that a user-defined function can be passed to parameter selDist to specify
    arbitrary distribution of fitness. A script-only feature is that a Python function
    can be provided through parameter postHook to process the population at each generation.
    '''
    #
    # convert regions to start/end positions
    ranges = []
    chromTypes = []
    for region in regions:
        start, end = [int(x) for x in region.split(':')[1].split('..')]
        ranges.append((start, end+1))
        if region.split(':')[0] == 'chrX':
            chromTypes.append(sim.CHROMOSOME_X)
            if len(regions) > 1:
                raise ValueError('The current implementation only allows one region if it is on chromosome X')
            logger.info('Chromosome {} is on chromosome X'.format(region))
        elif region.split(':')[0] == 'chrY':
            raise ValueError('The current implementation does not support chromosome Y')
            chromTypes.append(sim.CHROMOSOME_Y)
            logger.info('Chromosome {} is on chromosome Y'.format(region))
        else:
            chromTypes.append(sim.AUTOSOME)
    if logger:
        logger.info('%s regions with a total length of %d basepair.' % (len(ranges), sum([x[1]-x[0] for x in ranges])))
    #
    # set default parameter
    if selCoef is None:
        # set default parameters
        if selDist == 'mixed_gamma':
            selCoef = [0.0186, 0.0001, 0.184, 0.160*2, 0.5, 0.0001, 0.1]
        elif selDist == 'mixed_gamma1':
            selCoef = [0, -1, 0.562341, 0.01, 0.5, 0.00001, 0.1]
        elif selDist == 'gamma1':
            selCoef = [0.23, 0.185*2, 0.5]
        elif selDist == 'gamma2':
            selCoef = [0.184, 0.160*2, 0.5]
        elif selDist == 'gamma3':
            selCoef = [0.206, 0.146*2, 0.5]
        elif selDist == 'constant':
            selCoef = [0.01, 0.5]
        elif not callable(selDist):
            raise ValueError("Unsupported random distribution")
    else:
        # force to list type
        selCoef = list(selCoef)
    if len(steps) == 0:
        # at the end of each stage
        steps = G
    elif len(steps) == 1:
        # save step for each stage
        steps = steps * len(G)
    # use a right selection operator.
    collector = fitnessCollector()
    mode = {'multiplicative': sim.MULTIPLICATIVE,
        'additive': sim.ADDITIVE,
        'exponential': sim.EXPONENTIAL}[selModel]
    #
    if type(popFile) == str:
        popFile = [popFile, -1]
    #
    if callable(selDist):
        mySelector = MutSpaceSelector(selDist=selDist, mode=mode, output=collector.getCoef)
    elif selDist == 'mixed_gamma':
        mySelector = MutSpaceSelector(selDist=mixedGamma(selCoef), mode=mode, output=collector.getCoef)
    elif selDist == 'mixed_gamma1':
        mySelector = MutSpaceSelector(selDist=mixedGamma1(selCoef), mode=mode, output=collector.getCoef)
    elif selDist.startswith('gamma'):
        mySelector = MutSpaceSelector(selDist=[sim.GAMMA_DISTRIBUTION]+selCoef,
            mode=mode, output=collector.getCoef)
    elif selDist == 'constant':
        if selCoef == 0:
            mySelector = sim.NoneOp()
        else:
            mySelector = MutSpaceSelector(selDist=[sim.CONSTANT] + selCoef,
                mode=mode, output=collector.getCoef)
    #
    # Evolve
    if os.path.isfile(initPop):
        if logger:
            logger.info('Loading initial population %s...' % initPop)
        pop = sim.loadPopulation(initPop)
        if pop.numChrom() != len(regions):
            raise ValueError('Initial population %s does not have specified regions.' % initPop)
        for ch,reg in enumerate(regions):
            if pop.chromName(ch) != reg:
                raise ValueError('Initial population %s does not have region %s' % (initPop, reg))
        pop.addInfoFields(['fitness', 'migrate_to'])
    else:
        pop = sim.Population(size=N[0], loci=[10]*len(regions), chromNames=regions,
            infoFields=['fitness', 'migrate_to'], chromTypes=chromTypes)
    if logger:
        startTime = time.clock()
    #
    progGen = []
    # 0, G[0], G[0]+G[1], ..., sum(G)
    Gens = [sum(G[:i]) for i in range(len(G)+1)]
    for i in range(len(Gens)-1):
        progGen += range(Gens[i], Gens[i+1], steps[i])
    pop.evolve(
        initOps=sim.InitSex(),
        preOps=[
            sim.PyOutput('''Statistics outputted are
1. Generation number,
2. population size (a list),
3. number of segregation sites,
4. average number of segregation sites per individual
5. average allele frequency * 100
6. average fitness value
7. minimal fitness value of the parental population
''', at = 0)] + \
            [sim.PyOutput('Starting stage %d\n' % i, at = Gens[i]) for i in range(0, len(Gens))] + \
            # add alleles from an existing population 
            [sim.IfElse(extMutantFile != '',
                ifOps = [
                    sim.PyOutput('Loading and converting population %s' % extMutantFile),
                    sim.PyOperator(func=addMutantsFrom, param=(extMutantFile, regions, logger)),
                ], at = addMutantsAt),
            # revert alleles at fixed loci to wildtype
            RevertFixedSites(),
            # mutate in a region at rate mu, if verbose > 2, save mutation events to a file
            MutSpaceMutator(mu, ranges, {'finite_sites':1, 'infinite_sites':2}[mutationModel],
                output='' if verbose < 2 else '>>mutations.lst'),
            # selection on all loci
            mySelector,
            # output statistics in verbose mode
            sim.IfElse(verbose > 0, ifOps=[
                sim.Stat(popSize=True, meanOfInfo='fitness', minOfInfo='fitness'),
                NumSegregationSites(),
                sim.PyEval(r'"%5d %s %5d %.6f %.6f %.6f %.6f\n" '
                    '% (gen, subPopSize, numSites, avgSites, avgFreq*100, meanOfInfo["fitness"], minOfInfo["fitness"])',
                    output='>>' + statFile),
                ], at = progGen
            ),
            sim.IfElse(len(splitTo) > 1,
                sim.Migrator(rate=migrIslandRates(migrRate, len(splitTo)),
                    begin=splitAt + 1)
            ),
        ],
        matingScheme=sim.RandomMating(ops=MutSpaceRecombinator(recRate, ranges),
            subPopSize=multiStageDemoFunc(N, G, splitTo, splitAt)),
        postOps = [
            sim.NoneOp() if postHook is None else sim.PyOperator(func=postHook),
            sim.SavePopulation(popFile[0], at=popFile[1]),
        ],
        finalOps=[
            # revert fixed sites so that the final population does not have fixed sites
            RevertFixedSites(),
            sim.IfElse(verbose > 0, ifOps=[
                # statistics after evolution
                sim.Stat(popSize=True),
                NumSegregationSites(),
                sim.PyEval(r'"%5d %s %5d %.6f %.6f %.6f %.6f\n" '
                    '% (gen+1, subPopSize, numSites, avgSites, avgFreq*100, meanOfInfo["fitness"], minOfInfo["fitness"])',
                    output='>>' + statFile),
                sim.PyEval(r'"Simulated population has %d individuals, %d segregation sites.'
                           r'There are on average %.1f sites per individual. Mean allele frequency is %.4f%%.\n"'
                           r'% (popSize, numSites, avgSites, avgFreq*100)'),
            ]),
        ],
        gen = Gens[-1]
    )
    # record selection coefficients to population
    if len(collector.selCoef) == 0:
        # this must be the neutral case where a NonOp has been used.
        pop.dvars().selCoef = 0
    else:
        pop.dvars().selCoef = collector.selCoef
    #
    if logger:
        logger.info('Population simulation takes %.2f seconds' % (time.clock() - startTime))
    if markerFile or genotypeFile:
        if logger:
            logger.info('Saving marker information to file %s' % markerFile)
        mutants = saveMarkerInfoToFile(pop, markerFile, logger)
        if genotypeFile:
            if logger:
                logger.info('Saving genotype in .ped format to file %s' % genotypeFile)
            saveGenotypeToFile(pop, genotypeFile, mutants, logger)
    if mutantFile:
        if logger:
            logger.info('Saving mutants to file %s' % mutantFile)
        saveMutantsToFile(pop, mutantFile, logger=logger)
    return pop



def simulateArguments(parser):
    parser.add_argument('--regions', nargs='+',
        help='''One or more chromosome regions in the format of chr:start-end
        (e.g. chr21:33,031,597-33,041,570), or Field:Value from a region-based 
        annotation database (e.g. refGene.name2:TRIM2 or refGene_exon.name:NM_000947).
        Multiple chromosomal regions will be selected if the name matches more
        than one chromosomal regions. The regions will be marked by their indexes
        but values from another field will be used if the its name is appended
        (e.g. refGene.name2:BRCA2:name will mark each regions with name of isoforms).
        Chromosome positions are 1-based and are inclusive at both ends so the 
        chromosome region has a length of end-start+1 bp.''')
    parser.add_argument('--model', 
        help='''Simulation model, which defines the algorithm and default
            parameter to simulate data. A list of model-specific parameters
            could be specified to change the behavior of these models. Commands
            vtools show models and vtools show model MODEL can be used to list
            all available models and details of one model.''')
    

def expandRegions(arg_regions, proj, mergeRegions=True):
    regions = []
    for region in arg_regions:
        try:
            chr, location = region.split(':', 1)
            start, end = location.split('-')
            start = int(start.replace(',', ''))
            end = int(end.replace(',', ''))
            if start == 0 or end == 0:
                raise ValueError('0 is not allowed as starting or ending position')
            # start might be after end
            if start > end:
                regions.append((chr, start, end, '(reverse complementary)'))
            else:
                regions.append((chr, start, end, ''))
        except Exception as e:
            # this is not a format for chr:start-end, try field:name
            try:
                if region.count(':') == 1:
                    field, value = region.rsplit(':', 1)
                    comment_field = "''"
                else:
                    field, value, comment_field = region.rsplit(':', 2)
                # what field is this?
                query, fields = consolidateFieldName(proj, 'variant', field, False) 
                # query should be just one of the fields according to things that are passed
                if query.strip() not in fields:
                    raise ValueError('Could not identify field {} from the present project'.format(field))
                # now we have annotation database
                try:
                    annoDB = [x for x in proj.annoDB if x.linked_name.lower() == query.split('.')[0].lower()][0]
                except:
                    raise ValueError('Could not locate annotation database {} in the project'.format(query.split('.')[0]))
                #
                if annoDB.anno_type != 'range':
                    raise ValueError('{} is not linked as a range-based annotation database.'.format(annoDB.linked_name))
                # get the fields?
                chr_field, start_field, end_field = annoDB.build
                #
                # find the regions
                cur = proj.db.cursor()
                try:
                    cur.execute('SELECT {},{},{},{} FROM {}.{} WHERE {}="{}"'.format(
                        chr_field, start_field, end_field, comment_field, annoDB.linked_name, annoDB.name,
                        field.rsplit('.',1)[-1], value))
                except Exception as e:
                    raise ValueError('Failed to search range and comment field: {}'.format(e))
                for idx, (chr, start, end, comment) in enumerate(cur):
                    if start > end:
                        env.logger.warning('Ignoring unrecognized region chr{}:{}-{} from {}'
                            .format(chr, start + 1, end + 1, annoDB.linked_name))
                        continue
                    try:
                        if comment:
                            regions.append(('chr{}'.format(chr), int(start), int(end), comment))
                        else:
                            regions.append(('chr{}'.format(chr), int(start), int(end), '{} {}'.format(field, idx+1)))
                    except Exception as e:
                        env.logger.warning('Ignoring unrecognized region chr{}:{}-{} from {}'
                            .format(chr, start, end, annoDB.linked_name))
                if not regions:
                    env.logger.error('No valid chromosomal region is identified for {}'.format(region)) 
            except Exception as e:
                raise ValueError('Incorrect format for chromosomal region {}: {}'.format(region, e))
    # remove duplicates and merge ranges
    if not mergeRegions:
        return regions
    regions = list(set(regions))
    while True:
        merged = False
        for i in range(len(regions)):
            for j in range(i + 1, len(regions)):
                r1 = regions[i]
                r2 = regions[j]
                if r1 is None or r2 is None:
                    continue
                # reversed?
                reversed = r1[1] > r1[2] or r2[1] > r2[2]
                if reversed:
                    r1 = (r1[0], min(r1[1], r1[2]), max(r1[1], r1[2]), r1[3])
                    r2 = (r2[0], min(r2[1], r2[2]), max(r2[1], r2[2]), r2[3])
                if r1[0] == r2[0] and r1[2] >= r2[1] and r1[1] <= r2[2]:
                    env.logger.info('Merging regions {}:{}-{} ({}) and {}:{}-{} ({})'
                        .format(r2[0], r2[1], r2[2], r2[3], r1[0], r1[1], r1[2], r1[3]))
                    if reversed:
                        regions[i] = (r1[0], min(r1[1], r2[1]), max(r1[2], r2[2]), r1[3][:-1] + ', ' + r2[3][1:])
                    else:
                        regions[i] = (r1[0], max(r1[2], r2[2]), min(r1[1], r2[1]), r1[3][:-1] + ', ' + r2[3][1:])
                    regions[j] = None
                    merged = True
        if not merged:
            return sorted([x for x in regions if x is not None])


def extractFromVCF(filenameOrUrl, regions, output=''):
    # This is equivalent to 
    #
    # tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/
    #     ALL.2of4intersection.20100804.genotypes.vcf.gz 2:39967768-39967768
    #
    #
    tabixFetch(filenameOrUrl, regions)

def simulate(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # step 0: 
            # get the model of simulation
            print expandRegions(args.regions, proj)
            #extractFromVCF(os.path.expanduser('~/vtools/test/vcf/CEU.vcf.gz'),
            #    ['1:1-100000'])

    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


