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
    consolidateFieldName, RefGenome
from .importer import Importer
from .pipeline import Pipeline



import argparse
import time
import random
import subprocess

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



import simuOpt
simuOpt.setOptions(alleleType='mutant', optimized=False, quiet=True, version='1.1.2')
import simuPOP as sim


from srv import simuRareVariants2, getSelector
from simuPOP.demography import *

demo = MultiStageModel([
    InstantChangeModel(T=100, N0=200),
    ExponentialGrowthModel(T=100, NT=1000)
    ])
                      
                    

def simulateArguments(parser):
    parser.add_argument('model', 
        help='''Simulation model, which defines the algorithm and default
            parameter to simulate data. A list of model-specific parameters
            could be specified to change the behavior of these models. Commands
            vtools show simulations and vtools show simulation MODEL can be used to list
            all available models and details of one model.''')
    parser.add_argument('--seed', type=int,
        help='''Random seed for the simulation. By default, a random seed will be
            used to generate a random sample. A specific seed could be used to 
            reproduce a previously generated sample. The seed for a previous
            simulation run could be found from the simulation configuration
            file $model_$datetime.cfg''')
    



def simulate(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # step 1, create a simulation configuration file.
            model_name = os.path.basename(args.model).split('.', 1)[0]
            if args.seed is None:
                args.seed = random.randint(1, 2**32)
            cfg_file = '{}_{}.cfg'.format(model_name, args.seed)
                #time.strftime('%b%d_%H%M%S', time.localtime()),
            with open(cfg_file, 'w') as cfg:
                cfg.write('model={}\n'.format(args.model))
                cfg.write('seed={}\n'.format(args.seed))
                if '--seed' in sys.argv:
                    # skip the seed option so to stop pipeline from distinguishing the two commands
                    cmd_args = sys.argv[:sys.argv.index('--seed')] + sys.argv[sys.argv.index('--seed') + 2:]
                    cfg.write("command=vtools {}\n".format( 
                        subprocess.list2cmdline(cmd_args[1:])))
                else:
                    cfg.write("command={}\n".format(env.command_line))
            env.logger.info('Simulation configuration is saved to {}'.format(cfg_file))
            pipeline = Pipeline(proj, args.model, extra_args=args.unknown_args)
            pipeline.execute(None, [cfg_file], [], 1)
            #
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


