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
from .utils import ProgressBar, DatabaseEngine, delayedAction, env

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
def simulateArguments(parser):
    pass

def simulate(args):
    try:
        with Project(verbosity=args.verbosity) as proj:
            # step 0: 
            # get the model of simulation
            pass
    except Exception as e:
        env.logger.error(e)
        sys.exit(1)


