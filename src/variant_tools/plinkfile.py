#!/usr/bin/env
# -*- coding: utf8 -*-
# The libplinkio software library is distributed under the following terms:

# Copyright (c) 2012-2013, Matias Frånberg
# All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.

# * Neither the name of Mattias Frånberg nor the
#   names of its contributors may be used to endorse or promote products
#   derived from this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from . import cplinkio


class PlinkFile:
    ##
    # Opens the plink file at the given path.
    #
    # @param path The prefix for a .bed, .fam and .bim without
    #             the extension. E.g. for the files /plink/myfile.fam,
    #             /plink/myfile.bim, /plink/myfile.bed use the path
    #             /plink/myfile
    #
    def __init__(self, path):
        self.path = path
        self.handle = cplinkio.open(path)
        self.loci = cplinkio.get_loci(self.handle)
        self.samples = cplinkio.get_samples(self.handle)

    ##
    # Returns an iterator from the beginning of
    # the file.
    #
    def __iter__(self):
        cplinkio.reset_row(self.handle)

        return self

    ##
    # Returns the prefix path to the plink file, e.g.
    # without .bim, .bed or .fam.
    #
    def get_path(self):
        return self.path

    ##
    # Returns a list of the samples.
    #
    def get_samples(self):
        return self.samples

    ##
    # Returns a list of the loci.
    #
    def get_loci(self):
        return self.loci

    ##
    # Determines how the snps are stored. It will return
    # true if a row contains the genotypes of all individuals
    # from a single locus, false otherwise.
    #
    def one_locus_per_row(self):
        return cplinkio.one_locus_per_row(self.handle)

    ##
    # Goes to next row.
    #
    def __next__(self):
        row = cplinkio.next_row(self.handle)
        if not row:
            raise StopIteration
        return row

    ##
    # Closes the file.
    #
    def close(self):
        if self.handle:
            cplinkio.close(self.handle)
            self.handle = None

    ##
    # Transposes the file.
    #
    def transpose(self, new_path):
        return cplinkio.transpose(self.path, new_path)


class Sample:

    def __init__(self,
                 fid,
                 iid,
                 father_iid,
                 mother_iid,
                 sex,
                 affection,
                 phenotype=0.0):
        ##
        # Family id.
        #
        self.fid = fid

        ##
        # Individual id.
        #
        self.iid = iid

        ##
        # Individual id of father.
        #
        self.father_iid = father_iid

        ##
        # Individual id of mother.
        #
        self.mother_iid = mother_iid

        ##
        # Sex of individual.
        #
        self.sex = sex

        ##
        # Affection of individual, 0/1, case/control
        #
        self.affection = affection

        ##
        # Optional continuous phenotype
        #
        self.phenotype = phenotype

    def __str__(self):
        return "{0} {1} {2} {3}".format(self.fid, self.iid, self.sex,
                                        self.affection)


class Locus:

    def __init__(self, chromosome, name, position, bp_position, allele1,
                 allele2):
        ##
        # Chromosome number starting from 1
        #
        self.chromosome = chromosome

        ##
        # Name of the loci, usually RS number or
        # chrX:pos.
        #
        self.name = name

        ##
        # Position.
        #
        self.position = position

        ##
        # Base pair position.
        #
        self.bp_position = bp_position

        ##
        # First allele
        #
        self.allele1 = allele1

        ##
        # Second allele
        #
        self.allele2 = allele2

    def __str__(self):
        return "{0} {1}".format(self.chromosome, self.name)


##
# Opens the plink file at the given path.
#
# @param path The prefix for a .bed, .fam and .bim without
#             the extension. E.g. for the files /plink/myfile.fam,
#             /plink/myfile.bim, /plink/myfile.bed use the path
#             /plink/myfile
#
def open(path):
    return PlinkFile(path)
