#!/usr/bin/env python
#
# $File: preprocessor.py $
# $LastChangedDate: 2013-03-26 17:15:07 -0500 (Tue, 26 Mar 2013) $
# $Rev: 1775 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 Bo Peng (bpeng@mdanderson.org) and Gao Wang (wangow@gmail.com)
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


import sys, os
import itertools as it
from .utils import ProgressBar, RefGenome, env
from .plinkfile import PlinkFile

class BatchWriter:
    '''write text to file in batches (#lines)'''
    def __init__(self, fn, batch = 1000):
        self.batch = batch
        if os.path.exists(fn):
            os.remove(fn)
        self.fs = open(fn, 'a')
        self.counter = 0
        self.swap = ''

    def write(self, line):
        if line is None:
            self.fs.write(self.swap)
            self.fs.close()
        else:
            self.swap += line
            self.counter += 1
            if self.counter == self.batch:
                # time to write
                self.fs.write(self.swap)
                self.counter = 0
                self.swap = ''

class PlinkBinaryToVariants:
    """
    Class to write a PLINK BED format genotype dataset (.bed, .fam and .bim files)
    to vtools compatible variant format 
    
    c.f., http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

    @param    dataset: path + prefix for a .bed, .fam and .bim without extensions
    @param      build:
    @param     logger:

    Note: implementation of this class is based on libplinkio
    https://bitbucket.org/mattias_franberg/libplinkio

    The SNP will be encoded as follows:
    * 0 - Homozygous major
    * 1 - Heterozygous
    * 2 - Homozygous minor
    * 3 - Missing value

    Output of this program will retain such coding if the major allele is reference genotype
    in reference human genome; otherwise will re-code it by:
    * 0 --> 2
    * 2 --> 0

    WARNING
    -------
    The allele1 and allele2 in the *.bim file does not tell which is major and which is minor
    so treating allele1 as major and allele2 as minor will result in completely swapped
    genotype coding if it is not the case !!

    self.determineMajorAllele() Attempts to resolve this issue
    """
    def __init__(self, dataset, build, chrom_namemap = {}):
        # check file path
        for ext in ['.fam', '.bed', '.bim']:
            if not os.path.exists(dataset + ext):
                raise RuntimeError('Cannot find file {0}'.format(dataset + ext))
        self.dataset = dataset
        self.build = build
        self.cur = PlinkFile(self.dataset)
        # a list of sample names (sample ID's in .fam file)'''
        self.samples =  [x.iid for x in self.cur.get_samples()]
        if None in self.samples:
            raise ValueError("Cannot read sample ID from malformed '{0}.fam' file.".\
                             format(self.dataset))
        # iterator for variants info: chr, pos, allele1, allele2
        self.variants = it.chain(self.cur.get_loci())
        # reference genome object and ATCG dictionary
        self.hgref = RefGenome(build)
        self.CSTRANDS = {'A':'T',
                         'G':'C',
                         'T':'A',
                         'C':'G'}
        # status 0 for not flip, 1 for flip, -1 for bad match
        self.status = 0
        # chromosome naming convention map
        self.cmap = chrom_namemap
        # variant writer
        self.variant_writer = None
        self.data_writer = None

    def initWriter(self, ofile):
        self.variant_writer = BatchWriter(self.dataset + ".{0}adjusted".format(self.build), batch = 1000)
        self.data_writer = BatchWriter(ofile, batch = 5)

    def getValidatedLocusGenotype(self, chrom, pos, allele1, allele2, geno_cur):
        '''Use cgatools to obtain validated genotype for given locus.

        Assume allele1 is major allele and allele2 is minor allele
        - check which of the two input allele is reference allele
        - if none is reference, use the alternative strand and check again
        - if the locus is validated (self.status >= 0)
          - flip ref/alt coding if the minor allele is found to be reference allele
          - output its genotype
          otherwise return None
        @return a locus genotypes string
        '''
        if chrom in self.cmap:
            chrom = self.cmap[chrom]
        try:
            ref = self.hgref.getBase(chrom, pos)
        except Exception as e:
            env.logger.warning('Cannot find genomic coordinate {0}:{1} in reference genome {2}. '
                                'Input variant is ignored'.format(chrom, pos, self.build))
            self.status = -1
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(chrom, pos, 0, 0))
            return None
        self.status, strand, allele1, allele2 = self._matchref(ref, allele1, allele2)
        if self.status < 0:
            if self.status == -2:
                env.logger.warning('All genotypes for variant "{0}:{1}" are missing'.\
                                    format(chrom, pos))
            else:
                env.logger.warning('Variant "{0}:{1} {2} {3}" failed '
                                   'to match reference genome {4}/(A,T,C,G)'.\
                                    format(chrom, pos, allele1, allele2, ref))
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(chrom, pos, 0, 0))
            return None
        elif self.status == 0:
            if strand:
                env.logger.debug('Use alternative strand for {0}:{1}'.format(chrom, pos))
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(chrom, pos, allele1, allele2))
            return ','.join([chrom, str(pos), allele1, allele2]) + ',' + str(geno_cur)
        else:
            # have to flip the genotypes coding
            if strand:
                env.logger.debug('Use alternative strand for {0}:{1}'.format(chrom, pos))
            # env.logger.debug('Allele coding flipped for {0}:{1}'.format(chrom, pos))
            # Very time consuming compare to not flipping the genotype codes
            if self.variant_writer:
                self.variant_writer.write("{}\t{}\t{}\t{}\n".format(chrom, pos, allele2, allele1))
            return ','.join([chrom, str(pos), allele2, allele1]) + ',' + \
                ','.join([str(x) if x == 3 or x == 'E' else str(2 - x) for x in geno_cur])
            
    def getLociCounts(self):
        # FIXME: not efficient
        return len(self.cur.get_loci())
        
    def getHeader(self):
        '''a line of headers for the output text file'''
        return ','.join(
            ['#chr', 'pos', 'ref', 'alt'] + self.samples
            )
        
    def getLine(self, viter = None, giter = None, which_major = 1):
        '''a line of validated genotypes for a locus'''
        if viter is None: viter = self.variants
        if giter is None: giter = self.cur
        try:
            locus = next(viter)
            genotypes = giter.next()
        except StopIteration:
            self.cur.close()
            return False, None
        except:
            env.logger.error('Failed to retrieve locus {0}:{1} '
                                '(plinkio error)'.format(locus.chromosome, locus.bp_position))
            return True, None
        if which_major == 1:
            # allele 1 is the major allele
            return True, self.getValidatedLocusGenotype(str(locus.chromosome), int(locus.bp_position),
                                                        locus.allele1.upper(), locus.allele2.upper(),
                                                        genotypes)
        else:
            # allele 2 is the major allele
            return True, self.getValidatedLocusGenotype(str(locus.chromosome), int(locus.bp_position),
                                                        locus.allele2.upper(), locus.allele1.upper(),
                                                        genotypes)

    def determineMajorAllele(self, n=1000):
        '''The logic here is that for the first n loci we
        - assume the allele1 in *.bim is major allele
        - try to map to hg reference and see how many genotype coding have to be flipped
        - if too many (over half) have to be flipped we conclude that allele2 should be major allele
        @return: a guess of major allele based on the first n sample loci. -1 for bad matching
        '''
        # new temporary connection
        cur = PlinkFile(self.dataset)
        variants = it.chain(cur.get_loci())
        # counters
        m_ones = zeros = ones = 0
        i = 0
        while i < n:
            self.status = 0
            # self.status will be updated
            flag, line = self.getLine(viter = variants, giter = cur)
            if not flag:
                # end of line
                n = i
                break
            if self.status == 0:
                # not flipped
                zeros += 1
            elif self.status == 1:
                # flipped
                ones += 1
            elif self.status == -1:
                # bad match
                m_ones += 1
            else:
                # status = -2, monomorphic site
                # ignore this test
                if self.status == -2:
                    continue
            i += 1
        # check if so many are negative values
        if m_ones > float(n) / 2.0:
            return -9
        if ones > float(n - m_ones) / 2.0:
            # allele 2 seems to be major allele
            return 2
        else:
            return 1 
        

    def _matchref(self, ref, major, minor):
        '''try best to match reference allele
        @param    ref: reference base coding
        @param  major: major allele coding
        @param  minor: minor allele coding
        @return self.status (0 for no need to flip, 1 for having to flip)
        @return strand (0 for original, 1 for alternative),
        @return major and minor alleles (might be from alternative strand)

        Notice
        ------
        In *.bim file allele1 or allele2 (major/minor) can be 0 when only one allele is found
        in data. Such loci is not considered a variant site and will be ignored
        '''
        # monomorphic, missing, or invalid coding
        if major not in ['A','T','C','G','0'] or minor not in ['A','T','C','G','0']:
            return -1, 0, major, minor
        if major == '0' and minor == '0':
            return -2, 0, major, minor
        if major == '0' and minor != '0':
            major = minor
        if major != '0' and minor == '0':
            minor = major
        status = strand = 0
        if ref in [major, minor]:
            # allele found, determine flip status
            status = 0 if ref == major else 1
        else:
            # allele not found, use alternative strand
            strand = 1
            major = self.CSTRANDS[major]
            minor = self.CSTRANDS[minor]
            # allele still not found
            if ref not in [major, minor]:
                status = -1
            # allele found, determine flip status
            else:
                status = 0 if ref == major else 1
        return status, strand, major, minor

#
#
# Preprocessors of input files
# They will convert input files to intermediate variant based text files for easy import 
#
#

class Preprocessor:
    def __init__(self):
        '''Base preprocessor class that converts $files
        and write intermediate output files to $outdir/file.format'''
        pass

    def convert(self, files, output_files):
        for item, ofile in zip(files, output_files):
            env.logger.info('Convert {} to {}'.format(item, ofile))
        

class PlinkConverter(Preprocessor):
    def __init__(self, build, chrom_namemap = {}):
        Preprocessor.__init__(self)
        self.build = build
        self.cmap = chrom_namemap

    def convert(self, files, output_files):
        for item, ofile in zip(files, output_files):
            if os.path.exists(item + ".bed"):
                self.decode_plink(PlinkBinaryToVariants(item, self.build, self.cmap),
                                  ofile)
            else:
                import glob
                files = '/'.join([x for x in glob.glob(item + '*')])
                if files:
                    supported = ['*.bed/*.bim/*.fam']
                    raise ValueError("Unsupported input file '{}' (supported file types are {})".\
                                         format(files, ';'.join(supported)))
                else:
                    raise ValueError("Cannot find input files '{}'".format(item + '*'))
            
    def decode_plink(self, p2v, ofile, n = 1000):
        '''decode plink data from p2v object and output to ofile'''
        env.logger.info("Determining major/minor allele from data")
        # check major allele
        which_major = p2v.determineMajorAllele(n)
        # raise on bad match
        if which_major == -9:
            raise ValueError ('Invalid dataset {0}: too many unmatched loci to reference genome {1}. '
                              'Perhaps you specified the wrong build, or have too many unsupported allele '
                              'types (not A/T/C/G, e.g, indels I/D) in BED file which you have to remove before '
                              'import'.format(p2v.dataset, p2v.build))
        env.logger.debug("allele{} is major allele".format(which_major))
        # output
        nloci = p2v.getLociCounts()
        batch = int(nloci / 100)
        prog = ProgressBar('Decoding {0}'.format(p2v.dataset), nloci)
        p2v.initWriter(ofile)
        p2v.data_writer.write(p2v.getHeader() + '\n')
        p2v.variant_writer.write("#chr\tpos\tref\talt\n")
        count = 0
        while True:
            flag, line = p2v.getLine(which_major = which_major)
            count += 1
            if not flag:
                prog.done()
                p2v.data_writer.write(None)
                p2v.variant_writer.write(None)
                break
            else:
                if line is not None:
                    p2v.data_writer.write(line + '\n')
                if count % batch == 0 and count > batch:
                    prog.update(count)
