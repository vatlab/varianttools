import sys, os
import itertools as it
from variant_tools import cgatools_py3 as cga

import cplinkio

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
        self.handle = cplinkio.open( path )
        self.loci = cplinkio.get_loci( self.handle )
        self.samples = cplinkio.get_samples( self.handle )

    ##
    # Returns an iterator from the beginning of
    # the file.
    #
    def __iter__(self):
        cplinkio.reset_row( self.handle )

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
        return cplinkio.one_locus_per_row( self.handle )

    ##
    # Goes to next row.
    #
    def next(self):
        row = cplinkio.next_row( self.handle )
        if not row:
            raise StopIteration

        return row

    ##
    # For python 3.x.
    #
    def __next__(self):
        return self.next( )

    ##
    # Closes the file.
    #
    def close(self):
        if self.handle:
            cplinkio.close( self.handle )
            self.handle = None

    ##
    # Transposes the file.
    #
    def transpose(self, new_path):
        return cplinkio.transpose( self.path, new_path )



class PlinkBinaryToVariants:
    """
    Class to write a PLINK BED format genotype dataset (.bed, .fam and .bim files)
    to vtools compatible variant format 
    
    c.f., http://pngu.mgh.harvard.edu/~purcell/plink/binary.shtml

    @param    dataset: path + prefix for a .bed, .fam and .bim without extensions
    @param      hgref: path to human genome reference file, e.g., build37.crr

    Note: implementation of this class is based on libplinkio
    https://bitbucket.org/mattias_franberg/libplinkio

    The SNP will be encoded as follows:
    * 0 - Homozygous major
    * 1 - Hetrozygous
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
    def __init__(self, dataset, hgref, logger):
        # check file path
        for ext in ['.fam', '.bed', '.bim']:
            if not os.path.exists(dataset + ext):
                raise RuntimeError('Cannot find file {0}'.format(dataset + ext))
        self.dataset = dataset
        self.cur = pio.open(self.dataset)
        self.logger = logger
        # a list of sample names (sample ID's in .fam file)'''
        self.samples =  [x.iid for x in self.cur.get_samples()]
        # iterator for variants info: chr, pos, allele1, allele2
        self.variants = it.chain(self.cur.get_loci())
        # reference genome object and ATCG dictionary
        self.hgref = cga.CrrFile(hgref)
        self.CSTRANDS = {'A':'T',
                         'G':'C',
                         'T':'A',
                         'C':'G'}
        # status 0 for not flip, 1 for flip, -1 for bad match
        self.status = 0

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
        # 0 based chrom, 0 based position
        # FIXME is it correct to do chrom - 1 and pos - 1 for cga input?
        ref = self.hgref.getBase(cga.Location(chrom - 1, pos - 1))
        self.status, strand, allele1, allele2 = self._matchref(ref, allele1, allele2)
        if self.status < 0:
            self.logger.info('Invalid locus {0}:{1} (given allele1 is {2}<->{3}., allele2 is {4}<->{5} but reference is {6})\n'.format(chrom, pos, self.CSTRANDS[allele1], allele1, self.CSTRANDS[allele2], allele2, ref))
            return None
        elif self.status == 0:
            if strand: self.logger.info('Use alternative strand at {0}:{1}\n'.format(chrom, pos))
            return ', '.join([str(chrom), str(pos), allele1, allele2]) + ', ' + str(geno_cur)[1:-1]
        else:
            # have to flip the genotypes coding
            if strand: self.logger.info('Use alternative strand for {0}:{1}\n'.format(chrom, pos))
            # self.logger.info('Allele coding flipped for {0}:{1}\n'.format(chrom, pos))
            # Very time consuming compare to not flipping the genotype codes
            return ', '.join([str(chrom), str(pos), allele2, allele1]) + ', ' + \
                ', '.join([str(x) if x == 3 or x == 'E' else str(2 - x) for x in geno_cur])
            
    def getHeader(self):
        '''a line of headers for the output text file'''
        return ', '.join(
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
        else:
            if which_major == 1:
                # allele 1 is indeed the major allele
                return True, self.getValidatedLocusGenotype(int(locus.chromosome), int(locus.bp_position),
                                         locus.allele1.upper(), locus.allele2.upper(),
                                         genotypes)
            else:
                return True, self.getValidatedLocusGenotype(int(locus.chromosome), int(locus.bp_position),
                                         locus.allele2.upper(), locus.allele1.upper(),
                                         genotypes)


    def determineMajorAllele(self, n=1000):
        '''The logic here is that for the first n loci we
        - assume the allele1 in *.bim is major allele
        - try to map to hg reference and see how many genotype coding has to be flipped
        - if too many (over half) have to be flipped we conclude that allele2 should be major allele
        @return: a guess of major allele based on the first n samples. -1 for bad matching
        '''
        # new temporary connection
        cur = pio.open(self.dataset)
        variants = it.chain(cur.get_loci())
        # counters
        m_ones = zeros = ones = 0
        for i in range(n):
            self.status = 0
            # self.status will be updated
            self.getLine(viter = variants, giter = cur)
            if self.status == 0:
                # not flipped
                zeros += 1
            elif self.status == 1:
                # flipped
                ones += 1
            else:
                # bad match
                m_ones += 1
        # check if so many are negative values
        if m_ones > float(n) / 2.0:
            return -1
        if ones > float(n - m_ones) / 2.0:
            # allele 2 seems to be major allele
            return 2
        else:
            return 1 
        

    def _matchref(self, ref, major, minor):
        '''try best to match reference allele
        @return self.status (0 for no need to flip, 1 for having to flip)
        @return strand (0 for original, 1 for alternative),
        @return major and minor alleles (might be from alternative strand)
        '''
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

class Logger:
    def __init__(self):
        self.info = sys.stderr.write
        self.debug = sys.stderr.write

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='This is a standalone binary PLINK to text format converter. For specification of the output text file please read the *.fmt file distributed with this script.')
    parser.add_argument('data', help = 'prefix of the PLINK binary files data-set')
    parser.add_argument('-b', '--build', choices = ['hg19','hg18'], help = 'human genome build')
    args = parser.parse_args()
    logger = Logger()
    if args.build == 'hg19':
        hgref = '~/.variant_tools/ftp.completegenomics.com/ReferenceFiles/build37.crr'
    else:
        hgref = '~/.variant_tools/ftp.completegenomics.com/ReferenceFiles/build36.crr'
    hgref = os.path.expanduser(hgref)
    if not os.path.exists(hgref):
        sys.exit('Cannot find file {0}'.format(hgref))
    p2v = PlinkBinaryToVariants(os.path.expanduser(args.data), hgref, logger)
    # check major allele
    n = 1000
    logger.info('checking first {0} loci\n'.format(n))
    which_major = p2v.determineMajorAllele(1000)
    logger.info('check status {0}\n'.format(which_major))
    # quit on bad match
    if which_major == -1:
        sys.exit('Too many unmatched loci: perhaps wrong hg reference is used?')
    # output
    print(p2v.getHeader())
    while True:
        flag, line = p2v.getLine(which_major = which_major)
        if not flag:
            break
        else:
            print(line)
