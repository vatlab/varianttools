import sys, os
import itertools as it
from .utils import ProgressBar, RefGenome
from .plinkfile import PlinkFile

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
    def __init__(self, dataset, build, logger):
        # check file path
        for ext in ['.fam', '.bed', '.bim']:
            if not os.path.exists(dataset + ext):
                raise RuntimeError('Cannot find file {0}'.format(dataset + ext))
        self.dataset = dataset
        self.build = build
        self.cur = PlinkFile(self.dataset)
        self.logger = logger
        # a list of sample names (sample ID's in .fam file)'''
        self.samples =  [x.iid for x in self.cur.get_samples()]
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
        try:
            ref = self.hgref.getBase(chrom, pos)
        except Exception as e:
            self.logger.warning('Cannot find locus {0}:{1}. Input variant is ignored'.format(chrom, pos))
            self.status = -1
            self.cur.close()
            return None
        self.status, strand, allele1, allele2 = self._matchref(ref, allele1, allele2)
        if self.status < 0:
            self.logger.warning('Invalid locus {0}:{1} (given allele1 is {2}<->{3}., allele2 is {4}<->{5} but reference is {6})'.\
                                    format(chrom, pos, self.CSTRANDS[allele1], allele1, self.CSTRANDS[allele2], allele2, ref))
            return None
        elif self.status == 0:
            if strand:
                self.logger.debug('Use alternative strand for {0}:{1}'.format(chrom, pos))
            return ','.join([str(chrom), str(pos), allele1, allele2]) + ',' + str(geno_cur)[1:-1]
        else:
            # have to flip the genotypes coding
            if strand:
                self.logger.debug('Use alternative strand for {0}:{1}'.format(chrom, pos))
            # self.logger.debug('Allele coding flipped for {0}:{1}'.format(chrom, pos))
            # Very time consuming compare to not flipping the genotype codes
            return ','.join([str(chrom), str(pos), allele2, allele1]) + ',' + \
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
        else:
            if which_major == 1:
                # allele 1 is the major allele
                return True, self.getValidatedLocusGenotype(int(locus.chromosome), int(locus.bp_position),
                                         locus.allele1.upper(), locus.allele2.upper(),
                                         genotypes)
            else:
                # allele 2 is the major allele
                return True, self.getValidatedLocusGenotype(int(locus.chromosome), int(locus.bp_position),
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

    def convert(self, files, outdir_format, logger = None):
        for item in files:
            ofile = self._get_output_fn(item, outdir_format)
            if logger:
                logger.info('Convert {} to {}'.format(item, ofile))
        

class PlinkConverter(Preprocessor):
    def __init__(self, build):
        Preprocessor.__init__(self)
        self.build = build

    def convert(self, files, output_files, logger):
        for item, ofile in zip(files, output_files):
            if os.path.exists(item + ".bed"):
                self.decode_plink(PlinkBinaryToVariants(item, self.build, logger), ofile, logger=logger)
            else:
                import glob
                files = '/'.join([x for x in glob.glob(item + '*')])
                if files:
                    supported = ['*.bed/*.bim/*.fam']
                    raise ValueError("Unsupported input file '{}' (supported file types are {})".\
                                         format(files, ';'.join(supported)))
                else:
                    raise ValueError("Cannot find input files '{}'".format(item + '*'))
            
    def decode_plink(self, p2vObject, ofile, n = 1000, logger = None):
        '''decode plink data from p2vObject and output to ofile'''
        if logger: logger.info("Determining major/minor allele from data")
        # check major allele
        which_major = p2vObject.determineMajorAllele(n)
        # raise on bad match
        if which_major == -1:
            raise ValueError ('Invalid dataset {0}: too many unmatched loci to {1}. Perhaps wrong reference genome is used?'.\
                                  format(p2vObject.dataset, p2vObject.build))
        if logger: logger.debug("allele{} is major allele".format(which_major))
        # output
        nloci = p2vObject.getLociCounts()
        batch = int(nloci / 100)
        prog = ProgressBar('Decoding {0}'.format(p2vObject.dataset), nloci)
        if os.path.exists(ofile):
            os.remove(ofile)
        with open(ofile, 'a') as f:
            f.write(p2vObject.getHeader() + '\n')
            count = 0
            while True:
                flag, line = p2vObject.getLine(which_major = which_major)
                count += 1
                if not flag:
                    prog.done()
                    break
                else:
                    if line is not None:
                        f.write(line + '\n')
                    if count % batch == 0 and count > batch:
                        prog.update(count)
