
import os
from variant_tools.pipeline import PipelineAction
from variant_tools.utils import env, openFile
        
class CountMappedReads(PipelineAction):
    '''This action reads the input files in sam format and count the total
    number of reads and number of mapped reads. It raises a RuntimeError
    if the proportion of unmapped reads lower than the specified cutoff value
    (default to 0.8). This action writes a count file to specified output
    files and read from this file directly if the file already exists.'''
    def __init__(self, output, cutoff=0.2):
        self.cutoff = cutoff
        PipelineAction.__init__(self, cmd='CountMappedReads', output=output)

    def _execute(self, sam_file, pipeline=None):
        #
        # count total reads and unmapped reads
        #
        unmapped_count = 0
        total_count = 0
        with open(sam_file[0]) as sam:
           for line in sam:
               total_count += 1
               if 'XT:A:N' in line:
                   unmapped_count += 1
        mapped_count = total_count - unmapped_count
        with open(self.output[0], 'w') as target:
            target.write('{}\n{}\n'.format(mapped_count, total_count))
        if total_count == 0 or mapped_count * 1.0 / total_count < self.cutoff:
            raise RuntimeError('{}: {} out of {} reads ({:.2f}%) are mapped.'
                .format(sam_file[0], mapped_count, total_count,
                    0 if total_count == 0 else (100.*mapped_count/total_count)))
        else:
            env.logger.info('{}: {} out of {} reads ({:.2f}%) are mapped.'
                .format(sam_file[0], mapped_count, total_count,
                    0 if total_count == 0 else (100.*mapped_count/total_count)))
        return self.output


    
class GuessReadGroup(PipelineAction):
    def __init__(self, bamfile, rgfile):
        self.output_bam = bamfile
        PipelineAction.__init__(self, cmd='GuessReadGroup', output=rgfile)

    def __call__(self, fastq_filename, pipeline=None):
        '''Get read group information from names of fastq files.'''
        # Extract read group information from filename such as
        # GERALD_18-09-2011_p-illumina.8_s_8_1_sequence.txt. The files are named 
        # according to the lane that produced them and whether they
        # are paired or not: Single-end reads s_1_sequence.txt for lane 1;
        # s_2_sequence.txt for lane 2 Paired-end reads s_1_1_sequence.txt 
        # for lane 1, pair 1; s_1_2_sequence.txt for lane 1, pair 2
        #
        # This function return a read group string like '@RG\tID:foo\tSM:bar'
        #
        # ID* Read group identifier. Each @RG line must have a unique ID. The
        # value of ID is used in the RG
        #     tags of alignment records. Must be unique among all read groups
        #     in header section. Read group
        #     IDs may be modifid when merging SAM fies in order to handle collisions.
        # CN Name of sequencing center producing the read.
        # DS Description.
        # DT Date the run was produced (ISO8601 date or date/time).
        # FO Flow order. The array of nucleotide bases that correspond to the
        #     nucleotides used for each
        #     flow of each read. Multi-base flows are encoded in IUPAC format, 
        #     and non-nucleotide flows by
        #     various other characters. Format: /\*|[ACMGRSVTWYHKDBN]+/
        # KS The array of nucleotide bases that correspond to the key sequence
        #     of each read.
        # LB Library.
        # PG Programs used for processing the read group.
        # PI Predicted median insert size.
        # PL Platform/technology used to produce the reads. Valid values: 
        #     CAPILLARY, LS454, ILLUMINA,
        #     SOLID, HELICOS, IONTORRENT and PACBIO.
        # PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for
        #     SOLiD). Unique identifier.
        # SM Sample. Use pool name where a pool is being sequenced.
        #
        with open(self.output[0], 'w') as rg_output:
            filename = os.path.basename(fastq_filename[0])
            output = os.path.basename(self.output_bam)
            # sample name is obtained from output filename without file extension
            SM = output.split('.', 1)[0]
            # always assume ILLUMINA for this script and BWA for processing
            PL = 'ILLUMINA'  
            PG = 'BWA'
            #
            # PU is for flowcell and lane information, ID should be unique for each
            #     readgroup
            # ID is temporarily obtained from input filename without exteion
            ID = filename.split('.')[0]
            # try to get lan information from s_x_1/2 pattern
            try:
                PU = re.search('s_([^_]+)_', filename).group(1)
            except AttributeError:
                env.logger.warning('Failed to guess lane information from filename {}'
                    .format(filename))
                PU = 'NA'
            # try to get some better ID
            try:
                # match GERALD_18-09-2011_p-illumina.8_s_8_1_sequence.txt
                m = re.match('([^_]+)_([^_]+)_([^_]+)_s_([^_]+)_([^_]+)_sequence.txt', filename)
                ID = '{}.{}'.format(m.group(1), m.group(4))
            except AttributeError as e:
                env.logger.warning('Input fasta filename {} does not match a known'
                    ' pattern. ID is directly obtained from filename.'.format(filename))
            #
            rg = r'@RG\tID:{}\tPG:{}\tPL:{}\tPU:{}\tSM:{}'.format(ID, PG, PL, PU, SM)
            env.logger.info('Setting read group tag to {}'.format(rg))
            rg_output.write(rg)
        return self.output

 
class CheckFastqVersion(PipelineAction):
    def __init__(self, output):
        PipelineAction.__init__(self, 'CheckFastqVersion', output)

    def __call__(self, fastq_file, pipeline=None):
        '''Detect the version of input fastq file. This can be very inaccurate'''
        if not os.path.isfile(fastq_file[0]) and os.path.isfile(fastq_file[0] + '.file_info'):
            if os.path.isfile(self.output[0]):
                return self.output
            else:
                raise RuntimeError('A valid fastq file is needed to check version of fastq: .file_info detected')
        with open(self.output[0], 'w') as aln_param:
            #
            # This function assumes each read take 4 lines, and the last line contains
            # quality code. It collects about 1000 quality code and check their range,
            # and use it to determine if it is Illumina 1.3+
            #
            qual_scores = ''
            with openFile(fastq_file[0]) as fastq:
                while len(qual_scores) < 1000:
                    try:
                        line = fastq.readline().decode('utf-8')
                    except Exception as e:
                        raise RuntimeError('Failed to read fastq file {}: {}'
                            .format(fastq_file, e))
                    if not line.startswith('@'):
                        raise ValueError('Wrong FASTA file {}'.format(fastq_file))
                    line = fastq.readline().decode('utf-8')
                    line = fastq.readline().decode('utf-8')
                    if not line.startswith('+'):
                        env.logger.warning(
                            'Suspiciout FASTA file {}: third line does not start with "+".'
                            .foramt(fastq_file))
                        return 
                    line = fastq.readline().decode('utf-8')
                    qual_scores += line.strip()
            #
            min_qual = min([ord(x) for x in qual_scores])
            max_qual = max([ord(x) for x in qual_scores])
            env.logger.debug('FASTA file with quality score ranging {} to {}'
                .format(min_qual, max_qual))
            # Sanger qual score has range Phred+33, so 33, 73 with typical score range 0 - 40
            # Illumina qual scores has range Phred+64, which is 64 - 104 with typical score range 0 - 40
            if min_qual >= 64 or max_qual > 90:
                # option -I is needed for bwa if the input is Illumina 1.3+ read format (quliaty equals ASCII-64).
                aln_param.write('-I')
        return self.output



