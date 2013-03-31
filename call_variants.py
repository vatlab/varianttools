#!/usr/bin/env python3
#
import os
import sys
import subprocess
import glob
import argparse
import logging
import shutil
import tarfile
import copy
import gzip
import bz2
import zipfile
import time
from collections import defaultdict

#
# Runtime environment
#
class RuntimeEnvironment(object):
    '''Define the runtime environment of the pipeline, which presently provides
    maximum number of concurrent jobs, and a logger
    '''
    # the following make RuntimeEnvironment a singleton class
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(RuntimeEnvironment, cls).__new__(cls, *args, **kwargs)
        return cls._instance

    def __init__(self):
        self._max_jobs = 1
        self._logger = None
        #
        # additional parameters for args
        self.options = defaultdict(str)
        #
        # running_jobs implements a simple multi-processing queue system. Basically,
        # this variable holds the (popen, cmd, upon_succ) object for running jobs.
        #
        # Functions:
        #   BaseVariantCaller.call(cmd, upon_succ, wait=False)
        #      add an entry until there is less than self._max_jobs running jobs
        #   BaseVariantCaller.poll()
        #      check the number of running jobs
        #   BaseVariantCaller.wait()
        #      wait till all jobs are completed
        self.running_jobs = []
    
    #
    # max number of jobs
    #
    #
    def _setMaxJobs(self, x):
        try:
            self._max_jobs = int(x)
        except Exception as e:
            sys.exit('Failed to set max jobs: {}'.format(e))

    jobs = property(lambda self:self._max_jobs, _setMaxJobs)

    #
    # a global logger
    #
    class ColoredFormatter(logging.Formatter):
        ''' A logging format with colored output, which is copied from
        http://stackoverflow.com/questions/384076/how-can-i-make-the-python-logging-output-to-be-colored
        '''
        def __init__(self, msg):
            logging.Formatter.__init__(self, msg)
            self.LEVEL_COLOR = {
                'DEBUG': 'BLUE',
                'INFO': 'BLACK',
                'WARNING': 'PURPLE',
                'ERROR': 'RED',
                'CRITICAL': 'RED_BG',
                }
            self.COLOR_CODE={
                'ENDC':0,  # RESET COLOR
                'BOLD':1,
                'UNDERLINE':4,
                'BLINK':5,
                'INVERT':7,
                'CONCEALD':8,
                'STRIKE':9,
                'GREY30':90,
                'GREY40':2,
                'GREY65':37,
                'GREY70':97,
                'GREY20_BG':40,
                'GREY33_BG':100,
                'GREY80_BG':47,
                'GREY93_BG':107,
                'DARK_RED':31,
                'RED':91,
                'RED_BG':41,
                'LIGHT_RED_BG':101,
                'DARK_YELLOW':33,
                'YELLOW':93,
                'YELLOW_BG':43,
                'LIGHT_YELLOW_BG':103,
                'DARK_BLUE':34,
                'BLUE':94,
                'BLUE_BG':44,
                'LIGHT_BLUE_BG':104,
                'DARK_MAGENTA':35,
                'PURPLE':95,
                'MAGENTA_BG':45,
                'LIGHT_PURPLE_BG':105,
                'DARK_CYAN':36,
                'AUQA':96,
                'CYAN_BG':46,
                'LIGHT_AUQA_BG':106,
                'DARK_GREEN':32,
                'GREEN':92,
                'GREEN_BG':42,
                'LIGHT_GREEN_BG':102,
                'BLACK':30,
            }

        def colorstr(self, astr, color):
            return '\033[{}m{}\033[{}m'.format(self.COLOR_CODE[color], astr, self.COLOR_CODE['ENDC'])

        def format(self, record):
            record = copy.copy(record)
            levelname = record.levelname
            if levelname in self.LEVEL_COLOR:
                record.levelname = self.colorstr(levelname, self.LEVEL_COLOR[levelname])
                record.name = self.colorstr(record.name, 'BOLD')
                record.msg = self.colorstr(record.msg, self.LEVEL_COLOR[levelname])
            return logging.Formatter.format(self, record)

    def _setLogger(self, logfile=None):
        '''Create a logger with colored console output, and a log file if a
        filename is provided.'''
        # create a logger
        self._logger = logging.getLogger()
        self._logger.setLevel(logging.DEBUG)
        # output to standard output
        cout = logging.StreamHandler()
        cout.setLevel(logging.INFO)
        cout.setFormatter(RuntimeEnvironment.ColoredFormatter('%(levelname)s: %(message)s'))
        self._logger.addHandler(cout)
        if logfile is not None:
            # output to a log file
            ch = logging.FileHandler(logfile, 'a')
            ch.setLevel(logging.DEBUG)
            ch.setFormatter(logging.Formatter('%(asctime)s: %(levelname)s: %(message)s'))
            self._logger.addHandler(ch)

    logger = property(lambda self: self._logger, _setLogger)

# create a runtime environment object
env = RuntimeEnvironment()


#
# A simple job management scheme 
#
def run_command(cmd, upon_succ=None, wait=True):
    '''Call an external command, raise an error if it fails.
    If upon_succ is specified, the specified function and parameters will be
    evalulated after the job has been completed successfully.
    '''
    if wait or env.jobs == 1:
        try:
            env.logger.info('Running {}'.format(cmd))
            retcode = subprocess.call(cmd, shell=True)
            if retcode < 0:
                sys.exit("Command {} was terminated by signal {}".format(cmd, -retcode))
            elif retcode > 0:
                sys.exit("Command {} returned {}".format(cmd, retcode))
        except OSError as e:
            sys.exit("Execution of command {} failed: {}".format(cmd, e))
        # everything is OK
        if upon_succ:
            # call the function (upon_succ) using others as parameters.
            upon_succ[0](*(upon_succ[1:]))
    else:
        # if do not wait, look for running jobs
        while True:
            # if there are enough jobs running, wait
            if running_jobs() >= env.jobs:
                time.sleep(5)
            else:
                break
        # there is a slot, start running
        env.logger.info('Running {}'.format(cmd))
        proc = subprocess.Popen(cmd, shell=True)
        env.running_jobs.append([proc, cmd, upon_succ])

def running_jobs():
    '''check the number of running jobs'''
    count = 0
    for idx, job in enumerate(env.running_jobs):
        ret = job[0].poll()
        if ret is None:  # still running
            count += 1
        elif ret != 0:
            raise RuntimeError('Job {} failed.'.format(job[1]))
        else:
            # finish up
            if job[2]:
                # call the upon_succ function
                job[2][0](*(job[2][1:]))
            env.running_jobs[idx] = None
    # remove all completed jobs and exit
    env.running_jobs = [x for x in env.running_jobs if x]
    return count


def wait_all():
    '''Wait for all pending jobs to complete'''
    if not env.running_jobs:
        return
    for job in env.running_jobs:
        ret = job[0].wait()
        if ret != 0:
            raise RuntimeError('Job {} failed.'.format(job[1]))
        # run the upon_succ function
        if job[2]:
            job[2][0](*(job[2][1:]))
    # all jobs are completed
    env.running_jobs = []

#
# Check the existence of commands
#
def checkCmd(cmd):
    '''Check if a cmd exist'''
    if not hasattr(shutil, 'which'):
        env.logger.error('Please use Python 3.3 or higher for the use of shutil.which function')
        sys.exit(1)
    if shutil.which(cmd) is None:
        env.logger.error('Command {} does not exist. Please install it and try again.'.format(cmd))
        sys.exit(1)

def checkPicard():
    '''Check if picard is available, set PICARD_PATH if the path is specified in CLASSPATH'''
    if env.options['PICARD_PATH']:
        if not os.path.isfile(os.path.join(os.path.expanduser(env.options['PICARD_PATH']), 'SortSam.jar')):
            env.logger.error('Specified PICARD_PATH {} does not contain picard jar files.'.format(env.options['PICARD_PATH']))
            sys.exit(1)
    elif 'CLASSPATH' in os.environ:
        if not any([os.path.isfile(os.path.join(os.path.expanduser(x), 'SortSam.jar')) for x in os.environ['CLASSPATH'].split(':')]):
            env.logger.error('$CLASSPATH ({}) does not contain a path that contain picard jar files.'.format(os.environ['CLASSPATH']))
            sys.exit(1)
        for x in os.environ['CLASSPATH'].split(':'):
            if os.path.isfile(os.path.join(os.path.expanduser(x), 'SortSam.jar')):
                env.logger.info('Using picard under {}'.format(x))
                env.options['PICARD_PATH'] = os.path.expanduser(x)
                break
    else:
        env.logger.error('Please either specify path to picard using option PICARD_PATH=path, or set it in environment variable $CLASSPATH.')
        sys.exit(1)

def checkGATK():
    '''Check if GATK is available, set GATK_PATH from CLASSPATH if the path is specified in CLASSPATH'''
    if env.options['GATK_PATH']:
        if not os.path.isfile(os.path.join(os.path.expanduser(env.options['GATK_PATH']), 'GenomeAnalysisTK.jar')):
            env.logger.error('Specified GATK_PATH {} does not contain GATK jar files.'.format(env.options['GATK_PATH']))
            sys.exit(1)
    elif 'CLASSPATH' in os.environ:
        if not any([os.path.isfile(os.path.join(os.path.expanduser(x), 'GenomeAnalysisTK.jar')) for x in os.environ['CLASSPATH'].split(':')]):
            env.logger.error('$CLASSPATH ({}) does not contain a path that contain GATK jar files.'.format(os.environ['CLASSPATH']))
            sys.exit(1)
        else:
            for x in os.environ['CLASSPATH'].split(':'):
                if os.path.isfile(os.path.join(os.path.expanduser(x), 'GenomeAnalysisTK.jar')):
                    env.logger.info('Using GATK under {}'.format(x))
                    env.options['GATK_PATH'] = os.path.expanduser(x)
                    break
    else:
        env.logger.error('Please either specify path to GATK using option GATK_PATH=path, or set it in environment variable $CLASSPATH.')
        sys.exit(1) 

#
# Utility functions
# 
def downloadFile(URL, dest, quiet=False):
    '''Download a file from URL and save to dest'''
    # for some strange reason, passing wget without shell=True can fail silently.
    env.logger.info('Downloading {}'.format(URL))
    if os.path.isfile(dest):
        env.logger.warning('Using existing downloaded file {}.'.format(dest))
        return dest
    p = subprocess.Popen('wget {} -O {}_tmp {}'.format('-q' if quiet else '', dest, URL), shell=True)
    ret = p.wait()
    if ret == 0 and os.path.isfile(dest + '_tmp'):
        os.rename(dest + '_tmp', dest)
        return dest
    else:
        try:
            os.remove(dest + '_tmp')
        except OSError:
            pass
        raise RuntimeError('Failed to download {} using wget'.format(URL))

def fastqVersion(fastq_file):
    '''Detect the version of input fastq file. This can be very inaccurate'''
    #
    # This function assumes each read take 4 lines, and the last line contains
    # quality code. It collects about 1000 quality code and check their range,
    # and use it to determine if it is Illumina 1.3+
    #
    qual_scores = ''
    with open(fastq_file) as fastq:
        while len(qual_scores) < 1000:
            line = fastq.readline()
            if not line.startswith('@'):
                raise ValueError('Wrong FASTA file {}'.foramt(fastq_file))
            line = fastq.readline()
            line = fastq.readline()
            if not line.startswith('+'):
                env.logger.warning('Suspiciout FASTA file {}: third line does not start with "+".'.foramt(fastq_file))
                return 'Unknown'
            line = fastq.readline()
            qual_scores += line.strip()
    #
    min_qual = min([ord(x) for x in qual_scores])
    max_qual = max([ord(x) for x in qual_scores])
    env.logger.debug('FASTA file with quality score ranging {} to {}'.format(min_qual, max_qual))
    # Sanger qual score has range Phred+33, so 33, 73 with typical score range 0 - 40
    # Illumina qual scores has range Phred+64, which is 64 - 104 with typical score range 0 - 40
    if min_qual >= 64 or max_qual > 90:
        # option -I is needed for bwa if the input is Illumina 1.3+ read format (quliaty equals ASCII-64).
        return 'Illumina 1.3+'
    else:
        # no option is needed for bwa
        return 'Sanger'

def decompress(filename, dest_dir=None):
    '''If the file ends in .tar.gz, .tar.bz2, .bz2, .gz, .tgz, .tbz2, decompress it to
    dest_dir (current directory if unspecified), and return a list of files. Uncompressed
    files will be returned untouched.'''
    mode = None
    if filename.lower().endswith('.tar.gz') or filename.lower().endswith('.tar.bz2'):
        mode = 'r:gz'
    elif filename.lower().endswith('.tbz2') or filename.lower().endswith('.tgz'):
        mode = 'r:bz2'
    elif filename.lower().endswith('.tar'):
        mode = 'r'
    elif filename.lower().endswith('.gz'):
        dest_file = os.path.join('.' if dest_dir is None else dest_dir, os.path.basename(filename)[:-3])
        if os.path.isfile(dest_file):
            env.logger.warning('Using existing decompressed file {}'.format(dest_file))
        else:
            env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
            with gzip.open(filename, 'rb') as input, open(dest_file + '_tmp', 'wb') as output:
                buffer = input.read(10000000)
                while buffer:
                    output.write(buffer)
                    buffer = input.read(10000000)
            # only rename the temporary file to the right one after finishing everything
            # this avoids corrupted files
            os.rename(dest_file + '_tmp', dest_file)
        return [dest_file]
    elif filename.lower().endswith('.bz2'):
        dest_file = os.path.join('.' if dest_dir is None else dest_dir, os.path.basename(filename)[:-4])
        if os.path.isfile(dest_file):
            env.logger.warning('Using existing decompressed file {}'.format(dest_file))
        else:
            env.logger.info('Decompressing {} to {}'.format(filename, dest_file))
            with bz2.open(filename, 'rb') as input, open(dest_file + '_tmp', 'wb') as output:
                buffer = input.read(10000000)
                while buffer:
                    output.write(buffer)
                    buffer = input.read(10000000)
            # only rename the temporary file to the right one after finishing everything
            # this avoids corrupted files
            os.rename(dest_file + '_tmp', dest_file)
        return [dest_file]
    elif filename.lower().endswith('.zip'):
        bundle = zipfile.ZipFile(filename)
        dest_dir = '.' if dest_dir is None else dest_dir
        bundle.extractall(dest_dir)
        env.logger.info('Decompressing {} to {}'.format(filename, dest_dir))
        return [os.path.join(dest_dir, name) for name in bundle.namelist()]
    #
    # if it is a tar file
    if mode is not None:
        env.logger.info('Extracting fastq sequences from tar file {}'.format(filename))
        #
        # MOTE: open a compressed tar file can take a long time because it needs to scan
        # the whole file to determine its content. I am therefore creating a manifest
        # file for the tar file in the dest_dir, and avoid re-opening when the tar file
        # is processed again.
        manifest = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(filename) + '.manifest')
        all_extracted = False
        dest_files = []
        if os.path.isfile(manifest):
            all_extracted = True
            for f in [x.strip() for x in open(manifest).readlines()]:
                dest_file = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(f))
                if os.path.isfile(dest_file):
                    dest_files.append(dest_file)
                    env.logger.warning('Using existing extracted file {}'.format(dest_file))
                else:
                    all_extracted = False
        #
        if all_extracted:
            return dest_files
        #
        # create a temporary directory to avoid corrupted file due to interrupted decompress
        try:
            os.mkdir('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
        except:
            # directory might already exist
            pass
        #
        dest_files = []
        with tarfile.open(filename, mode) as tar:
            files = tar.getnames()
            # save content to a manifest
            with open(manifest, 'w') as manifest:
                for f in files:
                    manifest.write(f + '\n')
            for f in files:
                # if there is directory structure within tar file, decompress all to the current directory
                dest_file = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(f))
                dest_files.append(dest_file)
                if os.path.isfile(dest_file):
                    env.logger.warning('Using existing extracted file {}'.format(dest_file))
                else:
                    env.logger.info('Extracting {} to {}'.format(f, dest_file))
                    tar.extract(f, 'tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
                    # move to the top directory with the right name only after the file has been properly extracted
                    shutil.move(os.path.join('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'), f), dest_file)
        return dest_files
    # return source file if 
    return [filename]
   

#  Variant Caller
#
class BaseVariantCaller:
    '''A vase variant caller that is supposed to provide most of the utility functions
    and common operations that are not specific to any pipeline.
    '''
    def __init__(self, resource_dir, pipeline):
        self.resource_dir = os.path.join(os.path.expanduser(resource_dir), pipeline)
        if not os.path.isdir(self.resource_dir):
            env.logger.info('Creating resource directory {}'.format(self.resource_dir))
            os.makedirs(self.resource_dir)

    #
    # PREPARE RESOURCE
    #
    def downloadGATKResourceBundle(self, URL, files):
        '''Utility function to download GATK resource bundle. If files specified by
        files already exist, use the existing downloaded files.'''
        #
        if all([os.path.isfile(x) for x in files]):
            env.logger.warning('Using existing GATK resource')
        else:
            run_command('wget -r {}'.format(URL))
            # walk into the directory and get everything to the top directory
            # this is because wget -r saves files under URL/file
            for root, dirs, files in os.walk('.'):
                for name in files:
                    shutil.move(os.path.join(root, name), name)
        #
        # decompress all .gz files
        for gzipped_file in [x for x in os.listdir('.') if x.endswith('.gz') and not x.endswith('tar.gz')]:
            if os.path.isfile(gzipped_file[:-3]):
                env.logger.warning('Using existing decompressed file {}'.format(gzipped_file[:-3]))
            else:
                decompress(gzipped_file, '.')

    def buildBWARefIndex(self, ref_file):
        '''Create BWA index for reference genome file'''
        # bwa index -p hg19bwaidx -a bwtsw wg.fa
        if os.path.isfile('bwaidx.amb'):
            env.logger.warning('Using existing bwa indexed sequence bwaidx.amb')
        else:
            checkCmd('bwa')
            run_command('bwa index {} -p bwaidx -a bwtsw {}'.format(env.options['OPT_BWA_INDEX'], ref_file))

    def buildSamToolsRefIndex(self, ref_file):
        '''Create index for reference genome used by samtools'''
        if os.path.isfile('{}.fai'.format(ref_file)):
            env.logger.warning('Using existing samtools sequence index {}.fai'.format(ref_file))
        else:
            checkCmd('samtools')
            run_command('samtools faidx {} {}'.format(env.options['OPT_SAMTOOLS_FAIDX'], ref_file))

    # interface
    def checkResource(self):
        '''Check if needed resource is available.'''
        pass

    def prepareResourceIfNotExist(self):
        '''Prepare all resources for the pipeline. This is pipeline dependent.'''
        pass

    #
    # align and create bam file
    #
    def getFastaFiles(self, input_files, working_dir):
        '''Decompress input files to get a list of fastq files'''
        filenames = []
        for filename in input_files:
            filenames.extend(decompress(filename, working_dir))
        filenames.sort()
        return filenames

    def readGroup(self, filename, output):
        '''Get read group information from names of fastq files.'''
        # Extract read group information from filename such as
        # GERALD_18-09-2011_p-illumina.8_s_8_1_sequence.txt. The files are named 
        # according to the lane that produced them and whether they
        # are paired or not: Single-end reads s_1_sequence.txt for lane 1;
        # s_2_sequence.txt for lane 2 Paired-end reads s_1_1_sequence.txt 
        # for lane 1, pair 1; s_1_2_sequence.txt for lane 1, pair 2
        #
        # This function return a read group string like ‘@RG\tID:foo\tSM:bar’
        #
        # ID* Read group identifier. Each @RG line must have a unique ID. The value of ID is used in the RG
        #     tags of alignment records. Must be unique among all read groups in header section. Read group
        #     IDs may be modifid when merging SAM fies in order to handle collisions.
        # CN Name of sequencing center producing the read.
        # DS Description.
        # DT Date the run was produced (ISO8601 date or date/time).
        # FO Flow order. The array of nucleotide bases that correspond to the nucleotides used for each
        #     flow of each read. Multi-base flows are encoded in IUPAC format, and non-nucleotide flows by
        #     various other characters. Format: /\*|[ACMGRSVTWYHKDBN]+/
        # KS The array of nucleotide bases that correspond to the key sequence of each read.
        # LB Library.
        # PG Programs used for processing the read group.
        # PI Predicted median insert size.
        # PL Platform/technology used to produce the reads. Valid values: CAPILLARY, LS454, ILLUMINA,
        #     SOLID, HELICOS, IONTORRENT and PACBIO.
        # PU Platform unit (e.g. flowcell-barcode.lane for Illumina or slide for SOLiD). Unique identifier.
        # SM Sample. Use pool name where a pool is being sequenced.
        #
        # sample name is obtained from output filename
        SM = os.path.basename(output).split('.')[0]
        PL = 'ILLUMINA'  # always assume illumina in this script
        PG = 'BWA'
        # I do not know what this platform is, just return some info
        if not filename.endswith('_sequence.txt'):
            env.logger.error('Sequence filename {} does not ends with sequence.txt. Read group is set arbitrarily.'.format(filename))
            ID = filename.split('.')[0]
            PU = ID
        else:
            pieces = fastq_filename.split('_')
            if len(pieces) in [6, 7]:
                # GERALD  18-09-2011  p-illumina.8  s  8  1  sequence.txt
                ID = '{}.{}'.format(pieces[0], pieces[5])
                PU = pieces[5]
        #
        rg = r'@RG\tID:{}\tPG:{}\tPL:{}:PU:{}\tSM:{}'.format(ID, PG, PL, PU, SM)
        env.logger.info('Setting read group tag to {}'.format(rg))
        return rg

    def bwa_aln(self, fastq_files, working_dir):
        '''Use bwa aln to process fastq files'''
        for input_file in fastq_files:
            dest_file = '{}/{}.sai'.format(working_dir, os.path.basename(input_file))
            if os.path.isfile(dest_file):
                env.logger.warning('Using existing alignent index file {}'.format(dest_file))
            else:
                # input file should be in fastq format (-t 4 means 4 threads)
                opt = ' -I ' if fastqVersion(input_file) == 'Illumina 1.3+' else ''
                run_command('bwa aln {} {} -t 4 {}/bwaidx {} > {}_tmp'.format(opt, env.options['OPT_BWA_ALN'], self.resource_dir, 
                    input_file, dest_file), 
                    upon_succ=(os.rename, dest_file + '_tmp', dest_file), wait=False)
        # wait for all bwa aln jobs to be completed
        wait_all()

    def bwa_sampe(self, fastq_files, working_dir):
        '''Use bwa sampe to generate aligned sam files for paird end reads'''
        sam_files = []
        for idx in range(len(fastq_files)//2):
            f1 = fastq_files[2*idx]
            f2 = fastq_files[2*idx + 1]
            rg = self.getReadGroup(f1, working_dir)
            sam_file = '{}/{}_bwa.sam'.format(working_dir, os.path.basename(f1))
            if os.path.isfile(sam_file):
                env.logger.warning('Using existing sam file {}'.format(sam_file))
            else:
                run_command('bwa sampe {0} -r {1} {2}/bwaidx {3}/{4}.sai {3}/{5}.sai {6} {7} > {8}_tmp'.format(
                    env.options['OPT_BWA_SAMPE'], rg, self.resource_dir, 
                    working_dir, os.path.basename(f1), os.path.basename(f2), f1, f2, sam_file),
                    upon_succ=(os.rename, sam_file + '_tmp', sam_file), wait=False)
            sam_files.append(sam_file)
        # wait for all jobs to be completed
        wait_all()
        return sam_files

    def bwa_samse(self, fastq_files, working_dir):
        '''Use bwa sampe to generate aligned sam files'''
        sam_files = []
        for f in fastq_file:
            sam_file = '{}/{}_bwa.sam'.format(working_dir, os.path.basename(f))
            rg = self.getReadGroup(f, working_dir)
            if os.path.isfile(sam_file):
                env.logger.warning('Using existing sam file {}'.format(sam_file))
            else:
                run_command('bwa samse {0} -r {1} {2}/bwaidx {3}/{4}.sai {5} > {6}_tmp'.format(
                    env.options['OPT_BWA_SAMSE'], rg, self.resource_dir,
                    working_dir, os.path.basename(f), f, sam_file),
                    upon_succ=(os.rename, sam_file + '_tmp', sam_file), wait=False)
            sam_files.append(sam_file)
        # wait for all jobs to be completed
        wait_all()
        return sam_files        

    def sam2bam(self, sam_files):
        '''Convert sam file to sorted bam files.'''
        bam_files = []
        for sam_file in sam_files:
            bam_file = sam_file[:-4] + '.bam'
            if os.path.isfile(bam_file):
                env.logger.warning('Using existing bam file {}'.format(bam_file))
            else:
                run_command('samtools view {} -bt {}/ucsc.hg19.fasta.fai {} > {}_tmp'.format(
                    env.options['OPT_SAMTOOLS_VIEW'], self.resource_dir, sam_file, bam_file),
                    upon_succ=(os.rename, bam_file + '_tmp', bam_file), wait=False)
            bam_files.append(bam_file)
        # wait for all sam->bam jobs to be completed
        wait_all()
        #
        # sort bam files
        sorted_bam_files = []
        for bam_file in bam_files:
            sorted_bam_file = bam_file[:-4] + '_sorted.bam'
            if os.path.isfile(sorted_bam_file):
                env.logger.warning('Using existing sorted bam file {}'.format(sorted_bam_file))
            else:
                run_command('samtools sort {} {} {}_tmp'.format(
                    env.options['OPT_SAMTOOLS_SORT'], bam_file, sorted_bam_file[:-4]),
                    upon_succ=(os.rename, sorted_bam_file[:-4] + '_tmp.bam', sorted_bam_file), wait=False)
            sorted_bam_files.append(sorted_bam_file)
        wait_all()
        return sorted_bam_files

    def mergeBAMs(self, bam_files, output):
        '''merge sam files'''
        # use Picard merge, not samtools merge: 
        # Picard keeps RG information from all Bam files, whereas samtools uses only 
        # inf from the first bam file
        if os.path.isfile(output):
            env.logger.warning('Using existing merged bam file {}'.format(output))
        else:
            run_command('''java {} -jar {}/MergeSamFiles.jar {} {} USE_THREADING=true
    	        VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true
    	        OUTPUT={}'''.format(env.options['OPT_JAVA'],
                    env.options['PICARD_PATH'], env.options['OPT_PICARD_MERGESAMFILES'],
                    ' '.join(['INPUT={}'.format(x) for x in bam_files]), output[:-4] + '_tmp.bam').replace('\n', ' '),
                upon_succ=(os.rename, output[:-4] + '_tmp.bam', output))

    def indexBAMs(self, bam_files):
        '''Index the input bam file'''
        for bam_file in bam_files:
            if os.path.isfile('{}.bai'.format(bam_file)):
                env.logger.warning('Using existing bam index {}.bai'.format(bam_file))
            else:
                run_command('samtools index {0} {1} {1}_tmp.bai'.format(
                    env.options['OPT_SAMTOOLS_INDEX'], bam_file),
                    upon_succ=(os.rename, bam_file + '_tmp.bai', bam_file + '.bai'),
                    wait=False)
        #
        wait_all()

    def align(self, input_files, output):
        '''Align to the reference genome'''
        if not output.endswith('.bam'):
            env.logger.error('Plase specify a .bam file in the --output parameter')
            sys.exit(1)
        if os.path.isfile(output) and os.path.isfile(output + '.bai'):
            env.logger.warning('Using existing output file {}'.format(output))
            sys.exit(0)

    def realignIndels(self, input_files):
        '''Create realigner target and realign indels'''
        targets = []
        for bam_file in input_files:
            target = os.path.basename(bam_file)[:-4] + '.IndelRealignerTarget.intervals'
            if os.path.isfile(target):
                env.logger.warning('Using existing realigner target {}'.format(target))
            else:
                run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                    -R {4}/ucsc.hg19.fasta
                    -T RealignerTargetCreator
                    --mismatchFraction 0.0
                    -known {4}/dbsnp_137.hg19.vcf
                    -known {4}/hapmap_3.3.hg19.vcf
                    -known {4}/Mills_and_1000G_gold_standard.indels.hg19.vcf
                    -known {4}/1000G_phase1.indels.hg19.vcf
                    -o {5}_tmp '''.format(env.options['OPT_JAVA'], env.options['GATK_PATH'],
                    env.options['OPT_GATK_REALIGNERTARGETCREATOR'], bam_file, self.resource_dir,
                    target).replace('\n', ' '), upon_succ=(os.rename, target + '_tmp', target), wait=False)
            targets.append(target)
        # wait for the completion of gatk RealignerTargetCreator
        wait_all()
        # 
        # realign around known indels
        cleaned_bam_files = []
        for bam_file, target in zip(input_files, targets):
            cleaned_bam_file = os.path.basename(bam_file)[:-4] + '.clean.bam'
            if os.path.isfile(cleaned_bam_file):
                env.logger.warning('Using existing realigner bam file {}'.format(cleaned_bam_file))
            else:
                run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                    -R {4}/ucsc.hg19.fasta
                    -T IndelRealigner 
                    --targetIntervals {5}
                    -known {4}/dbsnp_137.hg19.vcf
                    -known {4}/hapmap_3.3.hg19.vcf
                    -known {4}/Mills_and_1000G_gold_standard.indels.hg19.vcf
                    -known {4}/1000G_phase1.indels.hg19.vcf
                    --consensusDeterminationModel USE_READS
                    -compress 0 -o {6}'''.format(env.options['OPT_JAVA'], env.options['GATK_PATH'],
                    env.options['OPT_GATK_REALIGNERTARGETCREATOR'], bam_file, self.resource_dir,
                    target, cleaned_bam_file[:-4] + '_tmp.bam').replace('\n', ' '),
                    upon_succ=(os.rename, cleaned_bam_file[:-4] + '_tmp.bam', cleaned_bam_file), wait=False)
            cleaned_bam_files.append(cleaned_bam_file)
        # wait for the completion of gatk IndelRealigner
        wait_all()
        # 
        return cleaned_bam_files

    def markDuplicates(self, input_files):
        '''Mark duplicate using picard'''
        dedup_bam_files = []
        for bam_file in input_files:
            dedup_bam_file = os.path.basename(bam_file)[:-4] + '.dedup.bam'
            metrics_file = os.path.basename(bam_file)[:-4] + '.metrics'
            if os.path.isfile(dedup_bam_file):
                env.logger.warning('Using existing realigner target {}'.format(target))
            else:
                run_command('''java {0} -jar {1}/MarkDuplicates.jar {2}
                    INPUT={3}
                    OUTPUT={4}
                    METRICS_FILE={5}
                    VALIDATION_STRINGENCY=LENIENT
                    '''.format(env.options['OPT_JAVA'], env.options['PICARD_PATH'],
                            env.options['OPT_PICARD_MARKDUPLICATES'], bam_file, dedup_bam_file[:-4] + '_tmp.bam',
                            metrics_file).replace('\n', ' '), 
                    upon_succ=(os.rename, dedup_bam_file[:-4] + '_tmp.bam', dedup_bam_file), wait=False)
            dedup_bam_files.append(dedup_bam_file)
        # wait for the completion of picard markduplicates
        wait_all()
        # 
        return dedup_bam_files

    def recalibrate(self, input_files):
        '''Create realigner target and realign indels'''
        targets = []
        for bam_file in input_files:
            target = os.path.basename(bam_file)[:-4] + '.grp'
            if os.path.isfile(target):
                env.logger.warning('Using existing base recalibrator target {}'.format(target))
            else:
                run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                    -T BaseRecalibrator
                    -R {4}/ucsc.hg19.fasta
                    -cov ReadGroupCovariate
                    -cov QualityScoreCovariate
                    -cov CycleCovariate
                    -cov ContextCovariate
                    -knownSites {4}/dbsnp_137.hg19.vcf
                    -knownSites {4}/hapmap_3.3.hg19.vcf
                    -knownSites {4}/1000G_omni2.5.hg19.vcf
                    -knownSites {4}/Mills_and_1000G_gold_standard.indels.hg19.vcf
                    -knownSites {4}/1000G_phase1.indels.hg19.vcf
                    -o {5}_tmp '''.format(env.options['OPT_JAVA'], env.options['GATK_PATH'],
                    env.options['OPT_GATK_BASERECALIBRATOR'], bam_file, self.resource_dir,
                    target).replace('\n', ' '), upon_succ=(os.rename, target + '_tmp', target), wait=False)
            targets.append(target)
        # wait for the completion of gatk BaseRecalibrator
        wait_all()
        #
        # recalibrate
        recal_bam_files = []
        for bam_file, target in zip(input_files, targets):
            recal_bam_file = os.path.basename(bam_file)[:-4] + '.recal.bam'
            if os.path.isfile(recal_bam_file):
                env.logger.warning('Using existing recalibrated bam file {}'.format(recal_bam_file))
            else:
                run_command('''java {0} -jar {1}/GenomeAnalysisTK.jar {2} -I {3} 
                    -R {4}/ucsc.hg19.fasta
                    -T PrintReads
                    -BQSR {5}
                    -baq=CalculationMode.CALCULATE_AS_NECESSARY
                    -o {6}'''.format(env.options['OPT_JAVA'], env.options['GATK_PATH'],
                    env.options['OPT_GATK_PRINTREADS'], bam_file, self.resource_dir,
                    target, recal_bam_file[:-4] + '_tmp.bam').replace('\n', ' '),
                    upon_succ=(os.rename, recal_bam_file[:-4] + '_tmp.bam', recal_bam_file), wait=False)
            recal_bam_files.append(recal_bam_file)
        # wait for the completion of gatk printreads
        wait_all()
        # 
        return recal_bam_files

    def callVariants(self, input_files, output):
        '''Call variants from a list of input files'''
        if not output.endswith('.vcf'):
           env.logger.error('Please specify a .vcf file in the --output parameter')
           sys.exit(1)
        if os.path.isfile(output):
            env.logger.warning('Using existing output file {}'.format(output))
            sys.exit(0)
        for bam_file in input_files:
            if not os.path.isfile(bam_file):
                env.logger.error('Input file {} does not exist'.format(bam_file))
                sys.exit(1)
            if not os.path.isfile(bam_file + '.bai'):
                env.logger.error('Input bam file {} is not indexed.'.format(bam_file))
                sys.exit(1)


class hg19_gatk_23(BaseVariantCaller):
    '''A variant caller that uses gatk resource package 2.3 to call variants
    from build hg19 of the human genome'''
    def __init__(self, resource_dir):
        self.pipeline = 'hg19_gatk_23'
        BaseVariantCaller.__init__(self, resource_dir, self.pipeline)

    def checkResource(self):
        '''Check if needed resource is available. This pipeline requires
        GATK resource bundle, commands wget, bwa, samtools, picard, and GATK. '''
        saved_dir = os.getcwd()
        os.chdir(self.resource_dir)
        files = ['ucsc.hg19.fasta.gz', 'bwaidx.amb', 'ucsc.hg19.fasta.fai']
        if not all([os.path.isfile(x) for x in files]):
            sys.exit('''GATK resource bundle does not exist in directory {}. Please run "call_variants.py prepare_resource" befor you execute this command.'''.format(self.resource_dir))
        #
        for cmd in ['wget',     # to download resource
                    'bwa',      # alignment
                    'samtools'  # merge sam files
                    ]:
            checkCmd(cmd)
        #
        checkPicard()
        checkGATK()
        #
        os.chdir(saved_dir)

    def prepareResourceIfNotExist(self):
        '''This function downloads the UCSC resource boundle for specified build and
        creates bwa and samtools indexed files from the whole genome sequence '''
        saved_dir = os.getcwd()
        os.chdir(self.resource_dir)
        #
        # these are pipeline specific
        self.downloadGATKResourceBundle('ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.3/hg19/*', files=['ucsc.hg19.fasta.gz'])
        self.buildBWARefIndex('ucsc.hg19.fasta')
        self.buildSamToolsRefIndex('ucsc.hg19.fasta')
        # 
        os.chdir(saved_dir)

    def align(self, input_files, output):
        '''Align reads to hg19 reference genome and return a sorted, indexed bam file.'''
        BaseVariantCaller.align(self, input_files, output)
        #
        # the working dir is a directory under output, the middle files are saved to this
        # directory to avoid name conflict
        working_dir = os.path.join(os.path.split(output)[0], os.path.basename(output) + '_align_cache')
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        #
        env.logger.info('Setting working directory to {}'.format(working_dir))
        #
        # step 1: decompress to get a list of fastq files
        fastq_files = self.getFastaFiles(input_files, working_dir)
        #
        # step 2: call bwa aln to produce .sai files
        self.bwa_aln(fastq_files, working_dir)
        #
        # step 3: generate .sam files for each pair of pairend reads, or reach file of unpaired reads
        paired = True
        if len(fastq_files) // 2 * 2 != len(fastq_files):
            env.logger.warning('Odd number of fastq files provided, not handled as paired end reads.')
            paired = False
        for idx in range(len(fastq_files)//2):
            f1 = fastq_files[2*idx]
            f2 = fastq_files[2*idx + 1]
            if len(f1) != len(f2):
                env.logger.warning('Filenames {}, {} are not paired, not handled as paired end reads.'.format(f1, f2))
                paired = False
                break
            diff = [ord(y)-ord(x) for x,y in zip(f1, f2) if x!=y]
            if diff != [1]:
                env.logger.warning('Filenames {}, {} are not paired, not handled as paired end reads.'.format(f1, f2))
                paired = False
                break
        #
        # sam files?
        if paired:
            sam_files = self.bwa_sampe(fastq_files, working_dir)
        else:
            sam_files = self.bwa_samse(fastq_files, working_dir)
        # 
        # step 4: convert sam to sorted bam files
        bam_files = self.sam2bam(sam_files)
        #
        # step 5: merge sorted bam files to output file
        if len(bam_files) > 1:
            self.mergeBAMs(bam_files, output)
        else:
            shutil.copy(bam_files[0], output)
        #
        # step 6: index the output bam file
        self.indexBAMs([output])

    def callVariants(self, input_files, output):
        '''Call variants from a list of input files'''
        BaseVariantCaller.callVariants(self, input_files, output)
        working_dir = os.path.join(os.path.split(output)[0], os.path.basename(output) + '_call_cache')
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        #
        env.logger.info('Setting working directory to {}'.format(working_dir))
        #
        # step 1: create indel realignment targets
        cleaned_bam_files = self.realignIndels(input_files)
        self.indexBAMs(cleaned_bam_files)
        #
        # step 2: Mark duplicate
        dedup_bam_files = self.markDuplicates(cleaned_bam_files)
        self.indexBAMs(dedup_bam_files)
        #
        # step 3: recalibration
        recal_bam_files = self.recalibrate(dedup_bam_files)
        self.indexBAMs(recal_bam_files)
        

if __name__ == '__main__':
    #
    options = [
        ('PICARD_PATH', ''),
        ('GATK_PATH', ''),
        ('OPT_JAVA', '-Xmx4g'),
        ('OPT_BWA_INDEX', ''),
        ('OPT_SAMTOOLS_FAIDX', ''),
        ('OPT_BWA_ALN', ''),
        ('OPT_BWA_SAMPE', ''),
        ('OPT_BWA_SAMSE', ''),
        ('OPT_SAMTOOLS_VIEW', ''),
        ('OPT_SAMTOOLS_SORT', ''),
        ('OPT_SAMTOOLS_INDEX', ''),
        ('OPT_PICARD_MERGESAMFILES', ''),
        ('OPT_GATK_REALIGNERTARGETCREATOR', ''),
        ('OPT_GATK_INDELREALIGNER', ''),
        ('OPT_PICARD_MARKDUPLICATES', ''),
        ('OPT_GATK_BASERECALIBRATOR', ''),
        ('OPT_GATK_PRINTREADS', ''),
        ]
    def addCommonArguments(parser, args):
        if 'pipeline' in args:
            parser.add_argument('--pipeline', nargs='?', default='hg19_gatk_23',
                choices=['hg19_gatk_23'],
                help='Name of the pipeline to be used to call variants.')
        if 'resource_dir' in args:
            parser.add_argument('--resource_dir', default='~/.variant_tools/var_caller', 
                help='''A directory for resources used by variant caller. Default to
                    ~/.variant_tools/var_caller.''')
        if 'set' in args:
            parser.add_argument('--set', nargs='*',
                help='''Set runtime variables in the format of NAME=value. NAME can be
                    PICARD_PATH (path to picard, should have a number of .jar files 
                    under it), GATK_PATH (path to gatk, should have GenomeAnalysisTK.jar
                    under it) for path to JAR files (can be ignored if the paths are
                    specified in environment variable $CLASSPATH), additional options
                    to command java (OPT_JAVA, parameter to the java command, default
                    to value "-Xmx4g"), and options to individual subcommands such as
                    OPT_BWA_INDEX (additional option to bwa index) and OPT_SAMTOOLS_FAIDX.
                    The following options are acceptable: {}'''.format(', '.join(
                    [x[0] for x in options])))
        if 'jobs' in args:
            parser.add_argument('-j', '--jobs', default=1, type=int,
                help='''Maximum number of concurrent jobs.''')
    #
    master_parser = argparse.ArgumentParser(description='''A pipeline to call variants
        from raw sequence files, or single-sample bam files. It works (tested) only
        for Illumina sequence data, and for human genome with build hg19 of the
        reference genome. This pipeline uses BWA for alignment and GATK for variant
        calling.''')

    subparsers = master_parser.add_subparsers(title='Available operations', dest='action')
    #
    # action prepare_resource
    #
    resource = subparsers.add_parser('prepare_resource',
        help='Prepare resources for subsequent variant calling operations.',
        description='''This operation downloads GATK resource bundle and creates
            indexed reference genomes to be used by other tools.''')
    addCommonArguments(resource, ['pipeline', 'resource_dir'])
    #
    # action align
    #
    align = subparsers.add_parser('align',
        help='''Align raw reads to reference genome and return a compressed BAM file.
            The input files should be reads for the same sample, which could be individual
            fastq files, a tar file with all fastq files, or their gziped or bzipped
            versions. Filenames ending with _1 _2 will be considered as paired end reads.''')
    align.add_argument('input_files', nargs='+',
        help='''One or more .txt, .fa, .fastq, .tar, .tar.gz, .tar.bz2, .tbz2, .tgz files
            that contain raw reads of a single sample.''')
    align.add_argument('--output', required=True,
        help='''Output aligned reads to a sorted BAM file $output.bam. ''')
    addCommonArguments(align, ['pipeline', 'resource_dir', 'set', 'jobs'])
    #
    # action call
    call = subparsers.add_parser('call',
        help='''Call variants from a list of BAM files.''')
    call.add_argument('input_files', nargs='+',
        help='''One or more BAM files.''')
    call.add_argument('--output', required=True,
        help='''Output called variants to the specified VCF file''')
    addCommonArguments(call, ['pipeline', 'resource_dir', 'set', 'jobs'])
    #
    args = master_parser.parse_args()
    #
    if hasattr(args, 'output'):
        working_dir = os.path.split(args.output)[0]
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        # screen + log file logging
        env.logger = os.path.join(working_dir, os.path.basename(args.output) + '.log')
    else:
        # screen only logging
        env.logger = None
    #
    # handling additional parameters
    # set default value
    for opt in options:
        if opt[1]:
            env.options[opt[0]] = opt[1]
    # override using command line values
    if hasattr(args, 'set') and args.set is not None:
        for arg in args.set:
            if '=' not in arg:
                sys.exit('Additional parameter should have form NAME=value')
            name, value = arg.split('=', 1)
            if name not in [x[0] for x in options]:
                env.logger.error('Unrecognized environment variable {}: {} are allowed.'.format(
                    name, ', '.join([x[0] for x in options])))
                sys.exit(1)
            env.options[name] = value
            env.logger.info('Environment variable {} is set to {}'.format(name, value))
    #
    # get a pipeline: args.pipeline is the name of the pipeline, also the name of the
    # class (subclass of VariantCaller) that implements the pipeline
    pipeline = eval(args.pipeline)(args.resource_dir)
    if args.action == 'prepare_resource':
        pipeline.prepareResourceIfNotExist()
    elif args.action == 'align':
        env.jobs = args.jobs
        checkPicard()
        pipeline.checkResource()
        pipeline.align(args.input_files, args.output)
    elif args.action == 'call':
        env.jobs = args.jobs
        checkPicard()
        checkGATK()
        pipeline.checkResource()
        pipeline.callVariants(args.input_files, args.output)

        
	
