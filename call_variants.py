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
import time

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
    # UTILITY FUNCTIONS
    #
    def checkCmd(self, cmd):
        '''Check if a cmd exist'''
        if not hasattr(shutil, 'which'):
            env.logger.error('Please use Python 3.3 or higher for the use of shutil.which function')
            sys.exit(1)
        if shutil.which(cmd) is None:
            env.logger.error('Command {} does not exist. Please install it and try again.'.format(cmd))
            sys.exit(1)
        
    def downloadFile(self, URL, dest, quiet=False):
        '''Download a file from URL and save to dest '''
        # for some strange reason, passing wget without shell=True can fail silently.
        env.logger.info('Downloading {}'.format(URL))
        if os.path.isfile(dest):
            env.logger.warning('Using existing downloaded file {}.'.format(dest))
            return dest
        p = subprocess.Popen('wget {} -O {} {}'.format('-q' if quiet else '', dest, URL), shell=True)
        ret = p.wait()
        if ret == 0 and os.path.isfile(dest):
            return dest
        else:
            try:
                os.remove(dest)
            except OSError:
                pass
            raise RuntimeError('Failed to download {} using wget'.format(URL))

    def call(self, cmd, upon_succ=None, wait=True):
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
                if self.poll() >= env.jobs:
                    time.sleep(5)
                else:
                    break
            # there is a slot, start running
            env.logger.info('Running {}'.format(cmd))
            proc = subprocess.Popen(cmd, shell=True)
            env.running_jobs.append([proc, cmd, upon_succ])

    def poll(self):
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


    def wait(self):
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
    
    def fastaVersion(self, fasta_file):
        '''Detect the version of input fasta file. This can be very inaccurate'''
        #
        # This function assumes each read take 4 lines, and the last line contains
        # quality code. It collects about 1000 quality code and check their range,
        # and use it to determine if it is Illumina 1.3+
        #
        qual_scores = ''
        with open(fasta_file) as fasta:
            while len(qual_scores) < 1000:
                line = fasta.readline()
                if not line.startswith('@'):
                    raise ValueError('Wrong FASTA file {}'.foramt(fasta_file))
                line = fasta.readline()
                line = fasta.readline()
                if not line.startswith('+'):
                    env.logger.warning('Suspiciout FASTA file {}: third line does not start with "+".'.foramt(fasta_file))
                    return 'Unknown'
                line = fasta.readline()
                qual_scores += line.strip()
        #
        min_qual = min([ord(x) for x in qual_scores])
        max_qual = max([ord(x) for x in qual_scores])
        env.logger.debug('FASTA file with quality score ranging {} to {}'.format(min_qual, max_qual))
        if min_qual >= 64:
            # option -I is needed for bwa if the input is Illumina 1.3+ read format (quliaty equals ASCII-64).
            return 'Illumina 1.3+'
        else:
            # no option is needed for bwa
            return 'Sanger'

    def decompress(self, filename, dest_dir=None):
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
        #
        # if it is a tar file
        if mode is not None:
            dest_files = []
            env.logger.info('Extracting fasta sequences from tar file {}'.format(filename))
            #
            # MOTE: open a compressed tar file can take a long time because it needs to scan
            # the whole file to determine its content
            #
            # create a temporary directory to avoid corrupted file due to interrupted decompress
            os.mkdir('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
            with tarfile.open(filename, mode) as tar:
                files = tar.getnames()
                for filename in files:
                    # if there is directory structure within tar file, decompress all to the current directory
                    dest_file = os.path.join( '.' if dest_dir is None else dest_dir, os.path.basename(filename))
                    dest_files.append(dest_file)
                    if os.path.isfile(dest_file):
                        env.logger.warning('Using existing extracted file {}'.format(dest_file))
                    else:
                        env.logger.info('Extracting {} to {}'.format(filename, dest_file))
                        tar.extract(filename, 'tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'))
                        # move to the top directory with the right name only after the file has been properly extracted
                        shutil.move(os.path.join('tmp' if dest_dir is None else os.path.join(dest_dir, 'tmp'), filename), dest_file)
            return dest_files
        # return source file if 
        return [filename]
       
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
            self.call('wget -r {}'.format(URL))
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
                self.decompress(gzipped_file, '.')
 
    def buildBWARefIndex(self, ref_file):
        '''Create BWA index for reference genome file'''
        # bwa index -p hg19bwaidx -a bwtsw wg.fa
        if os.path.isfile('bwaidx.amb'):
            env.logger.warning('Using existing bwa indexed sequence bwaidx.amb')
        else:
            self.checkCmd('bwa')
            self.call('bwa index -p bwaidx -a bwtsw {}'.format(ref_file))

    def buildSamToolsRefIndex(self, ref_file):
        '''Create index for reference genome used by samtools'''
        if os.path.isfile('{}.fai'.format(ref_file)):
            env.logger.warning('Using existing samtools sequence index {}.fai'.format(ref_file))
        else:
            self.checkCmd('samtools')
            self.call('samtools faidx {}'.format(ref_file))

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
        '''Decompress input files to get a list of fasta files'''
        filenames = []
        for filename in input_files:
            filenames.extend(self.decompress(filename, working_dir))
        return filenames

    def align(self):
        '''Align to the reference genome'''
        pass

    def callVariants(self, input_files, output):
        '''Call variants from a list of input files'''
        pass


class hg19_gatk_23(BaseVariantCaller):
    '''A variant caller that uses gatk resource package 2.3 to call variants
    from build hg19 of the human genome'''
    def __init__(self, resource_dir):
        self.pipeline = 'hg19_gatk_23'
        BaseVariantCaller.__init__(self, resource_dir, self.pipeline)

    def checkResource(self):
        '''Check if needed resource is available.'''
        saved_dir = os.getcwd()
        os.chdir(self.resource_dir)
        files = ['ucsc.hg19.fasta.gz', 'bwaidx.amb', 'ucsc.hg19.fasta.fai']
        if not all([os.path.isfile(x) for x in files]):
            sys.exit('''Resource does not exist. Please run "call_variants.py prepare_resource"
                befor you execute this command.''')
        for cmd in ['wget',     # to download resource
                    'bwa',      # alignment
                    'samtools'  # merge sam files
                    ]:
            self.checkCmd(cmd)
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
        if not output.endswith('.bam'):
            env.logger.error('Plase specify a .bam file in the --output parameter')
            sys.exit(1)
        if os.path.isfile(output):
            env.logger.warning('Using existing output file {}'.format(output))
            return
        working_dir = os.path.split(output)[0]
        env.logger.info('Setting working directory to {}'.format(working_dir))
        #
        fasta_files = self.getFastaFiles(input_files, working_dir)
        fasta_files.sort()
        for input_file in fasta_files:
            dest_file = '{}/{}.sai'.format(working_dir, os.path.basename(input_file))
            if os.path.isfile(dest_file):
                env.logger.warning('Using existing alignent index file {}'.format(dest_file))
            else:
                # input file should be in fasta format (-t 4 means 4 threads)
                fmt = self.fastaVersion(input_file)
                opt = ' -I ' if fmt == 'Illumina 1.3+' else ''
                self.call('bwa aln {} -t 4 {}/bwaidx {} > {}_tmp'.format(opt, self.resource_dir, 
                    input_file, dest_file), 
                    upon_succ=(os.rename, dest_file + '_tmp', dest_file), wait=False)
        # wait for all bwa aln jobs to be completed
        self.wait()
        #
        # generate .bam files for each pair of pairend reads
        paired = True
        if len(fasta_files) // 2 * 2 != len(fasta_files):
            env.logger.warning('Odd number of fasta files provided, not handled as paired end reads.')
            paired = False
        for idx in range(len(fasta_files)//2):
            f1 = fasta_files[2*idx]
            f2 = fasta_files[2*idx + 1]
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
        sam_files = []
        if paired:
            for idx in range(len(fasta_files)//2):
                f1 = fasta_files[2*idx]
                f2 = fasta_files[2*idx + 1]
                sam_file = '{}/{}_bwa.sam'.format(working_dir, os.path.basename(f1))
                if os.path.isfile(sam_file):
                    env.logger.warning('Using existing sam file {}'.format(sam_file))
                else:
                    self.call('bwa sampe {0}/bwaidx {1}/{2}.sai {1}/{3}.sai {4} {5} > {6}_tmp'.format(self.resource_dir, 
                        working_dir, os.path.basename(f1), os.path.basename(f2), f1, f2, sam_file),
                        upon_succ=(os.rename, sam_file + '_tmp', sam_file), wait=False)
                sam_files.append(sam_file)
            # wait for all jobs to be completed
            self.wait()
        else:
            for f in fasta_file:
                sam_file = '{}/{}_bwa.sam'.format(working_dir, os.path.basename(f))
                if os.path.isfile(sam_file):
                    env.logger.warning('Using existing sam file {}'.format(sam_file))
                else:
                    self.call('bwa samse {0}/bwaidx {1}/{2}.sai {3} > {4}_tmp'.format(self.resource_dir,
                        working_dir, os.path.basename(f), f, sam_file),
                        upon_succ=(os.rename, sam_file + '_tmp', sam_file), wait=False)
                sam_files.append(sam_file)
            # wait for all jobs to be completed
            self.wait()
        # 
        # convert sam to bam files
        for sam_file in sam_files:
            bam_file = sam_file[:-4] + '.bam'
            if os.path.isfile(bam_file):
                env.logger.warning('Using existing bam file {}'.format(bam_file))
            else:
                self.call('samtools view -bt {}/ucsc.hg19.fasta.fai {} > {}_tmp'.format(
                    self.resource_dir, sam_file, bam_file),
                    upon_succ=(os.rename, bam_file + '_tmp', bam_file), wait=False)
        # wait for all sam->bam jobs to be completed
        self.wait()
        #
        # merge sam files?
        if len(sam_files) > 1:
            self.call('samtools merge {} {}'.format(output, ' '.join([x[:-4] + '.bam' for x in sam_files]))) 
        else:
            shutil.copy(sam_files[0], output)

    def callVariants(self, input_files, output):
        '''Call variants from a list of input files'''
        if os.path.isfile(output):
            env.logger.warning('Using existing output file {}'.format(output))
            return
        working_dir = os.path.split(output)[0]
        env.logger.info('Setting working directory to {}'.format(working_dir))


if __name__ == '__main__':
    #
    # vars
    master_parser = argparse.ArgumentParser(description='''A pipeline to call variants
        from raw sequence files, or single-sample bam files. It works (tested) only
        for Illumina sequence data, and for human genome with build hg19 of the
        reference genome. This pipeline uses BWA for alignment and GATK for variant
        calling. ''')
    master_parser.add_argument('--pipeline', nargs='?', default='hg19_gatk_23',
        choices=['hg19_gatk_23'],
        help='Name of the pipeline to be used to call variants.')
    #
    subparsers = master_parser.add_subparsers(title='Available operations', dest='action')
    #
    # action prepare_resource
    #
    resource = subparsers.add_parser('prepare_resource',
        help='Prepare resources for subsequent variant calling operations.',
        description='''This operation downloads GATK resource bundle and creates
            indexed reference genomes to be used by other tools.''')
    resource.add_argument('--resource_dir', default='~/.variant_tools/var_caller', 
        help='A directory for resources used by variant caller')
    #
    # action align
    #
    align = subparsers.add_parser('align',
        help='''Align raw reads to reference genome and return a compressed BAM file.
            The input files should be reads for the same sample, which could be individual
            fasta files, a tar file with all fasta files, or their gziped or bzipped
            versions. Filenames ending with _1 _2 will be considered as paired end reads.''')
    align.add_argument('input_files', nargs='+',
        help='''One or more .txt, .fa, .fasta, .tar, .tar.gz, .tar.bz2, .tbz2, .tgz files
            that contain raw reads of a single sample.''')
    align.add_argument('--output', required=True,
        help='''Output aligned reads to a sorted BAM file $output.bam. ''')
    align.add_argument('--resource_dir', default='~/.variant_tools/var_caller', 
        help='A directory for resources used by variant caller')
    align.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs.''')
    #
    # action call
    call = subparsers.add_parser('call',
        help='''Call variants from a list of BAM files.''')
    call.add_argument('input_files', nargs='+',
        help='''One or more BAM files.''')
    call.add_argument('--output', required=True,
        help='''Output called variants to the specified VCF file''')
    call.add_argument('-j', '--jobs', default=1, type=int,
        help='''Maximum number of concurrent jobs.''')
    #
    args = master_parser.parse_args()
    #
    if hasattr(args, 'output'):
        working_dir = os.path.split(args.output)[0]
        if not os.path.isdir(working_dir):
            os.makedirs(working_dir)
        # screen + log file logging
        env.logger = os.path.join(working_dir, 'call_variant.log')
    else:
        # screen only logging
        env.logger = None
    #
    # get a pipeline: args.pipeline is the name of the pipeline, also the name of the
    # class (subclass of VariantCaller) that implements the pipeline
    pipeline = eval(args.pipeline)(args.resource_dir)
    if args.action == 'prepare_resource':
        pipeline.prepareResourceIfNotExist()
    elif args.action == 'align':
        env.jobs = args.jobs
        pipeline.checkResource()
        pipeline.align(args.input_files, args.output)
    elif args.action == 'call':
        env.jobs = args.jobs
        pipeline.checkResource()
        pipeline.callVariants(args.input_files, args.output)

        
	
