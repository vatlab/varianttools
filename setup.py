#!/usr/bin/env python2.7
#
# $File: setup.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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
#
from distutils.core import setup, Extension
from distutils.ccompiler import new_compiler

try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

# parallel compilation
import multiprocessing, multiprocessing.pool

def compile_parallel(
        self,
        sources,
        output_dir=None,
        macros=None,
        include_dirs=None,
        debug=0,
        extra_preargs=None,
        extra_postargs=None,
        depends=None):

    # Copied from distutils.ccompiler.CCompiler
    macros, objects, extra_postargs, pp_opts, build = self._setup_compile(
        output_dir, macros, include_dirs, sources, depends, extra_postargs)
    cc_args = self._get_cc_args(pp_opts, debug, extra_preargs)
    #
    def _single_compile(obj):

        try:
            src, ext = build[obj]
        except KeyError:
            return
        self._compile(obj, src, ext, cc_args, extra_postargs, pp_opts)
    # convert to list, imap is evaluated on-demand
    list(multiprocessing.pool.ThreadPool(multiprocessing.cpu_count()).imap(_single_compile, objects))
    return objects

import distutils.ccompiler
distutils.ccompiler.CCompiler.compile=compile_parallel
#
import sys, os, subprocess
# use ccache to speed up build
# try:
#     if subprocess.call(['ccache'], stderr = open(os.devnull, "w")):
#         os.environ['CC'] = 'ccache gcc'
# except OSError:
#     pass

#
# if building source package, we will need to have wrapper files for both
# versions of Python
#
SDIST = 'sdist' in sys.argv

try:
    import argparse
except ImportError:
    sys.exit('variant tools requires Python 2.7.2 or higher, or Python 3.2 or higher. Please upgrade your version (%s) of Python and try again.' % (sys.version.split()[0]))

# do not import variant_tools because __init__ might not be imported properly 
# before installation
with open('variant_tools/__init__.py') as init:
    for line in init:
        if line.startswith('VTOOLS_VERSION='):
            VTOOLS_VERSION = line[15:].strip().strip('"').strip("'")
            break

EMBEDDED_BOOST = os.path.isdir('boost_1_49_0')
SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword', '-w-511', '-w-509',
    '-outdir', 'variant_tools', '-py3']

SQLITE_FOLDER = 'sqlite/py3'
ASSO_WRAPPER_CPP_FILE = 'variant_tools/assoTests_wrap.cpp'
ASSO_WRAPPER_PY_FILE = 'variant_tools/assoTests.py'
ASSO_INTERFACE_FILE = 'variant_tools/assoTests.i'
CGATOOLS_WRAPPER_CPP_FILE = 'variant_tools/cgatools_wrap.cpp'
CGATOOLS_WRAPPER_PY_FILE = 'variant_tools/cgatools.py'
CGATOOLS_INTERFACE_FILE = 'variant_tools/cgatools.i'
UCSCTOOLS_WRAPPER_CPP_FILE = 'variant_tools/ucsctools_wrap.cpp'
UCSCTOOLS_WRAPPER_PY_FILE = 'variant_tools/ucsctools.py'
UCSCTOOLS_INTERFACE_FILE = 'variant_tools/ucsctools.i'
SQLITE_PY_FILE = 'variant_tools/vt_sqlite3.py'


ASSOC_HEADER = [
    'variant_tools/assoTests.i',
    'variant_tools/assoTests.h',
    'variant_tools/assoData.h',
    'variant_tools/action.h',
    'variant_tools/utils.h',
    'variant_tools/lm.h'
]

ASSOC_FILES = [
    'variant_tools/assoData.cpp',
    'variant_tools/action.cpp',
    'variant_tools/utils.cpp',
    'variant_tools/lm.cpp'
]

SQLITE_FILES = [os.path.join(SQLITE_FOLDER, x) for x in [
        'cache.c',
        'connection.c',
        'cursor.c',
        'microprotocols.c',
        'module.c',
        'prepare_protocol.c',
        'row.c',
        'statement.c',
        'util.c'
    ]
] + [ 'sqlite/sqlite3.c' ]

SQLITE_GSL = [
    'gsl/error.c',
    'gsl/sys/infnan.c',
    'gsl/sys/coerce.c',
    'gsl/sys/fdiv.c',
    'gsl/sys/pow_int.c',
    'gsl/sys/fcmp.c',
    'gsl/sys/log1p.c',
    'gsl/sys/invhyp.c',
    'gsl/complex/math.c',
    'gsl/specfunc/log.c',
    'gsl/specfunc/exp.c',
    'gsl/specfunc/expint.c',
    'gsl/specfunc/gamma.c',
    'gsl/specfunc/zeta.c',
    'gsl/specfunc/trig.c',
    'gsl/specfunc/elementary.c',
    'gsl/specfunc/psi.c'
    ]

LIB_GSL = [
   'gsl/error.c',
   'gsl/sys/infnan.c',
   'gsl/sys/coerce.c',
   'gsl/sys/fdiv.c',
   'gsl/sys/pow_int.c',
   'gsl/sys/fcmp.c',
   'gsl/sys/log1p.c',
   'gsl/sys/invhyp.c',
   'gsl/complex/math.c',
   'gsl/specfunc/beta.c',
   'gsl/specfunc/psi.c',
   'gsl/specfunc/trig.c',
   'gsl/specfunc/exp.c',
   'gsl/specfunc/expint.c',
   'gsl/specfunc/log.c',
   'gsl/specfunc/erfc.c',
   'gsl/specfunc/zeta.c',
   'gsl/specfunc/elementary.c',
   'gsl/specfunc/gamma.c',
   'gsl/specfunc/gamma_inc.c',
   'gsl/rng/rng.c',
   'gsl/rng/default.c',
   'gsl/rng/mt.c',
   'gsl/rng/types.c',
   'gsl/randist/binomial.c',
   'gsl/randist/binomial_tpe.c',
   'gsl/randist/beta.c',
   'gsl/randist/exponential.c',
   'gsl/randist/geometric.c',
   'gsl/randist/nbinomial.c',
   'gsl/randist/poisson.c',
   'gsl/randist/multinomial.c',
   'gsl/randist/chisq.c',
   'gsl/randist/gauss.c',
   'gsl/randist/tdist.c',
   'gsl/randist/gausszig.c',
   'gsl/randist/gamma.c',
   'gsl/randist/hyperg.c',
   'gsl/cdf/binomial.c',
   'gsl/cdf/beta.c',
   'gsl/cdf/betainv.c',
   'gsl/cdf/gauss.c',
   'gsl/cdf/gaussinv.c',
   'gsl/cdf/tdist.c',
   'gsl/cdf/tdistinv.c',
   'gsl/cdf/chisq.c',
   'gsl/cdf/chisqinv.c',
   'gsl/cdf/gamma.c',
   'gsl/cdf/gammainv.c',
   'gsl/cdf/hypergeometric.c',
   'gsl/cdf/poisson.c',
   #
   'gsl/blas/blas.c',
   'gsl/cblas/caxpy.c',
   'gsl/cblas/ccopy.c',
   'gsl/cblas/cdotc_sub.c',
   'gsl/cblas/cdotu_sub.c',
   'gsl/cblas/cgbmv.c',
   'gsl/cblas/cgemm.c',
   'gsl/cblas/cgemv.c',
   'gsl/cblas/cgerc.c',
   'gsl/cblas/cgeru.c',
   'gsl/cblas/chbmv.c',
   'gsl/cblas/chemm.c',
   'gsl/cblas/chemv.c',
   'gsl/cblas/cher2.c',
   'gsl/cblas/cher2k.c',
   'gsl/cblas/cher.c',
   'gsl/cblas/cherk.c',
   'gsl/cblas/chpmv.c',
   'gsl/cblas/chpr2.c',
   'gsl/cblas/chpr.c',
   'gsl/cblas/cscal.c',
   'gsl/cblas/csscal.c',
   'gsl/cblas/cswap.c',
   'gsl/cblas/csymm.c',
   'gsl/cblas/csyr2k.c',
   'gsl/cblas/csyrk.c',
   'gsl/cblas/ctbmv.c',
   'gsl/cblas/ctbsv.c',
   'gsl/cblas/ctpmv.c',
   'gsl/cblas/ctpsv.c',
   'gsl/cblas/ctrmm.c',
   'gsl/cblas/ctrmv.c',
   'gsl/cblas/ctrsm.c',
   'gsl/cblas/ctrsv.c',
   'gsl/cblas/dasum.c',
   'gsl/cblas/daxpy.c',
   'gsl/cblas/dcopy.c',
   'gsl/cblas/ddot.c',
   'gsl/cblas/dgbmv.c',
   'gsl/cblas/dgemm.c',
   'gsl/cblas/dgemv.c',
   'gsl/cblas/dger.c',
   'gsl/cblas/dnrm2.c',
   'gsl/cblas/drot.c',
   'gsl/cblas/drotg.c',
   'gsl/cblas/drotm.c',
   'gsl/cblas/drotmg.c',
   'gsl/cblas/dsbmv.c',
   'gsl/cblas/dscal.c',
   'gsl/cblas/dsdot.c',
   'gsl/cblas/dspmv.c',
   'gsl/cblas/dspr2.c',
   'gsl/cblas/dspr.c',
   'gsl/cblas/dswap.c',
   'gsl/cblas/dsymm.c',
   'gsl/cblas/dsymv.c',
   'gsl/cblas/dsyr2.c',
   'gsl/cblas/dsyr2k.c',
   'gsl/cblas/dsyr.c',
   'gsl/cblas/dsyrk.c',
   'gsl/cblas/dtbmv.c',
   'gsl/cblas/dtbsv.c',
   'gsl/cblas/dtpmv.c',
   'gsl/cblas/dtpsv.c',
   'gsl/cblas/dtrmm.c',
   'gsl/cblas/dtrmv.c',
   'gsl/cblas/dtrsm.c',
   'gsl/cblas/dtrsv.c',
   'gsl/cblas/dzasum.c',
   'gsl/cblas/dznrm2.c',
   'gsl/cblas/hypot.c',
   'gsl/cblas/icamax.c',
   'gsl/cblas/idamax.c',
   'gsl/cblas/isamax.c',
   'gsl/cblas/izamax.c',
   'gsl/cblas/sasum.c',
   'gsl/cblas/saxpy.c',
   'gsl/cblas/scasum.c',
   'gsl/cblas/scnrm2.c',
   'gsl/cblas/scopy.c',
   'gsl/cblas/sdot.c',
   'gsl/cblas/sdsdot.c',
   'gsl/cblas/sgbmv.c',
   'gsl/cblas/sgemm.c',
   'gsl/cblas/sgemv.c',
   'gsl/cblas/sger.c',
   'gsl/cblas/snrm2.c',
   'gsl/cblas/srot.c',
   'gsl/cblas/srotg.c',
   'gsl/cblas/srotm.c',
   'gsl/cblas/srotmg.c',
   'gsl/cblas/ssbmv.c',
   'gsl/cblas/sscal.c',
   'gsl/cblas/sspmv.c',
   'gsl/cblas/sspr2.c',
   'gsl/cblas/sspr.c',
   'gsl/cblas/sswap.c',
   'gsl/cblas/ssymm.c',
   'gsl/cblas/ssymv.c',
   'gsl/cblas/ssyr2.c',
   'gsl/cblas/ssyr2k.c',
   'gsl/cblas/ssyr.c',
   'gsl/cblas/ssyrk.c',
   'gsl/cblas/stbmv.c',
   'gsl/cblas/stbsv.c',
   'gsl/cblas/stpmv.c',
   'gsl/cblas/stpsv.c',
   'gsl/cblas/strmm.c',
   'gsl/cblas/strmv.c',
   'gsl/cblas/strsm.c',
   'gsl/cblas/strsv.c',
   'gsl/cblas/xerbla.c',
   'gsl/cblas/zaxpy.c',
   'gsl/cblas/zcopy.c',
   'gsl/cblas/zdotc_sub.c',
   'gsl/cblas/zdotu_sub.c',
   'gsl/cblas/zdscal.c',
   'gsl/cblas/zgbmv.c',
   'gsl/cblas/zgemm.c',
   'gsl/cblas/zgemv.c',
   'gsl/cblas/zgerc.c',
   'gsl/cblas/zgeru.c',
   'gsl/cblas/zhbmv.c',
   'gsl/cblas/zhemm.c',
   'gsl/cblas/zhemv.c',
   'gsl/cblas/zher2.c',
   'gsl/cblas/zher2k.c',
   'gsl/cblas/zher.c',
   'gsl/cblas/zherk.c',
   'gsl/cblas/zhpmv.c',
   'gsl/cblas/zhpr2.c',
   'gsl/cblas/zhpr.c',
   'gsl/cblas/zscal.c',
   'gsl/cblas/zswap.c',
   'gsl/cblas/zsymm.c',
   'gsl/cblas/zsyr2k.c',
   'gsl/cblas/zsyrk.c',
   'gsl/cblas/ztbmv.c',
   'gsl/cblas/ztbsv.c',
   'gsl/cblas/ztpmv.c',
   'gsl/cblas/ztpsv.c',
   'gsl/cblas/ztrmm.c',
   'gsl/cblas/ztrmv.c',
   'gsl/cblas/ztrsm.c',
   'gsl/cblas/ztrsv.c',
   #
   'gsl/linalg/svd.c',
   'gsl/linalg/lu.c',
   'gsl/linalg/bidiag.c',
   'gsl/linalg/householder.c',
   'gsl/matrix/matrix.c',
   'gsl/matrix/submatrix.c',
   'gsl/matrix/rowcol.c',
   'gsl/matrix/getset.c',
   'gsl/matrix/init.c',
   'gsl/vector/init.c',
   'gsl/matrix/swap.c',
   'gsl/vector/vector.c',
   'gsl/vector/copy.c',
   'gsl/vector/swap.c',
   'gsl/vector/subvector.c',
   'gsl/vector/oper.c',
   'gsl/matrix/oper.c',
   'gsl/matrix/copy.c',
   'gsl/block/init.c',
   #
   'gsl/permutation/init.c',
   'gsl/permutation/inline.c',
   'gsl/permutation/permutation.c',
   'gsl/permutation/permute.c'
    ]

LIB_BOOST = [
    'boost_1_49_0/libs/iostreams/src/bzip2.cpp',
    'boost_1_49_0/libs/iostreams/src/file_descriptor.cpp',
    'boost_1_49_0/libs/iostreams/src/gzip.cpp',
    'boost_1_49_0/libs/iostreams/src/mapped_file.cpp',
    'boost_1_49_0/libs/iostreams/src/zlib.cpp',
    'boost_1_49_0/libs/regex/src/c_regex_traits.cpp',
    'boost_1_49_0/libs/regex/src/cpp_regex_traits.cpp',
    'boost_1_49_0/libs/regex/src/cregex.cpp',
    'boost_1_49_0/libs/regex/src/fileiter.cpp',
    'boost_1_49_0/libs/regex/src/icu.cpp',
    'boost_1_49_0/libs/regex/src/instances.cpp',
    'boost_1_49_0/libs/regex/src/posix_api.cpp',
    'boost_1_49_0/libs/regex/src/regex.cpp',
    'boost_1_49_0/libs/regex/src/regex_debug.cpp',
    'boost_1_49_0/libs/regex/src/regex_raw_buffer.cpp',
    'boost_1_49_0/libs/regex/src/regex_traits_defaults.cpp',
    'boost_1_49_0/libs/regex/src/static_mutex.cpp',
    'boost_1_49_0/libs/regex/src/usinstances.cpp',
    'boost_1_49_0/libs/regex/src/w32_regex_traits.cpp',
    'boost_1_49_0/libs/regex/src/wc_regex_traits.cpp',
    'boost_1_49_0/libs/regex/src/wide_posix_api.cpp',
    'boost_1_49_0/libs/regex/src/winstances.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/codecvt_error_category.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/operations.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/path.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/path_traits.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/portability.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/unique_path.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/utf8_codecvt_facet.cpp',
    'boost_1_49_0/libs/filesystem/v3/src/windows_file_codecvt.cpp',
    'boost_1_49_0/libs/filesystem/v2/src/v2_operations.cpp',
    'boost_1_49_0/libs/filesystem/v2/src/v2_portability.cpp',
    'boost_1_49_0/libs/filesystem/v2/src/v2_path.cpp',
    'boost_1_49_0/libs/system/src/error_code.cpp'
]

LIB_PLINKIO = [
    'libplinkio/cplinkio.c',
    'libplinkio/snparray.c',
    'libplinkio/fam.c',
    'libplinkio/fam_parse.c',
    'libplinkio/bim.c',
    'libplinkio/bim_parse.c',
    'libplinkio/bed.c',
    'libplinkio/plinkio.c',
    'libplinkio/file.c',
    'libplinkio/bed_header.c',
    'libplinkio/libcsv.c', 
]

LIB_CGATOOLS = [
   'cgatools/util/BaseUtil.cpp',
   'cgatools/util/Md5.cpp',
   'cgatools/util/DelimitedFile.cpp',
   'cgatools/util/RangeSet.cpp',
   'cgatools/util/DelimitedLineParser.cpp',
   'cgatools/util/Streams.cpp',
   'cgatools/util/Exception.cpp',
   'cgatools/util/StringSet.cpp',
   'cgatools/util/Files.cpp',
   'cgatools/util/parse.cpp',
   #'cgatools/util/GenericHistogram.cpp',
   'cgatools/reference/ChromosomeIdField.cpp',
   'cgatools/reference/CompactDnaSequence.cpp',
   'cgatools/reference/CrrFile.cpp',
   'cgatools/reference/CrrFileWriter.cpp',
   'cgatools/reference/GeneDataStore.cpp'
]

# http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
LIB_UCSC_FILES = [
    'ucsc/lib/osunix.c',
    'ucsc/lib/asParse.c',
    'ucsc/lib/errabort.c',
    'ucsc/lib/regexHelper.c',
    'ucsc/lib/tokenizer.c',
    'ucsc/lib/dnautil.c',
    'ucsc/lib/sqlNum.c',
    'ucsc/lib/hash.c',
    'ucsc/lib/memalloc.c',
    'ucsc/lib/dlist.c',
    'ucsc/lib/linefile.c',
    'ucsc/lib/cheapcgi.c',
    'ucsc/lib/portimpl.c',
    'ucsc/lib/servBrcMcw.c',
    'ucsc/lib/servcl.c',
    'ucsc/lib/servcis.c',
    'ucsc/lib/servCrunx.c',
    'ucsc/lib/servmsII.c',
    'ucsc/lib/servpws.c',
    'ucsc/lib/mime.c',
    'ucsc/lib/intExp.c',
    'ucsc/lib/kxTok.c',
    'ucsc/lib/localmem.c',
    'ucsc/lib/pipeline.c',
    'ucsc/lib/rangeTree.c',
    'ucsc/lib/rbTree.c',
    'ucsc/lib/basicBed.c',
    'ucsc/lib/binRange.c',
    'ucsc/lib/psl.c',
    'ucsc/lib/ffAli.c',
    'ucsc/lib/aliType.c',
    'ucsc/lib/sqlList.c',
    'ucsc/lib/udc.c',
    'ucsc/lib/bits.c',
    'ucsc/lib/base64.c',
    'ucsc/lib/net.c',
    'ucsc/lib/internet.c',
    'ucsc/lib/https.c',
    'ucsc/lib/verbose.c',
    'ucsc/lib/wildcmp.c',
    'ucsc/lib/zlibFace.c',
    'ucsc/lib/obscure.c',
    'ucsc/lib/dystring.c',
    'ucsc/lib/common.c',
    'ucsc/lib/bPlusTree.c',
    'ucsc/lib/hmmstats.c',
    'ucsc/lib/cirTree.c',
    'ucsc/lib/bbiRead.c',
    'ucsc/lib/bigBed.c',
    'ucsc/lib/bwgQuery.c',
    'ucsc/lib/vcf.c',
    'ucsc/lib/bamFile.c',
    'ucsc/lib/htmshell.c',
    'ucsc/lib/filePath.c',
    'ucsc/tabix/bedidx.c',
    'ucsc/tabix/index.c',
    #'ucsc/tabix/bgzf.c',
    'ucsc/tabix/knetfile.c',
    'ucsc/tabix/kstring.c',
    'ucsc/samtools/sam.c',
    'ucsc/samtools/bam.c',
    'ucsc/samtools/bgzf.c',
    'ucsc/samtools/razf.c',
    'ucsc/samtools/faidx.c',
    'ucsc/samtools/bam_pileup.c',
    'ucsc/samtools/bam_index.c',
    'ucsc/samtools/bam_import.c',
    'ucsc/samtools/sam_header.c',
    'ucsc/samtools/bam_aux.c',
]
    
LIB_STAT = ['variant_tools/fisher2.c']

if not EMBEDDED_BOOST:
    def downloadProgress(count, blockSize, totalSize):
        perc = count * blockSize * 100 // totalSize
        if perc > downloadProgress.counter: 
            sys.stdout.write('.' * (perc - downloadProgress.counter))
            downloadProgress.counter = perc
        sys.stdout.flush()
    if sys.version_info.major == 2:
        from urllib import urlretrieve
    else:
        from urllib.request import urlretrieve
    import tarfile
    downloadProgress.counter = 0
    try:
        BOOST_URL = 'http://downloads.sourceforge.net/project/boost/boost/1.49.0/boost_1_49_0.tar.gz?r=&ts=1435893980&use_mirror=iweb'
        sys.stdout.write('Downloading boost C++ library 1.49.0 ')
        sys.stdout.flush()
        if not os.path.isfile('boost_1_49_0.tar.gz'):
            urlretrieve(BOOST_URL, 'boost_1_49_0.tar.gz', downloadProgress)
        sys.stdout.write('\n')
        # extract needed files
        with tarfile.open('boost_1_49_0.tar.gz', 'r:gz') as tar:
            files = [h for h in tar.getmembers() if h.name.startswith('boost_1_49_0/boost') \
                or h.name.startswith('boost_1_49_0/libs/iostreams') \
                or h.name.startswith('boost_1_49_0/libs/regex') \
                or h.name.startswith('boost_1_49_0/libs/filesystem') \
                or h.name.startswith('boost_1_49_0/libs/detail') \
                or h.name.startswith('boost_1_49_0/libs/system') ]
            sys.stdout.write('Extracting %d files\n' % len(files))
            tar.extractall(members=files)
        EMBEDDED_BOOST = True
    except Exception as e:
        print(e)
        print('The boost C++ library version 1.49.0 is not found under the current directory. Will try to use the system libraries.')


#
# Generate wrapper files (only in development mode)
#
#
if not os.path.isfile('variant_tools/swigpyrun.h'):
    try:
       ret = subprocess.call(['swig -python -external-runtime variant_tools/swigpyrun.h'], shell=True)
       if ret != 0: sys.exit('Failed to generate swig runtime header file.')
    except OSError as e:
        sys.exit('Failed to generate wrapper file. Please install swig (www.swig.org).')
#
# generate wrapper files for both versions of python. This will make sure sdist gets
# all files needed for the source package
#
# we re-generate wrapper files for all versions of python only for
# source distribution
if (not os.path.isfile(ASSO_WRAPPER_PY_FILE) or not os.path.isfile(ASSO_WRAPPER_CPP_FILE) \
  or os.path.getmtime(ASSO_WRAPPER_CPP_FILE) < max([os.path.getmtime(x) for x in ASSOC_HEADER + ASSOC_FILES])):
    print('Generating {}'.format(ASSO_WRAPPER_CPP_FILE))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', ASSO_WRAPPER_CPP_FILE, ASSO_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for association module.')
    os.rename('variant_tools/assoTests.py', ASSO_WRAPPER_PY_FILE)
#
if (not os.path.isfile(CGATOOLS_WRAPPER_PY_FILE)) or (not os.path.isfile(CGATOOLS_WRAPPER_CPP_FILE)) \
  or os.path.getmtime(CGATOOLS_WRAPPER_CPP_FILE) < os.path.getmtime(CGATOOLS_INTERFACE_FILE):
    print('Generating {}'.format(CGATOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', CGATOOLS_WRAPPER_CPP_FILE, CGATOOLS_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for cgatools.')
    os.rename('variant_tools/cgatools.py', CGATOOLS_WRAPPER_PY_FILE)
#
if (not os.path.isfile(UCSCTOOLS_WRAPPER_PY_FILE)) or (not os.path.isfile(UCSCTOOLS_WRAPPER_CPP_FILE)) \
  or os.path.getmtime(UCSCTOOLS_WRAPPER_CPP_FILE) < os.path.getmtime(UCSCTOOLS_INTERFACE_FILE):
    print('Generating {}'.format(UCSCTOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', UCSCTOOLS_WRAPPER_CPP_FILE, UCSCTOOLS_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for ucsctools.')
    os.rename('variant_tools/ucsctools.py', UCSCTOOLS_WRAPPER_PY_FILE)
         
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
if sys.platform in 'linux2':
    libs = ['stdc++']
    # gcc optimizations
    # gccargs = ["-O3", "-march=native"]
    gccargs = ["-O3", '-Wno-unused-local-typedef']
else:
    # mac system
    libs = []
    gccargs = ['-O3', '-Wno-unused-local-typedef', '-Wno-return-type']

ENV_INCLUDE_DIRS = [x for x in os.environ.get('LD_INCLUDE_PATH', '').split(os.pathsep) if x]
ENV_LIBRARY_DIRS = [x for x in os.environ.get('LD_LIBRARY_PATH', '').split(os.pathsep) if x]

if EMBEDDED_BOOST:
    try:
        c = new_compiler()
        if not os.path.isfile(os.path.join('build', c.static_lib_format % ('embedded_boost', c.static_lib_extension))):
            # -w suppress all warnings caused by the use of boost libraries
            objects = c.compile(LIB_BOOST,
                include_dirs=['boost_1_49_0'] + ENV_INCLUDE_DIRS,
                output_dir='build',
                extra_preargs = ['-w', '-fPIC'],
                macros = [('BOOST_ALL_NO_LIB', None)])
            c.create_static_lib(objects, "embedded_boost", output_dir='build')
    except Exception as e:
        sys.exit("Failed to build embedded boost library: {}".format(e))

# building other libraries
for files, incs, macs, libname in [
    (SQLITE_GSL, ['.', 'gsl'], [], 'sqlite_gsl'),
    (LIB_STAT, ['.'], [], 'stat'),
    (LIB_UCSC_FILES, ['ucsc/inc', 'ucsc/tabix', 'ucsc/samtools'], 
        [('USE_TABIX', '1'), ('_FILE_OFFSET_BITS', '64'), ('USE_BAM', '1'),
         ('_USE_KNETFILE', None), ('BGZF_CACHE', None)],
        'ucsc'),
    (LIB_CGATOOLS, ['.', 'cgatools', 'boost_1_49_0'],
        [('CGA_TOOLS_IS_PIPELINE', 0), ('CGA_TOOLS_VERSION', r'"1.6.0.43"')],
        'cgatools'),
    (LIB_GSL, ['.', 'gsl'], [], 'gsl')]:
    if os.path.isfile(os.path.join('build', c.static_lib_format % (libname, c.static_lib_extension))):
        continue
    try:
        print('Building embedded library {}'.format(libname))
        c = new_compiler()
        # -w suppress all warnings caused by the use of boost libraries
        objects = c.compile(files,
            include_dirs=incs + ENV_INCLUDE_DIRS,
            output_dir='build',
            extra_preargs = ['-w', '-fPIC'],
            macros = macs)
        c.create_static_lib(objects, libname, output_dir='build')
    except Exception as e:
        sys.exit("Failed to build embedded {} library: {}".format(libname, e))




setup(name = "variant_tools",
    version = VTOOLS_VERSION,
    description = "Variant tools: an integrated annotation and analysis package for next-generation sequencing data",
    author = 'Bo Peng',
    url = 'http://varianttools.sourceforge.net',
    author_email = 'bpeng@mdanderson.org',
    maintainer = 'Bo Peng',
    maintainer_email = 'varianttools-devel@lists.sourceforge.net',
    license = 'GPL3',
    classifiers = [
            'Development Status :: 5 - Production/Stable',
            'Environment :: Console',
            'Intended Audience :: Education',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Natural Language :: English',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS :: MacOS X',
            'Programming Language :: C++',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
        ],
    #packages = ['variant_tools'],
    scripts = [
        'vtools',
        'vtools_report'
    ],
    cmdclass = {'build_py': build_py },
    #package_dir = {'variant_tools': 'source'},
    py_modules = [
        # __init__ is automatically installed
        'variant_tools.annotation',
        'variant_tools.association',
        'variant_tools.assoTests',
        'variant_tools.cgatools',
        'variant_tools.compare',
        'variant_tools.exporter',
        'variant_tools.importer',
        'variant_tools.liftOver',
        'variant_tools.meta',
        'variant_tools.phenotype',
        'variant_tools.pipeline',
        'variant_tools.plinkfile',
        'variant_tools.plot',
        'variant_tools.project',
        'variant_tools.preprocessor',
        'variant_tools.rtester',
        'variant_tools.simulation',
        'variant_tools.site_options',
        'variant_tools.tester',
        'variant_tools.ucsctools',
        'variant_tools.update',
        'variant_tools.utils',
        'variant_tools.variant',
        'variant_tools.vt_sqlite3',
        ],
    ext_modules = [
        Extension('variant_tools._vt_sqlite3',
            # stop warning message for sqlite because it is written by us.
            extra_compile_args=['-w'],
            sources = SQLITE_FILES,
            define_macros = [('MODULE_NAME', '"vt_sqlite3"'), ('HAVE_USLEEP', None)],
            include_dirs = ['sqlite', SQLITE_FOLDER] + ENV_INCLUDE_DIRS,
        ),
        Extension('variant_tools._ucsctools',
            # stop warning message ucsctools because it is written by us.
            extra_compile_args=['-w'],
            sources = [UCSCTOOLS_WRAPPER_CPP_FILE],
            include_dirs = ['.', 'ucsc/inc', 'ucsc/tabix', 'ucsc/samtools'] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,
            define_macros =  [('USE_TABIX', '1'), ('_FILE_OFFSET_BITS', '64'), ('USE_BAM', '1'),
                ('_USE_KNETFILE', None), ('BGZF_CACHE', None)],
            libraries = libs + ['ucsc', 'z', 'bz2'],
        ),
        Extension('variant_tools.cplinkio',
            # stop warning message ucsctools because it is written by us.
            extra_compile_args=['-w'],
            sources = LIB_PLINKIO,
            include_dirs = ['libplinkio'] + ENV_INCLUDE_DIRS,
        ),
        Extension('variant_tools._vt_sqlite3_ext',
            # stop warning message for sqlite because it is written by us.
            sources = ['sqlite/vt_sqlite3_ext.cpp'],
            include_dirs = [".", 'ucsc/inc', 'ucsc/tabix', 'ucsc/samtools',
                'sqlite', "variant_tools", "gsl", "cgatools", "boost_1_49_0"] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"]  + ENV_LIBRARY_DIRS,
            libraries = ['sqlite_gsl', 'stat', 'ucsc', 'cgatools'] + \
                (['embedded_boost'] if EMBEDDED_BOOST else ['boost_iostreams', 'boost_regex', 'boost_filesystem']) + \
                ['z', 'bz2'],
            extra_compile_args = gccargs,
            define_macros = [
                ('MODULE_NAME', '"vt_sqlite3"'),
                ('BOOST_ALL_NO_LIB', None),  ('CGA_TOOLS_IS_PIPELINE', 0),
                ('CGA_TOOLS_VERSION', r'"1.6.0.43"'), ('USE_TABIX', '1'), ('USE_BAM', '1'),
                ('_FILE_OFFSET_BITS', '64'), ('_USE_KNETFILE', None),
                ('BGZF_CACHE', None)],
        ),
        Extension('variant_tools._cgatools',
            sources = [CGATOOLS_WRAPPER_CPP_FILE],
            libraries = ['cgatools'] + \
                (['embedded_boost'] if EMBEDDED_BOOST else ['boost_iostreams', 'boost_regex', 'boost_filesystem']) + \
                ['z', 'bz2'] + ENV_INCLUDE_DIRS,
            define_macros = [('BOOST_ALL_NO_LIB', None),  ('CGA_TOOLS_IS_PIPELINE', 0),
                ('CGA_TOOLS_VERSION', r'"1.6.0.43"')],
            extra_compile_args = gccargs,
            swig_opts = ['-O', '-shadow', '-c++', '-keyword'],
            include_dirs = [".", "cgatools", "boost_1_49_0"]  + ENV_LIBRARY_DIRS,
            library_dirs = ["build"],
        ),
        Extension('variant_tools._assoTests',
            sources = [ASSO_WRAPPER_CPP_FILE] + ASSOC_FILES,
            extra_compile_args = gccargs,
            libraries = libs + ['gsl', 'stat', 'blas'],
            library_dirs = ["build"]  + ENV_LIBRARY_DIRS,
            include_dirs = [".", "variant_tools", "gsl"] + ENV_INCLUDE_DIRS,
        )
      ] 
)

