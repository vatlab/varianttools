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
from distutils.ccompiler import new_compiler
from setuptools import find_packages, setup, Extension
import numpy as np
import sys, os, subprocess
from Cython.Build import cythonize
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
    sys.exit('variant tools requires Python 3.2 or higher. Please upgrade your version (%s) of Python and try again.' % (sys.version.split()[0]))

# do not import variant_tools because __init__ might not be imported properly 
# before installation
with open('src/variant_tools/_version.py') as init:
    for line in init:
        if line.startswith('VTOOLS_VERSION='):
            VTOOLS_VERSION = line[15:].strip().strip('"').strip("'")
            break

SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword', '-w-511', '-w-509',
    '-outdir', 'src/variant_tools', '-py3']

ASSO_WRAPPER_CPP_FILE = 'src/variant_tools/assoTests_wrap.cpp'
ASSO_WRAPPER_PY_FILE = 'src/variant_tools/assoTests.py'
ASSO_INTERFACE_FILE = 'src/variant_tools/assoTests.i'
CGATOOLS_WRAPPER_CPP_FILE = 'src/variant_tools/cgatools_wrap.cpp'
CGATOOLS_WRAPPER_PY_FILE = 'src/variant_tools/cgatools.py'
CGATOOLS_INTERFACE_FILE = 'src/variant_tools/cgatools.i'
UCSCTOOLS_WRAPPER_CPP_FILE = 'src/variant_tools/ucsctools_wrap.cpp'
UCSCTOOLS_WRAPPER_PY_FILE = 'src/variant_tools/ucsctools.py'
UCSCTOOLS_INTERFACE_FILE = 'src/variant_tools/ucsctools.i'
SQLITE_PY_FILE = 'src/variant_tools/vt_sqlite3.py'


ASSOC_HEADER = [
    'src/variant_tools/assoTests.i',
    'src/variant_tools/assoTests.h',
    'src/variant_tools/assoData.h',
    'src/variant_tools/action.h',
    'src/variant_tools/utils.h',
    'src/variant_tools/lm.h'
]

ASSOC_FILES = [
    'src/variant_tools/assoData.cpp',
    'src/variant_tools/action.cpp',
    'src/variant_tools/utils.cpp',
    'src/variant_tools/lm.cpp'
]

SQLITE_FOLDER = 'src/sqlite/py3'
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
] + [ 'src/sqlite/sqlite3.c' ]

SQLITE_GSL = [
    'src/gsl/error.c',
    'src/gsl/sys/infnan.c',
    'src/gsl/sys/coerce.c',
    'src/gsl/sys/fdiv.c',
    'src/gsl/sys/pow_int.c',
    'src/gsl/sys/fcmp.c',
    'src/gsl/sys/log1p.c',
    'src/gsl/sys/invhyp.c',
    'src/gsl/complex/math.c',
    'src/gsl/specfunc/log.c',
    'src/gsl/specfunc/exp.c',
    'src/gsl/specfunc/expint.c',
    'src/gsl/specfunc/gamma.c',
    'src/gsl/specfunc/zeta.c',
    'src/gsl/specfunc/trig.c',
    'src/gsl/specfunc/elementary.c',
    'src/gsl/specfunc/psi.c'
    ]

LIB_GSL = [
   'src/gsl/error.c',
   'src/gsl/sys/infnan.c',
   'src/gsl/sys/coerce.c',
   'src/gsl/sys/fdiv.c',
   'src/gsl/sys/pow_int.c',
   'src/gsl/sys/fcmp.c',
   'src/gsl/sys/log1p.c',
   'src/gsl/sys/invhyp.c',
   'src/gsl/complex/math.c',
   'src/gsl/specfunc/beta.c',
   'src/gsl/specfunc/psi.c',
   'src/gsl/specfunc/trig.c',
   'src/gsl/specfunc/exp.c',
   'src/gsl/specfunc/expint.c',
   'src/gsl/specfunc/log.c',
   'src/gsl/specfunc/erfc.c',
   'src/gsl/specfunc/zeta.c',
   'src/gsl/specfunc/elementary.c',
   'src/gsl/specfunc/gamma.c',
   'src/gsl/specfunc/gamma_inc.c',
   'src/gsl/rng/rng.c',
   'src/gsl/rng/default.c',
   'src/gsl/rng/mt.c',
   'src/gsl/rng/types.c',
   'src/gsl/randist/binomial.c',
   'src/gsl/randist/binomial_tpe.c',
   'src/gsl/randist/beta.c',
   'src/gsl/randist/exponential.c',
   'src/gsl/randist/geometric.c',
   'src/gsl/randist/nbinomial.c',
   'src/gsl/randist/poisson.c',
   'src/gsl/randist/multinomial.c',
   'src/gsl/randist/chisq.c',
   'src/gsl/randist/gauss.c',
   'src/gsl/randist/tdist.c',
   'src/gsl/randist/gausszig.c',
   'src/gsl/randist/gamma.c',
   'src/gsl/randist/hyperg.c',
   'src/gsl/cdf/binomial.c',
   'src/gsl/cdf/beta.c',
   'src/gsl/cdf/betainv.c',
   'src/gsl/cdf/gauss.c',
   'src/gsl/cdf/gaussinv.c',
   'src/gsl/cdf/tdist.c',
   'src/gsl/cdf/tdistinv.c',
   'src/gsl/cdf/chisq.c',
   'src/gsl/cdf/chisqinv.c',
   'src/gsl/cdf/gamma.c',
   'src/gsl/cdf/gammainv.c',
   'src/gsl/cdf/hypergeometric.c',
   'src/gsl/cdf/poisson.c',
   #
   'src/gsl/blas/blas.c',
   'src/gsl/cblas/caxpy.c',
   'src/gsl/cblas/ccopy.c',
   'src/gsl/cblas/cdotc_sub.c',
   'src/gsl/cblas/cdotu_sub.c',
   'src/gsl/cblas/cgbmv.c',
   'src/gsl/cblas/cgemm.c',
   'src/gsl/cblas/cgemv.c',
   'src/gsl/cblas/cgerc.c',
   'src/gsl/cblas/cgeru.c',
   'src/gsl/cblas/chbmv.c',
   'src/gsl/cblas/chemm.c',
   'src/gsl/cblas/chemv.c',
   'src/gsl/cblas/cher2.c',
   'src/gsl/cblas/cher2k.c',
   'src/gsl/cblas/cher.c',
   'src/gsl/cblas/cherk.c',
   'src/gsl/cblas/chpmv.c',
   'src/gsl/cblas/chpr2.c',
   'src/gsl/cblas/chpr.c',
   'src/gsl/cblas/cscal.c',
   'src/gsl/cblas/csscal.c',
   'src/gsl/cblas/cswap.c',
   'src/gsl/cblas/csymm.c',
   'src/gsl/cblas/csyr2k.c',
   'src/gsl/cblas/csyrk.c',
   'src/gsl/cblas/ctbmv.c',
   'src/gsl/cblas/ctbsv.c',
   'src/gsl/cblas/ctpmv.c',
   'src/gsl/cblas/ctpsv.c',
   'src/gsl/cblas/ctrmm.c',
   'src/gsl/cblas/ctrmv.c',
   'src/gsl/cblas/ctrsm.c',
   'src/gsl/cblas/ctrsv.c',
   'src/gsl/cblas/dasum.c',
   'src/gsl/cblas/daxpy.c',
   'src/gsl/cblas/dcopy.c',
   'src/gsl/cblas/ddot.c',
   'src/gsl/cblas/dgbmv.c',
   'src/gsl/cblas/dgemm.c',
   'src/gsl/cblas/dgemv.c',
   'src/gsl/cblas/dger.c',
   'src/gsl/cblas/dnrm2.c',
   'src/gsl/cblas/drot.c',
   'src/gsl/cblas/drotg.c',
   'src/gsl/cblas/drotm.c',
   'src/gsl/cblas/drotmg.c',
   'src/gsl/cblas/dsbmv.c',
   'src/gsl/cblas/dscal.c',
   'src/gsl/cblas/dsdot.c',
   'src/gsl/cblas/dspmv.c',
   'src/gsl/cblas/dspr2.c',
   'src/gsl/cblas/dspr.c',
   'src/gsl/cblas/dswap.c',
   'src/gsl/cblas/dsymm.c',
   'src/gsl/cblas/dsymv.c',
   'src/gsl/cblas/dsyr2.c',
   'src/gsl/cblas/dsyr2k.c',
   'src/gsl/cblas/dsyr.c',
   'src/gsl/cblas/dsyrk.c',
   'src/gsl/cblas/dtbmv.c',
   'src/gsl/cblas/dtbsv.c',
   'src/gsl/cblas/dtpmv.c',
   'src/gsl/cblas/dtpsv.c',
   'src/gsl/cblas/dtrmm.c',
   'src/gsl/cblas/dtrmv.c',
   'src/gsl/cblas/dtrsm.c',
   'src/gsl/cblas/dtrsv.c',
   'src/gsl/cblas/dzasum.c',
   'src/gsl/cblas/dznrm2.c',
   'src/gsl/cblas/hypot.c',
   'src/gsl/cblas/icamax.c',
   'src/gsl/cblas/idamax.c',
   'src/gsl/cblas/isamax.c',
   'src/gsl/cblas/izamax.c',
   'src/gsl/cblas/sasum.c',
   'src/gsl/cblas/saxpy.c',
   'src/gsl/cblas/scasum.c',
   'src/gsl/cblas/scnrm2.c',
   'src/gsl/cblas/scopy.c',
   'src/gsl/cblas/sdot.c',
   'src/gsl/cblas/sdsdot.c',
   'src/gsl/cblas/sgbmv.c',
   'src/gsl/cblas/sgemm.c',
   'src/gsl/cblas/sgemv.c',
   'src/gsl/cblas/sger.c',
   'src/gsl/cblas/snrm2.c',
   'src/gsl/cblas/srot.c',
   'src/gsl/cblas/srotg.c',
   'src/gsl/cblas/srotm.c',
   'src/gsl/cblas/srotmg.c',
   'src/gsl/cblas/ssbmv.c',
   'src/gsl/cblas/sscal.c',
   'src/gsl/cblas/sspmv.c',
   'src/gsl/cblas/sspr2.c',
   'src/gsl/cblas/sspr.c',
   'src/gsl/cblas/sswap.c',
   'src/gsl/cblas/ssymm.c',
   'src/gsl/cblas/ssymv.c',
   'src/gsl/cblas/ssyr2.c',
   'src/gsl/cblas/ssyr2k.c',
   'src/gsl/cblas/ssyr.c',
   'src/gsl/cblas/ssyrk.c',
   'src/gsl/cblas/stbmv.c',
   'src/gsl/cblas/stbsv.c',
   'src/gsl/cblas/stpmv.c',
   'src/gsl/cblas/stpsv.c',
   'src/gsl/cblas/strmm.c',
   'src/gsl/cblas/strmv.c',
   'src/gsl/cblas/strsm.c',
   'src/gsl/cblas/strsv.c',
   'src/gsl/cblas/xerbla.c',
   'src/gsl/cblas/zaxpy.c',
   'src/gsl/cblas/zcopy.c',
   'src/gsl/cblas/zdotc_sub.c',
   'src/gsl/cblas/zdotu_sub.c',
   'src/gsl/cblas/zdscal.c',
   'src/gsl/cblas/zgbmv.c',
   'src/gsl/cblas/zgemm.c',
   'src/gsl/cblas/zgemv.c',
   'src/gsl/cblas/zgerc.c',
   'src/gsl/cblas/zgeru.c',
   'src/gsl/cblas/zhbmv.c',
   'src/gsl/cblas/zhemm.c',
   'src/gsl/cblas/zhemv.c',
   'src/gsl/cblas/zher2.c',
   'src/gsl/cblas/zher2k.c',
   'src/gsl/cblas/zher.c',
   'src/gsl/cblas/zherk.c',
   'src/gsl/cblas/zhpmv.c',
   'src/gsl/cblas/zhpr2.c',
   'src/gsl/cblas/zhpr.c',
   'src/gsl/cblas/zscal.c',
   'src/gsl/cblas/zswap.c',
   'src/gsl/cblas/zsymm.c',
   'src/gsl/cblas/zsyr2k.c',
   'src/gsl/cblas/zsyrk.c',
   'src/gsl/cblas/ztbmv.c',
   'src/gsl/cblas/ztbsv.c',
   'src/gsl/cblas/ztpmv.c',
   'src/gsl/cblas/ztpsv.c',
   'src/gsl/cblas/ztrmm.c',
   'src/gsl/cblas/ztrmv.c',
   'src/gsl/cblas/ztrsm.c',
   'src/gsl/cblas/ztrsv.c',
   #
   'src/gsl/linalg/svd.c',
   'src/gsl/linalg/lu.c',
   'src/gsl/linalg/bidiag.c',
   'src/gsl/linalg/householder.c',
   'src/gsl/matrix/matrix.c',
   'src/gsl/matrix/submatrix.c',
   'src/gsl/matrix/rowcol.c',
   'src/gsl/matrix/getset.c',
   'src/gsl/matrix/init.c',
   'src/gsl/vector/init.c',
   'src/gsl/matrix/swap.c',
   'src/gsl/vector/vector.c',
   'src/gsl/vector/copy.c',
   'src/gsl/vector/swap.c',
   'src/gsl/vector/subvector.c',
   'src/gsl/vector/oper.c',
   'src/gsl/matrix/oper.c',
   'src/gsl/matrix/copy.c',
   'src/gsl/block/init.c',
   #
   'src/gsl/permutation/init.c',
   'src/gsl/permutation/inline.c',
   'src/gsl/permutation/permutation.c',
   'src/gsl/permutation/permute.c'
    ]

LIB_BOOST = [
    'src/boost_1_49_0/libs/iostreams/src/bzip2.cpp',
    'src/boost_1_49_0/libs/iostreams/src/file_descriptor.cpp',
    'src/boost_1_49_0/libs/iostreams/src/gzip.cpp',
    'src/boost_1_49_0/libs/iostreams/src/mapped_file.cpp',
    'src/boost_1_49_0/libs/iostreams/src/zlib.cpp',
    'src/boost_1_49_0/libs/regex/src/c_regex_traits.cpp',
    'src/boost_1_49_0/libs/regex/src/cpp_regex_traits.cpp',
    'src/boost_1_49_0/libs/regex/src/cregex.cpp',
    'src/boost_1_49_0/libs/regex/src/fileiter.cpp',
    'src/boost_1_49_0/libs/regex/src/icu.cpp',
    'src/boost_1_49_0/libs/regex/src/instances.cpp',
    'src/boost_1_49_0/libs/regex/src/posix_api.cpp',
    'src/boost_1_49_0/libs/regex/src/regex.cpp',
    'src/boost_1_49_0/libs/regex/src/regex_debug.cpp',
    'src/boost_1_49_0/libs/regex/src/regex_raw_buffer.cpp',
    'src/boost_1_49_0/libs/regex/src/regex_traits_defaults.cpp',
    'src/boost_1_49_0/libs/regex/src/static_mutex.cpp',
    'src/boost_1_49_0/libs/regex/src/usinstances.cpp',
    'src/boost_1_49_0/libs/regex/src/w32_regex_traits.cpp',
    'src/boost_1_49_0/libs/regex/src/wc_regex_traits.cpp',
    'src/boost_1_49_0/libs/regex/src/wide_posix_api.cpp',
    'src/boost_1_49_0/libs/regex/src/winstances.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/codecvt_error_category.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/operations.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/path.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/path_traits.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/portability.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/unique_path.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/utf8_codecvt_facet.cpp',
    'src/boost_1_49_0/libs/filesystem/v3/src/windows_file_codecvt.cpp',
    'src/boost_1_49_0/libs/filesystem/v2/src/v2_operations.cpp',
    'src/boost_1_49_0/libs/filesystem/v2/src/v2_portability.cpp',
    'src/boost_1_49_0/libs/filesystem/v2/src/v2_path.cpp',
    'src/boost_1_49_0/libs/system/src/error_code.cpp'
]

LIB_PLINKIO = [
    'src/libplinkio/cplinkio.c',
    'src/libplinkio/snparray.c',
    'src/libplinkio/fam.c',
    'src/libplinkio/fam_parse.c',
    'src/libplinkio/bim.c',
    'src/libplinkio/bim_parse.c',
    'src/libplinkio/bed.c',
    'src/libplinkio/plinkio.c',
    'src/libplinkio/file.c',
    'src/libplinkio/bed_header.c',
    'src/libplinkio/libcsv.c', 
]

LIB_CGATOOLS = [
   'src/cgatools/util/BaseUtil.cpp',
   'src/cgatools/util/Md5.cpp',
   'src/cgatools/util/DelimitedFile.cpp',
   'src/cgatools/util/RangeSet.cpp',
   'src/cgatools/util/DelimitedLineParser.cpp',
   'src/cgatools/util/Streams.cpp',
   'src/cgatools/util/Exception.cpp',
   'src/cgatools/util/StringSet.cpp',
   'src/cgatools/util/Files.cpp',
   'src/cgatools/util/parse.cpp',
   #'src/cgatools/util/GenericHistogram.cpp',
   'src/cgatools/reference/ChromosomeIdField.cpp',
   'src/cgatools/reference/CompactDnaSequence.cpp',
   'src/cgatools/reference/CrrFile.cpp',
   'src/cgatools/reference/CrrFileWriter.cpp',
   'src/cgatools/reference/GeneDataStore.cpp'
]

# http://hgdownload.cse.ucsc.edu/admin/jksrc.zip
LIB_UCSC_FILES = [
    'src/ucsc/lib/osunix.c',
    'src/ucsc/lib/asParse.c',
    'src/ucsc/lib/errabort.c',
    'src/ucsc/lib/regexHelper.c',
    'src/ucsc/lib/tokenizer.c',
    'src/ucsc/lib/dnautil.c',
    'src/ucsc/lib/sqlNum.c',
    'src/ucsc/lib/hash.c',
    'src/ucsc/lib/memalloc.c',
    'src/ucsc/lib/dlist.c',
    'src/ucsc/lib/linefile.c',
    'src/ucsc/lib/cheapcgi.c',
    'src/ucsc/lib/portimpl.c',
    'src/ucsc/lib/servBrcMcw.c',
    'src/ucsc/lib/servcl.c',
    'src/ucsc/lib/servcis.c',
    'src/ucsc/lib/servCrunx.c',
    'src/ucsc/lib/servmsII.c',
    'src/ucsc/lib/servpws.c',
    'src/ucsc/lib/mime.c',
    'src/ucsc/lib/intExp.c',
    'src/ucsc/lib/kxTok.c',
    'src/ucsc/lib/localmem.c',
    'src/ucsc/lib/pipeline.c',
    'src/ucsc/lib/rangeTree.c',
    'src/ucsc/lib/rbTree.c',
    'src/ucsc/lib/basicBed.c',
    'src/ucsc/lib/binRange.c',
    'src/ucsc/lib/psl.c',
    'src/ucsc/lib/ffAli.c',
    'src/ucsc/lib/aliType.c',
    'src/ucsc/lib/sqlList.c',
    'src/ucsc/lib/udc.c',
    'src/ucsc/lib/bits.c',
    'src/ucsc/lib/base64.c',
    'src/ucsc/lib/net.c',
    'src/ucsc/lib/internet.c',
    'src/ucsc/lib/https.c',
    'src/ucsc/lib/verbose.c',
    'src/ucsc/lib/wildcmp.c',
    'src/ucsc/lib/zlibFace.c',
    'src/ucsc/lib/obscure.c',
    'src/ucsc/lib/dystring.c',
    'src/ucsc/lib/common.c',
    'src/ucsc/lib/bPlusTree.c',
    'src/ucsc/lib/hmmstats.c',
    'src/ucsc/lib/cirTree.c',
    'src/ucsc/lib/bbiRead.c',
    'src/ucsc/lib/bigBed.c',
    'src/ucsc/lib/bwgQuery.c',
    'src/ucsc/lib/vcf.c',
    'src/ucsc/lib/bamFile.c',
    'src/ucsc/lib/htmshell.c',
    'src/ucsc/lib/filePath.c',
    'src/ucsc/tabix/bedidx.c',
    'src/ucsc/tabix/index.c',
    #'src/ucsc/tabix/bgzf.c',
    'src/ucsc/tabix/knetfile.c',
    'src/ucsc/tabix/kstring.c',
    'src/ucsc/samtools/sam.c',
    'src/ucsc/samtools/bam.c',
    'src/ucsc/samtools/bgzf.c',
    'src/ucsc/samtools/razf.c',
    'src/ucsc/samtools/faidx.c',
    'src/ucsc/samtools/bam_pileup.c',
    'src/ucsc/samtools/bam_index.c',
    'src/ucsc/samtools/bam_import.c',
    'src/ucsc/samtools/sam_header.c',
    'src/ucsc/samtools/bam_aux.c',
]
    
LIB_STAT = ['src/variant_tools/fisher2.c']

EMBEDDED_BOOST = os.path.isdir('src/boost_1_49_0')
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
        if not os.path.isfile('src/boost_1_49_0.tar.gz'):
            urlretrieve(BOOST_URL, 'src/boost_1_49_0.tar.gz', downloadProgress)
        sys.stdout.write('\n')
        # extract needed files
        with tarfile.open('src/boost_1_49_0.tar.gz', 'r:gz') as tar:
            files = [h for h in tar.getmembers() if h.name.startswith('boost_1_49_0/boost') \
                or h.name.startswith('boost_1_49_0/libs/iostreams') \
                or h.name.startswith('boost_1_49_0/libs/regex') \
                or h.name.startswith('boost_1_49_0/libs/filesystem') \
                or h.name.startswith('boost_1_49_0/libs/detail') \
                or h.name.startswith('boost_1_49_0/libs/system') ]
            sys.stdout.write('Extracting %d files\n' % len(files))
            tar.extractall(path='src', members=files)
        os.remove('src/boost_1_49_0.tar.gz')
        EMBEDDED_BOOST = True
    except Exception as e:
        print(e)
        print('The boost C++ library version 1.49.0 is not found under the current directory. Will try to use the system libraries.')


EMBEDDED_ZMQ= os.path.isdir('src/zeromq-4.0.3')
if not EMBEDDED_ZMQ:
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
        BOOST_URL = 'http://download.zeromq.org/zeromq-4.0.3.tar.gz'
        sys.stdout.write('Downloading ZeroMQ version 4.0.3')
        sys.stdout.flush()
        if not os.path.isfile('src/zeromq-4.0.3.tar.gz'):
            urlretrieve(BOOST_URL, 'src/zeromq-4.0.3.tar.gz', downloadProgress)
        sys.stdout.write('\n')
        # extract needed files
        with tarfile.open('src/zeromq-4.0.3.tar.gz', 'r:gz') as tar:
            tar.extractall(path='src')
        os.remove('src/zeromq-4.0.3.tar.gz')
        EMBEDDED_ZMQ = True
    except Exception as e:
        print(e)
        print('The ZeroMQ version 4.0.3 is not found under the current directory.')

#
# During development, if an interface file needs to be re-generated, please
# remove these files and they will be re-generated with SWIG
#
if not os.path.isfile('src/swigpyrun.h'):
    try:
       ret = subprocess.call(['swig -python -external-runtime src/swigpyrun.h -o build'], shell=True)
       if ret != 0: sys.exit('Failed to generate swig runtime header file.')
    except OSError as e:
        sys.exit('Failed to generate wrapper file. Please install swig (www.swig.org).')
#
# generate wrapper files for both versions of python. This will make sure sdist gets
# all files needed for the source package
#
# we re-generate wrapper files for all versions of python only for
# source distribution
if not os.path.isfile(ASSO_WRAPPER_PY_FILE) or not os.path.isfile(ASSO_WRAPPER_CPP_FILE):
    print('Generating {}'.format(ASSO_WRAPPER_CPP_FILE))
    print('swig ' + ' '.join(SWIG_OPTS + ['-o', ASSO_WRAPPER_CPP_FILE, ASSO_INTERFACE_FILE]))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', ASSO_WRAPPER_CPP_FILE, ASSO_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for association module.')
    os.rename('src/variant_tools/assoTests.py', ASSO_WRAPPER_PY_FILE)
#
if not os.path.isfile(CGATOOLS_WRAPPER_PY_FILE) or not os.path.isfile(CGATOOLS_WRAPPER_CPP_FILE):
    print('Generating {}'.format(CGATOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', CGATOOLS_WRAPPER_CPP_FILE, CGATOOLS_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for cgatools.')
    os.rename('src/variant_tools/cgatools.py', CGATOOLS_WRAPPER_PY_FILE)
#
if not os.path.isfile(UCSCTOOLS_WRAPPER_PY_FILE) or not os.path.isfile(UCSCTOOLS_WRAPPER_CPP_FILE):
    print('Generating {}'.format(UCSCTOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call('swig ' + ' '.join(SWIG_OPTS + ['-o', UCSCTOOLS_WRAPPER_CPP_FILE, UCSCTOOLS_INTERFACE_FILE]), shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for ucsctools.')
    os.rename('src/variant_tools/ucsctools.py', UCSCTOOLS_WRAPPER_PY_FILE)
         
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
            print('Building embedded boost library')
            objects = c.compile(LIB_BOOST,
                include_dirs=['src/boost_1_49_0'] + ENV_INCLUDE_DIRS,
                output_dir='build',
                extra_preargs = ['-w', '-fPIC'],
                macros = [('BOOST_ALL_NO_LIB', None)])
            c.create_static_lib(objects, "embedded_boost", output_dir='build')
    except Exception as e:
        sys.exit("Failed to build embedded boost library: {}".format(e))


if EMBEDDED_ZMQ:
  currentdir=os.getcwd()
  try:
    if not os.path.isfile("./src/zeromq-4.0.3/Makefile"):
      print('Building embedded zeromq library')
      os.chdir("./src/zeromq-4.0.3")
      ret=subprocess.call(["./configure", "--prefix="+currentdir+"/src/zeromq-4.0.3"],shell=True)
      if ret != 0: sys.exit('Failed to install zmq.')
      ret=subprocess.call(["make", "install"],shell=True)
      if ret != 0: sys.exit('Failed to install zmq.')
  except OSError as e:
        sys.exit('Failed to install zmq.')
  finally:
    os.chdir(currentdir)



# building other libraries
for files, incs, macs, libname in [
    (SQLITE_GSL, ['src', 'src/gsl'], [], 'sqlite_gsl'),
    (LIB_STAT, ['src'], [], 'stat'),
    (LIB_UCSC_FILES, ['src/ucsc/inc', 'src/ucsc/tabix', 'src/ucsc/samtools'], 
        [('USE_TABIX', '1'), ('_FILE_OFFSET_BITS', '64'), ('USE_BAM', '1'),
         ('_USE_KNETFILE', None), ('BGZF_CACHE', None)],
        'ucsc'),
    (LIB_CGATOOLS, ['src', 'src/cgatools', 'src/boost_1_49_0'],
        [('CGA_TOOLS_IS_PIPELINE', 0), ('CGA_TOOLS_VERSION', r'"1.6.0.43"')],
        'cgatools'),
    (LIB_GSL, ['src', 'src/gsl'], [], 'gsl')]:
    c = new_compiler()
    if os.path.isfile(os.path.join('build', c.static_lib_format % (libname, c.static_lib_extension))):
        continue
    try:
        print('Building embedded library {}'.format(libname))
        c = new_compiler()
        # -w suppress all warnings caused by the use of boost libraries
        objects = c.compile(files,
            include_dirs=incs  + ENV_INCLUDE_DIRS,
            output_dir='build',
            extra_preargs = ['-w', '-fPIC'],
            macros = macs)
        c.create_static_lib(objects, libname, output_dir='build')
    except Exception as e:
        sys.exit("Failed to build embedded {} library: {}".format(libname, e))

ext_modules=[
        Extension('variant_tools._vt_sqlite3',
            # stop warning message for sqlite because it is written by us.
            extra_compile_args=['-w'],
            sources = SQLITE_FILES,
            define_macros = [('MODULE_NAME', '"vt_sqlite3"'), ('HAVE_USLEEP', None)],
            include_dirs = ['src/sqlite', SQLITE_FOLDER] + ENV_INCLUDE_DIRS,
        ),
        Extension('variant_tools._ucsctools',
            # stop warning message ucsctools because it is written by us.
            extra_compile_args=['-w'],
            sources = [UCSCTOOLS_WRAPPER_CPP_FILE],
            include_dirs = ['.', 'src/ucsc/inc', 'src/ucsc/tabix', 'src/ucsc/samtools']  + ENV_INCLUDE_DIRS,
            library_dirs = ["build"]  + ENV_LIBRARY_DIRS,

            define_macros =  [('USE_TABIX', '1'), ('_FILE_OFFSET_BITS', '64'), ('USE_BAM', '1'),
                ('_USE_KNETFILE', None), ('BGZF_CACHE', None)],
            libraries = libs + ['ucsc', 'z', 'bz2'],
        ),
        Extension('variant_tools.cplinkio',
            # stop warning message ucsctools because it is written by us.
            extra_compile_args=['-w'],
            sources = LIB_PLINKIO,
            include_dirs = ['src/libplinkio']  + ENV_INCLUDE_DIRS,
        ),
        Extension('variant_tools._vt_sqlite3_ext',
            # stop warning message for sqlite because it is written by us.
            sources = ['src/sqlite/vt_sqlite3_ext.cpp'],
            include_dirs = ["src/", 'src/ucsc/inc', 'src/ucsc/tabix', 'src/ucsc/samtools',
                'src/sqlite', "src/variant_tools", "src/gsl", "src/cgatools", "src/boost_1_49_0"] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,
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

            include_dirs = ["src", "src/cgatools", "src/boost_1_49_0"] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,

        ),
        Extension('variant_tools._assoTests',
            sources = [ASSO_WRAPPER_CPP_FILE] + ASSOC_FILES,
            extra_compile_args = gccargs,
            libraries = libs + ['gsl', 'stat'], #, 'blas'],
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,
            include_dirs = ["src", "src/variant_tools", "src/gsl"] + ENV_INCLUDE_DIRS,

        )
      ]
ext_modules+=cythonize([Extension('variant_tools.io_vcf_read',
    # sources=['/Users/jma7/Development/VAT/VariantTools/src/variant_tools/io_vcf_read.c'],
    sources=['src/variant_tools/io_vcf_read.pyx'],
    include_dirs=[np.get_include()],
    library_dirs = ["build"])])


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
    install_requires=[
       'pyzmq'
    ],
    packages = find_packages('src'),
    package_dir = {'': 'src'},
    entry_points = '''
[console_scripts]
vtools = variant_tools.vtools:main
vtools_report = variant_tools.vtools_report:main
worker_run=variant_tools.worker_zmq:main
''',
    ext_modules = ext_modules
)

