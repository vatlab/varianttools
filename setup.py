#!/usr/bin/env python2.7
#
# $File: setup.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
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
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

import sys, os
try:
    import argparse
except ImportError:
    sys.exit('variant tools requires Python 2.7.2 or higher, or Python 3.2 or higher. Please upgrade your version (%s) of Python and try again.' % (sys.version.split()[0]))

from variant_tools import VTOOLS_VERSION

SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword', '-w-511',
    '-outdir', 'variant_tools']

if sys.version_info.major == 2:
    WRAPPER_CPP_FILE = 'variant_tools/assoTests_wrap_py2.cpp'
    WRAPPER_PY_FILE = 'variant_tools/assoTests_py2.py'
else:
    SWIG_OPTS.append('-py3')
    WRAPPER_CPP_FILE = 'variant_tools/assoTests_wrap_py3.cpp'
    WRAPPER_PY_FILE = 'variant_tools/assoTests_py3.py'


# generate wrapper files
if not os.path.isfile(WRAPPER_PY_FILE) or not os.path.isfile(WRAPPER_CPP_FILE) or os.path.getmtime(WRAPPER_CPP_FILE) < \
        max([os.path.getmtime(x) for x in ('variant_tools/assoTests.i', 'variant_tools/association.py')]):
    import subprocess
    print('Generating wrapper file')
    try:
        ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', WRAPPER_CPP_FILE, 'variant_tools/assoTests.i'], shell=False)
        if ret != 0:
            sys.exit('Failed to generate wrapper file. Please install swig.')
        os.rename('variant_tools/assoTests.py', WRAPPER_PY_FILE)
    except OSError as e:
        sys.exit('Failed to generate wrapper file. Please install swig: {}'.format(e))


LIB_GSL = [ 
   'variant_tools/gsl/error.c',
   'variant_tools/gsl/sys/infnan.c',
   'variant_tools/gsl/sys/coerce.c',
   'variant_tools/gsl/sys/fdiv.c',
   'variant_tools/gsl/sys/pow_int.c',
   'variant_tools/gsl/sys/fcmp.c',
   'variant_tools/gsl/sys/log1p.c',
   'variant_tools/gsl/sys/invhyp.c',
   'variant_tools/gsl/complex/math.c',
   'variant_tools/gsl/specfunc/beta.c',
   'variant_tools/gsl/specfunc/psi.c',
   'variant_tools/gsl/specfunc/trig.c',
   'variant_tools/gsl/specfunc/exp.c',
   'variant_tools/gsl/specfunc/expint.c',
   'variant_tools/gsl/specfunc/log.c',
   'variant_tools/gsl/specfunc/erfc.c',
   'variant_tools/gsl/specfunc/zeta.c',
   'variant_tools/gsl/specfunc/elementary.c',
   'variant_tools/gsl/specfunc/gamma.c',
   'variant_tools/gsl/specfunc/gamma_inc.c',
   'variant_tools/gsl/rng/rng.c',
   'variant_tools/gsl/rng/default.c',
   'variant_tools/gsl/rng/mt.c',
   'variant_tools/gsl/rng/types.c',
   'variant_tools/gsl/randist/binomial.c',
   'variant_tools/gsl/randist/binomial_tpe.c',
   'variant_tools/gsl/randist/beta.c',
   'variant_tools/gsl/randist/exponential.c',
   'variant_tools/gsl/randist/geometric.c',
   'variant_tools/gsl/randist/nbinomial.c',
   'variant_tools/gsl/randist/poisson.c',
   'variant_tools/gsl/randist/multinomial.c',
   'variant_tools/gsl/randist/chisq.c',
   'variant_tools/gsl/randist/gauss.c',
   'variant_tools/gsl/randist/tdist.c',
   'variant_tools/gsl/randist/gausszig.c',
   'variant_tools/gsl/randist/gamma.c',
   'variant_tools/gsl/randist/hyperg.c',
   'variant_tools/gsl/cdf/binomial.c',
   'variant_tools/gsl/cdf/beta.c',
   'variant_tools/gsl/cdf/betainv.c',
   'variant_tools/gsl/cdf/gauss.c',
   'variant_tools/gsl/cdf/gaussinv.c',
   'variant_tools/gsl/cdf/tdist.c',
   'variant_tools/gsl/cdf/tdistinv.c',
   'variant_tools/gsl/cdf/chisq.c',
   'variant_tools/gsl/cdf/chisqinv.c',
   'variant_tools/gsl/cdf/gamma.c',
   'variant_tools/gsl/cdf/gammainv.c',
   'variant_tools/gsl/cdf/hypergeometric.c',
   'variant_tools/gsl/cdf/poisson.c',
   #
   'variant_tools/gsl/blas/blas.c',
   'variant_tools/gsl/cblas/caxpy.c',
   'variant_tools/gsl/cblas/ccopy.c',
   'variant_tools/gsl/cblas/cdotc_sub.c',
   'variant_tools/gsl/cblas/cdotu_sub.c',
   'variant_tools/gsl/cblas/cgbmv.c',
   'variant_tools/gsl/cblas/cgemm.c',
   'variant_tools/gsl/cblas/cgemv.c',
   'variant_tools/gsl/cblas/cgerc.c',
   'variant_tools/gsl/cblas/cgeru.c',
   'variant_tools/gsl/cblas/chbmv.c',
   'variant_tools/gsl/cblas/chemm.c',
   'variant_tools/gsl/cblas/chemv.c',
   'variant_tools/gsl/cblas/cher2.c',
   'variant_tools/gsl/cblas/cher2k.c',
   'variant_tools/gsl/cblas/cher.c',
   'variant_tools/gsl/cblas/cherk.c',
   'variant_tools/gsl/cblas/chpmv.c',
   'variant_tools/gsl/cblas/chpr2.c',
   'variant_tools/gsl/cblas/chpr.c',
   'variant_tools/gsl/cblas/cscal.c',
   'variant_tools/gsl/cblas/csscal.c',
   'variant_tools/gsl/cblas/cswap.c',
   'variant_tools/gsl/cblas/csymm.c',
   'variant_tools/gsl/cblas/csyr2k.c',
   'variant_tools/gsl/cblas/csyrk.c',
   'variant_tools/gsl/cblas/ctbmv.c',
   'variant_tools/gsl/cblas/ctbsv.c',
   'variant_tools/gsl/cblas/ctpmv.c',
   'variant_tools/gsl/cblas/ctpsv.c',
   'variant_tools/gsl/cblas/ctrmm.c',
   'variant_tools/gsl/cblas/ctrmv.c',
   'variant_tools/gsl/cblas/ctrsm.c',
   'variant_tools/gsl/cblas/ctrsv.c',
   'variant_tools/gsl/cblas/dasum.c',
   'variant_tools/gsl/cblas/daxpy.c',
   'variant_tools/gsl/cblas/dcopy.c',
   'variant_tools/gsl/cblas/ddot.c',
   'variant_tools/gsl/cblas/dgbmv.c',
   'variant_tools/gsl/cblas/dgemm.c',
   'variant_tools/gsl/cblas/dgemv.c',
   'variant_tools/gsl/cblas/dger.c',
   'variant_tools/gsl/cblas/dnrm2.c',
   'variant_tools/gsl/cblas/drot.c',
   'variant_tools/gsl/cblas/drotg.c',
   'variant_tools/gsl/cblas/drotm.c',
   'variant_tools/gsl/cblas/drotmg.c',
   'variant_tools/gsl/cblas/dsbmv.c',
   'variant_tools/gsl/cblas/dscal.c',
   'variant_tools/gsl/cblas/dsdot.c',
   'variant_tools/gsl/cblas/dspmv.c',
   'variant_tools/gsl/cblas/dspr2.c',
   'variant_tools/gsl/cblas/dspr.c',
   'variant_tools/gsl/cblas/dswap.c',
   'variant_tools/gsl/cblas/dsymm.c',
   'variant_tools/gsl/cblas/dsymv.c',
   'variant_tools/gsl/cblas/dsyr2.c',
   'variant_tools/gsl/cblas/dsyr2k.c',
   'variant_tools/gsl/cblas/dsyr.c',
   'variant_tools/gsl/cblas/dsyrk.c',
   'variant_tools/gsl/cblas/dtbmv.c',
   'variant_tools/gsl/cblas/dtbsv.c',
   'variant_tools/gsl/cblas/dtpmv.c',
   'variant_tools/gsl/cblas/dtpsv.c',
   'variant_tools/gsl/cblas/dtrmm.c',
   'variant_tools/gsl/cblas/dtrmv.c',
   'variant_tools/gsl/cblas/dtrsm.c',
   'variant_tools/gsl/cblas/dtrsv.c',
   'variant_tools/gsl/cblas/dzasum.c',
   'variant_tools/gsl/cblas/dznrm2.c',
   'variant_tools/gsl/cblas/hypot.c',
   'variant_tools/gsl/cblas/icamax.c',
   'variant_tools/gsl/cblas/idamax.c',
   'variant_tools/gsl/cblas/isamax.c',
   'variant_tools/gsl/cblas/izamax.c',
   'variant_tools/gsl/cblas/sasum.c',
   'variant_tools/gsl/cblas/saxpy.c',
   'variant_tools/gsl/cblas/scasum.c',
   'variant_tools/gsl/cblas/scnrm2.c',
   'variant_tools/gsl/cblas/scopy.c',
   'variant_tools/gsl/cblas/sdot.c',
   'variant_tools/gsl/cblas/sdsdot.c',
   'variant_tools/gsl/cblas/sgbmv.c',
   'variant_tools/gsl/cblas/sgemm.c',
   'variant_tools/gsl/cblas/sgemv.c',
   'variant_tools/gsl/cblas/sger.c',
   'variant_tools/gsl/cblas/snrm2.c',
   'variant_tools/gsl/cblas/srot.c',
   'variant_tools/gsl/cblas/srotg.c',
   'variant_tools/gsl/cblas/srotm.c',
   'variant_tools/gsl/cblas/srotmg.c',
   'variant_tools/gsl/cblas/ssbmv.c',
   'variant_tools/gsl/cblas/sscal.c',
   'variant_tools/gsl/cblas/sspmv.c',
   'variant_tools/gsl/cblas/sspr2.c',
   'variant_tools/gsl/cblas/sspr.c',
   'variant_tools/gsl/cblas/sswap.c',
   'variant_tools/gsl/cblas/ssymm.c',
   'variant_tools/gsl/cblas/ssymv.c',
   'variant_tools/gsl/cblas/ssyr2.c',
   'variant_tools/gsl/cblas/ssyr2k.c',
   'variant_tools/gsl/cblas/ssyr.c',
   'variant_tools/gsl/cblas/ssyrk.c',
   'variant_tools/gsl/cblas/stbmv.c',
   'variant_tools/gsl/cblas/stbsv.c',
   'variant_tools/gsl/cblas/stpmv.c',
   'variant_tools/gsl/cblas/stpsv.c',
   'variant_tools/gsl/cblas/strmm.c',
   'variant_tools/gsl/cblas/strmv.c',
   'variant_tools/gsl/cblas/strsm.c',
   'variant_tools/gsl/cblas/strsv.c',
   'variant_tools/gsl/cblas/xerbla.c',
   'variant_tools/gsl/cblas/zaxpy.c',
   'variant_tools/gsl/cblas/zcopy.c',
   'variant_tools/gsl/cblas/zdotc_sub.c',
   'variant_tools/gsl/cblas/zdotu_sub.c',
   'variant_tools/gsl/cblas/zdscal.c',
   'variant_tools/gsl/cblas/zgbmv.c',
   'variant_tools/gsl/cblas/zgemm.c',
   'variant_tools/gsl/cblas/zgemv.c',
   'variant_tools/gsl/cblas/zgerc.c',
   'variant_tools/gsl/cblas/zgeru.c',
   'variant_tools/gsl/cblas/zhbmv.c',
   'variant_tools/gsl/cblas/zhemm.c',
   'variant_tools/gsl/cblas/zhemv.c',
   'variant_tools/gsl/cblas/zher2.c',
   'variant_tools/gsl/cblas/zher2k.c',
   'variant_tools/gsl/cblas/zher.c',
   'variant_tools/gsl/cblas/zherk.c',
   'variant_tools/gsl/cblas/zhpmv.c',
   'variant_tools/gsl/cblas/zhpr2.c',
   'variant_tools/gsl/cblas/zhpr.c',
   'variant_tools/gsl/cblas/zscal.c',
   'variant_tools/gsl/cblas/zswap.c',
   'variant_tools/gsl/cblas/zsymm.c',
   'variant_tools/gsl/cblas/zsyr2k.c',
   'variant_tools/gsl/cblas/zsyrk.c',
   'variant_tools/gsl/cblas/ztbmv.c',
   'variant_tools/gsl/cblas/ztbsv.c',
   'variant_tools/gsl/cblas/ztpmv.c',
   'variant_tools/gsl/cblas/ztpsv.c',
   'variant_tools/gsl/cblas/ztrmm.c',
   'variant_tools/gsl/cblas/ztrmv.c',
   'variant_tools/gsl/cblas/ztrsm.c',
   'variant_tools/gsl/cblas/ztrsv.c',
   #
   'variant_tools/gsl/linalg/svd.c',
   'variant_tools/gsl/linalg/bidiag.c',
   'variant_tools/gsl/linalg/householder.c',
   'variant_tools/gsl/matrix/matrix.c',
   'variant_tools/gsl/matrix/submatrix.c',
   'variant_tools/gsl/matrix/rowcol.c',
   'variant_tools/gsl/matrix/getset.c',
   'variant_tools/gsl/matrix/init.c',
   'variant_tools/gsl/vector/init.c',
   'variant_tools/gsl/matrix/swap.c',
   'variant_tools/gsl/vector/vector.c',
   'variant_tools/gsl/vector/copy.c',
   'variant_tools/gsl/vector/swap.c',
   'variant_tools/gsl/vector/subvector.c',
   'variant_tools/gsl/vector/oper.c',
   'variant_tools/gsl/matrix/oper.c',
   'variant_tools/gsl/matrix/copy.c',
   'variant_tools/gsl/block/init.c'
    ] 
    
# Under linux/gcc, lib stdc++ is needed for C++ based extension.
if sys.platform == 'linux2':
    libs = ['stdc++']
else:
    libs = []
  

setup(name = "variant_tools",
    version = VTOOLS_VERSION,
    description = "Variant tools: an integrated annotation and analysis package for next-gen sequencing data",
    author = 'Bo Peng',
    url = 'http://varianttools.sourceforge.net',
    author_email = 'bpeng@mdanderson.org',
    maintainer = 'Bo Peng',
    maintainer_email = 'varianttools-devel@lists.sourceforge.net',
    py_modules = [
        'variant_tools.__init__',
        'variant_tools.utils',
        'variant_tools.project',
        'variant_tools.importer',
        'variant_tools.exporter',
        'variant_tools.phenotype',
        'variant_tools.variant',
        'variant_tools.annotation',
        'variant_tools.liftOver',
        'variant_tools.association',
        WRAPPER_PY_FILE[:-3]          # assotests_pyX.py file without extension
    ],
    scripts = [
        'vtools',
        'vtools_report'
    ],
    cmdclass = {'build_py': build_py },
    package_dir = {'variant_tools': 'variant_tools'},
    ext_modules = [
        Extension('variant_tools/_assoTests',
            sources = [
                WRAPPER_CPP_FILE,
                'variant_tools/assoData.cpp',
                'variant_tools/utils.cpp'
                ] + LIB_GSL,
            library_dirs = [],
            libraries = libs,
            include_dirs = [".", "variant_tools", "variant_tools/gsl"],
        )
      ]
)

