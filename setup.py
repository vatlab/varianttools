#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 Bo Peng (bpeng@mdanderson.org)
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
import os
import subprocess
import sys

import numpy as np
from Cython.Build import cythonize
from setuptools import Extension, find_packages, setup

# do not import variant_tools because __init__ might not be imported properly
# before installation
with open('src/variant_tools/_version.py') as init:
    for line in init:
        if line.startswith('VTOOLS_VERSION='):
            VTOOLS_VERSION = line[15:].strip().strip('"').strip("'")
            break

SWIG_OPTS = [
    '-c++', '-python', '-O', '-shadow', '-keyword', '-w-511', '-w-509',
    '-outdir', 'src/variant_tools', '-py3'
]

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
    'src/variant_tools/assoTests.i', 'src/variant_tools/assoTests.h',
    'src/variant_tools/assoData.h', 'src/variant_tools/action.h',
    'src/variant_tools/utils.h', 'src/variant_tools/lm.h'
]

ASSOC_FILES = [
    'src/variant_tools/assoData.cpp', 'src/variant_tools/action.cpp',
    'src/variant_tools/utils.cpp', 'src/variant_tools/lm.cpp'
]

SQLITE_FOLDER = 'src/sqlite/py3'
SQLITE_FILES = [
    os.path.join(SQLITE_FOLDER, x) for x in [
        'cache.c', 'connection.c', 'cursor.c', 'microprotocols.c', 'module.c',
        'prepare_protocol.c', 'row.c', 'statement.c', 'util.c'
    ]
] + ['src/sqlite/sqlite3.c']

PLINKIO_FILES = [
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

CGATOOLS_FILES = [
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
UCSC_FILES = [
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

STAT_FILES = ['src/variant_tools/fisher2.c']

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
        print(
            'The boost C++ library version 1.49.0 is not found under the current directory. Will try to use the system libraries.'
        )

#
# During development, if an interface file needs to be re-generated, please
# remove these files and they will be re-generated with SWIG
#
if not os.path.isfile('src/swigpyrun.h'):
    try:
        ret = subprocess.call(
            ['swig -python -external-runtime src/swigpyrun.h -o build'],
            shell=True)
        if ret != 0:
            sys.exit('Failed to generate swig runtime header file.')
    except OSError:
        sys.exit(
            'Failed to generate wrapper file. Please install swig (www.swig.org).'
        )
#
# generate wrapper files for both versions of python. This will make sure sdist gets
# all files needed for the source package
#
# we re-generate wrapper files for all versions of python only for
# source distribution
if not os.path.isfile(ASSO_WRAPPER_PY_FILE) or not os.path.isfile(
        ASSO_WRAPPER_CPP_FILE):
    print('Generating {}'.format(ASSO_WRAPPER_CPP_FILE))
    print('swig ' +
          ' '.join(SWIG_OPTS +
                   ['-o', ASSO_WRAPPER_CPP_FILE, ASSO_INTERFACE_FILE]))
    ret = subprocess.call(
        'swig ' + ' '.join(SWIG_OPTS +
                           ['-o', ASSO_WRAPPER_CPP_FILE, ASSO_INTERFACE_FILE]),
        shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for association module.')
    os.rename('src/variant_tools/assoTests.py', ASSO_WRAPPER_PY_FILE)
#
if not os.path.isfile(CGATOOLS_WRAPPER_PY_FILE) or not os.path.isfile(
        CGATOOLS_WRAPPER_CPP_FILE):
    print('Generating {}'.format(CGATOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call(
        'swig ' +
        ' '.join(SWIG_OPTS +
                 ['-o', CGATOOLS_WRAPPER_CPP_FILE, CGATOOLS_INTERFACE_FILE]),
        shell=True)
    if ret != 0:
        sys.exit('Failed to generate wrapper file for cgatools.')
    os.rename('src/variant_tools/cgatools.py', CGATOOLS_WRAPPER_PY_FILE)
#
if not os.path.isfile(UCSCTOOLS_WRAPPER_PY_FILE) or not os.path.isfile(
        UCSCTOOLS_WRAPPER_CPP_FILE):
    print('Generating {}'.format(UCSCTOOLS_WRAPPER_CPP_FILE))
    ret = subprocess.call(
        'swig ' +
        ' '.join(SWIG_OPTS +
                 ['-o', UCSCTOOLS_WRAPPER_CPP_FILE, UCSCTOOLS_INTERFACE_FILE]),
        shell=True)
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

ENV_INCLUDE_DIRS = [os.path.join(os.environ['CONDA_PREFIX'], 'include')
                   ] if 'CONDA_PREFIX' in os.environ else []
ENV_INCLUDE_DIRS += [
    x for x in os.environ.get('LD_INCLUDE_PATH', '').split(os.pathsep) if x
]
ENV_LIBRARY_DIRS = [
    x for x in os.environ.get('LD_LIBRARY_PATH', '').split(os.pathsep) if x
]

if EMBEDDED_BOOST:
    try:
        from distutils.ccompiler import new_compiler

        c = new_compiler()
        if not os.path.isfile(
                os.path.join(
                    'build', c.static_lib_format %
                    ('embedded_boost', c.static_lib_extension))):
            # -w suppress all warnings caused by the use of boost libraries
            print('Building embedded boost library')
            objects = c.compile(
                LIB_BOOST,
                include_dirs=['src/boost_1_49_0'] + ENV_INCLUDE_DIRS,
                output_dir='build',
                extra_preargs=['-w', '-fPIC'],
                macros=[('BOOST_ALL_NO_LIB', None)])
            c.create_static_lib(objects, "embedded_boost", output_dir='build')
    except Exception as e:
        sys.exit("Failed to build embedded boost library: {}".format(e))


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
            sources = [UCSCTOOLS_WRAPPER_CPP_FILE] + UCSC_FILES,
            include_dirs = ['.', 'src/ucsc/inc', 'src/ucsc/tabix', 'src/ucsc/samtools']  + ENV_INCLUDE_DIRS,
            library_dirs = ["build"]  + ENV_LIBRARY_DIRS,

            define_macros =  [('USE_TABIX', '1'), ('_FILE_OFFSET_BITS', '64'), ('USE_BAM', '1'),
                ('_USE_KNETFILE', None), ('BGZF_CACHE', None)],
            libraries = libs + ['z', 'bz2'],
        ),
        Extension('variant_tools.cplinkio',
            # stop warning message ucsctools because it is written by us.
            extra_compile_args=['-w'],
            sources = PLINKIO_FILES,
            include_dirs = ['src/libplinkio']  + ENV_INCLUDE_DIRS,
        ),
        Extension('variant_tools._vt_sqlite3_ext',
            # stop warning message for sqlite because it is written by us.
            sources = ['src/sqlite/vt_sqlite3_ext.cpp'] + CGATOOLS_FILES + STAT_FILES + UCSC_FILES,
            include_dirs = ["src/", 'src/ucsc/inc', 'src/ucsc/tabix', 'src/ucsc/samtools',
                'src/sqlite', "src/variant_tools", "src/cgatools"] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,
            libraries = ['gsl'] + \
                (['embedded_boost'] if EMBEDDED_BOOST else ['boost_iostreams', 'boost_regex', 'boost_filesystem']) + \
                ['z', 'bz2'],
            extra_compile_args = gccargs,
            define_macros = [
                ('MODULE_NAME', '"vt_sqlite3"'),
                ('CGA_TOOLS_IS_PIPELINE', 0),
                ('CGA_TOOLS_VERSION', r'"1.6.0.43"'), ('USE_TABIX', '1'), ('USE_BAM', '1'),
                ('_FILE_OFFSET_BITS', '64'), ('_USE_KNETFILE', None),
                ('BGZF_CACHE', None)],
        ),
        Extension('variant_tools._cgatools',
            sources = [CGATOOLS_WRAPPER_CPP_FILE] + CGATOOLS_FILES,
            libraries = \
                (['embedded_boost'] if EMBEDDED_BOOST else ['boost_iostreams', 'boost_regex', 'boost_filesystem']) + \
                ['z', 'bz2'] + ENV_INCLUDE_DIRS,
            define_macros = [('CGA_TOOLS_IS_PIPELINE', 0),
                ('CGA_TOOLS_VERSION', r'"1.6.0.43"')],
            extra_compile_args = gccargs,
            swig_opts = ['-O', '-shadow', '-c++', '-keyword'],

            include_dirs = ["src", "src/cgatools"] + ENV_INCLUDE_DIRS,
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,

        ),
        Extension('variant_tools._assoTests',
            sources = [ASSO_WRAPPER_CPP_FILE] + STAT_FILES + ASSOC_FILES,
            extra_compile_args = gccargs,
            libraries = libs + ['gsl', 'blas'],
            library_dirs = ["build"] + ENV_LIBRARY_DIRS,
            include_dirs = ["src", "src/variant_tools"] + ENV_INCLUDE_DIRS,

        )
      ]
ext_modules += cythonize([
    Extension(
        'variant_tools.io_vcf_read',
        # sources=['/Users/jma7/Development/VAT/VariantTools/src/variant_tools/io_vcf_read.c'],
        sources=['src/variant_tools/io_vcf_read.pyx'],
        include_dirs=[np.get_include()],
        library_dirs=["build"])
])

CURRENT_DIR = os.path.abspath(os.path.dirname(__file__))


def get_long_description():
    with open(os.path.join(CURRENT_DIR, "README.md"), "r") as ld_file:
        return ld_file.read()


setup(
    name="variant_tools",
    version=VTOOLS_VERSION,
    description="Variant tools: an integrated annotation and analysis package for next-generation sequencing data",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author='Bo Peng',
    url='http://varianttools.sourceforge.net',
    author_email='bpeng@mdanderson.org',
    maintainer='Bo Peng',
    maintainer_email='varianttools-devel@lists.sourceforge.net',
    license='GPL3',
    classifiers=[
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
        'pyzmq',
        'tables',
        'numpy',
        'Cython',
        'scipy',
        'pyzmq',
        'pycurl',
        'tables',
    ],
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points='''
[console_scripts]
vtools = variant_tools.vtools:main
vtools_report = variant_tools.vtools_report:main
worker_run=variant_tools.worker_zmq:main
''',
    ext_modules=ext_modules)
