#!/usr/bin/env python
#
# $File: release.py $
# $LastChangedDate: 2013-01-30 18:29:35 -0600 (Wed, 30 Jan 2013) $
# $Rev: 1663 $
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

#
# use pyinstaller to build binaries. This script relies on other tools
# such as git and is only intended to be used by developers of variant
# tools.
#

import os
import sys
import subprocess
import shutil
import argparse
import zipfile
import platform
from variant_tools.utils import runCommand

if sys.version.startswith('2.7.4') or sys.version.startswith('3.3.1'):
    sys.exit('Python 2.7.4 and 3.3.1 cannot be used to release variant '
        'tools due to a regression bug in the gzip module.')

def modifyVersion(version):
    # modify source/__init__.py to write version string
    if version is not None:
        content = []
        rev = runCommand('svnversion').strip()
        with open('source/__init__.py', 'r') as init_file:
            for x in init_file.readlines():
                if x.startswith('VTOOLS_VERSION'):
                    content.append("VTOOLS_VERSION='{}'\n".format(version))
                elif x.startswith('# $Rev:'):
                    content.append("# $Rev: {} $\n".\
                                   format(rev))
                elif x.startswith("VTOOLS_REVISION='$Rev:"):
                    content.append("VTOOLS_REVISION='$Rev: {} $'\n".\
                                   format(rev))
                else:
                    content.append(x)
        with open('source/__init__.py', 'w') as init_file:
            init_file.write(''.join(content))
    # read version string
    from source import VTOOLS_VERSION
    return VTOOLS_VERSION

def setupEnvironment(version):
    #
    if version.endswith('svn'):
        print('WARNING: You are releasing a subversion version of variant tools.')
        print('To make a formal release, you will need to change version string in source/__init__.py')
    #
    if os.path.isdir('build'):
        print('Removing directory build')
        shutil.rmtree('build')
    #
    if os.path.isdir('variant_tools-{}'.format(version)):
        print('Removing directory variant_tools-{}'.format(version))
        shutil.rmtree('variant_tools-{}'.format(version))

def generateSWIGWrappers():
    # 
    # generate wrapper files to make sure they are up to date
    SWIG_OPTS = ['-c++', '-python', '-O', '-shadow', '-keyword', '-w-511',
        '-w-509', '-outdir', 'source']
    WRAPPER_CPP_FILE = 'source/assoTests_wrap_{}.cpp'
    WRAPPER_PY_FILE = 'source/assoTests_{}.py'
    CGATOOLS_WRAPPER_CPP_FILE = 'source/cgatools_wrap_{}.cpp'
    CGATOOLS_WRAPPER_PY_FILE = 'source/cgatools_{}.py'
    SQLITE_FOLDER = 'sqlite/{}'
    SQLITE_PY_FILE = 'source/vt_sqlite3_{}'
    with open(os.devnull, 'w') as fnull:
        try:
           print('Generating source/swigpyrun.h ...')
           ret = subprocess.call('swig -python -external-runtime source/swigpyrun.h', shell=True, stdout=fnull)
           if ret != 0:
               sys.exit('Failed to generate swig runtime header file.')
        except OSError as e:
            sys.exit('Failed to generate wrapper file. Please install swig (www.swig.org).')
        # 
        for ver in ['py2', 'py3']:
            print('Generating {} wrapper files for module assoTests ...'.format(ver))
            if ver == 'py3':
                SWIG_OPTS.append('-py3')
            ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', WRAPPER_CPP_FILE.format(ver), 'source/assoTests.i'], shell=False, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to generate wrapper file for association module.')
            os.rename('source/assoTests.py', WRAPPER_PY_FILE.format(ver))
            #
            print('Generating {} wrapper files for module cgatools ...'.format(ver))
            ret = subprocess.call(['swig'] + SWIG_OPTS + ['-o', CGATOOLS_WRAPPER_CPP_FILE.format(ver), 'source/cgatools.i'], shell=False, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to generate wrapper file for cgatools.')
            os.rename('source/cgatools.py', CGATOOLS_WRAPPER_PY_FILE.format(ver))
             

def buildVariantTools(extra_args):
    # build variant tools
    try:
        print('Building and installing variant tools ...')
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call('python setup.py install ' + ' '.join(extra_args), shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to build and install variant tools.')
    except Exception as e:
        sys.exit('Failed to build and install variant tools: {}'.format(e))


def buildSourcePackage():
    # build source package
    try:
        print('Building source package of variant tools {} ...'.format(version))
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call('python setup.py sdist', shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to build source package of variant tools.')
    except Exception as e:
        sys.exit('Failed to build source pacakge of variant tools: {}'.format(e))


def obtainPyInstaller(pyinstaller_dir):
    # check if pyinstaller is available
    #   if not, use git clone to get it
    #   if yes, try to update to the newest version
    pyinstaller_dir = os.path.expanduser(pyinstaller_dir.rstrip('/'))
    if pyinstaller_dir.endswith('pyinstaller'): pyinstaller_dir = pyinstaller_dir[:-12]
    if not os.path.isdir(pyinstaller_dir): pyinstaller_dir = os.getcwd()
    git_dir = os.path.join(pyinstaller_dir, 'pyinstaller')
    curdir = os.getcwd()
    if not os.path.isdir(git_dir):
        try:
            print('Downloading pyinstaller...')
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('git clone git://github.com/pyinstaller/pyinstaller.git {}'.format(git_dir), shell=True, stdout=fnull)
                if ret != 0:
                    sys.exit('Failed to clone pyinstaller. Please check if you have git installed.'
                        'You can also get pyinstaller manually anf decompress it under the pyinstaller directory.')
        except Exception as e:
            sys.exit('Failed to clone pyinstaller: {}'.format(e) + 
                'You can get pyinstaller manually anf decompress it under the pyinstaller directory '
                'if you are having trouble getting git installed.')
    else:
        os.chdir(git_dir)
        try:
            print('Updating pyinstaller ...')
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('git pull', shell=True, stdout=fnull)
                if ret != 0:
                    print('Failed to get latest version of pyinstaller. Using existing version.')
        except Exception as e:
            print('Failed to get latest version of pyinstaller ({}). Using existing version.'.format(e))
    os.chdir(curdir)
    return git_dir

def buildExecutables(git_dir, option):
    # use py installer to create executable
    for exe in ['vtools', 'vtools_report']:
        try:
            # remove existing files or directories
            dest = os.path.join('dist', exe)
            if os.path.isfile(dest):
                os.remove(dest)
            elif os.path.isdir(dest):
                shutil.rmtree(dest)
            if option == '--windowed':
                dest = os.path.join('dist', exe + '.app')
                if os.path.isdir(dest):
                    shutil.rmtree(dest)
            print('Building executable {} ...'.format(exe))
            with open(os.devnull, 'w') as fnull:
                ret = subprocess.call('python {} {} --log-level=ERROR {} '
                    .format(os.path.join(git_dir, 'pyinstaller.py'), option, exe),
                    shell=True, stdout=fnull)
                if ret != 0:
                    sys.exit('Failed to create executable for command {}'.format(exe))
        except Exception as e:
            sys.exit('Failed to create executable for command {}: {}'.format(exe, e))

def createZipPackage(version):
    # after the creation of commands, create a zip file with OS and version information
    zipfilename = os.path.join('dist', 'variant_tools-{}.{}.{}.zip'
        .format(version, platform.system(), platform.machine()))
    print('Adding executables to file {}'.format(zipfilename))
    with zipfile.ZipFile(zipfilename, 'w') as dist_file:
        vtools_cmd = 'vtools.exe' if os.name == 'win32' else 'vtools'
        dist_file.write(os.path.join('dist', vtools_cmd), vtools_cmd)
        vtools_report_cmd = 'vtools_report.exe' if os.name == 'win32' else 'vtools_report'
        dist_file.write(os.path.join('dist', vtools_report_cmd), vtools_report_cmd)

def createMacPackage(version):
    #
    # package source, required by PackageMaker
    src = os.path.join('dist', 'variant_tools')
    if os.path.isdir(src):
        shutil.rmtree(src)
    os.makedirs(src)
    # move files into src
    shutil.copy('README', src)
    shutil.copytree(os.path.join('dist', 'vtools.app'), os.path.join(src, 'variant_tools.app'))
    # merge vtools_report.app to variant_tools.app
    for root, dirs, files in os.walk(os.path.join('dist', 'vtools_report.app')):
        # root is dist/vtools_report.app/OTHERS
        # the destination should be dist/variant_tools/variant_tools.app/OTHERS
        new_root = root.replace('dist/vtools_report.app', os.path.join(src, 'variant_tools.app'))
        for f in files:
            if not os.path.isfile(os.path.join(new_root, f)):
                print('Copying {} to {}'.format(
                    os.path.join(root, f), os.path.join(new_root, f)))
                shutil.copy(os.path.join(root, f), os.path.join(new_root, f))
    #
    # package destination, within a dmg directory
    dest = os.path.join('dist', 'variant_tools-{}'.format(version))
    if os.path.isdir(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    # create a directory variant_tools and put everything inside it
    pkg = os.path.join(dest, 'variant_tools-{}.pkg'.format(version))
    # running packagemaker
    try:
        print('Building MacOSX package variant_tools-{}.pkg ...'.format(version))
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call(
                '/Applications/PackageMaker.app/Contents/MacOS/PackageMaker '
                '--verbose '
                '--doc development/MacOSX/variant_tools.pmdoc '
                '--out {}'.format(pkg),
                shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to create MacOSX package variant_tools-{}.pkg'
                    .format(version))
    except Exception as e:
        sys.exit('Failed to create MacOSX package variant_tools-{}.pkg: {}'
            .format(version, e))
    #
    # copy things into dmg
    shutil.copy('README', dest)
    shutil.copy('development/MacOSX/INSTALL', dest)
    for exe in [os.path.join('dist', x) for x in ['vtools', 'vtools_report']]:
        shutil.copy(exe, dest)
    dmg = os.path.join('dist', 'variant_tools-{}.dmg'.format(version))
    if os.path.isfile(dmg):
        os.remove(dmg)
    try:
        print('Building disk image variant_tools-{}.dmg ...'.format(version))
        with open(os.devnull, 'w') as fnull:
            ret = subprocess.call(
                'hdiutil create {} -volname variant_tools-{} -fs HFS+ -srcfolder {}'
                .format(dmg, version, dest), shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to create MacOSX disk image variant_tools-{}.dmg'.format(version))
    except Exception as e:
        sys.exit('Failed to create MacOSX disk image variant_tools-{}.dmg: {}'.format(version, e))

def createLinuxPackage(version):
    bundle = 'variant_tools-{}-{}-{}'.format(version, platform.system(), platform.machine())
    dest = os.path.join('dist', bundle)
    # merge the installations
    if os.path.isdir(dest):
        shutil.rmtree(dest)
    os.makedirs(dest)
    for folder in ['dist/vtools', 'dist/vtools_report']:
        os.system('rsync -auc {0}/* {1}'.format(folder, dest))
        shutil.rmtree(folder)
    with open('development/Linux/extractor-template.sh', 'r') as f:
        content = ["BUNDLE={0}\n".format(bundle) if x.startswith("BUNDLE=") else x for x in f.readlines()]
    with open(dest + '.sh', 'w') as f:
        f.write(''.join(content))
    shutil.copy('development/Linux/install.sh', dest)
    #
    os.system('cd dist; tar czf - {0} >> {0}.sh; cd -'.format(bundle))
    
def tagRelease(version):
    try:
        ret = subprocess.check_output(['svn', 'diff'], shell=True)
        if ret:
            sys.exit('Cannot tag release because there is uncommitted changes. '
                'Please commit the changes and try again.')
        with open(os.devnull, 'w') as fnull:
            print('Tagging release {}...'.format(version))
            ret = subprocess.call('svn copy svn+ssh://bpeng2000@svn.code.sf.net/p/varianttools/code/trunk '
                'svn+ssh://bpeng2000@svn.code.sf.net/p/varianttools/code/tag/v{} '
                ' -m "Version {} released at {}"'.format(version, version, time.asctime()),
                shell=True, stdout=fnull)
            if ret != 0:
                sys.exit('Failed to tag release {}.'.format(version))
    except Exception as e:
        sys.exit('Failed to tag release {}: {}'.format(version, e))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Create source distribution 
        and executables for a variant tools release. In addition to optional
        parameters version and tag, extra parameters would be specified and 
        will be passed directly to the 'python setup.py install' process. ''')
    parser.add_argument('--version',
        help='''Modify source/__init__.py to the specified version string and
            make the release.''')
    parser.add_argument('--tag', action='store_true',
        help='If specified, tag this release.')
    parser.add_argument('--pyinstaller_dir', default = '.',
        help='path to the directory where pyinstaller git clone is located.')
    parser.add_argument('--skip_rebuild', action='store_true',
        help=argparse.SUPPRESS)
    parser.add_argument('--skip_onefile', action='store_true',
        help=argparse.SUPPRESS)
    #
    # go to the root source directory
    os.chdir('..')
    # allow recognied parameters to be set to the build process
    args, argv = parser.parse_known_args()
    #
    version = modifyVersion(args.version)
    if not args.skip_rebuild:
        setupEnvironment(version)
        generateSWIGWrappers()
    buildVariantTools(argv)
    buildSourcePackage()
    git_dir = obtainPyInstaller(args.pyinstaller_dir)
    if platform.platform().startswith('Darwin'):
        # build mac mpkg package
        buildExecutables(git_dir, '--windowed')
        createMacPackage(version)
    else:
        # linux build self installation bundle
        buildExecutables(git_dir, '--onedir')
        createLinuxPackage(version)
    if not args.skip_onefile:
        buildExecutables(git_dir, '--onefile')
        createZipPackage(version)
    # if everything is OK, tag the release
    if args.tag:
        tagRelease(version)
    # if everything is done
    print('source packages and executables are successfully generated and '
        'saved to directory dist.')
