#!/usr/bin/env python
#
# $File: __init__.py $
# $LastChangedDate$
# $Rev$
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://varianttools.sourceforge.net for details.
#
# Copyright (C) 2011 - 2013 - 2013 Bo Peng (bpeng@mdanderson.org)
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

import sys

ver = sys.version_info
if (ver.major == 2 and (ver.minor, ver.micro) < (7, 1)) or (ver.major == 3 and (ver.minor, ver.micro) < (2, 0)):
    raise SystemError('variant tools requires Python 2.7.1, Python 3.2.1 or higher. Please upgrade your Python {}.{}.{}.'.format(ver.major, ver.minor, ver.micro))


#__all__ = ['utils', 'project', 'variant', 'annotation', 'assoTests', 'pipeline', 'simulation']

#
# The vtools version string consists of three digits where the REVISION will be
# changed during development, especially with each addition or change of feature.
# The development version usually has non-zero revision number where a formal 
# release has .0 revision number.
#
# MAJOR.MINOR.REVISION
#
VTOOLS_VERSION='2.7.2'
#
# Variant Tools no longer uses revision information after switching from SVN to GIT
# the version string will be more frequently updated.
#
# VTOOLS_REVISION='2961'
#
VTOOLS_FULL_VERSION='{} for Python {}.{}.{}'.format(VTOOLS_VERSION, ver.major, ver.minor, ver.micro)
VTOOLS_COPYRIGHT = '''variant tools {} : Copyright (c) 2011 - 2016 Bo Peng'''.format(VTOOLS_VERSION)
VTOOLS_CONTACT = '''Please visit http://varianttools.sourceforge.net for more information.'''

try:
    from .utils import *
    from .preprocessor import *
    from .pipeline import *
except ImportError as e:
    print('Failed to import variant_tools. Make sure that variant tools is properly installed '
        'and you are not importing variant_tools from its source directory.')
