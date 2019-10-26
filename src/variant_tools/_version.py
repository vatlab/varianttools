#!/usr/bin/env python
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit https://github.com/vatlab/varianttools for details.
#
# Copyright (C) 2011 - 2020 - 2013 Bo Peng (bpeng@mdanderson.org)
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

VTOOLS_VERSION='3.0.7'

pyver = sys.version_info
VTOOLS_FULL_VERSION='{} for Python {}.{}.{}'.format(VTOOLS_VERSION, pyver.major, pyver.minor, pyver.micro)
VTOOLS_COPYRIGHT = '''variant tools {} : Copyright (c) 2011 - 2016 Bo Peng'''.format(VTOOLS_VERSION)
VTOOLS_CONTACT = '''Please visit https://github.com/vatlab/varianttools for more information.'''
