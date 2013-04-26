#!/bin/sh
#
# $File: postinstall.sh $
# $LastChangedDate: 2013-04-16 13:32:03 -0500 (Tue, 16 Apr 2013) $
# $Rev: 1825 $
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

# install vtools and vtools_report to /usr/local/bin
#
if [ -f /usr/local/bin/vtools ]
then
    /bin/rm -f /usr/local/bin/vtools
fi

if [ -f /usr/local/bin/vtools_report ]
then
    /bin/rm -f /usr/local/bin/vtools_report
fi

/bin/ln -s /Applications/variant_tools/variant_tools.app/Contents/MacOS/vtools /usr/local/bin
/bin/ln -s /Applications/variant_tools/variant_tools.app/Contents/MacOS/vtools_report /usr/local/bin

exit 0
