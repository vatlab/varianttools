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

import math
import multiprocessing
import random
import sys
import time

import numpy as np
import tables as tb


def merge(*args):
    left, right = args[0] if len(args) == 1 else args
    left_length, right_length = len(left), len(right)
    left_index, right_index = 0, 0
    merged = []
    while left_index < left_length and right_index < right_length:
        if left[left_index] <= right[right_index]:
            merged.append(left[left_index])
            left_index += 1
        else:
            merged.append(right[right_index])
            right_index += 1
    if left_index == left_length:
        merged.extend(right[right_index:])
    else:
        merged.extend(left[left_index:])
    return merged


def merge_sort(data):
    length = len(data)
    if length <= 1:
        return data
    middle = int(length / 2)
    left = merge_sort(data[:middle])
    right = merge_sort(data[middle:])
    return merge(left, right)


def merge_sort_parallel(data):
    processes = multiprocessing.cpu_count()
    pool = multiprocessing.Pool(processes=processes)
    size = int(math.ceil(float(len(data)) / processes))
    data = [data[i * size:(i + 1) * size] for i in range(processes)]
    data = pool.map(merge_sort, data)

    while len(data) > 1:
        extra = data.pop() if len(data) % 2 == 1 else None
        data = [(data[i], data[i + 1]) for i in range(0, len(data), 2)]
        data = pool.map(merge, data) + ([extra] if extra else [])
    return data[0]


def index_HDF5_rowIDs(filePath):
    file = tb.open_file(filePath, "a")
    chrs = ["X", "Y"]
    chrs.extend(range(1, 23))
    rowIDs = []
    for chr in chrs:
        try:
            node = file.get_node("/chr" + str(chr))
            rownames = node.rownames[:].tolist()
            for idx, rowname in enumerate(rownames):
                rowIDs.append([rowname, idx])
        except tb.exceptions.NoSuchNodeError:
            pass
        except Exception as e:
            print(e)
            pass
    file.close()

    # for sort in merge_sort, merge_sort_parallel:
    start = time.time()
    rowIDs_index = merge_sort_parallel(rowIDs)
    end = time.time() - start
    # print(sort.__name__, end)
    return rowIDs_index


def binarySearch(sortedIDs, start, end, id):

    if end >= start:
        mid = start + int((end - start) / 2)
        if sortedIDs[mid][0] == id:
            return mid
        elif sortedIDs[mid][0] > id:
            return binarySearch(sortedIDs, start, mid - 1, id)
        else:
            return binarySearch(sortedIDs, mid + 1, end, id)
    else:
        return -1
