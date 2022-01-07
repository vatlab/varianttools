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
##

from heapq import heappush, heappop, heappushpop
from multiprocessing import Process, Pipe
from .utils import openFile, env


def TextReader(processor,
               input,
               varIdx,
               getNew,
               jobs,
               encoding,
               header,
               quiet=False):
    '''
    input: input file
    varIdx: variant index, if specified, only matching variants will be returned
        used by updater, or only new variants will be returned by an variant reader
    getNew: if getNew is true, only return variants that are NOT in varIdx. Otherwise,
        return variants that are in varIdx.
    jobs: number of jobs
    encoding: file encoding
    '''
    if jobs == 0:
        return EmbeddedTextReader(processor, input, varIdx, getNew, encoding,
                                  header, quiet)
    elif jobs == 1:
        return StandaloneTextReader(processor, input, varIdx, getNew, encoding,
                                    header, quiet)
    else:
        return MultiTextReader(processor, input, varIdx, getNew, jobs, encoding,
                               header, quiet)


class EmbeddedTextReader:
    #
    # This reader read the file from the main process. No separate process is spawned.
    #
    def __init__(self,
                 processor,
                 input,
                 varIdx,
                 getNew,
                 encoding,
                 header,
                 quiet=False):
        self.num_records = 0
        self.unprocessable_lines = 0
        self.processor = processor
        self.input = input
        self.varIdx = varIdx
        self.getNew = getNew
        self.encoding = encoding
        self.header = header
        self.quiet = quiet

    def records(self):
        first = True
        in_header = True
        line_no = 0
        with openFile(self.input) as input_file:
            for line in input_file:
                line_no += 1
                line = line.decode(self.encoding)
                try:
                    if in_header:
                        # default behavior
                        if self.header is None:
                            if not line.startswith('#'):
                                in_header = False
                        # skip a few lines
                        elif type(self.header) == int:
                            if line_no > self.header:
                                in_header = False
                        # a pattern
                        else:
                            if not self.header.match(line):
                                in_header = False
                        if in_header:
                            continue
                    for bins, rec in self.processor.process(line):
                        if first:
                            self.columnRange = self.processor.columnRange
                            first = False
                        self.num_records += 1
                        if self.varIdx is not None:
                            var_key = (rec[0], rec[2], rec[3])
                            if self.getNew:
                                # only need new variant, so continue if variant in varIdx.
                                if var_key in self.varIdx and rec[
                                        1] in self.varIdx[var_key]:
                                    continue
                            else:
                                # only need existing variant, continue if variant not in varIdx
                                if var_key not in self.varIdx or rec[
                                        1] not in self.varIdx[var_key]:
                                    continue
                        yield (line_no, bins, rec)
                except Exception as e:
                    if not self.quiet:
                        env.logger.debug(
                            'Failed to process line "{}...": {}'.format(
                                line[:20].strip(), e))
                    self.unprocessable_lines += 1


class ReaderWorker(Process):
    '''
    This class starts a process and use passed LineProcessor
    to process input line. If multiple works are started,
    they read lines while skipping lines (e.g. 1, 3, 5, 7, ...)
    '''

    def __init__(self,
                 processor,
                 input,
                 varIdx,
                 getNew,
                 output,
                 step,
                 index,
                 encoding,
                 header,
                 quiet=False):
        '''
        processor:  line processor
        input:      input filename
        output:     a pipe to write output
        varIdx:     a dictionary of variants, used to only return matching lines if specified
        step:       step between processing lines because multiple workers might read the same file
        index:      index of this worker in the worker group
        encoding:   file encoding
        '''
        Process.__init__(self, name='FileReader')
        self.daemon = True
        self.processor = processor
        self.input = input
        self.output = output
        self.step = step
        self.varIdx = varIdx
        self.getNew = getNew
        self.index = index
        self.encoding = encoding
        self.header = header
        self.quiet = quiet

    def run(self):
        first = True
        in_header = True
        num_records = 0
        unprocessable_lines = 0
        line_no = 0
        with openFile(self.input) as input_file:
            for line in input_file:
                line_no += 1
                if line_no % self.step != self.index:
                    continue
                line = line.decode(self.encoding)
                try:
                    if in_header:
                        # default behavior
                        if self.header is None:
                            if not line.startswith('#'):
                                in_header = False
                        # skip a few lines
                        elif type(self.header) == int:
                            if line_no > self.header:
                                in_header = False
                        # a pattern
                        else:
                            if not self.header.match(line):
                                in_header = False
                        if in_header:
                            continue
                    for bins, rec in self.processor.process(line):
                        if first:
                            self.output.send(self.processor.columnRange)
                            first = False
                        num_records += 1
                        if self.varIdx is not None:
                            var_key = (rec[0], rec[2], rec[3])
                            if self.getNew:
                                # only need new variant, so continue if variant in varIdx.
                                if var_key in self.varIdx and rec[
                                        1] in self.varIdx[var_key]:
                                    continue
                            else:
                                # only need existing variant, continue if variant not in varIdx
                                if var_key not in self.varIdx or rec[
                                        1] not in self.varIdx[var_key]:
                                    continue
                        self.output.send((line_no, bins, rec))
                except Exception as e:
                    if not self.quiet:
                        env.logger.debug(
                            'Failed to process line "{}...": {}'.format(
                                line[:20].strip(), e))
                    unprocessable_lines += 1
        # if still first (this thread has not read anything), still send the columnRange stuff
        if first:
            self.output.send(self.processor.columnRange)
        # everything is done, stop the pipe
        self.output.send(None)
        # and send the summary statistics
        self.output.send((num_records, unprocessable_lines))
        self.output.close()


class StandaloneTextReader:
    ''' This processor fire up 1 worker to read an input file
    and gather their outputs
    '''

    def __init__(self,
                 processor,
                 input,
                 varIdx,
                 getNew,
                 encoding,
                 header,
                 quiet=False):
        self.num_records = 0
        self.unprocessable_lines = 0
        #
        self.reader, w = Pipe(False)
        self.worker = ReaderWorker(processor, input, varIdx, getNew, w, 1, 0,
                                   encoding, header, quiet)
        self.worker.start()
        # the send value is columnRange
        self.columnRange = self.reader.recv()

    def records(self):
        while True:
            val = self.reader.recv()
            if val is None:
                self.num_records, self.unprocessable_lines = self.reader.recv()
                break
            else:
                yield val
        self.worker.join()


class MultiTextReader:
    '''This processor fire up num workers to read an input file
    and gather their outputs
    '''

    def __init__(self,
                 processor,
                 input,
                 varIdx,
                 getNew,
                 jobs,
                 encoding,
                 header,
                 quiet=False):
        self.readers = []
        self.workers = []
        self.num_records = 0
        self.unprocessable_lines = 0
        for i in range(jobs):
            r, w = Pipe(False)
            p = ReaderWorker(processor, input, varIdx, getNew, w, jobs, i,
                             encoding, header, quiet)
            self.readers.append(r)
            self.workers.append(p)
            p.start()
        # the send value is columnRange
        for reader in self.readers:
            self.columnRange = reader.recv()

    def records(self):
        all_workers = len(self.readers)
        still_working = len(self.readers)
        #
        # we need a heap to keep records read from multiple processes in order
        # we can not really guarantee this if there are large trunks of ignored
        # records but a heap size = 4 * number of readers should work in most cases
        #
        heap = []
        filled = False
        while True:
            for idx, reader in enumerate(self.readers):
                val = reader.recv()
                if val is None:
                    # one thread died ...
                    still_working -= 1
                    nr, sl = reader.recv()
                    self.num_records += nr
                    self.unprocessable_lines += sl
                    self.readers[idx] = None
                elif filled:
                    yield heappushpop(heap, val)
                else:
                    heappush(heap, val)
                    filled = len(heap) == 4 * len(self.readers)
            if still_working < all_workers:
                self.readers = [x for x in self.readers if x is not None]
                all_workers = len(self.readers)
                if all_workers == 0:
                    # return all things in the heap
                    for i in range(len(heap)):
                        yield heappop(heap)
                    break
        for p in self.workers:
            p.join()
