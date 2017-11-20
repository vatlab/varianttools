
# -*- coding: utf-8 -*-
"""
Extract data from VCF files.

This module contains Functions for extracting data from Variant Call Format (VCF) files and loading
into NumPy arrays, NumPy files, HDF5 files or Zarr array stores.

"""
from __future__ import absolute_import, print_function, division
import os
import sys
import re
import array
import time
from heapq import heappush, heappop, heappushpop
from multiprocessing import Process, Pipe, Value, Lock, Manager,Array
from multiprocessing import Queue as mpQueue
import queue
from itertools import repeat
from collections import defaultdict
from .project import Project, fileFMT
from .liftOver import LiftOverTool
from .utils import ProgressBar, lineCount, getMaxUcscBin, delayedAction, \
    openFile, DatabaseEngine, hasCommand, \
    downloadFile, env, RefGenome

import numpy as np
import tables as tb
#import HDF5_storage as storage
import glob
from shutil import copyfile
from .accessor import *

try:
    from variant_tools.cgatools import normalize_variant
except ImportError as e:
    sys.exit('Failed to import module ({})\n'
        'Please verify if you have installed variant tools successfully (using command '
        '"python setup.py install")'.format(e))


# preprocessors
from .preprocessor import *


import gzip
from collections import namedtuple
import warnings
import subprocess


from variant_tools.io_vcf_read import VCFChunkIterator, FileInputStream
# expose some names from cython extension
# from variant_tools.io_vcf_read import (  # noqa: F401
#     ANNTransformer, ANN_AA_LENGTH_FIELD, ANN_AA_POS_FIELD, ANN_ANNOTATION_FIELD,
#     ANN_ANNOTATION_IMPACT_FIELD, ANN_CDNA_LENGTH_FIELD, ANN_CDNA_POS_FIELD, ANN_CDS_LENGTH_FIELD,
#     ANN_CDS_POS_FIELD, ANN_DISTANCE_FIELD, ANN_FEATURE_ID_FIELD, ANN_FEATURE_TYPE_FIELD, ANN_FIELD,
#     ANN_FIELDS, ANN_GENE_ID_FIELD, ANN_GENE_NAME_FIELD, ANN_HGVS_C_FIELD, ANN_HGVS_P_FIELD,
#     ANN_RANK_FIELD, ANN_TRANSCRIPT_BIOTYPE_FIELD
# )


DEFAULT_BUFFER_SIZE = 2**14
DEFAULT_CHUNK_LENGTH = 2**13
DEFAULT_CHUNK_WIDTH = 2**6
DEFAULT_ALT_NUMBER = 3


def _prep_fields_param(fields):
    """Prepare the `fields` parameter, and determine whether or not to store samples."""

    store_samples = False

    if fields is None:
        # add samples by default
        return True, None

    if isinstance(fields, str):
        fields = [fields]
    else:
        fields = list(fields)

    if 'samples' in fields:
        fields.remove('samples')
        store_samples = True
    elif '*' in fields:
        store_samples = True

    return store_samples, fields


def _chunk_iter_progress(it, log, prefix):
    """Wrap a chunk iterator for progress logging."""
    n_variants = 0
    before_all = time.time()
    before_chunk = before_all
    for chunk, chunk_length, chrom, pos in it:
        after_chunk = time.time()
        elapsed_chunk = after_chunk - before_chunk
        elapsed = after_chunk - before_all
        n_variants += chunk_length
        chrom = str(chrom, 'ascii')
        message = (
            '%s %s rows in %.2fs; chunk in %.2fs (%s rows/s)' %
            (prefix, n_variants, elapsed, elapsed_chunk, int(chunk_length // elapsed_chunk))
        )
        if chrom:
            message += '; %s:%s' % (chrom, pos)
        print(message, file=log)
        log.flush()
        yield chunk, chunk_length, chrom, pos
        before_chunk = after_chunk
    after_all = time.time()
    elapsed = after_all - before_all
    print('%s all done (%s rows/s)' %
          (prefix, int(n_variants // elapsed)), file=log)
    log.flush()


def _chunk_iter_transform(it, transformers):
    for chunk, chunk_length, chrom, pos in it:
        for transformer in transformers:
            transformer.transform_chunk(chunk)
        yield chunk, chunk_length, chrom, pos


_doc_param_input = \
    """Path to VCF file on the local file system. May be uncompressed or gzip-compatible
        compressed file. May also be a file-like object (e.g., `io.BytesIO`)."""

_doc_param_fields = \
    """Fields to extract data for. Should be a list of strings, e.g., ``['variants/CHROM',
        'variants/POS', 'variants/DP', 'calldata/GT']``. If you are feeling lazy, you can drop
        the 'variants/' and 'calldata/' prefixes, in which case the fields will be matched against
        fields declared in the VCF header, with variants taking priority over calldata if a field
        with the same ID exists both in INFO and FORMAT headers. I.e., ``['CHROM', 'POS', 'DP',
        'GT']`` will work, although watch out for fields like 'DP' which can be both
        INFO and FORMAT. For convenience, some special string values are also recognized. To
        extract all fields, provide just the string ``'*'``. To extract all variants fields
        (including all INFO fields) provide ``'variants/*'``. To extract all calldata fields (i.e.,
        defined in FORMAT headers) provide ``'calldata/*'``."""

_doc_param_types = \
    """Overide data types. Should be a dictionary mapping field names to NumPy data types.
        E.g., providing the dictionary ``{'variants/DP': 'i8', 'calldata/GQ': 'i2'}`` will mean
        the 'variants/DP' field is stored in a 64-bit integer array, and the 'calldata/GQ' field
        is stored in a 16-bit integer array."""

_doc_param_numbers = \
    """Override the expected number of values. Should be a dictionary mapping field names to
        integers. E.g., providing the dictionary ``{'variants/ALT': 5, 'variants/AC': 5,
        'calldata/HQ': 2}`` will mean that, for each variant, 5 values are stored for the
        'variants/ALT' field, 5 values are stored for the 'variants/AC' field, and for each
        sample, 2 values are stored for the 'calldata/HQ' field."""

_doc_param_alt_number = \
    """Assume this number of alternate alleles and set expected number of values accordingly for
        any field declared with number 'A' or 'R' in the VCF meta-information."""

_doc_param_fills = \
    """Override the fill value used for empty values. Should be a dictionary mapping field names
        to fill values."""

_doc_param_region = \
    """Genomic region to extract variants for. If provided, should be a tabix-style region string,
        which can be either just a chromosome name (e.g., '2L'), or a chromosome name followed by
        1-based beginning and end coordinates (e.g., '2L:100000-200000'). Note that only variants
        whose start position (POS) is within the requested range will be included. This is slightly
        different from the default tabix behaviour, where a variant (e.g., deletion) may be included
        if its position (POS) occurs before the requested region but its reference allele overlaps
        the region - such a variant will *not* be included in the data returned by this function."""

_doc_param_tabix = \
    """Name or path to tabix executable. Only required if `region` is given. Setting `tabix` to
        `None` will cause a fall-back to scanning through the VCF file from the beginning, which
        may be much slower than tabix but the only option if tabix is not available on your system
        and/or the VCF file has not been tabix-indexed."""

_doc_param_samples = \
    """Selection of samples to extract calldata for. If provided, should be a list of strings
        giving sample identifiers. May also be a list of integers giving indices of selected
        samples."""

_doc_param_transformers = \
    """Transformers for post-processing data. If provided, should be a list of Transformer
        objects, each of which must implement a "transform()" method that accepts a dict
        containing the chunk of data to be transformed. See also the :class:`ANNTransformer`
        class which implements post-processing of data from SNPEFF."""

_doc_param_buffer_size = \
    """Size in bytes of the I/O buffer used when reading data from the underlying file or tabix
        stream."""

_doc_param_chunk_length = \
    """Length (number of variants) of chunks in which data are processed."""

_doc_param_log = \
    """A file-like object (e.g., `sys.stderr`) to print progress information."""





def _hdf5_setup_datasets(chunk, root, chunk_length, chunk_width, compression, compression_opts,
                         shuffle, overwrite, headers, vlen):
    import h5py

    # handle no input
    if chunk is None:
        raise RuntimeError('input file has no data?')

    # setup datasets
    # keys = sorted(chunk.keys())
    keys =["variants/ID","calldata/GT"]
    for k in keys:

        # obtain initial data
        data = chunk[k]

        # determine chunk shape
        if data.ndim == 1:
            chunk_shape = (chunk_length,)
        else:
            chunk_shape = (chunk_length, min(chunk_width, data.shape[1])) + data.shape[2:]

        # create dataset
        group, name = k.split('/')
        if name in root[group]:
            if overwrite:
                del root[group][name]
            else:
                raise ValueError('dataset exists at path %r; use overwrite=True to replace' % k)

        shape = (0,) + data.shape[1:]
        maxshape = (None,) + data.shape[1:]
        if data.dtype.kind == 'O':
            if vlen:
                dt = h5py.special_dtype(vlen=str)
            else:
                data = data.astype('S')
                dt = data.dtype
        else:
            dt = data.dtype
        ds = root[group].create_dataset(
            name, shape=shape, maxshape=maxshape, chunks=chunk_shape, dtype=dt,
            compression=compression, compression_opts=compression_opts, shuffle=shuffle
        )

        # copy metadata from VCF headers
        # meta = None
        # if group == 'variants' and name in headers.infos:
        #     meta = headers.infos[name]
        # elif group == 'calldata' and name in headers.formats:
        #     meta = headers.formats[name]
        # if meta is not None:
        #     ds.attrs['ID'] = meta['ID']
        #     ds.attrs['Number'] = meta['Number']
        #     ds.attrs['Type'] = meta['Type']
        #     ds.attrs['Description'] = meta['Description']

    return keys





_doc_param_chunk_width = \
    """Width (number of samples) to use when storing chunks in output."""



def iter_vcf_chunks(input,
                    fields=None,
                    types=None,
                    numbers=None,
                    alt_number=DEFAULT_ALT_NUMBER,
                    fills=None,
                    region=None,
                    tabix='tabix',
                    samples=None,
                    transformers=None,
                    buffer_size=DEFAULT_BUFFER_SIZE,
                    chunk_length=DEFAULT_CHUNK_LENGTH):
    """Iterate over chunks of data from a VCF file as NumPy arrays.

    Parameters
    ----------
    input : string
        {input}
    fields : list of strings, optional
        {fields}
    types : dict, optional
        {types}
    numbers : dict, optional
        {numbers}
    alt_number : int, optional
        {alt_number}
    fills : dict, optional
        {fills}
    region : string, optional
        {region}
    tabix : string, optional
        {tabix}
    samples : list of strings
        {samples}
    transformers : list of transformer objects, optional
        {transformers}
    buffer_size : int, optional
        {buffer_size}
    chunk_length : int, optional
        {chunk_length}

    Returns
    -------
    fields : list of strings
        Normalised names of fields that will be extracted.
    samples : ndarray
        Samples for which data will be extracted.
    headers : VCFHeaders
        Tuple of metadata extracted from VCF headers.
    it : iterator
        Chunk iterator.

    """

    # setup commmon keyword args
    kwds = dict(fields=fields, types=types, numbers=numbers, alt_number=alt_number,
                chunk_length=chunk_length, fills=fills, samples=samples)

    # obtain a file-like object
    close = False
    if isinstance(input, str) and input.endswith('gz'):

        if region and tabix and os.name != 'nt':

            try:
                # try tabix
                p = subprocess.Popen([tabix, '-h', input, region],
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.STDOUT,
                                     bufsize=0)

                # check if tabix exited early, look for tabix error
                time.sleep(.5)
                poll = p.poll()
                if poll is not None and poll > 0:
                    err = p.stdout.read()
                    err = str(err, 'ascii')
                    p.stdout.close()
                    raise RuntimeError(err.strip())
                fileobj = p.stdout
                close = True
                # N.B., still pass the region parameter through so we get strictly only
                # variants that start within the requested region. See also
                # https://github.com/alimanfoo/vcfnp/issues/54

            except FileNotFoundError:
                # no tabix, fall back to scanning
                warnings.warn('tabix not found, falling back to scanning to region')
                fileobj = gzip.open(input, mode='rb')
                close = True

            except Exception as e:
                warnings.warn('error occurred attempting tabix (%s); falling back to '
                              'scanning to region' % e)
                fileobj = gzip.open(input, mode='rb')
                close = True

        else:
            fileobj = gzip.open(input, mode='rb')
            close = True

    elif isinstance(input, str):
        # assume no compression
        fileobj = open(input, mode='rb', buffering=0)
        close = True

    elif hasattr(input, 'readinto'):
        fileobj = input

    else:
        raise ValueError('path must be string or file-like, found %r' % input)

    # setup input stream
    stream = FileInputStream(fileobj, buffer_size=buffer_size, close=close)

    # deal with region
    kwds['region'] = region

    # setup iterator
    fields, samples, headers, it = _iter_vcf_stream(stream, **kwds)

    # setup transformers
    if transformers is not None:
        # API flexibility
        if not isinstance(transformers, (list, tuple)):
            transformers = [transformers]
        for trans in transformers:
            fields = trans.transform_fields(fields)
        it = _chunk_iter_transform(it, transformers)

    return fields, samples, headers, it


iter_vcf_chunks.__doc__ = iter_vcf_chunks.__doc__.format(
    input=_doc_param_input,
    fields=_doc_param_fields,
    types=_doc_param_types,
    numbers=_doc_param_numbers,
    alt_number=_doc_param_alt_number,
    fills=_doc_param_fills,
    region=_doc_param_region,
    tabix=_doc_param_tabix,
    samples=_doc_param_samples,
    transformers=_doc_param_transformers,
    buffer_size=_doc_param_buffer_size,
    chunk_length=_doc_param_chunk_length,
    log=_doc_param_log,
)


FIXED_VARIANTS_FIELDS = (
    'CHROM',
    'POS',
    'ID',
    'REF',
    'ALT',
    'QUAL',
)


def _normalize_field_prefix(field, headers):

    # already contains prefix?
    if field.startswith('variants/') or field.startswith('calldata/'):
        return field

    # try to find in fixed fields
    elif field in FIXED_VARIANTS_FIELDS:
        return 'variants/' + field

    # try to find in FILTER
    elif field.startswith('FILTER_'):
        return 'variants/' + field

    # try to find in FILTER
    elif field in headers.filters:
        return 'variants/FILTER_' + field

    # try to find in INFO
    elif field in headers.infos:
        return 'variants/' + field

    # try to find in FORMAT
    elif field in headers.formats:
        return 'calldata/' + field

    else:
        # assume anything else in variants, even if header not found
        return 'variants/' + field


def _check_field(field, headers):

    # assume field is already normalized for prefix
    group, name = field.split('/')

    if group == 'variants':

        if name in FIXED_VARIANTS_FIELDS:
            pass

        elif name in ['numalt', 'svlen', 'is_snp']:
            # computed fields
            pass

        elif name.startswith('FILTER_'):
            filter_name = name[7:]
            if filter_name in headers.filters:
                pass
            else:
                warnings.warn('%r FILTER header not found' % filter_name)

        elif name in headers.infos:
            pass

        else:
            warnings.warn('%r INFO header not found' % name)

    elif group == 'calldata':

        if name in headers.formats:
            pass

        else:
            warnings.warn('%r FORMAT header not found' % name)

    else:
        # should never be reached
        raise ValueError('invalid field specification: %r' % field)


def _add_all_fields(fields, headers, samples):
    _add_all_variants_fields(fields, headers)
    if len(samples) > 0:
        _add_all_calldata_fields(fields, headers)


def _add_all_variants_fields(fields, headers):
    _add_all_fixed_variants_fields(fields)
    _add_all_info_fields(fields, headers)
    _add_all_filter_fields(fields, headers)
    # add in computed fields
    for f in 'variants/numalt', 'variants/svlen', 'variants/is_snp':
        if f not in fields:
            fields.append(f)


def _add_all_fixed_variants_fields(fields):
    for k in FIXED_VARIANTS_FIELDS:
        f = 'variants/' + k
        if f not in fields:
            fields.append(f)


def _add_all_info_fields(fields, headers):
    for k in headers.infos:
        f = 'variants/' + k
        if f not in fields:
            fields.append(f)


def _add_all_filter_fields(fields, headers):
    fields.append('variants/FILTER_PASS')
    for k in headers.filters:
        f = 'variants/FILTER_' + k
        if f not in fields:
            fields.append(f)


def _add_all_calldata_fields(fields, headers):
    # only add calldata fields if there are samples
    if headers.samples:
        for k in headers.formats:
            f = 'calldata/' + k
            if f not in fields:
                fields.append(f)


def _normalize_fields(fields, headers, samples):

    # setup normalized fields
    normed_fields = list()

    # special case, single field specification
    if isinstance(fields, str):
        fields = [fields]

    for f in fields:

        # special cases: be lenient about how to specify

        if f in ['*', 'kitchen sink']:
            _add_all_fields(normed_fields, headers, samples)

        elif f in ['variants', 'variants*', 'variants/*']:
            _add_all_variants_fields(normed_fields, headers)

        elif f in ['calldata', 'calldata*', 'calldata/*'] and len(samples) > 0:
            _add_all_calldata_fields(normed_fields, headers)

        elif f in ['INFO', 'INFO*', 'INFO/*', 'variants/INFO', 'variants/INFO*', 'variants/INFO/*']:
            _add_all_info_fields(normed_fields, headers)

        elif f in ['FILTER', 'FILTER*', 'FILTER/*', 'FILTER_*', 'variants/FILTER',
                   'variants/FILTER*', 'variants/FILTER/*', 'variants/FILTER_*']:
            _add_all_filter_fields(normed_fields, headers)

        # exact field specification

        else:

            # normalize field specification
            f = _normalize_field_prefix(f, headers)
            _check_field(f, headers)
            if f.startswith('calldata/') and len(samples) == 0:
                # only add calldata fields if there are samples
                pass
            elif f not in normed_fields:
                normed_fields.append(f)

    return normed_fields


default_integer_dtype = 'i4'
default_float_dtype = 'f4'
default_string_dtype = 'object'


def _normalize_type(t):
    if t == 'Integer':
        return np.dtype(default_integer_dtype)
    elif t == 'Float':
        return np.dtype(default_float_dtype)
    elif t == 'String':
        return np.dtype(default_string_dtype)
    elif t == 'Character':
        return np.dtype('S1')
    elif t == 'Flag':
        return np.dtype(bool)
    elif isinstance(t, str) and t.startswith('genotype/'):
        # custom genotype dtype
        return t
    elif isinstance(t, str) and t.startswith('genotype_ac/'):
        # custom genotype allele counts dtype
        return t
    else:
        return np.dtype(t)


default_types = {
    'variants/CHROM': 'object',
    'variants/POS': 'i4',
    'variants/ID': 'object',
    'variants/REF': 'object',
    'variants/ALT': 'object',
    'variants/QUAL': 'f4',
    'variants/DP': 'i4',
    'variants/AN': 'i4',
    'variants/AC': 'i4',
    'variants/AF': 'f4',
    'variants/MQ': 'f4',
    'variants/ANN': 'object',
    'calldata/GT': 'genotype/i1',
    'calldata/GQ': 'f4',
    'calldata/HQ': 'i1',
    'calldata/DP': 'i2',
    'calldata/AD': 'i2',
    'calldata/MQ0': 'i2',
    'calldata/MQ': 'f2',
}


def _normalize_types(types, fields, headers):

    # normalize user-provided types
    if types is None:
        types = dict()
    types = {_normalize_field_prefix(f, headers): _normalize_type(t)
             for f, t in types.items()}

    # setup output
    normed_types = dict()

    for f in fields:

        group, name = f.split('/')

        if f in types:
            normed_types[f] = types[f]

        elif f in default_types:
            normed_types[f] = _normalize_type(default_types[f])

        elif group == 'variants':

            if name in ['numalt', 'svlen', 'is_snp']:
                # computed fields, special case
                continue

            elif name.startswith('FILTER_'):
                normed_types[f] = np.dtype(bool)

            elif name in headers.infos:
                normed_types[f] = _normalize_type(headers.infos[name]['Type'])

            else:
                # fall back to string
                normed_types[f] = _normalize_type('String')
                warnings.warn('no type for field %r, assuming %s' %
                              (f, normed_types[f]))

        elif group == 'calldata':

            if name in headers.formats:
                normed_types[f] = _normalize_type(headers.formats[name]['Type'])

            else:
                # fall back to string
                normed_types[f] = _normalize_type('String')
                warnings.warn('no type for field %r, assuming %s' %
                              (f, normed_types[f]))

        else:
            raise RuntimeError('unpected field: %r' % f)

    return normed_types


default_numbers = {
    'variants/CHROM': 1,
    'variants/POS': 1,
    'variants/ID': 1,
    'variants/REF': 1,
    'variants/ALT': 'A',
    'variants/QUAL': 1,
    'variants/DP': 1,
    'variants/AN': 1,
    'variants/AC': 'A',
    'variants/AF': 'A',
    'variants/MQ': 1,
    'variants/ANN': 1,
    'calldata/DP': 1,
    'calldata/GT': 1,
    'calldata/GQ': 1,
    'calldata/HQ': 2,
    'calldata/AD': 'R',
    'calldata/MQ0': 1,
    'calldata/MQ': 1,
}


def _normalize_number(field, n, alt_number):
    if n == '.':
        return 1
    elif n == 'A':
        return alt_number
    elif n == 'R':
        return alt_number + 1
    elif n == 'G':
        return 3
    else:
        try:
            return int(n)
        except ValueError:
            warnings.warn('error parsing %r as number for field %r' % (n, field))
        return 1


def _normalize_numbers(numbers, fields, headers, alt_number):

    # normalize field prefixes
    if numbers is None:
        numbers = dict()
    numbers = {_normalize_field_prefix(f, headers): n for f, n in numbers.items()}

    # setup output
    normed_numbers = dict()

    for f in fields:

        group, name = f.split('/')

        if f in numbers:
            normed_numbers[f] = _normalize_number(f, numbers[f], alt_number)

        elif f in default_numbers:
            normed_numbers[f] = _normalize_number(f, default_numbers[f], alt_number)

        elif group == 'variants':

            if name in ['numalt', 'svlen', 'is_snp']:
                # computed fields, special case (for svlen, number depends on ALT)
                continue

            elif name.startswith('FILTER_'):
                normed_numbers[f] = 0

            elif name in headers.infos:
                normed_numbers[f] = _normalize_number(f, headers.infos[name]['Number'], alt_number)

            else:
                # fall back to 1
                normed_numbers[f] = 1
                warnings.warn('no number for field %r, assuming 1' % f)

        elif group == 'calldata':

            if name in headers.formats:
                normed_numbers[f] = _normalize_number(f, headers.formats[name]['Number'],
                                                      alt_number)

            else:
                # fall back to 1
                normed_numbers[f] = 1
                warnings.warn('no number for field %r, assuming 1' % f)

        else:
            raise RuntimeError('unexpected field: %r' % f)

    return normed_numbers


def _normalize_fills(fills, fields, headers):

    if fills is None:
        fills = dict()
    fills = {_normalize_field_prefix(f, headers): v
             for f, v in fills.items()}

    # setup output
    normed_fills = dict()

    for f in fields:

        if f in fills:
            normed_fills[f] = fills[f]

    return normed_fills


def _normalize_samples(samples, headers, types):
    loc_samples = np.zeros(len(headers.samples), dtype='u1')

    if samples is None:
        normed_samples = list(headers.samples)
        loc_samples.fill(1)

    else:
        samples = set(samples)
        normed_samples = []
        for i, s in enumerate(headers.samples):
            if i in samples:
                normed_samples.append(s)
                samples.remove(i)
                loc_samples[i] = 1
            elif s in samples:
                normed_samples.append(s)
                samples.remove(s)
                loc_samples[i] = 1
        if len(samples) > 0:
            warnings.warn('some samples not found, will be ignored: ' +
                          ', '.join(map(repr, sorted(samples))))

    t = default_string_dtype
    if types is not None:
        t = types.get('samples', t)
    normed_samples = np.array(normed_samples, dtype=t)

    return normed_samples, loc_samples


def _iter_vcf_stream(stream, fields, types, numbers, alt_number, chunk_length, fills, region,
                     samples):

    # read VCF headers
    headers = _read_vcf_headers(stream)

    # setup samples
    samples, loc_samples = _normalize_samples(samples=samples, headers=headers, types=types)

    # setup fields to read
    if fields is None:

        # choose default fields
        fields = list()
        _add_all_fixed_variants_fields(fields)
        fields.append('variants/FILTER_PASS')
        if len(samples) > 0 and 'GT' in headers.formats:
            fields.append('calldata/GT')

    else:
        fields = _normalize_fields(fields=fields, headers=headers, samples=samples)

    # setup data types
    types = _normalize_types(types=types, fields=fields, headers=headers)

    # setup numbers (a.k.a., arity)
    numbers = _normalize_numbers(numbers=numbers, fields=fields, headers=headers,
                                 alt_number=alt_number)

    # setup fills
    fills = _normalize_fills(fills=fills, fields=fields, headers=headers)

    # setup chunks iterator
    chunks = VCFChunkIterator(
        stream, chunk_length=chunk_length, headers=headers, fields=fields, types=types,
        numbers=numbers, fills=fills, region=region, loc_samples=loc_samples
    )

    return fields, samples, headers, chunks


# pre-compile some regular expressions
_re_filter_header = \
    re.compile('##FILTER=<ID=([^,]+),Description="([^"]+)">')
_re_info_header = \
    re.compile('##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)">')
_re_format_header = \
    re.compile('##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)">')


VCFHeaders = namedtuple('VCFHeaders', ['headers', 'filters', 'infos', 'formats', 'samples'])


def _read_vcf_headers(stream):

    # setup
    headers = []
    samples = None
    filters = dict()
    infos = dict()
    formats = dict()

    # read first header line
    header = stream.readline()
    header = str(header, 'ascii')

    while header and header[0] == '#':

        headers.append(header)

        if header.startswith('##FILTER'):

            match = _re_filter_header.match(header)
            if match is None:
                warnings.warn('invalid FILTER header: %r' % header)
            else:
                k, d = match.groups()
                if k in filters:
                    warnings.warn('multiple FILTER headers for %r' % k)
                filters[k] = {'ID': k, 'Description': d}

        elif header.startswith('##INFO'):

            match = _re_info_header.match(header)
            if match is None:
                warnings.warn('invalid INFO header: %r' % header)
            else:
                k, n, t, d = match.groups()
                if k in infos:
                    warnings.warn('multiple INFO headers for %r' % k)
                infos[k] = {'ID': k, 'Number': n, 'Type': t, 'Description': d}

        elif header.startswith('##FORMAT'):

            match = _re_format_header.match(header)
            if match is None:
                warnings.warn('invalid FORMAT header: %r' % header)
            else:
                k, n, t, d = match.groups()
                if k in formats:
                    warnings.warn('multiple FORMAT headers for %r' % k)
                formats[k] = {'ID': k, 'Number': n, 'Type': t, 'Description': d}

        elif header.startswith('#CHROM'):

            # parse out samples
            samples = header.strip().split('\t')[9:]
            break

        # read next header line
        header = stream.readline()
        header = str(header, 'ascii')

    # check if we saw the mandatory header line or not
    if samples is None:
        # can't warn about this, it's fatal
        raise RuntimeError('VCF file is missing mandatory header line ("#CHROM...")')

    return VCFHeaders(headers, filters, infos, formats, samples)





class HDF5GenotypeImportWorker(Process):
    '''This class starts a process, import genotype to a temporary HDF5 file.
        Args

            processor: a processor to parse the row 
            readQueue: a Queue of vcf rows
            variantIndex: a dictionary that returns ID for each variant.
            start_sample: the sample id of sample in the first column
            end_sample: the sample id of sample in the last column
            sample_ids: a list of sample IDS
            variant_count: count of variants
            proc_index: the index of process
            geno_info: genotype info other than GT
            dbLocation: the HDF5 file name 

    '''
    def __init__(self, chunk,variantIndex, start_sample,end_sample,sample_ids,variant_count, 
        proc_index,geno_info,dbLocation,build,fields):

        Process.__init__(self, name='GenotypeImporter')
        # self.daemon=True
        self.chunk=chunk
        self.variantIndex = variantIndex
        self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample=start_sample
        self.end_sample=end_sample
        self.geno_info=geno_info
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]
        self.sample_ids=sample_ids
        self.colnames=[sample_ids[i] for i in range(start_sample,end_sample)]
        self.genoCount=0
        self.dbLocation=dbLocation
        self.build=build
        self.fields=fields

        self.info={}
        self.rowData=[]
        if len(self.geno_info)>0:
            for info in self.geno_info:
                #indptr,indices,data,shape,rownames
                self.info[info.name]=[[],[],[],[],[]]
        

    # def writeIntoFile(self,chr):

    #     self.variant_count.value = self.variant_count.value + len(self.indptr)
    #     shape=(self.variant_count.value,self.end_sample-self.start_sample)

    #     storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
    #     # make a HMatrix object which is a matrix with rownames and colnames
    #     hmatrix=HMatrix(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames)
    #     # write GT into file
    #     storageEngine.store(hmatrix,chr)
    #     # write geno info into HDF5 if exists
    #     if len(self.geno_info)>0:
    #             for key,value in list(self.info.items()):
    #                 hmatrix=HMatrix(value[2],value[1],value[0],shape,value[3],self.colnames)
    #                 storageEngine.store(hmatrix,chr,key) 
        
    #     # clean up
    #     self.indptr=[]
    #     self.indices=[]
    #     self.data=[]
    #     self.rownames=[]

    #     if len(self.geno_info)>0:
    #         for info in self.geno_info:
    #             #indptr,indices,data,shape,rownames
    #             self.info[info.name]=[[],[],[],[],[]]

    def getGT(self,variant_id,GT,altIndex):
        for idx in range(self.start_sample,self.end_sample):
            if GT[idx] is not None:
                if altIndex==0:
                    if GT[idx]!=0:
                        if GT[idx]!=3 and GT[idx]!=4:
                            self.genoCount=self.genoCount+1
                            self.indices.append(idx-self.start_sample)
                            self.data.append(GT[idx])
                        else:
                            self.genoCount=self.genoCount+1
                            self.indices.append(idx-self.start_sample)
                            self.data.append(np.nan)  
                elif altIndex==1:
                        self.genoCount=self.genoCount+1
                        self.indices.append(idx-self.start_sample)   
                        if GT[idx]==3:
                            self.data.append(1)
                        elif GT[idx]==4:
                            self.data.append(2)
                        else:
                            self.data.append(np.nan)
            else:
                self.genoCount=self.genoCount+1
                self.indices.append(idx-self.start_sample)
                self.data.append(np.nan)
        self.indptr.append(self.genoCount)
        self.rownames.append(variant_id)

    #not used right now
    def getInfo(self,variant_id,infoDict,altIndex):
        genoDict={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                genoDict[info.name]=0
        for name,info in infoDict.items():
            if name=="GQ_geno" or name=="DP_geno":
                for idx in range(self.start_sample,self.end_sample):            
                    if info[idx] is not None:
                        if altIndex==0:
                            if info[idx]!='0':
                                self.info[name][2].append(info[idx])
                                genoDict[name]=genoDict[name]+1
                                self.info[name][1].append(idx-self.start_sample)
                        elif altIndex==1:   
                            if info[idx]!="0":
                                self.info[name][2].append(info[idx])
                                genoDict[name]=genoDict[name]+1
                                self.info[name][1].append(idx-self.start_sample)   
                    else:
                        genoDict[name]=genoDict[name]+1
                        self.info[name][1].append(idx-self.start_sample)   
                        self.info[name][2].append(np.nan)

                self.info[name][0].append(genoDict[name])
                self.info[name][3].append(variant_id)
            else:
                for idx in range(self.start_sample,self.end_sample):
                    infoString=",".join([str(x) for x in info[idx]])
                    self.info[name][2].append(infoString)
                    genoDict[name]=genoDict[name]+1
                    self.info[name][1].append(idx-self.start_sample)            
                self.info[name][0].append(genoDict[name])
                self.info[name][3].append(variant_id)

    def getInfoTable(self,variant_id,infoDict,altIndex):
        for idx in range(self.start_sample,self.end_sample):
            rowData=[variant_id,idx]
            samplePos=idx-self.start_sample
            for name in infoDict.keys():
                if name=="DP_geno" or name=="GQ_geno":
                    rowData.append(infoDict[name][samplePos])
                elif name=="AD_geno" and altIndex==0:
                    rowData.append(infoDict[name][samplePos][0])
                    rowData.append(infoDict[name][samplePos][1])
                elif name=="AD_geno" and altIndex==1:
                    rowData.append(infoDict[name][samplePos][0])
                    rowData.append(infoDict[name][samplePos][2])
                elif name=="PL_geno" and altIndex==0:
                    rowData.append(infoDict[name][samplePos][0])
                    rowData.append(infoDict[name][samplePos][1])
                    rowData.append(infoDict[name][samplePos][2])
                elif name=="PL_geno" and altIndex==1:
                    rowData.append(infoDict[name][samplePos][3])
                    rowData.append(infoDict[name][samplePos][4])
                    rowData.append(infoDict[name][samplePos][5])
            self.rowData.append(rowData)



    # # Used for iterator over set of samples
    # def getGT(self,variant_id,GT,altIndex):
    #     for idx in range(self.end_sample-self.start_sample):
    #         if GT[idx] is not None:
    #             if altIndex==0:
    #                 if GT[idx]!=0:
    #                     if GT[idx]!=3 and GT[idx]!=4:
    #                         self.genoCount=self.genoCount+1
    #                         self.indices.append(idx)
    #                         self.data.append(GT[idx])
    #                     else:
    #                         self.genoCount=self.genoCount+1
    #                         self.indices.append(idx)
    #                         self.data.append(np.nan)  
    #             elif altIndex==1:
    #                     self.genoCount=self.genoCount+1
    #                     self.indices.append(idx)   
    #                     if GT[idx]==3:
    #                         self.data.append(1)
    #                     elif GT[idx]==4:
    #                         self.data.append(2)
    #                     else:
    #                         self.data.append(np.nan)
    #         else:
    #             self.genoCount=self.genoCount+1
    #             self.indices.append(idx)
    #             self.data.append(np.nan)
    #     self.indptr.append(self.genoCount)
    #     self.rownames.append(variant_id)


    def writeIntoHDF(self,chr):
        storageEngine=Engine_Storage.choose_storage_engine(self.dbLocation)
        # if self.proc_index==1:
        #     print(self.proc_index,tolistTime,getGTtime,getInfotime)
        shape=(len(self.indptr),len(self.colnames))
        # make a HMatrix object which is a matrix with rownames and colnames
        hmatrix=HMatrix(self.data,self.indices,self.indptr,shape,self.rownames,self.colnames)
        # write GT into file
        # starttime=time.time()
        storageEngine.store(hmatrix,chr,"GT")
        # if self.proc_index==1:
        #     print(self.proc_index,"GTWrite",time.time()-starttime)
        # starttime=time.time()
        # if len(self.geno_info)>0:
        #     for key,value in list(self.info.items()):
        #         hmatrix=HMatrix(value[2],value[1],value[0],shape,value[3],self.colnames)
        #         storageEngine.store(hmatrix,chr,key) 
        #     for info in self.geno_info:
        #         #indptr,indices,data,shape,rownames
        #         self.info[info.name]=[[],[],[],[],[]]
        tableName="genoInfo"
        if len(self.geno_info)>0:
            # print(len(self.rowData),self.rowData[0])
            storageEngine.store_table(self.rowData,tableName,chr)
        # if self.proc_index==1:
        #     print(self.proc_index,"infoWrite",time.time()-starttime)
        storageEngine.close()
        self.indptr=[]
        self.indices=[]
        self.data=[]
        self.rownames=[]
        self.rowData=[]




    def processFile(self):
        prev_chr=self.chunk["variants/CHROM"][0].replace("chr","")
        tolistTime=0
        getGTtime=0
        getInfotime=0
        for i in range(len(self.chunk["variants/ID"])):
            infoDict={}
            chr=self.chunk["variants/CHROM"][i].replace("chr","")
            if chr!=prev_chr:
                self.writeIntoHDF(prev_chr)
                prev_chr=chr
                
            # if chr.startswith("chr"):
            #     chr=chr.replace("chr","")
            ref=self.chunk["variants/REF"][i]
            pos=self.chunk["variants/POS"][i]
            GT=self.chunk["calldata/GT"][i].tolist()
            # starttime=time.time()
            # if len(self.geno_info)>0:
            #     if "calldata/DP" in self.chunk:
            #         DP_geno=self.chunk["calldata/DP"][i].tolist()
            #         infoDict["DP_geno"]=DP_geno
            #     if "calldata/GQ" in self.chunk:
            #         GQ_geno=self.chunk["calldata/GQ"][i].tolist()
            #         infoDict["GQ_geno"]=GQ_geno
            #     if "calldata/AD" in self.chunk:
            #         AD_geno=self.chunk["calldata/AD"][i].tolist()
            #         infoDict["AD_geno"]=AD_geno
            #     if "calldata/PL" in self.chunk:
            #         PL_geno=self.chunk["calldata/PL"][i].tolist()
            #         infoDict["PL_geno"]=PL_geno
            # tolistTime+=time.time()-starttime
            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt=self.chunk["variants/ALT"][i][altIndex]
                if alt!="":
                    if tuple((chr, ref, alt)) in self.variantIndex:
                        variant_id  = self.variantIndex[tuple((chr, ref, alt))][pos][0]
                        # starttime=time.time()
                        self.getGT(variant_id,GT,altIndex)
                        # getGTtime+=time.time()-starttime

                        if len(self.geno_info)>0:
                            # starttime=time.time()
                            self.rowData.extend([[variant_id,idx,self.chunk["calldata/DP"][i][idx],self.chunk["calldata/GQ"][i][idx]] for idx in range(self.start_sample,self.end_sample)])
                            # self.rowData.extend([[variant_id,idx]+[self.chunk[field][i][idx] for field in self.fields] for idx in range(self.start_sample,self.end_sample)])
                            # self.getInfoTable(variant_id,infoDict,altIndex)
                            # getInfotime+=time.time()-starttime      
                    else:
                        rec=[str(chr),str(pos),ref,alt]  
                        msg=normalize_variant(RefGenome(self.build).crr, rec, 0, 1, 2, 3)
                        if tuple((rec[0], rec[2], rec[3])) in self.variantIndex:
                            variant_id  = self.variantIndex[tuple((rec[0], rec[2], rec[3]))][rec[1]][0]
                            starttime=time.time()
                            self.getGT(variant_id,GT,altIndex)
                            getGTtime+=time.time()-starttime
                            if len(self.geno_info)>0:
                                # starttime=time.time()
                                self.rowData.extend([[variant_id,idx,self.chunk["calldata/DP"][i][idx],self.chunk["calldata/GQ"][i][idx]] for idx in range(self.start_sample,self.end_sample)])
                                # self.rowData.extend([[variant_id,idx]+[self.chunk[field][i][idx] for field in self.fields] for idx in range(self.start_sample,self.end_sample)])
                                # self.getInfoTable(variant_id,infoDict,altIndex)
                                # getInfotime+=time.time()-starttime
        self.writeIntoHDF(chr)
   
        
      


    def run(self): 
        genoCount=0
        genoDict={}
        if len(self.geno_info)>0:
            for info in self.geno_info:
                genoDict[info.name]=0

        self.start_count = self.variant_count.value   
        self.processFile()
     

      

def updateSample(importer,start_sample,end_sample,sample_ids,names,allNames,HDF5fileName):
    cur=importer.db.cursor()
    
    for id in range(start_sample,end_sample):
        sql="UPDATE sample SET HDF5=? WHERE sample_id=? and sample_name=?"
        # task=(HDF5fileName,allNames[names[id]],sample_ids[id],names[id])
        task=(HDF5fileName,sample_ids[id],names[id])
        cur.execute(sql,task)
        
    


def manageHDF5(importer,allNames={}):

    cur=importer.db.cursor()
    sql="ALTER TABLE sample ADD COLUMN HDF5 CHAR(25)"
    try:
        cur.execute(sql)
    except:
        pass
   
    sql="SELECT sample_name,sample_id from sample"

  
    for rec in cur.execute(sql):
        if rec[0] not in allNames:
            allNames[rec[0]]=int(rec[1])
     
    return allNames




def importGenotypesInParallel(importer,num_sample=0):

    allNames=manageHDF5(importer)

    for count, input_filename in enumerate(importer.files):
        
        env.logger.info('{} variants from {} ({}/{})'.format('Importing', input_filename, count + 1, len(importer.files)))
        importer.importVariant(input_filename)

        env.logger.info('{:,} new variants {}{}{} from {:,} lines are imported.'\
            .format(importer.count[2], "(" if importer.count[2] else '', 
                ', '.join(['{:,} {}'.format(x, y) for x, y in \
                    zip(importer.count[3:8], ['SNVs', 'insertions', 'deletions', 'complex variants', 'unsupported']) if x > 0]),
                    ")" if importer.count[2] else '', importer.count[0]))
        # genotypes?
        if importer.genotype_field:
            importer.prober.reset()
        # if there are samples?
        sample_ids, genotype_status,names = importer.getSampleIDs(input_filename)
        if len(sample_ids) == 0:
            continue
        allNames=manageHDF5(importer,allNames)
        sample_ids=[int(allNames[name]) for name in names]
       
        workload=None
       
        #determine number of samples to be processed in each process
        if num_sample>0:
            #if number of samples in each HDF5 is determined
            num_files=int(len(sample_ids)/num_sample)
            workload=[num_sample]*num_files
            if len(sample_ids)%num_sample!=0:
                workload.append(len(sample_ids)%num_sample)
            if len(workload)<importer.jobs:
                importer.jobs=len(workload)
        else:
            if len(sample_ids)<50:
                importer.jobs=1
                workload=[len(sample_ids)]
            else:
                workload = [int(float(len(sample_ids)) / importer.jobs)] * importer.jobs
          
            unallocated = max(0, len(sample_ids) - sum(workload))
            for i in range(unallocated):
                workload[i % importer.jobs] += 1


        env.logger.debug("work load {}".format(workload))
        numTasks=len(workload)
        # importers = [None] * numProcess
        importers=[None]*importer.jobs
        variant_import_count = [Value('L', 0) for x in range(numTasks)]


        for i in range(len(importer.count)):
            importer.total_count[i] += importer.count[i]
            importer.count[i] = 0
        

        # env.logger.debug('Workload of processes: {}'.format(workload))

        line_no = 0
        num_lines=0    
        readQueue=[]

        taskQueue=queue.Queue()
        task=None


        compression='gzip'
        compression_opts=1
        shuffle=False
        overwrite=False
        vlen=True
        fields=None
        types=None
        numbers=None
      
        if (len(importer.genotype_info)>0):
            # fields=['variants/ID','variants/REF','variants/ALT','variants/POS','variants/CHROM','calldata/GT','calldata/GQ','calldata/DP','calldata/AD','calldata/PL']
            # types={'calldata/PL':'i1'}
            # numbers={'calldata/PL':6}
            fields=['variants/ID','variants/REF','variants/ALT','variants/POS','variants/CHROM','calldata/GT','calldata/GQ','calldata/DP']
        else:
            fields=['variants/ID','variants/REF','variants/ALT','variants/POS','variants/CHROM','calldata/GT']

        alt_number=DEFAULT_ALT_NUMBER
        fills=None
        region=None
        tabix='tabix'
        samples=None
        transformers=None
        buffer_size=DEFAULT_BUFFER_SIZE
        chunk_length=DEFAULT_CHUNK_LENGTH
        chunk_width=DEFAULT_CHUNK_WIDTH

        log=None
       
        _, samples, headers, it = iter_vcf_chunks(
                input_filename, fields=fields, types=types, numbers=numbers, alt_number=alt_number,
                buffer_size=buffer_size, chunk_length=chunk_length, fills=fills, region=region,
                tabix=tabix, samples=samples, transformers=transformers
            )

        #Put tasks in the queue first
        for job in range(numTasks):      
            if workload[job] == 0:
                continue
            # readQueue.append(mpQueue())

        prog = ProgressBar('Importing genotypes', importer.total_count[2],initCount=0)
        lines=0
        # print(fields[6:])
        for chunk, _, _, _ in it:     
            
            start_sample =0
            for job in range(numTasks):
                # readQueue[job].put(chunk)
                if workload[job] == 0:
                    continue
                end_sample = min(start_sample + workload[job], len(sample_ids))
                if end_sample <= start_sample:
                    continue
                HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
                updateSample(importer,start_sample,end_sample,sample_ids,names,allNames,HDFfile_Merge)
                taskQueue.put(HDF5GenotypeImportWorker(chunk, importer.variantIndex, start_sample, end_sample, 
                    sample_ids,variant_import_count[job], job, importer.genotype_info,HDFfile_Merge,importer.build,fields[6:]))
                start_sample = end_sample   
            while taskQueue.qsize()>0:
                for i in range(importer.jobs):    
                    if importers[i] is None or not importers[i].is_alive():     
                        task=taskQueue.get()
                        importers[i]=task
                        importers[i].start()           
                        break 
            for worker in importers:
                worker.join()
            lines+=chunk_length
            prog.update(lines)

        prog.done()

        # # Tried to generate iterator for sample list, but it seems not faster than parsing all samples
        # iterList=[]
        # dbList=[]
        # start_sample=0
        # for job in range(numTasks):
        #     if workload[job]==0:
        #         continue
        #     end_sample = min(start_sample + workload[job], len(sample_ids))
        #     if end_sample <= start_sample:
        #         continue
        #     samples=range(start_sample,end_sample)
        #     print(start_sample,end_sample)
        #     HDFfile_Merge="tmp_"+str(allNames[names[start_sample]])+"_"+str(allNames[names[end_sample-1]])+"_genotypes.h5"
        #     _, samples, headers, it = iter_vcf_chunks(
        #         input_filename, fields=fields, types=types, numbers=numbers, alt_number=alt_number,
        #         buffer_size=buffer_size, chunk_length=chunk_length, fills=fills, region=region,
        #         tabix=tabix, samples=samples, transformers=transformers)
        #     iterList.append(it)
        #     dbList.append(HDFfile_Merge)
        #     start_sample=end_sample

        # prog = ProgressBar('Importing genotypes', importer.total_count[2],initCount=0)

        # while True:
        #     for iterID,it in enumerate(iterList):
        #         try:
        #             chunk, _, _, _=next(it)
        #         except StopIteration:
        #             prog.done()
        #             sys.exit(1)
                
        #         start_sample=int(dbList[iterID].split("_")[1])-1
        #         end_sample=int(dbList[iterID].split("_")[2])
        #         print(start_sample,end_sample)
        #         taskQueue.put(HDF5GenotypeImportWorker(chunk, importer.variantIndex, start_sample, end_sample, 
        #                 sample_ids,variant_import_count[iterID], iterID, importer.genotype_info,dbList[iterID],importer.build))
        #     while taskQueue.qsize()>0:
        #         for i in range(importer.jobs):    
        #             if importers[i] is None or not importers[i].is_alive():     
        #                 task=taskQueue.get()
        #                 importers[i]=task
        #                 importers[i].start()           
        #                 break 
        #     for worker in importers:
        #         worker.join()

        #     prog.update(chunk_length)

   
