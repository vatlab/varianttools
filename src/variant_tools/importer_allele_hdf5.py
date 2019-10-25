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
"""
Extract data from VCF files.

This module contains Functions for extracting data from Variant Call Format (VCF) files and loading
into NumPy arrays, NumPy files, HDF5 files or Zarr array stores.

"""
import glob
import gzip
import os
import queue
import re
import subprocess
import time
import warnings
from collections import namedtuple
from multiprocessing import Process, Value
from shutil import copyfile

import numpy as np

from variant_tools.io_vcf_read import FileInputStream, VCFChunkIterator

from .accessor import Engine_Storage, Engine_Access
from .utils import DatabaseEngine, ProgressBar, env

# expose some names from cython extension
# from variant_tools.io_vcf_read import (  # noqa: F401
#     ANNTransformer, ANN_AA_LENGTH_FIELD, ANN_AA_POS_FIELD, ANN_ANNOTATION_FIELD,
#     ANN_ANNOTATION_IMPACT_FIELD, ANN_CDNA_LENGTH_FIELD, ANN_CDNA_POS_FIELD, ANN_CDS_LENGTH_FIELD,
#     ANN_CDS_POS_FIELD, ANN_DISTANCE_FIELD, ANN_FEATURE_ID_FIELD, ANN_FEATURE_TYPE_FIELD, ANN_FIELD,
#     ANN_FIELDS, ANN_GENE_ID_FIELD, ANN_GENE_NAME_FIELD, ANN_HGVS_C_FIELD, ANN_HGVS_P_FIELD,
#     ANN_RANK_FIELD, ANN_TRANSCRIPT_BIOTYPE_FIELD
# )

DEFAULT_BUFFER_SIZE = 2**14
DEFAULT_CHUNK_LENGTH = 2**11
DEFAULT_CHUNK_WIDTH = 2**6
DEFAULT_ALT_NUMBER = 3

READ_EXISTING_CHUNK_LENGTH = 2**10


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
        message = ('%s %s rows in %.2fs; chunk in %.2fs (%s rows/s)' %
                   (prefix, n_variants, elapsed, elapsed_chunk,
                    int(chunk_length // elapsed_chunk)))
        if chrom:
            message += '; %s:%s' % (chrom, pos)
        print(message, file=log)
        log.flush()
        yield chunk, chunk_length, chrom, pos
        before_chunk = after_chunk
    after_all = time.time()
    elapsed = after_all - before_all
    print(
        '%s all done (%s rows/s)' % (prefix, int(n_variants // elapsed)),
        file=log)
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


def _hdf5_setup_datasets(chunk, root, chunk_length, chunk_width, compression,
                         compression_opts, shuffle, overwrite, headers, vlen):
    import h5py

    # handle no input
    if chunk is None:
        raise RuntimeError('input file has no data?')

    # setup datasets
    # keys = sorted(chunk.keys())
    keys = ["variants/ID", "calldata/GT"]
    for k in keys:

        # obtain initial data
        data = chunk[k]

        # determine chunk shape
        if data.ndim == 1:
            chunk_shape = (chunk_length,)
        else:
            chunk_shape = (chunk_length, min(chunk_width,
                                             data.shape[1])) + data.shape[2:]

        # create dataset
        group, name = k.split('/')
        if name in root[group]:
            if overwrite:
                del root[group][name]
            else:
                raise ValueError(
                    'dataset exists at path %r; use overwrite=True to replace' %
                    k)

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
        root[group].create_dataset(
            name,
            shape=shape,
            maxshape=maxshape,
            chunks=chunk_shape,
            dtype=dt,
            compression=compression,
            compression_opts=compression_opts,
            shuffle=shuffle)

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
    kwds = dict(
        fields=fields,
        types=types,
        numbers=numbers,
        alt_number=alt_number,
        chunk_length=chunk_length,
        fills=fills,
        samples=samples)

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
                warnings.warn(
                    'tabix not found, falling back to scanning to region')
                fileobj = gzip.open(input, mode='rb')
                close = True

            except Exception as e:
                warnings.warn(
                    'error occurred attempting tabix (%s); falling back to '
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

        # else:
        #     warnings.warn('%r INFO header not found' % name)

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

        elif f in [
                'INFO', 'INFO*', 'INFO/*', 'variants/INFO', 'variants/INFO*',
                'variants/INFO/*'
        ]:
            _add_all_info_fields(normed_fields, headers)

        elif f in [
                'FILTER', 'FILTER*', 'FILTER/*', 'FILTER_*', 'variants/FILTER',
                'variants/FILTER*', 'variants/FILTER/*', 'variants/FILTER_*'
        ]:
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
    'variants/NS': 'i1',
    'calldata/GT': 'genotype/i1',
    'calldata/GQ': 'f4',
    'calldata/HQ': 'i1',
    'calldata/DP': 'i2',
    'calldata/GD': 'i2',
    'calldata/AD': 'i2',
    'calldata/MQ0': 'i2',
    'calldata/MQ': 'f2',
    'calldata/PL': 'i2',
    'calldata/NS': 'i1'
}


def _normalize_types(types, fields, headers):

    # normalize user-provided types
    if types is None:
        types = dict()
    types = {
        _normalize_field_prefix(f, headers): _normalize_type(t)
        for f, t in types.items()
    }

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
    'calldata/GD': 1,
    'calldata/GT': 1,
    'calldata/GQ': 1,
    'calldata/HQ': 2,
    'calldata/AD': 2,
    'calldata/PL': 4,
    'calldata/MQ0': 1,
    'calldata/MQ': 1,
    'calldata/NS': 1,
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
            warnings.warn('error parsing %r as number for field %r' %
                          (n, field))
        return 1


def _normalize_numbers(numbers, fields, headers, alt_number):

    # normalize field prefixes
    if numbers is None:
        numbers = dict()
    numbers = {
        _normalize_field_prefix(f, headers): n for f, n in numbers.items()
    }

    # setup output
    normed_numbers = dict()

    for f in fields:

        group, name = f.split('/')

        if f in numbers:
            normed_numbers[f] = _normalize_number(f, numbers[f], alt_number)

        elif f in default_numbers:
            normed_numbers[f] = _normalize_number(f, default_numbers[f],
                                                  alt_number)

        elif group == 'variants':

            if name in ['numalt', 'svlen', 'is_snp']:
                # computed fields, special case (for svlen, number depends on ALT)
                continue

            elif name.startswith('FILTER_'):
                normed_numbers[f] = 0

            elif name in headers.infos:
                normed_numbers[f] = _normalize_number(
                    f, headers.infos[name]['Number'], alt_number)

            else:
                # fall back to 1
                normed_numbers[f] = 1
                warnings.warn('no number for field %r, assuming 1' % f)

        elif group == 'calldata':

            if name in headers.formats:
                normed_numbers[f] = _normalize_number(
                    f, headers.formats[name]['Number'], alt_number)

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
    fills = {_normalize_field_prefix(f, headers): v for f, v in fills.items()}

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


def _iter_vcf_stream(stream, fields, types, numbers, alt_number, chunk_length,
                     fills, region, samples):

    # read VCF headers
    headers = _read_vcf_headers(stream)

    # setup samples
    samples, loc_samples = _normalize_samples(
        samples=samples, headers=headers, types=types)

    # setup fields to read
    if fields is None:

        # choose default fields
        fields = list()
        _add_all_fixed_variants_fields(fields)
        fields.append('variants/FILTER_PASS')
        if len(samples) > 0 and 'GT' in headers.formats:
            fields.append('calldata/GT')

    else:
        fields = _normalize_fields(
            fields=fields, headers=headers, samples=samples)

    # setup data types
    types = _normalize_types(types=types, fields=fields, headers=headers)

    # setup numbers (a.k.a., arity)
    numbers = _normalize_numbers(
        numbers=numbers, fields=fields, headers=headers, alt_number=alt_number)

    # setup fills
    fills = _normalize_fills(fills=fills, fields=fields, headers=headers)

    # setup chunks iterator
    chunks = VCFChunkIterator(
        stream,
        chunk_length=chunk_length,
        headers=headers,
        fields=fields,
        types=types,
        numbers=numbers,
        fills=fills,
        region=region,
        loc_samples=loc_samples)

    return fields, samples, headers, chunks


# pre-compile some regular expressions
_re_filter_header = \
    re.compile('##FILTER=<ID=([^,]+),Description="([^"]+)">')
_re_info_header = \
    re.compile('##INFO=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)">')
_re_format_header = \
    re.compile('##FORMAT=<ID=([^,]+),Number=([^,]+),Type=([^,]+),Description="([^"]+)">')

VCFHeaders = namedtuple('VCFHeaders',
                        ['headers', 'filters', 'infos', 'formats', 'samples'])


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
        raise RuntimeError(
            'VCF file is missing mandatory header line ("#CHROM...")')

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

    def __init__(self, chunk, variantIndex, start_sample, end_sample,
                 sample_ids, variant_count, proc_index, dbLocation,
                 genotype_info, build):

        Process.__init__(self, name='GenotypeImporter')
        # self.daemon=True
        self.chunk = chunk
        self.variantIndex = variantIndex
        self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample = start_sample
        self.end_sample = end_sample

        self.indptr = []
        self.indices = []
        self.data = []
        self.rownames = []
        self.sample_ids = sample_ids
        self.firstID = 0
        # if sample_ids[0]!=1 and start_sample!=0:
        #     self.firstID=sample_ids[0]
        #     self.start_sample=self.start_sample+1
        #     self.end_sample=self.end_sample+1
        # self.colnames=[self.sample_ids[i-self.firstID] for i in range(self.start_sample,self.end_sample)]
        self.colnames = [
            self.sample_ids[i]
            for i in range(self.start_sample, self.end_sample)
        ]
        self.genoCount = 0
        self.dbLocation = dbLocation
        self.build = build
        self.info = {}
        self.rowData = []
        self.info["GT"] = []
        self.info["Mask"] = []
        self.namedict = {}
        self.geno_info = []
        if "GT" in genotype_info:
            genotype_info.remove("GT")
        if len(genotype_info) > 0:
            for info in genotype_info:
                #indptr,indices,data,shape,rownames
                if not isinstance(info, str):
                    self.geno_info.append(info.name.replace("_geno", ""))
                else:
                    self.geno_info.append(info.replace("_geno", ""))
        for info in self.geno_info:
            self.info[info] = []
            if "calldata/" + info in self.chunk and np.nansum(
                    self.chunk["calldata/" + info][:10]) > 0:
                self.namedict[info] = "calldata/" + info
            elif "variants/" + info in self.chunk and np.nansum(
                    self.chunk["variants/" + info][:10]) > 0:
                self.namedict[info] = "variants/" + info

    # check io_vcf_read.pyx function vcf_genotype_parse to see the meaning of coding
    def get_geno(self, variant_id, pos, altIndex):
        self.rownames.append(variant_id)
        # print(self.dbLocation,self.start_sample,self.end_sample,self.firstID)

        if "calldata/GT" in self.chunk:

            GT = self.chunk["calldata/GT"][pos, self.start_sample -
                                           self.firstID:self.end_sample -
                                           self.firstID]
            GT = GT.astype(float)
            if altIndex == 0:
                GT[np.logical_or(GT == 3, GT == 4)] = np.nan
            elif altIndex == 1:
                # GT_geno[GT_geno==3]=1
                # GT_geno[GT_geno==4]=2
                GT[(GT != 3) & (GT != 4) & (GT != -1)] = np.nan
                # GT_geno[np.logical_and(GT_geno!=3, GT_geno!=4)]=np.nan
                GT[GT == 3] = 1
                GT[GT == 4] = 2
            GT[GT == -10] = np.nan
            self.info["GT"].append(GT)
            self.info["Mask"].append([1] * len(GT))
        else:
            # GT_geno=[np.nan]
            GT = [-1]
            self.info["GT"].append(GT)
            self.info["Mask"].append([1] * len(GT))

        if len(self.geno_info) > 0:
            # self.rowData.extend([[variant_id,idx,self.chunk["calldata/DP"][i][idx],self.chunk["calldata/GQ"][i][idx]] for idx in range(self.start_sample,self.end_sample)])
            # self.rowData.extend([[variant_id,idx]+[self.chunk[field][i][idx] for field in self.fields] for idx in range(self.start_sample,self.end_sample)])
            # self.getInfoTable(variant_id,infoDict,altIndex)
            for info in self.geno_info:
                if "variants" in self.namedict[info]:
                    self.info[info].append(
                        np.array([self.chunk[self.namedict[info]][pos]]))
                else:
                    self.info[info].append(self.chunk[self.namedict[info]]
                                           [pos, self.start_sample -
                                            self.firstID:self.end_sample -
                                            self.firstID])
                # print(self.namedict[info],info,self.start_sample,self.end_sample,pos,(self.chunk[self.namedict[info]][pos,self.start_sample-self.firstID:self.end_sample-self.firstID]))

    def writeIntoHDF(self, chr):
        storageEngine = Engine_Storage.choose_storage_engine(self.dbLocation)
        shape = np.array([len(self.rownames), len(self.colnames)])
        storageEngine.store(np.array(self.info["GT"]), chr, "GT")
        storageEngine.store(np.array(self.info["Mask"]), chr, "Mask")
        storageEngine.store(np.array(self.rownames), chr, "rownames")
        rowmask = np.zeros(len(self.rownames), dtype=np.bool)
        storageEngine.store(np.array(rowmask), chr, "rowmask")

        if not storageEngine.checkGroup(chr, "colnames"):
            storageEngine.store(np.array(self.colnames), chr, "colnames")
            colmask = np.zeros(len(self.colnames), dtype=np.bool)
            storageEngine.store(np.array(colmask), chr, "samplemask")

        storageEngine.store(shape, chr, "shape")

        self.info["GT"] = []
        self.info["Mask"] = []

        if len(self.geno_info) > 0:
            for info in self.geno_info:
                storageEngine.store_genoInfo(
                    np.array(self.info[info]), chr, info)
                self.info[info] = []

        storageEngine.close()
        self.rownames = []

    def run(self):

        prev_chr = self.chunk["variants/CHROM"][0].replace("chr", "")
        prev_variant_id = -1
        for i in range(len(self.chunk["variants/ID"])):
            chr = self.chunk["variants/CHROM"][i].replace("chr", "")
            if chr != prev_chr:
                self.writeIntoHDF(prev_chr)
                prev_chr = chr
            ref = self.chunk["variants/REF"][i]
            pos = self.chunk["variants/POS"][i]
            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt = self.chunk["variants/ALT"][i][altIndex]
                if alt != "":
                    if tuple((chr, ref, alt)) in self.variantIndex:
                        variant_id = self.variantIndex[tuple(
                            (chr, ref, alt))][pos][0]

                        if variant_id != prev_variant_id:
                            self.get_geno(variant_id, i, altIndex)
                            prev_variant_id = variant_id

                    else:
                        rec = [str(chr), str(pos), ref, alt]
                        #msg = normalize_variant(
                        #    RefGenome(self.build).crr, rec, 0, 1, 2, 3)
                        if tuple((rec[0], rec[2], rec[3])) in self.variantIndex:
                            variant_id = self.variantIndex[tuple(
                                (rec[0], rec[2], rec[3]))][rec[1]][0]
                            if variant_id != prev_variant_id:
                                self.get_geno(variant_id, i, altIndex)
                                prev_variant_id = variant_id
        self.writeIntoHDF(chr)


class HDF5GenotypeSortWorker(Process):

    def __init__(self, chunk, importer, start_sample, end_sample, sample_ids,
                 variant_count, proc_index, dbLocation, estart, eend, efirst):

        Process.__init__(self, name='GenotypeSorter')
        # self.daemon=True
        self.importer = importer
        self.chunk = chunk
        self.variantIndex = importer.variantIndex
        self.variant_count = variant_count
        self.proc_index = proc_index
        self.start_sample = start_sample
        self.end_sample = end_sample
        # self.geno_info=geno_info
        self.indptr = []
        self.indices = []
        self.data = []
        self.rownames = []
        self.sample_ids = sample_ids
        self.firstID = 0
        # if sample_ids[0]!=1 and start_sample!=0:
        #     self.firstID=sample_ids[0]
        #     self.start_sample=self.start_sample+1
        #     self.end_sample=self.end_sample+1
        # self.colnames=[self.sample_ids[i-self.firstID] for i in range(self.start_sample,self.end_sample)]
        self.colnames = [
            self.sample_ids[i]
            for i in range(self.start_sample, self.end_sample)
        ]
        self.genoCount = 0
        self.dbLocation = dbLocation

        self.originLocation = dbLocation.replace("_sort", "")
        self.build = importer.build
        self.info = {}
        self.rowData = []
        self.info["GT"] = []
        self.info["Mask"] = []
        self.namedict = {}
        self.geno_info = []
        self.estart = estart
        self.eend = eend
        self.efirst = efirst
        self.countcheck = 0
        self.readcheck = 0
        # if self.dbLocation==self.checkfileName:
        #     print("new process")
        if len(importer.genotype_info) > 0:
            for info in importer.genotype_info:
                #indptr,indices,data,shape,rownames
                if not isinstance(info, str):
                    self.geno_info.append(info.name)
                else:
                    self.geno_info.append(info)
        for info in self.geno_info:
            self.info[info] = []
            self.namedict[info] = "calldata/" + info

    # check io_vcf_read.pyx function vcf_genotype_parse to see the meaning of coding
    def get_geno(self, variant_id, pos, altIndex):
        self.readcheck += 1
        self.rownames.append(variant_id)

        GT = self.chunk["calldata/GT"][pos, self.start_sample -
                                       self.firstID:self.end_sample -
                                       self.firstID]
        GT = GT.astype(float)
        if altIndex == 0:
            GT[np.logical_or(GT == 3, GT == 4)] = np.nan
        elif altIndex == 1:
            GT[np.logical_and(GT != 3, GT != 4)] = np.nan
            GT[GT == 3] = 1
            GT[GT == 4] = 2
        # GT_geno[GT_geno==-10]=np.nan
        self.info["GT"].append(GT)
        self.info["Mask"].append([1] * len(GT))
        if len(self.geno_info) > 0:
            # self.rowData.extend([[variant_id,idx,self.chunk["calldata/DP"][i][idx],self.chunk["calldata/GQ"][i][idx]] for idx in range(self.start_sample,self.end_sample)])
            # self.rowData.extend([[variant_id,idx]+[self.chunk[field][i][idx] for field in self.fields] for idx in range(self.start_sample,self.end_sample)])
            # self.getInfoTable(variant_id,infoDict,altIndex)
            for info in self.geno_info:
                self.info[info].append(self.chunk[self.namedict[info]]
                                       [pos, self.start_sample -
                                        self.firstID:self.end_sample -
                                        self.firstID])

    def writeIntoHDF(self, chr):
        # if self.dbLocation==self.checkfileName:
        print("writing ", self.dbLocation, self.countcheck, self.readcheck, chr,
              len(self.rownames))

        storageEngine = Engine_Storage.choose_storage_engine(self.dbLocation)
        shape = np.array([len(self.rownames), len(self.colnames)])
        storageEngine.store(np.array(self.info["GT"]), chr, "GT")
        storageEngine.store(np.array(self.info["Mask"]), chr, "Mask")
        storageEngine.store(np.array(self.rownames), chr, "rownames")
        rowmask = np.zeros(len(self.rownames), dtype=np.bool)
        storageEngine.store(np.array(rowmask), chr, "rowmask")

        if not storageEngine.checkGroup(chr, "colnames"):
            storageEngine.store(np.array(self.colnames), chr, "colnames")
            colmask = np.zeros(len(self.colnames), dtype=np.bool)
            storageEngine.store(np.array(colmask), chr, "samplemask")

        storageEngine.store(shape, chr, "shape")

        self.info["GT"] = []
        self.info["Mask"] = []

        if len(self.geno_info) > 0:
            for info in self.geno_info:
                storageEngine.store_genoInfo(
                    np.array(self.info[info]), chr, info)
                self.info[info] = []

        storageEngine.close()
        self.rownames = []

    def get_existing_geno(self, chr, startPos, endPos):
        accessEngine = Engine_Access.choose_access_engine(self.originLocation)
        node = accessEngine.file.get_node("/chr" + str(chr))
        colnames = node.colnames[:].tolist()
        #shape = node.shape[:]
        colpos = list(map(lambda x: colnames.index(x), colnames))
        allgenoinfo = {}
        # if endPos>shape[0]:
        #     print(endPos,shape[0])
        #     endPos=shape[0]
        rownames, colnames, genotype = accessEngine.filter_on_genotypes(
            "", chr, node, "GT", startPos, endPos, colpos, [])
        allgenoinfo["GT"] = genotype
        if len(self.geno_info) > 0:
            for info in self.geno_info:
                _, _, genoinfo = accessEngine.filter_on_genotypes(
                    "", chr, node, info, startPos, endPos, colpos, [])
                allgenoinfo[info] = genoinfo

        # cur=self.importer.db.cursor()
        proj_file = self.importer.proj.name + '.proj'
        #
        # create a temporary directory
        db = DatabaseEngine()
        db.connect(proj_file)
        cur = db.cursor()

        sql = "SELECT variant_id,pos from variant where variant_id in (" + ",".join(
            [str(rowname) for rowname in rownames]) + ")"
        rowpos = []
        updatedrownames = []
        count = 0
        rows = cur.execute(sql)
        for row in rows:
            updatedrownames.append(row[0])
            rowpos.append(row[1])
            count += 1
        accessEngine.close()
        db.close()
        if count != len(rownames):
            print("alert", self.dbLocation, len(rowpos), len(updatedrownames),
                  len(rownames), startPos, endPos, count)

        return rowpos, updatedrownames, allgenoinfo

    def get_remaining(self, chr, startPos):
        accessEngine = Engine_Access.choose_access_engine(self.originLocation)
        node = accessEngine.file.get_node("/chr" + str(chr))
        rownames = node.rownames[:].tolist()
        colnames = node.colnames[:].tolist()
        colpos = list(map(lambda x: colnames.index(x), colnames))

        while startPos < len(rownames):
            endPos = startPos + READ_EXISTING_CHUNK_LENGTH
            if endPos >= len(rownames):
                endPos = len(rownames)
            chunk_rownames, colnames, genotype = accessEngine.filter_on_genotypes(
                "", chr, node, "GT", startPos, endPos, colpos, [])
            for rowname in chunk_rownames:
                self.rownames.append(rowname)
            for rec in genotype:
                self.info["GT"].append(rec)
                self.info["Mask"].append([1.0] * len(rec))
            if len(self.geno_info) > 0:
                for info in self.geno_info:
                    _, _, genoinfo = accessEngine.filter_on_genotypes(
                        "", chr, node, info, startPos, endPos, colpos, [])
                    for rec in genoinfo:
                        self.info[info].append(rec)
            self.writeIntoHDF(chr)
            # if self.dbLocation==self.checkfileName:
            #     print("remaining",startPos,endPos)
            startPos = endPos
        accessEngine.close()

    def run(self):
        prev_chr = self.chunk["variants/CHROM"][0].replace("chr", "")
        if self.efirst.value == True:
            storageEngine = Engine_Storage.choose_storage_engine(
                self.dbLocation)
            storageEngine.removeNode([prev_chr])
            # command = ["ptrepack", "-o", "--chunkshape=auto", "--propindexes", self.dbLocation, self.dbLocation+"test"]
            # call(command)
            storageEngine.close()
            self.efirst.value = False
        erowpos, erownames, egenoinfo = self.get_existing_geno(
            prev_chr, self.estart.value, self.eend.value)
        last = 0
        chunksize = len(self.chunk["variants/ID"])

        for i in range(len(self.chunk["variants/ID"])):
            chr = self.chunk["variants/CHROM"][i].replace("chr", "")
            if chr != prev_chr:
                if self.estart.value != self.eend.value:
                    self.estart.value = self.estart.value + last
                    self.get_remaining(prev_chr, self.estart.value)
                else:
                    self.writeIntoHDF(prev_chr)
                prev_chr = chr
                storageEngine = Engine_Storage.choose_storage_engine(
                    self.dbLocation)
                storageEngine.removeNode([prev_chr])
                storageEngine.close()
                self.estart.value = 0
                self.eend.value = self.estart.value + READ_EXISTING_CHUNK_LENGTH
                erowpos, erownames, egenoinfo = self.get_existing_geno(
                    prev_chr, self.estart.value, self.eend.value)
                last = 0
            ref = self.chunk["variants/REF"][i]
            pos = self.chunk["variants/POS"][i]

            while last < len(erowpos) and pos > erowpos[last]:
                self.countcheck += 1
                # if self.dbLocation==self.checkfileName:
                #     print("add old file",erowpos[last],last,i)
                self.rownames.append(erownames[last])
                for info in egenoinfo.keys():
                    self.info[info].append(egenoinfo[info][last, :])
                self.info["Mask"].append([1] * len(egenoinfo["GT"][last, :]))
                last += 1
            # if self.dbLocation==self.checkfileName:
            #         print("before",last,len(erowpos))

            # if last==len(erowpos) and last!=0:
            if last == READ_EXISTING_CHUNK_LENGTH:

                ##left over in new file
                self.estart.value = self.eend.value
                self.eend.value += READ_EXISTING_CHUNK_LENGTH
                # if self.dbLocation==self.checkfileName:
                #     print("read more from existing",self.estart.value,self.eend.value)
                erowpos, erownames, egenoinfo = self.get_existing_geno(
                    prev_chr, self.estart.value, self.eend.value)
                last = 0
                # if self.dbLocation==self.checkfileName:
                #     print("after",last,len(erowpos))
            elif last == len(erowpos) and last != 0:
                self.estart.value = self.eend.value

            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt = self.chunk["variants/ALT"][i][altIndex]
                if alt != "":
                    if tuple((chr, ref, alt)) in self.variantIndex:
                        variant_id = self.variantIndex[tuple(
                            (chr, ref, alt))][pos][0]
                        self.get_geno(variant_id, i, altIndex)

                    else:
                        rec = [str(chr), str(pos), ref, alt]
                        #msg = normalize_variant(
                        #    RefGenome(self.build).crr, rec, 0, 1, 2, 3)
                        if tuple((rec[0], rec[2], rec[3])) in self.variantIndex:
                            variant_id = self.variantIndex[tuple(
                                (rec[0], rec[2], rec[3]))][rec[1]][0]
                            self.get_geno(variant_id, i, altIndex)
                    # if self.dbLocation==self.checkfileName:
                    #     print("1",pos,variant_id)
        if last < len(erowpos):
            self.estart.value = self.estart.value + last
            self.eend.value = self.estart.value + READ_EXISTING_CHUNK_LENGTH
            # if self.dbLocation==self.checkfileName:
            #     print("process end, restart from new position",self.estart.value,self.eend.value)
        self.writeIntoHDF(chr)
        if chunksize < DEFAULT_CHUNK_LENGTH:
            self.get_remaining(prev_chr, self.estart.value)


def updateSample(cur, start_sample, end_sample, sample_ids, names, allNames,
                 HDF5fileName):
    firstID = 0
    # if sample_ids[0]!=1 and start_sample!=0:
    #     firstID=sample_ids[0]
    #     start_sample=start_sample+1
    #     end_sample=end_sample+1
    #     adjust=1
    #print(start_sample,end_sample,sample_ids,HDF5fileName)
    for id in range(start_sample, end_sample):
        try:
            sql = "UPDATE sample SET HDF5=? WHERE sample_id=? and sample_name=?"
            # task=(HDF5fileName,allNames[names[id]],sample_ids[id],names[id])
            task = (HDF5fileName, sample_ids[id - firstID], names[id - firstID])
            cur.execute(sql, task)
        except Exception as e:
            print(e)


def remove_duplicate_samples(cur, original_sample_ids):
    try:
        # sql="DELETE from SAMPLE where sample_id in (%s)" % ','.join(['?'] * len(original_sample_ids))
        # print(sql)
        # cur.execute(sql,original_sample_ids)
        for id in original_sample_ids:
            sql = "DELETE from SAMPLE where sample_id in (?)"
            cur.execute(sql, [str(id)])
        cur.execute('UPDATE project SET value=1 WHERE name="multiVCF"')
    except Exception as e:
        print(e)


def manageHDF5(cur, allNames={}):

    sql = "ALTER TABLE sample ADD COLUMN HDF5 CHAR(25)"
    try:
        cur.execute(sql)
    except:
        pass

    sql = "SELECT sample_name,sample_id from sample"

    for rec in cur.execute(sql):
        if rec[0] not in allNames:
            allNames[rec[0]] = int(rec[1])

    return allNames


def importGenotypesInParallel(importer, num_sample=0):
    cur = importer.db.cursor()
    allNames = manageHDF5(cur)

    for count, input_filename in enumerate(importer.files):

        env.logger.info('{} variants from {} ({}/{})'.format(
            'Importing', input_filename, count + 1, len(importer.files)))
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
        original_sample_ids, genotype_status, names = importer.getSampleIDs(
            input_filename)
        if len(original_sample_ids) == 0:
            continue
        allNames = manageHDF5(cur, allNames)
        sample_ids = [int(allNames[name]) for name in names]
        if original_sample_ids[0] != sample_ids[0]:
            remove_duplicate_samples(cur, original_sample_ids)
        workload = None

        #determine number of samples to be processed in each process
        if num_sample > 0:
            #if number of samples in each HDF5 is determined
            num_files = int(len(sample_ids) / num_sample)
            workload = [num_sample] * num_files
            if len(sample_ids) % num_sample != 0:
                workload.append(len(sample_ids) % num_sample)
            if len(workload) < importer.jobs:
                importer.jobs = len(workload)
        else:
            if len(sample_ids) < 50:
                importer.jobs = 1
                workload = [len(sample_ids)]
            else:
                workload = [int(float(len(sample_ids)) / importer.jobs)
                           ] * importer.jobs

            unallocated = max(0, len(sample_ids) - sum(workload))
            for i in range(unallocated):
                workload[i % importer.jobs] += 1

        env.logger.debug("work load {}".format(workload))
        numTasks = len(workload)
        # importers = [None] * numProcess
        importers = [None] * importer.jobs
        variant_import_count = [Value('L', 0) for x in range(numTasks)]

        for i in range(len(importer.count)):
            importer.total_count[i] += importer.count[i]
            importer.count[i] = 0

        # env.logger.debug('Workload of processes: {}'.format(workload))

        taskQueue = queue.Queue()
        task = None

        #compression = 'gzip'
        #compression_opts = 1
        #shuffle = False
        #overwrite = False
        #vlen = True
        fields = [
            'variants/ID', 'variants/REF', 'variants/ALT', 'variants/POS',
            'variants/CHROM'
        ]
        types = None
        numbers = None
        if (len(importer.genotype_field) > 0):
            fields.append('calldata/GT')
        if (len(importer.genotype_info) > 0):
            for genoinfo in importer.genotype_info:
                fields.append("calldata/" + genoinfo.name.replace("_geno", ""))
                fields.append("variants/" + genoinfo.name.replace("_geno", ""))
        alt_number = DEFAULT_ALT_NUMBER
        fills = None
        region = None
        tabix = 'tabix'
        samples = None
        transformers = None
        buffer_size = DEFAULT_BUFFER_SIZE
        chunk_length = DEFAULT_CHUNK_LENGTH
        #chunk_width = DEFAULT_CHUNK_WIDTH

        #log = None

        _, samples, headers, it = iter_vcf_chunks(
            input_filename,
            fields=fields,
            types=types,
            numbers=numbers,
            alt_number=alt_number,
            buffer_size=buffer_size,
            chunk_length=chunk_length,
            fills=fills,
            region=region,
            tabix=tabix,
            samples=samples,
            transformers=transformers)
        # print(importer.genotype_info)

        #Put tasks in the queue first
        for job in range(numTasks):
            if workload[job] == 0:
                continue
            # readQueue.append(mpQueue())
        prog = ProgressBar(
            'Importing genotypes', importer.total_count[2], initCount=0)
        lines = 0

        estart = [Value('L', 0) for x in range(numTasks)]
        eend = [Value('L', READ_EXISTING_CHUNK_LENGTH) for x in range(numTasks)]
        efirst = [Value('b', True) for x in range(numTasks)]

        for chunk, _, _, _ in it:
            start_sample = 0
            for job in range(numTasks):

                if workload[job] == 0:
                    continue
                end_sample = min(start_sample + workload[job], len(sample_ids))
                if end_sample <= start_sample:
                    continue
                if not importer.sort:
                    HDFfile_Merge = "tmp_" + str(
                        allNames[names[start_sample]]) + "_" + str(
                            allNames[names[end_sample - 1]]) + "_genotypes.h5"
                    updateSample(cur, start_sample, end_sample, sample_ids,
                                 names, allNames, HDFfile_Merge)
                    taskQueue.put(
                        HDF5GenotypeImportWorker(
                            chunk, importer.variantIndex, start_sample,
                            end_sample, sample_ids, variant_import_count[job],
                            job, HDFfile_Merge, importer.genotype_info,
                            importer.build))

                    # do_work.delay(chunk,importer.variantIndex,start_sample, end_sample,
                    #      sample_ids, job,HDFfile_Merge,importer.genotype_info,importer.build)

                else:
                    originalFile = "tmp_" + str(
                        allNames[names[start_sample]]) + "_" + str(
                            allNames[names[end_sample - 1]]) + "_genotypes.h5"
                    HDFfile_Merge = "tmp_" + str(
                        allNames[names[start_sample]]) + "_" + str(
                            allNames[names[end_sample -
                                           1]]) + "_sort_genotypes.h5"
                    if not os.path.isfile(HDFfile_Merge):
                        copyfile(originalFile, HDFfile_Merge)
                    updateSample(cur, start_sample, end_sample, sample_ids,
                                 names, allNames,
                                 HDFfile_Merge.replace("_sort", ""))

                    taskQueue.put(
                        HDF5GenotypeSortWorker(chunk, importer, start_sample,
                                               end_sample, sample_ids,
                                               variant_import_count[job], job,
                                               HDFfile_Merge, estart[job],
                                               eend[job], efirst[job]))
                start_sample = end_sample
            while taskQueue.qsize() > 0:
                for i in range(importer.jobs):
                    if importers[i] is None or not importers[i].is_alive():
                        task = taskQueue.get()
                        importers[i] = task
                        importers[i].start()
                        break
            for worker in importers:
                if worker is not None:
                    worker.join()
            lines += chunk_length
            prog.update(lines)
        for sortFile in glob.glob("tmp*_sort_genotypes.h5"):
            os.rename(sortFile, sortFile.replace("_sort", ""))

        prog.done()


class HDF5GenotypeUpdateWorker(Process):

    def __init__(self, chunk, variantIndex, start_sample, end_sample,
                 sample_ids, geno_info, HDFfile, build):
        Process.__init__(self, name='GenotypeUpdater')
        # self.daemon=True
        self.chunk = chunk
        self.variantIndex = variantIndex
        self.start_sample = start_sample
        self.end_sample = end_sample
        self.sample_ids = sample_ids
        self.geno_info = geno_info
        self.HDFfile = HDFfile
        self.build = build
        self.info = {}
        self.namedict = {}

        if len(self.geno_info) > 0:
            for info in self.geno_info:
                name = info.name.replace("_geno", "")
                self.info[name] = []
                # if "_geno" in info.name:
                #     self.namedict[info.name]="calldata/"+info.name.replace("_geno","")
                # else:
                #     self.namedict[info.name]="variants/"+info.name.replace("_geno","")
                if "calldata/" + name in self.chunk and np.nansum(
                        self.chunk["calldata/" + name][:10]) > 0:
                    self.namedict[name] = "calldata/" + name
                elif "variants/" + name in self.chunk and np.nansum(
                        self.chunk["variants/" + name][:10]) > 0:
                    self.namedict[name] = "variants/" + name

    def writeIntoHDF(self, chr):
        storageEngine = Engine_Storage.choose_storage_engine(self.HDFfile)
        if len(self.geno_info) > 0:
            for info in self.geno_info:
                # if "_geno" in info.name:
                name = info.name.replace("_geno", "")
                storageEngine.store_genoInfo(
                    np.array(self.info[name]), chr, name)
                # else:
                #     ar=np.array(self.info[info.name])
                #     storageEngine.store_genoInfo(np.reshape(ar,(len(ar),-1)),chr,info.name)
                self.info[name] = []
        storageEngine.close()

    def run(self):
        prev_chr = self.chunk["variants/CHROM"][0].replace("chr", "")
        for i in range(len(self.chunk["variants/ID"])):
            chr = self.chunk["variants/CHROM"][i].replace("chr", "")
            if chr != prev_chr:
                self.writeIntoHDF(prev_chr)
                prev_chr = chr
            #ref = self.chunk["variants/REF"][i]
            #pos = self.chunk["variants/POS"][i]

            for altIndex in range(len(self.chunk["variants/ALT"][i])):
                alt = self.chunk["variants/ALT"][i][altIndex]
                if alt != "":
                    for info in self.geno_info:
                        name = info.name.replace("_geno", "")
                        if "variants" in self.namedict[name]:
                            self.info[name].append(
                                np.array([self.chunk[self.namedict[name]][i]]))
                        else:
                            self.info[name].append(
                                self.chunk[self.namedict[name]]
                                [i, self.start_sample:self.end_sample])

        self.writeIntoHDF(chr)


def UpdateGenotypeInParallel(updater, input_filename, sample_ids, hdf5Files):

    jobs = updater.jobs
    #numTasks = jobs
    updaters = [None] * jobs
    taskQueue = queue.Queue()
    task = None
    variantIndex = updater.proj.createVariantMap('variant', False)

    fields = [
        'variants/ID', 'variants/REF', 'variants/ALT', 'variants/POS',
        'variants/CHROM'
    ]
    if (len(updater.genotype_info) > 0):
        for field in updater.genotype_info:
            fields.append("calldata/" + field.name.replace("_geno", ""))
            fields.append("variants/" + field.name.replace("_geno", ""))
    types = None
    numbers = None
    alt_number = DEFAULT_ALT_NUMBER
    fills = None
    region = None
    tabix = 'tabix'
    samples = None
    transformers = None
    buffer_size = DEFAULT_BUFFER_SIZE
    chunk_length = DEFAULT_CHUNK_LENGTH
    #chunk_width = DEFAULT_CHUNK_WIDTH

    #log = None
    _, samples, headers, it = iter_vcf_chunks(
        input_filename,
        fields=fields,
        types=types,
        numbers=numbers,
        alt_number=alt_number,
        buffer_size=buffer_size,
        chunk_length=chunk_length,
        fills=fills,
        region=region,
        tabix=tabix,
        samples=samples,
        transformers=transformers)

    for HDFfileName in hdf5Files:
        storageEngine = Engine_Storage.choose_storage_engine(HDFfileName)
        for info in updater.genotype_info:
            storageEngine.removeNode([], info.name)
            storageEngine.close()

    for chunk, _, _, _ in it:
        for HDFfileName in hdf5Files:
            # cur=updater.proj.db.cursor()
            # cur.execute('SELECT MIN(sample_id),Max(sample_id) FROM sample where HDF5="{}" '.format(HDFfileName))
            # result=cur.fetchone()
            # start_sample=result[0]-1
            # end_sample=result[1]
            cols = HDFfileName.split("_")
            start_sample = int(cols[1]) - 1
            end_sample = int(cols[2])
            taskQueue.put(
                HDF5GenotypeUpdateWorker(chunk, variantIndex, start_sample,
                                         end_sample, sample_ids,
                                         updater.genotype_info, HDFfileName,
                                         updater.proj.build))
        while taskQueue.qsize() > 0:
            for i in range(jobs):
                if updaters[i] is None or not updaters[i].is_alive():
                    task = taskQueue.get()
                    updaters[i] = task
                    updaters[i].start()
                    break
        for worker in updaters:
            if worker is not None:
                worker.join()
