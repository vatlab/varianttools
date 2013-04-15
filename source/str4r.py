#!/usr/bin/env python
# Copyright Gao Wang 2012
import sys, os, subprocess

if sys.version < '2.3': 
	set = frozenset = tuple
	basestring = str
elif sys.version < '2.4':
	from sets import Set as set, ImmutableSet as frozenset

if sys.version < '3.0':
	_mystr = _mybytes = lambda s:s
else:
	from functools import reduce
	long, basestring, unicode = int, str, str
	_mybytes = lambda s:bytes(s, 'utf8') #'ascii')
	_mystr = lambda s:str(s, 'utf8')
	StringType = str

from types import *

def NoneStr(obj):
    return 'NA'

def BoolStr(obj):
    return obj and 'TRUE' or 'FALSE'

def ReprStr(obj):
    return repr(obj)

def LongStr(obj):
    rtn = repr(obj)
    if rtn[-1] == 'L': rtn = rtn[:-1]
    return rtn

def ComplexStr(obj):
    return repr(obj).replace('j', 'i')

def SeqStr(obj, method='c(', tail=')'):
    if not obj:
        return method + tail
    # detect types
    if isinstance(obj, set):
        obj = list(obj)
    obj_not_none = [x for x in obj if x is not None]
    if len(obj_not_none) == 0:
        # all is None
        return (method + ','.join(list(map(Str4R, obj))) + tail)
    obj0 = obj_not_none[0]
    tp0 = type(obj0)
    simple_types = [str, bool, int, long, float, complex]
    num_types = [int, long, float, complex]
    is_int = tp0 in (int, long) 
    if tp0 not in simple_types and method == 'c(':
        method = 'list('
    else:
        tps = isinstance(obj0, basestring) and [StringType] or num_types
        for i in obj_not_none[1:]:
            tp = type(i)
            if tp not in tps:
                if method == 'c(':
                    method = 'list('
                is_int = False
                break
            elif is_int and tp not in (int, long):
                is_int = False
    # convert
    return (is_int and 'as.integer(' or '') + \
      method + ','.join(list(map(Str4R, obj))) + tail + \
      (is_int and ')' or '')

def DictStr(obj):
    return 'list(' + ','.join(['%s=%s' % (Str4R(a[0]), Str4R(a[1])) for a in list(obj.items())]) + ')'

def OtherStr(obj):
    if hasattr(obj, '__iter__'): # for iterators
        if hasattr(obj, '__len__') and len(obj) <= 10000:
            return SeqStr(list(obj))
        else: # waiting for better solution for huge-size containers
            return SeqStr(list(obj))
    return repr(obj)

base_tps = [type(None), bool, int, long, float, complex, str, unicode, list, tuple, set, frozenset, dict] # use type(None) instead of NoneType since the latter cannot be found in the types module in Python 3
base_tps.reverse()
str_func = {type(None):NoneStr, bool:BoolStr, long:LongStr, int:repr, float:repr, complex:ComplexStr, str:repr, unicode:repr, list:SeqStr, tuple:SeqStr, set:SeqStr, frozenset:SeqStr, dict:DictStr}

def Str4R(obj, method = 'c'):
    '''
    convert a Python basic object into an R object in the form of string.
    '''
    if method != 'c':
        return SeqStr(obj, method=method + '(')
    if type(obj) in str_func:
        return str_func[type(obj)](obj)
    for tp in base_tps:
        if isinstance(obj, tp):
            return str_func[tp](obj)
    return OtherStr(obj)
