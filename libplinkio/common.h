/* =====================================================================================
// 
//  This is a small C and Python library for reading Plink genotype files,
//  written by Mattias Franberg, version 0.2.2 
//  
//  https://bitbucket.org/mattias_franberg/libplinkio
//
//  This software is not licensed or copyrighted. The varianttools developers
//  have been contacting its author and will include the license information when we
//  hear from the author, or replace it with alternative implementation if the author
//  requests for a removal.
// 
 ===================================================================================== */



#ifndef _COMMON_H_
#define _COMMON_H_

#include <Python.h>

/**
 * Common integer conversion for python 3 and 2.x.
 */
#if PY_MAJOR_VERSION < 3
    #define PyLong_FromLong(x) ( (PyObject *) PyInt_FromLong( (long) ( x ) ) )
    #define PyLong_AsLong(x) ( (long) PyInt_AsLong( ( x ) ) )
#endif

#endif
