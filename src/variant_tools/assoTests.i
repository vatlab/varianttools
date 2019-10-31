//!/usr/bin/env python2.7
//
// $File: assoTests.i $
// $LastChangedDate$
// $Rev$
//
// This file is part of variant_tools, a software application to annotate,
// summarize, and filter variants for next-gen sequencing ananlysis.
// Please visit https://github.com/vatlab/varianttools for details.
//
// Copyright (C) 2011 Gao Wang (wangow@gmail.com) and Bo Peng (bpeng@mdanderson.org)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//


%module assoTests

%{
#include "utils.h"
#include "lm.h"
#include "assoData.h"
#include "action.h"
#include "assoTests.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_cdf.h"
%}

// gsl functions initialization
%inline %{
void my_error_handler(const char * reason, const char * file,
                       int line, int gsl_errno)
{
    fprintf(stderr, "GSL Error %d:\t%s", gsl_errno, reason);
}
int gsl_initialize()
{
    gsl_set_error_handler(&my_error_handler);
    return 1;
}
%} 

%init
%{
    vtools::initialize();
    gsl_initialize();
%}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(vtools::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(vtools::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(vtools::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(vtools::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%newobject *::clone;

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

namespace std
{
    %template(vectors)    vector<string>;
    %template(vectorf)    vector<double>; 
    %template(vectori)    vector<int>; 
    %template(matrixi)    vector<vector<int> >;
    %template(matrixf)    vector<vector<double> >;
    %template(vectora)    vector<vtools::BaseAction * >; 
}

%ignore vtools::PyAction::PyAction(const PyAction & rhs);
%ignore vtools::PyFunc;

%include "utils.h"
%include "lm.h"
%include "assoData.h"
%include "action.h"

%ignore vtools::Exception;
%ignore vtools::IndexError;
%ignore vtools::ValueError;
%ignore vtools::SystemError;
%ignore vtools::RuntimeError;

%include "assoTests.h"

// gsl functions
extern double gsl_cdf_gaussian_P(double x, double sigma); 
extern double gsl_cdf_gaussian_Q(double x, double sigma); 
extern double gsl_cdf_gaussian_Pinv(double P, double sigma); 
extern double gsl_cdf_gaussian_Qinv(double Q, double sigma); 
extern double gsl_cdf_ugaussian_P(double x); 
extern double gsl_cdf_ugaussian_Q(double x); 
extern double gsl_cdf_ugaussian_Pinv(double P); 
extern double gsl_cdf_ugaussian_Qinv(double Q); 
