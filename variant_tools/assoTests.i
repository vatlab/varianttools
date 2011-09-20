//!/usr/bin/env python2.7
//
// $File: assoTests.i $
// $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
// $Rev: 4234 $
//
// This file is part of variant_tools, a software application to annotate,
// summarize, and filter variants for next-gen sequencing ananlysis.
// Please visit http://variant_tools.sourceforge.net # for details.
//
// Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
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
#include "assoConfig.h"
#include "assoData.h"
#include "action.h"
#include "permutator.h"
#include "assoTests.h"
%}



%init
%{
    vtools::initialize();
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


%newobject vtools::AssoData::clone;
%newobject vtools::BaseAction::clone;
%newobject vtools::AssoTest::clone;

%include "std_vector.i"
%include "std_string.i"
%include "std_map.i"

namespace std
{
    /* %template(vectors)    vector<string>; */
    %template(vectorf)    vector<double>; 
    %template(vectori)    vector<int>; 
    %template(matrixi)    vector<vector<int> >;
    %template(vectora)    vector<vtools::BaseAction * >; 
}

%include "assoConfig.h"
%include "assoData.h"
%include "action.h"
%include "permutator.h"
%include "assoTests.h"
