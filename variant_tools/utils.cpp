/*
 *  $File: utils.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of variant_tools, a software application to annotate,
 *  summarize, and filter variants for next-gen sequencing ananlysis.
 *  Please visit http://varianttools.sourceforge.net for details.
 *
 *  Copyright (C) 2011 Gao Wang (wangow@gmail.com) and Bo Peng (bpeng@mdanderson.org)
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#include "utils.h"

#include "swigpyrun.h"

bool fEqual(double a, double b)
{
	return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}


void fRound(double & myValue, double PRECISION)
{
	double myRemainder = fmod(myValue, PRECISION);

	if (myRemainder > (0.5 * PRECISION)) {
		myValue = (myValue - myRemainder + PRECISION);
	}else {
		myValue = (myValue - myRemainder);
	}
	return;
}


namespace vtools {

swig_type_info * g_assoDataType = NULL;

PyObject * pyAssoDataObj(void * p)
{
	return SWIG_NewPointerObj(p, g_assoDataType, 0);
}


string PyObj_AsString(PyObject * str)
{
#if PY_VERSION_HEX >= 0x03000000
	char * cstr;
	char * newstr;
	Py_ssize_t len;
	str = PyUnicode_AsUTF8String(str);
	PyBytes_AsStringAndSize(str, &cstr, &len);
	newstr = (char *)malloc(len + 1);
	memcpy(newstr, cstr, len + 1);
	Py_XDECREF(str);
	string res(newstr);
	free(newstr);
	return res;
#else
	return string(PyString_AsString(str));
#endif
}


PyFunc::PyFunc(PyObject * func) : m_func(func), m_numArgs(0)
{
	Py_XINCREF(m_func);
	if (m_func == NULL)
		return;

	if (!PyCallable_Check(m_func))
		throw ValueError("Passed parameter should be None or a Python function");

	if (PyObject_HasAttrString(m_func, "__call__")) {
		if (PyObject_HasAttrString(m_func, "__args__")) {
			// in this case, a WithArgs m_funcect must have been passed.
			PyObject * args = PyObject_GetAttrString(m_func, "__args__");
			m_numArgs = PySequence_Size(args);
			for (size_t i = 0; i < m_numArgs; ++i) {
				PyObject * item = PySequence_GetItem(args, i);
				if (!PyString_Check(item))
					throw ValueError("Attribute args in a simuPOP WithArgs m_funcect should be a list of strings");
				m_args.push_back(PyObj_AsString(item));
				Py_DECREF(item);
			}
			Py_DECREF(args);
			// find its name.
			PyObject * func = PyObject_GetAttrString(m_func, "__call__");
			if (!PyCallable_Check(func))
				throw ValueError("The func attribute of the passed m_funcect should be callable.");
			if (PyObject_HasAttrString(func, "__name__")) {
				PyObject * name = PyObject_GetAttrString(func, "__name__");
				m_name = PyObj_AsString(name);
				Py_DECREF(name);
			} else if (PyObject_HasAttrString(func, "__call__")) {
				PyObject * func1 = PyObject_GetAttrString(func, "__call__");
				if (PyObject_HasAttrString(func1, "__name__")) {
					PyObject * name = PyObject_GetAttrString(func1, "__name__");
					m_name = PyObj_AsString(name);
					Py_DECREF(name);
				}
				Py_DECREF(func1);
			}
			Py_DECREF(func);
			return;
#if PY_VERSION_HEX < 0x03000000
		} else if (!PyObject_HasAttrString(m_func, "func_code")) {
#else
		} else if (!PyObject_HasAttrString(m_func, "__code__")) {
#endif
			// if there is no arg so it is a member function of a class
			m_func = PyObject_GetAttrString(m_func, "__call__");
		}
	}
	// is it unbounded?
#if PY_VERSION_HEX < 0x03000000
	int bounded = PyObject_HasAttrString(m_func, "im_self");
#else
	int bounded = PyObject_HasAttrString(m_func, "__self__");
#endif

	if (!PyObject_HasAttrString(m_func, "__name__"))
		throw ValueError("Cannot find name of the passed function.");
	// find its name.
	PyObject * name = PyObject_GetAttrString(m_func, "__name__");
	m_name = PyObj_AsString(name);
	Py_DECREF(name);

	// free python functions have a 'func_code' attribute
	// built-in functions might not have (e.g. random.random)
#if PY_VERSION_HEX < 0x03000000
	if (!PyObject_HasAttrString(m_func, "func_code"))
		return;
	PyObject * code = PyObject_GetAttrString(m_func, "func_code");
#else
	if (!PyObject_HasAttrString(m_func, "__code__"))
		return;
	PyObject * code = PyObject_GetAttrString(m_func, "__code__");
#endif

	if (!code)
		throw SystemError("Invalid attribute func_code or __code for a function m_funcect");
	// probe number of parameters
	PyObject * co_argcount = PyObject_GetAttr(code, PyString_FromString("co_argcount"));
	if (!co_argcount)
		throw RuntimeError("Invalid attribute co_argcount for a function m_funcect");
	// substract 1 if the method is bounded to remove the count for self.
	m_numArgs = PyInt_AsLong(co_argcount) - bounded;
	Py_DECREF(co_argcount);
	// probe parameter names
	PyObject * co_varnames = PyObject_GetAttr(code, PyString_FromString("co_varnames"));
	if (!co_varnames)
		RuntimeError("Invalid attribute co_varnames for a function m_funcect");
	for (size_t i = 0; i < m_numArgs; ++i) {
		PyObject * item = PyTuple_GetItem(co_varnames, i + bounded);
		m_args.push_back(PyObj_AsString(item));
	}
	Py_DECREF(co_varnames);
	// accepting arbitrary number of parameters?
	PyObject * co_flag = PyObject_GetAttrString(code, "co_flags");
	if (!co_flag)
		throw RuntimeError("Invalid attribute co_flags for a function m_funcect");
	m_flags = static_cast<unsigned char>(PyInt_AsLong(co_flag));
	Py_DECREF(co_flag);
	Py_DECREF(code);
}


// currently initialization does not do anything
void initialize()
{
	g_assoDataType = SWIG_TypeQuery("vtools::AssoData *");
	if (!g_assoDataType)
		throw RuntimeError("Cannot determine AssoData type from the python interface.");
}


}

