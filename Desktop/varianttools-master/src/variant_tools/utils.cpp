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
#include <cassert>
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


//!- Mann-Whitney rank test statistic "U"
// www.alglib.net
// http://en.wikipedia.org/wiki/Mann-Whitney_U
double Mann_Whitneyu(const std::vector<double> & x, const std::vector<double> & y)
{
	int i;
	int k;
	int t;
	double tmp;
	int tmpi;
	int ns;
	double u;
	int w;
	int n = (int)x.size();
	int m = (int)y.size();

	// Prepare

	if (n < 5 || m < 5) {
		std::clog << "Sample size too small" << std::endl;
		return -9;
	}
	ns = n + m;
	std::vector<double> r(ns);
	std::vector<int> c(ns);
	for (i = 0; i <= n - 1; i++) {
		r[i] = x[i];
		c[i] = 0;
	}
	for (i = 0; i <= m - 1; i++) {
		r[n + i] = y[i];
		c[n + i] = 1;
	}

	// sort {R, C}, QS: smaller scores ranks higher

	if (ns != 1) {
		i = 2;
		do {
			t = i;
			while (t != 1) {
				k = t / 2;
				if (r[k - 1] >= r[t - 1]) {
					t = 1;
				}else {
					tmp = r[k - 1];
					r[k - 1] = r[t - 1];
					r[t - 1] = tmp;
					tmpi = c[k - 1];
					c[k - 1] = c[t - 1];
					c[t - 1] = tmpi;
					t = k;
				}
			}
			i = i + 1;
		} while (i <= ns);
		i = ns - 1;
		do {
			tmp = r[i];
			r[i] = r[0];
			r[0] = tmp;
			tmpi = c[i];
			c[i] = c[0];
			c[0] = tmpi;
			t = 1;
			while (t != 0) {
				k = 2 * t;
				if (k > i) {
					t = 0;
				}else {
					if (k < i) {
						if (r[k] > r[k - 1]) {
							k = k + 1;
						}
					}
					if (r[t - 1] >= r[k - 1]) {
						t = 0;
					}else {
						tmp = r[k - 1];
						r[k - 1] = r[t - 1];
						r[t - 1] = tmp;
						tmpi = c[k - 1];
						c[k - 1] = c[t - 1];
						c[t - 1] = tmpi;
						t = k;
					}
				}
			}
			i = i - 1;
		} while (i >= 1);
	}

	// Compute U

	u = 0;
	w = 1;
	for (i = 0; i <= ns - 1; i++) {
		if (i == 0)
			w = 1;
		else {
			if (r[i] > r[i - 1])
				w = i + 1;
		}

		if (c[i] == 0) {
			//ranks (sum of) for cases
			//std::cout << w << " ";
			//std::cout << r[i] << " ";
			u = u + w;
		}
	}
	// comment these out since I am using the raw scores
	// rather than the actual U statistic
	// in which case ties not properly handled
	// this is based on communication with Dr. Browning (UW)
	// and Dr. YeeHim Cheung (Columbia)
	//u = n*m+m*(m+1)/2-u;
	//std::cout << u << " ";
	return u;
}


double chisq2X2stat(const std::vector<double> & regressors, const std::vector<double> & responses)
{
	//! - 2 by 2 Chisq test
	double A0 = 0.0, A1 = 0.0, U0 = 0.0, U1 = 0.0;

	for (size_t i = 0; i < regressors.size(); ++i) {
		if (fEqual(regressors[i], 0.0) && fEqual(responses[i], 0.0))
			U0 += 1.0;
		else if (fEqual(regressors[i], 0.0) && fEqual(responses[i], 1.0))
			A0 += 1.0;
		else if (fEqual(regressors[i], 1.0) && fEqual(responses[i], 0.0))
			U1 += 1.0;
		else if (fEqual(regressors[i], 1.0) && fEqual(responses[i], 1.0))
			A1 += 1.0;
		else {
			return -9;
		}
	}     // collect the contigency table

	double Aobs = A0 + A1;
	double Uobs = U0 + U1;
	double Tobs = Aobs + Uobs;
	double Obs0 = A0 + U0;
	double Obs1 = A1 + U1;

	double EA0 = (Aobs * Obs0) / Tobs; if (EA0 == 0) EA0 = 0.05;
	double EA1 = (Aobs * Obs1) / Tobs; if (EA1 == 0) EA1 = 0.05;
	double EU0 = (Uobs * Obs0) / Tobs; if (EU0 == 0) EU0 = 0.05;
	double EU1 = (Uobs * Obs1) / Tobs; if (EU1 == 0) EU1 = 0.05;

	double statistic = ( (A0 - EA0) * (A0 - EA0) ) / EA0
	                   + ( (A1 - EA1) * (A1 - EA1) ) / EA1
	                   + ( (U0 - EU0) * (U0 - EU0) ) / EU0
	                   + ( (U1 - EU1) * (U1 - EU1) ) / EU1;
	return statistic;
}


std::vector<double> dnhyper(double m, double n, double k, double ncp)
{
	double lo = std::max(0.0, k - n);
	double hi = std::min(k, m);

	std::vector<double> d(0);
	for (int i = lo; i <= hi; ++i) {
		d.push_back(log(gsl_ran_hypergeometric_pdf(i, m, n, k)) + log(ncp) * i);
	}
	double dm = *max_element(d.begin(), d.end());
	for (size_t i = 0; i < d.size(); ++i) {
		d[i] = exp(d[i] - dm);
	}
	std::transform(d.begin(), d.end(), d.begin(), std::bind2nd(std::divides<double>(), std::accumulate(d.begin(), d.end(), 0.0)));
	return d;
}


double mnhyper(double m, double n, double k, double ncp)
{
	double lo = std::max(0.0, k - n);
	double hi = std::min(k, m);

	if (ncp == 0.0) {
		return lo;
	}
	if (std::abs(ncp) >= DBL_MAX) {
		return hi;
	}
	std::vector<double> d = dnhyper(m, n, k, ncp);
	double s = 0.0;
	for (size_t i = 0; i < d.size(); ++i) {
		s += (lo + (double)i) * d[i];
	}
	return s;
}


double pnhyper(double q, double m, double n, double k, double ncp, bool upper_tail)
{
	double lo = std::max(0.0, k - n);
	double hi = std::min(k, m);

	if (ncp == 1.0) {
		if (upper_tail) {
			return gsl_cdf_hypergeometric_Q(q - 1, m, n, k);
		} else {
			return gsl_cdf_hypergeometric_P(q, m, n, k);
		}
	}
	if (ncp == 0.0) {
		if (upper_tail) {
			return (double)(q <= lo);
		} else {
			return (double)(q >= lo);
		}
	}
	if (std::abs(ncp) >= DBL_MAX) {
		if (upper_tail) {
			return (double)(q <= hi);
		} else {
			return (double)(q >= hi);
		}
	}
	double p = 0.0;
	std::vector<double> d = dnhyper(m, n, k, ncp);
	if (upper_tail) {
		for (size_t i = 0; i < d.size(); ++i) {
			if (lo + (double)i >= q) p += d[i];
		}
	} else {
		for (size_t i = 0; i < d.size(); ++i) {
			if (lo + (double)i <= q) p += d[i];
		}
	}
	return p;
}


double fexact2x2(std::vector<int> dat, std::string alternative, double ncp)
{
	/*
	   x[0] x[1]
	   x[2] x[3]
	 */
	double m = dat[0] + dat[2];
	double n = dat[1] + dat[3];
	double k = dat[0] + dat[1];
	double x = dat[0];
	double lo = std::max(0.0, k - n);
	double hi = std::min(k, m);
	double s = 0.0;

	if (alternative == "greater") {
		s = pnhyper(x, m, n, k, ncp, true);
	} else if (alternative == "less") {
		s = pnhyper(x, m, n, k, ncp, false);
	} else {
		if (ncp == 0) {
			return (double)(x == lo);
		} else if (std::abs(ncp) >= DBL_MAX) {
			return (double)(x == hi);
		} else {
			double relErr = 1.0 + 1E-7;
			std::vector<double> d = dnhyper(m, n, k, ncp);
			for (size_t i = 0; i < d.size(); ++i) {
				if (d[i] <= d[(size_t)(x - lo)] * relErr)
					s += d[i];
			}
		}
	}
	return s;
}


double calculateInbreedingCoef(const std::vector<int> & gt,
                               const std::vector<double> & maf)
{
	//P(Homo) = F + (1-F)P(Homo by chance)
	//P(Homo by chance) = p^2+q^2 for a biallelic locus.
	//For an individual with N genotyped loci, we
	//  1. count the total observed number of loci which are homozygous (O),
	//  2. calculate the total expected number of loci homozygous by chance (E)
	//Then, using the method of moments, we have
	//   O = NF + (1-F)E
	//Which rearranges to give
	//   F = (O-E)/(N-E)
	assert(gt.size() == maf.size());
	double O = 0.0, E = 0.0, N = 0.0;
	for (unsigned i = 0; i < maf.size(); ++i) {
		// skip missing sites, tri-allilic sites or uninformative marker
		if (maf[i] <= 1E-8 || gt[i] != gt[i] || gt[i] < 0)
			continue;
		// observed homozygous locus found
		if (gt[i] == 0 || gt[i] == 2) O++;
		// expected value under HWE
		E += 1 - 2 * maf[i] * (1 - maf[i]);
        N++;
	}
	// Finally compute F
	return (N - E > 0) ? (O - E) / (N - E) : std::numeric_limits<double>::quiet_NaN();
}


namespace vtools {

swig_type_info * g_assoDataType = NULL;

PyObject * pyAssoDataObj(void * p)
{
	return SWIG_NewPointerObj(p, g_assoDataType, 0);
}


std::string PyObj_AsString(PyObject * str)
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
	std::string res(newstr);
	free(newstr);
	return res;
#else
	return std::string(PyString_AsString(str));
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
