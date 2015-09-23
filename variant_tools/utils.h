/*
 *  $File: utils.h $
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
#ifndef _UTILS_H
#define _UTILS_H
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <functional>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <ctime>
#include <cfloat>
//#include <unistd.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "Python.h"
#include "assoTests.h"

// check if float numbers are equal
bool fEqual(double a, double b);

// round() float numbers
void fRound(double & myValue, double PRECISION);

// mannwhitney test
double Mann_Whitneyu(const std::vector<double> & x, const std::vector<double> & y);

// 2X2 table chisq statistic
double chisq2X2stat(const std::vector<double> & regressors, const std::vector<double> & responses);

// 2X2 table Fisher's test, translated from R's fisher.test
std::vector<double> dnhyper(double m, double n, double k, double ncp);

double mnhyper(double m, double n, double k, double ncp);

double pnhyper(double q, double m, double n, double k, double ncp = 1.0, bool upper_tail = false);

double fexact2x2(std::vector<int> dat, std::string alternative = "greater", double ncp = 1.0);

double calculateInbreedingCoef(const std::vector<int> & gt,
                               const std::vector<double> & maf);


// the "online variance" algorithm for calculating mean and variance
class RunningStat
{
public:
	RunningStat() : m_n(0) {}

	void Clear()
	{
		m_n = 0;
	}


	void Push(double x)
	{
		m_n++;

		// See Knuth TAOCP vol 2, 3rd edition, page 232
		if (m_n == 1) {
			m_oldM = m_newM = x;
			m_oldS = 0.0;
		}else {
			m_newM = m_oldM + (x - m_oldM) / m_n;
			m_newS = m_oldS + (x - m_oldM) * (x - m_newM);

			// set up for next iteration
			m_oldM = m_newM;
			m_oldS = m_newS;
		}
	}


	int NumDataValues() const
	{
		return m_n;
	}


	double Mean() const
	{
		return (m_n > 0) ? m_newM : 0.0;
	}


	double Variance() const
	{
		return ( (m_n > 1) ? m_newS / (m_n - 1) : 0.0);
	}


	double StandardDeviation() const
	{
		return sqrt(Variance() );
	}


private:
	int m_n;
	double m_oldM, m_newM, m_oldS, m_newS;
};


// pairwise sum of vector elements. can be used in std::accumulate for vector2f objects
struct VPlus
{
	template<typename T> std::vector<T> operator()(std::vector<T> x, std::vector<T> y)
	{
		//std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<T>());
		for (size_t i = 0; i < x.size(); ++i) {
			if (y[i] > 0) x[i] += y[i];
		}
		return x;
	}


};

template <typename T> std::string n2s(T Number)
{
	std::stringstream ss;

	ss << Number;
	return ss.str();
}


namespace std {
// order a vector by index specified in another vector
// Note that this is not exactly reordering, since the reorder vector contains the final location for each corresponding element in the data vector,
// which is different from a reorder vector that selects which initial data vector element should appear at each final location.
// this will not destruct the index vector
template< typename order_iterator, typename value_iterator >
void reorder(order_iterator order_begin, order_iterator order_end, value_iterator v)
{
	typedef typename iterator_traits< value_iterator >::value_type value_t;
	typedef typename iterator_traits< order_iterator >::value_type index_t;
	typedef typename iterator_traits< order_iterator >::difference_type diff_t;

	diff_t remaining = order_end - 1 - order_begin;
	for (index_t s = index_t(), d; remaining > 0; ++s) {
		for (d = order_begin[s]; d > s; d = order_begin[d]) ;
		if (d == s) {
			--remaining;
			value_t temp = v[s];
			while (d = order_begin[d], d != s) {
				swap(temp, v[d]);
				--remaining;
			}
			v[s] = temp;
		}
	}
}


// order a vector by index specified in another vector
// this is more efficient but will destruct the index vector
template< typename order_iterator, typename value_iterator >
void reorder_destructive(order_iterator order_begin, order_iterator order_end, value_iterator v)
{
	typedef typename iterator_traits< value_iterator >::value_type value_t;
	typedef typename iterator_traits< order_iterator >::value_type index_t;
	typedef typename iterator_traits< order_iterator >::difference_type diff_t;

	diff_t remaining = order_end - 1 - order_begin;
	for (index_t s = index_t(); remaining > 0; ++s) {
		index_t d = order_begin[s];
		if (d == (diff_t)-1) continue;
		--remaining;
		value_t temp = v[s];
		for (index_t d2; d != s; d = d2) {
			swap(temp, v[d]);
			swap(order_begin[d], d2 = (diff_t)-1);
			--remaining;
		}
		v[s] = temp;
	}
}


// stdout for vector
template<class T> ostream & operator<<(ostream & out, const vector<T> & vec)
{
	if (!vec.empty()) {
		typename vector<T>::const_iterator it = vec.begin();
		out << *it;
		for (++it; it != vec.end(); ++it)
			out << " " << *it ;
	}
	return out;
}


}

namespace vtools {

// random number generator from gsl
class RNG
{
public:
	// gsl_rng_alloc(gsl_rng_ranlxs2)
	RNG() : rng(gsl_rng_alloc(gsl_rng_mt19937) ) {};
	~RNG() { if (rng) gsl_rng_free(rng); }
	gsl_rng * get()
	{
		// time(NULL): number of seconds since 00:00:00 GMT Jan. 1, 1970
		__seed = static_cast<unsigned long>(time(NULL));
		gsl_rng_set(rng, __seed);
		return rng;
	}


	gsl_rng * get(const unsigned long seed)
	{
		__seed = seed;
		gsl_rng_set(rng, __seed);
		return rng;
	}


private:
	gsl_rng * rng;
	unsigned long __seed;
};

/** A wrapper to a python function
 *  CPPONLY
 */
class PyFunc
{
public:
	PyFunc(PyObject * func);

	~PyFunc()
	{
		Py_XDECREF(m_func);
	}


	PyFunc(const PyFunc & rhs) : m_func(rhs.m_func), m_name(rhs.m_name),
		m_numArgs(rhs.m_numArgs), m_args(rhs.m_args), m_flags(rhs.m_flags)
	{
		Py_XINCREF(m_func);
	}


	/// return number of arguments this function accepts.
	/// This function does not count tuple parameters.
	size_t numArgs() const
	{
		return m_numArgs;
	}


	std::string name() const
	{
		return m_name;
	}


	std::string arg(size_t arg) const
	{
		return m_args[arg];
	}


	bool isValid() const
	{
		return m_func != NULL;
	}


	PyObject * func() const
	{
		return m_func;
	}


	PyObject * operator()(const char * format, ...) const
	{
		va_list argptr;

		va_start(argptr, format);
		PyObject * arglist = Py_VaBuildValue(const_cast<char *>(format), argptr);
		va_end(argptr);
		PyObject * pyResult = PyEval_CallObject(m_func, arglist);

		Py_XDECREF(arglist);
		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed\n");
		}
		return pyResult;
	}


	PyObject * operator()(PyObject * args) const
	{
		PyObject * pyResult = PyEval_CallObject(m_func, args);

		if (pyResult == NULL) {
			PyErr_Print();
			PyErr_Clear();
			throw ValueError("Function call failed\n");
		}
		return pyResult;
	}


private:
	PyObject * m_func;

	std::string m_name;

	size_t m_numArgs;

	std::vector<std::string> m_args;

	int m_flags;
};


// initialize C++ module, currently does nothing
void initialize();

PyObject * pyAssoDataObj(void * p);

}
#endif
