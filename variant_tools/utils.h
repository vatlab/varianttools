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
#include <functional>
#include <vector>
#include <iostream>
#include <ctime>
#include <unistd.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_errno.h"

#include "assoTests.h"
// check if float numbers are equal
bool fEqual(double a, double b);

// round() float numbers
void fRound(double & myValue, double PRECISION);

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

namespace std {
// order a vector by index specified in another vector
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

// initialize C++ module, currently does nothing
void initialize();

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
		__seed = static_cast<unsigned long>(time(NULL) + getpid());
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


// linear model class
// members are
// m_ncol, m_nrow: dimensions
// m_x: predictor matrix
// m_y: respond vector
class BaseLm
{
public:
	BaseLm() : m_ncol(0), m_nrow(0), m_x(NULL), m_y(NULL)
	{
		gsl_set_error_handler_off();
	}


	~BaseLm()
	{
		if (m_x) {
			gsl_matrix_free(m_x);
		}
		if (m_y) {
			gsl_vector_free(m_y);
		}
	}


	BaseLm(const BaseLm & rhs) : m_ncol(rhs.m_ncol), m_nrow(rhs.m_nrow), m_x(NULL), m_y(NULL)
	{
		if (rhs.m_x) {
			m_x = gsl_matrix_alloc(m_nrow, m_ncol);
			gsl_matrix_memcpy(m_x, rhs.m_x);
		}
		if (rhs.m_y) {
			m_y = gsl_vector_alloc(m_nrow);
			gsl_vector_memcpy(m_y, rhs.m_y);
		}
	}


	virtual BaseLm * clone()
	{
		return new BaseLm(*this);
	}


	// initialize X, the predictor matrix
	void setX(const std::vector<std::vector<double> > & x)
	{
		if (x.size() == 0) {
			throw RuntimeError("No input data");
		}

		if (m_nrow == 0) m_nrow = x.front().size();
		else {
			if (m_nrow != x.front().size()) {
				throw ValueError("Dimension not match with input X");
			}
		}
		m_ncol = x.size();

		m_x = gsl_matrix_alloc(m_nrow, m_ncol);
		for (size_t j = 0; j < m_ncol; j++) {
			for (size_t i = 0; i < m_nrow; i++) {
				gsl_matrix_set(m_x, i, j, x[j][i]);
			}
		}
	}


	virtual void clear()
	{
		m_ncol = 0;
		m_nrow = 0;
		if (m_y) gsl_vector_free(m_y);
		if (m_x) gsl_matrix_free(m_x);
	}


	// initialize Y, the respond vector
	void setY(const std::vector<double> & y)
	{
		if (y.size() == 0) {
			throw RuntimeError("No input data");
		}
		if (m_nrow == 0) m_nrow = y.size();
		else {
			if (m_nrow != y.size()) {
				throw ValueError("Dimension not match with input Y");
			}
		}

		m_y = gsl_vector_alloc(m_nrow);
		for (size_t i = 0; i < m_nrow; i++) {
			gsl_vector_set(m_y, i, y[i]);
		}
	}


	// return X
	std::vector<std::vector<double> > getX()
	{
		std::vector<std::vector<double> > xout(m_ncol);
		for (size_t j = 0; j < m_ncol; j++) {
			for (size_t i = 0; i < m_nrow; i++) {
				xout[j].push_back(gsl_matrix_get(m_x, i, j));
			}
		}
		return xout;
	}


	// replace certain column in X, or replace the entire Y (when "which==0")
	void replaceCol(const std::vector<double> & col, int which)
	{
		if (which < 0 || which >= m_ncol) {
			throw IndexError("Invalid column index");
		}
		if (which != 0) {
			// set m_x
			if (!m_x) {
				throw RuntimeError("m_x matrix not initialized");
			}
			// will never replace the 0th col since it is (1...1)'
			for (size_t i = 0; i < m_nrow; ++i) {
				gsl_matrix_set(m_x, i, which, col[i]);
			}
		} else {
			// set m_y
			if (!m_y) {
				throw RuntimeError("m_y vector not initialized");
			}
			for (size_t i = 0; i < m_nrow; ++i) {
				gsl_vector_set(m_y, i, col[i]);
			}
		}
	}


protected:
	gsl_matrix * m_x;
	gsl_vector * m_y;
	int m_nrow;
	int m_ncol;
};


// fit the linear model
// via SVD technique
// members are
// m_beta: fitted model parameters (LSE)
// m_svd_: intermediate data from the SVD procedure
class LinearM : public BaseLm
{
public:
	LinearM() : BaseLm(), m_err(0), m_beta(NULL),
		m_svdS(NULL), m_svdV(NULL), m_svdU(NULL)
	{
	}


	~LinearM()
	{
		if (m_beta) {
			gsl_vector_free(m_beta);
		}
		if (m_svdS) {
			gsl_vector_free(m_svdS);
		}
		if (m_svdV) {
			gsl_matrix_free(m_svdV);
		}
		if (m_svdU) {
			gsl_matrix_free(m_svdU);
		}
	}


	LinearM(const LinearM & rhs) : BaseLm(rhs), m_err(rhs.m_err), m_beta(NULL),
		m_svdS(NULL), m_svdV(NULL), m_svdU(NULL)
	{
		if (rhs.m_beta) {
			m_beta = gsl_vector_alloc(m_ncol);
			gsl_vector_memcpy(m_beta, rhs.m_beta);
		}
		if (rhs.m_svdS) {
			m_svdS = gsl_vector_alloc(m_ncol);
			gsl_vector_memcpy(m_svdS, rhs.m_svdS);
		}
		if (rhs.m_svdU) {
			m_svdU = gsl_matrix_alloc(m_ncol, m_ncol);
			gsl_matrix_memcpy(m_svdU, rhs.m_svdU);
		}
		if (rhs.m_svdV) {
			m_svdV = gsl_matrix_alloc(m_ncol, m_ncol);
			gsl_matrix_memcpy(m_svdV, rhs.m_svdV);
		}
	}


	BaseLm * clone()
	{
		return new LinearM(*this);
	}


	void fit()
	{
		//fit data with a fresh start
		if (m_beta) {
			gsl_vector_free(m_beta);
		}
		if (m_svdS) {
			gsl_vector_free(m_svdS);
		}
		if (m_svdV) {
			gsl_matrix_free(m_svdV);
		}
		if (m_svdU) {
			gsl_matrix_free(m_svdU);
		}

		m_beta = gsl_vector_alloc(m_ncol);
		//compute X'Y
		gsl_vector * b = gsl_vector_alloc(m_ncol);
		m_err = gsl_blas_dgemv(CblasTrans, 1.0, m_x, m_y, 0.0, b);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemv(CblasTrans, 1.0, m_x, m_y, 0.0, b)");
		}
		//compute X'X
		//here m_svdU = X'X
		m_svdU = gsl_matrix_alloc(m_ncol, m_ncol);
		m_err = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, m_x, m_x, 0.0, m_svdU);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, m_x, m_x, 0.0, A)");
		}
		//svd for X'X
		//On output the matrix m_svdU will actually be computed.
		m_svdS = gsl_vector_alloc(m_ncol);
		m_svdV = gsl_matrix_alloc(m_ncol, m_ncol);
		gsl_vector * work = gsl_vector_alloc(m_ncol);
		m_err = gsl_linalg_SV_decomp(m_svdU, m_svdV, m_svdS, work);
		if (m_err != 0) {
			throw ValueError("Error in gsl_linalg_SV_decomp(A, V, s, work)");
		}


		//solve system Ax=b where x is beta
		m_err = gsl_linalg_SV_solve(m_svdU, m_svdV, m_svdS, b, m_beta);
		if (m_err != 0) {
			throw ValueError("Error in gsl_linalg_SV_solve(A, V, s, b, m_beta)");
		}


		gsl_vector_free(b);
		gsl_vector_free(work);
	}


	// obtain estimates for coefficients
	std::vector<double> getBeta()
	{
		std::vector<double> beta(m_ncol);
		for (size_t i = 0; i < m_ncol; ++i) {
			beta[i] = gsl_vector_get(m_beta, i);
		}
		return beta;
	}


	// obtain standard error for the LSE
	std::vector<double> getSEBeta()
	{
		if (!m_beta) {
			throw ValueError("Error in getSEBeta(): meed to fit the model first");
		}
		// compute (X'X)^-1 = V(diag(1/s))U'
		// diagnal matrix
		gsl_vector * oneovers = gsl_vector_alloc(m_ncol);
		for (size_t i = 0; i < m_ncol; ++i) {
			gsl_vector_set(oneovers, i, 1.0 / gsl_vector_get(m_svdS, i));
		}
		gsl_matrix * D = gsl_matrix_alloc(m_ncol, m_ncol);
		gsl_vector_view tmp = gsl_matrix_diagonal(D);
		gsl_matrix_set_zero(D);
		gsl_vector_memcpy(&tmp.vector, oneovers);
		gsl_vector_free(oneovers);

		gsl_matrix * V = gsl_matrix_alloc(m_ncol, m_ncol);
		m_err = gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m_svdV, D, 0.0, V);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m_svdV, D, 0.0, V)");
		}
		m_err = gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, m_svdU, 0.0, D);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, V, m_svdU, 0.0, D)");
		}
		gsl_matrix_free(V);

		// compute mse = (Y-Xb)'(Y-Xb) / (n-p)
		gsl_vector * fitted = gsl_vector_alloc(m_nrow);
		m_err = gsl_blas_dgemv(CblasNoTrans, 1.0, m_x, m_beta, 0.0, fitted);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemv(CblasNoTrans, 1.0, m_x, m_beta, 0.0, fitted)");
		}
		double mse = 0.0;
		for (size_t i = 0; i < m_nrow; ++i) {
			mse += pow(gsl_vector_get(m_y, i) - gsl_vector_get(fitted, i), 2.0);
		}
		gsl_vector_free(fitted);
		mse = mse / ((m_nrow - m_ncol) * 1.0);

		// s(b) = mse(X'X)^-1
		std::vector<double> seb(m_ncol);
		//tmp is the diagnal vector for D
		for (size_t i = 0; i < m_ncol; ++i) {
			seb[i] = sqrt(mse * gsl_vector_get(&tmp.vector, i));
		}
		gsl_matrix_free(D);
		return seb;
	}


private:
	gsl_vector * m_beta;
	gsl_vector * m_svdS;
	gsl_matrix * m_svdV;
	gsl_matrix * m_svdU;
	int m_err;
};

}
#endif
