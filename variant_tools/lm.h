/*
 *  $File: lm.h $
 *  $LastChangedDate: 2012-04-18 11:18:28 -0500 (Wed, 18 Apr 2012) $
 *  $Rev: 1093 $
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
#ifndef _LM_H
#define _LM_H
#include <vector>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_errno.h"

#include "assoTests.h"
#include "utils.h"

namespace vtools {
// linear model data class
// members are
// m_x: predictor matrix
// m_y: respond vector
// m_beta: fitted model parameters 
// m_sebeta: standard error for m_beta 
class LMData
{
public:
	LMData() : m_x(NULL), m_y(NULL),
		m_beta(0), m_sebeta(0)
	{
	}


	~LMData()
	{
		if (m_x) {
			gsl_matrix_free(m_x);
		}
		if (m_y) {
			gsl_vector_free(m_y);
		}
	}


	LMData(const LMData & rhs) : m_x(NULL), m_y(NULL),
		m_beta(rhs.m_beta), m_sebeta(rhs.m_sebeta)
	{
		if (rhs.m_x) {
			m_x = gsl_matrix_alloc(rhs.m_x->size1, rhs.m_x->size2);
			gsl_matrix_memcpy(m_x, rhs.m_x);
		}
		if (rhs.m_y) {
			m_y = gsl_vector_alloc(rhs.m_y->size);
			gsl_vector_memcpy(m_y, rhs.m_y);
		}
	}


	LMData * clone() const
	{
		return new LMData(*this);
	}


	void clear()
	{
		if (m_y) gsl_vector_free(m_y);
		if (m_x) gsl_matrix_free(m_x);
	}


	// initialize X, the predictor matrix
	bool setX(const std::vector<std::vector<double> > & x)
	{
		if (x.size() == 0) {
			throw RuntimeError("No input data");
		}

		if (m_y) {
			if (m_y->size != x.front().size()) {
				throw ValueError("Dimension not match with input X");
			}
		}

		m_x = gsl_matrix_alloc(x.front().size(), x.size());
		for (size_t j = 0; j < m_x->size2; j++) {
			for (size_t i = 0; i < m_x->size1; i++) {
				gsl_matrix_set(m_x, i, j, x[j][i]);
			}
		}
		return true;
	}


	// initialize Y, the respond vector
	bool setY(const std::vector<double> & y)
	{
		if (y.size() == 0) {
			throw RuntimeError("No input data");
		}
		if (m_x) {
			if (m_x->size1 != y.size()) {
				throw ValueError("Dimension not match with input Y");
			}
		}

		m_y = gsl_vector_alloc(y.size());
		for (size_t i = 0; i < m_y->size; i++) {
			gsl_vector_set(m_y, i, y[i]);
		}
		return true;
	}


	// return X
	std::vector<std::vector<double> > getX()
	{
		if (!m_x) {
			throw RuntimeError("m_x matrix not initialized");
		}
		std::vector<std::vector<double> > xout(m_x->size2);
		for (size_t j = 0; j < m_x->size2; j++) {
			for (size_t i = 0; i < m_x->size1; i++) {
				xout[j].push_back(gsl_matrix_get(m_x, i, j));
			}
		}
		return xout;
	}


	// replace certain column in X, or replace the entire Y (when "which==0")
	bool replaceColumn(const std::vector<double> & col, int which)
	{

		if (which != 0) {
			// set m_x
			if (!m_x) {
				throw RuntimeError("m_x matrix not initialized");
			}
			if (which < 0 || which >= m_x->size2) {
				throw IndexError("Invalid column index");
			}
			// will never replace the 0th col since it is (1...1)'
			for (size_t i = 0; i < m_x->size1; ++i) {
				gsl_matrix_set(m_x, i, which, col[i]);
			}
		} else {
			// set m_y
			if (!m_y) {
				throw RuntimeError("m_y vector not initialized");
			}
			for (size_t i = 0; i < m_y->size; ++i) {
				gsl_vector_set(m_y, i, col[i]);
			}
		}
		return true;
	}


	// obtain estimates for coefficients
	std::vector<double> getBeta() { return m_beta; }

	// obtain standard error
	std::vector<double> getSEBeta() { return m_sebeta; }


	bool setBeta(gsl_vector * rhs)
    {
        m_beta.resize(rhs->size);
        for (size_t i = 0; i < rhs->size; ++i) {
            m_beta[i] = gsl_vector_get(rhs, i);
        }
        return true;
    }


	bool setSEBeta(gsl_vector * rhs)
	{      
        m_sebeta.resize(rhs->size);
        for (size_t i = 0; i < rhs->size; ++i) {
            m_sebeta[i] = gsl_vector_get(rhs, i);
        }
		return true;
	}


private:
	gsl_matrix * m_x;
	gsl_vector * m_y;
    std::vector<double> m_beta;
	std::vector<double> m_sebeta;

public:
	gsl_matrix * x() { return m_x; }
	
    gsl_vector * y() { return m_y; }

};


class BaseLM
{
public:
	BaseLM() : m_err(0),
		m_beta(NULL), m_sebeta(NULL)
	{
		gsl_set_error_handler_off();
	}


	virtual ~BaseLM()
	{
	}


	BaseLM(const BaseLM & rhs) : m_err(rhs.m_err),
		m_beta(NULL), m_sebeta(NULL){
		if (rhs.m_beta) {
			m_beta = gsl_vector_alloc(rhs.m_beta->size);
			gsl_vector_memcpy(m_beta, rhs.m_beta);
		}
		if (rhs.m_sebeta) {
			m_sebeta = gsl_vector_alloc(rhs.m_sebeta->size);
			gsl_vector_memcpy(m_sebeta, rhs.m_sebeta);
		}
	}

	virtual BaseLM * clone() const
	{
		return new BaseLM(*this);
	}


	virtual bool fit(LMData & d)
	{
		throw RuntimeError("fit() method from the base linear model class should not be called");
		return true;
	}


	virtual bool evalSE(LMData & d)
	{
		throw RuntimeError("evalSE() method from the base linear model class should not be called");
		return true;
	}


protected:
	int m_err;
	gsl_vector * m_beta;
	gsl_vector * m_sebeta;
};


// fit the linear model
// via gsl
// members are
// m_svd_: intermediate data from the SVD procedure
class LinearM : public BaseLM
{
public:
	LinearM() : BaseLM(),
		m_svdS(NULL), m_svdV(NULL), m_svdU(NULL)
	{
	}


	~LinearM()
	{
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


	LinearM(const LinearM & rhs) : BaseLM(rhs),
		m_svdS(NULL), m_svdV(NULL), m_svdU(NULL)
	{
		if (rhs.m_svdS) {
			m_svdS = gsl_vector_alloc(rhs.m_svdS->size);
			gsl_vector_memcpy(m_svdS, rhs.m_svdS);
		}
		if (rhs.m_svdU) {
			m_svdU = gsl_matrix_alloc(rhs.m_svdU->size1, rhs.m_svdU->size2);
			gsl_matrix_memcpy(m_svdU, rhs.m_svdU);
		}
		if (rhs.m_svdV) {
			m_svdV = gsl_matrix_alloc(rhs.m_svdV->size1, rhs.m_svdV->size2);
			gsl_matrix_memcpy(m_svdV, rhs.m_svdV);
		}
	}


	BaseLM * clone() const
	{
		return new LinearM(*this);
	}


	bool fit(LMData & d)
	{
		gsl_matrix * x = d.x();
		size_t ncol = x->size2;
		gsl_vector * y = d.y();

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

		m_beta = gsl_vector_alloc(ncol);
		//compute X'Y
		gsl_vector * b = gsl_vector_alloc(ncol);
		m_err = gsl_blas_dgemv(CblasTrans, 1.0, x, y, 0.0, b);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemv(CblasTrans, 1.0, x, y, 0.0, b)");
		}
		//compute X'X
		//here m_svdU = X'X
		m_svdU = gsl_matrix_alloc(ncol, ncol);
		m_err = gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, x, x, 0.0, m_svdU);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, x, x, 0.0, A)");
		}
		//svd for X'X
		//On output the matrix m_svdU will actually be computed.
		m_svdS = gsl_vector_alloc(ncol);
		m_svdV = gsl_matrix_alloc(ncol, ncol);
		gsl_vector * work = gsl_vector_alloc(ncol);
		m_err = gsl_linalg_SV_decomp(m_svdU, m_svdV, m_svdS, work);
		if (m_err != 0) {
			throw ValueError("Error in gsl_linalg_SV_decomp(A, V, s, work)");
		}


		//solve system Ax=b where x is beta
		m_err = gsl_linalg_SV_solve(m_svdU, m_svdV, m_svdS, b, m_beta);
		if (m_err != 0) {
			throw ValueError("Error in gsl_linalg_SV_solve(A, V, s, b, m_beta)");
		}

		d.setBeta(m_beta);
		gsl_vector_free(b);
		gsl_vector_free(work);
		return true;
	}


	// obtain standard error for the LSE
	bool evalSE(LMData & d)
	{
		if (!m_beta) {
			throw ValueError("Error in evalSE(): need to fit the model first");
		}

		gsl_matrix * x = d.x();
		size_t nrow = x->size1;
		size_t ncol = x->size2;
		gsl_vector * y = d.y();

		if (m_beta->size != ncol) {
			throw ValueError("Error in evalSE(): fitted beta does not match input data dimension");
		}

		if (m_sebeta) {
			gsl_vector_free(m_sebeta);
		}

		// compute (X'X)^-1 = V(diag(1/s))U'
		// diagnal matrix
		gsl_vector * oneovers = gsl_vector_alloc(ncol);
		for (size_t i = 0; i < ncol; ++i) {
			gsl_vector_set(oneovers, i, 1.0 / gsl_vector_get(m_svdS, i));
		}
		gsl_matrix * D = gsl_matrix_alloc(ncol, ncol);
		gsl_vector_view tmp = gsl_matrix_diagonal(D);
		gsl_matrix_set_zero(D);
		gsl_vector_memcpy(&tmp.vector, oneovers);
		gsl_vector_free(oneovers);

		gsl_matrix * V = gsl_matrix_alloc(ncol, ncol);
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
		gsl_vector * fitted = gsl_vector_alloc(nrow);
		m_err = gsl_blas_dgemv(CblasNoTrans, 1.0, x, m_beta, 0.0, fitted);
		if (m_err != 0) {
			throw ValueError("Error in gsl_blas_dgemv(CblasNoTrans, 1.0, x, m_beta, 0.0, fitted)");
		}
		double mse = 0.0;
		for (size_t i = 0; i < nrow; ++i) {
			mse += pow(gsl_vector_get(y, i) - gsl_vector_get(fitted, i), 2.0);
		}
		gsl_vector_free(fitted);
		mse = mse / ((nrow - ncol) * 1.0);

		// s(b) = mse(X'X)^-1
		m_sebeta = gsl_vector_alloc(ncol);
		//tmp is the diagnal vector for D
		for (size_t i = 0; i < ncol; ++i) {
			gsl_vector_set(m_sebeta, i, sqrt(mse * gsl_vector_get(&tmp.vector, i)));
		}

		d.setSEBeta(m_sebeta);
		gsl_matrix_free(D);
		return true;
	}


private:
	gsl_vector * m_svdS;
	gsl_matrix * m_svdV;
	gsl_matrix * m_svdU;
};

}
#endif
