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
#include "gsl/gsl_permutation.h"
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
	bool setX(const std::vector<std::vector<double> > & x);

	// initialize Y, the respond vector
	bool setY(const std::vector<double> & y);

	// return X
	std::vector<std::vector<double> > getX();

	// replace certain column in X, or replace the entire Y (when "which==0")
	bool replaceColumn(const std::vector<double> & col, size_t which);

	// set estimates for coefficients
	bool setBeta(gsl_vector * rhs)
	{
		m_beta.resize(rhs->size);
		for (size_t i = 0; i < rhs->size; ++i) {
			m_beta[i] = gsl_vector_get(rhs, i);
		}
		return true;
	}


	// set standard error
	bool setSEBeta(gsl_vector * rhs)
	{
		m_sebeta.resize(rhs->size);
		for (size_t i = 0; i < rhs->size; ++i) {
			m_sebeta[i] = gsl_vector_get(rhs, i);
		}
		return true;
	}


	std::vector<double> getBeta() { return m_beta; }

	std::vector<double> getSEBeta() { return m_sebeta; }

	gsl_matrix * x() { return m_x; }

	gsl_vector * y() { return m_y; }

private:
	gsl_matrix * m_x;
	gsl_vector * m_y;
	std::vector<double> m_beta;
	std::vector<double> m_sebeta;

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
		if (m_beta) {
			gsl_vector_free(m_beta);
		}
		if (m_sebeta) {
			gsl_vector_free(m_sebeta);
		}
	}


	BaseLM(const BaseLM & rhs) : m_err(rhs.m_err),
		m_beta(NULL), m_sebeta(NULL)
	{
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


//!- multiple linear regression parameter estimate
//!- BETA= (X'X)^{-1}X'Y => (X'X)BETA = X'Y
//!- Solve the system via gsl_linalg_SV_solve()
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


	// fit linear model via LSE
	bool fit(LMData & d);

	// obtain standard error for the LSE
	bool evalSE(LMData & d);

private:
	gsl_vector * m_svdS;
	gsl_matrix * m_svdV;
	gsl_matrix * m_svdU;
};


// fit the binary logistic regression model
// via gsl
// members are
// m_HI: Hessian matrix
// m_pi: fitted probability
// m_iterations: number of iterations used in Newton method
class LogisticM : public BaseLM
{
public:
	LogisticM() : BaseLM(), m_iterations(0),
		m_HI(NULL), m_pi(NULL)
	{
	}


	~LogisticM()
	{
		if (m_pi) {
			gsl_vector_free(m_pi);
		}
		if (m_HI) {
			gsl_matrix_free(m_HI);
		}
	}


	LogisticM(const LogisticM & rhs) : BaseLM(rhs),
		m_iterations(rhs.m_iterations), m_HI(NULL), m_pi(NULL)
	{
		if (rhs.m_pi) {
			m_pi = gsl_vector_alloc(rhs.m_pi->size);
			gsl_vector_memcpy(m_pi, rhs.m_pi);
		}
		if (rhs.m_HI) {
			m_HI = gsl_matrix_alloc(rhs.m_HI->size1, rhs.m_HI->size2);
			gsl_matrix_memcpy(m_HI, rhs.m_HI);
		}
	}


	BaseLM * clone() const
	{
		return new LogisticM(*this);
	}


	// fit logisitic regression model via MLE
	bool fit(LMData & d);

	// obtain standard error for the MLE
	bool evalSE(LMData & d);

	// number of iterations
	bool niterations() { return m_iterations; }

private:
	unsigned m_iterations;
	gsl_matrix * m_HI;
	gsl_vector * m_pi;
};

}
#endif
