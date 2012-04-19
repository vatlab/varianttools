/*
 *  $File: lm.cpp $
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

#include "lm.h"
namespace vtools {

bool LMData::setX(const std::vector<std::vector<double> > & x)
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


bool LMData::setY(const std::vector<double> & y)
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


std::vector<std::vector<double> > LMData::getX()
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


bool LMData::replaceColumn(const std::vector<double> & col, int which)
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


bool LinearM::fit(LMData & d)
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


bool LinearM::evalSE(LMData & d)
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


}
