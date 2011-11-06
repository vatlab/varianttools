/*
 *  $File: utils.h $
 *  $LastChangedDate: 2011-07-06 23:27:10 -0500 (Wed, 06 Jul 2011) $
 *  $Rev: 4256 $
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
#include <cassert>

struct VPlus {
  template<typename T> std::vector<T> operator()(std::vector<T> x, std::vector<T> y) {
    //std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<T>());
    for (size_t i = 0; i < x.size(); ++i) {
      if (y[i] > 0) x[i] += y[i];
    }
    return x;
  }
};

bool fEqual(double a, double b);
void fRound(double& myValue, double PRECISION);

class BaseLm
{
  public:
    BaseLm() : m_ncol(0), m_nrow(0) 
    {
    }

    ~BaseLm()
    {
      gsl_matrix_free(m_x);
      gsl_vector_free(m_y);
    }

    virtual BaseLm * clone()
    {
      return new BaseLm(*this);
    }

    void setX(const std::vector<std::vector<double> > &x)
    { 
      assert(x.size() > 0);
      if (m_nrow == 0) m_nrow = x[0].size();
      else assert(m_nrow == x[0].size());
      m_ncol = x.size();

      m_x = gsl_matrix_alloc(m_nrow, m_ncol);
      for (size_t j = 0; j < m_ncol; j++) {
        for (size_t i = 0; i < m_nrow; i++) {
          gsl_matrix_set(m_x, i, j, x[j][i]);         
        }
      }
    }

    void setY(std::vector<double> &y)
    { 
      assert(y.size() > 0);
      if (m_nrow == 0) m_nrow = y.size();
      else assert(m_nrow == y.size());

      m_y = gsl_vector_alloc(m_nrow); 
      for (size_t i = 0; i < m_nrow; i++) {
        gsl_vector_set(m_y, i, y[i]);
      }
    }


    void replaceCol(const std::vector<double> &col, int which)
    {
      assert(which > 0);
      // will never replace the 0th col since it is (1...1)'
      for (size_t i = 0; i < m_nrow; ++i) {
        gsl_matrix_set(m_x, i, which, col[i]); 
      }
    }


  private:
    gsl_matrix *m_x;  
    gsl_vector *m_y;  
    int m_nrow;
    int m_ncol;
};

class LinearM : public BaseLm
{
  public:
    LinearM() : BaseLm(), m_err(0)
    {
    }

    ~LinearM()
    {
      gsl_vector_free(m_beta);
    }

    BaseLm * clone()
    {
      return new LinearM(*this);
    }

    void fit()
    {
      //compute X'Y
      gsl_vector *b = gsl_vector_alloc(m_ncol);
      m_err = gsl_blas_dgemv(CblasTrans, 1.0, m_x, m_y, 0.0, b);
      assert(m_err==0);
      //compute X'X
      gsl_matrix *A = gsl_matrix_alloc(m_ncol, m_ncol);
      m_err = gsl_blas_dgemv(CblasTrans, CblasNoTrans, 1.0, m_x, m_x, 0.0, A);
      assert(m_err==0);
      //svd for X'X
      //On output the matrix A is replaced by U
      gsl_vector *s = gsl_vector_alloc(m_ncol);
      gsl_matrix *V = gsl_matrix_alloc(m_ncol, m_ncol);
      gsl_vector *work = gsl_vector_alloc(m_ncol);
      m_err = gsl_linalg_SV_decomp(A, V, s, work);
      assert(m_err==0);
      //solve system Ax=b where x is beta
      m_err = gsl_linalg_SV_solve (A, V, s, b, m_beta);
      assert(m_err==0);
      //
      gsl_matrix_free(A);
      gsl_matrix_free(v);
      gsl_vector_free(b);
      gsl_vector_free(s);
      gsl_vector_free(work);
    }

    std::vector<double> getBeta()
    {
      std::vector<double> beta(m_ncol);
      for (size_t i = 0; i < m_ncol; ++i) beta[i] = gsl_vector_get(m_beta, i);
      return beta;
    }

  private:
    gsl_vector *m_beta;
    int m_err;

};
#endif
