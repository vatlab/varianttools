/*
 *  $File: assoData.h $
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

#ifndef _ASSODATA_H
#define _ASSODATA_H


#include <vector>
typedef std::vector<double> vectorf;
typedef std::vector<std::vector<double> > matrixf;
typedef std::vector<int> vectori;
typedef std::vector<std::vector<int> > matrixi;
#include <cassert>
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>
#include "assoConfig.h"
#include "utils.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"


namespace vtools {

  class AssoData
  {
    public:
      AssoData() : 
        m_phenotype(0), m_genotype(0), m_maf(0), m_X(0),
        m_statistic(0.0), m_pval(0.0), m_C(0), m_ncovar(0),
        m_sites(0), m_xbar(0.0), m_ybar(0.0), m_ncases(0),
        m_nctrls(0), m_model()
    {
    }


      ~AssoData(){}

      // make a copy of the data
      virtual AssoData * clone() const
      {
        return new AssoData(*this);
      }


      // set data
      double setGenotype(const matrixf & g)
      {
        //codings are 0, 1, 2, -9 and U(0,1) number for "expected" genotype
        m_genotype = g;
        return 0.0;
      }


      double setPhenotype(const vectorf & p)
      {
        //FIXME: have to consider missing phenotypes (especially when they are treated covariates) as well as multiple phenotypes
        m_phenotype = p;
        // set phenotype statistics
        m_ybar = std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0);
        m_ybar /= (1.0 * m_phenotype.size());
        m_ncases = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(), 
            std::bind2nd(std::equal_to<double>(),1.0));
        m_nctrls = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(), 
            std::bind2nd(std::equal_to<double>(),0.0));
        return m_ybar;
      }


      double setPhenotype(const vectorf & p, const matrixf & c)
      {
        m_phenotype = p;
        m_C = c;
        m_ncovar = c.size()-1;
        vectorf one(p.size());
        std::fill(one.begin(), one.end(), 1.0);
        m_C.push_back(one);
        m_model.clear();
        m_model.setY(m_phenotype);
        m_model.setX(m_C);
        // set phenotype statistics
        m_ybar = std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0);
        m_ybar /= (1.0 * m_phenotype.size());
        m_ncases = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(), 
            std::bind2nd(std::equal_to<double>(),1.0));
        m_nctrls = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(), 
            std::bind2nd(std::equal_to<double>(),0.0));
        return m_ybar;
      }


      void setMaf(const vectorf & maf)
      {
        //get this field directly from the variant table
        m_maf = maf;
      }

      void setMaf()
      {
        //sample based maf
        struct VPlus vplus;
        vectorf gx = m_genotype.front();
        for (size_t i = 0; i < gx.size(); ++i) {
          if (gx[i] < -1.0) {
            // missing
            gx[i] = 0.0;
            continue;
          } 
        }
        m_maf = std::accumulate(m_genotype.begin() + 1, m_genotype.end(), 
            gx, vplus);
        std::transform(m_maf.begin(), m_maf.end(), m_maf.begin(),
            std::bind2nd(std::divides<double>(), 2.0*m_genotype.size()));
      }

      void setMafWeight()
      {
        // compute w = 1 / sqrt(p*(1-p))
        if (m_maf.size()==0) {
          throw RuntimeError("MAF has not been calculated. Please calculate MAF prior to calculating weights.");
        }
        //
        m_weight.clear();
        for (size_t i = 0; i < m_maf.size(); ++i) {
          if (fEqual(m_maf[i], 0.0) || fEqual(m_maf[i], 1.0)) {
            m_weight.push_back(0.0);
          } else 
          {
            m_weight.push_back(1.0/sqrt(m_maf[i] * (1.0-m_maf[i])));
          }
        }
      }

      void setMafWeight(const std::vector<size_t> & idx) 
      {
        // FIXME weight by maf from selected samples (ctrls, low QT samples, etc)
      }


      double meanOfX()
      {
        m_xbar = std::accumulate(m_X.begin(), m_X.end(), 0.0);
        m_xbar /= (1.0 * m_X.size());
        return m_xbar;
      }


      // return some information
      vectorf phenotype()
      {
        return m_phenotype;
      }

      vectorf genotype()
      {
        return m_X;
      }

      matrixf raw_genotype()
      {
        return m_genotype;
      }

      matrixf covariates()
      {
        return m_model.getX();
      }

      vectorf maf()
      {
        return m_maf;
      }

      vectori sites()
      {
        return m_sites;
      }

      unsigned covarcounts()
      {
        return m_ncovar;
      }

      unsigned samplecounts()
      {
        return m_phenotype.size();
      }

      double pvalue()
      {
        return m_pval;
      }

      double setPvalue(double pval)
      {
        m_pval = pval;
        return 0.0;
      }      

      double setStatistic(double stat)
      {
        m_statistic = stat;
        return 0.0;
      }

      double statistic()
      {
        return m_statistic;
      }


    public:
      void permuteY()
      {
        random_shuffle(m_phenotype.begin(), m_phenotype.end());
      }

      void permuteX()
      {
        random_shuffle(m_X.begin(), m_X.end());
      }

      // manipulate data
      void sumToX()
      {
        m_X.resize(m_genotype.size());
        std::fill(m_X.begin(), m_X.end(), 0.0);
        for (size_t i = 0; i < m_genotype.size(); ++i) {
          //m_X[i] = std::accumulate(m_genotype[i].begin(), m_genotype[i].end(), 0.0);
          for (size_t j = 0; j < m_genotype[i].size(); ++j) {
            // check if we should skip this site
            if (m_sites.size() != 0) { 
              if (m_sites[j] == 0) continue;
            }

            if (m_genotype[i][j] > 0) m_X[i] += m_genotype[i][j]; 
          }
        }
      }

      void weightX()
      {
        if (m_weight.size() == 0) return;
        for (size_t i = 0; i < m_genotype.size(); ++i) {
          for (size_t j = 0; j < m_genotype[i].size(); ++j) {
            if (m_genotype[i][j] > 0) m_genotype[i][j] *= m_weight[j]; 
          }
        }
      }

      void binToX()
      {
        // binning the data with proper handling of missing genotype
        m_X.resize(m_genotype.size());
        std::fill(m_X.begin(), m_X.end(), 0.0);
        for (size_t i = 0; i < m_genotype.size(); ++i) {
          double pnovar = 1.0;
          for (size_t j = 0; j != m_genotype[i].size(); ++j) {  
            // check if we should skip this site
            if (m_sites.size() != 0) { 
              if (m_sites[j] == 0) continue;
            }

            if (m_genotype[i][j] >= 1.0) {  
              m_X[i] = 1.0;
              break;
            }
            else if (m_genotype[i][j] > 0.0) {
              pnovar *= (1.0 - m_genotype[i][j]);
            } 
            else;
          }
          if (pnovar < 1.0 && m_X[i] < 1.0) {
            m_X[i] = 1.0 - pnovar;
          }
        }
      }


      void setSitesByMaf(double upper, double lower)
      {
        //Will update m_sites with 0 = sites not to be analyzed; 1 = sites to be analyzed 
        m_sites.clear();
        if (upper > 1.0) {
          throw ValueError("Minor allele frequency value should not exceed 1");
        }
        if (lower < 0.0) {
          throw ValueError("Minor allele frequency should be a positive value");
        }

        if (fEqual(upper,1.0) && fEqual(lower,0.0)) return;

        for (size_t j = 0; j != m_maf.size(); ++j) {
          if (m_maf[j] <= lower || m_maf[j] > upper) { 
            m_sites.push_back(0);
          }
          else {
            m_sites.push_back(1);
          }
        }  
      }


      void simpleLinear()
      {
        // simple linear regression score test
        // FIXME: may later need other output fields such as beta, CI, etc	
        //!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation) 
        //!- See page 23 and 41 of Kutner's Applied Linear Stat. Model, 5th ed.
        //
        if (m_X.size() != m_phenotype.size()) {
          throw ValueError("Genotype/Phenotype length not equal!");
        }
        double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
        for (size_t i = 0; i != m_X.size(); ++i) {
          numerator += (m_X[i] - m_xbar) * m_phenotype[i];
          denominator += pow(m_X[i] - m_xbar, 2.0);
        }

        if (!fEqual(numerator, 0.0)) {  
          //!- Compute MSE and V[\hat{beta}]
          //!- V[\hat{beta}] = MSE / denominator
          double b1 = numerator / denominator;
          double b0 = m_ybar - b1*m_xbar;

          //SSE
          for (size_t i = 0; i != m_X.size(); ++i) {
            ysigma += pow(m_phenotype[i] - (b0+b1*m_X[i]), 2.0);
          }
          double varb = ysigma / (m_phenotype.size() - 2.0) / denominator;
          m_statistic = b1 / sqrt(varb); 
        }
        else m_statistic = 0.0;
      }


      void simpleLogit()
      {
        //!- Score test implementation for logistic regression model logit(p) = b0 + b1x 
        //!- labnotes vol.2 page 3
        //!- input phenotypes have to be binary values 0 or 1
        if (m_X.size() != m_phenotype.size()) {
          throw ValueError("Genotype/Phenotype length not equal!");
        }
        //double ebo = (1.0 * n1) / (1.0 * (m_phenotype.size()-n1));
        //double bo = log(ebo);

        double po = (1.0 * m_ncases) / (1.0 * m_phenotype.size());
        double ss = 0.0;
        // the score
        for (size_t i = 0; i != m_X.size(); ++i) { 
          ss += (m_X[i] - m_xbar) * (m_phenotype[i] - po);
        }
        double vm1 = 0.0;
        // variance of score, under the null
        for (size_t i = 0; i != m_X.size(); ++i) {
          vm1 += (m_X[i] - m_xbar) * (m_X[i] - m_xbar) * po * (1.0 - po);
        }

        m_statistic = ss / sqrt(vm1);

        //!-FIXME: (not sure why this happens)
        //!- w/ rounding to 0 I get strange number such as 3.72397e-35
        //!- this would lead to type I error problem
        fRound(m_statistic, 0.0001);
      }


      void multipleLinear()
      {
        //!- multiple linear regression parameter estimate
        //!- BETA= (X'X)^{-1}X'Y => (X'X)BETA = X'Y
        //!- Solve the system via gsl_linalg_SV_solve()
        if (m_X.size() != m_phenotype.size()) {
          throw ValueError("Genotype/Phenotype length not equal!");
        }
        m_model.replaceCol(m_X, m_C.size()-1);
        m_model.fit();
        vectorf beta = m_model.getBeta();
        vectorf seb = m_model.getSEBeta();
        m_statistic = beta[beta.size()-1] / seb[seb.size()-1]; 
      }


      void gaussianP(unsigned sided = 1)
      {
        if (sided == 1) {
          m_pval = gsl_cdf_ugaussian_Q(m_statistic);
        }
        else if (sided == 2) {
          m_pval = gsl_cdf_chisq_Q(m_statistic*m_statistic, 1.0);
        }
        else {
          throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
        }
      }

      void studentP(unsigned sided = 1)
      {
        // df = n - p where p = #covariates + 1 (for beta1) + 1 (for beta0) = m_ncovar+2
        if (sided == 1) {
          m_pval = gsl_cdf_tdist_Q(m_statistic, m_phenotype.size() - (m_ncovar+2.0));
        }
        else if (sided == 2) {
          double p = gsl_cdf_tdist_Q(m_statistic, m_phenotype.size() - (m_ncovar+2.0));
          m_pval = fmin(p, 1.0-p) * 2.0;
        }
        else {
          throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
        }
      }

      // permutation should use Clopper-Pearson 95% interval which is exact and conservative	

    private:
      /// raw phenotype and gneotype data
      vectorf m_phenotype;
      matrixf m_genotype;

      // covariates
      matrixf m_C;
      unsigned m_ncovar;

      // observed minor allele frequencies
      vectorf m_maf;
      /// translated genotype
      vectorf m_X;

      // sites that are involved in the association test
      vectori m_sites;

      // weights
      vectorf m_weight;

      double m_xbar;
      double m_ybar;
      unsigned m_ncases;
      unsigned m_nctrls;
      double m_pval;
      double m_statistic;
      LinearM m_model;
  };

}
#endif
