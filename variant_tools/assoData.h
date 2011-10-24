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
    m_statistic(0.0), m_pval(0.0)
	{

	}


	// make a copy of the data
	virtual AssoData * clone() const
	{
		return new AssoData(*this);
	}


	// set data
	void setGenotype(const matrixf & g)
	{
		//FIXME: codings are 0, 1, 2, -9 and float numbers for most likely genotypes
		m_genotype = g;
	}


	void setPhenotype(const vectorf & p)
	{
    //FIXME: have to consider missing phenotypes (especially when they are treated covariates) as well as multiple phenotypes
		m_phenotype = p;
	}

	void setMaf(const vectorf & maf)
	{
		//will get these info directly from the variant table
		//Or, alternatively, may want to get the freq info from external sources rather than sample-based
    m_maf = maf;
	}


	void setMaf()
	{
		//compute sample based maf
    m_maf = std::accumulate(m_genotype.begin() + 1, m_genotype.end(), m_genotype[0], vplus);
    std::transform(m_maf.begin(), m_maf.end(), m_maf.begin(),
std::bind2nd(std::divides<double>(), 2.0*m_genotype.size()));
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

	vectorf maf()
	{
		return m_maf;
	}

  double mean_phenotype()
  {
    double mean = std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0);
    return mean / (1.0 * m_phenotype.size());
  }

  double mean_genotype()
  {
    double mean = std::accumulate(m_X.begin(), m_X.end(), 0.0);
    return mean / (1.0 * m_X.size());
  }

  unsigned count_cases()
  {
    return (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(), std::bind2nd(std::equal_to<double>(),1.0));
  }

  double pvalue()
  {
    return m_pval;
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


	// manipulate data
	void sumToX()
	{
		m_X.resize(m_genotype.size());
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			m_X[i] = std::accumulate(m_genotype[i].begin(), m_genotype[i].end(), 0.0);
		}
	}

	void binToX()
	{
		// binning the data with proper handling of missing genotype
		m_X.resize(m_genotype.size());
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			m_X[i] = 0.0;
			double pnovar = 1.0;
			for (size_t j = 0; j != m_genotype[i].size(); ++j) {  
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


	void filterByMaf(double upper=0.01, double lower=0.0)
  {
    //FIXME: may want to do it smartly via sqlite. i.e., 
    // "select * from genotype_x where maf<upper and maf>lower"
    if (upper > 1.0 || lower < 0.0) {
      throw ValueError("Minor allele frequency value should fall between 0 and 1.");
    }
    if (fEqual(upper,1.0) && fEqual(lower,0.0)) return;

    matrixf xdat = m_genotype;
    m_genotype.clear();
    m_genotype.resize(xdat.size());

    for (size_t j = 0; j != m_maf.size(); ++j) {
      if (m_maf[j] <= lower || m_maf[j] > upper) 
        continue;
      else {
        for (size_t i = 0; i != xdat.size(); ++i)
          m_genotype[i].push_back(xdat[i][j]);
      }
    }  
  }


	void simpleLinear(double xbar, double ybar)
  {
    // simple linear regression score test
    // FIXME: may later need a struct of output such as beta, CI, etc	
    //!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation) 
    //!- See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
    //
    double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
    for (size_t i = 0; i != m_X.size(); ++i) {
      numerator += (m_X[i] - xbar) * m_phenotype[i];
      denominator += pow(m_X[i] - xbar, 2.0);
      //SSE
      ysigma += pow(m_phenotype[i] - ybar, 2.0);
    }

    if (!fEqual(numerator, 0.0)) {  
      //!- Compute MSE and V[\hat{beta}]
      //!- V[\hat{beta}] = MSE / denominator
      double varb = ysigma / (m_phenotype.size() * denominator);
      m_statistic = (numerator / denominator) / sqrt(varb); 
    }
    else m_statistic = 0.0;
  }


  void simpleLogit(double xbar, unsigned n1)
  {
    //!- Score test implementation for logistic regression model logit(p) = b0 + b1x 
    //!- labnotes vol.2 page 3
    //!- input phenotypes have to be binary values 0 or 1

    //double ebo = (1.0 * n1) / (1.0 * (m_phenotype.size()-n1));
    //double bo = log(ebo);

    double po = (1.0 * n1) / (1.0 * m_phenotype.size());
    double ss = 0.0;
    // the score
    for (size_t i = 0; i != m_X.size(); ++i) { 
      ss += (m_X[i] - xbar) * (m_phenotype[i] - po);
    }
    double vm1 = 0.0;
    // variance of score, under the null
    for (size_t i = 0; i != m_X.size(); ++i) {
      vm1 += (m_X[i] - xbar) * (m_X[i] - xbar) * po * (1.0 - po);
    }

    m_statistic = ss / sqrt(vm1);

    //!-FIXME: (not sure why this happens)
    //!- w/ rounding to 0 I get strange number such as 3.72397e-35
    //!- this would lead to type I error problem
    fRound(m_statistic, 0.0001);
  }


	void gaussianP(unsigned sided = 1)
  {
    if (sided == 1) {
      m_pval = gsl_cdf_ugaussian_Q(m_statistic);
    }
    else {
      m_pval = gsl_cdf_chisq_Q(m_statistic*m_statistic, 1.0);
    }
  }

// permutation should use Clopper-Pearson 95% interval which is exact and conservative	


private:
	/// raw phenotype and gneotype data
	vectorf m_phenotype;
	matrixf m_genotype;

  // observed minor allele frequencies
	vectorf m_maf;
	/// translated genotype
	vectorf m_X;

  double m_pval;
  double m_statistic;
};

}
#endif
