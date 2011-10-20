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
#include "assoConfig.h"
#include "utils.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

namespace vtools {


class AssoData
{
public:
	AssoData() : m_phenotype(0), m_genotype(0), m_maf(0), m_X(0)
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
		//FIXME: will get these info directly from the variant table
		//Or, alternatively, may want to get the freq info from external sources rather than sample-based
    m_maf = maf;
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
					pnovar *= m_genotype[i][j];
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

	double simpleLinear(double xbar, double ybar)
  {
    // simple linear regression score test
    // FIXME: may later need a struct of output such as beta, CI, etc	
    //!- Statistic: LSE (MLE) for beta, centered and scaled (bcz E[b] = 0 and sigma = 1 by simulation) 
    //!- See page 41 of Kutner's Applied Linear Stat. Model, 5th ed.
    //
    double statistic = 0.0;
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
      statistic = (numerator / denominator) / sqrt(varb); 
    } 
    
    return statistic;
  }

	double gaussianP(double statistic, int sided = 1)
  {
    double p = 1.0;
    if (sided == 1) {
      p = gsl_cdf_ugaussian_Q(statistic);
    }
    else {
      p = gsl_cdf_chisq_Q(statistic*statistic, 1.0);
    }
    return p;
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
};

}
#endif
