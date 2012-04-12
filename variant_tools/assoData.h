/*
 *  $File: assoData.h $
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
#include <limits>
#include "assoTests.h"
#include "utils.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"


namespace vtools {
/*  Association data object
 *  This class define memembers for all association data storage as well as simple methods on memeber data
 */
class AssoData
{
public:
	/*
	 *  members include
	 *  m_phenotype: phenotype vector (n X 1)
	 *  m_genotype: genotype matrix (n samples X p loci) -- raw genotypes
	 *  m_X: genotype vector (n X 1) -- genotype scores
	 *      Each element is a score combining information from an individual's p loci
	 *  m_maf: minor allele frequency vector (p X 1) for the p loci
	 *  m_statistic: statistics for all predictor variables from association tests.
	 *      For permutation based tests only genotype statistics is calculated (m_statistic.size() == 1)
	 *  m_se: standard error (or measures reflecting standard error) for m_statistic
	 *  m_pval: p-values for m_statistic
	 *  m_C: covariates matrix (n X m).
	 *      The first column of m_C is a vector of 1's (reserved for intercept estimate)
	 *      The last column of m_C will be filled with m_X
	 *  m_ncovar: number of covariates
	 *  m_xbar: mean(m_X)
	 *  m_ybar: mean(m_phenotype)
	 *  m_ncases: number of cases for case/ctrl data
	 *  m_nctrls: number of ctrls for case/ctrl data
	 *  m_model: model data object
	 *      Will be linear model for quantitative traits, logistic regression model for disease traits
	 */

	AssoData() :
		m_phenotype(0), m_genotype(0), m_maf(0), m_X(0),
		m_statistic(0), m_se(0), m_pval(0), m_C(0), m_ncovar(0),
		m_xbar(0.0), m_ybar(0.0), m_ncases(0),
		m_nctrls(0), m_model()
	{
	}


	~AssoData(){}

	// make a copy of the data
	virtual AssoData * clone() const
	{
		return new AssoData(*this);
	}


	// set raw genotypes
	double setGenotype(const matrixf & g)
	{
		//codings are 0, 1, 2, -9 and U(0,1) number for "expected" genotype
		m_genotype = g;
		return 0.0;
	}


	// set genotype scores
	double setX(const vectorf & g)
	{
		m_X = g;
		return 0.0;
	}


	// set phenotypes
	// input is a vector, i.e., no there is no other phenotype covariates and genotype will be the only predictor in the statistical model
	// As a result it re-sets the size of statistics, pvalues, etc, to 1
	double setPhenotype(const vectorf & p)
	{
		m_phenotype = p;
		// set phenotype mean
		m_ybar = std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0);
		m_ybar /= (1.0 * m_phenotype.size());
		// set case/ctrl counts
		//m_ncases = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(),
		//    std::bind2nd(std::equal_to<double>(),1.0));
		//m_nctrls = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(),
		//    std::bind2nd(std::equal_to<double>(),0.0));
		//
		// re-size statistics vector, etc, for there is only one predictor
		m_statistic.resize(1);
		m_se.resize(1);
		m_pval.resize(m_statistic.size());
		return m_ybar;
	}


	// set phenotypes and phenotype covariates
	// input is a phenotype vector and a covariates matrix
	double setPhenotype(const vectorf & p, const matrixf & c)
	{
		m_phenotype = p;
		m_C = c;
		// number of covariates. (need to subtract by 1 because first column of m_C is reserved with a vector of 1's)
		m_ncovar = c.size() - 1;
		vectorf one(p.size());
		std::fill(one.begin(), one.end(), 1.0);
		// reserve the last column to be a column of 1's, a placeholder for m_X which will be calculated and assigned to the last column of m_C
		m_C.push_back(one);
		// initialize the statistical model
		m_model.clear();
		m_model.setY(m_phenotype);
		m_model.setX(m_C);
		// set phenotype mean
		m_ybar = std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0);
		m_ybar /= (1.0 * m_phenotype.size());
		//m_ncases = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(),
		//    std::bind2nd(std::equal_to<double>(),1.0));
		//m_nctrls = (unsigned) std::count_if(m_phenotype.begin(), m_phenotype.end(),
		//    std::bind2nd(std::equal_to<double>(),0.0));
		//
		// re-size statistics vector, etc
		m_statistic.resize((m_ncovar + 1));
		m_se.resize((m_ncovar + 1));
		m_pval.resize(m_statistic.size());
		return m_ybar;
	}


	// set minor allele frequency from an external field
	// FIXME this function is not currently being used
	void setMaf(const vectorf & maf)
	{
		//get this field directly from the variant table
		m_maf = maf;

	}


	// set sample maf from data
	void setMaf()
	{
		//is problematic for variants on male chrX
		//but should be Ok if only use the relative mafs (e.g., weightings)

		m_maf.resize(m_genotype.front().size());
		std::fill(m_maf.begin(), m_maf.end(), 0.0);
		vectorf valid_all = m_maf;

		for (size_t j = 0; j < m_maf.size(); ++j) {
			// calc maf and loci counts for site j
			for (size_t i = 0; i < m_genotype.size(); ++i) {
				// genotype not missing
				if (!(m_genotype[i][j] < 0.0)) {
					valid_all[j] += 1.0;
					if (m_genotype[i][j] > 0.0) {
						m_maf[j] += m_genotype[i][j];
					}
				}
			}

			if (valid_all[j] > 0.0) {
				m_maf[j] = m_maf[j] / (valid_all[j] * 2.0);
			}
			//  FIXME : re-code genotype.  will be incorrect for male chrX
			if (m_maf[j] > 0.5) {
				m_maf[j] = 1.0 - m_maf[j];
				// recode genotypes
				for (size_t i = 0; i < m_genotype.size(); ++i) {
					// genotype not missing
					if (!(m_genotype[i][j] < 0.0)) {
						m_genotype[i][j] = 2.0 - m_genotype[i][j];
					}
				}
			}
		}
		/*
		   m_maf = std::accumulate(m_genotype.begin() + 1, m_genotype.end(),
		    m_genotype.front(), vplus);
		   std::transform(m_maf.begin(), m_maf.end(), m_maf.begin(),
		    std::bind2nd(std::divides<double>(), 2.0*m_genotype.size()));
		 */
	}


	// compute weight w = 1 / sqrt(p*(1-p))
	// using the entire sample
	void setMafWeight()
	{
		if (m_maf.size() == 0) {
			throw RuntimeError("MAF has not been calculated. Please calculate MAF prior to calculating weights.");
		}
		//
		m_weight.clear();
		for (size_t i = 0; i < m_maf.size(); ++i) {
			if (fEqual(m_maf[i], 0.0) || fEqual(m_maf[i], 1.0)) {
				m_weight.push_back(0.0);
			} else{
				m_weight.push_back(1.0 / sqrt(m_maf[i] * (1.0 - m_maf[i])));
			}
		}
	}


	// compute weight by maf from selected samples (ctrls, low QT samples, etc)
	// FIXME not yet implemented
	void setMafWeight(const std::vector<size_t> & idx)
	{
	}


	// calculate mean of m_X (genotype scores)
	double meanOfX()
	{
		m_xbar = std::accumulate(m_X.begin(), m_X.end(), 0.0);
		m_xbar /= (1.0 * m_X.size());
		return m_xbar;
	}


	/*
	 * get various data
	 *
	 * */
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


	unsigned covarcounts()
	{
		return m_ncovar;
	}


	unsigned samplecounts()
	{
		return m_phenotype.size();
	}


	vectorf pvalue()
	{
		return m_pval;
	}


	vectorf statistic()
	{
		return m_statistic;
	}


	vectorf se()
	{
		return m_se;
	}


	/*
	 * set p-value, etc
	 *
	 * */

	double setPvalue(vectorf pval)
	{
		m_pval = pval;
		return 0.0;
	}


	double setPvalue(double pval)
	{
		m_pval.resize(1);
		m_pval[0] = pval;
		return 0.0;
	}


	double setStatistic(vectorf stat)
	{
		m_statistic = stat;
		return 0.0;
	}


	double setSE(vectorf se)
	{
		m_se = se;
		return 0.0;
	}


	double setStatistic(double stat)
	{
		m_statistic.resize(1);
		m_statistic[0] = stat;
		return 0.0;
	}


	double setSE(double se)
	{
		m_se.resize(1);
		m_se[0] = se;
		return 0.0;
	}


public:
	// permute phenotype data
	void permuteY()
	{
		random_shuffle(m_phenotype.begin(), m_phenotype.end());
	}


	// permute genotype matrix by samples
	void permuteRawX()
	{
		random_shuffle(m_genotype.begin(), m_genotype.end());
	}


	// permute genotype scores
	void permuteX()
	{
		random_shuffle(m_X.begin(), m_X.end());
	}


	// calculate genotype scores from raw genotype
	// m_X = rowSum(m_genotype)
	void sumToX()
	{
		m_X.resize(m_genotype.size());
		std::fill(m_X.begin(), m_X.end(), 0.0);
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			//m_X[i] = std::accumulate(m_genotype[i].begin(), m_genotype[i].end(), 0.0);
			for (size_t j = 0; j < m_genotype[i].size(); ++j) {
				if (m_genotype[i][j] > 0) {
					m_X[i] += m_genotype[i][j];
				}
			}
		}
	}


	// calculate genotype scores from raw genotype
	// m_X = (rowSum(m_genotype) > 0)
	// m_X is binary: 0 for all wildtype, 1 for having at least one mutation
	void binToX()
	{
		m_X.resize(m_genotype.size());
		std::fill(m_X.begin(), m_X.end(), 0.0);
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			double pnovar = 1.0;
			for (size_t j = 0; j != m_genotype[i].size(); ++j) {
				if (m_genotype[i][j] >= 1.0) {
					m_X[i] = 1.0;
					break;
				} else if (m_genotype[i][j] > 0.0) {
					// binning the data with proper handling of missing genotype
					pnovar *= (1.0 - m_genotype[i][j]);
				} else ;
			}
			// all genotypes are missing: have to be represented as Pr(#mutation>=1)
			if (pnovar < 1.0 && m_X[i] < 1.0) {
				m_X[i] = 1.0 - pnovar;
			}
		}
	}


	// weight the raw genotypes by m_weight
	void weightX()
	{
		if (m_weight.size() == 0) {
			return;
		}
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			for (size_t j = 0; j < m_genotype[i].size(); ++j) {
				if (m_genotype[i][j] > 0) {
					m_genotype[i][j] *= m_weight[j];
				}
			}
		}
	}


	// remove variant sites having MAF <= lower_bound or MAF > upper_bound
	void setSitesByMaf(double upper, double lower)
	{
		if (upper > 1.0) {
			throw ValueError("Minor allele frequency value should not exceed 1");
		}
		if (lower < 0.0) {
			throw ValueError("Minor allele frequency should be a positive value");
		}

		if (fEqual(upper, 1.0) && fEqual(lower, 0.0)) return;

		for (size_t j = 0; j != m_maf.size(); ++j) {
			if (m_maf[j] <= lower || m_maf[j] > upper) {
				m_maf.erase(m_maf.begin() + j);
				for (size_t i = 0; i < m_genotype.size(); ++i) {
					m_genotype[i].erase(m_genotype[i].begin() + j);
				}
				--j;
			}
		}
		return;
	}


	// Wald's statistic for simple linear model Y = b0 + b1x
	void simpleLinear()
	{
		// simple linear regression score test
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
			double b0 = m_ybar - b1 * m_xbar;

			//SSE
			for (size_t i = 0; i != m_X.size(); ++i) {
				ysigma += pow(m_phenotype[i] - (b0 + b1 * m_X[i]), 2.0);
			}
			double varb = ysigma / (m_phenotype.size() - 2.0) / denominator;
			m_statistic[0] = b1 / sqrt(varb);
			m_se[0] = sqrt(varb);
		} else {
			m_statistic[0] = 0.0;
			m_se[0] = std::numeric_limits<double>::quiet_NaN();
		}
	}


	//!- Score test implementation for logistic regression model logit(p) = b0 + b1x
	void simpleLogit()
	{
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

		ss = ss / sqrt(vm1);

		//!-FIXME: (not sure why this happens)
		//!- w/ rounding to 0 I get strange number such as 3.72397e-35
		//!- this would lead to type I error problem
		fRound(ss, 0.0001);
		m_statistic[0] = ss;
		m_se[0] = sqrt(vm1);
	}


	// fitting / calculating wald's statistic for multiple linear regression model
	// Y = b0 + b1x1 + b2x2 + ... + bnxn
	void multipleLinear()
	{
		//!- multiple linear regression parameter estimate
		//!- BETA= (X'X)^{-1}X'Y => (X'X)BETA = X'Y
		//!- Solve the system via gsl_linalg_SV_solve()
		if (m_X.size() != m_phenotype.size()) {
			throw ValueError("Genotype/Phenotype length not equal!");
		}
		// reset phenotype data
		m_model.replaceCol(m_phenotype, 0);
		// reset genotype data
		m_model.replaceCol(m_X, m_C.size() - 1);
		// fit the linear regression model
		m_model.fit();
		// get statistics and set values to assoData object
		vectorf beta = m_model.getBeta();
		vectorf seb = m_model.getSEBeta();
		m_statistic[0] = beta.back() / seb.back();
		m_se[0] = seb.back();
		for (unsigned i = 1; i < beta.size() - 1; ++i) {
			m_statistic[i] = beta[i] / seb[i];
			m_se[i] = seb[i];
		}
	}


	// get p-value from statistic, assuming standard normal distribution
	void gaussianP(unsigned sided = 1)
	{
		if (sided == 1) {
			for (unsigned i = 0; i < m_statistic.size(); ++i) {
				m_pval[i] = gsl_cdf_ugaussian_Q(m_statistic[i]);
			}
		} else if (sided == 2) {
			for (unsigned i = 0; i < m_statistic.size(); ++i) {
				m_pval[i] = gsl_cdf_chisq_Q(m_statistic[i] * m_statistic[i], 1.0);
			}
		} else {
			throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
		}
	}


	// get p-value from statistic, assuming t distribution
	// degree of freedom will be calculated from data automatically
	void studentP(unsigned sided = 1)
	{
		// df = n - p where p = #covariates + 1 (for beta1) + 1 (for beta0) = m_ncovar+2
		if (sided == 1) {
			for (unsigned i = 0; i < m_statistic.size(); ++i) {
				m_pval[i] = gsl_cdf_tdist_Q(m_statistic[i], m_phenotype.size() - (m_ncovar + 2.0));
			}
		} else if (sided == 2) {
			for (unsigned i = 0; i < m_statistic.size(); ++i) {
				double p = gsl_cdf_tdist_Q(m_statistic[i], m_phenotype.size() - (m_ncovar + 2.0));
				m_pval[i] = fmin(p, 1.0 - p) * 2.0;
			}
		} else {
			throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
		}
	}


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

	// weights
	vectorf m_weight;
	// genotype/phenotype data summaries
	double m_xbar;
	double m_ybar;
	unsigned m_ncases;
	unsigned m_nctrls;
	// statistics from association tests
	vectorf m_pval;
	vectorf m_statistic;
	vectorf m_se;
	// statistical model object
	LinearM m_model;
};

}
#endif
