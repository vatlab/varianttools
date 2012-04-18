/**
 *  $File: assoData.cpp $
 *  $LastChangedDate$
 *  $Rev$
 *
 *  This file is part of variant_tools, a software application to annotate,
 *  summarize, and filter variants for next-gen sequencing ananlysis.
 *  Please visit http://varianttools.sourceforge.net for details.
 *
 *  Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
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


#include "assoData.h"

namespace vtools {


bool AssoData::setPhenotype(const vectorf & p)
{
	m_phenotype = p;
	// set phenotype mean
	setVar("ybar", (double)std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0) / (1.0 * m_phenotype.size()));
	// set case/ctrl counts
	//setVar("ncases", (int) std::count_if(m_phenotype.begin(), m_phenotype.end(),
	//    std::bind2nd(std::equal_to<double>(),1.0)));
	//setVar("nctrls", (int) std::count_if(m_phenotype.begin(), m_phenotype.end(),
	//    std::bind2nd(std::equal_to<double>(),0.0)));
	//
	// re-size statistics vector, etc, for there is only one predictor
	m_statistic.resize(1);
	m_se.resize(1);
	m_pval.resize(m_statistic.size());
	return true;
}


bool AssoData::setPhenotype(const vectorf & p, const matrixf & c)
{
	m_phenotype = p;
	m_C = c;
	// number of covariates. (need to subtract by 1 because first column of m_C is reserved with a vector of 1's)
	setVar("ncovar", int(c.size() - 1));
	vectorf one(p.size());
	std::fill(one.begin(), one.end(), 1.0);
	// reserve the last column to be a column of 1's, a placeholder for m_X which will be calculated and assigned to the last column of m_C
	m_C.push_back(one);
	// initialize the statistical model
	m_model.clear();
	m_model.setY(m_phenotype);
	m_model.setX(m_C);
	// set phenotype mean
	setVar("ybar", (double)std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0) / (1.0 * m_phenotype.size()));
	// set case/ctrl counts
	//setVar("ncases", (int) std::count_if(m_phenotype.begin(), m_phenotype.end(),
	//    std::bind2nd(std::equal_to<double>(),1.0)));
	//setVar("nctrls", (int) std::count_if(m_phenotype.begin(), m_phenotype.end(),
	//    std::bind2nd(std::equal_to<double>(),0.0)));
	// re-size statistics vector, etc
	m_statistic.resize(c.size());
	m_se.resize(m_statistic.size());
	m_pval.resize(m_statistic.size());
	return true;
}


void AssoData::sumToX()
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
	setVar("xbar", (double)std::accumulate(m_X.begin(), m_X.end(), 0.0) / (1.0 * m_X.size()));
}


void AssoData::binToX()
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
	setVar("xbar", (double)std::accumulate(m_X.begin(), m_X.end(), 0.0) / (1.0 * m_X.size()));
}


void AssoData::weightX(const vectorf & weight)
{
	if (weight.size() == 0) {
		return;
	}
	for (size_t i = 0; i < m_genotype.size(); ++i) {
		for (size_t j = 0; j < m_genotype[i].size(); ++j) {
			if (m_genotype[i][j] > 0) {
				m_genotype[i][j] *= weight[j];
			}
		}
	}
}


void AssoData::simpleLinear()
{
	// simple linear regression score test
	//!- See page 23 and 41 of Kutner's Applied Linear Stat. Model, 5th ed.
	//

	double xbar = getDoubleVar("xbar");
	double ybar = getDoubleVar("ybar");

	if (m_X.size() != m_phenotype.size()) {
		throw ValueError("Genotype/Phenotype length not equal!");
	}
	double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
	for (size_t i = 0; i != m_X.size(); ++i) {
		numerator += (m_X[i] - xbar) * m_phenotype[i];
		denominator += pow(m_X[i] - xbar, 2.0);
	}

	if (!fEqual(numerator, 0.0)) {
		//!- Compute MSE and V[\hat{beta}]
		//!- V[\hat{beta}] = MSE / denominator
		double b1 = numerator / denominator;
		double b0 = ybar - b1 * xbar;

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


void AssoData::simpleLogit()
{
	//!- labnotes vol.2 page 3
	//!- input phenotypes have to be binary values 0 or 1
	if (m_X.size() != m_phenotype.size()) {
		throw ValueError("Genotype/Phenotype length not equal!");
	}
	//double ebo = (1.0 * n1) / (1.0 * (m_phenotype.size()-n1));
	//double bo = log(ebo);
	double xbar = getDoubleVar("xbar");
	double po = (1.0 * getIntVar("ncases")) / (1.0 * m_phenotype.size());
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

	ss = ss / sqrt(vm1);

	//!-FIXME: (not sure why this happens)
	//!- w/ rounding to 0 I get strange number such as 3.72397e-35
	//!- this would lead to type I error problem
	fRound(ss, 0.0001);
	m_statistic[0] = ss;
	m_se[0] = sqrt(vm1);
}


void AssoData::multipleLinear()
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


void AssoData::gaussianP(unsigned sided)
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


void AssoData::studentP(unsigned sided)
{
	int ncovar = getIntVar("ncovar");

	// df = n - p where p = #covariates + 1 (for beta1) + 1 (for beta0) = ncovar+2
	if (sided == 1) {
		for (unsigned i = 0; i < m_statistic.size(); ++i) {
			m_pval[i] = gsl_cdf_tdist_Q(m_statistic[i], m_phenotype.size() - (ncovar + 2.0));
		}
	} else if (sided == 2) {
		for (unsigned i = 0; i < m_statistic.size(); ++i) {
			double p = gsl_cdf_tdist_Q(m_statistic[i], m_phenotype.size() - (ncovar + 2.0));
			m_pval[i] = fmin(p, 1.0 - p) * 2.0;
		}
	} else {
		throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
	}
}


}
