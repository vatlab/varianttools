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
#include <map>
#include <string>
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
private:
	typedef std::map<std::string, double> DoubleVars;
	typedef std::map<std::string, vectorf> ArrayVars;

public:
	/*
	 *  members include
	 *  m_phenotype: phenotype vector (n X 1)
	 *  m_genotype: genotype matrix (n samples X p loci) -- raw genotypes
	 *  m_X: genotype vector (n X 1) -- genotype scores
	 *      Each element is a score combining information from an individual's p loci
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
		m_phenotype(0), m_genotype(0), m_X(0),
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
	double setPhenotype(const vectorf & p);


	// set phenotypes and phenotype covariates
	// input is a phenotype vector and a covariates matrix
	double setPhenotype(const vectorf & p, const matrixf & c);


	// calculate mean of m_X (genotype scores)
	double meanOfX();

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


	matrixf & raw_genotype()
	{
		return m_genotype;
	}


	matrixf covariates()
	{
		return m_model.getX();
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
	void sumToX();

	// calculate genotype scores from raw genotype
	// m_X = (rowSum(m_genotype) > 0)
	// m_X is binary: 0 for all wildtype, 1 for having at least one mutation
	void binToX();

	// weight the raw genotypes by given weight
	void weightX(const vectorf & weight);

	// Wald's statistic for simple linear model Y = b0 + b1x
	void simpleLinear();

	//!- Score test implementation for logistic regression model logit(p) = b0 + b1x
	void simpleLogit();

	// fitting / calculating wald's statistic for multiple linear regression model
	// Y = b0 + b1x1 + b2x2 + ... + bnxn
	void multipleLinear();

	// get p-value from statistic, assuming standard normal distribution
	void gaussianP(unsigned sided = 1);

	// get p-value from statistic, assuming t distribution
	// degree of freedom will be calculated from data automatically
	void studentP(unsigned sided = 1);

	// store a double value with name 'name'
	void setVar(const string & name, const double value)
	{
		m_doubleVars[name] = value;
	}

	// store an array with name 'name'
	void setVar(const string & name, const vectorf & value)
	{
		m_arrayVars[name] = value;
	}

	bool hasVar(const string & name)
	{
		return (m_doubleVars.find(name) != m_doubleVars.end() ||
			m_arrayVars.find(name) != m_arrayVars.end());
	}

	double & getDoubleVar(const string & name)
	{
		DoubleVars::iterator it = m_doubleVars.find(name);
		if (it == m_doubleVars.end())
			throw ValueError("No double variable with name " + name + " can be found");
		return it->second;
	}

	vectorf & getArrayVar(const string & name)
	{
		ArrayVars::iterator it = m_arrayVars.find(name);
		if (it == m_arrayVars.end())
			throw ValueError("No array with name " + name + " can be found");
		return it->second;
	}

private:
	/// raw phenotype and gneotype data
	vectorf m_phenotype;
	matrixf m_genotype;

	// covariates
	matrixf m_C;
	unsigned m_ncovar;

	/// translated genotype
	vectorf m_X;

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

	// arbitrary double type of variables
	DoubleVars m_doubleVars;
	// arbitrary vectorf type of variables.
	ArrayVars m_arrayVars;
};

}
#endif
