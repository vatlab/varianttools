/**
 *  $File: action.cpp $
 *  $LastChangedDate: 2012-04-13 21:15:30 -0700 (Fri, 13 Apr 2012) $
 *  $Rev: 1075 $
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
#include <cfloat>
#include "action.h"
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"
#include "fisher2.h"
#include "time.h"

namespace vtools {

class Timer
{
public:
	Timer(int timeout) : m_total_time(timeout), m_calls_per_second(0)
	{
		m_start_time = time(NULL);
	}


	bool timeout(int caller_id = 0)
	{
		if (m_total_time == 0)
			return false;
		//count how many times the timeout() function is called
		//before a detectable time change happens (difftime > 0).
		if (m_calls_per_second <= 0) {
			// first detectable change
			int diff = static_cast<int>(difftime(time(NULL), m_start_time));
			if (diff)
				// m_calls_per_second might still be 0 when 1s is passed
				m_calls_per_second = std::max(1, -1 * m_calls_per_second);
			else
				m_calls_per_second--;
			return diff >= m_total_time;
		} else
			//check time() only every second
			return (caller_id % m_calls_per_second == 0) &&
			       static_cast<int>(difftime(time(NULL), m_start_time)) >= m_total_time;
	}


	int remaining()
	{
		// always return at least 1 second when timeout is set because 0 means no time limit is set
		return m_total_time > 0 ? std::max(1, m_total_time -
			static_cast<int>(difftime(time(NULL), m_start_time))) : 0;
	}


private:
	time_t m_start_time;
	int m_total_time;
	int m_calls_per_second;
};

bool SetMaf::apply(AssoData & d, int timeout)
{
	matrixf & genotype = d.raw_genotype();

	//is problematic for variants on male chrX
	//but should be Ok if only use the relative mafs (e.g., weightings)
	vectorf maf(genotype.front().size(), 0.0);
	// minor allele counts
	vectorf mac = maf;
	vectorf valid_all = maf;

	for (size_t j = 0; j < mac.size(); ++j) {
		// calc mac and loci counts for site j
		for (size_t i = 0; i < genotype.size(); ++i) {
			// genotype not missing
			if (genotype[i][j] == genotype[i][j]) {
				valid_all[j] += 1.0;
				if (genotype[i][j] >= 1.0) {
					mac[j] += genotype[i][j];
				}
			}
		}

		if (valid_all[j] > 0.0) {
			maf[j] = mac[j] / (valid_all[j] * 2.0);
		}
		//  FIXME : re-code genotype.  will be incorrect for male chrX
		if (maf[j] > 0.5) {
			maf[j] = 1.0 - maf[j];
			mac[j] = valid_all[j] * 2.0 - mac[j];
			// recode genotypes
			for (size_t i = 0; i < genotype.size(); ++i) {
				// genotype not missing
				if (genotype[i][j] == genotype[i][j]) {
					genotype[i][j] = 2.0 - genotype[i][j];
				}
			}
		}
	}
	d.setVar("mac", mac);
	d.setVar("maf", maf);
	// actual sample size
	d.setVar("maf_denominator", valid_all);
	return true;
}


bool FillGMissing::apply(AssoData & d, int timeout)
{
	if (!d.hasVar("maf")) {
		throw RuntimeError("Sample MAF is required for this operation but has not been calculated.");
	}

	RNG rng;
	gsl_rng * gslr = rng.get();
	double multiplier = (d.getIntVar("moi") > 0) ? d.getIntVar("moi") * 1.0 : 1.0;
	vectorf & maf = d.getArrayVar("maf");
	vectorf & mac = d.getArrayVar("mac");
	matrixf & genotype = d.raw_genotype();
	for (size_t j = 0; j < maf.size(); ++j) {
		// scan for site j
		for (size_t i = 0; i < genotype.size(); ++i) {
			// genotype missing; replace with maf
			if (genotype[i][j] != genotype[i][j]) {
				if (m_method == "maf") {
					genotype[i][j] = maf[j] * multiplier;
					// update the MAC
					// FIXME: not updating MAF here
					// although this would not make much difference
					mac[j] += genotype[i][j];
				}else if (m_method == "mlg") {
					// most likely genotype
					double p = gsl_rng_uniform(gslr);
					if (p > maf[j] * maf[j]) genotype[i][j] = 2.0;
					else if (p < maf[j]) genotype[i][j] = 0.0;
					else genotype[i][j] = 1.0;
					switch (d.getIntVar("moi")) {
					case 0:
						// recessive
					{
						genotype[i][j] = (fEqual(genotype[i][j], 2.0)) ? 1.0 : 0.0;
					}
					break;
					case 1:
						// dominant
					{
						genotype[i][j] = (!fEqual(genotype[i][j], 0.0)) ? 1.0 : 0.0;
					}
					break;
					default:
						break;
					}
				}
			}
		}
	}
	return true;
}


bool WeightByInfo::apply(AssoData & d, int timeout)
{
	for (size_t i = 0; i < m_info.size(); ++i) {
		if (d.hasVar("__var_" + m_info[i])) {
			// get var_info
			d.weightX(d.getArrayVar("__var_" + m_info[i]));
		} else if (d.hasVar("__geno_" + m_info[i])) {
			// get geno_info
			d.weightX(d.getMatrixVar("__geno_" + m_info[i]));
		} else if (d.hasVar(m_info[i])) {
			// get internal weight matrix, should be a vectorf, only one vector
			d.weightX(d.getMatrixVar(m_info[i])[0]);
		} else {
			throw ValueError("Cannot find genotype/variant information: " + m_info[i]);
		}
	}
	return true;
}


bool BrowningWeight::apply(AssoData & d, int timeout)
{
	double multiplier = (d.getIntVar("moi") > 0) ? d.getIntVar("moi") * 1.0 : 1.0;
	matrixf maf(m_model);
	matrixf valid_all(m_model);

	if (m_model == 0) {
		// Browing transformation by entire sample maf
		if (!d.hasVar("maf")) {
			throw RuntimeError("MAF has not been calculated. Please calculate it prior to calculating weights.");
		} else {
			maf.push_back(d.getArrayVar("maf"));
			valid_all.push_back(d.getArrayVar("maf_denominator"));
		}
	} else {
		// m_model == 1 : weight by "ctrl" only. The weight matrices will have one column only
		// m_model == 2 : weight by both models. The weight matrices will have 2 columns
		matrixf & genotype = d.raw_genotype();
		vectorf & phenotype = d.phenotype();
		double ybar = d.getDoubleVar("ybar");
		for (size_t s = 0; s < m_model; ++s) {
			maf[s].resize(genotype.front().size());
			std::fill(maf[s].begin(), maf[s].end(), 0.0);
			valid_all[s].resize(maf[s].size(), 0.0);
			std::fill(valid_all[s].begin(), valid_all[s].end(), 0.0);
		}
		for (size_t s = 0; s < m_model; ++s) {
			for (size_t j = 0; j < maf[s].size(); ++j) {
				// calc af and loci counts for site j
				for (size_t i = 0; i < genotype.size(); ++i) {
					if (s == 0) {
						// model 1: weight by ctrl
						if (phenotype[i] > ybar) continue;
					} else {
						// model 2: an additional weight by case
						if (phenotype[i] < ybar) continue;
					}
					// genotype not missing
					if (genotype[i][j] == genotype[i][j]) {
						valid_all[s][j] += 1.0;
						if (genotype[i][j] >= 1.0) {
							maf[s][j] += genotype[i][j];
						}
					}
				}

				if (valid_all[s][j] > 0.0) {
					maf[s][j] = (maf[s][j] + 1.0) / ((valid_all[s][j] + 1.0) * multiplier);
				}
			}
		}
	}
	// calculte weight. maf matrix will become the weight matrix
	for (size_t s = 0; s < maf.size(); ++s) {
		for (size_t i = 0; i < maf[s].size(); ++i) {
			if (fEqual(maf[s][i], 0.0) || fEqual(maf[s][i], 1.0) || fEqual(valid_all[s][i], 0.0)) {
				maf[s][i] = 0.0;
			} else {
				maf[s][i] = 1.0 / sqrt(maf[s][i] * (1.0 - maf[s][i]) * valid_all[s][i] * multiplier);
			}
		}
	}
	d.setVar("Browningweight", maf);
	return true;
}


bool SetSites::apply(AssoData & d, int timeout)
{
	if (!d.hasVar("maf") || !d.hasVar("mac")) {
		throw RuntimeError("MAF/MAC has not been calculated. Please calculate MAF/MAC prior to setting variant sites.");
	}
	//FIXME: all weights have to be trimed as well
	vectorf & maf = d.getArrayVar("maf");
	vectorf & mac = d.getArrayVar("mac");
	matrixf & genotype = d.raw_genotype();
	//
	if (!m_use_mac && m_upper > 1.0) {
		throw ValueError("Minor allele frequency value should not exceed 1");
	}
	if (!m_use_mac && m_lower < 0.0) {
		throw ValueError("Minor allele frequency should be a positive value");
	}

	if (!m_use_mac && fEqual(m_upper, 1.0) && fEqual(m_lower, 0.0))
		return true;

	for (size_t j = 0; j < maf.size(); ++j) {
		double im = m_use_mac ? mac[j] : maf[j];
		if (im <= m_lower || im > m_upper) {
			maf.erase(maf.begin() + j);
			mac.erase(mac.begin() + j);
			for (size_t i = 0; i < genotype.size(); ++i) {
				genotype[i].erase(genotype[i].begin() + j);
			}
			--j;
		}
	}
	return true;
}


bool CodeXByMOI::apply(AssoData & d, int timeout)
{
	int moi = d.getIntVar("moi");

	// additive mode, nothing to do
	if (moi == 2) return true;
	//
	matrixf & genotype = d.raw_genotype();
	// have to reset MAF
	vectorf valid_all = d.getArrayVar("maf_denominator");
	vectorf maf(genotype.front().size(), 0.0);
	vectorf mac = maf;
	//
	for (size_t j = 0; j < genotype.front().size(); ++j) {
		for (size_t i = 0; i < genotype.size(); ++i) {
			// missing genotype
			if (genotype[i][j] != genotype[i][j]) continue;
			// recessive
			if (moi == 0) genotype[i][j] = (fEqual(genotype[i][j], 2.0)) ? 1.0 : 0.0;
			// dominant
			if (moi == 1) genotype[i][j] = (!fEqual(genotype[i][j], 0.0)) ? 1.0 : 0.0;
			// re-calculate minor allele counts
			if (genotype[i][j] >= 1.0) mac[j] += genotype[i][j];
		}
		// denominator is total site counts itself for dominant and recessive
		if (valid_all[j] > 0.0) maf[j] = mac[j] / (valid_all[j] * 1.0);
	}
	//
	d.setVar("mac", mac);
	d.setVar("maf", maf);
	return true;
}


bool SumToX::apply(AssoData & d, int timeout)
{
	vectorf & X = d.genotype();
	matrixf & genotype = d.raw_genotype();

	X.resize(genotype.size());
	std::fill(X.begin(), X.end(), 0.0);
	for (size_t i = 0; i < genotype.size(); ++i) {
		//X[i] = std::accumulate(genotype[i].begin(), genotype[i].end(), 0.0);
		for (size_t j = 0; j < genotype[i].size(); ++j) {
			if (genotype[i][j] > 0.0) {
				X[i] += genotype[i][j];
			}
		}
	}
	d.setVar("xbar", (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size()));
	return true;
}


bool BinToX::apply(AssoData & d, int timeout)
{
	vectorf & X = d.genotype();
	matrixf & genotype = d.raw_genotype();

	X.resize(genotype.size());
	std::fill(X.begin(), X.end(), 0.0);
	for (size_t i = 0; i < genotype.size(); ++i) {
		double pnovar = 1.0;
		for (size_t j = 0; j != genotype[i].size(); ++j) {
			if (genotype[i][j] >= 1.0) {
				X[i] = 1.0;
				break;
			} else if (genotype[i][j] > 0.0) {
				// binning the data with proper handling of missing genotype
				pnovar *= (1.0 - genotype[i][j] / 2.0); // genotype[i][j]/2.0 would be the maf
			} else ;
		}
		// all genotypes are missing: have to be represented as Pr(#mutation>=1)
		if (pnovar < 1.0 && X[i] < 1.0) {
			X[i] = 1.0 - pnovar;
		}
	}
	d.setVar("xbar", (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size()));
	return true;
}


bool SimpleLinearRegression::apply(AssoData & d, int timeout)
{
	// simple linear regression score test
	//!- See page 23 and 41 of Kutner's Applied Linear Stat. Model, 5th ed.

	vectorf & X = d.genotype();
	vectorf & Y = d.phenotype();

	if (X.size() != Y.size()) {
		throw ValueError("Genotype/Phenotype length not equal!");
	}
	//
	double ybar = d.getDoubleVar("ybar");
	double xbar;
	try {
		xbar = d.getDoubleVar("xbar");
	} catch (...) {
		xbar = (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size());
	}
	//
	double numerator = 0.0, denominator = 0.0, ysigma = 0.0;
	for (size_t i = 0; i != X.size(); ++i) {
		numerator += (X[i] - xbar) * Y[i];
		denominator += pow(X[i] - xbar, 2.0);
	}

	if (!fEqual(numerator, 0.0)) {
		//!- Compute MSE and V[\hat{beta}]
		//!- V[\hat{beta}] = MSE / denominator
		double b1 = numerator / denominator;
		double b0 = ybar - b1 * xbar;

		//SSE
		for (size_t i = 0; i != X.size(); ++i) {
			ysigma += pow(Y[i] - (b0 + b1 * X[i]), 2.0);
		}
		double varb = ysigma / (Y.size() - 2.0) / denominator;
		d.setStatistic(b1);
		d.setSE(sqrt(varb));
	} else {
		d.setStatistic(0.0);
		d.setSE(std::numeric_limits<double>::quiet_NaN());
	}
	return true;
}


bool SimpleLogisticRegression::apply(AssoData & d, int timeout)
{
	//!- labnotes vol.2 page 3
	vectorf & X = d.genotype();
	vectorf & Y = d.phenotype();

	if (X.size() != Y.size()) {
		throw ValueError("Genotype/Phenotype length not equal");
	}
	//
	double xbar;
	try {
		xbar = d.getDoubleVar("xbar");
	} catch (...) {
		xbar = (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size());
	}
	//
	//double ebo = (1.0 * n1) / (1.0 * (Y.size()-n1));
	//double bo = log(ebo);
	double po = (1.0 * d.getIntVar("ncases")) / (1.0 * Y.size());
	double ss = 0.0, vm1 = 0.0;
	// the score, and variance of score under the null
	for (size_t i = 0; i != X.size(); ++i) {
		ss += (X[i] - xbar) * (Y[i] - po);
		vm1 += (X[i] - xbar) * (X[i] - xbar) * po * (1.0 - po);
	}

	if (!fEqual(ss, 0.0)) {
		d.setStatistic(ss);
		d.setSE(sqrt(vm1));
	} else {
		d.setStatistic(0.0);
		d.setSE(std::numeric_limits<double>::quiet_NaN());
	}
	return true;
}


bool MultipleRegression::apply(AssoData & d, int timeout)
{

	vectorf & X = d.genotype();
	vectorf & Y = d.phenotype();
	matrixf & C = d.covariates();
	vectorf & beta = d.statistic();

	if (X.size() != Y.size()) {
		throw ValueError("Genotype/Phenotype length not equal!");
	}
	LMData & mdata = d.modeldata();
	// reset phenotype data
	mdata.replaceColumn(Y, 0);
	// reset genotype data
	mdata.replaceColumn(X, C.size() - 1);
	// fit the multiple regression model
	BaseLM * model = m_getModel();
	model->fit(mdata);
	// get statistics
	beta = mdata.getBeta();
	beta.erase(beta.begin());
	// move the last element to first
	std::rotate(beta.begin(), beta.end() - 1, beta.end());

	// get/set standard error
	if (m_iSE) {
		vectorf & seb = d.se();
		model->evalSE(mdata);
		seb = mdata.getSEBeta();
		seb.erase(seb.begin());
		std::rotate(seb.begin(), seb.end() - 1, seb.end());
	}
	delete model;
	return true;
}


bool GaussianPval::apply(AssoData & d, int timeout)
{
	vectorf & statistic = d.statistic();
	vectorf & se = d.se();
	vectorf & pval = d.pvalue();

	if (m_tailed == 1) {
		for (unsigned i = 0; i < statistic.size(); ++i) {
			if (fEqual(se[i], 0.0)) {
				throw ValueError("Standard Error has not been calculated");
			}
			pval[i] = gsl_cdf_ugaussian_Q(statistic[i] / se[i]);
		}
	} else if (m_tailed == 2) {
		for (unsigned i = 0; i < statistic.size(); ++i) {
			if (fEqual(se[i], 0.0)) {
				throw ValueError("Standard Error has not been calculated");
			}
			pval[i] = gsl_cdf_chisq_Q(statistic[i] / se[i] * statistic[i] / se[i], 1.0);
		}
	} else {
		throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
	}
	return true;
}


bool StudentPval::apply(AssoData & d, int timeout)
{
	int ncovar = d.getIntVar("ncovar");
	vectorf & statistic = d.statistic();
	vectorf & se = d.se();
	vectorf & pval = d.pvalue();

	// df = n - p where p = #covariates + 1 (for beta1) + 1 (for beta0) = ncovar+2
	if (m_tailed == 1) {
		for (unsigned i = 0; i < statistic.size(); ++i) {
			if (fEqual(se[i], 0.0)) {
				throw ValueError("Standard Error has not been calculated");
			}
			pval[i] = gsl_cdf_tdist_Q(statistic[i] / se[i], d.samplecounts() - (ncovar + 2.0));
		}
	} else if (m_tailed == 2) {
		for (unsigned i = 0; i < statistic.size(); ++i) {
			if (fEqual(se[i], 0.0)) {
				throw ValueError("Standard Error has not been calculated");
			}
			double p = gsl_cdf_tdist_Q(statistic[i] / se[i], d.samplecounts() - (ncovar + 2.0));
			pval[i] = std::min(p, 1.0 - p) * 2.0;
		}
	} else {
		throw ValueError("Alternative hypothesis should be one-sided (1) or two-sided (2)");
	}
	return true;
}


bool Fisher2X2::apply(AssoData & d, int timeout)
{
	vectorf & genotype = d.genotype();
	vectorf & phenotype = d.phenotype();

	if (genotype.size() != phenotype.size()) {
		throw ValueError("genotype/phenotype do not match");
	}
	vectori twotwoTable(4, 0);
	for (size_t i = 0; i != genotype.size(); ++i) {

		if (genotype[i] != genotype[i]) {
			throw ValueError("Input genotype data have missing entries");
		}
		if (!(fEqual(phenotype[i], 1.0) || fEqual(phenotype[i], 0.0))) {
			throw ValueError("Input phenotype data not binary");
		}

		if (fEqual(phenotype[i], 1.0)) {
			if (genotype[i] > 0.0)
				twotwoTable[0] += 1;
			else
				twotwoTable[1] += 1;
		}else {
			if (genotype[i] > 0.0)
				twotwoTable[2] += 1;
			else
				twotwoTable[3] += 1;
		}
	}

	double pvalue = 0.0;
	if (m_tailed == 1) {
		pvalue = (m_midp)
		         ? (twotwoTable[0] > 0) * gsl_cdf_hypergeometric_P(
			(twotwoTable[0] - 1),
			(twotwoTable[0] + twotwoTable[2]),
			(twotwoTable[1] + twotwoTable[3]),
			(twotwoTable[0] + twotwoTable[1])
		    ) + 0.5 * gsl_ran_hypergeometric_pdf(
			twotwoTable[0],
			(twotwoTable[0] + twotwoTable[2]),
			(twotwoTable[1] + twotwoTable[3]),
			(twotwoTable[0] + twotwoTable[1])
		    ) : gsl_cdf_hypergeometric_P(
			(twotwoTable[0] - 1),
			(twotwoTable[0] + twotwoTable[2]),
			(twotwoTable[1] + twotwoTable[3]),
			(twotwoTable[0] + twotwoTable[1])
		    );
		pvalue = 1.0 - pvalue;
	} else {
		pvalue = fexact2x2(twotwoTable, "two");
	}
	d.setStatistic( (twotwoTable[1] * twotwoTable[2]) ? (double)twotwoTable[3] * (double)twotwoTable[0] / ((double)twotwoTable[1] * (double)twotwoTable[2]) : DBL_MAX);
	d.setPvalue(pvalue);
	return true;
}


bool MannWhitneyu::apply(AssoData & d, int timeout)
{
	vectorf & genotype = d.genotype();
	vectorf & phenotype = d.phenotype();
	double ybar = d.getDoubleVar("ybar");
	//
	vectorf caseScores(0), ctrlScores(0);

	for (size_t i = 0; i < phenotype.size(); ++i) {
		if (phenotype[i] > ybar) {
			caseScores.push_back(genotype[i]);
		}else {
			ctrlScores.push_back(genotype[i]);
		}
	}
	if (caseScores.size() < 5 || ctrlScores.size() < 5) {
		throw ValueError("Sample size too small to perform Mann-Whitney test.");
	}

	if (m_store) {
		if (!d.hasVar("RankStats")) {
			matrixf initstats(m_tailed);
			d.setVar("RankStats", initstats);
		}
		matrixf & wstats = d.getMatrixVar("RankStats");
		if (wstats.size() == 1) {
			// one-sided test
			wstats[0].push_back(Mann_Whitneyu(caseScores, ctrlScores));
			d.setStatistic(wstats[0]);

		} else {
			if (wstats[0].size() <= wstats[1].size()) {
				// should be testing model 1
				wstats[0].push_back(Mann_Whitneyu(caseScores, ctrlScores));
				d.setStatistic(wstats[0]);
			} else {
				// should be testing model 2
				// rank sum of scores in ctrls
				wstats[1].push_back(Mann_Whitneyu(ctrlScores, caseScores));
				d.setStatistic(wstats[1]);
			}
		}
	}else {
		/*
		   if (!d.hasVar("RankStatsOrig")) {
		    vectorf initstats(m_tailed);
		    d.setVar("RankStatsOrig", initstats);
		   }
		   vectorf & wstat = d.getArrayVar("RankStatsOrig");
		   if (wstat.size() == 1) {
		    // rank sum of scores in cases
		    wstat[0] = (Mann_Whitneyu(caseScores, ctrlScores));
		    d.setStatistic(wstat[0]);
		   } else {
		    wstat[0] = (Mann_Whitneyu(caseScores, ctrlScores));
		    wstat[1] = (Mann_Whitneyu(ctrlScores, caseScores));
		    d.setStatistic(wstats[1]);
		   }
		 */
		// rank sum of scores in cases
		d.setStatistic(Mann_Whitneyu(caseScores, ctrlScores));
	}
	return true;
}


bool WSSPvalue::apply(AssoData & d, int timeout)
{
	if (!d.hasVar("RankStats")) {
		throw ValueError("Cannot find Mann-Whitney test statistic");
	}

	matrixf & mwstats = d.getMatrixVar("RankStats");
	if (mwstats.size() != m_tailed) {
		throw ValueError("Rank Statistic does not match the alternative");
	}
	vectorf wsstat(0);
	//compute mean and se. skip the first element
	//which is the original statistic
	for (size_t s = 0; s < m_tailed; ++s) {
		double ntotal = mwstats[s].size() - 1.0;
		if (!(ntotal > 1.0)) {
			throw ValueError("Rank Statistic not saved");
		}
		double mean_stats = (double)std::accumulate(mwstats[s].begin() + 1, mwstats[s].end(), 0.0)
		                    / ntotal;
		double se_stats = 0.0;
		for (size_t i = 1; i < mwstats[s].size(); ++i) {
			se_stats += pow(mwstats[s][i] - mean_stats, 2.0);
		}
		se_stats = sqrt(se_stats / (ntotal - 1.0));
		// FIXME prevent extreme values: any better approach?
		if (fEqual(se_stats, 0.0)) se_stats = 1.0e-6;
		wsstat.push_back((mwstats[s][0] - mean_stats) / se_stats);
	}
	if (m_tailed == 1) {
		// one-sided
		d.setPvalue(gsl_cdf_ugaussian_Q(wsstat.front()));
	}else {
		// two-sided (?) FIXME
		double pval = 0.0;
		if (!fEqual(wsstat[0], wsstat[1])) {
			pval = std::min(1.0, std::min(gsl_cdf_ugaussian_Q(wsstat[0]), gsl_cdf_ugaussian_Q(wsstat[1])) * 2.0);
		} else {
			pval = gsl_cdf_ugaussian_Q(wsstat[0]);
		}
		d.setPvalue(pval);
	}
	return true;
}


bool FindGenotypePattern::apply(AssoData & d, int timeout)
{
	//!-Compute unique genotype patterns (string) as ID scores (double)
	d.setGenotypeId();
	vectorf & genotypeId = d.genotype_id();
	// unique genotype patterns
	vectorf uniquePattern = genotypeId;
	std::sort(uniquePattern.begin(), uniquePattern.end());
	std::vector<double>::iterator it = std::unique(uniquePattern.begin(), uniquePattern.end());
	uniquePattern.resize(it - uniquePattern.begin());
	if (fEqual(uniquePattern.front(), 0.0)) uniquePattern.erase(uniquePattern.begin());
	if (uniquePattern.size() == 0) {
		throw ValueError("Input genotype matrix does not have a variant");
	}
	// count number of sample individuals for each genotype pattern
	vectori uniquePatternCounts(uniquePattern.size(), 0);
	for (size_t i = 0; i < genotypeId.size(); ++i) {
		// for each sample, identify/count its genotype pattern
		for (size_t u = 0; u < uniquePattern.size(); ++u) {
			if (genotypeId[i] == uniquePattern[u]) {
				// genotype pattern identified
				++uniquePatternCounts[u];
				// count this genotype pattern
				break;
			} else ;
			// genotype pattern not found -- move on to next pattern
		}
	}
	d.setVar("uniqGPattern", uniquePattern);
	d.setVar("uniqGCounts", uniquePatternCounts);
	return true;
}


bool KBACtest::apply(AssoData & d, int timeout)
{
	vectorf & up = d.getArrayVar("uniqGPattern");
	vectori & upc = d.getIntArrayVar("uniqGCounts");
	vectorf & gId = d.genotype_id();
	vectorf & ydat = d.phenotype();
	unsigned nCases = d.getIntVar("ncases");
	unsigned nCtrls = d.getIntVar("nctrls");
	double ybar = d.getDoubleVar("ybar");
	//
	vectorf kbac(m_tailed);
	// KBAC weights (unique genotype pattern weights)
	matrixf upWeights(m_tailed);

	for (size_t s = 0; s < m_tailed; ++s) {
		// genotype pattern counts in cases (or in ctrls)
		vectori upcSub(up.size(), 0);
		for (size_t i = 0; i < ydat.size(); ++i) {
			if ((s) ? (ydat[i] <= ybar) : (ydat[i] > ybar)) {
				// identify/count genotype pattern in cases (or ctrls for two sided test)
				for (size_t u = 0; u < up.size(); ++u) {
					if (gId[i] == up[u]) {
						// genotype pattern identified/counted in cases (or ctrls for two sided test)
						++upcSub[u];
						break;
					}else ;
					// genotype pattern not found -- move on to next pattern
				}
			}else ;
		}


		// genotype pattern weights, the hypergeometric distribution cmf
		upWeights[s].resize(up.size());
		for (size_t u = 0; u < up.size(); ++u) {
			upWeights[s][u] = gsl_cdf_hypergeometric_P(
				upcSub[u],
				upc[u],
				ydat.size() - upc[u],
				(s) ? nCtrls : nCases
			    );
		}
		if (!m_weightOnly) {
			// KBAC statistic: sum of genotype pattern frequencies differences in cases vs. controls
			// weighted by the hypergeometric distribution kernel
			kbac[s] = 0.0;
			for (size_t u = 0; u < up.size(); ++u) {
				kbac[s] += upWeights[s][u] * ( (1.0 * upcSub[u]) / (1.0 * ((s) ? nCtrls : nCases))
				                              - (1.0 * (upc[u] - upcSub[u])) / (1.0 * ((s) ? nCases : nCtrls)) );
			}
		}
	}
	if (m_weightOnly) {
		d.setVar("KBACweight", upWeights);
	} else {
		if (m_tailed == 1) d.setStatistic(kbac[0]);
		else d.setStatistic(std::max(kbac[0], kbac[1]));
	}

	return true;
}


bool RBTtest::apply(AssoData & d, int timeout)
{
	vectorf & ydat = d.phenotype();
	matrixf & xdat = d.raw_genotype();
	unsigned nCases = d.getIntVar("ncases");
	unsigned nCtrls = d.getIntVar("nctrls");
	double ybar = d.getDoubleVar("ybar");
	matrixf RBTweights(m_tailed);

	for (size_t j = 0; j < xdat.front().size(); ++j) {

		//! - Count number of variants in cases/controls at a locus
		unsigned countcs = 0;
		unsigned countcn = 0;
		for (size_t i = 0; i < ydat.size(); ++i) {
			if (ydat[i] <= ybar) {
				countcn += (unsigned)xdat[i][j];
			} else {
				countcs += (unsigned)xdat[i][j];
			}
		}

		//float w;
		//k0=(int)(freq*2*nCtrls);
		if (countcs > 0 && (1.0 * countcs) / (1.0 * nCases) > (1.0 * countcn) / (1.0 * nCtrls)) {
			double f = gsl_cdf_poisson_P(countcn, nCtrls * (countcs + countcn) / (1.0 * ydat.size()))
			           * (1.0 - gsl_cdf_poisson_P(countcs - 1, nCases * (countcs + countcn) / (1.0 * ydat.size())));
			RBTweights[0].push_back(-log(f));
		} else {
			RBTweights[0].push_back(0.0);
		}
		if (m_tailed > 1) {
			//k0=(int)(freq*2*nCases);
			if (countcn > 0 && (1.0 * countcn) / (1.0 * nCtrls) > (1.0 * countcs) / (1.0 * nCases)) {
				double f = gsl_cdf_poisson_P(countcs, nCases * (countcs + countcn) / (1.0 * ydat.size()))
				           * (1.0 - gsl_cdf_poisson_P(countcn - 1, nCtrls * (countcs + countcn) / (1.0 * ydat.size())));
				RBTweights[1].push_back(-log(f));
			} else {
				RBTweights[1].push_back(0.0);
			}
		}
	}
	if (m_weightOnly) {
		d.setVar("RBTweight", RBTweights);
	} else {
		// RVP statistic: The max statistic in the manuscript
		// R - potentially risk and P - potentially protective
		double sumR = std::accumulate(RBTweights[0].begin(), RBTweights[0].end(), 0.0);
		if (m_tailed == 1) {
			d.setStatistic(sumR);
		} else {
			double sumP = std::accumulate(RBTweights[1].begin(), RBTweights[1].end(), 0.0);
			d.setStatistic(std::max(sumR, sumP));
		}
	}
	return true;
}


bool AdaptiveRvSum::apply(AssoData & d, int timeout)
{
	vectorf & ydat = d.phenotype();
	matrixf & xdat = d.raw_genotype();
	unsigned nCases = d.getIntVar("ncases");
	unsigned nCtrls = d.getIntVar("nctrls");
	double multiplier = (d.getIntVar("moi") > 0) ? d.getIntVar("moi") * 1.0 : 1.0;
	double ybar = d.getDoubleVar("ybar");
	//a binary sequence for whether or not there is an excess of rare variants in controls
	vectori cexcess(xdat.front().size(), 0);

	for (size_t j = 0; j < xdat.front().size(); ++j) {
		double tmp1 = 0.0;
		double tmp0 = 0.0;
		for (size_t i = 0; i < xdat.size(); ++i) {
			if (ydat[i] <= ybar) {
				//FIXME might want to use > 0.0 to handle imputed genotypes
				if (xdat[i][j] >= 1.0) {
					tmp0 += xdat[i][j];
				}
			}else {
				if (xdat[i][j] >= 1.0) {
					tmp1 += xdat[i][j];
				}
			}
		}
		if (tmp0 / (nCtrls * multiplier) > tmp1 / (nCases * multiplier)) {
			cexcess[j] = 1;
		}else {
			continue;
		}
	}

	//a binary sequence for whether or not a site should be recoded
	vectori rsites(xdat.front().size(), 0);

	for (size_t j = 0; j < xdat.front().size(); ++j) {
		if (cexcess[j] == 0) continue;

		// 2X2 Fisher's test onesided, no midp
		vectori twotwoTable(4, 0);
		for (size_t i = 0; i < xdat.size(); ++i) {
			double gtmp = (xdat[i][j] == xdat[i][j]) ? xdat[i][j] : 0.0;
			// note the difference here from Fisher2X2 class
			// since in this case we check if rare variants are siginificantly enriched in ctrls
			if (ydat[i] <= ybar) {
				if (gtmp > 0.0) twotwoTable[0] += 1;
				else twotwoTable[1] += 1;
			}else {
				if (gtmp > 0.0) twotwoTable[2] += 1;
				else twotwoTable[3] += 1;
			}
		}

		double pvalue = 1.0 - gsl_cdf_hypergeometric_P(
			(twotwoTable[0] - 1),
			(twotwoTable[0] + twotwoTable[2]),
			(twotwoTable[1] + twotwoTable[3]),
			(twotwoTable[0] + twotwoTable[1])
		    );

		if (pvalue < 0.1) {
			rsites[j] = 1;
		}else {
			continue;
		}
	}

	vectorf & X = d.genotype();
	X.resize(xdat.size());
	std::fill(X.begin(), X.end(), 0.0);
	//recode sites
	for (size_t i = 0; i < X.size(); ++i) {
		//  scan all sample individuals
		for (size_t j = 0; j < xdat.front().size(); ++j) {
			//!- Recode protective variants as in Pan 2010
			if (xdat[i][j] == xdat[i][j]) {
				X[i] += (rsites[j]) ? (multiplier - xdat[i][j]) : xdat[i][j];
			}
		}
	}

	d.setVar("xbar", (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size()));
	return true;
}


bool FindVariantPattern::apply(AssoData & d, int timeout)
{
	matrixf & xdat = d.raw_genotype();
	vectori vnum(xdat.front().size(), 0);

	for (size_t j = 0; j < xdat.front().size(); ++j) {
		for (size_t i = 0; i < xdat.size(); ++i) {
			if (xdat[i][j] == xdat[i][j]) {
				vnum[j] += (int)xdat[i][j];
			}
		}
	}
	d.setVar("allVPattern", vnum);
	// unique genotype patterns
	std::sort(vnum.begin(), vnum.end());
	std::vector<int>::iterator it = std::unique(vnum.begin(), vnum.end());
	vnum.resize(it - vnum.begin());
	if ((vnum.front() == 0 && vnum.size() == 1) || vnum.size() == 0) {
		throw ValueError("Input genotype matrix does not have a variant");
	}
	if (vnum.front() != 0) {
		vnum.insert(vnum.begin(), 0);
	}
	d.setVar("uniqVPattern", vnum);
	return true;
}


bool VTTest::apply(AssoData & d, int timeout)
{
	// VT method as in Price et al 2010
	//! - Define <b> 'allZs' </b>, a vector of the Z scores computed under different thresholds
	matrixf & xdat = d.raw_genotype();
	vectorf & ydat = d.phenotype();
	double ybar = d.getDoubleVar("ybar");
	vectori & vnum = d.getIntArrayVar("allVPattern");
	vectori & uniqv = d.getIntArrayVar("uniqVPattern");

	vectorf allZs(uniqv.size() - 1);
	vectorf zIa(xdat.size(), 0.0), zIb(xdat.size(), 0.0);

	//! - Iterate the following for all thresholds:
	for (size_t t = 1; t < uniqv.size(); ++t) {
		//! - - Record the index of loci that can be added at the new threshold criteria
		vectori idxesAdding(0);
		for (size_t s = 0; s < vnum.size(); ++s) {
			if (vnum[s] > uniqv[t - 1] && vnum[s] <= uniqv[t])
				idxesAdding.push_back(s);
		}

		//!- - For loci passing the threshold criteria, implement Price paper page 3 z(T) formula
		for (size_t i = 0; i < xdat.size(); ++i) {
			double zIai = 0.0, zIbi = 0.0;
			for (size_t j = 0; j < idxesAdding.size(); ++j) {
				size_t locIdx = (size_t)idxesAdding[j];
				zIai += (xdat[i][locIdx] > 0.0) ? (ydat[i] - ybar) * xdat[i][locIdx] : 0.0;
				zIbi += (xdat[i][locIdx] > 0.0) ? xdat[i][locIdx] * xdat[i][locIdx] : 0.0;
			}
			zIa[i] += zIai;
			zIb[i] += zIbi;
		}
		//! - Now each element in <b> 'allZs' </b> is z(T) for different thresholds
		double sb = std::accumulate(zIb.begin(), zIb.end(), 0.0);
		if (fEqual(sb, 0.0)) allZs[t - 1] = 0.0;
		else allZs[t - 1] = std::accumulate(zIa.begin(), zIa.end(), 0.0) / sqrt(sb);
	}
	//! - Compute zmax, the statistic; square it for two-sided test
	if (m_tailed == 2) {
		for (size_t i = 0; i < allZs.size(); ++i) {
			allZs[i] = allZs[i] * allZs[i];
		}
	}
	d.setStatistic(*max_element(allZs.begin(), allZs.end()));

	return true;
}


bool VTFisher::apply(AssoData & d, int timeout)
{
	// VT method using Fisher's test
	// will send out a stopping signal if the initial statistics are significant
	//! - Define <b> 'allPs' </b>, a vector of the p-values computed under different thresholds
	matrixf & xdat = d.raw_genotype();
	vectorf & ydat = d.phenotype();
	vectori & vnum = d.getIntArrayVar("allVPattern");
	vectori & uniqv = d.getIntArrayVar("uniqVPattern");

	vectorf allPs(uniqv.size() - 1);
	vectorf allStats(uniqv.size() - 1);
	vectorf XCurr(xdat.size(), 0.0);
	bool shouldstop = true;

	//! - Iterate the following for all thresholds:
	for (size_t t = 1; t < uniqv.size(); ++t) {
		//! - - Record the index of loci that can be added at the new threshold criteria
		vectori idxesAdding(0);
		for (size_t s = 0; s < vnum.size(); ++s) {
			if (vnum[s] > uniqv[t - 1] && vnum[s] <= uniqv[t])
				idxesAdding.push_back(s);
		}

		//!- - For loci passing the threshold, implement CMC Fisher
		vectorf X(xdat.size(), 0.0);
		for (size_t i = 0; i < X.size(); ++i) {
			for (size_t j = 0; j < idxesAdding.size(); ++j) {
				size_t locIdx = (size_t)idxesAdding[j];
				X[i] = (xdat[i][locIdx] > 0.0) ? 1.0 : 0.0;
			}
			XCurr[i] = (fEqual(X[i], 1.0)) ? 1.0 : XCurr[i];
		}
		// fisher's test
		vectori twotwoTable(4, 0);
		for (size_t i = 0; i < XCurr.size(); ++i) {

			if (!(fEqual(ydat[i], 1.0) || fEqual(ydat[i], 0.0))) {
				throw ValueError("Input ydat data not binary");
			}

			if (fEqual(ydat[i], 1.0)) {
				if (XCurr[i] > 0.0)
					twotwoTable[0] += 1;
				else
					twotwoTable[1] += 1;
			}else {
				if (XCurr[i] > 0.0)
					twotwoTable[2] += 1;
				else
					twotwoTable[3] += 1;
			}
		}


		double pvalue = 0.0;
		if (m_tailed == 1) {
			pvalue = (m_midp)
			         ? (twotwoTable[0] > 0) * gsl_cdf_hypergeometric_P(
				(twotwoTable[0] - 1),
				(twotwoTable[0] + twotwoTable[2]),
				(twotwoTable[1] + twotwoTable[3]),
				(twotwoTable[0] + twotwoTable[1])
			    ) + 0.5 * gsl_ran_hypergeometric_pdf(
				twotwoTable[0],
				(twotwoTable[0] + twotwoTable[2]),
				(twotwoTable[1] + twotwoTable[3]),
				(twotwoTable[0] + twotwoTable[1])
			    ) : gsl_cdf_hypergeometric_P(
				(twotwoTable[0] - 1),
				(twotwoTable[0] + twotwoTable[2]),
				(twotwoTable[1] + twotwoTable[3]),
				(twotwoTable[0] + twotwoTable[1])
			    );
			pvalue = 1.0 - pvalue;
		} else {
			pvalue = fexact2x2(twotwoTable, "two");
		}
		// end of fisher's test

		if (pvalue <= m_alpha) {
			// some tests can be siginficant.
			shouldstop = false;
		}
		allPs[t - 1] = pvalue;
		allStats[t - 1] = (twotwoTable[1] * twotwoTable[2]) ? (double)twotwoTable[3] * (double)twotwoTable[0] / ((double)twotwoTable[1] * (double)twotwoTable[2]) : DBL_MAX;
	}
	// with this adaptive Fisher's approach the resulting p-values will not be uniformly distributed
	vectorf validPs(0), validStats(0);
	for (size_t i = 0; i < allPs.size(); ++i) {
		if (allPs[i] > 0.0 && allPs[i] < 1.0) validPs.push_back(allPs[i]);
		if (std::abs(allStats[i]) < DBL_MAX) validStats.push_back(allStats[i]);
	}
	if (validStats.size() == 0) {
		d.setStatistic(0.0);
	} else {
		d.setStatistic(*max_element(validStats.begin(), validStats.end()));
	}
	if (shouldstop) {
		// set the final p-value. This will prevent carrying out permutation tests
		if (validPs.size() == 0) d.setPvalue(0.5);
		else d.setPvalue(*min_element(validPs.begin(), validPs.end()));
	}

	return true;
}


bool CalphaTest::apply(AssoData & d, int timeout)
{

	/*! * the c-alpha Statistic: sum of the std. error of variant counts in cases <br>
	 * *  Two-sided test, claims to be robust against protective mutations. <br>
	 * * See <em> Ben. Neale et al. 2011 PLoS Genet. </em> <br><br>
	 * * Implementation:
	 */
	vectorf & ydat = d.phenotype();
	matrixf & xdat = d.raw_genotype();
	double phat = d.getDoubleVar("ybar");
	//
	double calpT = 0.0;
	double calpV = 0.0;
	int singletonAll = 0;
	int singletonCases = 0;
	bool isEmptyData = true;

	for (size_t j = 0; j < xdat.front().size(); ++j) {

		//! - Count number of variants in cases/controls at a locus
		int countcs = 0;
		int countcn = 0;
		for (size_t i = 0; i < ydat.size(); ++i) {
			if (ydat[i] <= phat) countcn += (xdat[i][j] > 0.0) ? (int)xdat[i][j] : 0;
			else countcs += (xdat[i][j] > 0.0) ? (int)xdat[i][j] : 0;
		}
		//! - the c-alpha method implementation
		int ni = countcs + countcn;
		if (ni < 2) {
			singletonAll += ni;
			singletonCases += countcs;
			continue;
		}else {
			isEmptyData = false;
		}
		//!- * skip singletons
		double niv = ni * phat * (1 - phat);
		calpT += (countcs - ni * phat) * (countcs - ni * phat) - niv;
		for (int u = 0; u <= ni; ++u) {
			double tmess = (u - ni * phat) * (u - ni * phat) - niv;
			calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, ni);
		}
	}

	if (singletonAll >= 2) {
		isEmptyData = false;
		double niv = singletonAll * phat * (1 - phat);
		calpT += (singletonCases - singletonAll * phat) * (singletonCases - singletonAll * phat) - niv;
		for (int u = 0; u <= singletonAll; ++u) {
			double tmess = (u - singletonAll * phat) * (u - singletonAll * phat) - niv;
			calpV += tmess * tmess * gsl_ran_binomial_pdf(u, phat, singletonAll);
		}
	} else ;
	//!- * bin singletons
	if (isEmptyData) {
		throw ValueError("Cannot perform c-alpha test on data with all singletons");
	}
	//!- c-alpha statistic
	//reject the null hypothesis when Z is larger than expected
	//using a one-tailed standard normal distribution for reference
	if (m_permutation) {
		d.setStatistic(calpT / sqrt(calpV));
	} else {
		d.setStatistic(calpT);
		d.setSE(sqrt(calpV));
	}
	return true;
}


bool RareCoverTest::apply(AssoData & d, int timeout)
{
	// RareCover method, 2010 PLoS CompBio
	vectorf & ydat = d.phenotype();
	matrixf & xdat = d.raw_genotype();

	if (!d.hasVar("polymophic_index")) {
		vectorf & mafs = d.getArrayVar("maf");
		//!- Index of loci that are observed to be polymophic
		vectori vntVct(0);
		for (size_t j = 0; j < mafs.size(); ++j) {
			if (mafs[j] > 0.0) vntVct.push_back(j + 1);
		}
		if (vntVct.size() == 0) {
			throw ValueError("Input genotype matrix does not have a variant site");
		} else {
			d.setVar("polymophic_index", vntVct);
		}
	}
	vectori vntNow = d.getIntArrayVar("polymophic_index");
	//!- the current and the next statistics
	double sCurr = 0.0, sNext = 0.0;
	vectorf X(ydat.size(), 0.0);
	//!- the current and the next genotype coding
	vectorf XCurr(ydat.size(), 0.0);
	do {
		sCurr = sNext;
		unsigned rmIdx = 0;
		//!- the "test contributing" variant index, for the vntVct object
		bool rmIdxFlag = false;
		for (size_t t = 0; t < vntNow.size(); ++t) {

			if (vntNow[t] == 0) continue;
			size_t iIdx = vntNow[t] - 1;
			//!- the index of a variant site
			for (size_t i = 0; i < ydat.size(); ++i) {
				X[i] = (XCurr[i] + xdat[i][iIdx] > 0.0) ? 1.0 : 0.0;
			}
			double statistic = chisq2X2stat(X, ydat);

			if (statistic > sNext) {
				sNext = statistic;
				rmIdx = t;
				rmIdxFlag = true;
			} else continue;
		}
		//!- Now end up with a properly updated sNext statistic.
		if (rmIdxFlag == false) {
			// in this case, sNext is not updated, OK to break the loop.
			break;
		} else {
			size_t rmVnt = vntNow[rmIdx] - 1;
			//!- Update the genotype coding by adding in the contributing locus
			for (size_t i = 0; i < ydat.size(); ++i) {
				XCurr[i] = XCurr[i] + xdat[i][rmVnt];
			}
			//!- remove the contributing locus to avoid duplicated visit to it the next time.
			vntNow[rmIdx] = 0;
		}
	} while (sNext - sCurr > m_difQ);

	//!- Test statistic 'sNext'
	d.setStatistic(sNext);
	return true;

}


//////////////


bool PyAction::apply(AssoData & d, int timeout)
{
	// Passing d to the function
	PyObject * args = PyTuple_New(m_func.numArgs());

	for (size_t i = 0; i < m_func.numArgs(); ++i) {
		const string & arg = m_func.arg(i);
		if (arg == "data")
			PyTuple_SET_ITEM(args, i, pyAssoDataObj(&d));
		else
			throw ValueError("Callback function for action PyAction only accept parameter data:");
	}
	PyObject * res = m_func(args);
	return PyObject_IsTrue(res) == 1;
}


double BasePermutator::getP(unsigned pcount1, unsigned pcount2, size_t current, unsigned alt) const
{
	// calculate permutation based p-value
	double x;

	if (alt == 1) {
		x = 1.0 + pcount1;
	} else {
		x = std::min(pcount1 + 1.0, pcount2 + 1.0);
	}
	double pval = x / (current + 1.0);
	return (alt == 1) ? pval : pval * 2.0;
}


double BasePermutator::check(unsigned pcount1, unsigned pcount2, size_t current, unsigned alt, double sig) const
{
	// the adaptive p-value technique
	if (current % 1000 != 0 || current == 0) {
		return 9.0;
	}
	double x;
	if (alt == 1) {
		x = 1.0 + pcount1;
	} else {
		x = std::min(pcount1 + 1.0, pcount2 + 1.0);
	}


	double n = current + 1.0;
	double alpha = 0.05;

	/*
	 * use 95CI for the adaptive procedure
	 * There are many methods for computing the 95CI for binomial random variables
	 * Discussions on CI see Alan AGRESTI and Brent A. COULL, 1998
	 * ==== OPTION1: Clopper-Pearson interval, conservative ====
	 * bci <- function(n, x, alpha) {
	 *  lower <- (1+(n-x+1)/(x*qf(alpha/2, 2*x, 2*(n-x+1))))^(-1)
	 *  upper <- (1+(n-x+1)/((x+1)*qf(1-alpha/2, 2*(x+1), 2*(n-x))))^(-1)
	 *  return(c(lower, upper))
	 *  }
	 * ==== OPTION2: the exact procedure, not usable now due to a bug in gsl_cdf_fdist_Pinv ====
	 * double plw = 1.0 / (1.0+(n-x+1.0)/(x*gsl_cdf_fdist_Pinv(alpha/2.0, 2.0*x, 2.0*(n-x+1.0))));
	 * ==== OPTION3: Edwin B. Wilson interval, not very useful because it is overly stringent ====
	 * wci <- function(n, x, alpha) {
	 *  pval <- x/n
	 *  z <- qnorm(1.0-alpha/2.0)
	 *  zsq <- z*z
	 *  lower <- (pval + zsq / (2.0*n) - z * sqrt((pval*(1.0-pval)+zsq/(4.0*n))/(1.0*n))) / (1.0+zsq/(1.0*n))
	 *  upper <- (pval + zsq / (2.0*n) + z * sqrt((pval*(1.0-pval)+zsq/(4.0*n))/(1.0*n))) / (1.0+zsq/(1.0*n))
	 *  return(c(lower, upper))
	 *  }
	 *  ==== OPTION4: simple Normal approximation interval ====
	 *  will use this in the current implementation
	 */

	double pval = x / n;
	double z = gsl_cdf_gaussian_Pinv(1.0 - alpha / 2.0, 1.0);
	// OPTION3 implementation:
	//double zsq = z * z;
	//double plw = (pval + zsq / (2.0 * n) - z * sqrt((pval * (1.0 - pval) + zsq / (4.0 * n)) / (1.0 * n))) / (1.0 + zsq / (1.0 * n));
	// OPTION4 implementation:
	double plw = pval - z * sqrt(pval * (1.0 - pval) / n);
	plw = (alt == 1) ? plw : plw * 2.0;

	if (plw > sig) {
		return (alt == 1) ? pval : pval * 2.0;
	} else {
		return 9.0;
	}
}


bool AssoAlgorithm::apply(AssoData & d, int timeout)
{
	Timer timer(timeout);

	for (size_t j = 0; j < m_actions.size(); ++j) {
		try {
			// an action can throw StopIteration to stop the rest of actions to be applied
			if (!m_actions[j]->apply(d, timer.remaining()))
				break;
		} catch (RuntimeError & e) {
			std::string msg = "Operator " + m_actions[j]->name() + " raises an exception (" + e.message() + ")";
			throw RuntimeError(msg);
		} catch (ValueError & e) {
			std::string msg = "Operator " + m_actions[j]->name() + " raises an exception (" + e.message() + ")";
			throw ValueError(msg);
		}
	}
	return true;
}


bool FixedPermutator::apply(AssoData & d, int timeout)
{
	if (d.pvalue().size()) {
		// p-value has already been calculated
		// this is for tests such as VT-Fisher
		// where the tester is able to determine if there is need to pursue the permutation procedure
		// after an initial calculation (will set p-value and statistic then to stop doing any permutation)
		if (d.pvalue().front() > 0.0 && d.pvalue().front() < 1.0) {
			return true;
		}
	}

	RNG rng;
	gsl_rng * gslr = rng.get();

	unsigned permcount1 = 0, permcount2 = 0;
	double pvalue = 9.0;
	// statistics[0]: statistic
	// statistics[1]: actual number of permutations
	vectorf statistics(2);

	// permutation timer initialization
	Timer timer(timeout);
	// Running statistic calculator
	RunningStat rs;
	// permutation loop
	for (size_t i = 0; i < m_times; ++i) {
		// apply actions to data
		for (size_t j = 0; j < m_actions.size(); ++j) {
			m_actions[j]->apply(d);
		}
		double statistic = d.statistic()[0];
		// set statistic or count for "success"
		if (i == 0) {
			statistics[0] = statistic;
			if (statistics[0] != statistics[0]) {
				d.setStatistic(std::numeric_limits<double>::quiet_NaN());
				d.setSE(std::numeric_limits<double>::quiet_NaN());
				d.setPvalue(std::numeric_limits<double>::quiet_NaN());
				d.setVar("NPERM", 0);
				return true;
			}
		} else {
			if (statistic > statistics[0]) {
				++permcount1;
			} else if (statistic < statistics[0]) {
				++permcount2;
			} else {
				if (gsl_rng_uniform(gslr) > 0.5) ++permcount1;
				else ++permcount2;
			}
			rs.Push(statistic);
		}
		// adaptive p-value calculation checkpoint
		if (m_sig < 1.0) {
			pvalue = check(permcount1, permcount2, i, m_alternative, m_sig);
		}
		// check for timeout and properly return
		// this event will be flagged by a negative p-value
		// which will be caught by python codes
		// and log a warning message for the event
		if (timer.timeout(i)) {
			pvalue = -1.0 * getP(permcount1, permcount2, i, m_alternative);
		}
		if (pvalue <= 1.0) {
			statistics[1] = double(i);
			break;
		}
		m_permute->apply(d);
	}

	// Permutation finished. Set p-value, statistic, std error, actual number of permutations, etc
	if (pvalue <= 1.0) {
		d.setPvalue(pvalue);
	} else{
		d.setPvalue(getP(permcount1, permcount2, m_times, m_alternative));
	}

	statistics[1] = (statistics[1] > 0.0) ? statistics[1] : double(m_times);
	d.setStatistic(statistics[0]);
	d.setSE(rs.StandardDeviation());
	d.setVar("NPERM", int(statistics[1]));
	return true;
}


bool VariablePermutator::apply(AssoData & d, int timeout)
{
	if (d.pvalue().size()) {
		// p-value has already been calculated
		// this is for tests such as VT-Fisher
		// where the tester is able to determine if there is need to pursue the permutation procedure
		// after an initial calculation (will set p-value and statistic then to stop doing any permutation)
		if (d.pvalue().front() > 0.0 && d.pvalue().front() < 1.0) {
			return true;
		}
	}
	if (!d.hasVar("mac"))
		throw RuntimeError("MAC has not been calculated. Please calculate MAC prior to using variable thresholds method.");
	// MAC thresholds, uniq, sorted
	// each element in this vector of MAC thresholds will define one subset of variant sites
	vectorf umac = d.getArrayVar("mac");
	// make sure umac has no decimal values
	for (size_t i = 0; i < umac.size(); ++i) {
		umac[i] = double(int(umac[i] + 0.5));
	}
	try {
		std::sort(umac.begin(), umac.end());
		std::vector<double>::iterator it = std::unique(umac.begin(), umac.end());
		umac.resize(it - umac.begin());
		if (umac.size() == 0) throw "mac is empty";
		if (fEqual(umac.front(), 0.0)) {
			umac.erase(umac.begin());
		}
		if (umac.size() == 0) throw "mac is empty";
	} catch (...) {
		// nothing to do
		d.setPvalue(std::numeric_limits<double>::quiet_NaN());
		d.setStatistic(std::numeric_limits<double>::quiet_NaN());
		d.setSE(std::numeric_limits<double>::quiet_NaN());
		return true;
	}
	// now umac should be a sorted vector of unique MAC's from AssoData
	double maclower = umac.front() - std::numeric_limits<double>::epsilon();

	// ==== optimization for certain methods ====
	// determine whether to use a quicker permutation routine if the actions are simply "codeX + doRegression"
	// (we have a general framework for variable thresholds, or VT, procedure but lost efficiency for specific methods)
	// (this "quick VT" is to bypass the general framework for certain simple VT procedure and implemented optimization for them)
	// with the optimization computation time can be reduced by 21.7%
	//

	unsigned choice = 0;

	if (m_actions.size() == 2) {
		if ((m_actions[0]->name() == "BinToX" ||
		     m_actions[0]->name() == "SumToX") &&
		    (m_actions[1]->name().find("Regression") != std::string::npos)
		    ) {
			choice = 1;
		}
	}

	matrixf genotypes(0);
	std::vector<size_t> gindex(0);

	if (choice) {
		// obtain genotype scores for subsets of variants (determined by MAC cut-offs)
		// and store them in "genotypes"
		AssoData * dtmp = d.clone();
		for (size_t m = 0; m < umac.size(); ++m) {
			SetSites(umac[(umac.size() - m - 1)], maclower, true).apply(*dtmp);
			// m_actions[0] is some coding theme, which will generate genotype scores
			m_actions[0]->apply(*dtmp);
			genotypes.push_back(dtmp->genotype());
		}
		delete dtmp;
		// keep an index of individuals
		// This will be used for genotype permutations
		// when only this gindex will be permuted
		// and all vectors in "genotypes" will be re-ordered by gindex
		for (size_t i = 0; i < genotypes.front().size(); ++i) {
			gindex.push_back(i);
		}
	}
	//
	// ==== END optimization for certain methods ====
	//
	// apply variable thresholds w/i permutation test
	unsigned permcount1 = 0, permcount2 = 0;
	// max_ for maximum over all statistics (for one-sided test)
	// min_ for maximum over all statistics (with max_, for two-sided test)
	double max_obstatistic = 0.0, min_obstatistic = 0.0;
	vectorf vt_obstatistic(0);
	double pvalue = 9.0;
	// statistics[0]: statistic
	// statistics[1]: actual number of permutations
	vectorf statistics(2);

	// permutation timer initialization
	Timer timer(timeout);
	// Running statistic calculator
	RunningStat rs_max;
	RunningStat rs_min;

	RNG rng;
	gsl_rng * gslr = rng.get();

	// permutation loop
	for (size_t i = 0; i < m_times; ++i) {
		vectorf vt_statistic(0);
		// make a copy of data
		// since the VT procedure will eliminate variant sites at each subsetting
		AssoData * dtmp = d.clone();
		if (choice) {
			// quick VT method as is originally implemented
			// loop over each vector of genotype scores (in "genotype" vector)
			// and calculate statistics and store them in vt_statistic
			for (size_t m = 0; m < genotypes.size(); ++m) {
				// shuffle genotype scores by gindex
				if (m_permute->name() != "PermuteY") {
					reorder(gindex.begin(), gindex.end(), genotypes[m].begin());
				}
				dtmp->setX(genotypes[m]);
				// m_actions[1] is some regression model fitting
				m_actions[1]->apply(*dtmp);
				vt_statistic.push_back(dtmp->statistic()[0]);
			}
		} else {
			// regular VT method
			// for each MAC threshold,
			// eliminate sites that are not confined in the thresholds
			// and apply actions on them
			for (size_t m = 0; m < umac.size(); ++m) {
				SetSites(umac[(umac.size() - m - 1)], maclower, true).apply(*dtmp);
				for (size_t j = 0; j < m_actions.size(); ++j) {
					m_actions[j]->apply(*dtmp);
				}
				vt_statistic.push_back(dtmp->statistic()[0]);
			}
		}
		delete dtmp;
		double max_statistic = *max_element(vt_statistic.begin(), vt_statistic.end());
		double min_statistic = *min_element(vt_statistic.begin(), vt_statistic.end());
		if (i == 0) {
			max_obstatistic = max_statistic;
			min_obstatistic = min_statistic;
			vt_obstatistic = vt_statistic;
			// if max_ or min_ is NA
			if (max_obstatistic != max_obstatistic) {
				d.setStatistic(std::numeric_limits<double>::quiet_NaN());
				d.setSE(std::numeric_limits<double>::quiet_NaN());
				d.setPvalue(std::numeric_limits<double>::quiet_NaN());
				d.setVar("VT_MAF", std::numeric_limits<double>::quiet_NaN());
				d.setVar("NPERM", 0);
				return true;
			}

		} else {
			if (max_statistic >= max_obstatistic && min_statistic <= min_obstatistic) {
				// both the max_ and min_ are "successful"
				// then randomly count the success to either permcount1 or permcount2
				if (gsl_rng_uniform(gslr) > 0.5) ++permcount1;
				else ++permcount2;
			} else {
				// either max_ or min_ is successful
				if (max_statistic >= max_obstatistic) {
					++permcount1;
				}
				if (min_statistic <= min_obstatistic) {
					++permcount2;
				}
			}
			rs_max.Push(max_statistic);
			rs_min.Push(min_statistic);
		}

		// adaptive p-value calculation checkpoint
		if (m_sig < 1.0) {
			pvalue = check(permcount1, permcount2, i, m_alternative, m_sig);
		}
		// check for timeout and properly return
		// this event will be flagged by a negative p-value
		// which will be caught by python codes
		// and log a warning message for the event
		if (timer.timeout(i)) {
			pvalue = -1.0 * getP(permcount1, permcount2, i, m_alternative);
		}
		if (pvalue <= 1.0) {
			statistics[1] = double(i);
			break;
		}
		// permutation
		// shuffle gindex for the "quick VT" procedure
		if (choice && m_permute->name() != "PermuteY") {
			random_shuffle(gindex.begin(), gindex.end());
		} else {
			m_permute->apply(d);
		}
	}
	// set p-value
	if (pvalue <= 1.0) {
		d.setPvalue(pvalue);
	} else {
		d.setPvalue(getP(permcount1, permcount2, m_times, m_alternative));
	}

	// set statistic, standard deviation and the MAF that gives the statistic
	if (m_alternative == 1) {
		statistics[0] = max_obstatistic;
		d.setSE(rs_max.StandardDeviation());
	} else {
		statistics[0] = (permcount1 >= permcount2) ? min_obstatistic : max_obstatistic;
		d.setSE((permcount1 >= permcount2) ? rs_min.StandardDeviation() : rs_max.StandardDeviation());
	}
	double multiplier = (d.getIntVar("moi") > 0) ? d.getIntVar("moi") * 1.0 : 1.0;
	for (size_t m = 0; m < umac.size(); ++m) {
		if (fEqual(vt_obstatistic[m], statistics[0])) {
			d.setVar("VT_MAF", umac[(umac.size() - m - 1)] / (d.samplecounts() * multiplier));
			break;
		}
	}
	// set standard error and number of actual permutations
	statistics[1] = (statistics[1] > 0.0) ? statistics[1] : double(m_times);
	d.setStatistic(statistics[0]);
	d.setVar("NPERM", int(statistics[1]));
	return true;
}


// foreach model:
// 1. get weight from input info (via d.getVar)
// 2. generate weighted sum (for wss and RBT) or recode genotype patterns by corresponding weight
// (will modify d.genotype() but not d.raw_genotype())
// 3. apply a test statistic
// end
// max(sm1, sm2)
bool WeightedGenotypeTester::apply(AssoData & d, int timeout)
{
	// check input actions
	// the first should be a weighting theme
	// the 2nd should be a regression method, or the rank test
	if (m_actions.size() != 2) {
		throw ValueError("Exactly 2 actions is allowed (for now) for WeightedGenotypeTester (would someone want a combination of more than one weighting theme?)");
	}
	std::string wtheme = m_actions[0]->name();
	if (wtheme.find("weight") == std::string::npos) {
		throw ValueError("Invalid input action " + wtheme);
	}
	if (m_actions[1]->name().find("Regression") == std::string::npos && m_actions[1]->name() != "MannWhitneyu") {
		throw ValueError("Invalid input action " + m_actions[1]->name());
	}
	// apply a weighting theme first
	m_actions[0]->apply(d);

	// get weight generated by the action
	vectorf & X = d.genotype();
	matrixf & genotype = d.raw_genotype();
	matrixf & weights = d.getMatrixVar(wtheme);

	// computer X for 1 or 2 sided hypothesis
	for (size_t s = 0; s < m_model; ++s) {
		X.resize(genotype.size());
		std::fill(X.begin(), X.end(), 0.0);
		if (wtheme.find("KBAC") == std::string::npos) {
			// weighting theme is not KBAC, so let's do a regular weighted sum coding
			// get external information first, if any
			matrixf extern_weights(0);
			if (m_info.size() > 0) {
				for (size_t w = 0; w < m_info.size(); ++w) {
					if (d.hasVar("__var_" + m_info[w])) {
						// get var_info
						extern_weights.push_back(d.getArrayVar("__var_" + m_info[w]));
					} else {
						throw ValueError("Cannot find genotype/variant information: " + m_info[w]);
					}

				}
			}
			// apply weighted sum
			for (size_t i = 0; i < genotype.size(); ++i) {
				//X[i] = std::accumulate(genotype[i].begin(), genotype[i].end(), 0.0);
				for (size_t j = 0; j < genotype[i].size(); ++j) {
					if (genotype[i][j] > 0.0) {
						// internal weight
						double gtmp = genotype[i][j] * (weights[s][j] == weights[s][j]) ? weights[s][j] : 1.0;
						// external weight
						for (size_t w = 0; w < extern_weights.size(); ++w) {
							gtmp *= (extern_weights[w][j] == extern_weights[w][j]) ? extern_weights[w][j] : 1.0;
						}
						X[i] += gtmp;
					}
				}
			}
		} else {
			// now kbac scoring system
			// using KBAC weights directly as X
			vectorf & uniquePattern = d.getArrayVar("uniqGPattern");
			vectorf & genotypeId = d.genotype_id();
			matrixf & uniquePatternWeights = d.getMatrixVar("KBACweight");
			for (size_t i = 0; i < genotype.size(); ++i) {
				for (size_t u = 0; u < uniquePattern.size(); ++u) {
					if (fEqual(genotypeId[i], uniquePattern[u])) {
						X[i] = uniquePatternWeights[s][u];
						break;
					}else ;
				}
			}
		}
		d.setVar("xbar", (double)std::accumulate(X.begin(), X.end(), 0.0) / (1.0 * X.size()));
		// association testing
		m_actions[1]->apply(d);
		// set statistic to m_stat
		m_stats[s] = d.statistic();
	}
	if (m_stats[1].size() == 0) {
		// only exists one vector, just use it
		d.setStatistic(m_stats[0]);
	} else {
		for (size_t i = 0; i < 2; ++i) m_stats[i][0] = std::abs(m_stats[i][0]);
		// compare to models, choose the larger of the two
		(m_stats[0][0] > m_stats[1][0]) ? d.setStatistic(m_stats[0]) : d.setStatistic(m_stats[1]);
	}
	return true;
}


}
