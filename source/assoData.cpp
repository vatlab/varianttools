/**
 *  $File: assoData.cpp $
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


#include "assoData.h"
#include <cmath>
#include <algorithm>
namespace vtools {


bool AssoData::setPhenotype(const vectorf & p)
{
	m_phenotype = p;
	// set phenotype mean
	setVar("ybar", (double)std::accumulate(m_phenotype.begin(), m_phenotype.end(), 0.0) / (1.0 * m_phenotype.size()));
	setVar("ncovar", 0);
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
	// re-size statistics vector, etc
	m_statistic.resize(c.size());
	m_se.resize(m_statistic.size());
	m_pval.resize(m_statistic.size());
	return true;
}


void AssoData::weightX(const vectorf & weight)
{
	if (weight.size() == 0) {
		return;
	}
	for (size_t i = 0; i < m_genotype.size(); ++i) {
		for (size_t j = 0; j < m_genotype[i].size(); ++j) {
			if (m_genotype[i][j] > 0) {
				m_genotype[i][j] *= (weight[j] == weight[j]) ? weight[j] : 1.0;
			}
		}
	}
}


void AssoData::weightX(const matrixf & weight)
{
	if (weight.size() != m_genotype.size() || weight.front().size() != m_genotype.front().size()) {
		throw ValueError("Genotype and genotype information data do not match in scale");
	}
	for (size_t i = 0; i < m_genotype.size(); ++i) {
		for (size_t j = 0; j < m_genotype[i].size(); ++j) {
			if (m_genotype[i][j] > 0) {
				m_genotype[i][j] *= (weight[i][j] == weight[i][j]) ? weight[i][j] : 1.0;
			}
		}
	}
}


bool AssoData::setGenotypeId()
{
	// note: this will simply ignore missing data (treat as wildtype)

	m_genotype_id.resize(m_phenotype.size());

	for (size_t i = 0; i < m_phenotype.size(); ++i) {

		double vntIdL = 0.0;
		double vntIdR = 0.0;
		const double ixiix = pow(9.0, 10.0);
		unsigned lastCnt = 0;
		unsigned tmpCnt = 0;

		for (size_t j = 0; j < m_genotype.front().size(); ++j) {

			if (m_genotype[i][j] >= 1.0) {
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * m_genotype[i][j];
			}else {
				continue;
			}
			if (vntIdR >= ixiix) {
				vntIdL = vntIdL + 1.0;
				vntIdR = vntIdR - ixiix;
				lastCnt = lastCnt + tmpCnt + 1;
				tmpCnt = 0;
				continue;
			}else {
				++tmpCnt;
				continue;
			}
		}
		// one-to-one "ID number" for a genotype pattern
		m_genotype_id[i] = vntIdL + vntIdR * 1e-10;
	}

	vectorf v = m_genotype_id;
	std::sort(v.begin(), v.end());
	std::vector<double>::iterator it = std::unique(v.begin(), v.end());
	v.resize(it - v.begin());
	for (size_t i = 0; i < m_genotype_id.size(); ++i) {
		std::vector<double>::iterator low = std::lower_bound(v.begin(), v.end(), m_genotype_id[i]);
		m_genotype_id[i] = double(low - v.begin());
	}
	return true;
}


}
