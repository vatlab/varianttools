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
	if (weight.size() != m_genotype.size() || weight.front() != m_genotype.front()) {
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


}
