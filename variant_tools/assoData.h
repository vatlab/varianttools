/**
 *  $File: assoData.h $
 *  $LastChangedDate: 2011-07-06 23:27:10 -0500 (Wed, 06 Jul 2011) $
 *  $Rev: 4256 $
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

#ifndef _ASSODATA_H
#define _ASSODATA_H

#include <vector>
typedef std::vector<double> vectorf;
typedef std::vector<int> vectori;
typedef std::vector<std::vector<int> > matrixi;

#include <numeric>

#include "assoConfig.h"

namespace vtools {


class AssoData
{
public:
	AssoData() : m_phenotype(0), m_genotype(0), m_X(0)
	{
	}


	// make a copy of the data
	virtual AssoData * clone() const
	{
		return new AssoData(*this);
	}


	// set data
	void setGenotype(const matrixi & g)
	{
		m_genotype = g;
	}


	void setPhenotype(const vectorf & p)
	{
		m_phenotype = p;
	}


	// return some information
	vectorf phenotype()
	{
		return m_phenotype;
	}


	vectori genotype()
	{
		return m_X;
	}


public:
	void permuteY()
	{
	}


	// manipulate data
	void sumToX()
	{
		m_X.resize(m_genotype.size());
		for (size_t i = 0; i < m_genotype.size(); ++i) {
			m_X[i] = std::accumulate(m_genotype[i].begin(), m_genotype[i].end(), 1.0);
		}
	}


private:
	/// raw phenotype and gneotype data
	vectorf m_phenotype;
	matrixi m_genotype;

	/// translated genotype
	vectori m_X;
};

}

#endif
