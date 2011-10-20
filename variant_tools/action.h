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

#ifndef _ACTION_H
#define _ACTION_H

#include "assoData.h"

namespace vtools {


class BaseAction
{
public:
	BaseAction()
	{
	}


	virtual BaseAction * clone()
	{
		return new BaseAction(*this);
	}


	virtual double apply(AssoData & d)
	{
		return 0;
	}


};


class SumToX : public BaseAction
{
public:
	SumToX() : BaseAction()
	{
	}


	BaseAction * clone()
	{
		return new SumToX(*this);
	}


	double apply(AssoData & d)
	{
		d.sumToX();
		return 0;
	}


};


class SimpleLinearRegression : public BaseAction
{
public:
	SimpleLinearRegression() : BaseAction()
	{
	}


	BaseAction * clone()
	{
		return new SimpleLinearRegression(*this);
	}


	double apply(AssoData & d, int sided = 1)
  {    
    //assert(m_phenotype.size() == m_genotype.size());
    d.binToX();
    double xbar = d.mean_genotype();
    double ybar = d.mean_phenotype();
    double statistic = d.simpleLinear(xbar, ybar);
    return d.gaussianP(statistic, sided); 
  }


};


typedef std::vector<BaseAction *> vectora;


}

#endif
