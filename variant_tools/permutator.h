/**
 *  $File: permutator.h $
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

#ifndef _PERMUTATOR_H
#define _PERMUTATOR_H

#include "assoConfig.h"
#include "assoData.h"
#include "action.h"

namespace vtools {


class BasePermutator
{
public:
	BasePermutator(const vectora & actions)
	{
		for (size_t i = 0; i < actions.size(); ++i)
			m_actions.push_back(actions[i]->clone());
	}


	~BasePermutator()
	{
		for (size_t i = 0; i < m_actions.size(); ++i)
			delete m_actions[i];
	}


	BasePermutator(const BasePermutator & rhs)
	{
		for (size_t i = 0; i < m_actions.size(); ++i)
			delete m_actions[i];
		for (size_t i = 0; i < rhs.m_actions.size(); ++i)
			m_actions.push_back(rhs.m_actions[i]->clone());
	}


	virtual double permute(AssoData & d)
	{
		//
        throw RuntimeError("The base permutation class should not be called");
	}


protected:
	vectora m_actions;
};


class PhenoPermutator : public BasePermutator
{
public:
	PhenoPermutator(size_t times, const vectora & actions)
		: m_times(times), BasePermutator(actions)
	{
	}


	double permute(AssoData & d)
	{
		double all_p = 0;

		for (size_t i = 0; i < m_times; ++i) {
			d.permuteY();
			for (size_t j = 0; j < m_actions.size(); ++j) {
				all_p += m_actions[j]->apply(d);
			}
		}
		return all_p;
	}


private:
	size_t m_times;
};


}

#endif
