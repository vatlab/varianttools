/**
 *  $File: action.h $
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

#ifndef _ACTION_H
#define _ACTION_H

#include "assoData.h"

namespace vtools {

/*
 * action classes
 * wrappers of assoData member methods that can be used in association.py
 *
 */
class BaseAction
{
public:
	BaseAction()
	{
	}


	virtual ~BaseAction()
	{
	}


	virtual BaseAction * clone() const
	{
		return new BaseAction(*this);
	}


	virtual bool apply(AssoData & d)
	{
		throw RuntimeError("The base action class should not be called");
		return true;
	}


	// return a string of the class name
	virtual std::string name()
	{
		return "BASEACTION";
	}


};


typedef std::vector<BaseAction *> vectora;


class SumToX : public BaseAction
{
public:
	SumToX() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new SumToX(*this);
	}


	bool apply(AssoData & d)
	{
		d.sumToX();
		d.meanOfX();
		return true;
	}


	std::string name()
	{
		return "SumToX";
	}


};


class BinToX : public BaseAction
{
public:
	BinToX() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new BinToX(*this);
	}


	bool apply(AssoData & d)
	{
		d.binToX();
		d.meanOfX();
		return true;
	}


	std::string name()
	{
		return "BinToX";
	}


};


class PermuteX : public BaseAction
{
public:
	PermuteX() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new PermuteX(*this);
	}


	bool apply(AssoData & d)
	{
		d.permuteX();
		return true;
	}


	std::string name()
	{
		return "PermuteX";
	}


};

class PermuteRawX : public BaseAction
{
public:
	PermuteRawX() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new PermuteRawX(*this);
	}


	bool apply(AssoData & d)
	{
		d.permuteRawX();
		return true;
	}


	std::string name()
	{
		return "PermuteRawX";
	}


};

class PermuteY : public BaseAction
{
public:
	PermuteY() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new PermuteY(*this);
	}


	bool apply(AssoData & d)
	{
		d.permuteY();
		return true;
	}


	std::string name()
	{
		return "PermuteY";
	}


};


class SetMaf : public BaseAction
{
public:
	SetMaf() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new SetMaf(*this);
	}


	bool apply(AssoData & d)
	{
		d.setMaf();
		return true;
	}


	std::string name()
	{
		return "SetMaf";
	}


};


class WeightByAllMaf : public BaseAction
{
	// this will change the raw genotypes directly!

public:
	WeightByAllMaf() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new WeightByAllMaf(*this);
	}


	bool apply(AssoData & d)
	{
		d.setMafWeight();
		d.weightX();
		return true;
	}


	std::string name()
	{
		return "WeightByAllMaf";
	}


};


class SetSites : public BaseAction
{
public:
	SetSites(double upper = 1.0, double lower = 0.0) :
		BaseAction(), m_upper(upper), m_lower(lower)
	{
	}


	BaseAction * clone() const
	{
		return new SetSites(*this);
	}


	bool apply(AssoData & d)
	{
		d.setSitesByMaf(m_upper, m_lower);
		return true;
	}


	std::string name()
	{
		return "SetSites";
	}


private:
	double m_upper;
	double m_lower;
};


class SimpleLinearRegression : public BaseAction
{
public:
	SimpleLinearRegression() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new SimpleLinearRegression(*this);
	}


	bool apply(AssoData & d)
	{
		d.simpleLinear();
		return true;
	}


	std::string name()
	{
		return "LinearRegression";
	}


};


class MultipleLinearRegression : public BaseAction
{
public:
	MultipleLinearRegression() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new MultipleLinearRegression(*this);
	}


	bool apply(AssoData & d)
	{
		d.multipleLinear();
		return true;
	}


	std::string name()
	{
		return "LinearRegression";
	}


};


class SimpleLogisticRegression : public BaseAction
{
public:
	SimpleLogisticRegression() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new SimpleLogisticRegression(*this);
	}


	bool apply(AssoData & d)
	{
		//assert(m_phenotype.size() == m_X.size());
		//check for binary
		d.simpleLogit();
		return true;
	}


	std::string name()
	{
		return "LogisticRegression";
	}


};


class GaussianPval : public BaseAction
{
public:
	GaussianPval(unsigned sided = 1)
		: BaseAction(), m_sided(sided)
	{
	}


	BaseAction * clone() const
	{
		return new GaussianPval(*this);
	}


	bool apply(AssoData & d)
	{
		d.gaussianP(m_sided);
		return true;
	}


	std::string name()
	{
		return "GaussianPval";
	}


private:
	unsigned m_sided;
};


class StudentPval : public BaseAction
{
public:
	StudentPval(unsigned sided = 1)
		: BaseAction(), m_sided(sided)
	{
	}


	BaseAction * clone() const
	{
		return new StudentPval(*this);
	}


	bool apply(AssoData & d)
	{
		d.studentP(m_sided);
		return true;
	}


	std::string name()
	{
		return "StudentPval";
	}


private:
	unsigned m_sided;
};


/** Permutator class
 *
 **/
class BasePermutator : public BaseAction
{
public:
	BasePermutator(const vectora & actions = vectora()) : BaseAction(), m_actions()
	{
		for (size_t i = 0; i < actions.size(); ++i)
			m_actions.push_back(actions[i]->clone());
	}


	~BasePermutator()
	{
		for (size_t i = 0; i < m_actions.size(); ++i)
			delete m_actions[i];
		m_actions.clear();
	}


	BasePermutator(const BasePermutator & rhs) : m_actions()
	{
		for (size_t i = 0; i < rhs.m_actions.size(); ++i)
			m_actions.push_back(rhs.m_actions[i]->clone());
	}


	BaseAction * clone() const
	{
		return new BasePermutator(*this);
	}


	std::string name()
	{
		return "BASEPERMUTATOR";
	}


	virtual bool apply(AssoData & d)
	{
		throw RuntimeError("The base permutation class should not be called");
		return true;
	}


	// append action(s) to the end of the action list
	void append(const BaseAction & action)
	{
		m_actions.push_back(action.clone());
	}


	void extend(const vectora & actions)
	{
		for (size_t i = 0; i < actions.size(); ++i)
			m_actions.push_back(actions[i]->clone());
	}


	// implementation of adaptive permutation
	// for every 1000 permutations this function will calculate a p-value
	// and check if its 95% confidence interval would capture the required significance level "sig"
	// will continue permutation if the required "sig" is captured by the 95% CI
	// otherwise will quit permutation and use this p-value as the final p-value to report
	double check(unsigned pcount1, unsigned pcount2, size_t current, unsigned alt, double sig) const;

protected:
	vectora m_actions;

};

// Action executor. simply execute a sequence of actions one by on AssoData object
class AssoAlgorithm : public BasePermutator
{

public:
	// algorithm with a single action (most likely a permutator)
	// in which case this wrapper is actually not needed.
	AssoAlgorithm(const BaseAction & action)
		: BasePermutator()
	{
		append(action);
	}


	// algorithm with a series of actions
	AssoAlgorithm(const vectora & actions)
		: BasePermutator(actions)
	{
	}


	BaseAction * clone() const
	{
		return new AssoAlgorithm(*this);
	}


	std::string name()
	{
		return "AssoAlgorithm";
	}


	bool apply(AssoData & d);

};

/* permutator class
 * a "fixed" set of variant sites will be involved in permutation test
 * i.e., will apply actions / permutation test on all variant sites in AssoData
 *
 * data members
 * m_times: number of permutations
 * m_alternative: 1 or 2, for 1-sided or 2-sided tests
 * m_sig: required significance level "alpha"
 * m_actions: a sequence of actions to be applied to AssoData
 */
class FixedPermutator : public BasePermutator
{

public:
	FixedPermutator(char pm, unsigned alternative, size_t times, double sig, const vectora & actions)
		: m_times(times), m_alternative(alternative), m_sig(sig), BasePermutator(actions)
	{
		// permute phenotypes or permute genotype scores
		m_permute = pm == 'Y' ? (BaseAction *)(new PermuteY()) : (BaseAction *)(new PermuteX());
	}


	~FixedPermutator()
	{
		delete m_permute;
	}


	FixedPermutator(const FixedPermutator & rhs) :
		BasePermutator(rhs),
		m_times(rhs.m_times), m_permute(rhs.m_permute->clone()),
		m_alternative(rhs.m_alternative), m_sig(rhs.m_sig)
	{
	}


	BaseAction * clone() const
	{
		return new FixedPermutator(*this);
	}


	std::string name()
	{
		return "FixedPermutator";
	}


	bool apply(AssoData & d);

private:
	size_t m_times;
	BaseAction * m_permute;
	unsigned m_alternative;
	double m_sig;
};

/* permutator class
 * a "variable" set of variant sites will be involved in permutation test
 * i.e., will apply actions on a number of subsets of variant sites
 * and use minimized p-value from these multiple tests
 * family-wise error rate is properly controlled in permutation procedure
 * currently, subsets of variant sites are defined by MAF on the sites
 *
 * data members
 * m_times: number of permutations
 * m_alternative: 1 or 2, for 1-sided or 2-sided tests
 * m_sig: required significance level "alpha"
 * m_actions: a sequence of actions to be applied to AssoData
 */
class VariablePermutator : public BasePermutator
{
	// permutator for variable thresholds methods

public:
	VariablePermutator(char pm, unsigned alternative, size_t times, double sig, const vectora & actions)
		: m_times(times), m_alternative(alternative), m_sig(sig), BasePermutator(actions)
	{
		// permute phenotypes or permute genotype scores
		m_permute = pm == 'Y' ? (BaseAction *)(new PermuteY()) : (BaseAction *)(new PermuteRawX());
	}


	~VariablePermutator()
	{
		delete m_permute;
	}


	VariablePermutator(const VariablePermutator & rhs) :
		BasePermutator(rhs),
		m_times(rhs.m_times), m_permute(rhs.m_permute->clone()),
		m_alternative(rhs.m_alternative), m_sig(rhs.m_sig)
	{
	}


	VariablePermutator * clone()
	{
		return new VariablePermutator(*this);
	}


	std::string name()
	{
		return "VariablePermutator";
	}


	bool apply(AssoData & d);

private:
	size_t m_times;
	BaseAction * m_permute;
	unsigned m_alternative;
	double m_sig;
};

}

#endif
