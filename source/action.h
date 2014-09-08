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

#include <limits>
#include <numeric>
#include <functional>
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


	virtual bool apply(AssoData & d, int timeout = 0)
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

// calculate genotype scores from raw genotype
// m_X = rowSum(m_genotype)
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


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "SumToX";
	}


};


// calculate genotype scores from raw genotype
// m_X = (rowSum(m_genotype) > 0)
// m_X is binary: 0 for all wildtype, 1 for having at least one mutation
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


	bool apply(AssoData & d, int timeout = 0);

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


	bool apply(AssoData & d, int timeout = 0)
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


	bool apply(AssoData & d, int timeout = 0)
	{
		d.permuteRawX();
		return true;
	}


	std::string name()
	{
		return "PermuteRawX";
	}


};

class PermuteAllX : public BaseAction
{
public:
	PermuteAllX() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new PermuteAllX(*this);
	}


	bool apply(AssoData & d, int timeout = 0)
	{
		d.permuteAllX();
		return true;
	}


	std::string name()
	{
		return "PermuteAllX";
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


	bool apply(AssoData & d, int timeout = 0)
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


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "SetMaf";
	}


};


// m_model: weight by maf from all, or selected samples (ctrls, low QT samples, etc)
// m_model == 0 : set weight by maf
// m_model == 1 : set weight by maf of samples having phenotype < mean(phenotypes), a 1Xn matrix
// m_model == 2 : set weight by maf of samples having phenotype >/< mean(phenotypes), a 2Xn matrix
// compute weight w = 1 / sqrt(p*(1-p))
class BrowningWeight : public BaseAction
{
public:
	BrowningWeight(unsigned model) :
		BaseAction(), m_model(model)
	{
	}


	BaseAction * clone() const
	{
		return new BrowningWeight(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "Browningweight";
	}


private:
	unsigned m_model;

};


class FillGMissing : public BaseAction
{
public:
	FillGMissing(std::string method = "maf") :
		BaseAction(), m_method(method)
	{
	}


	BaseAction * clone() const
	{
		return new FillGMissing(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "FillGMissing";
	}


private:
	std::string m_method;

};


class WeightByInfo : public BaseAction
{
	// this will change the raw genotypes directly!

public:
	WeightByInfo(const vectors & info) : BaseAction(), m_info(info)
	{
	}


	BaseAction * clone() const
	{
		return new WeightByInfo(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "WeightByInfo";
	}


private:
	vectors m_info;
};


// remove variant sites having MAF <= lower_bound or MAF > upper_bound
// if use_mac, then remove variant sites having MAC <= lower_bound or MAC > upper_bound
class SetSites : public BaseAction
{
public:
	SetSites(double upper = 1.0, double lower = 0.0, bool use_mac = false) :
		BaseAction(), m_upper(upper), m_lower(lower), m_use_mac(use_mac)
	{
	}


	BaseAction * clone() const
	{
		return new SetSites(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "SetSites";
	}


private:
	double m_upper;
	double m_lower;
	bool m_use_mac;
};


//recode genotypes to 0,1 for dominant or recessive coding
//actually from simulation using kyrukov's model I only see 16 homozygote alleles out of 1000 sample * 1000 replicates
class CodeXByMOI : public BaseAction
{
public:
	CodeXByMOI() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new CodeXByMOI(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "CodeXByMOI";
	}


};
// Wald's statistic for simple linear model Y = b0 + b1x
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


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "SimpleLinearRegression";
	}


};

//!- Score test implementation for logistic regression model logit(p) = b0 + b1x
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


	bool apply(AssoData & d, int timeout = 0);


	std::string name()
	{
		return "SimpleLogisticRegression";
	}


};


// fitting / calculating wald's statistic for multiple linear/logistic regression model
// Y = b0 + b1x1 + b2x2 + ... + bnxn
class MultipleRegression : public BaseAction
{
public:
	MultipleRegression(bool iSE = true, unsigned method = 0) :
		BaseAction(), m_iSE(iSE), m_method(method)
	{
	}


	BaseAction * clone() const
	{
		return new MultipleRegression(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "MultipleRegression";
	}


private:
	bool m_iSE;
	unsigned m_method;

	BaseLM * m_getModel()
	{
		switch (m_method) {
		case 0: return new LinearM;
		case 1: return new LogisticM;
		default: return NULL;
		}
	}


};


// get p-value from statistic, assuming standard normal distribution
class GaussianPval : public BaseAction
{
public:
	GaussianPval(unsigned sided = 1)
		: BaseAction(), m_tailed(sided)
	{
	}


	BaseAction * clone() const
	{
		return new GaussianPval(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "GaussianPval";
	}


private:
	unsigned m_tailed;
};


// get p-value from statistic, assuming t distribution
// degree of freedom will be calculated from data automatically
class StudentPval : public BaseAction
{
public:
	StudentPval(unsigned sided = 1)
		: BaseAction(), m_tailed(sided)
	{
	}


	BaseAction * clone() const
	{
		return new StudentPval(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "StudentPval";
	}


private:
	unsigned m_tailed;
};

/* Fisher's 2x2 table test on m_X and m_phenotype:

   m_X\phenotype	phen=1	phen=0
    m_X = 0			n1		n2
    m_X > 0			n3		n4

   Result is the pvalue for m_pvalue; and n3 will be m_statistic
   If midp option is activated will then use midp correction for one-sided test.
 */
class Fisher2X2 : public BaseAction
{
public:
	Fisher2X2(unsigned alternative, bool midp = false)
		: BaseAction(), m_tailed(alternative), m_midp(midp)
	{
	}


	BaseAction * clone() const
	{
		return new Fisher2X2(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "Fisher2X2";
	}


private:
	unsigned m_tailed;
	bool m_midp;
};

/* MannWhitney rank test on m_X and m_phenotype.
   store = false for not storing the statistics from each permutation test; true is otherwise
   will be stored as d.setVar("RankStats", mwstats)
 */
class MannWhitneyu : public BaseAction
{
public:
	MannWhitneyu(unsigned alternative = 1, bool store = false)
		: BaseAction(), m_tailed(alternative), m_store(store)
	{
	}


	BaseAction * clone() const
	{
		return new MannWhitneyu(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "MannWhitneyu";
	}


private:
	unsigned m_tailed;
	bool m_store;
};


class WSSPvalue : public BaseAction
{
public:
	WSSPvalue(unsigned alternative) :
		BaseAction(), m_tailed(alternative)
	{
	}


	BaseAction * clone() const
	{
		return new WSSPvalue(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "WSSPvalue";
	}


private:
	unsigned m_tailed;
};


// find uniq genotype patterns with their counts
// will modify data by:
// 1. replace missing data with most likely genotype
// 2. set d.setVar("gPattern", genotypeId); d.setVar("uniqGPattern", uniquePattern); d.setVar("uniqGCounts", uniquePatternCounts);
class FindGenotypePattern : public BaseAction
{
public:
	FindGenotypePattern() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new FindGenotypePattern(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "FindGenotypePattern";
	}


};


class KBACtest : public BaseAction
{
public:
	KBACtest(unsigned alternative = 1, bool weightOnly = false) :
		BaseAction(), m_tailed(alternative), m_weightOnly(weightOnly)
	{
	}


	BaseAction * clone() const
	{
		return new KBACtest(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		// use the more approperate name here
		return "KBACweight";
	}


private:
	unsigned m_tailed;
	bool m_weightOnly;
};


class RBTtest : public BaseAction
{
public:
	RBTtest(unsigned alternative = 1, bool weightOnly = false) :
		BaseAction(), m_tailed(alternative), m_weightOnly(weightOnly)
	{
	}


	BaseAction * clone() const
	{
		return new RBTtest(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "RBTweight";
	}


private:
	unsigned m_tailed;
	bool m_weightOnly;
};

// recode protective rare variants and sum up the recoded genotype
// will not change d.raw_genotype; will change d.genotype()
// determine sites of rare variants having excess copies in ctrls
// and determine whether or not to recode that site via a Fisher's test
// evaluated at alpha = 0.1
class AdaptiveRvSum : public BaseAction
{
public:
	AdaptiveRvSum() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new AdaptiveRvSum(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "AdaptiveRvSum";
	}


};


class FindVariantPattern : public BaseAction
{
public:
	FindVariantPattern() : BaseAction()
	{
	}


	BaseAction * clone() const
	{
		return new FindVariantPattern(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "FindVariantPattern";
	}


};

// The original statistic in Price 2010 AJHG
class VTTest : public BaseAction
{
public:
	VTTest(unsigned alternative = 1) :
		BaseAction(), m_tailed(alternative)
	{
	}


	BaseAction * clone() const
	{
		return new VTTest(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "VTTest";
	}


private:
	unsigned m_tailed;
};


class VTFisher : public BaseAction
{
public:
	VTFisher(double alpha, unsigned alternative = 1, bool midp = false) :
		BaseAction(), m_tailed(alternative), m_midp(midp), m_alpha(alpha)
	{
	}


	BaseAction * clone() const
	{
		return new VTFisher(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "VTFisher";
	}


private:
	unsigned m_tailed;
	bool m_midp;
	double m_alpha;
};


class CalphaTest : public BaseAction
{
public:
	CalphaTest(bool permutation = false) :
		BaseAction(), m_permutation(permutation)
	{
	}


	BaseAction * clone() const
	{
		return new CalphaTest(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "CalphaTest";
	}


private:
	bool m_permutation;

};


class RareCoverTest : public BaseAction
{
public:
	RareCoverTest(double difQ = 0.5) :
		BaseAction(), m_difQ(difQ)
	{
	}


	BaseAction * clone() const
	{
		return new RareCoverTest(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "RareCoverTest";
	}


private:
	//!- the cut-off to use for the "heuristic greedy algorithm". = 0.5 as suggested by the paper
	double m_difQ;
};


////////////////////////////////
////////////////////////////////

/**
 * Permutator class
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


	virtual bool apply(AssoData & d, int timeout = 0)
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

	// calculate p-value in permutation tests
	double getP(unsigned pcount1, unsigned pcount2, size_t current, unsigned alt) const;

protected:
	vectora m_actions;

};

/** This action accepts a user-provided function and will pass the AssoData
 *  object to this function when the action is applied.
 */
class PyAction : public BaseAction
{
public:
	PyAction(PyObject * func) : BaseAction(), m_func(func)
	{
	}


	PyAction(const PyAction & rhs) : m_func(rhs.m_func)
	{
	}


	BaseAction * clone() const
	{
		return new PyAction(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "Calling python function";
	}


private:
	PyFunc m_func;
};


// Action executor. simply execute a sequence of actions one by one on AssoData object
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


	bool apply(AssoData & d, int timeout = 0);

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
		: BasePermutator(actions), m_times(times), m_alternative(alternative), m_sig(sig)
	{
		// permute phenotypes or permute genotype scores
		m_permute = pm == 'Y' ? (BaseAction *)(new PermuteY()) : (BaseAction *)(new PermuteAllX());
	}


	~FixedPermutator()
	{
		delete m_permute;
	}


	FixedPermutator(const FixedPermutator & rhs) :
		BasePermutator(rhs),
		m_times(rhs.m_times), m_alternative(rhs.m_alternative),
		m_sig(rhs.m_sig), m_permute(rhs.m_permute->clone())
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


	bool apply(AssoData & d, int timeout = 0);

private:
	size_t m_times;
	unsigned m_alternative;
	double m_sig;
	BaseAction * m_permute;
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
		: BasePermutator(actions), m_times(times), m_alternative(alternative), m_sig(sig)
	{
		// permute phenotypes or permute genotype scores
		m_permute = pm == 'Y' ? (BaseAction *)(new PermuteY()) : (BaseAction *)(new PermuteAllX());
	}


	~VariablePermutator()
	{
		delete m_permute;
	}


	VariablePermutator(const VariablePermutator & rhs) :
		BasePermutator(rhs),
		m_times(rhs.m_times), m_alternative(rhs.m_alternative),
		m_sig(rhs.m_sig), m_permute(rhs.m_permute->clone())
	{
	}


	BaseAction * clone() const
	{
		return new VariablePermutator(*this);
	}


	std::string name()
	{
		return "VariablePermutator";
	}


	bool apply(AssoData & d, int timeout = 0);

private:
	size_t m_times;
	unsigned m_alternative;
	double m_sig;
	BaseAction * m_permute;
};

/* VAT Stacking permutator
 * m_times: number of permutations
 * m_sig: required significance level "alpha"
 * m_actions: a sequence of actions to be applied to AssoData
 */

class StackingPermutator : public BasePermutator
{

public:
	StackingPermutator(const vectora & actions, const vectori & isVt,
		char pm, size_t times, double sig)
		: BasePermutator(actions), m_isVt(isVt), m_times(times), m_sig(sig)
	{
		// permute phenotypes or permute genotype scores
		m_permute = pm == 'Y' ? (BaseAction *)(new PermuteY()) : (BaseAction *)(new PermuteAllX());
	}


	~StackingPermutator()
	{
		delete m_permute;
	}


	StackingPermutator(const StackingPermutator & rhs) :
		BasePermutator(rhs),
		m_isVt(rhs.m_isVt),
		m_times(rhs.m_times),
		m_sig(rhs.m_sig), m_permute(rhs.m_permute->clone())
	{
	}


	BaseAction * clone() const
	{
		return new StackingPermutator(*this);
	}


	std::string name()
	{
		return "StackingPermutator";
	}


	bool apply(AssoData & d, int timeout = 0);

private:
	vectori m_isVt;
	size_t m_times;
	double m_sig;
	BaseAction * m_permute;
};


// this class works within a permutator, to perform weighted sum test based on different weighting themes
// this class exists to handle the chanllege of using a weighting theme to test for both one-sided and two-sided hypothesis within a permutation test
// while preserving the genotype codings
// Python level input should be WeightedGenotypeTester(self.alternative, self.weight, [t.weighter(), t.test()])
// lenght of vector of input actions should be exactly two for now
class WeightedGenotypeTester : public BasePermutator
{

public:
	WeightedGenotypeTester(unsigned alternative, const vectors & info, const vectora & actions) :
		BasePermutator(actions), m_model(alternative), m_info(info)
	{
		vectorf tmp(0);

		m_stats.resize(2, tmp);
	}


	BaseAction * clone() const
	{
		return new WeightedGenotypeTester(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "WeightedGenotypeTester";
	}


private:
	unsigned m_model;
	vectors m_info;
	matrixf m_stats;
};

// this class applies different weighting themes to the same data in one permutation test
// and use the weight that maximizes the statistic
// logic similar to WeightedGenotypeTester class
class OptimalWeightTester : public BasePermutator
{

public:
	OptimalWeightTester(const vectors & info, const vectora & actions) :
		BasePermutator(actions), m_info(info)
	{
		vectorf tmp(0);

		m_stats.resize(m_info.size(), tmp);
	}


	BaseAction * clone() const
	{
		return new OptimalWeightTester(*this);
	}


	bool apply(AssoData & d, int timeout = 0);

	std::string name()
	{
		return "OptimalWeightTester";
	}


private:
	vectors m_info;
	matrixf m_stats;
};

}
#endif
