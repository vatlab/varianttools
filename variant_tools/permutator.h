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


      virtual double apply(AssoData & d)
      {
        throw RuntimeError("The base permutation class should not be called");
        return 0;
      }


    protected:
      vectora m_actions;
  };


  class ActionExecuter : public BasePermutator
  {

    public:
      ActionExecuter(const vectora & actions)
        : BasePermutator(actions)
      {
      }


      double apply(AssoData & d)
      {
        for (size_t j = 0; j < m_actions.size(); ++j) {
          m_actions[j]->apply(d);
        }
        return 0;
      }
  };


  class FixedPermutator : public BasePermutator
  {

    public:
      FixedPermutator(size_t times, char pm, unsigned alternative, const vectora & actions)
        : m_times(times), m_alternative(alternative), BasePermutator(actions)
      {
        if (pm == 'Y') {
          PermuteY* permute;
          m_permute = permute->clone();
        }
        else {
          PermuteX* permute;
          m_permute = permute->clone();
        }
      }

      ~FixedPermutator()
      {
        delete m_permute;
      }
      
      double apply(AssoData & d)
      {
        unsigned permcount1 = 0, permcount2 = 0;
        double obstatistic = 0.0;
        double pvalue = 9.0;

        for (size_t i = 0; i < m_times; ++i) {
          for (size_t j = 0; j < m_actions.size(); ++j) {
            m_actions[j]->apply(d);
          }
          double statistic = d.statistic(); 
          if (i == 0) {
            obstatistic = statistic;
          }
          else {
            if (statistic >= obstatistic) 
              ++permcount1;
            if (statistic <= obstatistic)
              ++permcount2;
          }
          if (pvalue <= 1.0) {
            break;
          }
          m_permute->apply(d);
        }

        if (pvalue <= 1.0) {
          d.setPvalue(pvalue); 
        }
        else {
          if (m_alternative == 1) { 
            pvalue = (permcount1 + 1.0) / (m_times + 1.0);
          }
          else {
            double permcount = fmin(permcount1, permcount2);
            pvalue = 2.0 * (permcount + 1.0) / (m_times + 1.0);
          }
          d.setPvalue(pvalue);
        } 

        return 0.0;
        //return (double) std::count_if(all_statistic.begin(), all_statistic.end(), std::bind2nd(std::greater_equal<double>(),all_statistic[0]));
      }


    private:
      size_t m_times;
      BaseAction *m_permute;
      unsigned m_alternative;
  };


  class VariablePermutator : public BasePermutator
  {
    // permutator for variable thresholds methods
    public:
      VariablePermutator(size_t times, char pm, unsigned alternative, const vectora & actions)
        : m_times(times), m_alternative(alternative), BasePermutator(actions)
      {
        if (pm == 'Y') {
          PermuteY* permute;
          m_permute = permute->clone();
        }
        else {
          PermuteX* permute;
          m_permute = permute->clone();
        }
      }

      ~VariablePermutator()
      {
        delete m_permute;
      }

      double apply(AssoData & d)
      {

        // obtain proper thresholds cutoffs
        vectorf maf;
        if (d.sites().size() == 0) {
          maf = d.maf();
        }
        else {
          maf.resize(0);
          vectorf mafall = d.maf();
          vectori sites = d.sites();
          for (size_t i = 0; i < mafall.size(); ++i) {
            if (sites[i] == 1) maf.push_back(mafall[i]);
          }
        }
        std::sort(maf.begin(), maf.end());
        std::vector<double>::iterator it = std::unique(maf.begin(), maf.end());
        maf.resize(it - maf.begin()); 
        double maflower = maf[0] - std::numeric_limits<double>::epsilon();

        // apply variable thresholds w/i permutation test
        unsigned permcount1 = 0;
        double obstatistic = 0.0;
        double pvalue = 9.0;

        for (size_t i = 0; i < m_times; ++i) {
          vectorf vt_statistic(0);
          for (size_t m = 0; m < maf.size(); ++m) {
            d.setSitesByMaf(maflower, maf[m]);
            for (size_t j = 0; j < m_actions.size(); ++j) {
              m_actions[j]->apply(d);
            }
            if (m_alternative != 1) {
              vt_statistic.push_back(fabs(d.statistic()));
            }
            else {
              vt_statistic.push_back(d.statistic());
            }
          }
          double statistic = *max_element(vt_statistic.begin(), vt_statistic.end()); 
          if (i == 0) {
            obstatistic = statistic;
          }
          else {
            if (statistic >= obstatistic) 
              ++permcount1;
          }
          if (pvalue <= 1.0) {
            break;
          }
          m_permute->apply(d);
        }

        if (pvalue <= 1.0) {
          d.setPvalue(pvalue); 
        }
        else {
          pvalue = (1.0 + permcount1) / (1.0 + m_times);
          d.setPvalue(pvalue);
        } 
        return 0.0;
      }


    private:
      size_t m_times;
      BaseAction *m_permute;
      unsigned m_alternative;
  };

}
#endif
