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
#include <limits>

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

      double check(unsigned pcount1, unsigned pcount2, size_t current, unsigned alt, double sig) const
      {
        // the adaptive p-value technique
        if (current % 1000 != 0 || current == 0) {
          return 9.0;
        }
        double x;
        if (alt == 1) {
          x = 1.0 + pcount1;
        } else
        {
          x = fmin(pcount1 + 1.0, pcount2 + 1.0);
        }


        double n = current + 1.0;
        double alpha = 0.05;
        // use 95CI for the adaptive procedure 
        // Discussions on CI see Alan AGRESTI and Brent A. COULL, 1998
        // OPTION1: Clopper-Pearson interval, conservative
        //bci <- function(n, x, alpha) {
        //  lower <- (1+(n-x+1)/(x*qf(alpha/2, 2*x, 2*(n-x+1))))^(-1)
        //  upper <- (1+(n-x+1)/((x+1)*qf(1-alpha/2, 2*(x+1), 2*(n-x))))^(-1)
        //  return(c(lower, upper))
        //}
        // FIXME the exact procedure, not usable now due to a bug in gsl_cdf_fdist_Pinv
        //double plw = 1.0 / (1.0+(n-x+1.0)/(x*gsl_cdf_fdist_Pinv(alpha/2.0, 2.0*x, 2.0*(n-x+1.0))));
        // OPTION2: Edwin B. Wilson interval, recommanded
        //wci <- function(n, x, alpha) {
        // pval <- x/n
        // z <- qnorm(1.0-alpha/2.0) 
        // zsq <- z*z
        // lower <- (pval + zsq / (2.0*n) - z * sqrt((pval*(1.0-pval)+zsq/(4.0*n))/(1.0*n))) / (1.0+zsq/(1.0*n))
        // upper <- (pval + zsq / (2.0*n) + z * sqrt((pval*(1.0-pval)+zsq/(4.0*n))/(1.0*n))) / (1.0+zsq/(1.0*n))
        // return(c(lower, upper))
        //}
        double pval = x/n;
        double z = gsl_cdf_gaussian_Pinv(1.0-alpha/2.0, 1.0);
        double zsq = z*z;
        double plw = (pval + zsq / (2.0*n) - z * sqrt((pval*(1.0-pval)+zsq/(4.0*n))/(1.0*n))) / (1.0+zsq/(1.0*n));
        plw = (alt == 1) ? plw : plw*2.0;

        if (plw > sig) { 
          return (alt == 1) ? pval : pval*2.0;
        } else
        {
          return 9.0;
        }
      }

    protected:
      vectora m_actions;

  };


  class ActionExecutor : public BasePermutator
  {

    public:
      ActionExecutor(const vectora & actions)
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
      FixedPermutator(char pm, unsigned alternative, size_t times, double sig, const vectora & actions)
        : m_times(times), m_alternative(alternative), m_sig(sig), BasePermutator(actions)
      {
        if (pm == 'Y') {
          PermuteY* permute = new PermuteY();
          m_permute = permute->clone();
        } else
        {
          PermuteX* permute = new PermuteX();
          m_permute = permute->clone();
        }
      }

      ~FixedPermutator()
      {
        delete m_permute;
      }
      
      double apply(AssoData & d)
      {
        RNG rng;
        gsl_rng* gslr = rng.get();

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
          } else
          {
            if (statistic > obstatistic) {
              ++permcount1;
            } else if (statistic < obstatistic) 
            {
              ++permcount2;
            } else
            {
              if (gsl_rng_uniform(gslr) > 0.5) ++permcount1;
              else ++permcount2;
            }
          }
          // adaptive p-value calculation
          if (m_sig < 1.0) {
            pvalue = check(permcount1, permcount2, i, m_alternative, m_sig);
          }
          if (pvalue <= 1.0) {
            break;
          }
          m_permute->apply(d);
        }

        if (pvalue <= 1.0) {
          d.setPvalue(pvalue); 
        } else
        {
          if (m_alternative == 1) { 
            pvalue = (permcount1 + 1.0) / (m_times + 1.0);
          } else 
          {
            double permcount = fmin(permcount1, permcount2);
            pvalue = 2.0 * (permcount + 1.0) / (m_times + 1.0);
          }
          d.setPvalue(pvalue);
        } 

        d.setStatistic(obstatistic);
        return 0.0;
        //return (double) std::count_if(all_statistic.begin(), all_statistic.end(), std::bind2nd(std::greater_equal<double>(),all_statistic[0]));
      }


    private:
      size_t m_times;
      BaseAction *m_permute;
      unsigned m_alternative;
      double m_sig;
  };


  class VariablePermutator : public BasePermutator
  {
    // permutator for variable thresholds methods
    public:
      VariablePermutator(char pm, unsigned alternative, size_t times, double sig, const vectora & actions)
        : m_times(times), m_alternative(alternative), m_sig(sig), BasePermutator(actions)
      {
        if (pm == 'Y') {
          PermuteY* permute = new PermuteY();
          m_permute = permute->clone();
        } else 
        {
          PermuteX* permute = new PermuteX();
          m_permute = permute->clone();
        }
      }

      ~VariablePermutator()
      {
        delete m_permute;
      }

      double apply(AssoData & d)
      {
        if (d.maf().size()==0) {
          throw RuntimeError("MAF has not been calculated. Please calculate MAF prior to using variable thresholds method.");
        }
        
        RNG rng;
        gsl_rng* gslr = rng.get();

        // obtain proper thresholds cutoffs
        vectorf maf;
        if (d.sites().size() == 0) {
          maf = d.maf();
        } else
        {
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
        if (fEqual(maf.front(), 0.0)) {
          maf.erase(maf.begin());
        }
        if (fEqual(maf.back(), 1.0)) {
          maf.erase(maf.end());
        }
        if (maf.size() == 0) {
          // nothing to do
          d.setPvalue(std::numeric_limits<double>::quiet_NaN());
          d.setStatistic(std::numeric_limits<double>::quiet_NaN());
          return 0.0;
        }

        double maflower = maf.front() - std::numeric_limits<double>::epsilon();

        // apply variable thresholds w/i permutation test
        unsigned permcount1 = 0, permcount2 = 0;
        double max_obstatistic = 0.0, min_obstatistic = 0.0;
        double pvalue = 9.0;

        for (size_t i = 0; i < m_times; ++i) {
          vectorf vt_statistic(0);
          for (size_t m = 0; m < maf.size(); ++m) {
            d.setSitesByMaf(maf[m], maflower);
            for (size_t j = 0; j < m_actions.size(); ++j) {
              m_actions[j]->apply(d);
            }
            vt_statistic.push_back(d.statistic());
          }
          double max_statistic = *max_element(vt_statistic.begin(), vt_statistic.end()); 
          double min_statistic = *min_element(vt_statistic.begin(), vt_statistic.end()); 
          if (i == 0) {
            max_obstatistic = max_statistic;
            min_obstatistic = min_statistic;
          } else          
          {
            if (max_statistic>= max_obstatistic && min_statistic <= min_obstatistic) {
              if (gsl_rng_uniform(gslr) > 0.5) ++permcount1;
              else ++permcount2;
            } else
            {
              if (max_statistic >= max_obstatistic) {
                ++permcount1;
              }
              if (min_statistic <= min_obstatistic) {
                ++permcount2;
              }
            }
          }

          // adaptive p-value calculation
          if (m_sig < 1.0) {
            pvalue = check(permcount1, permcount2, i, m_alternative, m_sig);
          }
          if (pvalue <= 1.0) {
            break;
          }
          m_permute->apply(d);
        }

        if (pvalue <= 1.0) {
          d.setPvalue(pvalue); 
        } else
        {
          if (m_alternative == 1) { 
            pvalue = (permcount1 + 1.0) / (m_times + 1.0);
          } else
          {
            double permcount = fmin(permcount1, permcount2);
            pvalue = 2.0 * (permcount + 1.0) / (m_times + 1.0);
          }
          d.setPvalue(pvalue);
        } 

        // set statistic, a bit involved
        if (m_alternative == 1) { 
          d.setStatistic(max_obstatistic);
        } else 
        {
          (permcount1 >= permcount2)?d.setStatistic(max_obstatistic):d.setStatistic(min_obstatistic);
        }

        return 0.0;
      }

    private:
      size_t m_times;
      BaseAction *m_permute;
      unsigned m_alternative;
      double m_sig;
  };

}
#endif
