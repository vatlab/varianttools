/**
 *  $File: action.h $
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
        throw RuntimeError("The base action class should not be called");
        return 0;
      }
  };


  typedef std::vector<BaseAction *> vectora;


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
        d.meanOfX();
        return 0;
      }
  };


  class BinToX : public BaseAction
  {
    public:
      BinToX() : BaseAction()
    {
    }

      BaseAction * clone()
      {
        return new BinToX(*this);
      }

      double apply(AssoData & d)
      {
        d.binToX();
        d.meanOfX();
        return 0;
      }
  };


  class SetMaf : public BaseAction
  {
    public:
      SetMaf() : BaseAction()
    {
    }

      BaseAction * clone()
      {
        return new SetMaf(*this);
      }

      double apply(AssoData & d)
      {
        d.setMaf();
        return 0;
      }
  };


  class SetSites : public BaseAction
  {
    public:
      SetSites(double upper=1.0, double lower=0.0) : 
        BaseAction(), m_upper(upper), m_lower(lower)
    {
    }

      BaseAction * clone()
      {
        return new SetSites(*this);
      }

      double apply(AssoData & d)
      {
        d.setSitesByMaf(m_upper, m_lower);
        return 0;
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

      BaseAction * clone()
      {
        return new SimpleLinearRegression(*this);
      }

      double apply(AssoData & d)
      {    
        d.simpleLinear();
        return 0;
      }
  };


  class MultipleLinearRegression : public BaseAction
  {
    public:
      MultipleLinearRegression() : BaseAction()
    {
    }

      BaseAction * clone()
      {
        return new MultipleLinearRegression(*this);
      }

      double apply(AssoData & d)
      {    
        d.multipleLinear();
        return 0;
      }
  };


  class SimpleLogisticRegression : public BaseAction
  {
    public:
      SimpleLogisticRegression() : BaseAction()
    {
    }

      BaseAction * clone()
      {
        return new SimpleLogisticRegression(*this);
      }

      double apply(AssoData & d)
      {    
        //assert(m_phenotype.size() == m_X.size());
        //check for binary
        d.simpleLogit();
        return 0;
      }
  };


  class GaussianPval : public BaseAction
  {
    public:
      GaussianPval(unsigned sided=1) 
        : BaseAction(), m_sided(sided)
      {
      }

      BaseAction * clone()
      {
        return new GaussianPval(*this);
      }

      double apply(AssoData & d)
      {
        d.gaussianP(m_sided); 
        return 0;
      }

    private:
      unsigned m_sided;
  };


  class StudentPval : public BaseAction
  {
    public:
      StudentPval(unsigned sided=1) 
        : BaseAction(), m_sided(sided)
      {
      }

      BaseAction * clone()
      {
        return new StudentPval(*this);
      }

      double apply(AssoData & d)
      {
        d.studentP(m_sided); 
        return 0;
      }

    private:
      unsigned m_sided;
  };

}

#endif
