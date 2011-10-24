/*
 *  $File: utils.h $
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
#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <limits>
#include <algorithm>
#include <functional>

extern struct VPlus {
  template<typename T> std::vector<T> operator()(std::vector<T> x, std::vector<T> y) {
    std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<T>());
    return x;
  }
} vplus;

inline bool fEqual(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

inline void fRound(double& myValue, double PRECISION) 
{
  double myRemainder = fmod(myValue, PRECISION);
  if (myRemainder > ( 0.5 * PRECISION ))
  {
    myValue = (myValue - myRemainder + PRECISION);
  }
  else
  {
    myValue = (myValue  - myRemainder);
  }
  return;
}
#endif
