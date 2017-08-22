// Copyright 2010 Complete Genomics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License. You
// may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

#ifndef CGATOOLS_UTIL_INDIRECTCOMPARATOR_HPP_
#define CGATOOLS_UTIL_INDIRECTCOMPARATOR_HPP_ 1

//! @file IndirectComparator.hpp

#include "cgatools/core.hpp"

namespace cgatools { namespace util {

    template <class Container>
    class IndirectComparator
    {
    public:
        IndirectComparator(const Container& cc)
            : cc_(cc)
        {
        }

        template <class Index>
        bool operator()(const Index& lhs, const Index& rhs) const
        {
            return cc_[lhs] < cc_[rhs];
        }

    private:
        const Container& cc_;
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_INDIRECTCOMPARATOR_HPP_
