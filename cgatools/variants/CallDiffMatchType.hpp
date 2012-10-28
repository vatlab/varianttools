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

#ifndef CGATOOLS_VARIANTS_CALLDIFFMATCHTYPE_HPP_
#define CGATOOLS_VARIANTS_CALLDIFFMATCHTYPE_HPP_ 1

//! @file CallDiffMatchType.hpp

#include "cgatools/core.hpp"

namespace cgatools { namespace variants { namespace cdmt {

    //! Enumeration describing the type of match found one one
    //! allele.
    enum MatchType
    {
        REF_IDENTICAL   = 0, //!< files identical, ref-consistent
        ALT_IDENTICAL   = 1, //!< files identical, ref-inconsistent
        REF_CONSISTENT  = 2, //!< files consistent, ref-consistent
        ALT_CONSISTENT  = 3, //!< files consistent, ref-inconsistent
        ONLY_A          = 4, //!< files inconsistent, A is ref-inconsistent
        ONLY_B          = 5, //!< files inconsistent, B is ref-inconsistent
        MISMATCH        = 6, //!< files inconsistent, both ref-inconsistent
        PHASE_MISMATCH  = 7, //!< files inconsistent, consistent if no hapLink
        PLOIDY_MISMATCH = 8, //!< superlocus inconsistent ploidy
        MATCH_TYPE_LAST = 9
    };

} } } // cgatools::variants::cdmt

#endif // CGATOOLS_VARIANTS_CALLDIFFMATCHTYPE_HPP_
