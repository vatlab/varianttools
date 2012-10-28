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

#ifndef CGATOOLS_VARIANTS_ALLELEDIFFSEGMENT_HPP_
#define CGATOOLS_VARIANTS_ALLELEDIFFSEGMENT_HPP_ 1

//! @file AlleleDiffSegment.hpp

#include "cgatools/core.hpp"

namespace cgatools { namespace variants {

    //! Structure to keep track of a segment the reference, as well as
    //! calls aligned to that segment for two calldiff files. This is
    //! used for narrowing the range of non-identical PhasedAllele
    //! instances to specific calls.
    struct AlleleDiffSegment
    {
        AlleleDiffSegment(const std::pair<uint32_t,uint32_t>& lhs,
                          const std::pair<uint32_t,uint32_t>& rhs,
                          const std::pair<uint32_t,uint32_t>& ref,
                          const std::pair<uint32_t,uint32_t>& lhsCalls,
                          const std::pair<uint32_t,uint32_t>& rhsCalls)
            : lhs_(lhs),
              rhs_(rhs),
              ref_(ref),
              lhsCalls_(lhsCalls),
              rhsCalls_(rhsCalls)
        {
        }

        std::pair<uint32_t,uint32_t> lhs_; //!< pos in PhasedAllele::allele_
        std::pair<uint32_t,uint32_t> rhs_; //!< pos in PhasedAllele::allele_
        std::pair<uint32_t,uint32_t> ref_; //!< pos in reference sequence
        std::pair<uint32_t,uint32_t> lhsCalls_; //!< pos in PhasedAllele::calls_
        std::pair<uint32_t,uint32_t> rhsCalls_; //!< pos in PhasedAllele::calls_
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_ALLELEDIFFSEGMENT_HPP_
