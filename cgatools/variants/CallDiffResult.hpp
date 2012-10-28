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

#ifndef CGATOOLS_VARIANTS_CALLDIFFRESULT_HPP_
#define CGATOOLS_VARIANTS_CALLDIFFRESULT_HPP_ 1

//! @file CallDiffResult.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/Call.hpp"
#include "cgatools/variants/CallDiffMatchType.hpp"
#include "cgatools/variants/PhasedHypothesis.hpp"

namespace cgatools { namespace variants {

    //! Class to hold the result of diffing two superloci.
    class CallDiffResult
    {
    public:
        void setHypotheses(const PhasedHypothesis& lhs, const PhasedHypothesis& rhs,
                           const std::vector<size_t>& perm);

        //! For calls participating in two CallDiffSegment instances or
        //! alleles having multiple calls, this returns the merged match
        //! type associated with the call or allele as a whole.
        static cdmt::MatchType mergeMatchTypes(cdmt::MatchType mtA, cdmt::MatchType mtB);

        //! Returns string identifier associated with the given match type.
        static const char* getMatchTypeString(cdmt::MatchType mt);

        //! Returns ;-separated string identifier for the given list of
        //! cdmt::MatchType instances.
        static std::string getMatchTypeString(const std::vector<cdmt::MatchType>& mt);

        //! Returns string identifier that is the simplified version of
        //! the match type identifier (i.e. ref-identical -> identical,
        //! alt-identical -> identical).
        static std::string getSimpleMatchTypeString(const std::vector<cdmt::MatchType>& mt);

        //! Returns true for the identical and compatible match types.
        static bool isConsistent(cdmt::MatchType mt);

        //! Returns true for the identical match types.
        static bool isIdentical(cdmt::MatchType mt);

        //! For each allele, the match type.
        std::vector<cdmt::MatchType> matchType_;

        //! For each allele, the match type if the broken-down
        //! AlleleDiffSegments are used.
        std::vector<cdmt::MatchType> matchTypeBySegment_;

        //! For each file, the PhasedHypothesis for this CallDiffResult.
        boost::array<PhasedHypothesis, 2> hyp_;

        //! For each file, for each allele, for each call, the cdmt::MatchType
        //! and Call.
        boost::array<std::vector< std::vector< std::pair<cdmt::MatchType, const Call*> > >, 2> callClass_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_CALLDIFFRESULT_HPP_
