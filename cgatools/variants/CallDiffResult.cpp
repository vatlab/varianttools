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

#include "cgatools/core.hpp"
#include "cgatools/variants/CallDiffResult.hpp"
#include "cgatools/util/IndirectComparator.hpp"

namespace cgatools { namespace variants {

    using namespace std;
    using namespace cgatools::util;
    using cdmt::MatchType;

    void CallDiffResult::setHypotheses(const PhasedHypothesis& lhs,
                                       const PhasedHypothesis& rhs,
                                       const std::vector<size_t>& perm)
    {
        CGA_ASSERT(lhs.size() == rhs.size());
        CGA_ASSERT(lhs.size() == perm.size());
        CGA_ASSERT(lhs.size() == matchType_.size());

        vector<size_t> matchTypeOrder;
        for(size_t ii=0; ii<matchType_.size(); ii++)
            matchTypeOrder.push_back(ii);
        std::sort(matchTypeOrder.begin(), matchTypeOrder.end(),
                  IndirectComparator< vector<MatchType> >(matchType_));

        hyp_[0] = PhasedHypothesis(lhs.getRange(), lhs.size());
        hyp_[1] = PhasedHypothesis(rhs.getRange(), rhs.size());
        for(size_t ii=0; ii<matchType_.size(); ii++)
        {
            hyp_[0][ii] = lhs[perm[matchTypeOrder[ii]]];
            hyp_[1][ii] = rhs[matchTypeOrder[ii]];
        }

        std::sort(matchType_.begin(), matchType_.end());
    }

    cdmt::MatchType CallDiffResult::mergeMatchTypes(
        cdmt::MatchType mtA, cdmt::MatchType mtB)
    {
        if (mtB < mtA)
            std::swap(mtA, mtB);
        if (cdmt::ALT_IDENTICAL == mtA && cdmt::REF_CONSISTENT == mtB)
            return cdmt::ALT_CONSISTENT;
        if ( (cdmt::ALT_IDENTICAL == mtA || cdmt::ALT_CONSISTENT == mtA) &&
             cdmt::ONLY_A <= mtB )
            return cdmt::MISMATCH;
        if (cdmt::ONLY_A == mtA && cdmt::ONLY_B == mtB)
            return cdmt::MISMATCH;
        return mtB;
    }

    const char* CallDiffResult::getMatchTypeString(cdmt::MatchType mt)
    {
        switch(mt)
        {
        case cdmt::REF_IDENTICAL:
            return "ref-identical";
        case cdmt::REF_CONSISTENT:
            return "ref-consistent";
        case cdmt::ALT_IDENTICAL:
            return "alt-identical";
        case cdmt::ALT_CONSISTENT:
            return "alt-consistent";
        case cdmt::ONLY_A:
            return "onlyA";
        case cdmt::ONLY_B:
            return "onlyB";
        case cdmt::MISMATCH:
            return "mismatch";
        case cdmt::PHASE_MISMATCH:
            return "phase-mismatch";
        case cdmt::PLOIDY_MISMATCH:
            return "ploidy-mismatch";
        default:
            CGA_ASSERT(false);
            return "unknown";
        }
    }

    std::string CallDiffResult::getMatchTypeString(const std::vector<cdmt::MatchType>& mt)
    {
        vector<MatchType> mtSorted(mt);
        std::sort(mtSorted.begin(), mtSorted.end());
        string result;
        for(size_t ii=0; ii<mtSorted.size(); ii++)
        {
            if (ii>0)
                result += ";";
            result += getMatchTypeString(mtSorted[ii]);
        }
        return result;
    }

    std::string CallDiffResult::getSimpleMatchTypeString(const std::vector<cdmt::MatchType>& mt)
    {
        CGA_ASSERT(mt.size() > 0);
        MatchType mm = mt[0];
        for(size_t ii=1; ii<mt.size(); ii++)
            mm = mergeMatchTypes(mt[ii], mm);
        switch(mm)
        {
        case cdmt::REF_IDENTICAL:
            return "identical";
        case cdmt::ALT_IDENTICAL:
            return "identical";
        case cdmt::REF_CONSISTENT:
            return "consistent";
        case cdmt::ALT_CONSISTENT:
            return "consistent";
        case cdmt::ONLY_A:
            return "onlyA";
        case cdmt::ONLY_B:
            return "onlyB";
        case cdmt::MISMATCH:
            return "mismatch";
        case cdmt::PHASE_MISMATCH:
            return "phase-mismatch";
        case cdmt::PLOIDY_MISMATCH:
            return "ploidy-mismatch";
        default:
            CGA_ASSERT(false);
            return "unknown";
        }
    }

    bool CallDiffResult::isConsistent(cdmt::MatchType mt)
    {
        return mt <= cdmt::ALT_CONSISTENT;
    }

    bool CallDiffResult::isIdentical(cdmt::MatchType mt)
    {
        return mt <= cdmt::ALT_IDENTICAL;
    }

} } // cgatools::variants
