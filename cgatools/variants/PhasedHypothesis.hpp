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

#ifndef CGATOOLS_VARIANTS_PHASEDHYPOTHESIS_HPP_
#define CGATOOLS_VARIANTS_PHASEDHYPOTHESIS_HPP_ 1

//! @file PhasedHypothesis.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/PhasedAllele.hpp"
#include "cgatools/variants/CallDiffMatchType.hpp"
#include "cgatools/variants/AlleleDiffSegment.hpp"

namespace cgatools { namespace variants {

    class Superlocus;
    class CallDiffResult;

    //! Class to hold a single possible phasing of loci in a Superlocus.
    class PhasedHypothesis
    {
    public:
        //! Constructs an empty PhasedHypothesis.
        PhasedHypothesis();

        //! Constructs the PhasedHypothesis.
        PhasedHypothesis(const reference::Range& range, size_t ploidy);

        //! Returns the ploidy of the PhasedHypothesis.
        size_t size() const
        {
            return alleles_.size();
        }

        //! Returns the reference range of the PhasedHypothesis.
        const reference::Range& getRange() const
        {
            return range_;
        }

        //! Returns the given PhasedAllele.
        const PhasedAllele& operator[](size_t offset) const
        {
            return alleles_[offset];
        }

        //! Returns the given PhasedAllele.
        PhasedAllele& operator[](size_t offset)
        {
            return alleles_[offset];
        }

        void swap(PhasedHypothesis& other);

        //! Returns the best calldiff result, given the list of
        //! PhasedHypothesis for two files.
        static void findBestDiff(const Superlocus& sl,
                                 size_t maxHypothesisCount,
                                 const std::vector<PhasedHypothesis>& lhs,
                                 const std::vector<PhasedHypothesis>& rhs,
                                 const reference::CrrFile& crr,
                                 CallDiffResult& result);

        static void testVariant(const Superlocus& sl,
                                const std::vector<PhasedHypothesis>& hyps,
                                const reference::Range& varRange,
                                const std::string& varAllele,
                                const reference::CrrFile& crr,
                                std::string& result);

    private:
        static void diffHypotheses(const PhasedHypothesis& lhs,
                                   const PhasedHypothesis& rhs,
                                   const std::string& refSequence,
                                   const std::vector<size_t>& perm,
                                   std::vector<cdmt::MatchType>& result);
        static void testVariantOneHypothesis(
            const PhasedHypothesis& hyp,
            const PhasedAllele& rhs,
            const reference::Range& variantRange,
            const reference::Range& extendedVariantRange,
            const std::string& refSequence,
            std::string& result);
        static char testPhasedAllele(const PhasedAllele& lhs,
                                     const PhasedAllele& rhs,
                                     const std::vector<AlleleDiffSegment>& segs,
                                     int lSeg,
                                     int rSeg);

        static void classifyCalls(
            const PhasedAllele& lhs,
            const PhasedAllele& rhs,
            const std::string& refSequence,
            const reference::CrrFile& crr,
            const reference::Range& range,
            cdmt::MatchType mt,
            std::vector< std::pair<cdmt::MatchType, const Call*> >& lhsCallClass,
            std::vector< std::pair<cdmt::MatchType, const Call*> >& rhsCallClass,
            cdmt::MatchType& mtBySeg);

        static cdmt::MatchType matchTypeForCall(
            cdmt::MatchType mt,
            const reference::CrrFile& crr,
            const Call* call);

        static void splitIntoSegments(const PhasedAllele& lhs,
                                      const PhasedAllele& rhs,
                                      const std::string& refSequence,
                                      const reference::Range& range,
                                      std::vector<AlleleDiffSegment>& segs);
        static bool segCanSplit(const Call& call, uint32_t pos);

        reference::Range range_;
        std::vector<PhasedAllele> alleles_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_PHASEDHYPOTHESIS_HPP_
