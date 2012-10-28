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

#ifndef CGATOOLS_VARIANTS_SUPERLOCUS_HPP_
#define CGATOOLS_VARIANTS_SUPERLOCUS_HPP_ 1

//! @file Superlocus.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/PhasedHypothesis.hpp"
#include "cgatools/variants/Locus.hpp"

#include <deque>

namespace cgatools { namespace variants {

    class SuperlocusIterator;

    //! The Superlocus class. Describes a Locus list for a set of files
    //! that must be considered together in context.
    class Superlocus
    {
    public:
        //! Constructs an empty Superlocus.
        Superlocus();

        //! Returns all the hapLink-consistent phasings of this Superlocus.
        void buildPhasedHypotheses(std::vector< std::vector<PhasedHypothesis> >& hypotheses,
                                   size_t maxHypothesisCount,
                                   bool useHapLinks) const;

        //! Returns the reference range for this Superlocus.
        const reference::Range& getRange() const
        {
            return range_;
        }

        //! Returns a numeric identifier for the Superlocus.
        uint32_t getId() const
        {
            return id_;
        }

        //! Returns the set of Locus for one file of this
        //! Superlocus.
        std::pair<std::deque<Locus>::const_iterator, std::deque<Locus>::const_iterator>
        getLoci(size_t fileOffset) const;

        size_t getGenomeCount() const
        {
            return queues_.size();
        }

        static uint32_t extendRightByPrefixMatching(
            const reference::CrrFile& crr, const reference::Location& loc, const std::string& sequence);
        static uint32_t extendLeftBySuffixMatching(
            const reference::CrrFile& crr, const reference::Location& loc, const std::string& sequence);

    private:
        void buildPhasedHypotheses(std::vector<PhasedHypothesis>& hypotheses,
                                   size_t maxHypothesisCount,
                                   bool useHapLinks,
                                   size_t fileOffset) const;
        bool hypothesisPermutationsAreEqual(const PhasedHypothesis& hypothesis) const;
        bool callPermutationsAreEqual(const Call& lhs, const Call& rhs) const;
        bool allelePermutationsAreEqual(const Locus& locus) const;
        void addCalls(PhasedHypothesis& hypothesis,
                      const Locus& locus,
                      const std::vector<size_t>& perm) const;
        void addSequence(PhasedHypothesis& hypothesis,
                         const Locus& locus,
                         const std::vector<size_t>& perm,
                         const std::string& sequence) const;
        bool areHapLinksConsistent(const PhasedHypothesis& hypothesis,
                                   const Locus& locus,
                                   const std::vector<size_t>& perm) const;
        void downSample(std::vector<PhasedHypothesis>& hypotheses,
                        size_t maxHypothesisCount) const;
        uint32_t hashHypothesis(const PhasedHypothesis& hyp) const;

        std::vector< std::deque<Locus> > queues_;
        reference::Range range_;
        reference::Range searchRange_;
        uint32_t id_;

        friend class SuperlocusIterator;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_SUPERLOCUS_HPP_
