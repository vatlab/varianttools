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

#ifndef CGATOOLS_VARIANTS_PHASEDALLELE_HPP_
#define CGATOOLS_VARIANTS_PHASEDALLELE_HPP_ 1

//! @file PhasedAllele.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/Call.hpp"
#include "cgatools/variants/Locus.hpp"

namespace cgatools { namespace variants {

    //! Class to hold allele-specific information for a PhasedHypothesis.
    class PhasedAllele
    {
    public:
        //! Creates an empty PhasedAllele.
        PhasedAllele()
            : pos_(1, 0)
        {
        }

        //! Adds the given call to the PhasedAllele.
        void addCall(const Call& call, const Locus& locus, const reference::CrrFile& crr);

        //! Adds the given call to the allele sequence, overriding the
        //! sequence for the call. This is used for no-call and ref
        //! calls whenever the locus ends in the middle of the no-call
        //! or ref call.
        void addSequence(const Call& call, const Locus& locus, const std::string& sequence);

        //! Returns true if any call in this PhasedAllele has the hapLink.
        bool hasHapLink(const std::string& hapLink) const;

        //! Returns the full sequence of this PhasedAllele.
        const std::string& allele() const
        {
            return allele_;
        }

        //! Returns the list of positions in allele_ of each call in
        //! calls_, plus allele_.size().
        const std::vector<uint32_t>& pos() const
        {
            return pos_;
        }

        //! Returns the list of calls comprising this PhasedAllele.
        const std::vector<const Call*>& calls() const
        {
            return calls_;
        }

        //! Returns the list of loci that own the calls comprising this PhasedAllele.
        const std::vector<const Locus*>& loci() const
        {
            return loci_;
        }

    private:
        std::string allele_;
        std::vector<uint32_t> pos_; //! pos in allele_ of calls_[ii]
        std::vector<const Call*> calls_;
        std::vector<const Locus*> loci_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_PHASEDALLELE_HPP_
