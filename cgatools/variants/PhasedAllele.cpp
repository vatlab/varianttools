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
#include "cgatools/variants/PhasedAllele.hpp"

namespace cgatools { namespace variants {

    void PhasedAllele::addCall(const Call& call, const Locus& locus, const reference::CrrFile& crr)
    {
        allele_.append(call.calledSequence(crr));
        pos_.push_back(allele_.size());
        calls_.push_back(&call);
        loci_.push_back(&locus);
    }

    void PhasedAllele::addSequence(const Call& call, const Locus& locus, const std::string& sequence)
    {
        allele_.append(sequence);
        pos_.push_back(allele_.size());
        calls_.push_back(&call);
        loci_.push_back(&locus);
    }

    bool PhasedAllele::hasHapLink(const std::string& hapLink) const
    {
        CGA_ASSERT(hapLink.size() != 0);
        for(size_t ii=0; ii<calls_.size(); ii++)
        {
            if (calls_[ii]->hapLink_ == hapLink)
                return true;
        }
        return false;
    }

} } // cgatools::variants
