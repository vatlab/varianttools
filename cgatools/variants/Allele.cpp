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
#include "cgatools/variants/Allele.hpp"
#include "cgatools/variants/Locus.hpp"
#include "cgatools/variants/Call.hpp"

#include <boost/algorithm/string.hpp>

#include <set>

using namespace std;

namespace cgatools { namespace variants {

    void Allele::clearCalls()
    {
        callOffsets_.clear();
    }

    void Allele::addCallOffset(size_t offset)
    {
        callOffsets_.push_back(offset);
    }

    void Allele::setLocus(const Locus* locus)
    {
        locus_ = locus;
    }

    std::string Allele::getAlleleSequence() const
    {
        const vector<Call>& calls = locus_->getCalls();
        string result;
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
            result += calls[callOffsets_[ii]].calledSequence(locus_->getReference());
        return result;
    }

    int32_t Allele::getMinVarScoreVAF() const
    {
        const vector<Call>& calls = locus_->getCalls();
        int32_t s = std::numeric_limits<int32_t>::max();
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
            s = std::min(s, calls[callOffsets_[ii]].varScoreVAF_);
        return s;
    }

    int32_t Allele::getMinVarScoreEAF() const
    {
        const vector<Call>& calls = locus_->getCalls();
        int32_t s = std::numeric_limits<int32_t>::max();
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
            s = std::min(s, calls[callOffsets_[ii]].varScoreEAF_);
        return s;
    }

    Call::VarQuality Allele::getMinVarQuality() const
    {
        const vector<Call>& calls = locus_->getCalls();
        int s = Call::VAR_QUALITY_HIGH;
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
            s = std::min(s, int(calls[callOffsets_[ii]].varQuality_));
        return Call::VarQuality(s);
    }

    const Call& Allele::getCall(size_t index) const
    {
        CGA_ASSERT(index < callOffsets_.size());
        return locus_->getCalls()[callOffsets_[index]];
    }

    std::string Allele::getRefSequence() const
    {
        const vector<Call>& calls = locus_->getCalls();
        string result;
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
            result += calls[callOffsets_[ii]].refSequence(locus_->getReference());
        return result;
    }

    const std::string& Allele::getHapLink() const
    {
        const vector<Call>& calls = locus_->getCalls();
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
        {
            if (calls[callOffsets_[ii]].hapLink_.size() > 0)
                return calls[callOffsets_[ii]].hapLink_;
        }
        return calls[callOffsets_[0]].hapLink_;
    }

    std::string Allele::getXRef() const
    {
        set<string> xrefs;
        const vector<Call>& calls = locus_->getCalls();
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
        {
            const Call& c = calls[callOffsets_[ii]];
            if (!c.xRef_.empty())
            {
                string::const_iterator first = c.xRef_.begin();
                size_t pos = 0;
                while (true)
                {
                    size_t next = c.xRef_.find_first_of(';', pos);
                    if (string::npos == next)
                    {
                        xrefs.insert(string(first+pos, c.xRef_.end()));
                        break;
                    }
                    xrefs.insert(string(first+pos, first+next));
                    pos = next+1;
                }
            }
        }
        return boost::join(xrefs, ";");
    }

    bool Allele::hasNoCalls() const
    {
        const vector<Call>& calls = locus_->getCalls();
        for(size_t ii=0; ii<callOffsets_.size(); ii++)
        {
            if (calls[callOffsets_[ii]].hasNoCalls())
                return true;
        }
        return false;
    }

} } // cgatools::variants
