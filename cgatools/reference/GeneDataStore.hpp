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

#ifndef CGATOOLS_REFERENCE_GENEDATASTORE_HPP_
#define CGATOOLS_REFERENCE_GENEDATASTORE_HPP_ 1

//! @file GeneDataStore.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/RangeIntersector.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/reference/CrrFile.hpp"

namespace cgatools { namespace reference {

struct GeneDescriptionData
{
    std::string mrnaAcc_, proteinAcc_, geneId_, symbol_;
    bool strand_;
    std::vector<Range> exonRanges_;
    Range genomeCodingRange_;
};

class GeneDataStore
{
private:
    typedef util::IntervalTree<Range,
                               Location,
                               GeneDescriptionData,
                               RangeOverlap,
                               GetRangeBoundary > DataStore;
public:
    typedef DataStore::QueryResultType QueryResultType;

    GeneDataStore(const reference::CrrFile& crr, const std::string& fn)
    {
        load(crr, fn);
    }

    const util::DelimitedFileMetadata& getMetadata() const
    {
        return metadata_;
    }

    void intersect(const reference::Range& range,
                   std::vector<QueryResultType>& result) const
    {
        tree_.intersect(range, result);
    }

private:
    util::DelimitedFileMetadata metadata_;
    DataStore tree_;

    void load(const reference::CrrFile& crr, const std::string& fn);
};

}}

#endif 
