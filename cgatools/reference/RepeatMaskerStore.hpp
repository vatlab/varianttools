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

#ifndef CGATOOLS_REFERENCE_REPEATMASKERSTORE_HPP_
#define CGATOOLS_REFERENCE_REPEATMASKERSTORE_HPP_ 1

//! @file RepeatMaskerStore.hpp

#include "cgatools/core.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"

namespace cgatools { namespace reference {

struct RepeatMaskerAnnotation
{
    std::string name_, family_;
    double divergence_;
    bool strand_;
};

class RepeatMaskerStore :
    public RangeAnnotationStore<RepeatMaskerStore, RepeatMaskerAnnotation>
{
public:

    RepeatMaskerStore(const reference::CrrFile& crr, const std::string& fn)
        : Base(crr)
    {
        load(fn);
    }

    void bindColumns(util::DelimitedFile& df,
                     reference::Range& range,
                     RepeatMaskerAnnotation& data)
    {
        using namespace util;
        bindRangeColumns(df, range);
        df.addField(StringField("repName", &data.name_));
        df.addField(StringField("repFamily", &data.family_));
        df.addField(ValueField<double>("divergence", &data.divergence_));
        df.addField(StrandField("strand", &data.strand_));
    }
};

}}

#endif 
