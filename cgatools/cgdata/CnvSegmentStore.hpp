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

#ifndef CGATOOLS_CGDATA_CNVSEGMENTSTORE_HPP_
#define CGATOOLS_CGDATA_CNVSEGMENTSTORE_HPP_ 1

//! @file CnvSegmentStore.hpp

#include "cgatools/core.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"

namespace cgatools { namespace cgdata {

//! Per-segment CNV data.
//! Relative coverage will be negative if it is no-called ("N") in the file.
struct CnvSegmentData
{
    double relativeCoverage_;
    std::string calledPloidy_;
    std::string calledLevel_;
};

//! Loads CNV segment data from a file and caches it in memory.
class CnvSegmentStore :
    public reference::RangeAnnotationStore<CnvSegmentStore, CnvSegmentData>
{
public:
    //! Loads all CNV segment data from a given export package.
    CnvSegmentStore(const reference::CrrFile& crr, const GenomeMetadata& exp,
                    bool isDiploid)
        : Base(crr)
    {
        if (isDiploid)
            load( exp.getCnvSegmentsDiploidFileName() );
        else
            load( exp.getCnvSegmentsNondiploidFileName() );
    }

    //! True if the ploidy calls are present in the export package.
    bool hasCalledPloidy() const { return hasCalledPloidy_; }

    //! True if the ploidy calls are present in the export package.
    bool hasCalledLevel() const { return hasCalledLevel_; }

    //! Returns the CNV segment data for the segment with the longest overlap
    //! with the given range. In case of a tie, returns the data for the
    //! first segment in the reference order.
    const CnvSegmentData* getBestOverlappingSegment(const reference::Range& r) const;

public:
    //! This function is not a part of the true public interface of this class.
    void bindColumns(util::DelimitedFile& df, reference::Range& range, CnvSegmentData& data);

private:
    bool hasCalledPloidy_;
    bool hasCalledLevel_;
};

}}

#endif 
