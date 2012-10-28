// Copyright 2011 Complete Genomics, Inc.
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

#ifndef CGATOOLS_CGDATA_CNVDETAILSTORE_HPP_
#define CGATOOLS_CGDATA_CNVDETAILSTORE_HPP_ 1

//! @file CnvDetailStore.hpp

#include "cgatools/core.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"

namespace cgatools { namespace cgdata {

//! Per-window CNV data.
//! Relative coverage will be negative if it is no-called ("N") in the file.
struct CnvDetailData
{
    double relativeCoverage_;
    std::string calledPloidy_;
    std::string calledLevel_;
    std::string bestLAF_;
    std::string lowLAF_;
    std::string highLAF_;
};

//! Loads CNV window data from a file and caches it in memory.
class CnvDetailStore :
    public reference::RangeAnnotationStore<CnvDetailStore, CnvDetailData>
{
public:
    //! Loads all CNV window data from a given export package.
    CnvDetailStore(const reference::CrrFile& crr, const GenomeMetadata& exp,
                   bool isDiploid, bool isSomatic)
        : Base(crr)
    {
        if (isDiploid){
            if (isSomatic){
                load( exp.getCnvDetailsSomaticDiploidFileName() );
            } else {
                load( exp.getCnvDetailsDiploidFileName() );
            }
        } else {
            if (isSomatic){
                load( exp.getCnvDetailsSomaticNondiploidFileName() );
            } else {
                load( exp.getCnvDetailsNondiploidFileName() );
            }
        }
    }

    //! True if the ploidy calls are present in the export package.
    bool hasCalledPloidy() const { return hasCalledPloidy_; }

    //! True if the level calls are present in the export package.
    bool hasCalledLevel() const { return hasCalledLevel_; }

    //! True if LAF information is present in the Details file.
    bool hasLAF() const { return hasLAF_; }

    //! Returns the CNV window data for the window with the longest overlap
    //! with the given range. In case of a tie, returns the data for the
    //! first window in the reference order.
    const CnvDetailData* getBestOverlappingDetail(const reference::Range& r) const;

public:
    //! This function is not a part of the true public interface of this class.
    void bindColumns(util::DelimitedFile& df, reference::Range& range, CnvDetailData& data);

private:
    bool hasCalledPloidy_;
    bool hasCalledLevel_;
    bool hasLAF_;
};

}}

#endif 
