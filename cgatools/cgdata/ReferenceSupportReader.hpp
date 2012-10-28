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

#ifndef CGATOOLS_CGDATA_REFERENCESUPPORTREADER_HPP_
#define CGATOOLS_CGDATA_REFERENCESUPPORTREADER_HPP_ 1

//! @file ReferenceSupportReader.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"

#include <limits>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace cgdata {

    //! Sequential reader for CGI reference support files provided
    //! as a part of the genome data package.
    class ReferenceSupportReader : boost::noncopyable
    {
    public:
        ReferenceSupportReader(const reference::CrrFile& crr,
                               const GenomeMetadata& exp);

        //! Buffer reference support data for the specified range of the genome.
        //! Switching between chromosomes is expensive. Within chromosome, only forward
        //! movement is allowed: the file is always read sequentially. In other words,
        //! if the currently buffered range is OR, and the new range is NR, then
        //! (OR.begin_ <= NR.begin_ && OR.end_ <= NR.end_) must hold.
        void seek(const reference::Range& r);

        //! Returns the minimum reference score for the specified range.
        //! The range must be within the currently buffered segment. If no data
        //! is in the file for the specified range, returns noDataValue.
        int32_t getMinScore(const reference::Range& r,
                            int32_t noDataValue =
                                std::numeric_limits<int32_t>::min() ) const;

        //! Returns the maximum sum-of-weights coverage number over all bases of
        //! the specified range.
        int32_t getMaxWeightedCoverage(const reference::Range& r,
                                       int32_t noDataValue = 0) const;

        //! Returns the maximum unique sequence coverage (fully mapped)
        //! over all bases of the specified range.
        int32_t getMaxUniqueSequenceCoverage(const reference::Range& r,
                                             int32_t noDataValue = 0) const;

    private:
        struct File
        {
            std::string filename_;
            boost::shared_ptr<std::istream> f_;
            util::DelimitedFile parser_;
            bool hasWeightedCoverage_;
            uint32_t offset_;
            int32_t score_;
            int32_t uniqueCoverage_;
            int32_t weightedCoverage_;

            File(const std::string& fn);
        };

        struct DataItem
        {
            bool hasData_;
            int32_t score_;
            int32_t uniqueCoverage_;
            int32_t weightedCoverage_;

            DataItem() : hasData_(false), score_(0) {}
        };

        const reference::CrrFile& crr_;
        const GenomeMetadata& exp_;

        uint16_t chromosome_;
        uint32_t bufferStartOffset_;
        std::vector<DataItem> buffer_;
        boost::scoped_ptr<File> file_;

        void checkRange(const reference::Range& r) const;
    };

} } // cgatools::cgdata

#endif // CGATOOLS_CGDATA_REFERENCESUPPORTREADER_HPP_
