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
#include "cgatools/util/Streams.hpp"
#include "cgatools/cgdata/ReferenceSupportReader.hpp"

#include <boost/format.hpp>

namespace cgatools { namespace cgdata {

    ReferenceSupportReader::ReferenceSupportReader(
        const reference::CrrFile& crr, const GenomeMetadata& exp)
        :   crr_(crr), exp_(exp), chromosome_(65535)
    {
    }

    ReferenceSupportReader::File::File(const std::string& fn)
        :   filename_(fn),
            f_(util::InputStream::openCompressedInputStreamByExtension(fn)),
            parser_(*f_, fn)
    {
        int version = parser_.getMetadata().getFormatVersion();

        parser_.addField(util::ValueField<uint32_t>("offset", &offset_));
        parser_.addField(util::ValueField<int32_t>("refScore", &score_));

        if (version >= 1001)
        {
            hasWeightedCoverage_ = true;
            parser_.addField(util::ValueField<int32_t>("uniqueSequenceCoverage", &uniqueCoverage_));
            parser_.addField(util::ValueField<int32_t>("weightSumSequenceCoverage", &weightedCoverage_));
        }
        else
        {
            hasWeightedCoverage_ = false;
            parser_.addField(util::ValueField<int32_t>("coverage", &uniqueCoverage_));
        }
    }

    void ReferenceSupportReader::seek(const reference::Range& r)
    {
        if (r.chromosome_ != chromosome_)
        {
            std::string fn =
                exp_.getReferenceSupportFileName(
                    crr_.listChromosomes()[r.chromosome_].getName());
            file_.reset(new File(fn));
            chromosome_ = r.chromosome_;
            if (!file_->parser_.next())
                throw util::Exception("no data in file: " + file_->filename_);
            buffer_.clear();
            bufferStartOffset_ = 0;
        }

        if (r.begin_ < bufferStartOffset_)
            throw util::Exception("reference support reader: cannot seek backward");
        if (bufferStartOffset_ + buffer_.size() > r.begin_)
            // Part of the buffer is still useful, move it to the beginning
            buffer_.erase(buffer_.begin(),
                          buffer_.begin() + r.begin_ - bufferStartOffset_);
        else
            buffer_.clear();
        bufferStartOffset_ = r.begin_;
        if (r.length() < buffer_.size())
            throw util::Exception("reference support reader: cannot seek backward");
        buffer_.resize(r.length(), DataItem());

        while (file_->offset_ < r.end_)
        {
            if (file_->offset_ >= bufferStartOffset_)
            {
                DataItem& d = buffer_[file_->offset_ - bufferStartOffset_];
                d.hasData_ = true;
                d.score_ = file_->score_;
                d.uniqueCoverage_ = file_->uniqueCoverage_;
                d.weightedCoverage_ = file_->weightedCoverage_;
            }
            if (!file_->parser_.next())
            {
                file_->offset_ = std::numeric_limits<uint32_t>::max();
                break;
            }
        }
    }

    void ReferenceSupportReader::checkRange(const reference::Range& r) const
    {
        if (r.chromosome_ != chromosome_ || r.begin_ < bufferStartOffset_ ||
            r.end_ > bufferStartOffset_ + buffer_.size())
        {
            throw util::Exception("specified data range is not in the buffer");
        }
    }

    int32_t ReferenceSupportReader::getMinScore(const reference::Range& r,
                                                int32_t noDataValue) const
    {
        checkRange(r);

        int32_t val = std::numeric_limits<int32_t>::max();
        bool hasData = false;
        for (size_t ii = r.begin_; ii < r.end_; ++ii)
        {
            const DataItem& d = buffer_[ii - bufferStartOffset_];
            if (d.hasData_)
            {
                val = std::min(val, d.score_);
                hasData = true;
            }
        }
        if (hasData)
            return val;
        else
            return noDataValue;
    }

    int32_t ReferenceSupportReader::getMaxWeightedCoverage(const reference::Range& r,
                                                           int32_t noDataValue) const
    {
        checkRange(r);

        if (!file_->hasWeightedCoverage_)
            return noDataValue;

        int32_t val = std::numeric_limits<int32_t>::min();
        bool hasData = false;
        for (size_t ii = r.begin_; ii < r.end_; ++ii)
        {
            const DataItem& d = buffer_[ii - bufferStartOffset_];
            if (d.hasData_)
            {
                val = std::max(val, d.weightedCoverage_);
                hasData = true;
            }
        }
        if (hasData)
            return val;
        else
            return noDataValue;
    }

    int32_t ReferenceSupportReader::getMaxUniqueSequenceCoverage(
        const reference::Range& r, int32_t noDataValue) const
    {
        checkRange(r);

        int32_t val = std::numeric_limits<int32_t>::min();
        bool hasData = false;
        for (size_t ii = r.begin_; ii < r.end_; ++ii)
        {
            const DataItem& d = buffer_[ii - bufferStartOffset_];
            if (d.hasData_)
            {
                val = std::max(val, d.uniqueCoverage_);
                hasData = true;
            }
        }
        if (hasData)
            return val;
        else
            return noDataValue;
    }


} } // cgatools::cgdata
