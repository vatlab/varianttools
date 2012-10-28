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

#ifndef CGA_TOOLS_LANE_BATCH_CACHE_HPP_
#define CGA_TOOLS_LANE_BATCH_CACHE_HPP_ 1

//! @file LaneBatchCache.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"

#include <boost/filesystem/path.hpp>
#include <boost/ptr_container/ptr_map.hpp>

#include <vector>

namespace cgatools { namespace util {
    class DelimitedFileMetadata;
}}

namespace cgatools { namespace mapping {

class BaseLaneBatchStreams
{
public:
    static std::string getBatchStreamKey(const std::string& slide, const std::string& lane, size_t batchNo);
};

//! creates a family of lane/batch files and generates standard metadata header
class OutLaneBatchStreams : public BaseLaneBatchStreams
{
public:
    typedef boost::ptr_map<std::string,util::CompressedOutputStream> OutputLaneBatches;

    OutLaneBatchStreams(const boost::filesystem::path& outputPrefix, const std::string & fileType)
        :outputPrefix_(outputPrefix), fileType_(fileType)
    {}

    OutLaneBatchStreams(
        const boost::filesystem::path& outputPrefix, 
        const std::string & fileType,
        const std::string & headerLine
    )
        :outputPrefix_(outputPrefix), fileType_(fileType), headerLine_(headerLine)
    {}

    virtual ~OutLaneBatchStreams() {}

    std::ostream& getBatchStream(const std::string& slide, const std::string& lane, size_t batchNo);

protected:
    virtual void writeMetadata(std::ostream &output,
                    const std::string& slide, const std::string& lane, size_t batchNo) const;

    virtual void addMetadata(util::DelimitedFileMetadata &meta) const {}

    boost::filesystem::path outputPrefix_;
    const std::string       fileType_;
    const std::string       headerLine_;
    OutputLaneBatches       outputLanes_;
};

class InLaneBatchStreams : public BaseLaneBatchStreams
{
public:
    typedef std::vector<boost::filesystem::path>        InputBatchFiles;
    typedef boost::ptr_map<std::string,InputBatchFiles>    InputLaneBatches;

    InLaneBatchStreams(const boost::filesystem::path& inputRoot)
        :inputRoot_(inputRoot)
    {}

    //! collects the files found in or below inputRoot and splits each batch into a separate bucket
    void collectFiles();

    const InputBatchFiles& getBatchFiles(
        const std::string& slide, const std::string& lane, size_t batchNo) const;

protected:

    boost::filesystem::path inputRoot_;
    InputLaneBatches        inputBatches_;
    InputBatchFiles         emptyBatch_;
};

} } // cgatools::mapping

#endif // CGA_TOOLS_LANE_BATCH_CACHE_HPP_
