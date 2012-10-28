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

#ifndef CGA_TOOLS_EVIDENCECACHE_HPP_
#define CGA_TOOLS_EVIDENCECACHE_HPP_ 1

//! @file EvidenceCache.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/RangeSet.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"
#include "cgatools/mapping/LaneBatchCache.hpp"

namespace cgatools { namespace cgdata {
    class GenomeMetadata;
}}

namespace cgatools { namespace reference {
    class CrrFile;
}}

namespace cgatools { namespace util {
    class DelimitedFile;
}}


namespace cgatools { namespace mapping {

class CacheOutStreams : public OutLaneBatchStreams
{
public:
    CacheOutStreams(const boost::filesystem::path& outputPrefix)
        :OutLaneBatchStreams(outputPrefix, "EVIDENCE_CACHE")
    {}

protected:
    virtual void writeMetadata(std::ostream &output,
        const std::string& slide, const std::string& lane, size_t batchNo) const;
};


class LibraryMetadataContainer;

class EvidenceCacheBuilder 
{
public:
    EvidenceCacheBuilder(const cgdata::GenomeMetadata& genomeMetadata,const reference::CrrFile& reference);

    void processChrData(uint16_t chr, const boost::filesystem::path& outputPrefix, 
        const util::FastRangeSet::RangeSet &exportRanges);

    void exportRanges(const boost::filesystem::path& outputPrefix, const util::FastRangeSet& ranges);

protected:
    boost::shared_ptr<LibraryMetadataContainer> libraries_;

    const cgdata::GenomeMetadata&           genomeMetadata_;
    const reference::CrrFile&               reference_;
};

class EvidenceCacheDnbRecord : public cgdata::EvidenceReader::DnbRecord
{
public:
    double                   alleleConcordance_;
};

typedef std::multimap<uint64_t,EvidenceCacheDnbRecord> BatchRecords;

class EvidenceCacheReader
{
public:
    typedef std::vector<boost::filesystem::path>    InputBatchFiles;
    typedef boost::ptr_map<uint64_t,InputBatchFiles>  InputLaneBatches;

    EvidenceCacheReader(const boost::filesystem::path& rootDir)
        : inputBatches_(rootDir)
    {
        inputBatches_.collectFiles();
    }

    void readBatchRecords(const std::string& slide, const std::string& lane, size_t batchNo, 
                                            BatchRecords& result, const reference::CrrFile& crr);

protected:

    static void initCacheRecordParser(
        util::DelimitedFile& delimitedFile, EvidenceCacheDnbRecord& record, const reference::CrrFile& crr);

    InLaneBatchStreams      inputBatches_;
};


} } // cgatools::mapping

#endif // CGA_TOOLS_EVIDENCECACHE_HPP_
