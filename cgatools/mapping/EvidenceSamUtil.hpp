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

#ifndef CGA_TOOLS_EVIDENCE_SAM_UTIL_HPP_
#define CGA_TOOLS_EVIDENCE_SAM_UTIL_HPP_ 1

//! @file EvidenceSamUtil.hpp

#include "cgatools/core.hpp"
#include "cgatools/mapping/Cigar.hpp"
#include "cgatools/mapping/SamRecord.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"

#include <vector>
#include <sstream>
#include <boost/array.hpp>

namespace cgatools { namespace util {
    class DelimitedFile;
}};

namespace cgatools { namespace mapping {

    class SamSequenceSplitter {
    public:
        SamSequenceSplitter(const Cigar& extCigar, const std::string& dnbSequence, 
            const std::string& dnbScore, bool negativeStrand);

        Cigar               srcExtCigar_;
        std::string         srcSequence_;
        std::string         srcScores_;

        Cigar               cigar_;
        std::stringstream   sequence_;
        std::stringstream   scores_;

        Cigar               tagCigar_;
        std::stringstream   tagSequence_;
        std::stringstream   tagScores_;
    protected:
        void convert();
    };

    class EvidenceDnbRecord : public cgdata::EvidenceReader::DnbRecord 
    {
        friend std::ostream& operator<< (std::ostream& out, const EvidenceDnbRecord& r);
    public:
        EvidenceDnbRecord() : mateMappingQuality_(-1)
        {}
        void initParser(util::DelimitedFile &delimitedFile, size_t formatVersion, 
                                const reference::CrrFile& crr);
        //! The function is relevant only for circular chromosomes and 
        //! is used to set the offset within the chromosome dimensions
        void adjustOffset(int chrLength);

        int16_t             mateMappingQuality_;
    };

    class EvidenceCacheDnbRecord;

    class EvidenceSamRecord : public SamRecord
    {
    public:
        EvidenceSamRecord(
            const EvidenceDnbRecord& evidenceRecord,
            bool isPrimary,
            uint8_t mate,
            size_t halfDnbSize,
            const reference::CrrFile& reference
            );

        EvidenceSamRecord(
            const EvidenceCacheDnbRecord& evidenceRecord,
            const std::string& sequence,
            const std::string& scores,
            bool isPrimary,
            uint8_t mate,
            size_t halfDnbSize,
            const reference::CrrFile& reference
            );

        size_t  intervalId_;
        int     alleleIndex_;
        double  armWeight_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_EVIDENCE_SAM_UTIL_HPP_
