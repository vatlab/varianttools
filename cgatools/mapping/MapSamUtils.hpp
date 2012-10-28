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

#ifndef CGA_TOOLS_COMMAND_MAP_SAM_UTILS_HPP_
#define CGA_TOOLS_COMMAND_MAP_SAM_UTILS_HPP_ 1

//! @file MapSamUtils.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/mapping/SamRecord.hpp"
#include "cgatools/cgdata/Dnb.hpp"

#include <vector>
#include <boost/array.hpp>

namespace cgatools { namespace util {
    class DelimitedFile;
}};

namespace cgatools { namespace mapping {

    class ReadsFlagParser {
    public:
        ReadsFlagParser():flags_(0) {}
        ReadsFlagParser(size_t flags):flags_(flags) {}
        bool HalfDnbNoMatches(size_t side) const {return side==0 ? (flags_ & 0x1)!=0:(flags_ & 0x4)!=0;}
        bool HalfDnbNoMapOverflow(size_t side) const {return side==0 ? (flags_ & 0x2)!=0:(flags_&0x8)!=0;}
        bool NoHalfDnbMappings(size_t side) const {return side==0 ? (flags_ & 0x3)!=0:(flags_ & 0xC)!=0;}
        bool NoMappings() const {return (flags_ & 0x3)!=0 && (flags_ & 0xC)!=0;}

        size_t flags_;
    };

    class ReadsRecord {
        friend std::ostream& operator<< (std::ostream& out, const ReadsRecord& r);
    public:
        ReadsRecord():recordIndex_(-1) {}

        void initParser(util::DelimitedFile &delimitedFile);

        std::string     reads_;
        std::string     scores_;
        ReadsFlagParser flags_;

        //additional fields
        int             recordIndex_;
    };

    class MappingsFlagParser {
    public:
        MappingsFlagParser():flags_(0) {}
        MappingsFlagParser(size_t flags):flags_(flags) {}
        bool    LastDnbRecord() const {return (flags_ & 0x1) != 0;}
        void    setLastDnbRecord() {flags_ |= 0x1;}
        size_t  getSide() const {return (flags_ & 0x2) > 0;}
        void    setSide(size_t side) {CGA_ASSERT(side<2); flags_ |= side<<1;}
        size_t  getStrand() const {return (flags_ & 0x4) > 0;}
        void    setStrand(size_t strand) {CGA_ASSERT(strand<2); flags_ |= strand<<2;}

        size_t flags_;
    };

    class MappingsRecord {
        friend std::ostream& operator<< (std::ostream& out, const MappingsRecord& r);
    public:
        MappingsRecord()
            : flags_(0), chr_(""), offsetInChr_(-1), weightChar_(33),
            bestMate_(-1), isPrimary_(false), recordIndex_(-1)
        {
            gaps_.assign(0);
        }

        void initParser(util::DelimitedFile &delimitedFile);

        std::string createCigar(const cgdata::HalfDnbStructure::Reads &readLengths) const;

        MappingsFlagParser  flags_;
        std::string         chr_;
        int                 offsetInChr_;
        boost::array<int,3> gaps_;
        uint8_t             weightChar_;
        int                 bestMate_;

        //additional fields
        bool                isPrimary_;
        int                 recordIndex_;
    };

    typedef std::vector<MappingsRecord> MappingsRecords;

    class BaseMappingSamRecord : public SamRecord
    {
    public:
        BaseMappingSamRecord(
            const std::string& readName,
            const MappingsRecord &mapping,
            const cgdata::DnbStructure& dnbStructure,
            const std::string& fullReadSequence,
            const std::string& fullReadScores,
            const reference::CrrFile& reference
            );
    };


} } // cgatools::mapping

#endif // CGA_TOOLS_COMMAND_MAP_SAM_UTILS_HPP_
