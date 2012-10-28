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

#ifndef CGA_TOOLS_COMMAND_MAP2SAM_CONVERTER_HPP_
#define CGA_TOOLS_COMMAND_MAP2SAM_CONVERTER_HPP_ 1

//! @file Map2SamConverter.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/RangeSet.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "MapSamUtils.hpp"

#include <boost/array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <string>
#include <map>

namespace cgatools { namespace mapping {

    class Map2SamConfig {
    public:
        Map2SamConfig() 
            :recordsFrom_(0)
            ,recordsTo_(std::numeric_limits<size_t>::max())
            ,dumpDebugInfo_(true)
            ,caveat_(false)
        {}

        size_t  recordsFrom_;
        size_t  recordsTo_;
        bool    dumpDebugInfo_;
        bool caveat_;

        util::StringVector exportRegions_;

        std::string exportRootDirName_;
        std::string inputReadsBatchId_;
        std::string inputReadsFileName_;
        std::string inputMappingsFileName_;
        std::string referenceFileName_;

        std::string outputFileName_;
        util::StringVector outputStreamNames_;

        std::string commandLine_;

        SamGeneratorConfig  samGeneratorConfig_;
    };

    class LibraryData;

    class Map2SamConverter 
    {
    public:
        typedef boost::shared_ptr<std::istream> InStream;
        typedef boost::ptr_vector<SamRecord> SamRecordArray;
        static const size_t MAX_SIDES = 2; //max number of mates in a DNB

        Map2SamConverter(const Map2SamConfig &config, std::ostream &outSamFile);
        virtual ~Map2SamConverter() {}

        virtual void init();
        void run();

    protected:
        SamFileHeaderBlock createHeader();

        //Export mapping record in SAM format
        virtual void writeMappingRecord(const SamRecord &m) const;

        void processDnbRecord(const ReadsRecord& readsRecord, 
            const mapping::MappingsRecords& mappingsRecords) const;

        void outputSamMappings(const mapping::ReadsRecord& readsRecord, 
            SamRecordArray& samMappings) const;

        std::string generateDnbId(const mapping::ReadsRecord& readsRecord) const;

        void convertBaseMappingsIntoSamMappings(const mapping::ReadsRecord& readsRecord, 
            SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const;

        //! a template function allowing functionality extension by a derived class, 
        //! such as merge additional mappings
        //! returns true if the processing has replaced the standard processing 
        virtual bool processMappings(const mapping::ReadsRecord& readsRecord, 
            SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const = 0;

        reference::Range getMappingRange(const SamRecord &mapping) const;

        //! detect the chunk number for the older versions
        size_t getChunkNumber(const std::string &fileName, size_t formatVersion) const;

        size_t batchNumber_;
        size_t formatVersion_;
        std::string slide_;
        std::string lane_;
        std::string laneId_;

        InStream readsFileStream_;
        InStream mappingsFileStream_;

        boost::shared_ptr<util::DelimitedFile> readsFile_;
        boost::shared_ptr<util::DelimitedFile> mappingsFile_;

        boost::shared_ptr<LibraryData>              library_;

        reference::CrrFile      reference_;

        std::ostream &          outSamFile_;

        const Map2SamConfig &   config_;

        boost::scoped_ptr<util::FastRangeSet>           exportRegions_;
        boost::scoped_ptr<mapping::SamRecordGenerator>  mappingSamRecordGenerator_;
    };

    class BaseMap2SamConverter : public Map2SamConverter
    {
    public:
        BaseMap2SamConverter(const Map2SamConfig &config, std::ostream &outSamFile)
            : Map2SamConverter(config, outSamFile) {}

    protected:
        //! a template function allowing functionality extension by a derived class, 
        //! such as merge additional mappings
        //! returns true if the processing has replaced the standard processing 
        virtual bool processMappings(const mapping::ReadsRecord& readsRecord, 
            SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const;

        size_t detectPrimaryMapping(const MappingsRecords &mappingsRecords, 
            bool oneArmOnly, SamRecordArray& samMappings) const;
    };


} } // cgatools::mapping

#endif // CGA_TOOLS_COMMAND_MAP2SAM_CONVERTER_HPP_
