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

#ifndef CGA_TOOLS_SAM_RECORD_HPP_
#define CGA_TOOLS_SAM_RECORD_HPP_ 1

//! @file SamRecord.hpp

#include "cgatools/core.hpp"
#include "SamOptions.hpp"

#include <vector>
#include <string>
#include <map>
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace reference {
    class CrrFile;
}}

namespace cgatools { namespace mapping {

    typedef std::pair<uint16_t,uint16_t> UInt16Pair;

    class SamRecord
    {
        friend std::ostream& operator<< (std::ostream& ost, const SamRecord& r);
    public:
        enum RecordType
        {
            DEFAULT = 0,
            BASE_MAPPING,
            EVIDENCE,
            EVIDENCE_CACHE
        };

        typedef std::vector<SamRecord *> SamRecords;
        typedef std::vector<const SamRecord *> ConstSamRecords;

        SamRecord(
            const std::string& readName, bool isMapped, bool onNegativeStrand, bool isPrimary,
            uint8_t side, uint16_t chr, int32_t position,
            const std::string& extCigar,
            uint8_t mappingQuality,
            bool    isConsistentMapQ, 
            const std::string& fullReadSequence,
            const std::string& fullReadScores,
            UInt16Pair sequenceStartAndLength
            );

        virtual ~SamRecord() {}

        bool correctPosition(const reference::CrrFile& reference);

        RecordType  typeId_;

        std::string readName_;
        bool        isMapped_;
        bool        onNegativeStrand_;
        bool        isPrimary_;
        bool        isGroupPrimary_;
        uint8_t     side_;
        uint16_t    chr_;
        int32_t     position_;
        std::string extCigar_;

        std::string fullReadSequence_;
        std::string fullReadScores_;
        bool        isSvCandidate_;

        // start and length of the sequence corresponding to this mapping(0) or the mate (1)
        UInt16Pair sequenceStartAndLength_[2];

        SamRecords  mates_;
        SamRecords  alternatives_;

        uint8_t     getMappingQuality() const {return mappingQuality_;}
        void        setMappingQuality(uint8_t value, bool isConsistent) {
            mappingQuality_ = value; isConsistentMapQ_ = isConsistent;
        }

        //the flag is used for correct mate output in the case only one arm is mapped
        bool        isArtificialMateReported() const {return isArtificialMateReported_;}
        void        setArtificialMateReported(bool isReported) {isArtificialMateReported_ = true;}

        bool        isConsistent() const {return isConsistentMapQ_;}
    protected:
        uint8_t     mappingQuality_;
        bool        isConsistentMapQ_;
        bool        isArtificialMateReported_;
    };
    
    //! Used to manage sections of the SAM file header
    //! use "operator <<" to serialize the content into a SAM file
    class SamFileHeaderBlock
    {
        friend std::ostream& operator<< (std::ostream& ostr,const SamFileHeaderBlock& block);
    public:
        typedef std::vector<SamFileHeaderBlock> Children;

        SamFileHeaderBlock(std::string id = "None")
            :id_(id),type_(id),separator_("\t")
        {}
        SamFileHeaderBlock(const std::string &id, 
            const std::string& type, const std::string& separator, const std::string& value)
            :id_(id),type_(type),separator_(separator),value_(value)
        {}
        SamFileHeaderBlock(
            const std::string &id, const std::string& separator, const std::string& value)
            :id_(id),type_(id),separator_(separator),value_(value)
        {}

        const SamFileHeaderBlock& getChild(const std::string &child) const;

        //removes the children of the same type but with different ids
        void  removeNotMatchingChildren(const SamFileHeaderBlock &child);

        bool operator== (const SamFileHeaderBlock& other) const {
            return id_ == other.id_ && type_ == other.type_;
        }

        SamFileHeaderBlock& get(const SamFileHeaderBlock& b);
        SamFileHeaderBlock& add(const SamFileHeaderBlock& b);

        std::string id_;
        std::string type_;
        std::string separator_;
        std::string value_;
        Children children_;
    };

    class EvidenceSamRecord;
    class SamSequenceSplitter;

    //! Generates a SAM file record based on a given mapping record
    //! The final logic defining what goes into the output file is defined here
    class SamRecordGenerator
    {
    public:
        class OutputFileDescriptor
        {
        public:
            OutputFileDescriptor(std::ostream &  outStream)
                :hasHeader_(false), hasRecords_(false), outStream_(outStream)
            {}
            ~OutputFileDescriptor();

            void writeHeader(const SamFileHeaderBlock& header);

            bool            hasHeader_;
            bool            hasRecords_;
            std::ostream &  outStream_;
        };

        static const char SAM_SEPARATOR = '\t';

        SamRecordGenerator(
            std::ostream& outSamFile, 
            const reference::CrrFile& reference,
            const SamGeneratorConfig& config,
            const std::vector<std::string>& outStreams
            );

        ~SamRecordGenerator();

        void mappingRecordToSam(const SamRecord& record);

        void setHeader(const SamFileHeaderBlock& header);
    protected:
        //obsolete function to detect consistency
        bool isConsistent(const SamRecord &r) const;

        OutputFileDescriptor& getOutputStream(const std::string &id);

        //! Generate additional tags: mate sequence
        void printMateSequence(OutputFileDescriptor& out, 
            const std::string &mateSeq, const std::string &mateScore);

        void printReadGroup(OutputFileDescriptor& out);
        void printNegativeGapTag(OutputFileDescriptor& out, const mapping::SamSequenceSplitter &splitter);
    
        //! Generate inconsistent mapping tag
        void flagAsSVCandidate(OutputFileDescriptor& out);
        //! print allele info for EvidenceRecords
        void printAlleleInfoTag(OutputFileDescriptor& out, const mapping::EvidenceSamRecord& evidenceRecord);
        //! print alternative mappings or mates
        void printAlternatives(OutputFileDescriptor& out, const std::string& tag, 
            const mapping::SamRecord::SamRecords& samRecords, size_t startFrom);

        uint32_t getAdjustedSamPosition(const SamRecord& record ) const;
        std::string getSamChr(uint16_t chr) const;

        typedef std::map<std::string,boost::shared_ptr<OutputFileDescriptor> > OutStreamMap;
        OutStreamMap outFiles_;

        typedef boost::shared_ptr<std::ostream> OpenStreamPtr;
        typedef std::vector<OpenStreamPtr> OutStreams;
        OutStreams  openStreams_;

        const reference::CrrFile&   reference_;
        SamFileHeaderBlock          header_;
        std::string                 readGroup_;
        const SamGeneratorConfig &  samGeneratorConfig_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_SAM_RECORD_HPP_
