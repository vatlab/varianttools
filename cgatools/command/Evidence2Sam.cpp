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
#include "cgatools/command/Evidence2Sam.hpp"
#include "cgatools/mapping/EvidenceSamUtil.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/CrrFile.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace cgatools { namespace command {

    //! Contains an "LRU cache"-like structure: a list of records + an indexing multimap
    class EvidenceDnbBuffer {
    public:
        static const int32_t LOOK_AHEAD_DISTANCE = 1000;

        typedef std::list<mapping::EvidenceDnbRecord>        EvidenceDnbs;
        typedef std::vector<EvidenceDnbs::iterator> EvidenceDnbSet;

        class IntervalDescriptor {
        public:
            IntervalDescriptor(size_t id = 0)
                : id_(id)
                ,minPosition_(std::numeric_limits<int>::max())
                ,maxPosition_(std::numeric_limits<int>::min())
            {}
            int32_t id_;
            int     minPosition_;
            int     maxPosition_;
        };

        typedef std::list<IntervalDescriptor>   Intervals;
        typedef std::multimap<uint64_t, EvidenceDnbs::iterator> EvidenceDnbIndex;

        uint64_t getDnbId(const mapping::EvidenceDnbRecord &evidenceRecord) { 
            //slideId - 5 decimal digits
            //laneId - 2 decimal digits
            //lanePartId - 2 decimal digits
            //dnbIdInAPart - 10 decimal digits
            //--------------------------------
            // = 19 decimal digits < log10(2^64)
            return  (util::parseValue<int>(evidenceRecord.slide_.substr(2,5)
                + evidenceRecord.lane_.substr(1,2)) * 100
                + evidenceRecord.fileNumInLane_) * 10000000000
                + evidenceRecord.dnbOffsetInLaneFile_;
        }

        //! loads DNBs interval by interval into the buffer until
        //! the desired window size is loaded or EOF is reached
        //! returns true if EOF has been reached otherwise false
        bool updateBuffer(util::DelimitedFile& file, mapping::EvidenceDnbRecord &record) {
            while (intervals_.empty() ||
                intervals_.back().minPosition_ - intervals_.front().maxPosition_ < LOOK_AHEAD_DISTANCE) 
            {
                intervals_.push_back(IntervalDescriptor(record.intervalId_));
                if (loadInterval(file,record,intervals_.back()))
                    return true;
            }
            return false;
        }

        //! loads DNBs of the next interval into the buffer
        //! returns true if EOF has been reached otherwise false
        bool loadInterval(util::DelimitedFile& file, 
            mapping::EvidenceDnbRecord &evidenceRecord, 
            IntervalDescriptor &interval) 
        {
            bool result = false;

            EvidenceDnbs::iterator it;

            while(interval.id_ == evidenceRecord.intervalId_) {
                uint64_t dnbId = getDnbId(evidenceRecord);

                // update interval bounds
                if (evidenceRecord.offsetInReference_[0] < interval.minPosition_)
                    interval.minPosition_ = evidenceRecord.offsetInReference_[0];
                if (evidenceRecord.offsetInReference_[0] > interval.maxPosition_)
                    interval.maxPosition_ = evidenceRecord.offsetInReference_[0];

                it = evidenceDnbs_.insert(evidenceDnbs_.end(),evidenceRecord);
                evidenceDnbIndex_.insert(std::make_pair(dnbId,it));

                if (!file.next()) {
                    result = true;
                    break;
                }
            }

            return result;
        }

        void removeLastInterval(int32_t intervalId) {
            CGA_ASSERT(!intervals_.empty());
            CGA_ASSERT_EQ(intervals_.front().id_,intervalId);
            while (!evidenceDnbs_.empty() && 
                evidenceDnbs_.front().intervalId_ == intervalId) 
            {
                removeMappingRecord(evidenceDnbs_.begin());
            }
            intervals_.pop_front();
        }

        EvidenceDnbs::iterator removeMappingRecord(EvidenceDnbs::iterator record) {
            uint64_t dnbId = getDnbId(*record);
            EvidenceDnbIndex::iterator it = evidenceDnbIndex_.find(dnbId);
            CGA_ASSERT_MSG(it!=evidenceDnbIndex_.end(),"corrupted evidenceDnbIndex_");
            for (;it->second != record; ++it) {
                CGA_ASSERT_MSG(it!=evidenceDnbIndex_.end() && it->first == dnbId,
                    "corrupted evidenceDnbIndex_: " << CGA_VOUT(it->first) << CGA_VOUT(dnbId));
            }
            evidenceDnbIndex_.erase(it);
            EvidenceDnbs::iterator result = evidenceDnbs_.erase(record);
            return result;
        }

        void findDnbRecordsOfTheSameDnb(EvidenceDnbs::iterator& record, EvidenceDnbSet& result) {
            uint64_t dnbId = getDnbId(*record);
            EvidenceDnbIndex::iterator it = evidenceDnbIndex_.find(dnbId);
            CGA_ASSERT_MSG(it!=evidenceDnbIndex_.end(),"corrupted evidenceDnbIndex_");
            for (;it!=evidenceDnbIndex_.end() && it->first == dnbId; ++it) {
                result.push_back(it->second);
            }
        }

        bool empty() const {
            CGA_ASSERT_EQ(evidenceDnbs_.empty(),evidenceDnbIndex_.empty());
            return evidenceDnbs_.empty();
        }

        const IntervalDescriptor &  getLastInterval() {return intervals_.front();}
        EvidenceDnbs::iterator      getEvidenceDnbsBegin() {return evidenceDnbs_.begin();}

        const Intervals &       getIntervals() const {return intervals_;}
        const EvidenceDnbs &    getEvidenceDnbs() const {return evidenceDnbs_;}
        const EvidenceDnbIndex& getEvidenceDnbIndex() const {return evidenceDnbIndex_;}

    protected:

        Intervals           intervals_;
        EvidenceDnbs        evidenceDnbs_;
        EvidenceDnbIndex    evidenceDnbIndex_;
    };

    class FullMappingRange {
    public:
        static const uint32_t ARM_OVERLAP_LENGTH = 45;

        FullMappingRange(mapping::EvidenceDnbRecord &evidenceRecord)
        {
            //0 - left Arm, 2 - right Arm
            int index = ((evidenceRecord.side_) ^ uint8_t(evidenceRecord.strand_)) * 2; 

            armRanges_[index] = evidenceRecord.offsetInReference_[0];
            armRanges_[index+1] = evidenceRecord.offsetInReference_[0]+ARM_OVERLAP_LENGTH;
            armRanges_[2-index] = evidenceRecord.offsetInReference_[1];
            armRanges_[2-index+1] = evidenceRecord.offsetInReference_[1]+ARM_OVERLAP_LENGTH;

            CGA_ASSERT_MSG(armRanges_[0]<=armRanges_[1],
                CGA_VOUT(armRanges_[0])<<CGA_VOUT(armRanges_[1])<<CGA_VOUT(evidenceRecord));
            CGA_ASSERT_MSG(armRanges_[2]<=armRanges_[3],
                CGA_VOUT(armRanges_[2])<<CGA_VOUT(armRanges_[3])<<CGA_VOUT(evidenceRecord));
        }

        //returns true if ether left or right arms overlap
        bool correspondingArmsOverlap(const FullMappingRange &r) const {
            for(size_t i=0; i<2; ++i) {
                if (r.armRanges_[i]>=armRanges_[0] && r.armRanges_[i]<=armRanges_[1])
                    return true;
                if (r.armRanges_[i+2]>=armRanges_[2] && r.armRanges_[i+2]<=armRanges_[3])
                    return true;
            }
            return false;
        }

        boost::array<int,4> armRanges_;
    };

    class Evidence2SamConverter {
    public:
        static const char SAM_SEPARATOR = '\t';

        typedef boost::shared_ptr<std::istream> InStream;

        Evidence2SamConverter(const Evidence2SamConfig &config, std::ostream &outSamFile)
          : formatVersion_(0), 
            halfDnbSize_(35), 
            chrSequence_(NULL),
            outSamFile_(outSamFile), 
            config_(config)
        {
            init();
        }

        struct SortEvidenceDnbIteratorsByInterval {
            bool operator() (const EvidenceDnbBuffer::EvidenceDnbs::const_iterator &o1,
                             const EvidenceDnbBuffer::EvidenceDnbs::const_iterator &o2) 
            { return o1->intervalId_ < o2->intervalId_; }
        };

        struct SortEvidenceDnbIteratorsByStrand {
            bool operator() (const EvidenceDnbBuffer::EvidenceDnbs::const_iterator &o1,
                const EvidenceDnbBuffer::EvidenceDnbs::const_iterator &o2) 
            { return o1->side_ == o2->side_ ? o1->strand_ < o2->strand_ : o1->side_ < o2->side_; }
        };

        EvidenceDnbBuffer::EvidenceDnbs::iterator removeFromBuffer(
            EvidenceDnbBuffer& buffer,
            EvidenceDnbBuffer::EvidenceDnbs::iterator record)
        {
            if (config_.verboseOutput_) {
                (*deletedEvidenceRecords_) << *record << std::endl;
            }
            return buffer.removeMappingRecord(record);
        }

        void getBetterMapping(mapping::EvidenceDnbRecord &targetRecord,
                                const mapping::EvidenceDnbRecord &anotherRecord) 
        {
            if (targetRecord.mappingQuality_ > anotherRecord.mappingQuality_)
                targetRecord = anotherRecord;
        }

        void getCompositeMapping(mapping::EvidenceDnbRecord &targetRecord,
                                 const mapping::EvidenceDnbRecord &anotherRecord) 
        {
            targetRecord.offsetInReference_[1] = anotherRecord.offsetInReference_[0];
            targetRecord.referenceAlignment_[1] = anotherRecord.referenceAlignment_[0];
            targetRecord.mateMappingQuality_ = anotherRecord.mappingQuality_;
        }

        void processIntersectingMappings(
            EvidenceDnbBuffer &recordBuffer,
            EvidenceDnbBuffer::EvidenceDnbSet &intersectingDnbSet) 
        {

            if (intersectingDnbSet.size()>2) 
            {
                // get the best mapping from the same interval, remove the other
                // eliminates DNBs supporting different alleles of the same interval
                std::sort(intersectingDnbSet.begin(),intersectingDnbSet.end(),
                    SortEvidenceDnbIteratorsByInterval());
                for(int prevI=intersectingDnbSet.size()-1, i=prevI-1; i>=0; prevI=i, --i) {
                    if (intersectingDnbSet[prevI]->intervalId_==intersectingDnbSet[i]->intervalId_) 
                    {
                        getBetterMapping(*intersectingDnbSet[i],*intersectingDnbSet[prevI]);
                        removeFromBuffer(recordBuffer,intersectingDnbSet[prevI]);
                        intersectingDnbSet.erase(intersectingDnbSet.begin()+prevI);
                    }
                }
            }

            if (intersectingDnbSet.size()>=2) 
            {
                // get the best mapping from neighbor intervals, remove the other
                std::sort(intersectingDnbSet.begin(),intersectingDnbSet.end(),
                    SortEvidenceDnbIteratorsByStrand());
                for(int prevI=intersectingDnbSet.size()-1, i=prevI-1; i>=0; prevI=i, --i) {
                    if (intersectingDnbSet[prevI]->side_==intersectingDnbSet[i]->side_ &&
                        intersectingDnbSet[prevI]->strand_==intersectingDnbSet[i]->strand_) 
                    {
                        getBetterMapping(*intersectingDnbSet[i],*intersectingDnbSet[prevI]);
                        removeFromBuffer(recordBuffer,intersectingDnbSet[prevI]);
                        intersectingDnbSet.erase(intersectingDnbSet.begin()+prevI);
                    }
                }
            }

            if (intersectingDnbSet.size()==2) {
                mapping::EvidenceDnbRecord &record0 = *intersectingDnbSet.front();
                mapping::EvidenceDnbRecord &record1 = *intersectingDnbSet.back();
                if (record0.intervalId_ != record1.intervalId_ && record0.side_!=record1.side_) 
                {
                    //combine a pair of mappings that belong to different intervals into one mapping
                    getCompositeMapping(record0,record1);
                    removeFromBuffer(recordBuffer,intersectingDnbSet.back());
                    intersectingDnbSet.pop_back();
                    return;
                }
            }

            if (intersectingDnbSet.size()>1) 
            {
                if (config_.verboseOutput_) {
                    std::cout << std::endl << "*************" 
                                << intersectingDnbSet.size() 
                                << "*************" << std::endl;
                    for (size_t i=0; i<intersectingDnbSet.size(); ++i) {
                        std::cout << i << ": " << *intersectingDnbSet[i] << std::endl;
                    }
                }
                for(int i=intersectingDnbSet.size()-2; i>=0; --i) {
                    getBetterMapping(*intersectingDnbSet[i],*intersectingDnbSet.back());
                    removeFromBuffer(recordBuffer,intersectingDnbSet.back());
                    intersectingDnbSet.pop_back();
                }
                if (config_.verboseOutput_) {
                    std::cout << "-- leaving: " << *intersectingDnbSet.front() << std::endl;
                }
            }
        }

        //! collects all the intervals in the maximum mate length distance
        void processIntervals(mapping::EvidenceDnbRecord &evidenceRecord) 
        {
            EvidenceDnbBuffer recordBuffer;
            EvidenceDnbBuffer::EvidenceDnbSet dnbSet;
            EvidenceDnbBuffer::EvidenceDnbSet intersectingDnbSet;
            bool eof = recordBuffer.updateBuffer(*evidenceDnbsFile_,evidenceRecord);
            while (!recordBuffer.empty()) {
                const EvidenceDnbBuffer::IntervalDescriptor &interval = recordBuffer.getLastInterval();
                EvidenceDnbBuffer::EvidenceDnbs::iterator it=recordBuffer.getEvidenceDnbsBegin();
                CGA_ASSERT_MSG(it->intervalId_ >= interval.id_,
                    CGA_VOUT(it->intervalId_)<<CGA_VOUT(interval.id_));
                while (it->intervalId_==interval.id_) {
                    dnbSet.clear();
                    recordBuffer.findDnbRecordsOfTheSameDnb(it,dnbSet);
                    CGA_ASSERT(!dnbSet.empty());
                    if (dnbSet.size()>1) {
                        intersectingDnbSet.clear();
                        FullMappingRange currentRecordRange(*it);
                        BOOST_FOREACH(EvidenceDnbBuffer::EvidenceDnbs::iterator &mapIt, dnbSet)
                        {
                            if (mapIt!=it) 
                            {
                                if (currentRecordRange.correspondingArmsOverlap(FullMappingRange(*mapIt))
                                    || mapIt->intervalId_ == interval.id_)
                                    intersectingDnbSet.push_back(mapIt);
                            }
                        }
                        if (!intersectingDnbSet.empty()) 
                        {
                            intersectingDnbSet.push_back(it);
                            processIntersectingMappings(recordBuffer,intersectingDnbSet);
                            it = intersectingDnbSet.front();
                        }
                    }
                    generateSamRecord(*it);
                    removeFromBuffer(recordBuffer,it);
                    it = recordBuffer.getEvidenceDnbsBegin();
                    if (recordBuffer.empty())
                        break;
                }
                recordBuffer.removeLastInterval(interval.id_);
                if (!eof)
                    eof = recordBuffer.updateBuffer(*evidenceDnbsFile_,evidenceRecord);
            }

        }


        void processOneIntervalRemoveIntervalDuplicates(int32_t intervalId_,
                                mapping::EvidenceDnbRecord &evidenceRecord)
        {
            typedef std::map<uint64_t, mapping::EvidenceDnbRecord> DuplicatedDnbs;
            DuplicatedDnbs duplicatedDnbs;


            while(intervalId_ == evidenceRecord.intervalId_) {
                //slideId - 5 decimal digits
                //laneId - 2 decimal digits
                //lanePartId - 2 decimal digits
                //dnbIdInAPart - 10 decimal digits
                //--------------------------------
                // = 19 decimal digits < log10(2^64)
                uint64_t dnbId = 
                    (util::parseValue<int>(evidenceRecord.slide_.substr(2,5)
                        + evidenceRecord.lane_.substr(1,2)) * 100
                    + evidenceRecord.fileNumInLane_) * 10000000000
                    + evidenceRecord.dnbOffsetInLaneFile_;

                DuplicatedDnbs::iterator it = duplicatedDnbs.find(dnbId);
                if (it!=duplicatedDnbs.end()) { 
                    if (it->second.mappingQuality_<evidenceRecord.mappingQuality_)
                        it->second = evidenceRecord;
                } else
                    duplicatedDnbs[dnbId] = evidenceRecord;

                if (!evidenceDnbsFile_->next())
                    break;
            }

            BOOST_FOREACH(DuplicatedDnbs::value_type &p, duplicatedDnbs) {
                generateSamRecord(p.second);
            }

        }

        bool processOneInterval(int32_t intervalId_, mapping::EvidenceDnbRecord &evidenceRecord) 
        {

            while(intervalId_ == evidenceRecord.intervalId_) {
                generateSamRecord(evidenceRecord);
                if (!evidenceDnbsFile_->next())
                    return false;
            }
            return true;
        }

        void run() {
            samRecordGenerator_->setHeader(createHeader());

            mapping::EvidenceDnbRecord evidenceRecord;
            evidenceRecord.initParser(*evidenceDnbsFile_, formatVersion_, reference_);

            if (config_.keepDuplicates_) {
                for (bool notEof = evidenceDnbsFile_->next(); notEof;
                    notEof = processOneInterval(evidenceRecord.intervalId_,evidenceRecord)) 
                    ;
            } else {
                if (evidenceDnbsFile_->next()) // input file is not empty
                    processIntervals(evidenceRecord);
            }
        }


    protected:
        mapping::SamFileHeaderBlock createHeader() {
            const util::DelimitedFile::Metadata& readsHeader = evidenceDnbsFile_->getMetadata();

            mapping::SamFileHeaderBlock h("");

            h.get(mapping::SamFileHeaderBlock("@HD","@HD","",""))
                .add(mapping::SamFileHeaderBlock("VN","VN","\t","1.4"));

            BOOST_FOREACH(const reference::CompactDnaSequence &s, reference_.listChromosomes()) 
            {
                h.get(mapping::SamFileHeaderBlock(s.getName(),"@SQ","\n",""))
                    .add(mapping::SamFileHeaderBlock("SN","\t",s.getName()))
                    .add(mapping::SamFileHeaderBlock("LN","\t",boost::lexical_cast<std::string>(s.length())))
                    .add(mapping::SamFileHeaderBlock("AS","\t",assemblyId_))
                    .add(mapping::SamFileHeaderBlock("UR","\t",config_.referenceFileName_));
            }

            std::string date = readsHeader.get("GENERATED_AT");
            date = date.substr(0,date.find(' '));
            boost::gregorian::date dd(boost::gregorian::from_simple_string(date));

            h.get(mapping::SamFileHeaderBlock("@RG","@RG","\n",""))
                .add(mapping::SamFileHeaderBlock("ID","\t",evidenceFileId_))
                .add(mapping::SamFileHeaderBlock("SM","\t",
                    (readsHeader.hasKey("SAMPLE") ? readsHeader.get("SAMPLE"):"not-provided") ))
                //.add(SamFileHeaderBlock("LB","\t",readsHeader.get("LIBRARY"))) - may be from different libs
                .add(mapping::SamFileHeaderBlock("PU","\t",evidenceFileId_))
                .add(mapping::SamFileHeaderBlock("CN","\t","\"Complete Genomics\""))
                .add(mapping::SamFileHeaderBlock("DT","\t",boost::gregorian::to_iso_extended_string(dd)))
                .add(mapping::SamFileHeaderBlock("PL","\t","\"Complete Genomics\""));

            h.get(mapping::SamFileHeaderBlock("@PG","@PG","\n",""))
                .add(mapping::SamFileHeaderBlock("ID","\t","cgatools"))
                .add(mapping::SamFileHeaderBlock("VN","\t",CGA_TOOLS_VERSION))
                .add(mapping::SamFileHeaderBlock("CL","\t",config_.commandLine_));

            return h;
        }

        void generateSamRecord(mapping::EvidenceDnbRecord &evidenceRecord) 
        {
            if (config_.verboseOutput_) {
                (*selectedEvidenceRecords_) << evidenceRecord << std::endl;
            }
            const size_t sides = 2;

            evidenceRecord.adjustOffset(
                reference_.listChromosomes()[evidenceRecord.chromosome_].length());

            boost::ptr_vector<mapping::SamRecord> records;
            for (size_t i=0; i<sides; ++i) 
                records.push_back(new mapping::EvidenceSamRecord(
                    evidenceRecord, !config_.keepDuplicates_, i, halfDnbSize_, reference_));

            for (size_t i=0; i<sides; ++i)
                records[i].mates_.push_back(&records[1-i]);

            for (size_t i=0; i<sides; ++i) {
                int startPos = evidenceRecord.offsetInReference_[i];
                mapping::Cigar c(evidenceRecord.referenceAlignment_[i]);
                reference::Range range(evidenceRecord.chromosome_,
                    startPos,startPos+c.getReferenceLength());
                if (!exportRegions_->intersects(range))
                    continue;

                samRecordGenerator_->mappingRecordToSam(records[i]);
            }
        }


        void init() {
            evidenceDnbsFileStream_ = util::InputStream::
                openCompressedInputStreamByExtension(config_.evidenceDnbsFileName_);
            evidenceDnbsFile_.reset(new util::DelimitedFile(*evidenceDnbsFileStream_,
                                                            config_.evidenceDnbsFileName_));

            const util::DelimitedFile::Metadata& evidenceFileHeader = evidenceDnbsFile_->getMetadata();

            formatVersion_ = evidenceFileHeader.getFormatVersion();
            std::string chrName = evidenceFileHeader.get("CHROMOSOME");
            if (evidenceFileHeader.hasKey("ASSEMBLY_ID")) {
                assemblyId_ = evidenceFileHeader.get("ASSEMBLY_ID");
            } else {
                std::string fname = boost::filesystem::path(config_.evidenceDnbsFileName_).leaf();
                std::string fprefix = "evidenceDnbs-"+chrName+"-";
                size_t pos=fname.find(fprefix);
                if (pos==std::string::npos)
                    throw util::Exception("Cannot extract assembly ID from file: "+fname);
                pos+=fprefix.length();
                size_t endpos = fname.find("-", pos);
                if (endpos==std::string::npos)
                    throw util::Exception("Cannot extract assembly ID from file: "+fname);
                endpos = fname.find("-", endpos);
                if (endpos==std::string::npos)
                    throw util::Exception("Cannot extract assembly ID from file: "+fname);
                assemblyId_ = fname.substr(pos,endpos-pos);
            }
            evidenceFileId_ = assemblyId_+"-"+chrName;

            if (!config_.referenceFileName_.empty())
                reference_.open(config_.referenceFileName_);

            //define genomic ranges. If no ranges were provided - use whole region
            exportRegions_.reset(new util::FastRangeSet(reference_));
            if (config_.exportRegions_.empty())
                exportRegions_->addWholeReference();
            else
                exportRegions_->add(config_.exportRegions_);

            std::cerr << "Export the sequence: " << chrName 
                << ", length: " << reference_.listChromosomes()[reference_.getChromosomeId(chrName)].length()
                << std::endl;

            if (config_.verboseOutput_) {
                CGA_ASSERT_MSG(!config_.outputFileName_.empty(),"The output file name must be provided");
                selectedEvidenceRecords_.reset(new util::OutputStream(config_.outputFileName_+".sel"));
                deletedEvidenceRecords_.reset(new util::OutputStream(config_.outputFileName_+".del"));
            }

            samRecordGenerator_.reset(
                new mapping::SamRecordGenerator(outSamFile_,reference_,
                    config_.samGeneratorConfig_,std::vector<std::string>()));
        }

        const reference::CompactDnaSequence& getChrSequence( const std::string& chrName )
        {
            if (chrSequence_==NULL)
                chrSequence_ = &reference_.listChromosomes()[reference_.getChromosomeId(chrName)];
            CGA_ASSERT_EQ(chrSequence_->getName(),chrName);
            return *chrSequence_;
        }


        size_t      formatVersion_;
        std::string assemblyId_;
        std::string evidenceFileId_;

        size_t   halfDnbSize_;

        InStream evidenceDnbsFileStream_;

        boost::shared_ptr<util::DelimitedFile>          evidenceDnbsFile_;

        reference::CrrFile                              reference_;
        const reference::CompactDnaSequence *           chrSequence_;

        std::ostream &                                  outSamFile_;
        boost::shared_ptr<std::ostream>                 selectedEvidenceRecords_;
        boost::shared_ptr<std::ostream>                 deletedEvidenceRecords_;

        const Evidence2SamConfig &                      config_;

        boost::scoped_ptr<util::FastRangeSet>               exportRegions_;
        boost::scoped_ptr<mapping::SamRecordGenerator>  samRecordGenerator_;
    };

    Evidence2Sam::Evidence2Sam(const std::string& name)
        : Command(name,
                  "Converts CGI variant evidence data into SAM format.",
                  "0.3 or later",
        "The evidence2sam converter takes as input evidence mapping files (evidenceDnbs-*) "
        "and generates one SAM file as an output. The output is sent into stdout by default. "
        "By default, all the evidence mapping records from the input are "
        "converted into a pair of corresponding SAM records - one record for each HalfDNB. "
        "The negative gaps in CGI mappings are represented using GS/GQ/GC tags."
        )
    {
        options_.add_options()
            ("evidence-dnbs,e", po::value<std::string>(&config_.evidenceDnbsFileName_),
             "Input evidence dnbs file.")
            ("reference,s", po::value<std::string>(&config_.referenceFileName_),
             "Reference file.")
            ("output,o", po::value<std::string>(&config_.outputFileName_)->default_value("STDOUT"),
             "The output SAM file (may be omitted for stdout).")
            ("extract-genomic-region,r", po::value<util::StringVector>(&config_.exportRegions_),
             "defines a region as a half-open interval 'chr,from,to'.")
            ("keep-duplicates", po::bool_switch(&config_.keepDuplicates_)->default_value(false),
             "Keep local duplicates of DNB mappings."
             "All the output SAM records will be marked as not primary if this option is used.")
            ("add-allele-id", po::bool_switch(
                &config_.samGeneratorConfig_.printAlleleInfo_)->default_value(false),
             "Generate interval id and allele id tags.")
             ;

        addSamConfigOptions(config_.samGeneratorConfig_, options_);

        hiddenOptions_.add_options()
            ("debug-output,v", po::bool_switch(&config_.verboseOutput_)->default_value(false),
            "Generate verbose debug output. Please don't rely on this option in production.")
            ;
    }

    int Evidence2Sam::run(po::variables_map& vm)
    {
        //requireParam(vm, "evidence-dnbs");
        requireParam(vm, "reference");

        config_.commandLine_ = getCommandLine();

        Evidence2SamConverter convertor(config_, openStdout(config_.outputFileName_));
        convertor.run();

        return 0;
    }


} } // cgatools::command
