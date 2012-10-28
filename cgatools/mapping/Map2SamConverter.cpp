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
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "Map2SamConverter.hpp"
#include "LibraryData.hpp"
#include "MapSamUtils.hpp"
#include "Cigar.hpp"

#include <boost/date_time.hpp>
#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#include <numeric>


namespace cgatools { namespace mapping {

    Map2SamConverter::Map2SamConverter( const Map2SamConfig &config, std::ostream &outSamFile ) 
        : batchNumber_(0), formatVersion_(0), outSamFile_(outSamFile), config_(config)
    {
    }

    void Map2SamConverter::run()
    {
        CGA_ASSERT_MSG(readsFileStream_.get()!=NULL,"Use Init to initialize the class");
        
        mappingSamRecordGenerator_->setHeader(createHeader());

        mapping::ReadsRecord readsRecord;
        readsRecord.initParser(*readsFile_);
        mapping::MappingsRecord mappingsRecord;
        mappingsRecord.initParser(*mappingsFile_);
        mapping::MappingsRecords baseMappings;

        int readRecordIndex = -1;
        while (readsFile_->next()) 
        {
            readsRecord.recordIndex_ = ++readRecordIndex;
            baseMappings.clear();
            if (!readsRecord.flags_.NoMappings()) 
            {
                int mappingRecordIndex = -1;
                while (mappingsFile_->next()) {
                    mappingsRecord.recordIndex_ = ++mappingRecordIndex;
                    baseMappings.push_back(mappingsRecord);
                    if (mappingsRecord.flags_.LastDnbRecord())
                        break;
                }
                CGA_ASSERT_MSG(!baseMappings.empty(),
                    "The mappings and reads files are not in sync. Excess of reads.");
            } else {
                // Don't output SAM records for the DNBs that have no mappings
                if (config_.samGeneratorConfig_.skipNotMapped_)
                    continue;
            }

            if (size_t(readRecordIndex) < config_.recordsFrom_)
                continue;
            else if (size_t(readRecordIndex) >= config_.recordsTo_)
                break;

            processDnbRecord(readsRecord, baseMappings);
        }
        CGA_ASSERT_MSG(config_.recordsTo_ != std::numeric_limits<size_t>::max() || !mappingsFile_->next(),
            "The mappings and reads files are not in sync. Excess of mappings.");
    }

    namespace {
        void mateMappings(mapping::SamRecord &m0,mapping::SamRecord &m1) 
        {
            CGA_ASSERT_NEQ(m0.side_,m1.side_);
            m0.mates_.push_back(&m1);
            m1.mates_.push_back(&m0);
        }
    }

    std::string Map2SamConverter::generateDnbId(const mapping::ReadsRecord& readsRecord) const
    {
        return laneId_ + "-" + boost::lexical_cast<std::string>(batchNumber_) 
            + ":" + boost::lexical_cast<std::string>(readsRecord.recordIndex_);
    }

    void Map2SamConverter::convertBaseMappingsIntoSamMappings(const mapping::ReadsRecord& readsRecord, 
        SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const
    {
        std::string dnbId = generateDnbId(readsRecord);

        samMappings.reserve(baseMappingRecords.size());

        //generate records for the input mappings
        BOOST_FOREACH(const mapping::MappingsRecord &m,baseMappingRecords) 
        {
            std::auto_ptr<BaseMappingSamRecord> samRecord(
                new mapping::BaseMappingSamRecord(
                    dnbId, m, library_->dnbStructure_, 
                    readsRecord.reads_, readsRecord.scores_, 
                    reference_)
                );
            if (samRecord->correctPosition(reference_))
                samMappings.push_back(samRecord);
        }
    }

    void Map2SamConverter::processDnbRecord( const mapping::ReadsRecord& readsRecord, 
        const mapping::MappingsRecords& baseMappingRecords ) const
    {
        SamRecordArray samMappings;

        if (!baseMappingRecords.empty())
            processMappings(readsRecord, samMappings, baseMappingRecords);

        outputSamMappings(readsRecord, samMappings);
    }

    //! The function does several modifications to the given samRecords
    void Map2SamConverter::outputSamMappings(const mapping::ReadsRecord& readsRecord, 
        SamRecordArray& samMappings) const
    {
        boost::array<size_t,MAX_SIDES> countBySide;
        countBySide.assign(0);

        boost::array<size_t,MAX_SIDES> primaryBySide; //sanity check
        primaryBySide.assign(0);

        bool primaryHasNoMates = false;
        bool notMatedMappingsFound = false;
        BOOST_FOREACH(mapping::SamRecord &m, samMappings) 
        {
            ++countBySide[m.side_];
            if (m.isPrimary_)
            {
                ++primaryBySide[m.side_];
                primaryHasNoMates = primaryHasNoMates || m.mates_.empty();
            }
            notMatedMappingsFound = notMatedMappingsFound || m.mates_.empty();
            BOOST_FOREACH(mapping::SamRecord &altM, samMappings)
                if (&m != &altM && m.side_ == altM.side_)
                    m.alternatives_.push_back(&altM);
        }

        size_t primaryBySideTotal = 0;
        for (size_t i=0; i<MAX_SIDES; ++i)
        {
            primaryBySideTotal += primaryBySide[i];
            CGA_ASSERT_MSG(primaryBySide[i]<=1, 
                CGA_VOUT(i)
                << CGA_VOUT(countBySide[i])
                << CGA_VOUT(primaryBySide[i])
                << CGA_VOUT(readsRecord)
            );
        }

        CGA_ASSERT_MSG(
            (primaryBySideTotal==1 && primaryHasNoMates)
            || (primaryBySideTotal==2 && !primaryHasNoMates)
            || (primaryBySideTotal==0 && samMappings.empty()), 
            "Error in primary detection for " << readsRecord
            );

        bool isSvCandidateFlag = 
            config_.samGeneratorConfig_.mateSvCandidates_ && 
            (samMappings.size()==MAX_SIDES) 
            && samMappings.front().side_!=samMappings.back().side_
            && samMappings.front().mates_.empty()
            && samMappings.back().mates_.empty()
            ;

        SamRecordArray notMappedRecords;

        if (isSvCandidateFlag) 
        {
            BOOST_FOREACH(mapping::SamRecord &m, samMappings)
                m.isPrimary_ = true;
            CGA_ASSERT(samMappings.front().mates_.empty());
            CGA_ASSERT(samMappings.back().mates_.empty());
            mateMappings(samMappings.front(),samMappings.back());
        } else if (notMatedMappingsFound || samMappings.empty()) 
        {
            //generate zero records for the missing arm mappings
            std::string dnbId = generateDnbId(readsRecord);
            mapping::MappingsRecord zeroMapping;
            zeroMapping.isPrimary_ = true;
            for (size_t i=0; i<MAX_SIDES; ++i)
            {
                zeroMapping.flags_.setSide(i);
                notMappedRecords.push_back(
                    new mapping::BaseMappingSamRecord(
                        dnbId, zeroMapping, library_->dnbStructure_, 
                        readsRecord.reads_, readsRecord.scores_, 
                        reference_)
                    );
            }

            if (samMappings.empty()) 
            {
                //mate empty records
                if (!config_.samGeneratorConfig_.skipNotMapped_)
                    mateMappings(notMappedRecords.front(),notMappedRecords.back());
            } else 
            {
                //add not mapped mates to single arm mappings
                BOOST_FOREACH(mapping::SamRecord &m, samMappings) {
                    if (m.mates_.empty())
                    {
                        mateMappings(m,notMappedRecords[1-m.side_]);
                        if (m.isPrimary_)
                            std::swap(notMappedRecords[1-m.side_].mates_.front(),
                                    notMappedRecords[1-m.side_].mates_.back());
                    }
                }
            }
        }

        if (countBySide[0]+countBySide[1] == 0)
        {
            if (!config_.samGeneratorConfig_.skipNotMapped_)
            {
                writeMappingRecord(notMappedRecords[0]);
                writeMappingRecord(notMappedRecords[1]);
            }
        } else
        {
            BOOST_FOREACH(mapping::SamRecord &m, samMappings) 
            {
                if (exportRegions_->intersects(getMappingRange(m)))
                {
                    m.isSvCandidate_ = isSvCandidateFlag;
                    if (m.isPrimary_ 
                        //&& countBySide[1-m.side_] == 0 
                        && m.mates_.size()==1 && (!m.mates_[0]->isMapped_)
                        && (!config_.samGeneratorConfig_.skipNotMapped_))
                    {
                        m.setArtificialMateReported(true);
                        writeMappingRecord(m);
                        writeMappingRecord(notMappedRecords[1-m.side_]);
                    } else 
                    {
                        writeMappingRecord(m);
                    }
                }
            }
        }
    }

    void Map2SamConverter::writeMappingRecord(const SamRecord &m) const
    {
        mappingSamRecordGenerator_->mappingRecordToSam(m);
    }


    reference::Range Map2SamConverter::getMappingRange( const mapping::SamRecord &mapping ) const
    {
        int mappingLen = mapping::Cigar(mapping.extCigar_).getReferenceLength();

        return reference::Range(mapping.chr_, mapping.position_,mapping.position_+mappingLen);
    }

    size_t Map2SamConverter::getChunkNumber( const std::string &fileName, size_t formatVersion ) const
    {
        if (formatVersion < 4)
            return 1;
        static const std::string fnamePrefix = "reads_";
        size_t lastFoundPos = std::string::npos;
        for (size_t pos=0;
            pos!=std::string::npos;
            pos=config_.inputReadsFileName_.find(fnamePrefix,pos+1))
        {
            lastFoundPos = pos;
        }
        CGA_ASSERT_MSG(lastFoundPos!=std::string::npos,
            "don't recognize the read file name format: can't find '"
            << fnamePrefix << "' in "<<fileName);
        if (config_.inputReadsFileName_.find("_",lastFoundPos+fnamePrefix.size()) != std::string::npos ||
            4 != formatVersion)
        {
            lastFoundPos = config_.inputReadsFileName_.find("_",lastFoundPos+fnamePrefix.size());
        }
        else
        {
            lastFoundPos += 5;
        }
        CGA_ASSERT_MSG(lastFoundPos!=std::string::npos,
            "don't recognize the read file name format: can't find "
            "chunk number separator '_' after '" << fnamePrefix << "' in " << fileName);
        ++lastFoundPos;
        size_t chunkNumEnd = config_.inputReadsFileName_.find(".",lastFoundPos);
        return util::parseValue<size_t>(fileName.substr(lastFoundPos,chunkNumEnd-lastFoundPos));
    }

    void Map2SamConverter::init()
    {
        readsFileStream_ = util::InputStream::
                openCompressedInputStreamByExtension(config_.inputReadsFileName_);

        mappingsFileStream_ = util::InputStream::
                openCompressedInputStreamByExtension(config_.inputMappingsFileName_);

        readsFile_.reset(new util::DelimitedFile(*readsFileStream_, config_.inputReadsFileName_));
        mappingsFile_.reset(new util::DelimitedFile(*mappingsFileStream_, config_.inputMappingsFileName_));

        const util::DelimitedFile::Metadata& readsHeader = readsFile_->getMetadata();

        formatVersion_ = readsHeader.getFormatVersion();
        if (formatVersion_ >= 1001) {
            batchNumber_ = 
                util::parseValue<size_t>(readsHeader.get("BATCH_FILE_NUMBER").c_str());
        } else {
            batchNumber_ = getChunkNumber(config_.inputReadsFileName_, formatVersion_);
            std::cerr << "for the EXP format 1.0 and older the following batch number is detected: " 
                << batchNumber_ << std::endl;
        }

        slide_ = readsHeader.get("SLIDE");
        lane_  = readsHeader.get("LANE");
        laneId_ = slide_+"-"+lane_;

        if (!config_.referenceFileName_.empty())
            reference_.open(config_.referenceFileName_);

        //define genomic ranges. If no ranges were provided - use whole region
        exportRegions_.reset(new util::FastRangeSet(reference_));
        if (config_.exportRegions_.empty())
            exportRegions_->addWholeReference();
        else
            exportRegions_->add(config_.exportRegions_);

        mappingSamRecordGenerator_.reset(
            new mapping::SamRecordGenerator(outSamFile_,reference_,
                config_.samGeneratorConfig_, config_.outputStreamNames_));

        cgdata::GenomeMetadata gm(config_.exportRootDirName_);
        library_.reset(new LibraryData(gm.getLibraryMetadata(gm.getLibraryName(laneId_)),true));
    }

    SamFileHeaderBlock Map2SamConverter::createHeader()
    {
        const util::DelimitedFile::Metadata& readsHeader = readsFile_->getMetadata();
        std::string assemblyId = (readsHeader.hasKey("ASSEMBLY_ID") 
            ? readsHeader.get("ASSEMBLY_ID") : "not defined");
        std::string fieldSize = (readsHeader.hasKey("FIELD_SIZE") 
            ? readsHeader.get("FIELD_SIZE") : "0");

        SamFileHeaderBlock h("");

        std::string refId = "notDefined";
        std::string sp = "notDefined";
        if (reference_.listChromosomes()[0].getMd5Digest().hex()=="9ebc6df9496613f373e73396d5b3b6b6")
        {
            refId = "NCBI36";
            sp = "Homo sapiens";
        } else if (reference_.listChromosomes()[0].getMd5Digest().hex()=="1b22b98cdeb4a9304cb5d48026a85128")
        {
            refId = "GRCh37";
            sp = "Homo sapiens";
        }

        h.get(SamFileHeaderBlock("@HD","@HD","","")).add(SamFileHeaderBlock("VN","VN","\t","1.3"));
        BOOST_FOREACH(const reference::CompactDnaSequence &s, reference_.listChromosomes()) 
        {
            h.get(SamFileHeaderBlock(s.getName(),"@SQ","\n",""))
                .add(SamFileHeaderBlock("SN","\t",s.getName()))
                .add(SamFileHeaderBlock("LN","\t",boost::lexical_cast<std::string>(s.length())))
                .add(SamFileHeaderBlock("AS","\t",refId))
                .add(SamFileHeaderBlock("UR","\t",config_.referenceFileName_))
                .add(SamFileHeaderBlock("M5","\t",s.getMd5Digest().hex()))
                .add(SamFileHeaderBlock("SP","\t",sp));
        }

        std::string date = readsHeader.get("GENERATED_AT");
        date = date.substr(0,date.find(' '));
        boost::gregorian::date dd(boost::gregorian::from_simple_string(date));

        h.get(SamFileHeaderBlock("@RG","@RG","\n",""))
            .add(SamFileHeaderBlock("ID","\t",laneId_))
            .add(SamFileHeaderBlock("DS","\t",fieldSize))
            .add(SamFileHeaderBlock("DT","\t",boost::gregorian::to_iso_extended_string(dd)))
            .add(SamFileHeaderBlock("LB","\t",readsHeader.get("LIBRARY")))
            .add(SamFileHeaderBlock("PU","\t",laneId_))
            .add(SamFileHeaderBlock("CN","\t","COMPLETEGENOMICS"))
            .add(SamFileHeaderBlock("PL","\t","COMPLETEGENOMICS"))
            .add(SamFileHeaderBlock("SM","\t",readsHeader.get("SAMPLE")));

        h.get(SamFileHeaderBlock("@PG","@PG","\n",""))
            .add(SamFileHeaderBlock("ID","\t","cgatools"))
            .add(SamFileHeaderBlock("VN","\t",CGA_TOOLS_VERSION))
            .add(SamFileHeaderBlock("CL","\t",config_.commandLine_));

        return h;
    }

    bool BaseMap2SamConverter::processMappings(const mapping::ReadsRecord& readsRecord, 
        SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const
    {
        convertBaseMappingsIntoSamMappings(readsRecord, samMappings, baseMappingRecords);

        if (!readsRecord.flags_.NoMappings()) 
        {
            bool isSingleArm = readsRecord.flags_.HalfDnbNoMatches(0) 
                || readsRecord.flags_.HalfDnbNoMapOverflow(0)
                || readsRecord.flags_.HalfDnbNoMatches(1)
                || readsRecord.flags_.HalfDnbNoMapOverflow(1);
            size_t primaryRecords = detectPrimaryMapping(baseMappingRecords, isSingleArm, samMappings);
            CGA_ASSERT_MSG(primaryRecords <= 2,
                CGA_VOUT(primaryRecords)
                <<CGA_VOUT(readsRecord.flags_.HalfDnbNoMatches(0))
                <<CGA_VOUT(readsRecord.flags_.HalfDnbNoMatches(1))
                );
        }

        for (size_t i=0; i<samMappings.size(); ++i)
        {
            mapping::SamRecord &m = samMappings[i];
            if (m.mates_.empty())
            {
                int bestMate = baseMappingRecords[i].bestMate_;
                if (bestMate>=0 && size_t(bestMate)!=i)
                {
                    m.mates_.push_back(&samMappings[bestMate]);
                    if (int(i)==baseMappingRecords[bestMate].bestMate_)
                        samMappings[bestMate].mates_.insert(samMappings[bestMate].mates_.begin(),&m);
                    else
                        samMappings[bestMate].mates_.push_back(&m);
                }
            }
        }
        return true;
    }

    namespace {
        class PrimaryMappingInfo {
        public:
            PrimaryMappingInfo(int maxWeight,size_t maxWeightInd,size_t bestMappingsCount)
                : maxWeight_(maxWeight), maxWeightInd_(maxWeightInd), bestMappingsCount_(bestMappingsCount)
            {}
            size_t  getBestMappingsCount() {return bestMappingsCount_;}
            bool    bestMappingsFound() const {return bestMappingsCount_>0;}

            int     maxWeight_;
            size_t  maxWeightInd_;
            size_t  bestMappingsCount_;
        };
    }

    size_t BaseMap2SamConverter::detectPrimaryMapping(
        const MappingsRecords &mappingsRecords, bool oneArmOnly, SamRecordArray& samMappings) const
    {
        //this happens if some mappings lie outside the chromosome boundaries
        //it's an error situation and it was most likely reported upstream -- ignore detection of primary
        if (mappingsRecords.size()!=samMappings.size())
            return 0;

        std::vector<PrimaryMappingInfo> 
            mapInfo(MAX_SIDES,PrimaryMappingInfo(0,mappingsRecords.size(),0)),
            singleMapInfo(MAX_SIDES,PrimaryMappingInfo(0,mappingsRecords.size(),0));

        for (size_t i=0; i<mappingsRecords.size(); ++i) 
        {
            const MappingsRecord &m = mappingsRecords[i];
            PrimaryMappingInfo &mi = (m.bestMate_==int(i))
                ? singleMapInfo[m.flags_.getSide()]
                : mapInfo[m.flags_.getSide()];
            if (m.weightChar_>mi.maxWeight_) {
                mi.maxWeight_ = m.weightChar_;
                mi.maxWeightInd_ = i;
                mi.bestMappingsCount_ = 1;
            } else if (m.weightChar_==mi.maxWeight_) {
                ++mi.bestMappingsCount_;
            }
        }
        using namespace boost::lambda;

        //int bestCount = std::accumulate(mapInfo.begin(), mapInfo.end(), 0,
        //    _1 += bind(&PrimaryMappingInfo::getBestMappingsCount, _2));

        if (oneArmOnly)
        {
            int bestSingleFound = std::accumulate(singleMapInfo.begin(), singleMapInfo.end(), 0, 
                _1 += bind(&PrimaryMappingInfo::bestMappingsFound, _2));
            CGA_ASSERT_EQ(bestSingleFound,1);

            BOOST_FOREACH(const PrimaryMappingInfo& i,singleMapInfo)
                if (i.bestMappingsFound())
                {
                    samMappings[i.maxWeightInd_].isPrimary_ = true;
                    return 1;
                }
            CGA_ERROR_EX("No single arm best record found in " << mappingsRecords.size() <<" records");
        } else {
            size_t bestFound = std::accumulate(mapInfo.begin(), mapInfo.end(), 0,
                _1 += bind(&PrimaryMappingInfo::bestMappingsFound, _2));
            if (bestFound==MAX_SIDES) 
            {
                int ind = mapInfo[0].maxWeightInd_;
                int bestMate = mappingsRecords[ind].bestMate_;
                while (bestMate != ind &&
                       mappingsRecords[bestMate].bestMate_ != bestMate &&
                       mappingsRecords[bestMate].bestMate_ != ind &&
                       mappingsRecords[bestMate].weightChar_ >= mappingsRecords[ind].weightChar_)
                {
                    ind = bestMate;
                    bestMate = mappingsRecords[ind].bestMate_;
                }
                const mapping::MappingsRecord &m = mappingsRecords[ind];
                samMappings[m.recordIndex_].isPrimary_ = true;
                samMappings[m.bestMate_].isPrimary_ = true;
                return 2;
            } else {
                size_t result = 0;

                BOOST_FOREACH(const PrimaryMappingInfo& i, singleMapInfo)
                {
                    if (i.bestMappingsFound()) 
                    {
                        samMappings[i.maxWeightInd_].isPrimary_ = true;
                        ++result;
                        break;
                    }
                }
                return result;
            }
        }
    }

} } // cgatools::mapping
