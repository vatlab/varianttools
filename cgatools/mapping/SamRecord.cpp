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
#include "cgatools/mapping/SamRecord.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/mapping/EvidenceSamUtil.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "Cigar.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <cmath>

namespace cgatools { namespace mapping {

    namespace {
        class SamRecordSequenceSplitter : public mapping::SamSequenceSplitter
        {
        public:
            SamRecordSequenceSplitter(const SamRecord& record)
                : mapping::SamSequenceSplitter(record.extCigar_,
                record.fullReadSequence_.substr(record.sequenceStartAndLength_[0].first,
                record.sequenceStartAndLength_[0].second),
                record.fullReadScores_.substr(record.sequenceStartAndLength_[0].first,
                record.sequenceStartAndLength_[0].second),
                record.onNegativeStrand_)
            {}
        };

    }

    SamRecord::SamRecord( 
        const std::string& readName, 
        bool isMapped, 
        bool onNegativeStrand, 
        bool isPrimary, 
        uint8_t side, 
        uint16_t chr, 
        int32_t position, 
        const std::string& extCigar, 
        uint8_t mappingQuality, 
        bool    isConsistentMapQ, 
        const std::string& fullReadSequence, 
        const std::string& fullReadScores, 
        UInt16Pair sequenceStartAndLength) 
    :
        typeId_(DEFAULT)
        ,readName_(readName)
        ,isMapped_(isMapped)
        ,onNegativeStrand_(onNegativeStrand)
        ,isPrimary_(isPrimary)
        ,isGroupPrimary_(false)
        ,side_(side)
        ,chr_(chr)
        ,position_(position)
        ,extCigar_(extCigar)
        ,fullReadSequence_(fullReadSequence)
        ,fullReadScores_(fullReadScores)
        ,isSvCandidate_(false)
        ,mappingQuality_(mappingQuality)
        ,isConsistentMapQ_(isConsistentMapQ)
        ,isArtificialMateReported_(false)
    {
        sequenceStartAndLength_[0] = sequenceStartAndLength;

        if (sequenceStartAndLength.first==0)
            sequenceStartAndLength_[1] = UInt16Pair(sequenceStartAndLength.second,
                                            fullReadSequence.length()-sequenceStartAndLength.second);
        else
            sequenceStartAndLength_[1] = UInt16Pair(0,sequenceStartAndLength.first);
    }

    bool SamRecord::correctPosition(const reference::CrrFile& reference)
    {
        if ((!isMapped_) && position_<0)
            position_ = 0;
        size_t chrLength = reference.listChromosomes()[chr_].length();
        if (position_<0) {
            if(!reference.listChromosomes()[chr_].isCircular())
                return false;
            position_+=chrLength;
            if (position_<0)
                return false;
        } else if (position_>=int(chrLength)) {
            if(!reference.listChromosomes()[chr_].isCircular())
                return false;
            position_-=chrLength;
            if(position_>int(chrLength))
                return false;
        }
        return true;
    }

    bool SamRecordGenerator::isConsistent(const SamRecord &r) const
    {
        CGA_ASSERT(!r.mates_.empty());

        if ((!r.isMapped_) ||(!r.mates_[0]->isMapped_))
            return false;
        const SamRecord &mate = *r.mates_[0];
        if (r.chr_!=mate.chr_ || r.onNegativeStrand_ != mate.onNegativeStrand_)
            return false;

        if (std::abs(int(r.position_ - mate.position_)) > samGeneratorConfig_.maxConsistentRange_)
        {
            if(reference_.listChromosomes()[r.chr_].isCircular())
            {
                return int(reference_.listChromosomes()[r.chr_].length()) -
                        std::abs(r.position_-mate.position_) < samGeneratorConfig_.maxConsistentRange_;
            } else
               return false;
        }
        return true;
    }

    void SamRecordGenerator::OutputFileDescriptor::writeHeader(const SamFileHeaderBlock& header)
    {
        CGA_ASSERT(!hasHeader_);
        outStream_ << header << std::endl;
        hasHeader_ = true;
    }

    SamRecordGenerator::OutputFileDescriptor::~OutputFileDescriptor()
    {
        if (hasHeader_ && !hasRecords_)
            outStream_ << "empty\t0\t*\t0\t0\t0M\t*\t0\t0\t\t" << std::endl;
    }

    void SamRecordGenerator::printMateSequence(OutputFileDescriptor& out, 
        const std::string &mateSeq, const std::string &mateScore )
    {
        out.outStream_  << SAM_SEPARATOR << "R2:Z:" << mateSeq;         // TAG R2:Z:Mate sequence
        out.outStream_  << SAM_SEPARATOR << "Q2:Z:" << mateScore;       // TAG Q2:Z:Mate scores
    }

    void SamRecordGenerator::flagAsSVCandidate(OutputFileDescriptor& out)
    {
        out.outStream_  << SAM_SEPARATOR << "XS:I:1";                   // TAG XS:i:SVCandidate
    }

    void SamRecordGenerator::printNegativeGapTag(OutputFileDescriptor& out, 
        const mapping::SamSequenceSplitter &splitter)
    {
        std::string tagReadsStr(splitter.tagSequence_.str());
        if (!tagReadsStr.empty()) {
            out.outStream_ << SAM_SEPARATOR << "GC:Z:" << splitter.tagCigar_;       // TAG GC:Z:NEGATIVE GAPS
            out.outStream_ << SAM_SEPARATOR << "GS:Z:" << tagReadsStr;              // TAG GS:Z:NEGATIVE GAPS
            out.outStream_ << SAM_SEPARATOR << "GQ:Z:" << splitter.tagScores_.str();// TAG GQ:Z:NEGATIVE GAPS
        }
    }

    void SamRecordGenerator::printAlleleInfoTag(OutputFileDescriptor& out, 
        const mapping::EvidenceSamRecord& evidenceRecord)
    {
        out.outStream_ << SAM_SEPARATOR << "ZI:I:" << evidenceRecord.intervalId_;
        out.outStream_ << SAM_SEPARATOR << "ZA:I:" << evidenceRecord.alleleIndex_;
    }

    void SamRecordGenerator::printReadGroup(OutputFileDescriptor& out)
    {
        out.outStream_ << SAM_SEPARATOR << "RG:Z:" << readGroup_;
    }

    void SamRecordGenerator::printAlternatives(
        OutputFileDescriptor& out,
        const std::string& tag, 
        const mapping::SamRecord::SamRecords& samRecords,
        size_t startFrom
    )
    {
        if (samRecords.size() <= startFrom)
            return;

        out.outStream_ << SAM_SEPARATOR << tag << ":Z:";
        for (size_t i=startFrom; i<samRecords.size(); ++i)
        {
            const mapping::SamRecord& m = *samRecords[i];
            if (!m.isMapped_)
                continue;

            out.outStream_ << getSamChr(m.chr_) << ",";
            out.outStream_ << getAdjustedSamPosition(m) << ",";

            SamRecordSequenceSplitter s(m);

            out.outStream_ << s.cigar_ << ",";
            out.outStream_ << int(m.getMappingQuality()) << ";";
        }
    }

    std::string SamRecordGenerator::getSamChr(uint16_t chr) const
    {
        return reference_.listChromosomes()[chr].getName();
    }

    SamRecordGenerator::OutputFileDescriptor& SamRecordGenerator::getOutputStream(const std::string &id)
    {
        OutStreamMap::iterator it = outFiles_.find(id);
        if (it==outFiles_.end())
            it = outFiles_.find("all");
        CGA_ASSERT_MSG(it!=outFiles_.end(),
            "Unknown stream name: " << id << ". was not able to substitute by 'all'");

        OutputFileDescriptor& out = *it->second;
        CGA_ASSERT(out.hasHeader_);

        return out;
    }

    void SamRecordGenerator::mappingRecordToSam(const SamRecord& record)
    {
        // apply primary mappings only filter
        if (samGeneratorConfig_.primaryMappingsOnly_ && !record.isPrimary_)
            return;

        CGA_ASSERT(!record.mates_.empty());
        const SamRecord* mate = record.mates_[0];

        CGA_ASSERT_MSG( (record.isPrimary_ && mate->isPrimary_) || !record.isPrimary_, 
            "Error in primary flags in " << record << " Primary mappings have to have primary mates");

        std::string chr = (record.isMapped_ ? getSamChr(record.chr_) 
                                            : (mate->isMapped_ ? getSamChr(mate->chr_) : "*"));

        CGA_ASSERT_L(0,outFiles_.size());
        OutputFileDescriptor& out = getOutputStream(chr);
        out.hasRecords_ = true;

        mapping::SamRecordSequenceSplitter splitter(record);

        bool consistent = record.isConsistent();

        size_t flag = 0x0001;                       //FLAG: the read is paired in sequencing

        if (consistent)
            flag |= 0x0002;                         //FLAG: both mates are mapped consistently

        if (!record.isMapped_)
            flag |= 0x0004;                         //FLAG: the query sequence itself is unmapped
        if (!mate->isMapped_)
            flag |= 0x0008;                         //FLAG: the mate is unmapped

        if (mate->onNegativeStrand_)
        {
            CGA_ASSERT_MSG(mate->isMapped_, "unmapped on a negative strand");
            flag |= 0x0020;                         //FLAG: the mate is on the reverse strand
        }

        if (record.onNegativeStrand_)
        {
            CGA_ASSERT_MSG(record.isMapped_, "unmapped on a negative strand");
            flag = flag | 0x0010;                   //FLAG: on the reverse strand
        }

        if (record.side_ == 0)
            flag |= 0x0040;                         //FLAG: the read is the first read in a pair
        else
            flag |= 0x0080;                         //FLAG: the read is the second read in a pair

        if (!record.isPrimary_)
            flag = flag | 0x0100;                   //FLAG: the alignment is not primary


        out.outStream_ << record.readName_;                                      // QNAME
        out.outStream_ << SAM_SEPARATOR << flag;                                 // FLAG

        out.outStream_ << SAM_SEPARATOR << chr;                                  // RNAME

        out.outStream_ << SAM_SEPARATOR << (record.isMapped_ ?
            getAdjustedSamPosition(record) : getAdjustedSamPosition(*mate));     // POS (1-based)

        out.outStream_ << SAM_SEPARATOR << (record.isMapped_ ?
            int(record.getMappingQuality()) : 0);                                // MAPQ

        if (record.isMapped_)                                                    // CIGAR
        {
            Cigar outCigar = samGeneratorConfig_.packCigar_ 
                ? splitter.cigar_.pack() 
                : splitter.cigar_;

            if (samGeneratorConfig_.removePaddingAtCigarEnds_)
                outCigar.trancatePaddings();

            out.outStream_ 
                << SAM_SEPARATOR 
                << (outCigar);
        } else {
            out.outStream_ << SAM_SEPARATOR << "*";  
        }

        if (record.isMapped_ 
            && record.isArtificialMateReported()
            && (!mate->isMapped_) 
            && mate->mates_[0]==&record 
            )
        {
            mate = &record; //print record coordinates for unmapped primary mate
        }

        bool printMatePosition = (mate->isMapped_ && !record.isMapped_)
            || consistent 
            || record.isSvCandidate_
            || mate == &record
            ;

        out.outStream_ << SAM_SEPARATOR << ((!printMatePosition)
            ? "*" :
              (mate->chr_==record.chr_ || record.chr_==0) ?
                                                      "=" :
                                                      getSamChr(mate->chr_));     // MRNM

        out.outStream_ << SAM_SEPARATOR << 
            (printMatePosition ? getAdjustedSamPosition(*mate) : 0);              // MPOS (1-based)

        out.outStream_ << SAM_SEPARATOR << 
            ((mate->chr_!=record.chr_ || (!consistent) || (!record.isMapped_))
            ? 0 : int(getAdjustedSamPosition(*mate))
                    -int(getAdjustedSamPosition(record)));                        // ISIZE

        mate = NULL; //the mate could be overridden above to the record in a case it wasn't mapped

        out.outStream_ << SAM_SEPARATOR << splitter.sequence_.str();              // SEQ
        out.outStream_ << SAM_SEPARATOR << splitter.scores_.str();                // QUAL

        if (!readGroup_.empty())
            printReadGroup(out);

        printNegativeGapTag(out,splitter);

        CGA_ASSERT_MSG(record.mates_[0]->isMapped_ || record.mates_.size()==1, 
            CGA_VOUT(record.mates_[0]->isMapped_)<<CGA_VOUT(record.mates_.size()));

        if (samGeneratorConfig_.addMateSequenceAndScore_
            || (samGeneratorConfig_.addMateSequenceAndScoreIfNotLocallyAvailable_ 
                && (!consistent)
                && (samGeneratorConfig_.skipNotMapped_ || !record.isPrimary_)
               )
            ) 
        {
            printMateSequence(
                out,
                record.fullReadSequence_.substr(record.sequenceStartAndLength_[1].first,
                record.sequenceStartAndLength_[1].second),
                record.fullReadScores_.substr(record.sequenceStartAndLength_[1].first,
                record.sequenceStartAndLength_[1].second)
                );
        }

        if (record.isSvCandidate_)
            flagAsSVCandidate(out);

        if (record.typeId_==SamRecord::EVIDENCE || record.typeId_==SamRecord::EVIDENCE_CACHE)
        {
            const EvidenceSamRecord& evidenceRecord = static_cast<const EvidenceSamRecord &>(record);
            if (samGeneratorConfig_.printAlleleInfo_)
                printAlleleInfoTag(out,evidenceRecord);
        }

        if (samGeneratorConfig_.addAlternativeMappings_)
            printAlternatives(out,"XA",record.alternatives_,0);

        if (samGeneratorConfig_.addAlternativeMates_)
            printAlternatives(out,"ZM",record.mates_,1);

        out.outStream_ << std::endl;
    }

    uint32_t SamRecordGenerator::getAdjustedSamPosition(const SamRecord& record) const
    {
        if (!record.isMapped_)
        {
            if (record.position_<0)
                return 0;
            else
                return record.position_+1;
        }
        int result = record.position_;
        size_t chrLength = reference_.listChromosomes()[record.chr_].length();
        if (result < 0) {
            result+=chrLength;
            CGA_ASSERT_MSG(reference_.listChromosomes()[record.chr_].isCircular(),
                reference_.listChromosomes()[record.chr_].getName());
            CGA_ASSERT_LE(0,result);
        } else if (record.position_ >= int(chrLength)) {
            result-=chrLength;
            CGA_ASSERT_MSG(reference_.listChromosomes()[record.chr_].isCircular(),
                reference_.listChromosomes()[record.chr_].getName());
            CGA_ASSERT_L(result,int(chrLength));
        }
        return result+1;
    }

    void SamRecordGenerator::setHeader(const SamFileHeaderBlock& header)
    {
        header_ = header;
        readGroup_ = header.getChild("@RG").getChild("ID").value_;
        BOOST_FOREACH(OutStreamMap::value_type &os, outFiles_)
        {
            SamFileHeaderBlock h = header;
            
            //remove chromosomes that are not present in the SAM records
            if ((!samGeneratorConfig_.mateSvCandidates_) && (os.first!="all"))
            {
                h.removeNotMatchingChildren(SamFileHeaderBlock(os.first, "@SQ", "", ""));
            }

            os.second->writeHeader(h);
        }
    }

    SamRecordGenerator::SamRecordGenerator( 
        std::ostream& outSamFile, 
        const reference::CrrFile& reference, 
        const SamGeneratorConfig& config, 
        const std::vector<std::string>& outStreams 
        ) 
    :   reference_(reference), samGeneratorConfig_(config)
    {
        bool allStreamIsDefined = false;
        if (!outStreams.empty())
        {
            BOOST_FOREACH(const std::string s, outStreams)
            {
                std::vector<std::string> parsedStr;
                boost::split(parsedStr,s,boost::is_any_of(","));
                if (parsedStr.size()!=2)
                    CGA_ERROR_EX("The output stream format should be 'chrName,filePath', got instead:" << s);
                const std::string& chr = parsedStr[0];
                const std::string& path = parsedStr[1];
                if (chr=="all")
                    allStreamIsDefined = true;

                OpenStreamPtr ostr = util::OutputStream::openCompressedOutputStreamByExtension(path);
                openStreams_.push_back(ostr);

                outFiles_.insert(std::make_pair(chr, new OutputFileDescriptor(*ostr)));
            }
        }

        if (!allStreamIsDefined)
            outFiles_.insert(std::make_pair("all", new OutputFileDescriptor(outSamFile)));
    }

    SamRecordGenerator::~SamRecordGenerator()
    {
        outFiles_.clear();
    }

    std::ostream& operator<< (std::ostream& ost, const SamRecord& r)
    {
        char sep='\t';
        ost         << r.typeId_
            << sep  << r.readName_
            << sep  << (r.isMapped_?'M':'m')
            << sep  << (r.onNegativeStrand_?'-':'+')
            << sep  << (r.isPrimary_?'P':'S')
            << sep  << int(r.side_)
            << sep  << r.chr_
            << sep  << r.position_
            << sep  << r.extCigar_
            << sep  << int(r.mappingQuality_) << '=' << std::pow(double(10),double(r.mappingQuality_)/(-10))
            << sep  << r.sequenceStartAndLength_[0].first << ';' << r.sequenceStartAndLength_[0].second
            << sep  << r.sequenceStartAndLength_[1].first << ';' << r.sequenceStartAndLength_[1].second

            ;

        ost << sep << "m:";
        for(SamRecord::SamRecords::const_iterator it=r.mates_.begin();it!=r.mates_.end();++it)
        {
            if (it!=r.mates_.begin())
                ost << ';';
            ost << int((**it).chr_) << ((**it).onNegativeStrand_?'-':'+') << (**it).position_;
        }
        return ost;
    }

    std::ostream& operator<< (std::ostream& ostr,const SamFileHeaderBlock& block) 
    {
        ostr << block.type_;
        if (block.children_.empty()) 
        {
            ostr << ':' << block.value_;
        } else 
        {
            BOOST_FOREACH(const SamFileHeaderBlock& b, block.children_)
                ostr << b.separator_ << b;
        }
        return ostr;
    }


    SamFileHeaderBlock& SamFileHeaderBlock::get(const SamFileHeaderBlock& b)
    {
        Children::iterator it = std::find(children_.begin(),children_.end(),b);
        if (it==children_.end())
        {
            children_.push_back(b);
            it = children_.end()-1;
        }
        return *it;
    }

    SamFileHeaderBlock& SamFileHeaderBlock::add(const SamFileHeaderBlock& b)
    {
        Children::iterator it = std::find(children_.begin(), children_.end() ,b);
        if (it == children_.end())
            children_.push_back(b);
        else
            CGA_ERROR_EX("The child with the given id already exists:'" << b << "'.");
        return *this;
    }

    const SamFileHeaderBlock& SamFileHeaderBlock::getChild(const std::string &child) const
    {
        Children::const_iterator it = std::find(children_.begin(), children_.end(), child);
        if (it != children_.end())
            return *it;
        else
            CGA_ERROR_EX("Invalid child id:'" << child << "'.");
    }

    namespace {
        
        class IsNotMatching {
        public:
            IsNotMatching(const SamFileHeaderBlock &b) : b_(b) {}

            bool operator() (const SamFileHeaderBlock &other) {
                return other.type_ == b_.type_ && other.id_ != b_.id_;
            }
            const SamFileHeaderBlock &b_;
        };

    }

    void SamFileHeaderBlock::removeNotMatchingChildren(const SamFileHeaderBlock &child)
    {
        Children::iterator it = std::remove_if(children_.begin(), children_.end(), IsNotMatching(child));
        children_.erase(it, children_.end());
    }

} } // cgatools::mapping
