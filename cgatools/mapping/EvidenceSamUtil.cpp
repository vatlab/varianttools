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
#include "cgatools/mapping/EvidenceSamUtil.hpp"
#include "cgatools/mapping/EvidenceCache.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/CrrFile.hpp"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

namespace cgatools { namespace mapping {

    SamSequenceSplitter::SamSequenceSplitter( const Cigar& extCigar, 
        const std::string& dnbSequence, 
        const std::string& dnbScore, 
        bool negativeStrand ) 
        :   srcExtCigar_(extCigar),
        srcSequence_(dnbSequence),
        srcScores_(dnbScore),
        tagCigar_(true,true)
    {
        if (negativeStrand) {
            std::reverse(srcSequence_.begin(),srcSequence_.end());
            BOOST_FOREACH(char &ch,srcSequence_)
                ch = util::baseutil::complement(ch);
            std::reverse(srcScores_.begin(),srcScores_.end());
        }
        convert();
    }

    void SamSequenceSplitter::convert()
    {
        size_t sequencePosition = 0;
        size_t skipSequenceLen = 0;
        BOOST_FOREACH(Cigar::CigarElement e, srcExtCigar_.getParsedCigar()) {
            switch (e.type_) {
                case '=':
                case 'X':
                case 'M':
                case 'I':
                    if (skipSequenceLen>0) {
                        if (e.length_ <= skipSequenceLen) {
                            skipSequenceLen -= e.length_;
                            sequencePosition += e.length_;
                            break;
                        } else {
                            tagCigar_.add(Cigar::CigarElement('S',e.length_-skipSequenceLen));
                            e.length_ -= skipSequenceLen;
                            sequencePosition += skipSequenceLen;
                            skipSequenceLen = 0;
                        }
                    } else {
                        tagCigar_.add(Cigar::CigarElement('S',e.length_));
                    }
                    sequence_<<srcSequence_.substr(sequencePosition,e.length_);
                    scores_<<srcScores_.substr(sequencePosition,e.length_);
                    sequencePosition+=e.length_;

                case 'N':
                    cigar_.add(e);
                    break;

                case 'P':
                case 'D':
                    if (skipSequenceLen==0)
                        cigar_.add(e);
                    break;

                case 'B':
                    CGA_ASSERT_MSG(!cigar_.getParsedCigar().empty(),"a cigar can't start from xxB");
                    sequence_.seekp(-int(e.length_),std::ios_base::cur);
                    scores_.seekp(-int(e.length_),std::ios_base::cur);
                    size_t startPos = sequencePosition-e.length_;
                    for (size_t i=0; i<e.length_; ++i) {
                        char scoreL = srcScores_[startPos+i];
                        char scoreR = srcScores_[sequencePosition+i];
                        if (scoreL>scoreR) {
                            sequence_.seekp(1,std::ios_base::cur);
                            scores_.seekp(1,std::ios_base::cur);
                        } else {
                            sequence_<<srcSequence_[sequencePosition+i];
                            scores_<<scoreR;
                        }
                    }
                    tagSequence_<<srcSequence_.substr(startPos,e.length_*2);
                    tagScores_<<srcScores_.substr(startPos,e.length_*2);
                    
                    CGA_ASSERT_MSG(tagCigar_.size()>0 && tagCigar_.back().type_=='S' 
                        && tagCigar_.back().length_>=e.length_,
                        "Wrong negative gap position: " << srcExtCigar_);
                    tagCigar_.back().length_-=e.length_;
                    tagCigar_.add(Cigar::CigarElement('G',e.length_));
                    skipSequenceLen = e.length_;
                    break;
            }
        }
    }


    void EvidenceDnbRecord::initParser( util::DelimitedFile &delimitedFile, 
        size_t formatVersion, const reference::CrrFile& crr)
    {
        cgdata::EvidenceReader::DnbRecord::bindToParser(delimitedFile,*this, crr);
    }

    void EvidenceDnbRecord::adjustOffset( int chrLength )
    {
        for (size_t i=0; i<2; ++i)
            if (offsetInReference_[i]<0) {
                offsetInReference_[i]+=chrLength;
                CGA_ASSERT_LE(0,offsetInReference_[i]);
            } else if (offsetInReference_[i]>=chrLength) {
                offsetInReference_[i]-=chrLength;
                CGA_ASSERT_L(offsetInReference_[i],chrLength);
            }
    }

    EvidenceSamRecord::EvidenceSamRecord( 
        const EvidenceDnbRecord& evidenceRecord, 
        bool isPrimary, 
        uint8_t mate, 
        size_t halfDnbSize, 
        const reference::CrrFile& reference
        ) 
    :   
        SamRecord(
        evidenceRecord.getId()
        ,true
        ,evidenceRecord.strand_
        ,isPrimary
        ,(evidenceRecord.side_==mate) ? 0:1
        ,evidenceRecord.chromosome_
        ,evidenceRecord.offsetInReference_[mate]
        ,evidenceRecord.referenceAlignment_[mate]
        ,(mate>0 && evidenceRecord.mateMappingQuality_>=0) 
            ? evidenceRecord.mateMappingQuality_-33:evidenceRecord.mappingQuality_-33
        ,true
        ,evidenceRecord.sequence_
        ,evidenceRecord.scores_
        ,std::make_pair((evidenceRecord.side_==mate) ? 0:halfDnbSize,halfDnbSize)
        )
        ,intervalId_(evidenceRecord.intervalId_)
        ,alleleIndex_(evidenceRecord.alleleIndex_)
        ,armWeight_(-1)
    {
        typeId_ = EVIDENCE;
    }

    EvidenceSamRecord::EvidenceSamRecord( 
        const EvidenceCacheDnbRecord& evidenceRecord, 
        const std::string& sequence,
        const std::string& scores,
        bool isPrimary, 
        uint8_t mate, 
        size_t halfDnbSize, 
        const reference::CrrFile& reference
        ) 
    :   
        SamRecord(
        evidenceRecord.getId()
        ,true
        ,evidenceRecord.strand_
        ,isPrimary
        ,(evidenceRecord.side_==mate) ? 0:1
        ,evidenceRecord.chromosome_
        ,evidenceRecord.offsetInReference_[mate]
        ,evidenceRecord.referenceAlignment_[mate]
        ,evidenceRecord.mappingQuality_-33
        ,true
        ,sequence
        ,scores
        ,std::make_pair((evidenceRecord.side_==mate) ? 0:halfDnbSize,halfDnbSize)
        )
        ,intervalId_(evidenceRecord.intervalId_)
        ,alleleIndex_(evidenceRecord.alleleIndex_)
        ,armWeight_(evidenceRecord.alleleConcordance_)
    {
        typeId_ = EVIDENCE_CACHE;
    }

    std::ostream& operator<< (std::ostream& out, const EvidenceDnbRecord& r)
    {
        out << cgdata::EvidenceReader::DnbRecord(r)
            << '\t' << r.mateMappingQuality_;
        return out;
    }
} } // cgatools::mapping
