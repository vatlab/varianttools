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
#include "cgatools/mapping/AlleleAlignment.hpp"
#include "cgatools/mapping/EvidenceSamUtil.hpp"
#include "cgatools/mapping/GapsEstimator.hpp"
#include <math.h>

namespace cgatools { namespace mapping {

    void AlleleSequenceAligner::setInterval( 
        uint16_t chr, size_t begin, size_t length, const std::string* allele )
    {
        chr_ = chr;
        begin_ = begin;
        length_ = length;
        allele_ = allele;
        if (allele!=NULL)
            alleleLength_ = allele_->length();
        else
            alleleLength_ = 0;
        chrSequence_ = &reference_.listChromosomes()[chr_];
        CGA_ASSERT_LE(begin_+length_,chrSequence_->length());
    }

    void GapAndConcordanceExtractor::run( uint8_t side, bool strand, int32_t offsetInAllele,
        const std::string & alignmentCigar, const std::string & readSequence, const std::string & readScores)
    {
        gaps_.clear();
        concordance_ = 1;
        mismatches_ = 0;

        size_t startPos=0;
        int strandDirection = 1;

        if (side!=0)
            startPos+=halfDnbSize_;
        if (strand) {
            startPos+=halfDnbSize_-1;
            strandDirection=-1;
        }

        DnbSequenceIterator reads(readSequence,startPos,strandDirection,true);
        DnbSequenceIterator scores(readScores,startPos,strandDirection);

        Cigar c(alignmentCigar);
        endOffsetInAllele_ = offsetInAllele;
        char lastOp = '0';
        for (size_t i=0; i<c.getParsedCigar().size(); ++i)
        {
            const Cigar::CigarElement& e = c[i];
            if (e.type_==lastOp)
                gaps_.push_back(0);
            switch (e.type_)
            {
            case 'B' : 
                gaps_.push_back(-int8_t(e.length_));
                endOffsetInAllele_-=e.length_;
                break;
            case 'N' : 
                gaps_.push_back(e.length_);
                endOffsetInAllele_+=e.length_;
                break;
            case 'D':
            case 'I':
                CGA_ASSERT_MSG(false,"No insertions or deletions are expected in the input!");
                break;
            case 'M' : 
                for (int j=e.length_; j>0; --j)
                {
                    if (*reads==alleleAligner_[endOffsetInAllele_]) {
                        concordance_*=1.0-pow(double(10),double(*scores-33)/(-10));
                    } else {
                        concordance_*=pow(double(10),double(*scores-33)/(-10))/3;
                        ++mismatches_;
                    }

                    ++reads;
                    ++scores;
                    ++endOffsetInAllele_;
                }
                break;
            default:
                CGA_ASSERT_MSG(false,"The operation is not supported: " << e.type_ << " in " << c);
                break;
            }
            lastOp = e.type_;
        }
        if (strand)
            std::reverse(gaps_.begin(),gaps_.end());
    }


    void GapProbabilityAndConcordanceExtractor::run( 
        uint8_t side, bool strand, int32_t offsetInAllele, 
        const std::string & alignmentCigar, 
        const std::string & readSequence, 
        const std::string & readScores )
    {
        GapAndConcordanceExtractor::run(side,strand,offsetInAllele,alignmentCigar,readSequence,readScores);

        cgdata::SmallGapTuple gapTuple;
        CGA_ASSERT_EQ(gapTuple.size(),gaps_.size());
        for (size_t i=0; i<gapTuple.size();++i)
            gapTuple[i] = gaps_[i];

        int endPosition = (side==int8_t(strand)) ? offsetInAllele : endOffsetInAllele_;
        mapping::GapEst::TemplateSequenceRetriever<AlleleSequenceAligner> 
            retriever(alleleAligner_,alleleAligner_.getChrSequence().isCircular(),endPosition,side,strand);
        gapProbability_ = gapsEstimators_[side]->getProbability(gapTuple,retriever);
    }

} } // cgatools::mapping
