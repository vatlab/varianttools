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

#ifndef CGA_TOOLS_ALLELE_ALIGNMENT_HPP_
#define CGA_TOOLS_ALLELE_ALIGNMENT_HPP_ 1

//! @file AlleleAlignment.hpp

#include "cgatools/core.hpp"

#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"
#include "cgatools/util/BaseUtil.hpp"

namespace cgatools { namespace mapping {

//! The class is used to iterate through the bases of a variation allele from an evidenceInterval file
//! operator[] returns either a base from the allele or a reference base if the index points outside 
//! of the allele range
class AlleleSequenceAligner
{
public:
    AlleleSequenceAligner(const reference::CrrFile& reference)
        :reference_(reference)
        ,chrSequence_(NULL)
        ,allele_(NULL)
        ,alleleLength_(0)
    {}

    //! set current allele
    void setInterval(uint16_t chr, size_t begin, size_t length, const std::string* allele);

    //! returns the i-th base in the current allele 
    char operator[] (int i) const
    {
        if (i<0)
            return chrSequence_->getUnambiguousBase(begin_+i);
        if (i >= int(alleleLength_))
            return chrSequence_->getUnambiguousBase(begin_+length_-alleleLength_+i);
        return (*allele_)[i];
    }

    //! the position after the last available base
    int endSequencePosition() const {
        return chrSequence_->length()-length_+alleleLength_-begin_;
    }

    int startSequencePosition() const {
        return -int(begin_);
    }

    const reference::CompactDnaSequence& getChrSequence() const {return *chrSequence_;}
protected:
    const reference::CrrFile& reference_;
    const reference::CompactDnaSequence* chrSequence_;
    const std::string * allele_;
    size_t              alleleLength_;
    uint16_t            chr_;
    int64_t             begin_;
    size_t              length_;
};

//! An iterator through an aligned read sequence
//! The iterator is able to return bases in reverse order or/and complimented
class DnbSequenceIterator
{
public:
    //@param direction - 1 normal order, -1 reversed order
    //@param needsComplementing - defines if the sequence needs complementing when reversed
    DnbSequenceIterator(const std::string & sequence, 
        size_t startFrom, int direction=1, bool needsComplementing=false)
        :sequence_(sequence)
        ,position_(startFrom)
        ,direction_(direction)
        ,needsComplementing_(needsComplementing && direction_<0)
    {
        CGA_ASSERT_LE(0,position_);
        CGA_ASSERT_L(position_,sequence_.size());
    }

    DnbSequenceIterator& operator++ () 
    {
        position_+=direction_; 
        return *this;
    }

    //! return the current base
    char operator* () const {
        CGA_ASSERT_LE(0,position_);
        CGA_ASSERT_L(position_,sequence_.size());
        if (needsComplementing_)
            return util::baseutil::complement(sequence_[position_]);
        else
            return sequence_[position_];
    }

protected:
    std::string sequence_;
    size_t      position_;
    int         direction_;
    bool        needsComplementing_;
};

//! The class scans a given alignment for a given allele, provided by AlleleSequenceAligner 
//! in the constructor. While scanning it detects the gaps and computes the product of
//!  probabilities of bases (concordance), the number of mismatches and end position of the alignment
class GapAndConcordanceExtractor
{
public:
    GapAndConcordanceExtractor(const AlleleSequenceAligner& alleleAligner)
        :
        concordance_(-1)
        ,endOffsetInAllele_(0)
        ,mismatches_(0)
        ,alleleAligner_(alleleAligner)
        ,halfDnbSize_(35)
    {}

    void run( uint8_t side, bool strand, int32_t offsetInAllele, const std::string & alignmentCigar, 
                    const std::string & readSequence, const std::string & readScores);

    std::vector<int8_t> gaps_;
    double  concordance_;
    int32_t endOffsetInAllele_;
    size_t  mismatches_;

protected:
    const AlleleSequenceAligner& alleleAligner_;
    size_t  halfDnbSize_;
};

namespace GapEst {
    class GapsEstimator;
};

//! The class extends GapAndConcordanceExtractor by computing the probablility of sequence-dependent gaps
//! in addition to the metrics computed by the parent class
class GapProbabilityAndConcordanceExtractor : public GapAndConcordanceExtractor
{
public:

    //! takes an allele aligned and a sequence-dependant gap estimator
    GapProbabilityAndConcordanceExtractor(const AlleleSequenceAligner& alleleAligner,
        boost::array<boost::shared_ptr<GapEst::GapsEstimator>,2> gapsEstimators)
        : GapAndConcordanceExtractor(alleleAligner),gapsEstimators_(gapsEstimators)
    {}
    void run( uint8_t side, bool strand, int32_t offsetInAllele, const std::string & alignmentCigar, 
        const std::string & readSequence, const std::string & readScores);

    double  gapProbability_;
protected:
    boost::array<boost::shared_ptr<GapEst::GapsEstimator>,2> gapsEstimators_;
};

} } // cgatools::mapping

#endif // CGA_TOOLS_ALLELE_ALIGNMENT_HPP_
