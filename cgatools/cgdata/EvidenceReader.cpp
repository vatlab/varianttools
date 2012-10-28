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
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"

#include <boost/foreach.hpp>
#include <sstream>

namespace cgatools { namespace cgdata {

    using namespace cgatools::util;
    using namespace cgatools::reference;
    using namespace cgatools::variants;
    namespace bu = cgatools::util::baseutil;

    using namespace std;

    class EvidenceScoreField : public util::DelimitedFieldParser
    {
    public:
        EvidenceScoreField(
            const std::string& name,
            int32_t* evidenceScoreVAF,
            int32_t* evidenceScoreEAF)
            : DelimitedFieldParser(name),
              evidenceScoreVAF_(evidenceScoreVAF),
              evidenceScoreEAF_(evidenceScoreEAF)
        {
        }

        void parse(const char* first, const char* last)
        {
            *evidenceScoreVAF_ = parseValue<int32_t>(first, last);
            *evidenceScoreEAF_ = *evidenceScoreVAF_;
        }

    private:
        int32_t* evidenceScoreVAF_;
        int32_t* evidenceScoreEAF_;
    };

    class AlleleIndexesField : public DelimitedFieldParser
    {
    public:
        AlleleIndexesField(const std::string& name, std::vector<uint16_t>* alleleIndexes)
            : DelimitedFieldParser(name),
              alleleIndexes_(alleleIndexes)
        {
        }

        void parse(const char* first, const char* last)
        {
            alleleIndexes_->clear();
            while (first < last)
            {
                const char* next = first;
                while (next != last && *next != ';')
                    ++next;

                alleleIndexes_->push_back(parseValue<uint16_t>(first, next));
                first = next+1;
            }
        }

    private:
        std::vector<uint16_t>* alleleIndexes_;
    };

    EvidenceReader::IntervalRecord::IntervalRecord()
        :intervalId_(-1)
        ,chromosome_(std::numeric_limits<uint16_t>::max())
        ,offset_(std::numeric_limits<uint32_t>::max())
        ,length_(std::numeric_limits<uint32_t>::max())
        ,ploidy_(std::numeric_limits<uint16_t>::max())
        ,evidenceScoreVAF_(-1)
        ,evidenceScoreEAF_(-1)
    {
    }

    reference::Range EvidenceReader::IntervalRecord::getRange() const
    {
        return Range(chromosome_, offset_, offset_+length_);
    }

    namespace {
        //! Returns true if a given sequence is fully called and base-to-base
        //! compatible to reference sequence (which may contain ambiguous bases).
        bool calledReferenceCompatible(const std::string& seq, const std::string& refSeq)
        {
            if (seq.size() != refSeq.size())
                return false;
            for (size_t ii = 0; ii < seq.size(); ++ii)
            {
                if (!bu::isValidBase(seq[ii]))
                    return false;
                if (!bu::isConsistent(seq[ii], refSeq[ii]))
                    return false;
            }
            return true;
        }
    }

    bool EvidenceReader::IntervalRecord::isCompatible(
        uint16_t alleleIndex, const variants::Call& call, const reference::CrrFile& crr) const
    {
        string seq = call.calledSequence(crr);
        bool isRef = !seq.empty() &&
                     calledReferenceCompatible(seq, call.refSequence(crr));
        // For ref calls, skip I's. Assuming all other kinds of
        // calls are split by non-zero-length ref calls, for all other
        // calls we wish to explicitly include I's at the beginning and
        // end of the call of interest.
        return isCompatible(alleleIndex, call.range_, seq, isRef);
    }

    bool EvidenceReader::IntervalRecord::isCompatible(
        uint16_t alleleIndex,
        const reference::Range& range,
        const string& seq, bool ignoreAdjacentInsertions) const
    {
        if (range.chromosome_ != chromosome_)
            return false;

        if (range.begin_ < offset_ || range.end_ > offset_+length_)
            return false;

        // Skip sequence in the allele alignment until we get to the
        // call of interest.
        string alignment = getAlignment(alleleIndex);
        size_t alleleOffset = 0, refOffset = offset_, alignOffset = 0;
        while (refOffset < range.begin_)
        {
            CGA_ASSERT(alignOffset < alignment.size());
            switch (alignment[alignOffset])
            {
            case 'M':
                refOffset++;
                alleleOffset++;
                break;
            case 'I':
                alleleOffset++;
                break;
            case 'D':
                refOffset++;
                break;
            default:
                CGA_ASSERT(false);
                break;
            }
            alignOffset++;
        }

        if (ignoreAdjacentInsertions)
        {
            while (alignOffset < alignment.size() && 'I' == alignment[alignOffset])
            {
                alignOffset++;
                alleleOffset++;
            }
        }

        size_t origAlleleOffset = alleleOffset;
        while (alignOffset < alignment.size() &&
               ( refOffset < range.end_ ||
                 (refOffset == range.end_ && 'I' == alignment[alignOffset]) ))
        {
            if (ignoreAdjacentInsertions && refOffset == range.end_)
                break;

            switch (alignment[alignOffset])
            {
            case 'M':
                refOffset++;
                alleleOffset++;
                break;
            case 'I':
                alleleOffset++;
                break;
            case 'D':
                refOffset++;
                break;
            default:
                CGA_ASSERT(false);
                break;
            }
            alignOffset++;
        }

        string matchSeq = alleles_[alleleIndex].substr(origAlleleOffset, alleleOffset-origAlleleOffset);
        return seq == matchSeq;
    }

    std::string EvidenceReader::IntervalRecord::getAlignment(size_t alleleIndex) const
    {
        CGA_ASSERT(alleleIndex <= alleleAlignments_.size());

        if (0 == alleleIndex)
            return string(length_, 'M');

        return cigarToAlignment(alleleAlignments_[alleleIndex-1]);
    }

    std::string EvidenceReader::IntervalRecord::cigarToAlignment(const std::string& cigar)
    {
        string result;

        size_t offset = 0;
        while (offset < cigar.size())
        {
            size_t digitEnd = offset;
            while (digitEnd < cigar.size() && isdigit(cigar[digitEnd]))
                digitEnd++;

            if (offset == digitEnd || digitEnd == cigar.size())
                throw Exception("invalid cigar: "+cigar);

            size_t count = parseValue<size_t>(&cigar[offset], &cigar[digitEnd]);
            result.append(count, cigar[digitEnd]);

            offset = digitEnd+1;
        }

        return result;
    }

    void EvidenceReader::IntervalRecord::bindToParser( 
        util::DelimitedFile &parser, IntervalRecord& rec, const reference::CrrFile& crr )
    {
        parser.addField(ValueField<int32_t>("IntervalId", &rec.intervalId_));
        parser.addField(ChromosomeIdField("Chromosome", &rec.chromosome_, crr));
        parser.addField(ValueField<uint32_t>("OffsetInChromosome", &rec.offset_));
        parser.addField(ValueField<uint32_t>("Length", &rec.length_));
        parser.addField(ValueField<uint16_t>("Ploidy", &rec.ploidy_));
        if ( parser.hasField("Score") &&
             !(parser.hasField("EvidenceScoreVAF") || parser.hasField("EvidenceScoreEAF")) )
        {
            parser.addField(EvidenceScoreField("Score", &rec.evidenceScoreVAF_,
                                               &rec.evidenceScoreEAF_));
        }
        else
        {
            parser.addField(ValueField<int32_t>("EvidenceScoreVAF", &rec.evidenceScoreVAF_));
            parser.addField(ValueField<int32_t>("EvidenceScoreEAF", &rec.evidenceScoreEAF_));
        }
        parser.addField(AlleleIndexesField("AlleleIndexes", &rec.alleleIndexes_));
        parser.addField(StringField("Allele0", &rec.alleles_[0]));
        parser.addField(StringField("Allele1", &rec.alleles_[1]));
        parser.addField(StringField("Allele2", &rec.alleles_[2]));
        parser.addField(StringField("Allele3", &rec.alleles_[3]), DelimitedFile::FPT_OPTIONAL);
        parser.addField(StringField("Allele1Alignment", &rec.alleleAlignments_[0]));
        parser.addField(StringField("Allele2Alignment", &rec.alleleAlignments_[1]));
        parser.addField(StringField("Allele3Alignment", &rec.alleleAlignments_[2]),
                        DelimitedFile::FPT_OPTIONAL);
    }

    std::string EvidenceReader::DnbRecord::getId() const
    {
        ostringstream result;
        result << slide_ << "-" << lane_ << "-" << fileNumInLane_ << ":" << dnbOffsetInLaneFile_;
        return result.str();
    }

    bool EvidenceReader::DnbRecord::hasOverlap(const reference::Range& r, bool allowGap) const
    {
        std::string align = EvidenceReader::IntervalRecord::cigarToAlignment(referenceAlignment_[0]);

        int32_t referenceOffset = offsetInReference_[0];
        int32_t rangeStart = r.begin_, rangeEnd = r.end_;
        if (referenceOffset > rangeEnd)
        {
            return false;
        }
        else if (referenceOffset == rangeEnd)
        {
            char code = align[0];
            return code == 'I' || code == 'P';
        }
        uint32_t backtracks = 0;
        BOOST_FOREACH(char code, align)
        {
            // If the DNB spans the range of a length-changing variation,
            // this DNB is certainly interesting. If we walked out of the
            // variation interval of a length-preserving variation, then
            // by now we know that the overlap was only with DNB gap, and
            // the DNB is not interesting.
            if (referenceOffset == rangeEnd)
                return allowGap;
            if (referenceOffset >= rangeStart && (code != 'N' || allowGap))
                return true;
            if (code == 'B')
            {
                ++backtracks;
            }
            else if (backtracks != 0)
            {
                if (code != 'D')
                    --backtracks;
            }
            else
            {
                if (code != 'P' && code != 'I')
                    ++referenceOffset;
            }
        }
        return false;
    }

    void EvidenceReader::DnbRecord::bindToParser(
        util::DelimitedFile &parser, EvidenceReader::DnbRecord& rec, const reference::CrrFile& crr )
    {
        parser.addField(ValueField<int32_t>("IntervalId", &rec.intervalId_));
        parser.addField(ChromosomeIdField("Chromosome", &rec.chromosome_,crr));
        parser.addField(StringField("Slide", &rec.slide_));
        parser.addField(StringField("Lane", &rec.lane_));
        if (parser.hasField("FileNumInLane"))
            parser.addField(ValueField<uint16_t>("FileNumInLane", &rec.fileNumInLane_));
        else
            rec.fileNumInLane_ = 0; // old exports had one file per lane
        if (parser.hasField("DnbOffsetInLaneFile"))
            parser.addField(ValueField<uint32_t>("DnbOffsetInLaneFile", &rec.dnbOffsetInLaneFile_));
        else
            parser.addField(ValueField<uint32_t>("OffsetInLane", &rec.dnbOffsetInLaneFile_));
        parser.addField(util::ValueField<uint8_t>("AlleleIndex", &rec.alleleIndex_));

        parser.addField(util::SideField("Side",&rec.side_));
        parser.addField(util::StrandField("Strand",&(rec.strand_)));
        parser.addField(util::ValueField<int32_t>("OffsetInAllele",&rec.offsetInAllele_));
        parser.addField(util::StringField("AlleleAlignment",&rec.alleleAlignment_));

        parser.addField(ValueField<int32_t>("OffsetInReference", &rec.offsetInReference_[0]));
        parser.addField(StringField("ReferenceAlignment", &rec.referenceAlignment_[0]));
        parser.addField(ValueField<int32_t>("MateOffsetInReference", &rec.offsetInReference_[1]));
        parser.addField(StringField("MateReferenceAlignment", &rec.referenceAlignment_[1]));
        parser.addField(ValueField<int32_t>("ScoreAllele0", &rec.scoreAllele_[0]));
        parser.addField(ValueField<int32_t>("ScoreAllele1", &rec.scoreAllele_[1]));
        parser.addField(ValueField<int32_t>("ScoreAllele2", &rec.scoreAllele_[2]).
                        exception("", std::numeric_limits<int32_t>::min()));
        parser.addField(ValueField<int32_t>("ScoreAllele3", &rec.scoreAllele_[3]).
                        exception("", std::numeric_limits<int32_t>::min()), DelimitedFile::FPT_OPTIONAL);
        parser.addField(util::CharField("MappingQuality",(char *)&rec.mappingQuality_));
        parser.addField(util::StringField("Sequence",&rec.sequence_));
        parser.addField(util::StringField("Scores",&rec.scores_));
    }

    EvidenceReader::DnbRecord::DnbRecord()
        :intervalId_(-1)
        ,fileNumInLane_(std::numeric_limits<uint16_t>::max())
        ,dnbOffsetInLaneFile_(std::numeric_limits<uint32_t>::max())
        ,alleleIndex_(std::numeric_limits<uint8_t>::max())
        ,side_(std::numeric_limits<uint8_t>::max())
        ,strand_(true)
        ,offsetInAllele_(std::numeric_limits<int32_t>::min())
        ,mappingQuality_(std::numeric_limits<uint8_t>::max())
    {
        offsetInReference_.assign(-1);
        scoreAllele_.assign(-1);
    }


    EvidenceReader::IntervalsFile::IntervalsFile(const std::string& fn, const CrrFile& crr)
        :   filename_(fn),
            f_(util::InputStream::openCompressedInputStreamByExtension(fn)),
            parser_(*f_, fn)
    {
        IntervalRecord::bindToParser(parser_,rec_,crr);
    }

    EvidenceReader::DnbsFile::DnbsFile(const std::string& fn, const reference::CrrFile& crr)
        : filename_(fn),
          f_(util::InputStream::openCompressedInputStreamByExtension(fn)),
          parser_(*f_, fn)
    {
        DnbRecord::bindToParser(parser_,rec_, crr);
    }

    EvidenceReader::EvidenceReader(
        const reference::CrrFile& crr, const GenomeMetadata& exp)
        :   crr_(crr), exp_(exp), dnbsChromosome_(std::numeric_limits<uint16_t>::max()), inInterval_(false)
    {
    }

    void EvidenceReader::openIntervals(uint16_t chr)
    {
        if (0 == intervalsFile_.get() || chr != intervalsFile_->rec_.chromosome_)
        {
            std::string fn =
                exp_.getEvidenceIntervalsFileName(
                crr_.listChromosomes()[chr].getName());
            intervalsFile_.reset(new IntervalsFile(fn, crr_));
            intervalsFile_->rec_.chromosome_ = chr;
            nextInChr();
        }
    }

    void EvidenceReader::seek(const reference::Range& r)
    {
        openIntervals(r.chromosome_);

        while (intervalsFile_->rec_.offset_ + intervalsFile_->rec_.length_ < r.begin_)
            nextInChr();
        inInterval_ = (intervalsFile_->rec_.offset_ + intervalsFile_->rec_.length_ >= r.end_ &&
                       intervalsFile_->rec_.offset_ <= r.begin_);
    }

    void EvidenceReader::seekToChr( uint16_t chr )
    {
        intervalsFile_.reset(NULL);
        openIntervals(chr);
    }

    void EvidenceReader::nextInChr()
    {
        if (!intervalsFile_->parser_.next())
        {
            intervalsFile_->rec_.offset_ = std::numeric_limits<uint32_t>::max();
            intervalsFile_->rec_.length_ = 0;
            inInterval_ = false;
        } else 
        {
            inInterval_ = true;
        }
    }

    uint32_t EvidenceReader::countSupportingDnbs(uint16_t alleleIndex, int32_t scoreThreshold)
    {
        CGA_ASSERT(inInterval_);
        if (intervalsFile_->rec_.alleleIndexes_.size() != 2)
            return 0;
        if (intervalsFile_->rec_.alleleIndexes_[0] == intervalsFile_->rec_.alleleIndexes_[1])
            return 0;
        if (intervalsFile_->rec_.alleleIndexes_[0] != alleleIndex &&
            intervalsFile_->rec_.alleleIndexes_[1] != alleleIndex)
            return 0;

        seekDnbs();
        uint16_t altIndex = intervalsFile_->rec_.alleleIndexes_[0];
        if (altIndex == alleleIndex)
            altIndex = intervalsFile_->rec_.alleleIndexes_[1];

        uint32_t result = 0;
        BOOST_FOREACH(const DnbRecord& rec, dnbs_)
        {
            if (rec.scoreAllele_[alleleIndex] - rec.scoreAllele_[altIndex] >= scoreThreshold)
                result++;
        }
        return result;
    }

    void EvidenceReader::seekDnbs()
    {
        openDnbs();

        if (dnbsFile_->rec_.intervalId_ <= intervalsFile_->rec_.intervalId_)
            dnbs_.clear();
        while (dnbsFile_->rec_.intervalId_ <= intervalsFile_->rec_.intervalId_) 
        {
            if (dnbsFile_->rec_.intervalId_ == intervalsFile_->rec_.intervalId_)
            {
                dnbs_.push_back(dnbsFile_->rec_);
            }
            if (!nextDnbs())
                break;
        } 
        
    }

    void EvidenceReader::openDnbs()
    {
        CGA_ASSERT(inInterval());
        if (dnbsChromosome_ != intervalsFile_->rec_.chromosome_)
        {
            std::string fn =
                exp_.getEvidenceDnbsFileName(
                crr_.listChromosomes()[intervalsFile_->rec_.chromosome_].getName());
            dnbsFile_.reset(new DnbsFile(fn, crr_));
            dnbsChromosome_ = intervalsFile_->rec_.chromosome_;
            nextDnbs();
        }
    }

    bool EvidenceReader::nextDnbs()
    {
        if (!dnbsFile_->parser_.next())
        {
            dnbsFile_->rec_.intervalId_ = std::numeric_limits<uint32_t>::max();
            dnbsFile_->rec_.scoreAllele_.assign(std::numeric_limits<int32_t>::min());
            return false;
        }
        return true;
    }

    std::ostream& operator <<(std::ostream& out, const EvidenceReader::DnbRecord& r) 
    {
        out << r.intervalId_
            << '\t' << r.chromosome_
            << '\t' << r.slide_
            << '\t' << r.lane_
            << '\t' << r.fileNumInLane_
            << '\t' << r.dnbOffsetInLaneFile_
            << '\t' << r.alleleIndex_
            << '\t' << r.side_
            << '\t' << r.strand_
            << '\t' << r.offsetInAllele_
            << '\t' << r.alleleAlignment_
            << '\t' << r.offsetInReference_[0]
            << '\t' << r.referenceAlignment_[0]
            << '\t' << r.offsetInReference_[1]
            << '\t' << r.referenceAlignment_[1]
            << '\t' << r.mappingQuality_
            << '\t' << r.scoreAllele_[0]
            << '\t' << r.scoreAllele_[1]
            << '\t' << r.scoreAllele_[2]
            << '\t' << r.sequence_
            << '\t' << r.scores_
            ;
        return out;
    }

} } // cgatools::cgdata
