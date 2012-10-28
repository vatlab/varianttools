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

#ifndef CGATOOLS_CGDATA_EVIDENCEREADER_HPP_
#define CGATOOLS_CGDATA_EVIDENCEREADER_HPP_ 1

//! @file EvidenceReader.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/variants/Call.hpp"

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace cgdata {

    class EvidenceReader : boost::noncopyable
    {
    public:
        struct IntervalRecord
        {
            IntervalRecord();
            static void bindToParser(
                util::DelimitedFile &parser, IntervalRecord& rec, const reference::CrrFile& crr);

            int32_t intervalId_;
            uint16_t chromosome_;
            uint32_t offset_;
            uint32_t length_;
            uint16_t ploidy_;
            int32_t evidenceScoreVAF_;
            int32_t evidenceScoreEAF_;
            std::vector<uint16_t> alleleIndexes_;
            boost::array<std::string, 4> alleles_;
            boost::array<std::string, 3> alleleAlignments_;

            reference::Range getRange() const;

            //! Returns true if the given call is compatible to the
            //! evidence allele at the specified index.
            bool isCompatible(
                uint16_t alleleIndex, const variants::Call& call,
                const reference::CrrFile& crr) const;

            //! Returns true if the sequence is compatible to the
            //! evidence allele at the specified index.
            //! ignoreAdjacentInsertions controls whether the comparison
            //! to the allele includes or skips the inserted sequence at
            //! the boundaries of the given range.
            bool isCompatible(
                uint16_t alleleIndex,
                const reference::Range& range,
                const std::string& seq,
                bool ignoreAdjacentInsertions) const;

            std::string getAlignment(size_t alleleIndex) const;
            static std::string cigarToAlignment(const std::string& cigar);
        };

        struct DnbRecord
        {
            friend std::ostream& operator<< (std::ostream& out, const DnbRecord& r);

            DnbRecord();
            std::string getId() const;

            //! Returns true of the DNB overlaps the specified range
            //! (overlap by gap is counted if allowGap is true)
            bool hasOverlap(const reference::Range& r, bool allowGap) const;

            static void bindToParser(util::DelimitedFile &parser, 
                                        DnbRecord& rec, const reference::CrrFile& crr);

            int32_t     intervalId_;
            uint16_t chromosome_;
            std::string slide_;
            std::string lane_;
            uint16_t    fileNumInLane_;
            uint32_t    dnbOffsetInLaneFile_;
            uint8_t     alleleIndex_;
            uint8_t     side_;
            bool        strand_;
            int32_t                     offsetInAllele_;
            std::string                 alleleAlignment_;
            boost::array<int32_t,2>     offsetInReference_;
            boost::array<std::string,2> referenceAlignment_;

            boost::array<int32_t,4>     scoreAllele_;
            uint8_t                     mappingQuality_;

            std::string sequence_;
            std::string scores_;
        };

        EvidenceReader(const reference::CrrFile& crr,
                       const GenomeMetadata& exp);

        void seek(const reference::Range& r);

        //! moves the cursor to the first interval of the given chromosome
        //! sets inInterval() into false if the cursor moves behind the end of the chromosome
        void seekToChr(uint16_t chr);

        //! moves the cursor to the next interval in the current chromosome
        //! sets inInterval() into false if the cursor moves behind the end of the chromosome
        void nextInChr();

        bool inInterval() const
        {
            return inInterval_;
        }

        const IntervalRecord& getInterval() const
        {
            CGA_ASSERT(inInterval_);
            return intervalsFile_->rec_;
        }

        const std::vector<DnbRecord>& getDnbs()
        {
            seekDnbs();
            return dnbs_;
        }

        //! Returns the evidenceScoreVAF of the interval at the current
        //! position or zero if the last seek didn't find an overlapping
        //! interval.
        int32_t getEvidenceScoreVAF() const
        {
            if (!inInterval_)
                return 0;
            else
                return intervalsFile_->rec_.evidenceScoreVAF_;
        }

        //! Returns the evidenceScoreEAF of the interval at the current
        //! position or zero if the last seek didn't find an overlapping
        //! interval.
        int32_t getEvidenceScoreEAF() const
        {
            if (!inInterval_)
                return 0;
            else
                return intervalsFile_->rec_.evidenceScoreEAF_;
        }

        //! Returns the number of DNBs supporting the given allele of
        //! the top hypothesis, such that the
        //! ScoreAllele[alleleIndex]-ScoreAllele[alternate] >
        //! scoreThreshold.
        uint32_t countSupportingDnbs(uint16_t alleleIndex, int32_t scoreThreshold);

    private:
        struct IntervalsFile
        {
            std::string filename_;
            boost::shared_ptr<std::istream> f_;
            util::DelimitedFile parser_;
            IntervalRecord rec_;

            IntervalsFile(const std::string& fn, const reference::CrrFile& crr);
        };

        struct DnbsFile
        {
            std::string filename_;
            boost::shared_ptr<std::istream> f_;
            util::DelimitedFile parser_;
            DnbRecord rec_;

            DnbsFile(const std::string& fn, const reference::CrrFile& crr);
        };

        //! opens the interval file of the given chromosome and points cursor to the first interval
        void openIntervals(uint16_t chr);

        void openDnbs();
        void seekDnbs();
        bool nextDnbs();

        const reference::CrrFile& crr_;
        const GenomeMetadata& exp_;

        boost::scoped_ptr<IntervalsFile> intervalsFile_;

        uint16_t dnbsChromosome_;
        int32_t dnbsIntervalId_;
        boost::scoped_ptr<DnbsFile> dnbsFile_;
        std::vector<DnbRecord> dnbs_;

        bool inInterval_;
    };

} } // cgatools::cgdata

#endif // CGATOOLS_CGDATA_EVIDENCEREADER_HPP_
