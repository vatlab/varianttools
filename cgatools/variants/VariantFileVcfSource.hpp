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

#ifndef CGATOOLS_VARIANTS_VARIANTFILEVCFSOURCE_HPP_
#define CGATOOLS_VARIANTS_VARIANTFILEVCFSOURCE_HPP_ 1

//! @file VariantFileVcfSource.hpp

#include "cgatools/core.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/variants/SuperlocusIterator.hpp"

#include <deque>
#include <map>
#include <set>
#include <vector>

// Forward declarations
namespace cgatools { namespace cgdata {
    class GenomeMetadata;
} } // cgatools::cgdata

namespace cgatools { namespace variants {

    // Forward declarations
    class SubFieldAnnotation;
    struct SubFieldCall;
    namespace calib {
        class CalibratedScorer;
    }

    class VariantFileVcfRecordWriter : public cgatools::conv::VcfRecordWriter
    {
    public:
        VariantFileVcfRecordWriter(
            const std::vector< boost::shared_ptr<VariantFileIterator> >& var,
            const std::vector< std::string> fieldNames,
            const cgatools::reference::CrrFile& crr,
            const std::string& calibPrefix);

        // Overridden methods of VcfRecordWriter.
        cgatools::reference::Location getLocation() const;
        void writeRef(std::ostream& out) const;
        void writeAlt(std::ostream& out) const;
        void writeInfo(std::ostream& out) const;
        void writeFormat(std::ostream& out) const;
        void writeSample(std::ostream& out, size_t gIdx) const;

        // Other public methods etc.

        struct VcfData
        {
            reference::Range range_;

            // alleles_[alleleIdx_[0][hap]] <= sequence of haplotype hap of file 0
            // alleles_[alleleIdx_[1][hap]] <= sequence of haplotype hap of file 1
            std::vector<std::string> alleles_;
            std::vector< std::vector<size_t> > alleleIdx_;
            std::vector<uint32_t> ps_; // ps_[0] <= PS for file 0
            std::vector< std::vector<int32_t> > hq_;
            std::vector< std::vector<int32_t> > ehq_;
            std::vector< std::vector<int32_t> > cehq_; // calibrated ehq
        };

        const VcfData& getRecord() const;

        void setLocus(
            const cgatools::reference::Range& range,
            const Superlocus& sl,
            std::vector<PhasedHypothesis>& hyp);

        size_t getSampleCount() const;

        const char* getSomaticStatusAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls) const;

        void addSubFieldHeaderRecords(
            std::vector<cgatools::conv::VcfSubFieldHeaderRecord>& result) const;

        size_t getAlleleCount(size_t alleleIdx) const;

    private:
        struct HaploEntry
        {
            HaploEntry()
                : pos_(0),
                  lastSeenLocusCount_(0),
                  allele1_(0)
            {
            }

            HaploEntry(uint32_t pos, uint32_t locusCount, bool allele1)
                : pos_(pos),
                  lastSeenLocusCount_(locusCount),
                  allele1_(allele1)
            {
            }

            uint32_t pos_;
            uint32_t lastSeenLocusCount_;
            bool allele1_;
        };

        void addAnnotation(std::vector< boost::shared_ptr<SubFieldAnnotation> >& anns,
                           SubFieldAnnotation* pAnn);

        std::string getRefColumn() const;

        void fixN(std::string& result) const;

        bool isEmptyAnn(const std::string& ann) const;

        uint32_t swapIfNeeded(
            size_t gIdx,
            PhasedHypothesis& hyp,
            std::map< std::string, HaploEntry >& haploMap) const;

        bool callBelongsInRecord(const Call& call) const;

        size_t addVcfAllele(const std::string& allele);

        void appendScores(const PhasedAllele& allele, size_t gIdx);

        double getCalibratedScore(size_t gIdx, const Locus* locus, const Call* call, bool eaf);

        bool isHomozygousAlt(const Locus* locus) const;

        std::string getLocusVarType(const Locus* locus);

        std::string sanitizeVarType(const std::string& varType);

        std::string getVcfAllele(const Superlocus& sl,
                                 const PhasedAllele& allele);

    private:
        const cgatools::reference::CrrFile& crr_;
        const std::vector< boost::shared_ptr<VariantFileIterator> >& var_;
        std::string calibPrefix_;
        std::vector<bool> hasCalibration_;
        std::map< std::string, boost::shared_ptr<calib::CalibratedScorer> > calibrations_;
        const Superlocus* sl_;
        VcfData rec_;
        std::vector< std::map< std::string, HaploEntry > > haploMap_;
        std::vector< std::vector<SubFieldCall> > infoCalls_;
        std::vector< std::vector< std::vector<SubFieldCall> > > gtCalls_;
        std::vector< boost::shared_ptr<SubFieldAnnotation> > infoAnn_;
        std::vector< boost::shared_ptr<SubFieldAnnotation> > gtAnn_;
        std::vector< std::vector<std::string> > gtAnnData_;
        std::vector<bool> gtAnnFlag_;
        std::vector< std::string > asmIds_;
        std::vector< std::string > asmIdSuffixes_;
        std::vector< bool > hasSomatic_;
        std::vector< std::vector<size_t> > otherIdx_;
        uint32_t locusCount_;
    };

    class VariantFileVcfRecordSource : public cgatools::conv::VcfRecordSource
    {
    public:
        VariantFileVcfRecordSource(
            const std::vector< boost::shared_ptr<VariantFileIterator> >& var,
            const std::vector<std::string> fieldNames,
            const cgatools::reference::CrrFile& crr,
            const std::string& calibPrefix,
            bool includeNoCalls);

        // Overridden methods of VcfRecordSource.

        std::vector<cgatools::conv::VcfSubFieldHeaderRecord> getSubFieldHeaderRecords() const;
        std::string getSource(size_t idxGenome) const;
        std::vector<cgatools::conv::VcfKvHeaderRecord> getKeyValueHeaderRecords(size_t idxGenome) const;
        std::string getAssemblyId(size_t idxGenome) const;

        bool eof() const;
        cgatools::conv::VcfRecordSource& operator++();
        const cgatools::conv::VcfRecordWriter& operator*() const;
        const cgatools::conv::VcfRecordWriter* operator->() const;

    private:
        class OrderBySizeAscending
        {
        public:
            OrderBySizeAscending(
                const std::vector< std::vector<PhasedHypothesis> >& hypAll)
                : hypAll_(hypAll)
            {
            }

            bool operator()(size_t lhs, size_t rhs) const
            {
                if (hypAll_[lhs].size() != hypAll_[rhs].size())
                    return hypAll_[lhs].size() < hypAll_[rhs].size();
                return lhs < rhs;
            }

        private:
            const std::vector< std::vector<PhasedHypothesis> >& hypAll_;
        };

        void initSuperlocus();

        void splitForVcf(
            const cgatools::variants::Superlocus& sl,
            std::vector<cgatools::reference::Range>& ranges) const;
        void addVariantRanges(
            std::vector<cgatools::reference::Range>& ranges,
            const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
            std::deque<cgatools::variants::Locus>::const_iterator>& loci,
            std::vector<uint32_t>& insPositions) const;
        bool isInsertionAtLeft(
            const std::string& rr,
            const std::string& aa) const;
        bool isInsertionAtRight(
            const std::string& rr,
            const std::string& aa) const;
        void addNoCallPositions(
            std::vector<cgatools::reference::Range>& ranges,
            const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
            std::deque<cgatools::variants::Locus>::const_iterator>& loci,
            const std::vector<uint32_t>& insPositions) const;
        void addNoCallSplits(
            std::vector<uint32_t>& ncSplits,
            const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
            std::deque<cgatools::variants::Locus>::const_iterator>& loci) const;
        void splitAndAddRange(
            std::vector<cgatools::reference::Range>& ranges,
            const cgatools::reference::Range& range,
            const std::vector<uint32_t>& ncSplits) const;

        void phaseHaplotypes(
            const cgatools::variants::Superlocus& sl,
            std::vector<cgatools::variants::PhasedHypothesis>& hyp) const;
        int scoreHypothesis(
            const cgatools::variants::PhasedHypothesis& hyp,
            const std::set<std::string>& popAlleles) const;
        void addPopAlleles(
            const cgatools::variants::PhasedHypothesis& hyp,
            std::set<std::string>& popAlleles) const;

        const cgatools::reference::CrrFile& crr_;
        std::vector< boost::shared_ptr<VariantFileIterator> > var_;
        SuperlocusIterator slIt_;
        boost::shared_ptr<VariantFileVcfRecordWriter> writer_;
        std::vector<cgatools::reference::Range> ranges_;
        size_t rangesIdx_;
        std::vector<cgatools::variants::PhasedHypothesis> hyp_;
        bool includeNoCalls_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_VARIANTFILEVCFSOURCE_HPP_
