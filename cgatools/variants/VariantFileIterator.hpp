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

#ifndef CGATOOLS_VARIANTS_VARIANTFILEITERATOR_HPP_
#define CGATOOLS_VARIANTS_VARIANTFILEITERATOR_HPP_ 1

//! @file VariantFileIterator.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/variants/Locus.hpp"

namespace cgatools { namespace variants {

    class VariantFileIterator : boost::noncopyable
    {
        typedef reference::CrrFile CrrFile;

    public:
        VariantFileIterator(const CrrFile& crr);
        void open(const std::string& fn);
        void close();

        bool eof() const
        {
            return eof_;
        }

        const Locus& operator*() const
        {
            return locus_;
        }

        const Locus* operator->() const
        {
            return &locus_;
        }

        VariantFileIterator& operator++()
        {
            readLocus();
            return *this;
        }

        const std::string& getFileName()
        {
            return name_;
        }

        void error(const std::string& msg) const;
        void locusCallError(const std::string& msg, const Locus& locus, const Call& call) const;
        void locusError(const std::string& msg, const Locus& locus) const;

        void setReferenceCoverValidation(bool validate)
        {
            referenceCoverValidation_ = validate;
        }

        const util::DelimitedFile::Metadata& getMetadata() const
        {
            return df_->getMetadata();
        }

        const std::vector<std::string>& getColumnHeaders() const
        {
            return df_->getColumnHeaders();
        }

        const std::vector<std::string>& getAnnotationColumnHeaders() const
        {
            return annotationColumnHeaders_;
        }

        bool xRefIsAlleleSpecific() const
        {
            return locus_.xRefIsAlleleSpecific();
        }

        const std::string& getFilterString() const
        {
            return filterStr_;
        }

        //! Returns the requested read counts, as defined in the
        //! masterVar file, or -1 if this read count is empty as
        //! annotated in a masterVar file. If the input file does not
        //! have read count annotations, the read counts are computed
        //! from the evidence.
        void getReadCounts(int32_t* allele1ReadCount,
                           int32_t* allele2ReadCount,
                           int32_t* referenceAlleleReadCount,
                           int32_t* totalReadCount,
                           cgdata::EvidenceReader& evidence,
                           const Locus& locus) const;

        bool hasAnnotation(const std::string& name) const;

        //! Returns the index of a given annotation column in the Locus::extras_
        //! vector, or Locus::extras_.size() if the annotation is not present.
        size_t getAnnotationIndex(const std::string& name) const;

        void fillOlplFileMetadata(util::DelimitedFile::Metadata& meta) const;
        void writeOlplFileHeader(std::ostream& out) const;

        bool isOlpl() const
        {
            return oneLinePerLocus_;
        }

        class CallFilter
        {
        public:
            virtual ~CallFilter() { }
            virtual bool passesFilter(const Locus& locus, const Call& call) const = 0;
        };

    private:
        //! Locus data in OLPL file
        struct LocusData
        {
            LocusData()
                : locus_(0),
                  ploidy_(0),
                  range_(0,0,0)
            {
                alleleVarScoreVAF_.assign(Call::EMPTY_SCORE);
                alleleVarScoreEAF_.assign(Call::EMPTY_SCORE);
                alleleVarQuality_ .assign(Call::VAR_QUALITY_EMPTY);
            }

            uint32_t locus_;
            uint16_t ploidy_;
            reference::Range range_;
            std::string zygosity_;
            std::string varType_;
            std::string reference_;
            boost::array<std::string, 2> alleleSeq_;
            boost::array<int32_t, 2> alleleVarScoreVAF_;
            boost::array<int32_t, 2> alleleVarScoreEAF_;
            boost::array<Call::VarQuality, 2> alleleVarQuality_;
            boost::array<std::string, 2> hapLink_;
            boost::array<std::string, 2> xRef_;
            std::vector<std::string> extras_;
        };

        void applyFilters();
        void applyNoCall(Call& call);
        void parseFilters(const std::string& fnWithFilters);
        boost::shared_ptr<VariantFileIterator::CallFilter>
        parseFilter(const std::string& andPart) const;
        void readPending();
        void readLocus();
        void readMultilineLocus();
        void readOlplLocus();
        void fillOlplLocus();
        void detectFileType();
        void addOlplParsers();
        void createNoCallLocus(const reference::Range& r);
        bool isPureNoCall(const LocusData& d);

        const CrrFile* crr_;
        boost::shared_ptr<std::istream> istream_;
        boost::shared_ptr<util::DelimitedFile> df_;
        std::string name_;
        std::string filterStr_;
        Locus locus_;
        Call pendingCall_;
        LocusData pendingLocusLine_;
        bool emptyFile_;
        bool eof_;
        bool hasPending_;
        bool referenceCoverValidation_;
        boost::uint64_t headerEnd_;
        boost::uint64_t fileEnd_;

        std::vector< std::vector< boost::shared_ptr<CallFilter> > > filters_;

        bool oneLinePerLocus_;
        std::vector<std::string> annotationColumnHeaders_;
        size_t allele1ReadCountOffset_;
        size_t allele2ReadCountOffset_;
        size_t referenceAlleleReadCountOffset_;
        size_t totalReadCountOffset_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_VARIANTFILEITERATOR_HPP_
