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

#ifndef CGATOOLS_VARIANTS_CALL_HPP_
#define CGATOOLS_VARIANTS_CALL_HPP_ 1

//! @file Call.hpp

#include "cgatools/core.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include <string>

namespace cgatools { namespace variants {

    //! A struct that corresponds to a single line of a Complete
    //! Genomics variant file.
    struct Call
    {
    public:
        enum VarQuality
        {
            VAR_QUALITY_EMPTY = 0,
            VAR_QUALITY_LOW   = 1,
            VAR_QUALITY_HIGH  = 2
        };

        class VarQualityField : public util::DelimitedFieldParser
        {
        public:
            VarQualityField(const std::string& name, VarQuality* varQuality);

            void parse(const char* first, const char* last);

        private:
            VarQuality* varQuality_;
        };

        class TotalScoreField : public util::DelimitedFieldParser
        {
        public:
            TotalScoreField(
                const std::string& name,
                int32_t* varScoreVAF,
                int32_t* varScoreEAF);

            void parse(const char* first, const char* last);

        private:
            int32_t* varScoreVAF_;
            int32_t* varScoreEAF_;
        };

        //! The "all" haplotype == [1..ploidy_] haplotypes.
        static const uint16_t ALL_HAPLOTYPES;

        //! Unknown ploidy, normally in no-ref regions on genome.
        static const uint16_t UNKNOWN_PLOIDY;

        //! The in-memory score corresponding to empty score on disk.
        static const int32_t EMPTY_SCORE;

        //! Construct an empty call.
        Call ();

        //! Returns the called sequence as a string, replacing "=" with
        //! the corresponding reference sequence.
        std::string calledSequence(const reference::CrrFile& crr) const;

        //! Returns the reference sequence as a string, replacing "="
        //! with the corresponding reference sequence.
        std::string refSequence(const reference::CrrFile& crr) const;

        //! Returns true iff the alleleSeq_ is consistent with the
        //! reference.
        bool isRefConsistent(const reference::CrrFile& crr) const;

        //! Returns true iff this call's sequence sequence contains N or
        //! ? characters.
        bool hasNoCalls() const;

        //! Add fields to parser so that it can parse lines, putting
        //! data in this call.
        void addFieldParsers(util::DelimitedFile& df, const reference::CrrFile& crr);

        //! Writes the call out, using the given separator charactor.
        std::ostream& write(std::ostream& out,
                            const reference::CrrFile& crr,
                            const char sep = '\t') const;

        //! Get the string form of varQuality.
        static const char* varQualityToString(const VarQuality varQuality);

        //! Parse a varQuality string.
        static VarQuality parseVarQuality(const char* first, const char* last);

        //! Returns header for a call output.
        static std::string getHeader (const char separator = '\t');

        //! The locus column in the variant file.
        uint32_t locus_;

        //! The ploidy column in the variant file. The value is
        //! Call::UNKNOWN_PLOIDY if the variant file specifies "?".
        uint16_t ploidy_;

        //! The haplotype column in the variant file. The value is
        //! Call::ALL_HAPLOTYPES if the variant file specifies "all".
        uint16_t haplotype_;

        //! The chromosome, begin, and end columns in the variant file.
        reference::Range range_;

        //! The varType column in the variant file.
        std::string varType_;

        //! The reference column in the variant file.
        std::string reference_;

        //! The alleleSeq column in the variant file.
        std::string alleleSeq_;

        //! The varScoreVAF column in the variant file.
        int32_t varScoreVAF_;

        //! The varScoreEAF column in the variant file. For older var
        //! files that have totalScore instead of varScoreVAF and
        //! varScoreEAF, the totalScore is filled into the varScoreEAF
        //! column.
        int32_t varScoreEAF_;

        //! The varQuality column in the variant file. For old var files
        //! that have no varQuality, this field gets VAR_QUALITY_EMPTY.
        VarQuality varQuality_;

        //! The hapLink column in the variant file.
        std::string hapLink_;

        //! The xRef column in the variant file.
        std::string xRef_;
    };

    //! Convert some field values that exist only in old variation files
    //! to their contemporary equivalents
    void upgradeOldFieldValues(Call& c);

    std::ostream& operator<<(std::ostream& out, const Call::VarQuality& varQuality);

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_CALL_HPP_
