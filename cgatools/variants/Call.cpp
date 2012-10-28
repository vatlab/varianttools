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
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/variants/Call.hpp"

#include <sstream>
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace variants {

    using namespace cgatools::util;
    namespace bu = baseutil;
    using reference::ChromosomeIdField;
    using namespace std;

    const uint16_t Call::ALL_HAPLOTYPES = 0;
    const uint16_t Call::UNKNOWN_PLOIDY = 0;
    const int32_t  Call::EMPTY_SCORE    = 0;

    Call::VarQualityField::VarQualityField(const std::string& name, VarQuality* varQuality)
        : DelimitedFieldParser(name),
          varQuality_(varQuality)
    {
    }

    void Call::VarQualityField::parse(const char* first, const char* last)
    {
        *varQuality_ = Call::parseVarQuality(first, last);
    }

    Call::TotalScoreField::TotalScoreField(
        const std::string& name,
        int32_t* varScoreVAF,
        int32_t* varScoreEAF)
        : DelimitedFieldParser(name),
          varScoreVAF_(varScoreVAF),
          varScoreEAF_(varScoreEAF)
    {
    }

    void Call::TotalScoreField::parse(const char* first, const char* last)
    {
        if (first == last)
            *varScoreVAF_ = Call::EMPTY_SCORE;
        else
            *varScoreVAF_ = parseValue<int32_t>(first, last);
        *varScoreEAF_ = *varScoreVAF_;
    }

    Call::Call()
        : locus_(0),
          ploidy_(UNKNOWN_PLOIDY),
          haplotype_(ALL_HAPLOTYPES),
          range_(0,0,0),
          varScoreVAF_(EMPTY_SCORE),
          varScoreEAF_(EMPTY_SCORE),
          varQuality_(VAR_QUALITY_EMPTY)
    {
    }

    std::string Call::calledSequence(const reference::CrrFile& crr) const
    {
        if ("=" == alleleSeq_)
            return crr.getSequence(range_);

        return alleleSeq_;
    }

    std::string Call::refSequence(const reference::CrrFile& crr) const
    {
        if ("=" == reference_)
            return crr.getSequence(range_);

        return reference_;
    }

    bool Call::isRefConsistent(const reference::CrrFile& crr) const
    {
        if ("=" == alleleSeq_)
            return true;

        if ("?" == alleleSeq_)
            return true;

        if ("=" == reference_)
            return bu::isConsistent(alleleSeq_, crr.getSequence(range_));

        return bu::isConsistent(alleleSeq_, reference_);
    }

    bool Call::hasNoCalls() const
    {
        return alleleSeq_.find('?') != string::npos || alleleSeq_.find('N') != string::npos;
    }

    void Call::addFieldParsers(util::DelimitedFile& df, const reference::CrrFile& crr)
    {
        df.addField(ValueField<uint32_t>("locus", &locus_));
        df.addField(ValueField<uint16_t>("ploidy", &ploidy_).
                    exception("?", UNKNOWN_PLOIDY));
        if (df.hasField("haplotype") && !df.hasField("allele"))
            df.addField(ValueField<uint16_t>("haplotype", &haplotype_).
                        exception("all", ALL_HAPLOTYPES));
        else
            df.addField(ValueField<uint16_t>("allele", &haplotype_).
                        exception("all", ALL_HAPLOTYPES));
        df.addField(ChromosomeIdField("chromosome", &range_.chromosome_, crr));
        df.addField(ValueField<uint32_t>("begin", &range_.begin_));
        df.addField(ValueField<uint32_t>("end", &range_.end_));
        df.addField(StringField("varType", &varType_));
        df.addField(StringField("reference", &reference_));
        df.addField(StringField("alleleSeq", &alleleSeq_));
        if ( df.hasField("totalScore") &&
             !(df.hasField("varScoreVAF") || df.hasField("varScoreEAF")) )
        {
            df.addField(TotalScoreField("totalScore", &varScoreVAF_, &varScoreEAF_));
        }
        else
        {
            df.addField(ValueField<int32_t>("varScoreVAF", &varScoreVAF_).
                        exception("", EMPTY_SCORE));
            df.addField(ValueField<int32_t>("varScoreEAF", &varScoreEAF_).
                        exception("", EMPTY_SCORE));
        }
        df.addField(VarQualityField("varQuality", &varQuality_), DelimitedFile::FPT_OPTIONAL);
        df.addField(StringField("hapLink", &hapLink_), DelimitedFile::FPT_OPTIONAL);
        df.addField(StringField("xRef", &xRef_), DelimitedFile::FPT_OPTIONAL);
    }

    std::ostream& Call::write(std::ostream& out,
                              const reference::CrrFile& crr,
                              const char sep) const
    {
        out << locus_ << sep;
        if (UNKNOWN_PLOIDY == ploidy_)
            out << "?" << sep;
        else
            out << ploidy_ << sep;
        if (Call::ALL_HAPLOTYPES == haplotype_)
            out << "all" << sep;
        else
            out << haplotype_ << sep;
        out << crr.listChromosomes()[range_.chromosome_].getName() << sep
            << range_.begin_ << sep
            << range_.end_ << sep
            << varType_ << sep
            << reference_ << sep
            << alleleSeq_ << sep;
        // A varScoreVAF of 0 indicates this call wasn't
        // scored. Otherwise, it was scored, both for VAF and EAF.
        if (varScoreVAF_ != EMPTY_SCORE)
            out << varScoreVAF_;
        out << sep;
        if (varScoreVAF_ != EMPTY_SCORE)
            out << varScoreEAF_;
        out << sep
            << varQualityToString(varQuality_) << sep
            << hapLink_ << sep
            << xRef_;
        return out;
    }

    const char* Call::varQualityToString(const VarQuality varQuality)
    {
        if (VAR_QUALITY_EMPTY == varQuality)
            return "";
        else if (VAR_QUALITY_LOW == varQuality)
            return "VQLOW";
        else if (VAR_QUALITY_HIGH == varQuality)
            return "VQHIGH";
        else
        {
            // Unrecognized varQuality.
            CGA_ASSERT(false);
            return "";
        }
    }

    Call::VarQuality Call::parseVarQuality(const char* first, const char* last)
    {
        if (first == last)
            return VAR_QUALITY_EMPTY;
        else if (first+5 == last &&
                 boost::iequals(string(first,last), "VQLOW"))
            return VAR_QUALITY_LOW;
        else if (first+6 == last &&
                 boost::iequals(string(first,last), "VQHIGH"))
            return VAR_QUALITY_HIGH;
        else
            throw Exception("unrecognized varQuality: "+string(first, last));
    }

    std::string Call::getHeader(const char separator)
    {
        std::ostringstream out;
        out << "locus"
            << separator <<"ploidy"
            << separator <<"allele"
            << separator << "chromosome"
            << separator << "begin"
            << separator << "end"
            << separator << "varType"
            << separator << "reference"
            << separator << "alleleSeq"
            << separator << "varScoreVAF"
            << separator << "varScoreEAF"
            << separator << "varQuality"
            << separator << "hapLink"
            << separator << "xRef";
        return out.str();
    }

    void upgradeOldFieldValues(Call& c)
    {
        if (c.varType_ == "=")
            c.varType_ = "ref";
        else if (c.varType_ == "delins")
            c.varType_ = "sub";
    }

    std::ostream& operator<<(std::ostream& out, const Call::VarQuality& varQuality)
    {
        out << Call::varQualityToString(varQuality);
        return out;
    }

} } // cgatools::variants
