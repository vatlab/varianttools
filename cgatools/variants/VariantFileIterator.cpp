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
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"

#include <sstream>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace variants {

    using std::ostringstream;
    using std::string;
    using std::set;
    using std::vector;
    namespace ba = boost::algorithm;
    using namespace util;
    using namespace cgatools::cgdata;

    class XRefField : public DelimitedFieldParser
    {
    public:
        XRefField(const std::string& name, std::string* val1, std::string* val2)
            : DelimitedFieldParser(name),
              val1_(val1),
              val2_(val2)
        {
        }

        void parse(const char* first, const char* last)
        {
            val1_->assign(first, last);
            *val2_ = *val1_;
        }

    private:
        std::string* val1_;
        std::string* val2_;
    };

    class HomCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        bool passesFilter(const Locus& locus, const Call& call) const
        {
            if (locus.getCalls().size() != 2 || locus.getPloidy() != 2 ||
                locus.getCalls()[0].alleleSeq_ != locus.getCalls()[1].alleleSeq_)
                return false;
            return true;
        }
    };

    class HetCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        bool passesFilter(const Locus& locus, const Call& call) const
        {
            return !homFilter_.passesFilter(locus, call);
        }

    private:
        HomCallFilter homFilter_;
    };

    class VarTypeCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        VarTypeCallFilter(const std::string& varType)
            : varType_(varType)
        {
        }

        bool passesFilter(const Locus& locus, const Call& call) const
        {
            return varType_ == call.varType_;
        }

    private:
        std::string varType_;
    };

    class VarScoreVafCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        VarScoreVafCallFilter(int32_t minScore)
            : minScore_(minScore)
        {
        }

        bool passesFilter(const Locus& locus, const Call& call) const
        {
            return call.varScoreVAF_ < minScore_;
        }

    private:
        int32_t minScore_;
    };

    class VarScoreEafCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        VarScoreEafCallFilter(int32_t minScore)
            : minScore_(minScore)
        {
        }

        bool passesFilter(const Locus& locus, const Call& call) const
        {
            return call.varScoreEAF_ < minScore_;
        }

    private:
        int32_t minScore_;
    };

    class VarQualityCallFilter : public VariantFileIterator::CallFilter
    {
    public:
        VarQualityCallFilter(Call::VarQuality excludeQuality)
            : excludeQuality_(excludeQuality)
        {
        }

        bool passesFilter(const Locus& locus, const Call& call) const
        {
            return call.varQuality_ != excludeQuality_;
        }

    private:
        Call::VarQuality excludeQuality_;
    };

    VariantFileIterator::VariantFileIterator(const CrrFile& crr)
        : crr_(&crr),
          locus_(crr),
          emptyFile_(true),
          eof_(true),
          hasPending_(false),
          referenceCoverValidation_(true),
          headerEnd_(0),
          fileEnd_(0),
          oneLinePerLocus_(false)
    {
    }

    void VariantFileIterator::open(const std::string& fnWithFilters)
    {
        // Parse filters.
        parseFilters(fnWithFilters);

        istream_ = InputStream::openCompressedInputStreamByExtension(name_);
        df_.reset(new DelimitedFile(*istream_, name_));
        detectFileType();
        if (oneLinePerLocus_)
        {
            addOlplParsers();
            locus_.setLocusAnnotations(
                LocusAnnotations(annotationColumnHeaders_,
                                 ! (df_->hasField("xRef") &&
                                    !(df_->hasField("allele1XRef") || df_->hasField("allele2XRef")))));
        }
        else
        {
            pendingCall_.addFieldParsers(*df_, *crr_);
            allele1ReadCountOffset_         = locus_.extras_.size();
            allele2ReadCountOffset_         = locus_.extras_.size();
            referenceAlleleReadCountOffset_ = locus_.extras_.size();
            totalReadCountOffset_           = locus_.extras_.size();
            locus_.setLocusAnnotations(LocusAnnotations(annotationColumnHeaders_, true));
        }
        eof_ = false;

        readPending();
        emptyFile_ = !hasPending_;

        readLocus();
    }

    void VariantFileIterator::parseFilters(const std::string& fnWithFilters)
    {
        try
        {
            size_t poundIdx = fnWithFilters.find('#');
            if (string::npos == poundIdx)
            {
                name_ = fnWithFilters;
                filterStr_ = "";
                return;
            }

            name_ = fnWithFilters.substr(0, poundIdx);
            filterStr_ = fnWithFilters.substr(poundIdx+1);

            vector<string> orParts;
            boost::split(orParts, filterStr_, boost::is_any_of(","));
            BOOST_FOREACH(const string& orPart, orParts)
            {
                vector<string> andParts;
                boost::split(andParts, orPart, boost::is_any_of(":"));
                filters_.push_back(vector< boost::shared_ptr<CallFilter> >());
                BOOST_FOREACH(const string& andPart, andParts)
                    filters_.back().push_back(parseFilter(andPart));
            }
        }
        catch(std::exception& ee)
        {
            throw Exception("failed to parse filters: "+fnWithFilters+": "+ee.what());
        }
    }

    boost::shared_ptr<VariantFileIterator::CallFilter>
    VariantFileIterator::parseFilter(const std::string& andPart) const
    {
        boost::shared_ptr<CallFilter> result;
        if ("hom" == andPart)
            result.reset(new HomCallFilter());
        else if ("het" == andPart)
            result.reset(new HetCallFilter());
        else if (boost::starts_with(andPart, "varType="))
            result.reset(new VarTypeCallFilter(andPart.substr(8)));
        else if (boost::starts_with(andPart, "varScoreVAF<"))
            result.reset(new VarScoreVafCallFilter(
                             boost::lexical_cast<int32_t>(andPart.substr(12))));
        else if (boost::starts_with(andPart, "varScoreEAF<"))
            result.reset(new VarScoreEafCallFilter(
                             boost::lexical_cast<int32_t>(andPart.substr(12))));
        else if (boost::starts_with(andPart, "varQuality!="))
            result.reset(new VarQualityCallFilter(
                             Call::parseVarQuality(&andPart[12], &andPart[andPart.size()])));
        if (0 == result.get())
            throw Exception("unrecognized filter: "+andPart);
        return result;
    }

    void VariantFileIterator::close()
    {
        df_.reset(static_cast<DelimitedFile*>(0));
        istream_.reset(static_cast<std::istream*>(0));
    }

    void VariantFileIterator::readPending()
    {
        hasPending_ = df_->next();
        if (hasPending_ && !oneLinePerLocus_)
            upgradeOldFieldValues(pendingCall_);
    }

    void VariantFileIterator::readLocus()
    {
        if (oneLinePerLocus_)
            readOlplLocus();
        else
            readMultilineLocus();
        applyFilters();
    }

    void VariantFileIterator::applyFilters()
    {
        for(bool noCallApplied=true; noCallApplied; )
        {
            noCallApplied = false;
            BOOST_FOREACH(Call& call, locus_.getCalls())
            {
                // Don't filter unscored calls.
                if (Call::EMPTY_SCORE == call.varScoreVAF_)
                    continue;

                for(size_t ii=0; ii<filters_.size(); ii++)
                {
                    bool pass = true;
                    for(size_t jj=0; jj<filters_[ii].size(); jj++)
                    {
                        const CallFilter& filt = *(filters_[ii][jj]);
                        if (!filt.passesFilter(locus_, call))
                        {
                            pass = false;
                            break;
                        }
                    }
                    if (pass)
                    {
                        applyNoCall(call);
                        locus_.setType("");
                        noCallApplied = true;
                        break;
                    }
                }
            }
        }
    }

    void VariantFileIterator::applyNoCall(Call& call)
    {
        call.varType_ = "no-call";
        call.alleleSeq_ = "?";
        call.varScoreVAF_ = Call::EMPTY_SCORE;
        call.varScoreEAF_ = Call::EMPTY_SCORE;
        call.varQuality_  = Call::VAR_QUALITY_EMPTY;
        call.hapLink_ = "";
        call.xRef_ = "";
    }

    void VariantFileIterator::readMultilineLocus()
    {
        if (!hasPending_)
        {
            eof_ = true;
            return;
        }

        locus_.clearCalls();
        locus_.addCall(pendingCall_);
        for(;;)
        {
            readPending();
            if ( (!hasPending_) || locus_.getId() != pendingCall_.locus_ )
                break;
            locus_.addCall(pendingCall_);
        }

        try
        {
            locus_.initFromCalls();
        }
        catch(std::exception& ee)
        {
            error(ee.what());
        }

        typedef reference::Location Location;
        if (hasPending_)
        {
            if (pendingCall_.range_.beginLocation() < locus_.getRange().endLocation())
            {
                if (pendingCall_.range_.chromosome_ != locus_.getRange().chromosome_)
                    locusCallError("Variant file chromosome order differs from reference chromosome order",
                                   locus_, pendingCall_);
                else
                    locusCallError("Variant file loci not well ordered", locus_, pendingCall_);
            }
            if (referenceCoverValidation_ &&
                locus_.getRange().endLocation() != pendingCall_.range_.beginLocation())
            {
                // This should only happen at chromosome ends.
                Location locusEnd = locus_.getRange().endLocation();
                Location pendingBegin = pendingCall_.range_.beginLocation();
                if (locusEnd.chromosome_ >= pendingBegin.chromosome_ ||
                    pendingBegin.offset_ != 0)
                {
                    locusCallError("Not all reference bases covered by calls", locus_, pendingCall_);
                }

                // Allow calls to be omitted that would cover no-calls
                // at the end of a chromosome. Otherwise, calls need to
                // cover each base of the reference.
                if (crr_->listChromosomes()[locusEnd.chromosome_].length() != locusEnd.offset_)
                {
                    const reference::CompactDnaSequence& ch = crr_->listChromosomes()[locusEnd.chromosome_];
                    for(size_t ii=locusEnd.offset_; ii<ch.length(); ii++)
                    {
                        if (ch.getBase(ii) != 'N')
                            locusCallError("Not all reference bases covered by calls", locus_, pendingCall_);
                    }
                }
            }
        }
        else
        {
            Location locusEnd = locus_.getRange().endLocation();
            if ( referenceCoverValidation_ &&
                 locusEnd.offset_ != crr_->listChromosomes()[locusEnd.chromosome_].length() )
            {
                // Allow calls to be omitted that would cover no-calls
                // at the end of a chromosome. Otherwise, calls need to
                // cover each base of the reference.
                const reference::CompactDnaSequence& ch = crr_->listChromosomes()[locusEnd.chromosome_];
                for(size_t ii=locusEnd.offset_; ii<ch.length(); ii++)
                {
                    if (ch.getBase(ii) != 'N')
                        locusError("Variant file ends before chromosome ends", locus_);
                }
            }
        }
    }

    void VariantFileIterator::error(const std::string& msg) const
    {
        throw Exception(name_ + ": " + msg);
    }

    void VariantFileIterator::locusCallError(
        const std::string& msg, const Locus& locus, const Call& call) const
    {
        ostringstream out;
        out << locus;
        call.write(out, *crr_);
        throw Exception(name_ + ": " + msg + ":\n" + out.str());
    }

    void VariantFileIterator::locusError(const std::string& msg, const Locus& locus) const
    {
        ostringstream out;
        out << locus;
        throw Exception(name_ + ": " + msg + ":\n" + out.str());
    }

    void VariantFileIterator::getReadCounts(int32_t* allele1ReadCount,
                                            int32_t* allele2ReadCount,
                                            int32_t* referenceAlleleReadCount,
                                            int32_t* totalReadCount,
                                            cgdata::EvidenceReader& evidence,
                                            const Locus& locus) const
    {
        *allele1ReadCount = -1;
        *allele2ReadCount = -1;
        *referenceAlleleReadCount = -1;
        *totalReadCount = -1;
        if (allele1ReadCountOffset_ < locus.extras_.size() &&
            allele2ReadCountOffset_ < locus.extras_.size() &&
            referenceAlleleReadCountOffset_ < locus.extras_.size() &&
            totalReadCountOffset_ < locus.extras_.size())
        {
            if ("" != locus.extras_[allele1ReadCountOffset_])
                *allele1ReadCount = parseValue<uint32_t>(locus.extras_[allele1ReadCountOffset_]);
            if ("" != locus.extras_[allele2ReadCountOffset_])
                *allele2ReadCount = parseValue<uint32_t>(locus.extras_[allele2ReadCountOffset_]);
            if ("" != locus.extras_[referenceAlleleReadCountOffset_])
                *referenceAlleleReadCount =
                    parseValue<uint32_t>(locus.extras_[referenceAlleleReadCountOffset_]);
            if ("" != locus.extras_[totalReadCountOffset_])
                *totalReadCount = parseValue<uint32_t>(locus.extras_[totalReadCountOffset_]);
            return;
        }

        locus.computeReadCounts(allele1ReadCount,
                                allele2ReadCount,
                                referenceAlleleReadCount,
                                totalReadCount,
                                evidence);
    }

    bool VariantFileIterator::hasAnnotation(const std::string& name) const
    {
        return annotationColumnHeaders_.size() != getAnnotationIndex(name);
    }

    size_t VariantFileIterator::getAnnotationIndex(const std::string& name) const
    {
        std::vector<std::string>::const_iterator p =
                std::find(annotationColumnHeaders_.begin(),
                          annotationColumnHeaders_.end(),
                          name);
        if (p == annotationColumnHeaders_.end())
            return annotationColumnHeaders_.size();
        else
            return static_cast<size_t>(&(*p) - &(annotationColumnHeaders_[0]));
    }

    void VariantFileIterator::detectFileType()
    {
        const DelimitedFileMetadata& m = df_->getMetadata();
        oneLinePerLocus_ = m.hasKey("TYPE") && (m.get("TYPE") == "VAR-OLPL");
    }

    //-------- OLPL specific code -------------------------------------------

    void VariantFileIterator::readOlplLocus()
    {
        const reference::Range& lr = locus_.getRange();

        // Handle the very end of the file
        if (!hasPending_)
        {
            if (emptyFile_)
            {
                // Translate empty files to empty files
                eof_ = true;
            }
            else
            {
                uint32_t chrEnd = crr_->listChromosomes()[lr.chromosome_].length();
                if (lr.end_ < chrEnd)
                    // Finish the last chromosome
                    createNoCallLocus(reference::Range(lr.chromosome_, lr.end_, chrEnd));
                else
                    eof_ = true;
            }
            return;
        }

        const reference::Range& nextr = pendingLocusLine_.range_;

        // Handle chromosome switches between current locus and the next line
        if (lr.chromosome_ != nextr.chromosome_)
        {
            if (lr.chromosome_ > nextr.chromosome_)
                error("Variant file chromosome order differs from reference chromosome order");

            uint32_t chrEnd = crr_->listChromosomes()[lr.chromosome_].length();
            if (lr.end_ < chrEnd)
            {
                // Finish the previous chromosome
                createNoCallLocus(reference::Range(lr.chromosome_, lr.end_, chrEnd));
                return;
            }
            else if (nextr.begin_ > 0)
            {
                // Fill the gap before the first locus in the new chromosome
                createNoCallLocus(reference::Range(nextr.chromosome_, 0, nextr.begin_));
                return;
            }
        }

        // Handle the gap between the last locus and the next line
        if (lr.chromosome_ == nextr.chromosome_ && lr.end_ != nextr.begin_)
        {

            if (lr.end_ > nextr.begin_)
                error("Variant file is out of order");
            createNoCallLocus(reference::Range(lr.chromosome_, lr.end_, nextr.begin_));
            return;
        }

        // Copy the line contents to the locus_
        fillOlplLocus();

        // Parse the next line of data
        readPending();
    }

    void VariantFileIterator::fillOlplLocus()
    {
        const LocusData& d = pendingLocusLine_;

        Call c;
        c.locus_ = d.locus_;
        c.ploidy_ = d.ploidy_;
        c.range_ = d.range_;
        c.reference_ = d.reference_;
        c.xRef_ = "";

        locus_.clearCalls();
        locus_.setType(d.varType_);

        if (isPureNoCall(d))
        {
            c.haplotype_ = Call::ALL_HAPLOTYPES;
            c.alleleSeq_ = "?";
            c.varScoreVAF_ = Call::EMPTY_SCORE;
            c.varScoreEAF_ = Call::EMPTY_SCORE;
            c.varQuality_  = Call::VAR_QUALITY_EMPTY;
            c.hapLink_.clear();
            c.varType_ = d.varType_ == "complex" ? "no-call" : d.varType_;
            locus_.addCall(c);
        }
        else if (d.varType_ == "ref" && (d.zygosity_ == "hom" || d.zygosity_ == "hap"))
        {
            c.haplotype_ = Call::ALL_HAPLOTYPES;
            c.alleleSeq_ = "=";
            c.varScoreVAF_ = Call::EMPTY_SCORE;
            c.varScoreEAF_ = Call::EMPTY_SCORE;
            c.varQuality_  = Call::VAR_QUALITY_EMPTY;
            c.hapLink_.clear();
            c.varType_ = d.varType_;
            locus_.addCall(c);
        }
        else
        {
            for (uint16_t hap = 0; hap < d.ploidy_; ++hap)
            {
                c.haplotype_ = hap + 1;
                c.alleleSeq_ = d.alleleSeq_[hap];
                c.varScoreVAF_ = d.alleleVarScoreVAF_[hap];
                c.varScoreEAF_ = d.alleleVarScoreEAF_[hap];
                c.varQuality_  = d.alleleVarQuality_ [hap];
                c.hapLink_ = d.hapLink_[hap];
                c.xRef_    = d.xRef_[hap];

                if (c.alleleSeq_.find('N') != std::string::npos ||
                    c.alleleSeq_.find('?') != std::string::npos )
                    c.varType_ = "no-call";
                else if (c.alleleSeq_ == c.reference_)
                    c.varType_ = "ref";
                else if (d.varType_ != "complex")
                    c.varType_ = d.varType_;
                else if (c.alleleSeq_.empty())
                    c.varType_ = "del";
                else if (c.reference_.empty())
                    c.varType_ = "ins";
                else if (1 == c.range_.length() && 1 == c.alleleSeq_.size())
                    c.varType_ = "snp";
                else
                    c.varType_ = "sub";
                locus_.addCall(c);
            }
        }

        try
        {
            locus_.initFromCalls();
        }
        catch(std::exception& ee)
        {
            error(ee.what());
        }

        CGA_ASSERT(locus_.extras_.size() == d.extras_.size());
        locus_.extras_ = d.extras_;
    }

    bool VariantFileIterator::isPureNoCall(const LocusData& d)
    {
        for (size_t hap = 0; hap < d.ploidy_; ++hap)
        {
            if (d.alleleSeq_[hap] != "?")
                return false;
        }
        return true;
    }

    void VariantFileIterator::createNoCallLocus(const reference::Range& r)
    {
        Call c;
        c.locus_ = 0;
        c.ploidy_ = 1;
        c.haplotype_ = Call::ALL_HAPLOTYPES;
        c.range_ = r;
        c.varType_ = "no-call";
        c.reference_ = "=";
        c.alleleSeq_ = "?";
        c.varScoreVAF_ = Call::EMPTY_SCORE;
        c.varScoreEAF_ = Call::EMPTY_SCORE;
        c.varQuality_  = Call::VAR_QUALITY_EMPTY;
        c.hapLink_.clear();
        c.xRef_.clear();

        locus_.clearCalls();
        locus_.addCall(c);
        locus_.initFromCalls();
        locus_.setType("");

        BOOST_FOREACH(std::string& ex, locus_.extras_)
            ex.clear();
    }

    void VariantFileIterator::addOlplParsers()
    {
        DelimitedFile& df = *df_;
        LocusData& d = pendingLocusLine_;

        df.addField(ValueField<uint32_t>("locus", &d.locus_));
        df.addField(ValueField<uint16_t>("ploidy", &d.ploidy_));
        df.addField(reference::ChromosomeIdField("chromosome", &d.range_.chromosome_, *crr_));
        df.addField(ValueField<uint32_t>("begin", &d.range_.begin_));
        df.addField(ValueField<uint32_t>("end", &d.range_.end_));
        df.addField(StringField("zygosity", &d.zygosity_));
        df.addField(StringField("varType", &d.varType_));
        df.addField(StringField("reference", &d.reference_));
        df.addField(StringField("allele1Seq", &d.alleleSeq_[0]));
        df.addField(StringField("allele2Seq", &d.alleleSeq_[1]));
        if ( df.hasField("allele1Score") &&
             !(df.hasField("allele1VarScoreEAF") || df.hasField("allele1VarScoreVAF")) )
        {
            df.addField(Call::TotalScoreField("allele1Score",
                                              &d.alleleVarScoreVAF_[0],
                                              &d.alleleVarScoreEAF_[0]));
        }
        else
        {
            df.addField(ValueField<int32_t>("allele1VarScoreVAF", &d.alleleVarScoreVAF_[0]).
                        exception("", Call::EMPTY_SCORE));
            df.addField(ValueField<int32_t>("allele1VarScoreEAF", &d.alleleVarScoreEAF_[0]).
                        exception("", Call::EMPTY_SCORE));
        }
        if ( df.hasField("allele2Score") &&
             !(df.hasField("allele2VarScoreEAF") || df.hasField("allele2VarScoreVAF")) )
        {
            df.addField(Call::TotalScoreField("allele2Score",
                                              &d.alleleVarScoreVAF_[1],
                                              &d.alleleVarScoreEAF_[1]));
        }
        else
        {
            df.addField(ValueField<int32_t>("allele2VarScoreVAF", &d.alleleVarScoreVAF_[1]).
                        exception("", Call::EMPTY_SCORE));
            df.addField(ValueField<int32_t>("allele2VarScoreEAF", &d.alleleVarScoreEAF_[1]).
                        exception("", Call::EMPTY_SCORE));
        }
        df.addField(Call::VarQualityField("allele1VarQuality", &d.alleleVarQuality_[0]),
                    DelimitedFile::FPT_OPTIONAL);
        df.addField(Call::VarQualityField("allele2VarQuality", &d.alleleVarQuality_[1]),
                    DelimitedFile::FPT_OPTIONAL);
        df.addField(StringField("allele1HapLink", &d.hapLink_[0]));
        df.addField(StringField("allele2HapLink", &d.hapLink_[1]));

        StringSet callColumns("locus,ploidy,chromosome,begin,end,zygosity,"
                              "varType,reference,allele1Seq,allele2Seq,"
                              "allele1Score,allele2Score,"
                              "allele1VarScoreVAF,allele2VarScoreVAF,"
                              "allele1VarScoreEAF,allele2VarScoreEAF,"
                              "allele1VarQuality,allele2VarQuality,"
                              "allele1HapLink,allele2HapLink", "", "");
        if (df.hasField("xRef") && !(df.hasField("allele1XRef") || df.hasField("allele2XRef")))
        {
            df.addField(XRefField("xRef", &d.xRef_[0], &d.xRef_[1]));
            callColumns.insert("xRef");
        }
        else
        {
            df.addField(StringField("allele1XRef", &d.xRef_[0]));
            df.addField(StringField("allele2XRef", &d.xRef_[1]));
            callColumns.insert("allele1XRef");
            callColumns.insert("allele2XRef");
        }

        BOOST_FOREACH(const std::string& hdr, df.getColumnHeaders())
        {
            if (callColumns.contains(hdr))
                continue;
            annotationColumnHeaders_.push_back(hdr);
        }

        d.extras_.resize(annotationColumnHeaders_.size());
        for (size_t ii = 0; ii < annotationColumnHeaders_.size(); ++ii)
            df.addField(StringField(annotationColumnHeaders_[ii], &d.extras_[ii]));

        allele1ReadCountOffset_ = getAnnotationIndex("allele1ReadCount");
        allele2ReadCountOffset_ = getAnnotationIndex("allele2ReadCount");
        referenceAlleleReadCountOffset_ = getAnnotationIndex("referenceAlleleReadCount");
        totalReadCountOffset_ = getAnnotationIndex("totalReadCount");
    }

    void VariantFileIterator::fillOlplFileMetadata(util::DelimitedFile::Metadata& meta) const
    {
        const util::DelimitedFile::Metadata& vfMeta = df_->getMetadata();

        // Transfer all keys over from the original file
        meta = vfMeta;

        // Replace the producer information
        meta.initDefaults();
        meta.set("TYPE", "VAR-OLPL");
    }

    void VariantFileIterator::writeOlplFileHeader(std::ostream& out) const
    {
        Locus::writeOneLineFileHeader(out, xRefIsAlleleSpecific());
        BOOST_FOREACH(const std::string& exh, annotationColumnHeaders_)
            out << '\t' << exh;
    }

} } // cgatools::variants
