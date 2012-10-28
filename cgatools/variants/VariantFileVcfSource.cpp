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
#include "cgatools/variants/VariantFileVcfSource.hpp"
#include "cgatools/variants/calib/CalibratedScorer.hpp"
#include "cgatools/util/parse.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

#define CGATOOLS_HAPLOID_GL_IS_TWO_VALUES 1

using namespace std;
using namespace cgatools::util;
using namespace cgatools::reference;
using namespace cgatools::conv;
using cgatools::variants::calib::CalibratedScorer;

using boost::lexical_cast;
using boost::shared_ptr;

namespace cgatools { namespace variants {

    //////////////////////////////////////////////////////////////////
    // AnnotationRetriever and derived classes.
    //////////////////////////////////////////////////////////////////

    class AnnotationRetriever
    {
    public:
        virtual ~AnnotationRetriever() { }
        virtual std::string getAnnotation(const Call& call, const Locus& locus) = 0;
    };

    class LocusAnnotationRetriever : public AnnotationRetriever
    {
    public:
        LocusAnnotationRetriever(const std::string& mvName)
            : mvNames_(1, mvName)
        {
        }

        LocusAnnotationRetriever(const std::string& mvName1, const std::string& mvName2)
        {
            mvNames_.push_back(mvName1);
            mvNames_.push_back(mvName2);
        }

        LocusAnnotationRetriever(const std::vector<std::string>& mvNames)
            : mvNames_(mvNames)
        {
        }

        std::string getAnnotation(const Call& call, const Locus& locus)
        {
            BOOST_FOREACH(const string& name, mvNames_)
            {
                if (locus.hasAnnotation(name))
                {
                    const string& ann = locus.getAnnotation(name);
                    if (ann.size() > 0)
                        return ann;
                }
            }

            return string();
        }

    private:
        std::vector< std::string > mvNames_;
    };

    class AlleleAnnotationRetriever : public AnnotationRetriever
    {
    public:
        AlleleAnnotationRetriever(const std::string& mvName)
        {
            mvName_[0] = vector<string>(1, "allele1"+mvName);
            mvName_[1] = vector<string>(1, "allele2"+mvName);
        }

        AlleleAnnotationRetriever(const std::string& mvName1, const std::string& mvName2)
        {
            mvName_[0].push_back("allele1"+mvName1);
            mvName_[0].push_back("allele1"+mvName2);
            mvName_[1].push_back("allele2"+mvName1);
            mvName_[1].push_back("allele2"+mvName2);
        }

        AlleleAnnotationRetriever(const std::vector<std::string>& mvNames)
        {
            BOOST_FOREACH(const string& mvName, mvNames)
            {
                mvName_[0].push_back("allele1"+mvName);
                mvName_[1].push_back("allele2"+mvName);
            }
        }

        std::string getAnnotation(const Call& call, const Locus& locus)
        {
            size_t hap = 0;
            if (1 == call.haplotype_)
                hap = 0;
            else if (2 == call.haplotype_)
                hap = 1;
            else
                return string();

            const vector<string>& mvNames = mvName_[hap];
            BOOST_FOREACH(const string& name, mvNames)
            {
                if (locus.hasAnnotation(name))
                {
                    const string& ann = locus.getAnnotation(name);
                    if (ann.size() > 0)
                        return ann;
                }
            }

            return string();
        }

    private:
        boost::array<std::vector<std::string>, 2> mvName_;
    };

    class XRefAnnotationRetriever : public AnnotationRetriever
    {
    public:
        std::string getAnnotation(const Call& call, const Locus& locus)
        {
            return call.xRef_;
        }
    };

    //////////////////////////////////////////////////////////////////
    // AnnotationAccumulator and derived classes.
    //////////////////////////////////////////////////////////////////

    class AnnotationAccumulator
    {
    public:
        virtual ~AnnotationAccumulator() { }
        virtual void init() = 0;
        virtual void addAnnotation(const std::string& ann) = 0;
        virtual std::string getAnnotation() = 0;
    };

    class JoinAnnotationAccumulator : public AnnotationAccumulator
    {
    public:
        JoinAnnotationAccumulator(const std::string& outSep)
            : outSep_(outSep)
        {
        }

        void init()
        {
            ann_.clear();
        }

        void addAnnotation(const std::string& ann)
        {
            vector<string> parts;
            boost::split(parts, ann, boost::is_any_of(";"));
            BOOST_FOREACH(string& part, parts)
            {
                boost::replace_all(part, ":", "|");
                boost::replace_all(part, " ", "_");
                ann_.insert(part);
            }
        }

        std::string getAnnotation()
        {
            if (0 == ann_.size())
                return ".";
            return boost::join(ann_, outSep_);
        }

    private:
        std::set<std::string> ann_;
        std::string outSep_;
    };

    class FirstNonEmptyAccumulator : public AnnotationAccumulator
    {
    public:
        void init()
        {
            ann_.clear();
        }

        void addAnnotation(const std::string& ann)
        {
            if (ann_.size() > 0)
                return;

            ann_ = ann;
        }

        std::string getAnnotation()
        {
            if (0 == ann_.size())
                return ".";
            return ann_;
        }

    private:
        std::string ann_;
    };

    template <typename TT>
    class MaxAnnotationAccumulator : public AnnotationAccumulator
    {
    public:
        MaxAnnotationAccumulator()
        {
            init();
        }

        void init()
        {
            val_ = std::numeric_limits<TT>::min();
        }

        void addAnnotation(const std::string& ann)
        {
            val_ = std::max(val_, parseValue<TT>(ann));
        }

        std::string getAnnotation()
        {
            if (std::numeric_limits<TT>::min() != val_)
                return lexical_cast<string>(val_);
            return ".";
        }

    private:
        TT val_;
    };

    //////////////////////////////////////////////////////////////////
    // SubFieldCall
    //////////////////////////////////////////////////////////////////

    struct SubFieldCall
    {
        SubFieldCall(
            const cgatools::variants::Locus* locus,
            const cgatools::variants::Call*  call)
            : locus_(locus),
              call_(call)
        {
        }

        const cgatools::variants::Locus* locus_;
        const cgatools::variants::Call*  call_;
    };

    //////////////////////////////////////////////////////////////////
    // SubFieldAnnotation and derived classes.
    //////////////////////////////////////////////////////////////////

    class SubFieldAnnotation
    {
    public:
        typedef VariantFileVcfRecordWriter::VcfData VcfData;

        virtual ~SubFieldAnnotation() { }

        SubFieldAnnotation(
            VcfSubFieldHeaderRecord::Key key,
            const std::string& id,
            const std::string& number,
            const std::string& type,
            const std::string& description)
            : key_(key),
              id_(id),
              number_(number),
              type_(type),
              description_(description)
        {
        }

        virtual cgatools::conv::VcfSubFieldHeaderRecord
        getHeaderRecord() const
        {
            return VcfSubFieldHeaderRecord(key_, id_, number_, type_, description_);
        }

        virtual void
        writeHeader(std::ostream& out) const
        {
            out << "<ID=" << id_
                << ",Number=" << number_
                << ",Type=" << type_
                << ",Description=\"" << description_
                << "\">";
        }

        virtual void
        writeId(std::ostream& out)
        {
            out << id_;
        }

        //! Gets the annotation for a subfield. A return value of "."
        //! indicates no information was found.
        //! @gIdx  The sample or genome offset within the list of
        //!        samples, or numeric_limits<size_t>::max() if this is
        //!        an INFO annotation.
        //! @calls The set of calls for this SubFieldAnnotation,
        //!        organized by allele then by call within the
        //!        allele. For INFO subfields, the allele refers to an
        //!        ALT allele. (The calls for that ALT allele, for all
        //!        genomes, are merged into the calls vector.) For
        //!        genotype subfields, the allele refers to an allele
        //!        within the genome.
        virtual std::string
        getAnnotation(size_t gIdx,
                      const std::vector< std::vector<SubFieldCall> >& calls) = 0;

        std::string getId() const
        {
            return id_;
        }

    private:
        VcfSubFieldHeaderRecord::Key key_;
        std::string id_;
        std::string number_;
        std::string type_;
        std::string description_;
    };

    class GlobalSubField : public SubFieldAnnotation
    {
    public:
        GlobalSubField(boost::shared_ptr<AnnotationRetriever> retriever,
                       boost::shared_ptr<AnnotationAccumulator> acc,
                       const std::vector<bool>& allowedGenomes,
                       VcfSubFieldHeaderRecord::Key key,
                       const std::string& id,
                       const std::string& number,
                       const std::string& type,
                       const std::string& description)
            : SubFieldAnnotation(key, id, number, type, description),
              allowedGenomes_(allowedGenomes),
              retriever_(retriever),
              acc_(acc)
        {
        }

        std::string
        getAnnotation(size_t gIdx,
                      const std::vector< std::vector<SubFieldCall> >& calls)
        {
            acc_->init();
            if (0 == allowedGenomes_.size() || allowedGenomes_[gIdx])
            {
                for(size_t ii=0; ii<calls.size(); ii++)
                {
                    for(size_t jj=0; jj<calls[ii].size(); jj++)
                    {
                        const Call&  call  = *(calls[ii][jj].call_);
                        const Locus& locus = *(calls[ii][jj].locus_);

                        string ann(retriever_->getAnnotation(call, locus));
                        if (!ann.empty())
                            acc_->addAnnotation(ann);
                    }
                }
            }

            return acc_->getAnnotation();
        }

    private:
        std::vector<bool> allowedGenomes_;
        boost::shared_ptr<AnnotationRetriever> retriever_;
        boost::shared_ptr<AnnotationAccumulator> acc_;
    };

    class AlleleSubField : public SubFieldAnnotation
    {
    public:
        AlleleSubField(boost::shared_ptr<AnnotationRetriever> retriever,
                       boost::shared_ptr<AnnotationAccumulator> acc,
                       const std::vector<bool>& allowedGenomes,
                       VcfSubFieldHeaderRecord::Key key,
                       const std::string& id,
                       const std::string& number,
                       const std::string& type,
                       const std::string& description,
                       const std::string& outAlleleSep = ",")
            : SubFieldAnnotation(key, id, number, type, description),
              allowedGenomes_(allowedGenomes),
              retriever_(retriever),
              acc_(acc),
              outAlleleSep_(outAlleleSep)
        {
        }

        std::string
        getAnnotation(size_t gIdx,
                      const std::vector< std::vector<SubFieldCall> >& calls)
        {
            bool found = false;
            vector<string> annByAllele;
            for(size_t ii=0; ii<calls.size(); ii++)
            {
                acc_->init();
                if (0 == allowedGenomes_.size() || allowedGenomes_[gIdx])
                {
                    for(size_t jj=0; jj<calls[ii].size(); jj++)
                    {
                        const Call&  call  = *(calls[ii][jj].call_);
                        const Locus& locus = *(calls[ii][jj].locus_);

                        string ann(retriever_->getAnnotation(call, locus));
                        if (!ann.empty())
                            acc_->addAnnotation(ann);
                    }
                }
                annByAllele.push_back(acc_->getAnnotation());
                found |= "." != annByAllele.back();
            }
            return boost::join(annByAllele, outAlleleSep_);
        }

    private:
        std::vector<bool> allowedGenomes_;
        boost::shared_ptr<AnnotationRetriever> retriever_;
        boost::shared_ptr<AnnotationAccumulator> acc_;
        std::string outAlleleSep_;
    };

    class FTSubField : public SubFieldAnnotation
    {
    public:
        FTSubField()
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT,
                                 "FT", "1", "String", "Genotype filters")
        {
        }

        std::string
        getAnnotation(size_t gIdx,
                      const std::vector< std::vector<SubFieldCall> >& calls)
        {
            bool found = false;
            bool vqLow = false;
            for(size_t ii=0; ii<calls.size(); ii++)
            {
                for(size_t jj=0; jj<calls[ii].size(); jj++)
                {
                    const Call& call = *(calls[ii][jj].call_);

                    if (Call::VAR_QUALITY_EMPTY == call.varQuality_)
                        continue;

                    found = true;
                    if (Call::VAR_QUALITY_HIGH != call.varQuality_)
                        vqLow = true;
                }
            }

            bool sqLow = false;
            for(size_t ii=0; ii<calls.size(); ii++)
            {
                for(size_t jj=0; jj<calls[ii].size(); jj++)
                {
                    const Locus& locus = *(calls[ii][jj].locus_);

                    if ( !(locus.hasAnnotation("somaticScore") &&
                           "" != locus.getAnnotation("somaticScore") &&
                           locus.hasAnnotation("somaticQuality")) )
                        continue;

                    found = true;
                    if ("SQHIGH" != locus.getAnnotation("somaticQuality"))
                        sqLow = true;
                }
            }

            if (!found)
                return ".";
            else if (vqLow && sqLow)
                return "VQLOW;SQLOW";
            else if (vqLow)
                return "VQLOW";
            else if (sqLow)
                return "SQLOW";
            else
                return "PASS";
        }
    };

    class NsSubField : public SubFieldAnnotation
    {
    public:
        NsSubField(const std::string& id,
                   const std::string& description,
                   VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_INFO, id, "1", "Integer", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            return lexical_cast<string>(conv_->getSampleCount());
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class AnSubField : public SubFieldAnnotation
    {
    public:
        AnSubField(const std::string& id,
                   const std::string& description,
                   VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_INFO, id, "1", "Integer", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            size_t count = 0;
            for(size_t ii=0; ii<=calls.size(); ii++)
                count += conv_->getAlleleCount(ii);
            return lexical_cast<string>(count);
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class AcSubField : public SubFieldAnnotation
    {
    public:
        AcSubField(const std::string& id,
                   const std::string& description,
                   VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_INFO, id, "A", "Integer", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            if (calls.size() < 1)
                return ".";

            vector<string> ac;
            for(size_t ii=1; ii<=calls.size(); ii++)
                ac.push_back(lexical_cast<string>(conv_->getAlleleCount(ii)));
            return boost::join(ac, ",");
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class SomaticStatusSubField : public SubFieldAnnotation
    {
    public:
        SomaticStatusSubField(const std::string& id,
                              const std::string& description,
                              VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "1", "String", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            return conv_->getSomaticStatusAnnotation(gIdx, calls);
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class GenotypeSubField : public SubFieldAnnotation
    {
    public:
        GenotypeSubField(const std::string& id,
                         const std::string& description,
                         VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "1", "String", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            const VcfData& rec = conv_->getRecord();
            vector<string> gt;
            for(size_t jj=0; jj<rec.alleleIdx_[gIdx].size(); jj++)
            {
                if (rec.alleleIdx_[gIdx][jj] >= rec.alleles_.size())
                    gt.push_back(".");
                else
                    gt.push_back(boost::lexical_cast<string>(rec.alleleIdx_[gIdx][jj]));
            }
            return boost::join(gt, rec.ps_[gIdx] != std::numeric_limits<uint32_t>::max() ? "|" : "/");
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class PhaseSetSubField : public SubFieldAnnotation
    {
    public:
        PhaseSetSubField(const std::string& id,
                         const std::string& description,
                         VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "1", "Integer", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            const VcfData& rec = conv_->getRecord();
            if (std::numeric_limits<uint32_t>::max() != rec.ps_[gIdx])
                return boost::lexical_cast<string>(rec.ps_[gIdx]);
            else
                return ".";
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class GenotypeQualitySubField : public SubFieldAnnotation
    {
    public:
        GenotypeQualitySubField(const std::string& id,
                                const std::string& description,
                                VariantFileVcfRecordWriter* conv)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "1", "Integer", description),
              conv_(conv)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            const VcfData& rec = conv_->getRecord();
            int32_t gq = std::numeric_limits<int32_t>::max();
            for(size_t jj=0; jj<rec.hq_[gIdx].size(); jj++)
                gq = std::min(gq, rec.hq_[gIdx][jj]);
            if (gq > 0 && gq < std::numeric_limits<int32_t>::max())
                return lexical_cast<string>(gq);
            else
                return ".";
        }

    private:
        VariantFileVcfRecordWriter* conv_;
    };

    class GenotypeLikelihoodSubField : public SubFieldAnnotation
    {
    public:
        GenotypeLikelihoodSubField(const std::string& id,
                                   const std::string& description,
                                   VariantFileVcfRecordWriter* conv,
                                   const std::vector< std::vector<int32_t> >& hq,
                                   const std::vector<bool>* hasScore = 0)
#if CGATOOLS_HAPLOID_GL_IS_TWO_VALUES
            // VCF validator appears to assume "G" is for diploid
            // genotypes, even though spec says GL haploid sites use
            // haploid genotypes. Should be fixed in future version of
            // VCF validator, in which case we can change "." to "G"
            // below.
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, ".", "Integer", description),
#else
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "G", "Integer", description),
#endif
              conv_(conv),
              hq_(hq),
              hasScore_(hasScore)
        {
        }

        struct ScoredPremise
        {
            ScoredPremise(size_t alleleIdx, size_t count, int32_t score)
                : alleleIdx_(alleleIdx),
                  count_(count),
                  score_(score)
            {
            }

            size_t alleleIdx_;
            size_t count_;
            int32_t score_;
        };

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            if (0 != hasScore_)
            {
                if ( !(*hasScore_)[gIdx] )
                    return ".";
            }

            const VcfData& rec = conv_->getRecord();
            const vector<int32_t>& hq       = hq_[gIdx];
            const vector<size_t>& alleleIdx = rec.alleleIdx_[gIdx];
            const vector<string>& alleles   = rec.alleles_;
            vector<ScoredPremise> premises;

#if CGATOOLS_HAPLOID_GL_IS_TWO_VALUES
            // VCF spec says "In presence of the GT field the same ploidy is
            // expected and the canonical order is used". This seems to
            // indicate two comma-separated values are expected for
            // haploid. But in that case, if we use a "G" as the Number
            // entry in the GL header, VCF validation fails. In that case, I
            // guess we had better say Number=".".
            if ( 1 == alleleIdx.size() )
            {
                // Haploid
                if ( alleleIdx[0] > alleles.size() || 0 == alleleIdx[0] )
                {
                    // Not called alt
                    return ".";
                }

                // VCF spec says "In presence of the GT field the same
                // ploidy is expected and the canonical order is used".
                CGA_ASSERT(1 == hq.size());
                vector<string> gl;
                for(size_t ii=0; ii<alleles.size(); ii++)
                    gl.push_back(ii == alleleIdx[0] ? "0" : lexical_cast<string>(-hq[0]));
                return boost::join(gl, ",");
            }
#else
            if ( 1 == alleleIdx.size() )
            {
                // Haploid
                if ( alleleIdx[0] > alleles.size() || 0 == alleleIdx[0] )
                {
                    // Not called alt
                    return ".";
                }

                CGA_ASSERT(1 == hq.size());
                premises.push_back(ScoredPremise(alleleIdx[0], 1, hq[0]));
            }
#endif

            else
            {
                // Diploid
                CGA_ASSERT( 2 == alleleIdx.size() );

                // Full no-call -> "."
                if ( alleleIdx[0] >= alleles.size() && alleleIdx[1] >= alleles.size() )
                    return ".";

                CGA_ASSERT( alleleIdx.size() == hq.size() );

                if ( 0 == rec.hq_[gIdx][0] && 0 == rec.hq_[gIdx][1] )
                    return "."; // Ref call, or score information otherwise unavailable.

                if ( alleleIdx[0] == alleleIdx[1] )
                {
                    // full call, hom
                    if ( ! (0 != alleleIdx[0] && alleleIdx[0] < alleles.size() ) )
                    {
                        std::cerr << "WARNING: homozygous ref call with scores! at "
                                  << (rec.range_.chromosome_+1)
                                  << ":" << (rec.range_.begin_+1)
                                  << " allele " << rec.alleles_[alleleIdx[0]] << std::endl;
                    }

                    int32_t minHq = std::min(hq[0], hq[1]);
                    int32_t maxHq = std::max(hq[0], hq[1]);

                    // At least one allele scored at maxHq
                    premises.push_back(ScoredPremise(alleleIdx[0], 1, maxHq));
                    // Two alleles scored at minHq
                    premises.push_back(ScoredPremise(alleleIdx[0], 2, minHq));
                }
                else
                {
                    // Het or call opposite no-call
                    for(size_t ii=0; ii<2; ii++)
                    {
                        if ( alleleIdx[ii] < alleles.size() )
                            premises.push_back(ScoredPremise(alleleIdx[ii], 1, hq[ii]));
                    }
                }
            }

            vector<string> gl;
            for(size_t gt1=0; gt1<alleles.size(); gt1++)
            {
                for(size_t gt0=0; gt0<=gt1; gt0++)
                {
                    boost::array<size_t, 2> gt;
                    gt[0] = gt0;
                    gt[1] = gt1;
                    int32_t gli = 0;
                    BOOST_FOREACH(const ScoredPremise& pp, premises)
                    {
                        size_t matchCount = std::count(gt.begin(), gt.end(), pp.alleleIdx_);
                        if (matchCount < pp.count_)
                            gli = std::min(gli, -pp.score_);
                    }
                    if ( 1 == alleleIdx.size() && gt[0] != gt[1] )
                        gl.push_back(".");
                    else
                        gl.push_back(lexical_cast<string>(gli));
                }
            }
            return boost::join(gl, ",");
        }

    private:
        VariantFileVcfRecordWriter* conv_;
        const std::vector< std::vector<int32_t> >& hq_;
        const std::vector<bool>* hasScore_;
    };

    class HaplotypeQualitySubField : public SubFieldAnnotation
    {
    public:
        HaplotypeQualitySubField(const std::string& id,
                                 const std::string& description,
                                 VariantFileVcfRecordWriter* conv,
                                 const std::vector< std::vector<int32_t> >& hq,
                                 const std::vector<bool>* hasScore = 0)
            : SubFieldAnnotation(VcfSubFieldHeaderRecord::VCF_FORMAT, id, "2", "Integer", description),
              conv_(conv),
              hq_(hq),
              hasScore_(hasScore)
        {
        }

        std::string
        getAnnotation(
            size_t gIdx,
            const std::vector< std::vector<SubFieldCall> >& calls)
        {
            if (0 != hasScore_)
            {
                if ( !(*hasScore_)[gIdx] )
                    return ".";
            }

            const VcfData& rec = conv_->getRecord();
            const vector<int32_t>& hqi = hq_[gIdx];
            vector<string> hq(2, ".");
            for(size_t jj=0; jj<2; jj++)
            {
                if (jj < rec.hq_[gIdx].size() && rec.hq_[gIdx][jj] > 0)
                    hq[jj] = lexical_cast<string>(hqi[jj]);
            }
            return boost::join(hq, ",");
        }

    private:
        VariantFileVcfRecordWriter* conv_;
        const std::vector< std::vector<int32_t> >& hq_;
        const std::vector<bool>* hasScore_;
    };


    //////////////////////////////////////////////////////////////////
    // VariantFileVcfRecordWriter
    //////////////////////////////////////////////////////////////////

    VariantFileVcfRecordWriter::VariantFileVcfRecordWriter(
        const std::vector< boost::shared_ptr<cgatools::variants::VariantFileIterator> >& var,
        const std::vector< std::string> fieldNames,
        const cgatools::reference::CrrFile& crr,
        const std::string& calibPrefix)
        : crr_(crr),
          var_(var),
          calibPrefix_(calibPrefix),
          hasCalibration_(var.size(), false),
          sl_(0),
          haploMap_(var.size()),
          locusCount_(0)
    {
        set<string> ach;
        bool hasLocusSpecificXRef = false;
        for(size_t ii=0; ii<var_.size(); ii++)
        {
            VariantFileIterator& vit = *(var_[ii]);
            const vector<string>& iiach = vit.getAnnotationColumnHeaders();
            BOOST_FOREACH(const string& ss, iiach)
                ach.insert(ss);
            if (!vit.xRefIsAlleleSpecific())
                hasLocusSpecificXRef = true;
        }

        typedef shared_ptr<AnnotationRetriever>   PAR;
        typedef shared_ptr<AnnotationAccumulator> PAA;
        const vector<string>& fn = fieldNames;
        if (std::find(fn.begin(), fn.end(), "GT") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new GenotypeSubField(
                    "GT",
                    "Genotype",
                    this));
        }
        if (std::find(fn.begin(), fn.end(), "PS") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new PhaseSetSubField(
                    "PS",
                    "Phase Set",
                    this));
        }
        if (std::find(fn.begin(), fn.end(), "NS") != fn.end())
        {
            addAnnotation(infoAnn_,
                          new NsSubField(
                              "NS",
                              "Number of Samples With Data",
                              this));
        }
        if (std::find(fn.begin(), fn.end(), "AN") != fn.end())
        {
            addAnnotation(infoAnn_,
                          new AnSubField(
                              "AN",
                              "Total number of alleles in called genotypes",
                              this));
        }
        if (std::find(fn.begin(), fn.end(), "AC") != fn.end())
        {
            addAnnotation(infoAnn_,
                          new AcSubField(
                              "AC",
                              "Allele count in genotypes, for each ALT allele",
                              this));
        }
        if (std::find(fn.begin(), fn.end(), "SS") != fn.end())
        {
            addAnnotation(gtAnn_,
                          new SomaticStatusSubField(
                              "SS",
                              "Somatic Status: Germline, Somatic, LOH, or . (Unknown)",
                              this));
        }
        if (std::find(fn.begin(), fn.end(), "FT") != fn.end())
        {
            addAnnotation(gtAnn_, new FTSubField());
        }
        if ( std::find(fn.begin(), fn.end(), "CGA_SXR") != fn.end() ||
             std::find(fn.begin(), fn.end(), "CGA_XR") != fn.end() )
        {
            if (hasLocusSpecificXRef)
            {
                addAnnotation(
                    gtAnn_,
                    new GlobalSubField(
                        PAR(new XRefAnnotationRetriever()),
                        PAA(new JoinAnnotationAccumulator(",")),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "CGA_SXR",
                        ".",
                        "String",
                        "Per-sample external database reference (dbSNP, COSMIC, etc)"));
            }
            else
            {
                addAnnotation(
                    infoAnn_,
                    new AlleleSubField(
                        PAR(new XRefAnnotationRetriever()),
                        PAA(new JoinAnnotationAccumulator("&")),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_INFO,
                        "CGA_XR",
                        "A",
                        "String",
                        "Per-ALT external database reference (dbSNP, COSMIC, etc)"));
            }
        }
        if (std::find(fn.begin(), fn.end(), "CGA_FI") != fn.end())
        {
            if (ach.find("allele1Gene") != ach.end() ||
                ach.find("allele2Gene") != ach.end())
            {
                addAnnotation(
                    infoAnn_,
                    new AlleleSubField(
                        PAR(new AlleleAnnotationRetriever("Gene")),
                        PAA(new JoinAnnotationAccumulator("&")),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_INFO,
                        "CGA_FI",
                        "A",
                        "String",
                        "Functional impact annotation"));
            }
        }
        if (std::find(fn.begin(), fn.end(), "GQ") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new GenotypeQualitySubField(
                    "GQ",
                    "Genotype Quality",
                    this));
        }
        if (std::find(fn.begin(), fn.end(), "HQ") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new HaplotypeQualitySubField(
                    "HQ",
                    "Haplotype Quality",
                    this,
                    rec_.hq_));
        }
        if (std::find(fn.begin(), fn.end(), "EHQ") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new HaplotypeQualitySubField(
                    "EHQ",
                    "Haplotype Quality, Equal Allele Fraction Assumption",
                    this,
                    rec_.ehq_));
        }
        if (!calibPrefix_.empty())
        {
            for(size_t ii=0; ii<var_.size(); ii++)
            {
                if (var_[ii]->hasAnnotation("totalReadCount"))
                    hasCalibration_[ii] = true;
            }
        }
        if (std::find(fn.begin(), fn.end(), "CGA_CEHQ") != fn.end())
        {
            if (calibPrefix_.empty())
                throw Exception("calibration-root must be specified for CGA_CEHQ annotation");
            addAnnotation(
                gtAnn_,
                new HaplotypeQualitySubField(
                    "CGA_CEHQ",
                    "Calibrated Haplotype Quality, Equal Allele Fraction Assumption",
                    this,
                    rec_.cehq_,
                    &hasCalibration_));
        }
        if (std::find(fn.begin(), fn.end(), "GL") != fn.end())
        {
            addAnnotation(
                gtAnn_,
                new GenotypeLikelihoodSubField(
                    "GL",
                    "Genotype Likelihood",
                    this,
                    rec_.hq_));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_CEGL") != fn.end())
        {
            if (calibPrefix_.empty())
                throw Exception("calibration-root must be specified for CGA_CEHQ annotation");
            addAnnotation(
                gtAnn_,
                new GenotypeLikelihoodSubField(
                    "CGA_CEGL",
                    "Calibrated Genotype Likelihood, Equal Allele Fraction Asssumption",
                    this,
                    rec_.cehq_,
                    &hasCalibration_));
        }
        if (std::find(fn.begin(), fn.end(), "DP") != fn.end())
        {
            if (ach.find("totalReadCount") != ach.end())
            {
                addAnnotation(
                    gtAnn_,
                    new GlobalSubField(
                        PAR(new LocusAnnotationRetriever("totalReadCount")),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "DP",
                        "1",
                        "Integer",
                        "Total Read Depth"));
            }
        }
        if (std::find(fn.begin(), fn.end(), "AD") != fn.end())
        {
            if (ach.find("allele1ReadCount") != ach.end() ||
                ach.find("allele2ReadCount") != ach.end())
            {
                addAnnotation(
                    gtAnn_,
                    new AlleleSubField(
                        PAR(new AlleleAnnotationRetriever("ReadCount")),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "AD",
                        "2",
                        "Integer",
                        "Allelic depths (number of reads in each observed allele)"));
            }
        }
        if (std::find(fn.begin(), fn.end(), "CGA_RDP") != fn.end())
        {
            if (ach.find("referenceAlleleReadCount") != ach.end())
            {
                addAnnotation(
                    gtAnn_,
                    new GlobalSubField(
                        PAR(new LocusAnnotationRetriever("referenceAlleleReadCount")),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        vector<bool>(),
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "CGA_RDP",
                        "1",
                        "Integer",
                        "Number of reads observed supporting the reference allele"));
            }
        }
        // Determine asm ids, asm id suffixes, and whether each
        // sample has somatic data.
        vector<string> allSuffixes;
        for(size_t ii=0; ii<var_.size(); ii++)
        {
            asmIds_.push_back(var_[ii]->getMetadata().get("ASSEMBLY_ID"));
            asmIdSuffixes_.push_back("");
            if ( (!boost::ends_with(asmIds_[ii], "-ASM")) &&
                 string::npos != asmIds_[ii].rfind('-') )
                asmIdSuffixes_.back() = asmIds_[ii].substr(asmIds_[ii].rfind('-'));
            if ( "" != asmIdSuffixes_[ii] &&
                 std::find(allSuffixes.begin(), allSuffixes.end(), asmIdSuffixes_[ii]) ==
                 allSuffixes.end() )
                allSuffixes.push_back(asmIdSuffixes_[ii]);
            const vector<string>& headers = var_[ii]->getColumnHeaders();
            hasSomatic_.push_back(
                std::find(headers.begin(), headers.end(), "somaticScore") != headers.end());
        }
        // Detect paired samples by ASM ID.
        otherIdx_.resize(var_.size());
        for(size_t ii=0; ii<otherIdx_.size(); ii++)
        {
            if ( asmIdSuffixes_[ii] != "" )
            {
                string asmIdPrefix = asmIds_[ii].substr(0, asmIds_[ii].size()-asmIdSuffixes_[ii].size());
                for(size_t jj=0; jj<otherIdx_.size(); jj++)
                {
                    if (hasSomatic_[ii] != hasSomatic_[jj] &&
                        asmIdSuffixes_[jj] != "" &&
                        asmIdPrefix == asmIds_[jj].substr(
                            0, asmIds_[jj].size()-asmIdSuffixes_[jj].size()))
                    {
                        otherIdx_[ii].push_back(jj);
                    }
                }
            }
        }
        // "Other" read counts.
        bool otherCountsExist = false;
        BOOST_FOREACH(const string& suffix, allSuffixes) {
            const char* rcAnn[] =
                {
                    "totalReadCount",
                    "allele1ReadCount",
                    "allele2ReadCount",
                    "referenceAlleleReadCount"
                };
            for(size_t jj=0; jj<sizeof(rcAnn)/sizeof(rcAnn[0]); jj++)
            {
                if (ach.find(rcAnn[jj]+suffix) != ach.end())
                    otherCountsExist = true;
            }
        }
        if (otherCountsExist)
        {
            vector<bool> allowedGenomes(var_.size(), false);
            for(size_t ii=0; ii<allowedGenomes.size(); ii++)
                allowedGenomes[ii] = 1 == otherIdx_[ii].size();
            vector<string> odpNames, oadNames, ordpNames;
            BOOST_FOREACH(const string& suffix, allSuffixes)
            {
                odpNames.push_back("totalReadCount"+suffix);
                oadNames.push_back("ReadCount"+suffix);
                ordpNames.push_back("referenceAlleleReadCount"+suffix);
            }
            if (std::find(fn.begin(), fn.end(), "CGA_ODP") != fn.end())
            {
                addAnnotation(
                    gtAnn_,
                    new GlobalSubField(
                        PAR(new LocusAnnotationRetriever(odpNames)),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        allowedGenomes,
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "CGA_ODP",
                        "1",
                        "Integer",
                        "Other Total Read Depth (read depth of other sample "
                        "in somatic comparison)"));
            }
            if (std::find(fn.begin(), fn.end(), "CGA_OAD") != fn.end())
            {
                addAnnotation(
                    gtAnn_,
                    new AlleleSubField(
                        PAR(new AlleleAnnotationRetriever(oadNames)),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        allowedGenomes,
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "CGA_OAD",
                        "2",
                        "Integer",
                        "Other Allelic depths (number of reads in other sample "
                        "in each observed allele of this sample)"));
            }
            if (std::find(fn.begin(), fn.end(), "CGA_ORDP") != fn.end())
            {
                addAnnotation(
                    gtAnn_,
                    new GlobalSubField(
                        PAR(new LocusAnnotationRetriever(ordpNames)),
                        PAA(new MaxAnnotationAccumulator<int32_t>()),
                        allowedGenomes,
                        VcfSubFieldHeaderRecord::VCF_FORMAT,
                        "CGA_ORDP",
                        "1",
                        "Integer",
                        "Other Reference Depth (number of reads in other sample "
                        "observed supporting the reference allele)"));
            }
        }
        if (std::find(fn.begin(), fn.end(), "CGA_PFAM") != fn.end() &&
            ach.find("pfam") != ach.end())
        {
            addAnnotation(
                infoAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("pfam")),
                    PAA(new JoinAnnotationAccumulator(",")),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_INFO,
                    "CGA_PFAM",
                    ".",
                    "String",
                    "PFAM Domain"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_MIRB") != fn.end() &&
            ach.find("miRBaseId") != ach.end())
        {
            addAnnotation(
                infoAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("miRBaseId")),
                    PAA(new JoinAnnotationAccumulator(",")),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_INFO,
                    "CGA_MIRB",
                    ".",
                    "String",
                    "miRBaseId"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_RPT") != fn.end() &&
            ach.find("repeatMasker") != ach.end())
        {
            addAnnotation(
                infoAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("repeatMasker")),
                    PAA(new JoinAnnotationAccumulator(",")),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_INFO,
                    "CGA_RPT",
                    ".",
                    "String",
                    "repeatMasker overlap information"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_SDO") != fn.end() &&
            ach.find("segDupOverlap") != ach.end())
        {
            addAnnotation(
                infoAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("segDupOverlap")),
                    PAA(new MaxAnnotationAccumulator<int32_t>()),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_INFO,
                    "CGA_SDO",
                    "1",
                    "Integer",
                    "Number of distinct segmental duplications that overlap this locus"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_SOMC") != fn.end() &&
            ach.find("somaticCategory") != ach.end())
        {
            addAnnotation(
                gtAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("somaticCategory")),
                    PAA(new FirstNonEmptyAccumulator()),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_FORMAT,
                    "CGA_SOMC",
                    "1",
                    "String",
                    "somaticCategory"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_SOMR") != fn.end() &&
            ach.find("somaticRank") != ach.end())
        {
            addAnnotation(
                gtAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("somaticRank")),
                    PAA(new FirstNonEmptyAccumulator()),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_FORMAT,
                    "CGA_SOMR",
                    "1",
                    "String",
                    "somaticRank"));
        }
        if (std::find(fn.begin(), fn.end(), "CGA_SOMS") != fn.end() &&
            ach.find("somaticScore") != ach.end())
        {
            addAnnotation(
                gtAnn_,
                new GlobalSubField(
                    PAR(new LocusAnnotationRetriever("somaticScore")),
                    PAA(new FirstNonEmptyAccumulator()),
                    vector<bool>(),
                    VcfSubFieldHeaderRecord::VCF_FORMAT,
                    "CGA_SOMS",
                    "1",
                    "Float",
                    "somaticScore"));
        }
    }

    Location VariantFileVcfRecordWriter::getLocation() const
    {
        return rec_.range_.beginLocation();
    }

    void VariantFileVcfRecordWriter::writeRef(std::ostream& out) const
    {
        if (rec_.range_.length() >= 200 && 1 == rec_.alleles_.size())
        {
            // No-call locus.
            out << rec_.alleles_[0].substr(0, 1);
            return;
        }
        out << rec_.alleles_[0];
    }

    void VariantFileVcfRecordWriter::writeAlt(std::ostream& out) const
    {
        if (rec_.range_.length() >= 200 && 1 == rec_.alleles_.size())
        {
            // No-call locus.
            out << "<CGA_NOCALL>";
            return;
        }
        if (rec_.alleles_.size() < 2)
        {
            out << ".";
            return;
        }
        for(size_t ii=1; ii<rec_.alleles_.size(); ii++)
        {
            if (ii > 1)
                out << ",";

            out << rec_.alleles_[ii];
        }
    }

    void VariantFileVcfRecordWriter::writeInfo(std::ostream& out) const
    {
        size_t countWritten = 0;
        if (rec_.range_.length() >= 200 && 1 == rec_.alleles_.size())
        {
            // No-call locus.  1-based position of last base of
            // no-call. So no need to increment position from our
            // 0-based position of the base after the no-call.
            out << "END=" << rec_.range_.end_;
            countWritten++;
        }
        if (0 == infoAnn_.size())
        {
            out << ".";
            return;
        }

        for(size_t ii=0; ii<infoAnn_.size(); ii++)
        {
            string ann = infoAnn_[ii]->getAnnotation(
                std::numeric_limits<size_t>::max(), infoCalls_);
            if (isEmptyAnn(ann))
                continue;

            if (countWritten > 0)
                out << ";";

            infoAnn_[ii]->writeId(out);
            out << "=";
            out << ann;
            countWritten++;
        }
        if (0 == countWritten)
            out << ".";
    }

    void VariantFileVcfRecordWriter::writeFormat(std::ostream& out) const
    {
        bool wroteOne = false;
        for(size_t ii=0; ii<gtAnn_.size(); ii++)
        {
            if (gtAnnFlag_[ii])
            {
                if (wroteOne)
                    out << ":";
                gtAnn_[ii]->writeId(out);
                wroteOne = true;
            }
        }
    }

    void VariantFileVcfRecordWriter::writeSample(std::ostream& out, size_t gIdx) const
    {
        CGA_ASSERT(gIdx < gtCalls_.size());
        bool wroteOne = false;
        for(size_t ii=0; ii<gtAnn_.size(); ii++)
        {
            if (gtAnnFlag_[ii])
            {
                if (wroteOne)
                    out << ":";
                out << gtAnnData_[gIdx][ii];
                wroteOne = true;
            }
        }
    }

    const VariantFileVcfRecordWriter::VcfData& VariantFileVcfRecordWriter::getRecord() const
    {
        return rec_;
    }

    void VariantFileVcfRecordWriter::setLocus(
        const cgatools::reference::Range& range,
        const cgatools::variants::Superlocus& sl,
        std::vector<cgatools::variants::PhasedHypothesis>& hyp)
    {
        locusCount_++;
        sl_ = &sl;
        rec_ = VcfData();
        rec_.range_ = range;
        rec_.alleleIdx_.resize(hyp.size());
        rec_.ps_ .resize(hyp.size(), std::numeric_limits<uint32_t>::max());
        rec_.hq_ .resize(hyp.size());
        rec_.ehq_.resize(hyp.size());
        rec_.cehq_.resize(hyp.size());

        infoCalls_.clear();
        gtCalls_.clear();
        gtCalls_.resize(hyp.size(), vector< vector<SubFieldCall> >(2));

        // Check on phasing.
        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            if (2 != hyp[ii].size())
                continue;

            rec_.ps_[ii] = swapIfNeeded(ii, hyp[ii], haploMap_[ii]);
        }

        addVcfAllele(getRefColumn());
        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            if (0 == hyp[ii].size())
            {
                // Multi-ploidy locus. No-call.
                rec_.alleleIdx_[ii].push_back(std::numeric_limits<size_t>::max());
                continue;
            }

            for(size_t jj=0; jj<hyp[ii].size(); jj++)
            {
                rec_.alleleIdx_[ii].push_back(
                    addVcfAllele(getVcfAllele(sl, hyp[ii][jj])));

                // HACK -- VCF validator appears to disallow any
                // locus where the length of the ALT does not equal
                // the length of the REF and the first REF base
                // doesn't equal first base of ALT, for each ALT. To
                // circumvent this, at this point we explicitly
                // check this. If the check fails, we must either
                // extend the VCF locus to the left by one base or
                // give up (i.e. report "."). We can only extend the
                // VCF locus to the left by one base without
                // stomping on the previous locus if this is the
                // beginning of the superlocus. This is rare (10's
                // of records per genome).
                if (rec_.alleleIdx_[ii].back() == rec_.alleles_.size()-1)
                {
                    string& rr  = rec_.alleles_[0];
                    string& alt = rec_.alleles_[rec_.alleleIdx_[ii].back()];
                    CGA_ASSERT(rr.size() > 0 && alt.size() > 0);

                    if (rr.size() != alt.size() && rr[0] != alt[0])
                    {
                        if (0 != range.begin_ &&
                            range.beginLocation() == sl.getRange().beginLocation())
                        {
                            Range newRange(range);
                            newRange.begin_--;
                            setLocus(newRange, sl, hyp);
                            return;
                        }
                        else
                        {
                            // Give up. No-call the allele.
                            rec_.alleles_.resize(rec_.alleles_.size()-1);
                            rec_.alleleIdx_[ii].back() = std::numeric_limits<size_t>::max();
                        }
                    }
                }

                infoCalls_.resize(rec_.alleles_.size()-1);

                for(size_t kk=0; kk<hyp[ii][jj].calls().size(); kk++)
                {
                    const Call*  call  = hyp[ii][jj].calls()[kk];
                    const Locus* locus = hyp[ii][jj].loci() [kk];
                    if ( !callBelongsInRecord(*call) )
                        continue;

                    SubFieldCall sfc(locus, call);
                    if (jj < gtCalls_[ii].size())
                        gtCalls_[ii][jj].push_back(sfc);
                    if (rec_.alleleIdx_[ii].back() < rec_.alleles_.size() &&
                        rec_.alleleIdx_[ii].back() != 0)
                        infoCalls_[rec_.alleleIdx_[ii].back()-1].push_back(sfc);
                }
                appendScores(hyp[ii][jj], ii);
            }
        }

        gtAnnData_.resize(gtCalls_.size());
        gtAnnFlag_.clear();
        gtAnnFlag_.resize(gtAnn_.size(), false);
        for(size_t gIdx=0; gIdx<gtCalls_.size(); gIdx++)
        {
            gtAnnData_[gIdx].resize(gtAnn_.size());
            for(size_t ii=0; ii<gtAnn_.size(); ii++)
            {
                gtAnnData_[gIdx][ii] = gtAnn_[ii]->getAnnotation(gIdx, gtCalls_[gIdx]);
                if ( (!isEmptyAnn(gtAnnData_[gIdx][ii])) ||
                     "GT" == gtAnn_[ii]->getId() ||
                     "PS" == gtAnn_[ii]->getId() )
                    gtAnnFlag_[ii] = true;
            }
        }
    }

    size_t VariantFileVcfRecordWriter::getSampleCount() const
    {
        return var_.size();
    }

    const char* VariantFileVcfRecordWriter::getSomaticStatusAnnotation(
        size_t gIdx,
        const std::vector< std::vector<SubFieldCall> >& calls) const
    {
        if (!hasSomatic_[gIdx])
            return ".";

        size_t otherIdx = var_.size();
        BOOST_FOREACH(size_t oId, otherIdx_[gIdx])
        {
            if (hasSomatic_[oId])
                continue;

            if (otherIdx != var_.size())
                return "."; // multiple possible baseline genomes
            otherIdx = oId;
        }
        if (var_.size() == otherIdx)
            return ".";

        for(size_t ii=0; ii<calls.size(); ii++)
        {
            BOOST_FOREACH(const SubFieldCall& sfc, calls[ii])
            {
                const Locus& locus = *sfc.locus_;
                if (locus.hasAnnotation("somaticScore") &&
                    locus.getAnnotation("somaticScore") != "")
                    return "Somatic";
            }
        }

        vector<size_t> aT(rec_.alleleIdx_[gIdx]);
        vector<size_t> aN(rec_.alleleIdx_[otherIdx]);

        bool hasNoCalls = false;
        for(size_t jj=0; jj<aT.size(); jj++)
        {
            if (aT[jj] >= rec_.alleles_.size())
            {
                hasNoCalls = true;
                break;
            }
        }
        for(size_t jj=0; jj<aN.size(); jj++)
        {
            if (aN[jj] >= rec_.alleles_.size())
            {
                hasNoCalls = true;
                break;
            }
        }

        if (hasNoCalls)
            return ".";

        std::sort(aT.begin(), aT.end());
        std::sort(aN.begin(), aN.end());

        if (aT == aN)
        {
            if (aT.size() > 0 && aT.back() != 0)
                return "Germline";
            return ".";
        }

        // If tumor is homozygous, normal is heterozygous, and one
        // allele is identical, return LOH.
        if ( 2 == aT.size() && 2 == aN.size() &&
             aN[0] != aN[1] && aT[0] == aT[1] &&
             ( aN[0] == aT[0] || aN[0] == aT[1] ||
               aN[1] == aT[0] || aN[1] == aT[1] ) )
            return "LOH";

        return ".";
    }

    void VariantFileVcfRecordWriter::addSubFieldHeaderRecords(
        std::vector<cgatools::conv::VcfSubFieldHeaderRecord>& result) const
    {
        for(size_t ii=0; ii<infoAnn_.size(); ii++)
            result.push_back(infoAnn_[ii]->getHeaderRecord());

        for(size_t ii=0; ii<gtAnn_.size(); ii++)
            result.push_back(gtAnn_[ii]->getHeaderRecord());
    }

    size_t VariantFileVcfRecordWriter::getAlleleCount(size_t alleleIdx) const
    {
        CGA_ASSERT(alleleIdx < rec_.alleles_.size());
        size_t count = 0;
        for(size_t gIdx=0; gIdx<rec_.alleleIdx_.size(); gIdx++)
        {
            for(size_t jj=0; jj<rec_.alleleIdx_[gIdx].size(); jj++)
            {
                if ( alleleIdx == rec_.alleleIdx_[gIdx][jj] )
                    count++;
            }
        }
        return count;
    }

    void VariantFileVcfRecordWriter::addAnnotation(std::vector< boost::shared_ptr<SubFieldAnnotation> >& anns,
                                                   SubFieldAnnotation* pAnn)
    {
        shared_ptr<SubFieldAnnotation> ann(pAnn);
        anns.push_back(ann);
    }

    std::string VariantFileVcfRecordWriter::getRefColumn() const
    {
        string result = crr_.getSequence(rec_.range_);
        CGA_ASSERT(0 != result.size());

        fixN(result);

        return result;
    }

    void VariantFileVcfRecordWriter::fixN(std::string& result) const
    {
        // VCF alleles and REF column must only contain A,C,G,T, or
        // N chars.
        // http://sourceforge.net/mailarchive/message.php?msg_id=27619958
        for(size_t ii=0; ii<result.size(); ii++)
        {
            char ch = result[ii];
            if ('A' != ch && 'C' != ch && 'G' != ch && 'T' != ch && 'N' != ch)
                result[ii] = 'N';
        }
    }

    bool VariantFileVcfRecordWriter::isEmptyAnn(const std::string& ann) const
    {
        return ann.size() ==
            size_t( std::count(ann.begin(), ann.end(), '.') +
                    std::count(ann.begin(), ann.end(), ',') );
    }

    uint32_t VariantFileVcfRecordWriter::swapIfNeeded(
        size_t gIdx,
        PhasedHypothesis& hyp,
        std::map< std::string, HaploEntry >& haploMap) const
    {
        CGA_ASSERT(2 == hyp.size());

        // Check haploMap to see if we need to be swapped.
        uint32_t ps = std::numeric_limits<uint32_t>::max();
        for(size_t jj=0; jj<hyp.size(); jj++)
        {
            const PhasedAllele& allele = hyp[jj];
            bool hapFound = false;
            BOOST_FOREACH(const Call* call, allele.calls())
            {
                if ( !callBelongsInRecord(*call) )
                    continue;

                if ( "" != call->hapLink_ &&
                     haploMap.find(call->hapLink_) != haploMap.end() )
                {
                    haploMap[call->hapLink_].lastSeenLocusCount_ = locusCount_;
                    ps = haploMap[call->hapLink_].pos_+1;
                    size_t prevAllele = haploMap[call->hapLink_].allele1_ ? 1 : 0;
                    if (prevAllele != jj)
                    {
//                             cerr << "swap "
//                                  << var_[gIdx]->getMetadata().get("ASSEMBLY_ID") << " "
//                                  << crr_.listChromosomes()[rec_.range_.chromosome_].getName()
//                                  << "," << (rec_.range_.begin_+1)
//                                  << endl;
                        std::swap(hyp[0], hyp[1]);
                    }
                    hapFound = true;
                    break;
                }
            }
            if (hapFound)
                break;
        }

        // Add all the call information to haploMap.
        for(size_t jj=0; jj<hyp.size(); jj++)
        {
            const PhasedAllele& allele = hyp[jj];
            BOOST_FOREACH(const Call* call, allele.calls())
            {
                if ( !callBelongsInRecord(*call) )
                    continue;

                if ( "" != call->hapLink_ &&
                     haploMap.find(call->hapLink_) == haploMap.end() )
                {
                    haploMap[call->hapLink_] = HaploEntry(rec_.range_.begin_, locusCount_, 0 != jj);
                    if (std::numeric_limits<uint32_t>::max() == ps)
                        ps = rec_.range_.begin_+1;
                }
            }
        }

        // Make sure haploMap doesn't grow too much.
        const size_t MAX_HAPLO_SIZE = 1000;
        if (haploMap.size() > MAX_HAPLO_SIZE)
        {
            vector<uint32_t> lastSeens;
            map< std::string, HaploEntry >::iterator first, last=haploMap.end();
            for(first=haploMap.begin(); first!=last; ++first)
            {
                lastSeens.push_back(first->second.lastSeenLocusCount_);
            }

            std::sort(lastSeens.begin(), lastSeens.end());

            uint32_t minLastSeen = lastSeens[MAX_HAPLO_SIZE / 2];
            for(first=haploMap.begin(); first!=last; )
            {
                if (first->second.lastSeenLocusCount_ < minLastSeen)
                {
                    map< std::string, HaploEntry >::iterator next = first;
                    ++next;
                    haploMap.erase(first);
                    first = next;
                }
                else
                {
                    ++first;
                }
            }
        }

        return ps;
    }

    bool VariantFileVcfRecordWriter::callBelongsInRecord(const Call& call) const
    {
        if (0 == call.range_.length() &&
            call.range_.beginLocation() != sl_->getRange().beginLocation())
        {
            // Inserts are added to record that includes preceding
            // reference base, except when they are at the beginning
            // of the superlocus, in which case they must be
            // attached to the existing record.
            if (call.range_.beginLocation() <= rec_.range_.beginLocation())
                return false;
        }
        if (rec_.range_.beginLocation() <= call.range_.beginLocation() &&
            call.range_.endLocation() <= rec_.range_.endLocation())
        {
            return true;
        }
        return false;
    }

    size_t VariantFileVcfRecordWriter::addVcfAllele(const std::string& allele)
    {
        if ("." == allele)
            return std::numeric_limits<size_t>::max();

        vector<string>::const_iterator iter =
            std::find(rec_.alleles_.begin(), rec_.alleles_.end(), allele);
        if (iter == rec_.alleles_.end())
        {
            rec_.alleles_.push_back(allele);
            return rec_.alleles_.size() - 1;
        }
        else
        {
            return iter - rec_.alleles_.begin();
        }
    }

    void VariantFileVcfRecordWriter::appendScores(const PhasedAllele& allele, size_t gIdx)
    {
        if (rec_.alleleIdx_[gIdx].back() >= rec_.alleles_.size())
        {
            rec_.hq_ [gIdx].push_back(0);
            rec_.ehq_[gIdx].push_back(0);
            rec_.cehq_[gIdx].push_back(0);
            return;
        }

        int32_t hq   = 0;
        int32_t ehq  = 0;
        double cehq = 0;

        for(size_t ii=0; ii<allele.calls().size(); ii++)
        {
            const Call* call   = allele.calls()[ii];
            const Locus* locus = allele.loci() [ii];
            if ( !callBelongsInRecord(*call) )
                continue;

            if (Call::EMPTY_SCORE != call->varScoreVAF_)
            {
                hq  = std::max(hq,  call->varScoreVAF_);
                ehq = std::max(ehq, call->varScoreEAF_);
                if (!calibPrefix_.empty())
                    cehq = std::max(cehq, getCalibratedScore(gIdx, locus, call, true));
            }
        }

        rec_.hq_  [gIdx].push_back(hq);
        rec_.ehq_ [gIdx].push_back(ehq);
        rec_.cehq_[gIdx].push_back(cehq);
    }

    double VariantFileVcfRecordWriter::getCalibratedScore(
        size_t gIdx, const Locus* locus, const Call* call, bool eaf)
    {
        CGA_ASSERT(!calibPrefix_.empty());
        CGA_ASSERT(eaf);

        if (0 == call->varScoreVAF_)
            return 0.0; // No score info.

        const string& totalReadCount = locus->getAnnotation("totalReadCount");
        if ("" == totalReadCount)
            return 0; // No coverage info.

        const char* scoreType = "fp";
        string varType = call->varType_;
        if (isHomozygousAlt(locus))
        {
            int32_t maxScore = std::numeric_limits<int32_t>::min();
            int32_t minScore = std::numeric_limits<int32_t>::max();

            BOOST_FOREACH(const Call& lc, locus->getCalls())
            {
                maxScore = std::max(maxScore, lc.varScoreEAF_);
                minScore = std::min(minScore, lc.varScoreEAF_);
            }

            CGA_ASSERT(minScore <= maxScore);

            if ( call->varScoreEAF_ > minScore || (maxScore == minScore && 1 == call->haplotype_) )
                scoreType = "fp";
            else
                scoreType = "oc";
        }
        else
        {
            if (call->isRefConsistent(crr_))
            {
                scoreType = "uc";
                varType = getLocusVarType(locus);
            }
            else
                scoreType = "fp";
        }
        varType = sanitizeVarType(varType);

        const string& softwareVersion = var_[gIdx]->getMetadata().getSoftwareVersionString();
        string calibrationKey = string(varType) + "-" + scoreType + "-" + softwareVersion;
        if (calibrations_.end() == calibrations_.find(calibrationKey))
        {
            calibrations_[calibrationKey].reset(
                new CalibratedScorer(varType, scoreType, eaf, calibPrefix_, softwareVersion));
        }
        int32_t cvg = parseValue<int32_t>(totalReadCount);
        return calibrations_[calibrationKey]->getCalibratedScore(cvg, call->varScoreEAF_);
    }

    bool VariantFileVcfRecordWriter::isHomozygousAlt(const Locus* locus) const
    {
        if (2 != locus->getPloidy())
            return false;

        if (2 != locus->getCalls().size())
            return false;

        if (locus->getCalls()[0].alleleSeq_ != locus->getCalls()[1].alleleSeq_)
            return false;

        if (locus->getCalls()[0].hasNoCalls() || locus->getCalls()[0].isRefConsistent(crr_))
            return false;

        return true;
    }

    std::string VariantFileVcfRecordWriter::getLocusVarType(const Locus* locus)
    {
        string varType;
        BOOST_FOREACH(const Call& lc, locus->getCalls())
        {
            if ( "snp" == lc.varType_ ||
                 "del" == lc.varType_ ||
                 "ins" == lc.varType_ ||
                 "sub" == lc.varType_ )
            {
                if ( (!varType.empty()) && varType != lc.varType_ )
                    return "sub";
                varType = lc.varType_;
            }
        }
        if ( varType.empty() )
            varType = "sub";
        return varType;
    }

    std::string VariantFileVcfRecordWriter::sanitizeVarType(const std::string& varType)
    {
        if ( "snp" == varType ||
             "del" == varType ||
             "ins" == varType ||
             "sub" == varType )
        {
            return varType;
        }
        return "sub";
    }

    std::string VariantFileVcfRecordWriter::getVcfAllele(const Superlocus& sl,
                                                         const PhasedAllele& allele)
    {
        string result;

        // Reference filler as necessary.
        if (rec_.range_.beginLocation() < sl.getRange().beginLocation())
        {
            result += crr_.getSequence(Range(rec_.range_.beginLocation(),
                                             sl.getRange().beginLocation()));
        }

        Range cRange(rec_.range_);
        cRange.begin_ = std::max(rec_.range_.begin_, sl.getRange().begin_);
        cRange.end_ = std::min(rec_.range_.end_, sl.getRange().end_);
        BOOST_FOREACH(const Call* call, allele.calls())
        {
            // Inserts must be combined with preceding base, which
            // is not possible for inserts before the first base.
            if (0 == call->range_.length() && 0 == call->range_.begin_)
                return ".";
            else if (callBelongsInRecord(*call) &&
                     cRange.beginLocation() <= call->range_.beginLocation() &&
                     call->range_.endLocation() <= cRange.endLocation())
            {
                result += call->calledSequence(crr_);
            }
            else if (cRange.beginLocation() < call->range_.endLocation() &&
                     call->range_.beginLocation() < cRange.endLocation())
            {
                // Even if record isn't fully contained by VCF
                // record but still overlaps VCF record, we need to
                // figure out the appropriate sequence for the ALT column.
                if (string::npos != call->alleleSeq_.find('?') ||
                    string::npos != call->alleleSeq_.find('N'))
                    result += "?";
                else
                {
                    CGA_ASSERT(call->isRefConsistent(crr_));
                    Range subRange(cRange);
                    subRange.begin_ = std::max(cRange.begin_, call->range_.begin_);
                    subRange.end_ = std::min(cRange.end_, call->range_.end_);
                    if ("?" != call->alleleSeq_ && "=" != call->alleleSeq_)
                    {
                        result += call->alleleSeq_.substr(
                            subRange.begin_-call->range_.begin_,
                            subRange.length());
                    }
                    else
                    {
                        result += crr_.getSequence(subRange);
                    }
                }
            }
        }

        // Reference filler as necessary.
        if (sl.getRange().endLocation() < rec_.range_.endLocation())
        {
            result += crr_.getSequence(Range(sl.getRange().endLocation(),
                                             rec_.range_.endLocation()));
        }

        if (string::npos != result.find('?'))
            return ".";
        fixN(result);
        return result;
    }

    //////////////////////////////////////////////////////////////////
    // VariantFileVcfRecordSource
    //////////////////////////////////////////////////////////////////

    VariantFileVcfRecordSource::VariantFileVcfRecordSource(
        const std::vector< boost::shared_ptr<VariantFileIterator> >& var,
        const std::vector<std::string> fieldNames,
        const CrrFile& crr,
        const std::string& calibPrefix,
        bool includeNoCalls)
        : crr_(crr),
          var_(var),
          slIt_(0, 0, 0),
          includeNoCalls_(includeNoCalls)
    {
        vector<VariantFileIterator*> varp;
        for(size_t ii=0; ii<var_.size(); ii++)
            varp.push_back(var_[ii].get());

        slIt_.setVariantFiles(varp);
        writer_.reset(new VariantFileVcfRecordWriter(var_, fieldNames, crr, calibPrefix));

        slIt_.seekFirst();
        initSuperlocus();

        if (!eof())
            writer_->setLocus(ranges_[rangesIdx_], *slIt_, hyp_);
    }

    std::vector<cgatools::conv::VcfSubFieldHeaderRecord>
    VariantFileVcfRecordSource::getSubFieldHeaderRecords() const
    {
        vector<VcfSubFieldHeaderRecord> result;
        writer_->addSubFieldHeaderRecords(result);
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "END", "1", "Integer",
                             "End position of the variant described in this record"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "CGA_NOCALL", "", "",
                             "No-called record"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FILTER, "VQLOW", "", "",
                             "Quality not VQHIGH"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FILTER, "SQLOW", "", "",
                             "Somatic quality not SQHIGH"));
        return result;
    }

    std::string
    VariantFileVcfRecordSource::getSource(size_t idxGenome) const
    {
        return var_[idxGenome]->getMetadata().getSoftwareVersionString();
    }

    std::vector<cgatools::conv::VcfKvHeaderRecord>
    VariantFileVcfRecordSource::getKeyValueHeaderRecords(size_t idxGenome) const
    {
        vector<VcfKvHeaderRecord> result;

        const char* transferKeys[] =
        {
            "GENOME_REFERENCE",
            "GENE_ANNOTATIONS",
            "DBSNP_BUILD",
            "COSMIC",
            "DGV_VERSION",
            "MIRBASE_VERSION",
            "PFAM_DATE",
            "REPMASK_GENERATED_AT",
            "SEGDUP_GENERATED_AT"
        };

        BOOST_FOREACH(const char* key, transferKeys)
        {
            if (var_[idxGenome]->getMetadata().hasKey(key))
                result.push_back(VcfKvHeaderRecord(string("source_")+key,
                                                   var_[idxGenome]->getMetadata().get(key)));
        }

        result.push_back(VcfKvHeaderRecord("phasing", "partial"));
        return result;
    }

    std::string
    VariantFileVcfRecordSource::getAssemblyId(size_t idxGenome) const
    {
        return var_[idxGenome]->getMetadata().get("ASSEMBLY_ID");
    }

    bool VariantFileVcfRecordSource::eof() const
    {
        return rangesIdx_ >= ranges_.size();
    }

    cgatools::conv::VcfRecordSource& VariantFileVcfRecordSource::operator++()
    {
        if (rangesIdx_+1 < ranges_.size())
        {
            rangesIdx_++;
        }
        else
        {
            ++slIt_;
            initSuperlocus();
        }

        if (!eof())
            writer_->setLocus(ranges_[rangesIdx_], *slIt_, hyp_);

        return *this;
    }

    const cgatools::conv::VcfRecordWriter& VariantFileVcfRecordSource::operator*() const
    {
        return *writer_;
    }

    const cgatools::conv::VcfRecordWriter* VariantFileVcfRecordSource::operator->() const
    {
        return writer_.get();
    }

    void VariantFileVcfRecordSource::initSuperlocus()
    {
        ranges_.clear();
        rangesIdx_ = 0;
        hyp_.clear();

        for (;;)
        {
            if (slIt_.eof())
                return;

            const Superlocus& sl = *slIt_;
            splitForVcf(sl, ranges_);
            phaseHaplotypes(sl, hyp_);

            if (0 != ranges_.size())
                return;

            ++slIt_;
        }
    }

    void VariantFileVcfRecordSource::splitForVcf(
        const cgatools::variants::Superlocus& sl,
        std::vector<cgatools::reference::Range>& ranges) const
    {
        ranges.clear();

        vector<uint32_t> insPositions;
        for(size_t ii=0; ii<sl.getGenomeCount(); ii++)
            addVariantRanges(ranges, sl.getLoci(ii), insPositions);
        for(size_t ii=0; ii<sl.getGenomeCount(); ii++)
            addNoCallPositions(ranges, sl.getLoci(ii), insPositions);

        if (0 != ranges.size())
        {
            std::sort(ranges.begin(), ranges.end());

            size_t to = 0;
            bool sticky = 0 == ranges[to].length();
            for(size_t from=1; from<ranges.size(); from++)
            {
                sticky |= std::find(insPositions.begin(), insPositions.end(), ranges[to].end_)
                    != insPositions.end();
                if (ranges[to].endLocation() < ranges[from].beginLocation())
                {
                    sticky = 0 == ranges[from].length();
                    ranges[++to] = ranges[from];
                    continue;
                }
                if ( ranges[to].endLocation() == ranges[from].beginLocation() &&
                     0 != ranges[from].length() &&
                     !sticky )
                {
                    sticky = false;
                    ranges[++to] = ranges[from];
                    continue;
                }
                CGA_ASSERT(ranges[to].chromosome_ == ranges[from].chromosome_);
                if (ranges[to].end_ < ranges[from].end_)
                {
                    sticky = 0 == ranges[from].length();
                    ranges[to].end_ = ranges[from].end_;
                }
            }

            ranges.erase(ranges.begin()+to+1, ranges.end());
        }

        if (includeNoCalls_)
        {
            // Add in the no-call ranges in between ranges. Split the
            // between-variant ranges by no-call boundaries.
            vector<uint32_t> ncSplits;
            for(size_t ii=0; ii<sl.getGenomeCount(); ii++)
                addNoCallSplits(ncSplits, sl.getLoci(ii));

            size_t origSize = ranges.size();
            for(size_t ii=0; ii<origSize; ii++)
            {
                // Split and add range before ranges[ii].
                Range rangeBefore(sl.getRange());
                rangeBefore.end_ = ranges[ii].begin_;
                if (0 != ii)
                    rangeBefore.begin_ = ranges[ii-1].end_;
                if (rangeBefore.begin_ < rangeBefore.end_)
                    splitAndAddRange(ranges, rangeBefore, ncSplits);
            }

            // Split and add range after last range in ranges.
            Range rangeAfter(sl.getRange());
            if (0 != origSize)
                rangeAfter.begin_ = ranges[origSize-1].end_;
            else if (0 == sl.getRange().length() && 0 != rangeAfter.begin_)
                rangeAfter.begin_--;
            if (rangeAfter.begin_ < rangeAfter.end_)
                splitAndAddRange(ranges, rangeAfter, ncSplits);

            std::sort(ranges.begin(), ranges.end());

            // Make sure we've covered all the bases of the superlocus,
            // and none of our ranges overlap.
            uint32_t begin = sl.getRange().begin_;
            for(size_t ii=0; ii<ranges.size(); ii++)
            {
                if (0 == ii)
                {
                    CGA_ASSERT(ranges[ii].begin_ <= begin);
                }
                else
                {
                    CGA_ASSERT(ranges[ii].begin_ == begin);
                }
                begin = ranges[ii].end_;
            }
            CGA_ASSERT(begin >= sl.getRange().end_);
        }
    }

    void VariantFileVcfRecordSource::addVariantRanges(
        std::vector<cgatools::reference::Range>& ranges,
        const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
        std::deque<cgatools::variants::Locus>::const_iterator>& loci,
        std::vector<uint32_t>& insPositions) const
    {
        BOOST_FOREACH(const Locus& locus, loci)
        {
            BOOST_FOREACH(const Call& call, locus.getCalls())
            {
                if ( (!call.isRefConsistent(crr_)) &&
                     string::npos == call.alleleSeq_.find('?') &&
                     string::npos == call.alleleSeq_.find('N') )
                {
                    // Need to add the variant's range. For VCF, indels
                    // and subs are described with the adjacent
                    // reference base.
                    Range range(call.range_);
                    if (call.range_.length() != call.alleleSeq_.size())
                    {
                        if (0 == call.reference_.size() || 0 == call.alleleSeq_.size() ||
                            call.reference_[0] != call.alleleSeq_[0] ||
                            isInsertionAtLeft(call.reference_, call.alleleSeq_))
                        {
                            if (isInsertionAtRight(call.reference_, call.alleleSeq_))
                                insPositions.push_back(call.range_.end_);
                            if (range.begin_ > 0)
                                range.begin_--;
                            else
                                range.end_++;
                        }
                    }
                    ranges.push_back(range);
                }
            }
        }
    }

    bool VariantFileVcfRecordSource::isInsertionAtLeft(
        const std::string& rr,
        const std::string& aa) const
    {
        if (rr.size() >= aa.size())
            return false;
        
        size_t sm;
        for(sm=0; sm<rr.size(); sm++)
        {
            if (rr[rr.size()-sm-1] != aa[aa.size()-sm-1])
                break;
        }
        if (sm == rr.size())
            return true;

        // If the input variant aligns with two bases equal to ref
        // that's consistent with insertion at the left, consider this
        // an insertion at left and prepend a base in the VCF output.
        size_t szDiff = aa.size()-rr.size();
        if (rr.size() > 2 && rr[0] == aa[szDiff] && rr[1] == aa[szDiff+1])
            return true;

        return false;
    }

    bool VariantFileVcfRecordSource::isInsertionAtRight(
        const std::string& rr,
        const std::string& aa) const
    {
        if (rr.size() >= aa.size())
            return false;
        
        size_t pm;
        for(pm=0; pm<rr.size(); pm++)
        {
            if (rr[pm] != aa[pm])
                break;
        }
        if (pm == rr.size())
            return true;

        // If the input variant aligns with two bases equal to ref
        // that's consistent with insertion at the right, consider this
        // an insertion at right.
        size_t szDiff = aa.size()-rr.size();
        if (rr.size() > 2 && rr[rr.size()-1] == aa[aa.size()-1-szDiff] &&
            rr[rr.size()-2] == aa[aa.size()-2-szDiff])
            return true;

        return false;
    }

    void VariantFileVcfRecordSource::addNoCallPositions(
        std::vector<cgatools::reference::Range>& ranges,
        const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
        std::deque<cgatools::variants::Locus>::const_iterator>& loci,
        const std::vector<uint32_t>& insPositions) const
    {
        if (0 == insPositions.size())
            return;

        BOOST_FOREACH(const Locus& locus, loci)
        {
            BOOST_FOREACH(const Call& call, locus.getCalls())
            {
                if (0 == call.range_.length())
                    continue;

                if ( string::npos == call.alleleSeq_.find('?') )
                    continue;

                BOOST_FOREACH(uint32_t pos, insPositions)
                {
                    if (call.range_.begin_ <= pos && pos < call.range_.end_)
                    {
                        uint32_t begin = pos;
                        if (0 != begin)
                            begin--;
                        ranges.push_back(Range(call.range_.chromosome_, begin, pos+1));
                    }
                }
            }
        }
    }

    void VariantFileVcfRecordSource::addNoCallSplits(
        std::vector<uint32_t>& ncSplits,
        const std::pair<std::deque<cgatools::variants::Locus>::const_iterator,
        std::deque<cgatools::variants::Locus>::const_iterator>& loci) const
    {
        BOOST_FOREACH(const Locus& locus, loci)
        {
            BOOST_FOREACH(const Call& call, locus.getCalls())
            {
                ncSplits.push_back(call.range_.begin_);
                ncSplits.push_back(call.range_.end_);

                // 0-length calls are attached to locus of previous base.
                if (0 == call.range_.length() && call.range_.begin_ > 0)
                    ncSplits.push_back(call.range_.begin_-1);
            }
        }

        // Sort and uniquify.
        std::sort(ncSplits.begin(), ncSplits.end());
        ncSplits.erase(std::unique(ncSplits.begin(), ncSplits.end()), ncSplits.end());
    }

    void VariantFileVcfRecordSource::splitAndAddRange(
        std::vector<cgatools::reference::Range>& ranges,
        const cgatools::reference::Range& range,
        const std::vector<uint32_t>& ncSplits) const
    {
        uint32_t begin = range.begin_;
        BOOST_FOREACH(uint32_t split, ncSplits)
        {
            if (begin < split && split < range.end_)
            {
                Range newRange(range);
                newRange.begin_ = begin;
                newRange.end_ = split;
                ranges.push_back(newRange);
                begin = split;
            }
        }
        if (begin < range.end_)
        {
            Range newRange(range);
            newRange.begin_ = begin;
            ranges.push_back(newRange);
        }
    }

    void VariantFileVcfRecordSource::phaseHaplotypes(
        const cgatools::variants::Superlocus& sl,
        std::vector<cgatools::variants::PhasedHypothesis>& hyp) const
    {
        size_t maxHypothesisCount = 32;
        vector< vector<PhasedHypothesis> > hypAll;
        sl.buildPhasedHypotheses(hypAll, maxHypothesisCount, true);
        hyp.clear();
        hyp.resize(hypAll.size());

        // Order hypAll by count ascending.
        vector<size_t> order;
        for(size_t ii=0; ii<hypAll.size(); ii++)
            order.push_back(ii);
        std::sort(order.begin(), order.end(), OrderBySizeAscending(hypAll));

        // For each genome, pick the phasing that maximizes count of
        // fully called haplotypes. If multiple best phasings exist by
        // that metric, pick the phasing the minimizes the count of
        // distinct fully called haplotypes in the population.
        set<string> popAlleles;
        for(size_t ii=0; ii<order.size(); ii++)
        {
            size_t gIdx = order[ii];
            vector<PhasedHypothesis>& gHyp = hypAll[gIdx];
            if (0 == gHyp.size())
            {
                // Multi-ploidy locus. Return empty hyp.
                continue;
            }

            size_t bestIdx = 0;
            int bestScore = scoreHypothesis(gHyp[0], popAlleles);
            for(size_t hIdx=1; hIdx<gHyp.size(); hIdx++)
            {
                int score = scoreHypothesis(gHyp[hIdx], popAlleles);
                if (score > bestScore)
                {
                    bestIdx = hIdx;
                    bestScore = score;
                }
            }

            hyp[gIdx] = gHyp[bestIdx];
            addPopAlleles(hyp[gIdx], popAlleles);
        }
    }

    int VariantFileVcfRecordSource::scoreHypothesis(
        const cgatools::variants::PhasedHypothesis& hyp,
        const std::set<std::string>& popAlleles) const
    {
        set<string> popInc;
        int score = 0;
        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            const PhasedAllele& pa = hyp[ii];
            const string& seq = pa.allele();
            if ( string::npos == seq.find('?') &&
                 string::npos == seq.find('N') )
            {
                score += 1000;
                if (popAlleles.find(seq) == popAlleles.end())
                    popInc.insert(seq);
            }
        }
        score -= popInc.size();
        return score;
    }

    void VariantFileVcfRecordSource::addPopAlleles(
        const cgatools::variants::PhasedHypothesis& hyp,
        std::set<std::string>& popAlleles) const
    {
        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            const PhasedAllele& pa = hyp[ii];
            const string& seq = pa.allele();
            if ( string::npos == seq.find('?') &&
                 string::npos == seq.find('N') )
                popAlleles.insert(seq);
        }
    }

} } // cgatools::variants
