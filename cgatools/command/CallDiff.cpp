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
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/command/CallDiff.hpp"
#include "cgatools/variants/SuperlocusIterator.hpp"
#include "cgatools/variants/calib/CalibratedScorer.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"
#include "cgatools/cgdata/ReferenceSupportReader.hpp"
#include "cgatools/util/parse.hpp"

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/algorithm/string.hpp>

#include <cmath>

namespace cgatools { namespace command {

    using namespace cgatools::util;
    using namespace cgatools::reference;
    using namespace cgatools::variants;
    using namespace cgatools::variants::calib;
    using namespace cgatools::cgdata;
    namespace bu = cgatools::util::baseutil;

    using namespace std;

    using boost::array;
    using boost::shared_ptr;
    namespace ba = boost::algorithm;

    // SQHIGH threshold.
    const double SOMATIC_QUAL_THRESH = -10.0;

    CallDiff::CallDiff(const std::string& name)
        : Command(name,
                  "Compares two Complete Genomics variant files.",
                  "0.3 or later",
                  "Compares two Complete Genomics variant files. Divides the genome "
                  "up into superloci of nearby variants, then compares the superloci. Also "
                  "refines the comparison to determine per-call or per-locus comparison "
                  "results.\n"
                  "\n"
                  "Comparison results are usually described by a semi-colon separated string, "
                  "one per allele. Each allele's comparison result is one of the following "
                  "classifications:\n\n"
                  "    ref-identical   \tThe alleles of the two variant files are identical, "
                  "and they are consistent with the reference.\n"
                  "    alt-identical   \tThe alleles of the two variant files are identical, "
                  "and they are inconsistent with the reference.\n"
                  "    ref-consistent  \tThe alleles of the two variant files are consistent, "
                  "and they are consistent with the reference.\n"
                  "    alt-consistent  \tThe alleles of the two variant files are consistent, "
                  "and they are inconsistent with the reference.\n"
                  "    onlyA           \tThe alleles of the two variant files are inconsistent, "
                  "and only file A is inconsistent with the reference.\n"
                  "    onlyB           \tThe alleles of the two variant files are inconsistent, "
                  "and only file B is inconsistent with the reference.\n"
                  "    mismatch        \tThe alleles of the two variant files are inconsistent, "
                  "and they are both inconsistent with the reference.\n"
                  "    phase-mismatch  \tThe two variant files would be consistent if the hapLink "
                  "field had been empty, but they are inconsistent.\n"
                  "    ploidy-mismatch \tThe superlocus did not have uniform ploidy.\n"
                  "\n"
                  "In some contexts, this classification is rolled up into a simplified "
                  "classification, which is one of \"identical\", \"consistent\", \"onlyA\", "
                  "\"onlyB\", or \"mismatch\".\n"
                  "\n"
                  "A good place to start looking at the results is the superlocus-output "
                  "file. It has columns defined as follows:\n\n"
                  "    SuperlocusId   \tAn identifier given to the superlocus.\n"
                  "    Chromosome     \tThe name of the chromosome.\n"
                  "    Begin          \tThe 0-based offset of the start of the superlocus.\n"
                  "    End            \tThe 0-based offset of the base one past the end of the superlocus.\n"
                  "    Classification \tThe match classification of the superlocus.\n"
                  "    Reference      \tThe reference sequence.\n"
                  "    AllelesA       \tA semicolon-separated list of the alleles (one per haplotype) "
                  "for variant file A, for the phasing with the best comparison result.\n"
                  "    AllelesB       \tA semicolon-separated list of the alleles (one per haplotype) "
                  "for variant file B, for the phasing with the best comparison result.\n"
                  "\n"
                  "The locus-output file contains, for each locus in file A and file B that is not "
                  "consistent with the reference, an annotated set of calls for the locus. The calls "
                  "are annotated with the following columns:\n\n"
                  "    SuperlocusId            \tThe id of the superlocus containing the locus.\n"
                  "    File                    \tThe variant file (A or B).\n"
                  "    LocusClassification     \tThe locus classification is determined by the "
                  "varType column of the call that is inconsistent with the reference, concatenated "
                  "with a modifier that describes whether the locus is heterozygous, homozygous, "
                  "or contains no-calls. If there is no one variant in the locus (i.e., it is "
                  "heterozygous alt-alt), the locus classification begins with \"other\".\n"
                  "    LocusDiffClassification \tThe match classification for the locus. This is "
                  "defined to be the best of the comparison of the locus to the same region in the "
                  "other file, or the comparison of the superlocus.\n\n"
                  "The somatic output file contains a list of putative somatic variations of genome A. "
                  "The output includes only those loci that can be classified as snp, del, ins or sub "
                  "in file A, and are called reference in the file B. "
                  "Every locus is annotated with the following columns:\n\n"
                  "    VarCvgA                 \tThe totalReadCount from file A for this locus (computed "
                  "on the fly if file A is not a masterVar file).\n"
                  "    VarScoreA               \tThe varScoreVAF from file A, or varScoreEAF if the "
                  "\"--diploid\" option is used.\n"
                  "    RefCvgB                 \tThe maximum of the uniqueSequenceCoverage values for "
                  "the locus in genome B.\n"
                  "    RefScoreB               \tMinimum of the reference scores of the locus in "
                  "genome B.\n"
                  "    SomaticCategory         \tThe category used for determining the calibrated scores "
                  "and the SomaticRank.\n"
                  "    VarScoreACalib          \tThe calibrated variant score of file A, under the "
                  "model selected by using or not using the \"--diploid\" option, and corrected for "
                  "the count of heterozygous variants observed in this genome. See user guide for more "
                  "information.\n"
                  "    VarScoreBCalib          \tThe calibrated reference score of file B, under the "
                  "model selected by using or not using the \"--diploid\" option, and corrected for "
                  "the count of heterozygous variants observed in this genome. See user guide for more "
                  "information.\n"
                  "    SomaticRank             \tThe estimated rank of this somatic mutation, amongst "
                  "all true somatic mutations within this SomaticCategory. The value is a number "
                  "between 0 and 1; a value of 0.012 means, for example, that an estimated 1.2% of "
                  "the true somatic mutations in this somaticCategory have a somaticScore less than "
                  "the somaticScore for this mutation. See user guide for more information.\n"
                  "    SomaticScore            \tAn integer that provides a total order on quality "
                  "for all somatic mutations. It is equal to -10*log10( P(false)/P(true) ), under "
                  "the assumption that this genome has a rate of somatic mutation equal to 1/Mb for "
                  "SomaticCategory snp, 1/10Mb for SomaticCategory ins, 1/10Mb for SomaticCategory "
                  "del, and 1/20Mb for SomaticCategory sub. The computation is based on the assumptions "
                  "described in the user guide, and is affected by choice of variant model selected by "
                  "using or not using the \"--diploid\" option.\n"
                  "    SomaticQuality          \tEqual to VQHIGH for all somatic mutations where "
                  "SomaticScore &gt;= -10. Otherwise, this column is empty.\n"
            ),
          extend3Mers_(4),
          extendBases_(0)
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The input crr file.")
            ("variantsA", po::value<string>(&variantFileNameA_),
             "The \"A\" input variant file.")
            ("variantsB", po::value<string>(&variantFileNameB_),
             "The \"B\" input variant file.")
            ("output-prefix", po::value<string>(&oPrefix_),
             "The path prefix for all output reports.")
            ("reports", po::value<string>(&reports_)->
             default_value("SuperlocusOutput,SuperlocusStats,LocusOutput,LocusStats"),
             "Comma-separated list of reports to generate. (Beware any reports whose "
             "name begins with \"Debug\".) A report is one of:\n"
             "    SuperlocusOutput      \tReport for superlocus classification.\n"
             "    SuperlocusStats       \tReport for superlocus classification stats.\n"
             "    LocusOutput           \tReport for locus classification.\n"
             "    LocusStats            \tReport for locus stats.\n"
             "    VariantOutput         \tBoth variant files annotated by comparison results."
             "If the somatic output report is requested, file A is also annotated with the "
             "same score ranks as produced in that report.\n"
             "    SomaticOutput         \tReport for the list of simple variations "
             "that are present only in file \"A\", annotated with the score that "
             "indicates the probability of the variation being truly somatic. "
             "Requires beta, genome-rootA, and genome-rootB options to be provided as well. "
             "Note: generating this report slows calldiff by 10x-20x.\n"
             "    DebugCallOutput       \tReport for call classification.\n"
             "    DebugSuperlocusOutput \tReport for debug superlocus information.\n"
             "    DebugSomaticOutput    \tReport for distribution estimates used for "
             "somatic rescoring. Only produced if SomaticOutput is also turned on.\n")
            ("diploid", po::bool_switch(&eaf_)->default_value(false),
             "Uses varScoreEAF instead of varScoreVAF in somatic score computations. "
             "Also, uses diploid variant model instead of variable allele mixture model.\n")
            ("locus-stats-column-count", po::value<size_t>(&statsColumnCount_)->default_value(15),
             "The number of columns for locus compare classification in the locus stats file.")
            ("max-hypothesis-count", po::value<size_t>(&maxHypothesisCount_)->default_value(32),
             "The maximum number of possible phasings to consider for a superlocus.")
            ("no-reference-cover-validation", po::bool_switch(&noReferenceCoverValidation_)
             ->default_value(false),
             "Turns off validation that all bases of a chromosome are covered by calls "
             "of the variant file.")
            ("genome-rootA", po::value<string>(&exportRootA_),
             "The \"A\" genome directory, for example /data/GS00118-DNA_A01; "
             "this directory is expected to contain ASM/REF and ASM/EVIDENCE subdirectories.")
            ("genome-rootB", po::value<string>(&exportRootB_),
             "The \"B\" genome directory.")
            ("calibration-root", po::value<string>(&calibPrefix_),
             "The directory containing calibration data. For example, there should "
             "exist a file calibration-root/0.0.0/metrics.tsv.")
            ("beta", "This flag enables the SomaticOutput report, which is beta functionality.");
            ;

        hiddenOptions_.add_options()
            ("somatic-input", po::value<string>(&somaticScoreDebugInput_)->default_value(""),
             "Prefix to previously created SomaticOutput and DebugSomaticOutput files. Used "
             "for debugging SomaticOutput report.")
            ;

        positionalOptions_.add("variantsA", 1);
        positionalOptions_.add("variantsB", 2);
    }

    typedef std::map< std::pair<int32_t,int32_t>, size_t> ScoreHist;

    struct VarData
    {
        VarData(uint32_t slid, const std::string& locusClass,
                int32_t varCvgA, int32_t varScoreA,
                int32_t refCvgB, int32_t refScoreB, const Call& call,
                const std::string& somaticCategory)
            :   superlocusId_(slid),
                somaticCategory_(somaticCategory),
                locusClass_(locusClass),
                cvgA_(varCvgA),
                varScoreA_(varScoreA),
                cvgB_(refCvgB),
                refScoreB_(refScoreB),
                varScoreACalib_(-1000),
                refScoreBCalib_(1000),
                varSensitivityA_(0),
                refSensitivityB_(0),
                somaticSensitivity_(0),
                somaticScore_(0),
                errorEstimate_(0),
                call_(call)
        {}

        uint32_t superlocusId_;
        std::string somaticCategory_;
        std::string locusClass_;
        int32_t cvgA_;
        int32_t varScoreA_;
        int32_t cvgB_;
        int32_t refScoreB_;

        double varScoreACalib_;
        double refScoreBCalib_;
        double varSensitivityA_;
        double refSensitivityB_;
        double somaticSensitivity_;
        double somaticScore_;
        double errorEstimate_;

        Call call_;

        bool operator<(const VarData& x) const
        {
            return call_.range_ < x.call_.range_;
        }
    };

    //! Class that maintains a mapping between score and cumulative
    //! count or fraction of some metric.
    class ScoreToCumulative
    {
    public:
        void addMapping(int score, double value)
        {
            data_[score] += value;
        }

        void finishInit(bool normalize)
        {
            if (normalize)
            {
                double cumulative = 0.0;
                map<int,double>::iterator first, last=data_.end();
                for(first=data_.begin(); first!=last; ++first)
                    cumulative += first->second;
                for(first=data_.begin(); first!=last; ++first)
                    first->second /= cumulative;
            }
            double cumulative = 0.0;
            map<int,double>::reverse_iterator first, last=data_.rend();
            for(first=data_.rbegin(); first!=last; ++first)
            {
                cumulative += first->second;
                first->second = cumulative;
            }
            if (normalize && data_.size() > 0)
            {
                // Sanity check the final value is close to 1.
                CGA_ASSERT(data_.begin()->second >  .999);
                CGA_ASSERT(data_.begin()->second < 1.001);
            }
        }

        // Returns cumulative sum of values for scores at or above score.
        double getCumulativeAtOrAbove(int score) const
        {
            map<int,double>::const_iterator iter = data_.lower_bound(score+1);
            if (iter == data_.end())
                return 0.0;
            return iter->second;
        }

        // Returns cumulative sum of values for scores below score.
        double getCumulativeBelow(int score) const
        {
            if (0 == data_.size())
                return 0.0;
            return data_.begin()->second - getCumulativeAtOrAbove(score);
        }

        void dump() const
        {
            cout << "DUMP" << endl;
            map<int,double>::const_iterator first, last=data_.end();
            for(first=data_.begin(); first!=last; ++first)
                cout << first->first << " " << first->second << endl;
            cout << "ENDDUMP" << endl;
        }

    private:
        std::map< int, double > data_;
    };

    class RecalibratedSomaticScoreCalc
    {
    public:
        class CalibratedScoreHist
        {
        public:
            CalibratedScoreHist(const CalibratedScorer& calib,
                                const ScoreHist& hist,
                                double trueCount,
                                double trueSomaticCount,
                                const CalibratedScoreHist* rs = 0)
                : calib_(calib),
                  rs_(rs),
                  trueCount_(trueCount),
                  trueSomaticCount_(trueSomaticCount)
            {
                counts_.resize(calib.getBinner().getBinCount(),
                               vector<uint32_t>(calib.getMaxScore()-calib.getMinScore()+1, 0));

                for(ScoreHist::const_iterator iter=hist.begin(); iter!=hist.end(); ++iter)
                    getCount(iter->first.first, iter->first.second) += iter->second;

                // Init calibScoreToSensitivity_.
                for(size_t ii=0; ii<counts_.size(); ii++)
                {
                    int32_t cvg = calib_.getBinner().getMinCvg(ii);
                    for(size_t jj=0; jj<counts_[ii].size(); jj++)
                    {
                        int32_t score = calib_.getMinScore() + jj;
                        double trueCount = getPTrue(cvg, score) * counts_[ii][jj];
                        calibScoreToSensitivity_.addMapping(
                            static_cast<int>(getS(cvg, score)), trueCount);
                    }
                }
                calibScoreToSensitivity_.finishInit(true);

                if (rs_)
                {
                    // Compute estimated count of true positive and true
                    // negative. In reality, these numbers are scaled
                    // versions of their selves. (True negatives are
                    // based on reference sampling, true positives are
                    // based on het variants.) These will be useful in
                    // the fraction computations below.
                    double trueVarCount = countTrue();
                    double trueRefCount = rs_->countTrue();

                    // Init somaticScoreToSensitivity_,
                    // somaticScoreToErrorFraction_.
                    for(size_t iiV=0; iiV<counts_.size(); iiV++)
                    {
                        int32_t varCvg = calib_.getBinner().getMinCvg(iiV);
                        for(size_t jjV=0; jjV<counts_[iiV].size(); jjV++)
                        {
                            int32_t varScore = calib_.getMinScore() + jjV;
                            double varFraction =
                                getPTrue(varCvg, varScore) * counts_[iiV][jjV] / trueVarCount;

                            for(size_t iiR=0; iiR<rs_->counts_.size(); iiR++)
                            {
                                int32_t refCvg = rs_->calib_.getBinner().getMinCvg(iiR);
                                for(size_t jjR=0; jjR<rs_->counts_[iiR].size(); jjR++)
                                {
                                    int32_t refScore = rs_->calib_.getMinScore() + jjR;
                                    double refFraction =
                                        rs_->getPTrue(refCvg, refScore) *
                                        rs_->counts_[iiR][jjR] / trueRefCount;

                                    double somScore =
                                        getSomaticScore(varCvg, varScore, refCvg, refScore);
                                    double fraction = varFraction * refFraction;

                                    somaticScoreToSensitivity_.addMapping(
                                        static_cast<int>(somScore), fraction);

                                    // Estimate false positive count for the set
                                    // of tumor variants with this
                                    // refScore,varScore combination.
                                    double fpEst =
                                        (1.0-getPTrue(varCvg, varScore)) * counts_[iiV][jjV] *
                                        refFraction;

                                    // Estimate false negative count for the set
                                    // of normal reference positions with this
                                    // refScore,varScore combination. Important
                                    // to use rs_->trueCount_ here, as it
                                    // reflects the count of bases in the
                                    // genome, rather than refCountTotal, which
                                    // is based on reference sampling.
                                    double fnEst =
                                        (1.0-rs_->getPTrue(refCvg, refScore)) * varFraction *
                                        rs_->counts_[iiR][jjR];

                                    somaticScoreToErrorFraction_.addMapping(
                                        static_cast<int>(somScore), fnEst+fpEst);
                                }
                            }
                        }
                    }
                    somaticScoreToSensitivity_   .finishInit(true);
                    somaticScoreToErrorFraction_ .finishInit(true);

//                     somaticScoreToErrorFraction_ .dump();
                }
            }

            uint32_t& getCount(int32_t cvg, int32_t score)
            {
                return counts_[calib_.getBinner().getBin(cvg)][scoreToOffset(score)];
            }

            const uint32_t& getCount(int32_t cvg, int32_t score) const
            {
                return counts_[calib_.getBinner().getBin(cvg)][scoreToOffset(score)];
            }

            double getS(int32_t cvg, int32_t score) const
            {
                return calib_.getCalibratedScore(cvg, score);
            }

            double getL(int32_t cvg, int32_t score) const
            {
                return calib_.getCalibratedLikelihoodRatio(cvg, score);
            }

            double getPTrue(int32_t cvg, int32_t score) const
            {
                return calib_.getPTrue(cvg, score);
            }

            double getSomaticScore(int32_t varCvgA, int32_t varScoreA,
                                   int32_t refCvgB, int32_t refScoreB) const
            {
                CGA_ASSERT(0 != rs_);

                double La = getL(varCvgA, varScoreA);
                double Lb = rs_->getL(refCvgB, refScoreB);
                double LSom = ( La *      trueCount_ / trueSomaticCount_ +
                                Lb * rs_->trueCount_ / trueSomaticCount_ );
                return std::floor(-10.0 * std::log10(LSom) + 0.5);
            }

            // Returns sensitivity of applying a calibrated score
            // threshold of "score", relative to sensitivity of not
            // applying a calibrated score threshold of "score".
            double calibScoreToSensitivity(double calibScore) const
            {
                return calibScoreToSensitivity_.getCumulativeAtOrAbove(
                    static_cast<int>(calibScore));
            }

            double getVarSensitivity(int32_t cvg, int32_t score) const
            {
                return calibScoreToSensitivity(getS(cvg, score));
            }

            double getSensitivity(int32_t varCvgA, int32_t varScoreA,
                                  int32_t refCvgB, int32_t refScoreB) const
            {
                double somScore = getSomaticScore(varCvgA, varScoreA, refCvgB, refScoreB);
                return somaticScoreToSensitivity_.getCumulativeAtOrAbove(
                    static_cast<int>(somScore));
            }

            double getErrorFraction(double somScore) const
            {
                return somaticScoreToErrorFraction_.getCumulativeAtOrAbove(
                    static_cast<int>(somScore));
            }

            double countTrue() const
            {
                double result = 0.0;
                for(size_t ii=0; ii<counts_.size(); ii++)
                {
                    int32_t varCvg = calib_.getBinner().getMinCvg(ii);
                    for(size_t jj=0; jj<counts_[ii].size(); jj++)
                    {
                        int32_t varScore = calib_.getMinScore() + jj;
                        result += getPTrue(varCvg, varScore) * counts_[ii][jj];
                    }
                }
                return result;
            }

        private:
            size_t scoreToOffset(int32_t score) const
            {
                score = std::max(score, calib_.getMinScore());
                score = std::min(score, calib_.getMaxScore());
                size_t offset = score-calib_.getMinScore();
                return offset;
            }

            CalibratedScorer calib_;
            const CalibratedScoreHist* rs_;
            std::vector< std::vector<uint32_t> > counts_;
            double trueCount_;
            double trueSomaticCount_;
            ScoreToCumulative calibScoreToSensitivity_;
            ScoreToCumulative somaticScoreToSensitivity_;
            ScoreToCumulative somaticScoreToErrorFraction_;
        };

        RecalibratedSomaticScoreCalc(const std::string& subTypeName,
                                     const std::string& calibPrefix,
                                     const std::string& softwareVersion,
                                     bool isEaf,
                                     const ScoreHist& rsHist,
                                     const ScoreHist& varHist,
                                     const CrrFile& crr)
        {
            // Estimate count of bases matching reference.
            uint64_t estimatedRefBases = 0;
            BOOST_FOREACH(const Range& range, crr.listContigs())
            {
                estimatedRefBases += range.length();
            }

            // Estimate count of variants.
            size_t estimatedTrueCount = 0;
            for(ScoreHist::const_iterator iter=varHist.begin(); iter!=varHist.end(); ++iter)
            {
                // If varScoreVAF is at least 30, consider it to be true.
                if (iter->first.second >= 30)
                    estimatedTrueCount += iter->second;
            }
            double basesPerVariant = double(estimatedRefBases) / double(estimatedTrueCount);

            double somCount;
            if ("snp" == subTypeName)
                somCount = estimatedRefBases / 1000000.0;
            else if ("ins" == subTypeName)
                somCount = estimatedRefBases / 10000000.0;
            else if ("del" == subTypeName)
                somCount = estimatedRefBases / 10000000.0;
            else if ("sub" == subTypeName)
                somCount = estimatedRefBases / 20000000.0;
            else
            {
                somCount = estimatedRefBases / 1000000.0;
                CGA_ASSERT(false);
            }

            CalibratedScorer refScorer(subTypeName, "fn", isEaf,
                                       calibPrefix, softwareVersion, 0.0, basesPerVariant);
            CalibratedScorer varScorer(subTypeName, "fp", isEaf,
                                       calibPrefix, softwareVersion,
                                       isEaf ? 0.0 : 0.5,
                                       basesPerVariant);
            rs_ .reset(new CalibratedScoreHist(
                           refScorer, rsHist,  estimatedRefBases, somCount));
            var_.reset(new CalibratedScoreHist(
                           varScorer, varHist, estimatedTrueCount, somCount, rs_.get()));
        }

        const CalibratedScoreHist& getRs() const
        {
            return *rs_;
        }

        const CalibratedScoreHist& getVar() const
        {
            return *var_;
        }

    private:
        boost::shared_ptr<CalibratedScoreHist> rs_;
        boost::shared_ptr<CalibratedScoreHist> var_;
    };

    class SomaticScoreCalc
    {
    public:
        SomaticScoreCalc(const std::string& rootA, const std::string& rootB,
                         const reference::CrrFile& crr, bool eaf,
                         const std::string& debugFileName,
                         const std::string& calibPrefix,
                         const std::string& softwareVersion,
                         const std::string& inputPrefix)
            : crr_(crr), scoreSrcA_(rootA, crr), scoreSrcB_(rootB, crr),
              eaf_(eaf), debugFileName_(debugFileName),
              calibPrefix_(calibPrefix),
              softwareVersion_(softwareVersion)
        {
            subtypes_.insert("snp");
            subtypes_.insert("del");
            subtypes_.insert("ins");
            subtypes_.insert("sub");

            snpSubtypes_.insert("snp");

            // Load score distributions from DebugSomaticOutput.csv.
            loadScoreDistributions(inputPrefix+"DebugSomaticOutput.tsv");

            // Load VarData from SomaticOutput.tsv.
            loadVarData(inputPrefix+"SomaticOutput.tsv", crr);

            validateCalibration();
        }

        SomaticScoreCalc(const std::string& rootA, const std::string& rootB,
                         const reference::CrrFile& crr, bool eaf,
                         const std::string& debugFileName,
                         const std::string& calibPrefix,
                         const std::string& softwareVersion)
            : crr_(crr), scoreSrcA_(rootA, crr), scoreSrcB_(rootB, crr),
              eaf_(eaf), debugFileName_(debugFileName),
              calibPrefix_(calibPrefix),
              softwareVersion_(softwareVersion)
        {
            subtypes_.insert("snp");
            subtypes_.insert("del");
            subtypes_.insert("ins");
            subtypes_.insert("sub");

            snpSubtypes_.insert("snp");

            validateCalibration();
        }

        void validateCalibration()
        {
            BOOST_FOREACH(const string& subtype, subtypes_)
            {
                CalibratedScorer refScorer(
                    subtype, "fn", eaf_, calibPrefix_, softwareVersion_, 0.0, 0.0);
                CalibratedScorer varScorer(
                    subtype, "fp", eaf_, calibPrefix_, softwareVersion_, eaf_ ? 0.0 : 0.5, 0.0);
            }
        }

        void loadScoreDistributions(const std::string& fileName)
        {
            InputStream in(fileName);
            string line;

            // Search for header.
            string header;
            while (std::getline(in, header))
            {
                if ( boost::ends_with(header, "-VarScoreA") ||
                     boost::ends_with(header, "-RefScoreB") )
                {
                    std::getline(in, line);
                    if (!boost::starts_with(line, "Score\tCoverage\tCount"))
                        continue;

                    ScoreHist hist;
                    while (std::getline(in, line))
                    {
                        if ("" == line)
                            break;

                        vector<string> fields;
                        boost::split(fields, line, boost::is_any_of("\t"));
                        int32_t score = boost::lexical_cast<int>(fields[0]);
                        int32_t coverage = boost::lexical_cast<int>(fields[1]);
                        int32_t count = boost::lexical_cast<int>(fields[2]);
                        hist[make_pair(score,coverage)] = count;
                    }

                    string locusClassId = header.substr(0, header.size()-10);

                    if (boost::ends_with(header, "-RefScoreB"))
                        refScores_ = hist;
                    else if (boost::ends_with(header, "-VarScoreA"))
                        variations_[locusClassId].varADistribution_ = hist;
                }
            }
        }

        void loadVarData(const std::string& fileName, const CrrFile& crr)
        {
            InputStream in(fileName);
            DelimitedFile df(in, fileName);

            Call call;
            uint32_t superlocusId;
            int32_t varScore, refScore, cvgA, cvgB;
            string locusClass, locusClassId;

            df.addField(ValueField<uint32_t>("SuperlocusId", &superlocusId));
            df.addField(StringField("LocusClassification", &locusClass));
            call.addFieldParsers(df, crr);
            df.addField(ValueField<int32_t>("VarScoreA", &varScore));
            df.addField(ValueField<int32_t>("RefScoreB", &refScore));
            df.addField(ValueField<int32_t>("VarCvgA", &cvgA));
            df.addField(ValueField<int32_t>("RefCvgB", &cvgB));
            df.addField(StringField("SomaticCategory", &locusClassId));

            while (df.next())
            {
                locusClassId = locusClassId.substr(locusClassId.size()-3);

                Subtype& st = variations_[locusClassId];

                st.items_.push_back(
                    VarData(superlocusId, locusClass, cvgA, varScore,
                            cvgB, refScore, call, locusClassId));
            }
        }

        void moveToSuperlocus(Range xrg, const SuperlocusIterator& slIt, bool haploid)
        {
            extend(xrg);

            Location refStart = slIt.getPrecedingRefStart();
            if (refStart < xrg.beginLocation() &&
                xrg.begin_ - refStart.offset_ >= MIN_REFERENCE_STRETCH &&
                refStart.chromosome_ == xrg.chromosome_)
            {
                Range allRg(refStart.chromosome_, refStart.offset_, xrg.begin_);
                scoreSrcB_.ref_.seek(allRg);
                for(uint32_t offset=refStart.offset_; offset<xrg.begin_; offset++)
                {
                    Range sampleRg(refStart.chromosome_, offset, offset+1);
                    int32_t cvg = scoreSrcB_.ref_.getMaxUniqueSequenceCoverage(sampleRg);
                    int32_t score = scoreSrcB_.ref_.getMinScore(sampleRg);
                    if (!haploid)
                        ++refScores_[make_pair(cvg,score)];
                }
            }

            scoreSrcA_.ref_.seek(xrg);
            scoreSrcB_.ref_.seek(xrg);
        }

        bool isVariantType(const std::string& locusClass,
                           const std::set<std::string>& types)
        {
            set<string>::const_iterator first, last=types.end();
            for(first=types.begin(); first!=last; ++first)
            {
                if (ba::ends_with(locusClass, *first))
                    return true;
            }
            return false;
        }

        int32_t getMaxAltVarScoreVAF(const variants::Locus& locus)
        {
            int32_t result = -1;
            BOOST_FOREACH(const variants::Call& call, locus.getCalls())
            {
                if ( ! ("snp" == call.varType_ ||
                        "ins" == call.varType_ ||
                        "del" == call.varType_ ||
                        "sub" == call.varType_) )
                    continue;

                if (int32_t(call.varScoreVAF_) > result)
                    result = call.varScoreVAF_;
            }
            return result;
        }

        int32_t getMaxAltVarScoreEAF(const variants::Locus& locus)
        {
            int32_t result = -1;
            BOOST_FOREACH(const variants::Call& call, locus.getCalls())
            {
                if ( ! ("snp" == call.varType_ ||
                        "ins" == call.varType_ ||
                        "del" == call.varType_ ||
                        "sub" == call.varType_) )
                    continue;

                if (int32_t(call.varScoreEAF_) > result)
                    result = call.varScoreEAF_;
            }
            return result;
        }

        void processLocus(uint32_t superlocusId,
                          const variants::Locus& locus,
                          const variants::VariantFileIterator& varFile,
                          const std::string& locusClass,
                          const vector<cdmt::MatchType>& mt,
                          size_t refBAlleleCount)
        {
            if (!isVariantType(locusClass, subtypes_))
                return;

            // Skip chromosome M: it has unusual score distribution and too little
            // data to empirically assess it
            uint16_t chrId = locus.getRange().chromosome_;
            if (crr_.listChromosomes()[chrId].getName() == "chrM")
            {
                return;
            }

            Range xrg = locus.getRange();
            extend(xrg);

            int32_t varScore = eaf_ ? getMaxAltVarScoreEAF(locus) : getMaxAltVarScoreVAF(locus);
            int32_t allele1ReadCount, allele2ReadCount, referenceAlleleReadCount, varCvg;
            varFile.getReadCounts(&allele1ReadCount,
                                  &allele2ReadCount,
                                  &referenceAlleleReadCount,
                                  &varCvg,
                                  scoreSrcA_.evidence_,
                                  locus);

            string locusClassId = locusClass.substr(locusClass.size()-3);
            Subtype& st = variations_[locusClassId];

            if (boost::starts_with(locusClass, "het-"))
            {
                // Include this locus in *ScoreA distributions.
                ++st.varADistribution_[make_pair(varCvg, varScore)];
            }

            // Ignore this locus if it's not somatic.
            if (std::find(mt.begin(), mt.end(), cdmt::ONLY_A) == mt.end())
                return;
            if (std::count(mt.begin(), mt.end(), cdmt::ONLY_A) +
                std::count(mt.begin(), mt.end(), cdmt::REF_IDENTICAL) !=
                locus.getPloidy() &&
                refBAlleleCount != locus.getPloidy())
                return;

            int32_t refcvg   = scoreSrcB_.ref_.getMaxUniqueSequenceCoverage(xrg);
            int32_t refscore = scoreSrcB_.ref_.getMinScore(xrg);

            // Since we only consider simple loci, the very first non-reference call
            // is the one we are interested in
            const Call* call = 0;
            BOOST_FOREACH(const Call& c, locus.getCalls())
            {
                if (subtypes_.count(c.varType_) != 0)
                {
                    call = &c;
                    break;
                }
            }
            CGA_ASSERT(0 != call);
            st.items_.push_back(
                VarData(superlocusId, locusClass, varCvg, varScore,
                        refcvg, refscore, *call, locusClassId));
        }

        void dump(std::ostream& out)
        {

            out << "SuperlocusId\tLocusClassification\t"
                << Call::getHeader() <<
                "\tVarCvgA\tVarScoreA\tRefCvgB\tRefScoreB\t"
                "SomaticCategory\t"
                "VarScoreACalib\t"
                "RefScoreBCalib\t"
                "SomaticRank\t"
                "SomaticScore\t"
                "SomaticQuality\n";

            const std::vector<VarData>& final = getFinalData();
            BOOST_FOREACH(const VarData& v, final)
            {
                out << v.superlocusId_ << "\t";
                out << v.locusClass_ << "\t";
                v.call_.write(out, crr_);
                out << "\t" << v.cvgA_
                    << "\t" << v.varScoreA_
                    << "\t" << v.cvgB_
                    << "\t" << v.refScoreB_
                    << "\t" << v.somaticCategory_
                    << "\t" << boost::format("%.0f") % v.varScoreACalib_
                    << "\t" << boost::format("%.0f") % v.refScoreBCalib_
                    << "\t" << boost::format("%.3f") % (1.0-v.somaticSensitivity_)
                    << "\t" << boost::format("%.0f") % v.somaticScore_
                    << "\t" << ( v.somaticScore_ >= SOMATIC_QUAL_THRESH ? "SQHIGH" : "" )
                    << "\n";
            }
        }

        class RescoreOrder
        {
        public:
            bool operator()(const std::string& lhs, const std::string& rhs) const
            {
                if (lhs == rhs)
                    return false;
                if ("snp" == lhs)
                    return true;
                if ("snp" == rhs)
                    return false;
                if ("ins" == lhs)
                    return true;
                if ("ins" == rhs)
                    return false;
                return lhs < rhs;
            }
        };

        //! Retrieves all scored somatic variations, must be called
        //! after processing all superloci.
        const std::vector<VarData>& getFinalData()
        {
            if (finalData_.size() != 0)
                return finalData_;

            shared_ptr<ostream> dbg;
            if (debugFileName_ != "")
                dbg = OutputStream::openCompressedOutputStreamByExtension(debugFileName_);

            vector<string> keys;
            for( SomaticVariations::const_iterator first=variations_.begin();
                 first!=variations_.end();
                 ++first)
                keys.push_back(first->first);
            std::sort(keys.begin(), keys.end(), RescoreOrder());
            BOOST_FOREACH(string& key, keys)
            {
                rescoreSubtype(key, variations_[key], finalData_, dbg);
            }
            std::sort(finalData_.begin(), finalData_.end());

            return finalData_;
        }

    private:
        static const uint32_t MIN_REFERENCE_STRETCH = 1;

        struct ScoreSource
        {
            std::string root_;
            cgdata::GenomeMetadata exp_;
            cgdata::EvidenceReader evidence_;
            cgdata::ReferenceSupportReader ref_;

            ScoreSource(const std::string& root, const reference::CrrFile& crr)
                : root_(root), exp_(root), evidence_(crr, exp_), ref_(crr, exp_)
            {}
        };

        struct Subtype
        {
            ScoreHist varADistribution_;
            std::vector<VarData> items_;
        };

        typedef std::map<std::string, Subtype> SomaticVariations;

        const reference::CrrFile& crr_;
        ScoreSource scoreSrcA_, scoreSrcB_;
        bool eaf_;
        std::string debugFileName_;
        std::string calibPrefix_;
        std::string softwareVersion_;
        int32_t refScoreB_, evidenceScore_;

        std::set<std::string> subtypes_, snpSubtypes_;
        ScoreHist refScores_;
        SomaticVariations variations_;
        std::vector<VarData> finalData_;

        void extend(Range& xrg)
        {
            if (xrg.begin_ > 0)
                --xrg.begin_;
            ++xrg.end_;
        }

        void rescoreSubtype(
            std::string subTypeName,
            const Subtype& st,
            std::vector<VarData>& out,
            boost::shared_ptr<ostream>& dbg) const
        {
            if (st.items_.empty())
                return;

            if (dbg)
            {
                dumpCounts(st.varADistribution_, subTypeName + "-VarScoreA", *dbg);
                dumpCounts(refScores_, subTypeName + "-RefScoreB", *dbg);
            }

            RecalibratedSomaticScoreCalc calc(
                subTypeName, calibPrefix_, softwareVersion_, eaf_, refScores_,
                st.varADistribution_, crr_);
            size_t startingCount = out.size();
            BOOST_FOREACH(const VarData& v, st.items_)
            {
                out.push_back(v);

                VarData& o = out.back();
                const RecalibratedSomaticScoreCalc::CalibratedScoreHist& vsh = calc.getVar();
                const RecalibratedSomaticScoreCalc::CalibratedScoreHist& rsh = calc.getRs();

                o.varSensitivityA_    = vsh.getVarSensitivity(v.cvgA_, v.varScoreA_);
                o.refSensitivityB_    = rsh.getVarSensitivity(v.cvgB_, v.refScoreB_);
                o.somaticSensitivity_ = vsh.getSensitivity(v.cvgA_, v.varScoreA_,
                                                           v.cvgB_, v.refScoreB_);

                o.varScoreACalib_ = vsh.getS(v.cvgA_, v.varScoreA_);
                o.refScoreBCalib_ = rsh.getS(v.cvgB_, v.refScoreB_);
                o.somaticScore_  = vsh.getSomaticScore(v.cvgA_, v.varScoreA_,
                                                       v.cvgB_, v.refScoreB_);
                o.errorEstimate_ = vsh.getErrorFraction(o.somaticScore_);
            }

            if (dbg)
            {
                dumpReplicateRoc(subTypeName + "-ReplicateRoc",
                                 calc,
                                 out.begin()+startingCount,
                                 out.end(),
                                 *dbg);
            }
        }

        void dumpCounts(
            const ScoreHist& hist, const std::string& title, std::ostream& out) const
        {
            out << title << "\n";
            out << "Score\tCoverage\tCount\n";
            ScoreHist::const_iterator first, last=hist.end();
            for(first=hist.begin(); first!=last; ++first)
            {
                out << first->first.first << "\t"
                    << first->first.second << "\t"
                    << first->second << "\n";
            }
            out << endl;
        }

        void dumpReplicateRoc(
            const std::string& title,
            const RecalibratedSomaticScoreCalc& calc,
            std::vector<VarData>::const_iterator vFirst,
            std::vector<VarData>::const_iterator vLast,
            std::ostream& out
            ) const
        {
            out << title << "\n";
            out << "SomaticScore\tDiscordantCount\tScaledSensitivity\n";

            int MAX_SSCORE = 100;
            double tcount = calc.getVar().countTrue();
            map<int, uint32_t> disc, sens;

            for(vector<VarData>::const_iterator vIter=vFirst; vIter!=vLast; ++vIter)
            {
                const VarData& v = *vIter;
                int score = static_cast<int>(v.somaticScore_);
                score = std::min(std::max(score, -MAX_SSCORE), MAX_SSCORE);
                sens[score] = static_cast<uint32_t>(tcount * v.somaticSensitivity_);
                disc[score]++;
            }

            for(int score=MAX_SSCORE; score>=-MAX_SSCORE; score--)
            {
                sens[score] = std::max(sens[score], sens[score+1]);
                disc[score] += disc[score+1];
            }

            for(int score=-MAX_SSCORE; score<=MAX_SSCORE; score++)
            {
                out << score
                    << "\t" << disc[score]
                    << "\t" << sens[score]
                    << endl;
            }
            out << endl;
        }
    };

    class MafReport
    {
    public:
        enum VariantClassification
        {
            VC_SILENT,
            VC_UNKNOWN,
            VC_UPSTREAM,
            VC_UTR,
            VC_INTRON,
            VC_SPLICE_SITE,
            VC_MISSENSE,
            VC_IN_FRAME,
            VC_NONSTOP,
            VC_FRAME_SHIFT,
            VC_MISSTART,
            VC_NONSENSE
        };

        struct MafRecord
        {
            MafRecord()
                : varScoreARank_(-1.0),
                  refScoreBRank_(-1.0)
            {
                totalReadCounts_.assign(-1);
                altReadCounts_.assign(-1);
            }

            std::string getVariantType() const
            {
                const char* result = "Unknown";
                result = getVariantType(result, alleles_[0][0]);
                result = getVariantType(result, alleles_[0][1]);
                result = getVariantType(result, alleles_[1][0]);
                result = getVariantType(result, alleles_[1][1]);
                return result;
            }

            const char* getMutationStatus() const
            {
                if (!bu::isCalledSequence(alleles_[1][0])) return "Unknown";
                if (!bu::isCalledSequence(alleles_[1][1])) return "Unknown";

                if ( alleles_[1][0] == refAllele_ && alleles_[1][1] == refAllele_)
                {
                    if ( bu::isCalledSequence(alleles_[0][0]) && alleles_[0][0] != refAllele_)
                        return "Somatic";
                    if ( bu::isCalledSequence(alleles_[0][1]) && alleles_[0][1] != refAllele_)
                        return "Somatic";
                }

                if (!bu::isCalledSequence(alleles_[0][0])) return "Unknown";
                if (!bu::isCalledSequence(alleles_[0][1])) return "Unknown";

                if (alleles_[0][0] == alleles_[1][0] &&
                    alleles_[0][1] == alleles_[1][1])
                    return "Germline";
                if ( alleles_[0][0] == alleles_[0][1] &&
                     alleles_[1][0] != alleles_[1][1] &&
                     ( alleles_[0][0] == alleles_[1][0] || alleles_[0][1] == alleles_[1][1] ) )
                    return "LOH";
                return "Unknown";
            }

            std::vector<std::string> hugoSymbol_;
            std::vector<std::string> entrezGeneId_;
            reference::Range range_; //
            std::vector<VariantClassification> variantClassification_;
            std::string refAllele_;
            // alleles_[0] <- tumor, alleles_[1] <- normal
            boost::array< boost::array<std::string, 2>, 2 > alleles_;
            boost::array< int, 2 > totalReadCounts_;
            boost::array< int, 2 > altReadCounts_;
            std::string dbSnpRs_;
            std::vector<std::string> pfamDomain_;
            std::vector<std::string> miRna_;
            double varScoreARank_, refScoreBRank_;

        private:
            const char* getVariantType(const char* orig, const std::string& allele) const
            {
                if (allele == refAllele_)
                    return orig;
                if (!bu::isCalledSequence(allele))
                {
                    if ( 0 == refAllele_.size() && ba::starts_with(orig, "Unknown") )
                        return "Unknown-Ins";
                    return orig;
                }
                if (bu::isConsistent(refAllele_, allele))
                    return orig;
                const char* vt = "Sub";
                if (allele.size() == 1 && refAllele_.size() == 1)
                    vt = "SNP";
                else if (0 == allele.size() && 0 != refAllele_.size())
                    vt = "Del";
                else if (0 == refAllele_.size() && 0 != allele.size())
                    vt = "Ins";
                if ( (!ba::starts_with(orig, "Unknown")) && 0 != strcmp(orig, vt) )
                    return "Unknown";
                return vt;
            }
        };

        MafReport(const std::string& fn,
                  const std::string& exportRootA,
                  const std::string& exportRootB,
                  SomaticScoreCalc& somatic,
                  const CrrFile& crr)
            : fn_(fn),
              exportRootA_(exportRootA),
              exportRootB_(exportRootB),
              metaA_(exportRootA),
              metaB_(exportRootB),
              somatic_(somatic),
              crr_(crr),
              gscCenter_("COMPLETE_GENOMICS"),
              ncbiBuild_(metaA_.getGenomeReference()),
              tumorBarCode_("TUMOR_BAR_CODE"),
              normalBarCode_("NORMAL_BAR_CODE")
        {
            evReader_[0].reset(new EvidenceReader(crr, metaA_));
            evReader_[1].reset(new EvidenceReader(crr, metaB_));
        }

        void addMafRecords(const CallDiffResult& dr, const Superlocus& sl)
        {
            array<PhasedHypothesis, 2> hyp;
            if (1 == dr.hyp_[0].size())
            {
                hyp[0] = PhasedHypothesis(dr.hyp_[0].getRange(), 2);
                hyp[0][0] = dr.hyp_[0][0];
                hyp[0][1] = dr.hyp_[0][0];

                hyp[1] = PhasedHypothesis(dr.hyp_[1].getRange(), 2);
                hyp[1][0] = dr.hyp_[1][0];
                hyp[1][1] = dr.hyp_[1][0];
            }
            else if (2 == dr.hyp_[0].size())
            {
                hyp = dr.hyp_;
            }
            else
                return; // not haploid or diploid -- give up

            // If call-specific comparison is compatible with superlocus
            // comparison, break superlocus up by ref-inconsistent calls
            // of both files.
            if (canSplitByCallsOfBothFiles(dr))
            {
                vector<Range> ranges;
                addVariantRanges(ranges, sl.getLoci(0));
                if (dr.hyp_[0].size() > 1)
                    addVariantRanges(ranges, sl.getLoci(1));
                orderRangesForMaf(ranges);
                addMafRecords(ranges, hyp);
            }

            // If identical, break superlocus up by ref-inconsistent
            // calls of file A.
            else if (allAllelesAreIdentical(dr))
            {
                vector<Range> ranges;
                addVariantRanges(ranges, sl.getLoci(0));
                orderRangesForMaf(ranges);
                hyp[1] = hyp[0];
                addMafRecords(ranges, hyp);
            }

            // Use superlocus.
            else
            {
                recs_.push_back(MafRecord());

                MafRecord& rec = recs_.back();
                rec.range_ = sl.getRange();
                rec.alleles_[0][0] = hyp[0][0].allele();
                rec.alleles_[0][1] = hyp[0][1].allele();
                rec.alleles_[1][0] = hyp[1][0].allele();
                rec.alleles_[1][1] = hyp[1][1].allele();

                addAnnotation(rec, hyp);
            }
        }

        void finish()
        {
            annotateByGeneFile(exportRootA_);
            annotateByGeneFile(exportRootB_);
            annotateByMirna(exportRootA_);
            annotateByMirna(exportRootB_);
            annotateBySomaticScores();

            shared_ptr<ostream> out = OutputStream::openCompressedOutputStreamByExtension(fn_);
            writeHeader(*out);
            BOOST_FOREACH(const MafRecord& rec, recs_)
            {
                writeRecord(*out, rec);
            }
        }

    private:
        class MafRecordLtRange
        {
        public:
            bool operator()(const MafRecord& rec, const Range& range) const
            {
                return rec.range_ < range;
            }

            bool operator()(const MafRecord& rec0, const MafRecord& rec1) const
            {
                return rec0.range_ < rec1.range_;
            }

            bool operator()(const Range& range, const MafRecord& rec) const
            {
                return range < rec.range_;
            }
        };

        void annotateByGeneFile(const std::string exportRoot)
        {
            cgdata::GenomeMetadata meta(exportRoot);
            string fn = meta.getGeneFileName();
            boost::shared_ptr<istream> in = InputStream::openCompressedInputStreamByExtension(fn);
            DelimitedFile df(*in, fn);

            Range gRange;
            string geneId, symbol, component, impact, pfam;

            df.addField(ChromosomeIdField("chromosome", &gRange.chromosome_, crr_));
            df.addField(ValueField<uint32_t>("begin", &gRange.begin_));
            df.addField(ValueField<uint32_t>("end", &gRange.end_));
            df.addField(StringField("geneId", &geneId));
            df.addField(StringField("symbol", &symbol));
            df.addField(StringField("component", &component));
            df.addField(StringField("impact", &impact));
            df.addField(StringField("pfam", &pfam));

            while (df.next())
            {
                vector<MafRecord>::iterator iter =
                    std::lower_bound(recs_.begin(), recs_.end(), gRange, MafRecordLtRange());
                if (iter != recs_.begin())
                    --iter;
                for(; iter != recs_.end() && iter->range_.beginLocation() <= gRange.endLocation(); ++iter)
                {
                    MafRecord& rec = *iter;
                    if (gRange.beginLocation() < rec.range_.beginLocation())
                        continue;
                    if (rec.range_.endLocation() < gRange.endLocation())
                        continue;

                    // Apply gene record to MafRecord.
                    if (std::find(rec.entrezGeneId_.begin(), rec.entrezGeneId_.end(), geneId) ==
                        rec.entrezGeneId_.end())
                    {
                        rec.entrezGeneId_.push_back(geneId);
                        rec.hugoSymbol_.push_back(symbol);
                        rec.variantClassification_.push_back(getClassification(component, impact));
                        rec.pfamDomain_.push_back(pfam);
                    }
                    else
                    {
                        vector<string>::const_iterator itGid =
                            std::find(rec.entrezGeneId_.begin(), rec.entrezGeneId_.end(), geneId);
                        CGA_ASSERT(itGid != rec.entrezGeneId_.end());
                        size_t ii = itGid - rec.entrezGeneId_.begin();
                        rec.variantClassification_[ii] =
                            mergeClassification(getClassification(component, impact),
                                                rec.variantClassification_[ii]);
                        if (0 == rec.pfamDomain_[ii].size())
                            rec.pfamDomain_[ii] = pfam;
                    }
                }
            }
        }

        void annotateByMirna(const std::string exportRoot)
        {
            cgdata::GenomeMetadata meta(exportRoot);
            string fn = meta.getNcRNAFileName();
            boost::shared_ptr<istream> in = InputStream::openCompressedInputStreamByExtension(fn);
            DelimitedFile df(*in, fn);

            Range mRange;
            string miRBaseId;

            df.addField(ChromosomeIdField("chr", &mRange.chromosome_, crr_));
            df.addField(ValueField<uint32_t>("begin", &mRange.begin_));
            df.addField(ValueField<uint32_t>("end", &mRange.end_));
            df.addField(StringField("miRBaseId", &miRBaseId));

            while (df.next())
            {
                vector<string> ids;
                ba::split(ids, miRBaseId, ba::is_any_of(";"));

                vector<MafRecord>::iterator iter =
                    std::lower_bound(recs_.begin(), recs_.end(), mRange, MafRecordLtRange());
                if (iter != recs_.begin())
                    --iter;
                for(; iter != recs_.end() && iter->range_.beginLocation() <= mRange.endLocation(); ++iter)
                {
                    MafRecord& rec = *iter;
                    if (mRange.beginLocation() < rec.range_.beginLocation())
                        continue;
                    if (rec.range_.endLocation() < mRange.endLocation())
                        continue;

                    // Apply ncRNA record to MafRecord.
                    BOOST_FOREACH(const string& id, ids)
                    {
                        if (std::find(rec.miRna_.begin(), rec.miRna_.end(), id) ==
                            rec.miRna_.end())
                        {
                            rec.miRna_.push_back(id);
                        }
                    }
                }
            }
        }

        void annotateBySomaticScores()
        {
            const vector<VarData>& final = somatic_.getFinalData();
            BOOST_FOREACH(const VarData& vv, final)
            {
                const Range& vRange = vv.call_.range_;

                vector<MafRecord>::iterator iter =
                    std::lower_bound(recs_.begin(), recs_.end(), vRange, MafRecordLtRange());
                if (iter != recs_.begin())
                    --iter;
                for(; iter != recs_.end() && iter->range_.beginLocation() <= vRange.endLocation(); ++iter)
                {
                    MafRecord& rec = *iter;
                    if (vRange.beginLocation() < rec.range_.beginLocation())
                        continue;
                    if (rec.range_.endLocation() < vRange.endLocation())
                        continue;

                    // Apply somatic score to MafRecord.
                    if (rec.refScoreBRank_ < 0.0 ||
                        std::min(rec.varScoreARank_, rec.refScoreBRank_) <
                        std::min(1.0-vv.varSensitivityA_, 1.0-vv.refSensitivityB_))
                    {
                        rec.varScoreARank_ = 1.0 - vv.varSensitivityA_;
                        rec.refScoreBRank_ = 1.0 - vv.refSensitivityB_;
                    }
                }
            }
        }

        VariantClassification getClassification(const std::string& component,
                                                const std::string& impact)
        {
            if ("TSS-UPSTREAM" == component)
                return VC_UPSTREAM;
            if ("INTRON" == component)
                return VC_INTRON;
            if (ba::starts_with(component, "UTR"))
                return VC_UTR;
            if ("DONOR" == component || "ACCEPTOR" == component || ba::starts_with(component, "SPAN"))
                return VC_SPLICE_SITE;
            if ("NO-CHANGE" == impact || "COMPATIBLE" == impact || "SYNONYMOUS" == impact)
                return VC_SILENT;
            if ("MISSENSE" == impact)
                return VC_MISSENSE;
            if ("NONSTOP" == impact)
                return VC_NONSTOP;
            if ("MISSTART" == impact)
                return VC_MISSTART;
            if ("NONSENSE" == impact)
                return VC_NONSENSE;
            if ("DELETE" == impact || "DELETE+" == impact || "INSERT" == impact || "INSERT+" == impact)
                return VC_IN_FRAME;
            if ("FRAMESHIFT" == impact)
                return VC_FRAME_SHIFT;
            return VC_UNKNOWN;
        }

        VariantClassification mergeClassification(VariantClassification vc1, VariantClassification vc2)
        {
            return std::max(vc1, vc2);
        }

        const char* vcToString(VariantClassification vc) const
        {
            switch(vc)
            {
            case VC_SILENT:
                return "Silent";
            case VC_UNKNOWN:
                return "Unknown";
            case VC_UPSTREAM:
                return "Tss_Upstream";
            case VC_UTR:
                return "Utr";
            case VC_INTRON:
                return "Intron";
            case VC_SPLICE_SITE:
                return "Splice_Site";
            case VC_MISSENSE:
                return "Missense";
            case VC_IN_FRAME:
                return "In_Frame";
            case VC_NONSTOP:
                return "Nonstop";
            case VC_FRAME_SHIFT:
                return "Frame_Shift";
            case VC_MISSTART:
                return "Misstart";
            case VC_NONSENSE:
                return "Nonsense";
            default:
                CGA_ASSERT(false);
                return "Unknown";
            }
        }

        bool canSplitByCallsOfBothFiles(const CallDiffResult& dr)
        {
            for(size_t ii=0; ii<dr.matchType_.size(); ii++)
            {
                if (dr.matchType_[ii] < dr.matchTypeBySegment_[ii])
                    return false;
//                 if (dr.isMismatch(ii) && dr.inconsistentSegmentCount_[ii] > 1)
//                     return false;
            }
            return true;
        }

        bool allAllelesAreIdentical(const CallDiffResult& dr)
        {
            for(size_t ii=0; ii<dr.matchType_.size(); ii++)
            {
                if (!CallDiffResult::isIdentical(dr.matchType_[ii]))
                    return false;
            }
            return true;
        }

        void addVariantRanges(std::vector<Range>& ranges,
                              const std::pair<std::deque<Locus>::const_iterator,
                              std::deque<Locus>::const_iterator>& loci)
        {
            BOOST_FOREACH(const Locus& locus, loci)
            {
                BOOST_FOREACH(const Call& call, locus.getCalls())
                {
                    if (!call.isRefConsistent(crr_))
                        ranges.push_back(call.range_);
                }
            }
        }

        void orderRangesForMaf(std::vector<Range>& ranges)
        {
            if (ranges.size() <= 1)
                return;

            std::sort(ranges.begin(), ranges.end());

            size_t to = 0;
            bool sticky = 0 == ranges[to].length();
            for(size_t from=1; from<ranges.size(); from++)
            {
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

        void addMafRecords(const std::vector<Range>& ranges,
                           const boost::array<PhasedHypothesis, 2>& hyp)
        {
            BOOST_FOREACH(const Range& range, ranges)
            {
                recs_.push_back(MafRecord());

                MafRecord& rec = recs_.back();
                rec.range_ = range;
                rec.alleles_[0][0]  = getAllele(hyp[0][0], range);
                rec.alleles_[0][1]  = getAllele(hyp[0][1], range);
                rec.alleles_[1][0] = getAllele(hyp[1][0], range);
                rec.alleles_[1][1] = getAllele(hyp[1][1], range);
                rec.dbSnpRs_ = getDbSnpRs(hyp, range);

                addAnnotation(rec, hyp);
            }
        }

        std::string getAllele(const PhasedAllele& allele, const Range& range)
        {
            string result;
            BOOST_FOREACH(const Call* call, allele.calls())
            {
                if (range.beginLocation() <= call->range_.beginLocation() &&
                    call->range_.endLocation() <= range.endLocation())
                {
                    result += call->calledSequence(crr_);
                }
                else if (range.beginLocation() < call->range_.endLocation() &&
                         call->range_.beginLocation() < range.endLocation())
                {
                    if (string::npos != call->alleleSeq_.find('?') ||
                        string::npos != call->alleleSeq_.find('N'))
                        result += "?";
                    else
                    {
                        Range subRange(range);
                        subRange.begin_ = std::max(range.begin_, call->range_.begin_);
                        subRange.end_ = std::min(range.end_, call->range_.end_);
                        result += crr_.getSequence(subRange);
                    }
                }
            }
            return result;
        }

        std::string getDbSnpRs(const boost::array<PhasedHypothesis, 2>& hyps,
                               const Range& range)
        {
            set<int64_t> rsSet;
            BOOST_FOREACH(const PhasedHypothesis& hyp, hyps)
            {
                for(size_t ii=0; ii<hyp.size(); ii++)
                {
                    const PhasedAllele& allele = hyp[ii];

                    BOOST_FOREACH(const Call* call, allele.calls())
                    {
                        if ( ( range.beginLocation() < call->range_.endLocation() &&
                               call->range_.beginLocation() < range.endLocation() ) ||
                             ( 0 == call->range_.length() &&
                               range.beginLocation() <= call->range_.endLocation() &&
                               call->range_.beginLocation() <= range.endLocation() ) )
                        {
                            if (0 != call->xRef_.size())
                                addDbSnpRs(rsSet, call->xRef_);
                        }
                    }
                }
            }
            string result;
            BOOST_FOREACH(int64_t rs, rsSet)
            {
                if (result.size() != 0)
                    result += ";";
                result += "rs";
                result += boost::lexical_cast<string>(rs);
            }
            return result;
        }

        void addDbSnpRs(std::set<int64_t>& rsSet, const std::string& xRef)
        {
            vector<string> parts;
            ba::split(parts, xRef, ba::is_any_of(";"));
            BOOST_FOREACH(const string& part, parts)
            {
                if (ba::starts_with(part, "dbsnp"))
                {
                    size_t pos = part.find(":rs");
                    if (string::npos != pos)
                    {
                        try
                        {
                            rsSet.insert(util::parseValue<int64_t>(part.substr(pos+3)));
                        }
                        catch (const std::exception&)
                        {
                            // ignore parse failure
                        }
                    }
                }
            }
        }

        void addAnnotation(MafRecord& rec, const boost::array<PhasedHypothesis, 2>& hyp)
        {
            rec.refAllele_ = crr_.getSequence(rec.range_);
            for(size_t ii=0; ii<2; ii++)
            {
                if (bu::isCalledSequence(rec.alleles_[ii][0]) && bu::isCalledSequence(rec.alleles_[ii][1]) &&
                    (rec.alleles_[ii][0] != rec.refAllele_ || rec.alleles_[ii][1] != rec.refAllele_))
                {
                    addAlleleCountAnnotation(rec, hyp[ii], ii);
                }
            }
        }

        class ByCallRangeLt
        {
        public:
            bool operator()(const std::pair<const Call*,size_t>& lhs,
                            const std::pair<const Call*,size_t>& rhs) const
            {
                return lhs.first->range_ < rhs.first->range_;
            }
        };

        void addAlleleCountAnnotation(MafRecord& rec,
                                      const PhasedHypothesis& hyp,
                                      size_t fileIdx)
        {
            const Range& range = rec.range_;
            std::vector< std::pair<const Call*,size_t> > orderedCalls;
            for(size_t ii=0; ii<hyp.size(); ii++)
            {
                const PhasedAllele& allele = hyp[ii];

                BOOST_FOREACH(const Call* call, allele.calls())
                {
                    if (!call->isRefConsistent(crr_) && bu::isCalledSequence(call->calledSequence(crr_)) &&
                        range.beginLocation() <= call->range_.beginLocation() &&
                        call->range_.endLocation() <= range.endLocation())
                    {
                        orderedCalls.push_back(make_pair(call, ii));
                    }
                }

            }

            std::sort(orderedCalls.begin(), orderedCalls.end(), ByCallRangeLt());

            for(size_t ii=0; ii<orderedCalls.size(); ii++)
            {
                const Call& call = *orderedCalls[ii].first;
                size_t alleleIdx = orderedCalls[ii].second;
                if (addAlleleCountAnnotation(call, rec, fileIdx, alleleIdx))
                    return;
            }
        }

        bool addAlleleCountAnnotation(const Call& call,
                                      MafRecord& rec,
                                      size_t fileIdx,
                                      size_t alleleIdx)
        {
            evReader_[fileIdx]->seek(call.range_);

            if (!evReader_[fileIdx]->inInterval())
                return false;

            const EvidenceReader::IntervalRecord interval = evReader_[fileIdx]->getInterval();
            bool hetAltAltInterval = std::find(interval.alleleIndexes_.begin(),
                                               interval.alleleIndexes_.end(),
                                               2) != interval.alleleIndexes_.end();
            bool hetRec = rec.alleles_[fileIdx][0] != rec.alleles_[fileIdx][1];
            bool altAltRec = rec.alleles_[fileIdx][0] != rec.refAllele_ &&
                rec.alleles_[fileIdx][1] != rec.refAllele_;

            boost::array<bool, 3> aiIsRef; // allele index supports ref for this MAF locus?
            aiIsRef[0] = true;
            aiIsRef[1] = false;
            aiIsRef[2] = false;
            if (hetAltAltInterval && hetRec && !altAltRec)
            {
                // The variation interval is het alt-alt, but this MAF
                // record is for a subset of that variation interval
                // that is het ref-alt. Determine if allele 1 or allele
                // 2 is ref for the region defined by the MAF record.
                bool compat0 = interval.isCompatible(interval.alleleIndexes_[0], call, crr_);
                bool compat1 = interval.isCompatible(interval.alleleIndexes_[1], call, crr_);
                if (compat0 == compat1)
                    return false; // This call cannot distinguish.

                if (compat0)
                    aiIsRef[interval.alleleIndexes_[1]] = true;
                else
                    aiIsRef[interval.alleleIndexes_[0]] = true;
            }

            rec.totalReadCounts_[fileIdx] = 0;
            rec.altReadCounts_[fileIdx] = 0;
            set<string> dnbsSeen;
            BOOST_FOREACH(const EvidenceReader::DnbRecord& dnb, evReader_[fileIdx]->getDnbs())
            {
                string id = dnb.getId();
                if (dnbsSeen.count(id) != 0)
                    continue;
                dnbsSeen.insert(id);

                int32_t maxRef = dnb.scoreAllele_[0];
                int32_t maxAlt = dnb.scoreAllele_[1];
                if (hetAltAltInterval)
                {
                    if (hetRec && !altAltRec)
                    {
                        if (aiIsRef[1])
                        {
                            CGA_ASSERT(!aiIsRef[2]);
                            maxRef = std::max(maxRef, dnb.scoreAllele_[1]);
                            maxAlt = dnb.scoreAllele_[2];
                        }
                        else
                        {
                            CGA_ASSERT(aiIsRef[2]);
                            maxRef = std::max(maxRef, dnb.scoreAllele_[2]);
                        }
                    }
                    else
                        maxAlt = std::max(maxAlt, dnb.scoreAllele_[2]);
                }

                const int32_t SCORE_THRESHOLD = 3;
                if (maxRef + SCORE_THRESHOLD <= maxAlt)
                {
                    rec.altReadCounts_[fileIdx]++;
                    rec.totalReadCounts_[fileIdx]++;
                }
                else if (maxAlt + SCORE_THRESHOLD <= maxRef)
                {
                    rec.totalReadCounts_[fileIdx]++;
                }
            }

            return true;
        }

        void writeHeader(std::ostream& out) const
        {
            out << "Hugo_Symbol\tEntrez_Gene_Id\tGSC Center\tNCBI_Build\tChromosome\t"
                "Start_position\tEnd_position\tStrand\tVariant_Classification\tVariant_Type\t"
                "Reference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tdbSNP_RS\t"
                "dbSNP_Val_Status\tTumor_Sample_Barcode\tMatched_Norm_Sample_Barcode\t"
                "Match_Norm_Seq_Allele1\tMatch_Norm_Seq_Allele2\tTumor_Validation_Allele1\t"
                "Tumor_Validation_Allele2\tMatch_Norm_Validation_Allele1\t"
                "Match_Norm_Validation_Allele2\tVerification_Status\tValidation_Status\t"
                "Mutation_Status\tValidation_Method\t"
                "PFAM_DOMAIN\tmiRNA\tTumor_VarScore_Rank\t"
                "Match_Norm_RefScore_Rank\tTumor_ReadCount_Alt\t"
                "Tumor_ReadCount_Total\tNormal_ReadCount_Alt\t"
                "Normal_ReadCount_Total" << endl;
        }

        void writeRecord(std::ostream& out, const MafRecord& rec) const
        {
            string chrom = crr_.listChromosomes()[rec.range_.chromosome_].getName();
            if (ba::starts_with(chrom, "chr"))
                chrom = chrom.substr(3);
            writeSemiSeparated(out, rec.hugoSymbol_);
            out << "\t";
            writeSemiSeparated(out, rec.entrezGeneId_);
            out << "\t";
            uint32_t begin = rec.range_.begin_;
            uint32_t end = rec.range_.end_;
            if (begin != end)
                begin++;
            out
                << gscCenter_ << "\t"
                << ncbiBuild_ << "\t"
                << chrom << "\t"
                << begin << "\t"
                << end << "\t"
                << "+" << "\t";
            writeSemiSeparated(out, rec.variantClassification_);
            out << "\t";
            const char* mutationStatus = rec.getMutationStatus();
            const char* validationMethod = "";
            if (0 == strcmp(mutationStatus, "Somatic") &&
                rec.varScoreARank_ >= 0.025 &&
                rec.refScoreBRank_ >= 0.025)
                validationMethod = "bioinformatic";
            out
                << rec.getVariantType() << "\t"
                << rec.refAllele_ << "\t"
                << rec.alleles_[0][0] << "\t"
                << rec.alleles_[0][1] << "\t"
                << rec.dbSnpRs_ << "\t"
                << "" << "\t"
                << tumorBarCode_ << "\t"
                << normalBarCode_ << "\t"
                << rec.alleles_[1][0] << "\t"
                << rec.alleles_[1][1] << "\t"
                << "?" << "\t" // Tumor_Validation_Allele1
                << "?" << "\t" // Tumor_Validation_Allele2
                << "?" << "\t" // Match_Norm_Validation_Allele1
                << "?" << "\t" // Match_Norm_Validation_Allele2
                << "Unknown" << "\t" // Verification_Status
                << "Unknown" << "\t" // Validation_Status
                << mutationStatus << "\t"
                << validationMethod << "\t" // Validation_Method
                ;
            writeSemiSeparated(out, rec.pfamDomain_);
            out << "\t";
            writeSemiSeparated(out, rec.miRna_);
            out << "\t";
            if (0 != strcmp(mutationStatus, "Somatic"))
                out << "\t";
            else if (rec.refScoreBRank_ < 0.0)
                out << "0.000\t0.000";
            else
                out << boost::format("%.3f") % rec.varScoreARank_
                    << "\t" << boost::format("%.3f") % rec.refScoreBRank_;
            for(size_t ii=0; ii<2; ii++)
            {
                if (rec.altReadCounts_[ii] >= 0)
                    out << "\t" << rec.altReadCounts_[ii];
                else
                    out << "\t";
                if (rec.totalReadCounts_[ii] >= 0)
                    out << "\t" << rec.totalReadCounts_[ii];
                else
                    out << "\t";
            }
            out << "\n";
        }

        template <typename TT>
        void writeSemiSeparated(std::ostream& out, const std::vector<TT>& vals) const
        {
            for(size_t ii=0; ii<vals.size(); ii++)
            {
                if (0 != ii)
                    out << ";";
                out << vals[ii];
            }
        }

        void writeSemiSeparated(std::ostream& out, const std::vector<VariantClassification>& vals) const
        {
            for(size_t ii=0; ii<vals.size(); ii++)
            {
                if (0 != ii)
                    out << ";";
                out << vcToString(vals[ii]);
            }
        }

        bool hasRefInconsistentCalls(const Superlocus& sl, size_t fileId)
        {
            BOOST_FOREACH(const Locus& locus, sl.getLoci(fileId))
            {
                if (!locus.isRefConsistent())
                    return true;
            }
            return false;
        }

        bool isFullyCalled(const Superlocus& sl, size_t fileId)
        {
            BOOST_FOREACH(const Locus& locus, sl.getLoci(fileId))
            {
                if (locus.hasNoCalls())
                    return false;
            }
            return true;
        }

        std::string fn_;
        std::string exportRootA_;
        std::string exportRootB_;
        GenomeMetadata metaA_, metaB_;
        SomaticScoreCalc& somatic_;
        const CrrFile& crr_;
        std::vector<MafRecord> recs_;
        std::string gscCenter_;
        std::string ncbiBuild_;
        std::string tumorBarCode_;
        std::string normalBarCode_;
        boost::array< boost::shared_ptr<EvidenceReader>, 2> evReader_;
    };

    void CallDiff::writeVariantOutput(const CrrFile& crr,
                                      const string& inpFn, const string& outFn,
                                      const LocusDiffClassStore& ldc,
                                      SomaticScoreCalc* somatic) const
    {
        const char SEP = '\t';

        VariantFileIterator locIt(crr);
        locIt.open(inpFn);
        if (noReferenceCoverValidation_)
            locIt.setReferenceCoverValidation(false);
        shared_ptr<std::ostream> out =
                OutputStream::openCompressedOutputStreamByExtension(outFn);

        // Prepare somatic data
        map<uint32_t, const VarData*> somaticData;
        if (somatic)
        {
            const vector<VarData>& final = somatic->getFinalData();
            BOOST_FOREACH(const VarData& vd, final)
                somaticData[vd.call_.locus_] = &vd;
        }

        // Write header
        DelimitedFileMetadata meta;
        locIt.fillOlplFileMetadata(meta);
        meta.sort();
        *out << meta;
        locIt.writeOlplFileHeader(*out);
        *out << SEP << "locusDiffClassification";
        if (somatic)
            *out << SEP << "somaticCategory"
                 << SEP << "somaticRank"
                 << SEP << "somaticScore"
                 << SEP << "somaticQuality"
                ;
        *out << '\n';

        // Copy input file to output file with extgra annotations
        while (!locIt.eof())
        {
            locIt->writeAsOneLine(*out, true);
            *out << SEP;
            LocusDiffClassStore::const_iterator ii = ldc.find(locIt->getId());
            if (ii != ldc.end())
                *out << CallDiffResult::getMatchTypeString(ii->second);
            if (somatic)
            {
                if (0 == somaticData.count(locIt->getId()))
                {
                    *out << SEP //<< "somaticCategory"
                         << SEP //<< "somaticRank"
                         << SEP //<< "somaticScore"
                         << SEP //<< "somaticQuality"
                        ;
                }
                else
                {
                    const VarData& v = *somaticData[locIt->getId()];
                    *out << SEP << v.somaticCategory_
                         << "\t" << boost::format("%.3f") % (1.0-v.somaticSensitivity_)
                         << "\t" << boost::format("%.0f") % v.somaticScore_
                         << "\t" << ( v.somaticScore_ >= SOMATIC_QUAL_THRESH ? "SQHIGH" : "" )
                        ;
                }
            }
            *out << '\n';
            ++locIt;
        }
    }

    namespace {
        void checkAsmId(
                const DelimitedFileMetadata& varMeta,
                const std::string& root)
        {
            if (root.empty())
                return;

            const string KEY = "ASSEMBLY_ID";
            if (!varMeta.hasKey(KEY)) // Old files didn't have assembly ID
                return;

            GenomeMetadata mtd(root);

            if (varMeta.get(KEY) != mtd.getAsmId())
            {
                throw Exception("variation file assembly ID '" + varMeta.get(KEY) +
                                "' doesn't match data package assembly ID '" +
                                mtd.getAsmId() + "'");
            }

        }
    }

    int CallDiff::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "variantsA");
        requireParam(vm, "variantsB");

        StringSet reports(reports_,
                          "SuperlocusOutput,SuperlocusStats,LocusOutput,LocusStats,VariantOutput,"
                          "SomaticOutput,DebugCallOutput,DebugSuperlocusOutput,DebugSomaticOutput,Maf",
                          "unknown report type");
        StringSet betaReports("VariantOutput,SomaticOutput,DebugSomaticOutput,Maf", "", "");

        shared_ptr<std::ostream> superOut, superDebug, locusOut, callOut;
        if (0 != reports.count("SuperlocusOutput"))
        {
            superOut = OutputStream::openCompressedOutputStreamByExtension(
                oPrefix_ + "SuperlocusOutput.tsv");
            *superOut << "SuperlocusId\tChromosome\tBegin\tEnd\tClassification\tReference\t"
                "AllelesA\tAllelesB\n";
        }
        if (0 != reports.count("DebugSuperlocusOutput"))
        {
            superDebug = OutputStream::openCompressedOutputStreamByExtension(
                oPrefix_ + "DebugSuperlocusOutput.tsv");
        }
        if (0 != reports.count("LocusOutput"))
        {
            locusOut = OutputStream::openCompressedOutputStreamByExtension(
                oPrefix_ + "LocusOutput.tsv");
            *locusOut << "SuperlocusId\tFile\tLocusClassification\tLocusDiffClassification\t"
                      << Call::getHeader();
            *locusOut << "\n";
        }
        if (0 != reports.count("DebugCallOutput"))
        {
            callOut = OutputStream::openCompressedOutputStreamByExtension(
                oPrefix_ + "DebugCallOutput.tsv");
            *callOut << "SuperlocusId\tFile\tSuperlocusClassification\tCallClassification\t"
                     << Call::getHeader() << "\n";
        }

        reference::CrrFile crr(referenceFileName_);

        VariantFileIterator locItA(crr), locItB(crr);
        locItA.open(variantFileNameA_);
        locItB.open(variantFileNameB_);

        checkAsmId(locItA.getMetadata(), exportRootA_);
        checkAsmId(locItB.getMetadata(), exportRootB_);

        if (noReferenceCoverValidation_)
        {
            locItA.setReferenceCoverValidation(false);
            locItB.setReferenceCoverValidation(false);
        }

        boost::array<LocusDiffClassStore, 2> locusDiffClassCache;
        bool needVariantOutput = reports.contains("VariantOutput");

        LocusStats locusStats;
        SuperlocusStats superlocusStats;

        SuperlocusIterator slIt(extend3Mers_, extendBases_);
        slIt.setVariantFiles(locItA, locItB);

        if (0 == vm.count("beta"))
        {
            BOOST_FOREACH(const string& rep, reports)
                if (betaReports.contains(rep))
                    throw util::Exception("--beta flag required to enable " + rep);
        }

        boost::scoped_ptr<SomaticScoreCalc> somaticScoreCalc;
        if (reports.contains("SomaticOutput") || reports.contains("Maf"))
        {
            if (vm.count("genome-rootA") != 1 || vm.count("genome-rootB") != 1)
            {
                throw util::Exception("for SomaticOutput, genome-rootA and genome-rootB must be specified");
            }
            if (vm.count("calibration-root") != 1)
            {
                throw util::Exception("for SomaticOutput, calibration-root must be specified");
            }
            string debugFileName;
            if (0 != reports.count("DebugSomaticOutput"))
            {
                debugFileName = oPrefix_ + "DebugSomaticOutput.tsv";
            }
            string swVersion = locItA.getMetadata().getSoftwareVersionString();
            if ("" == swVersion)
                throw Exception("no SOFTWARE_VERSION defined in "+variantFileNameA_);

            if (somaticScoreDebugInput_ != "")
            {
                somaticScoreCalc.reset(
                    new SomaticScoreCalc(
                        exportRootA_, exportRootB_, crr, eaf_, debugFileName,
                        calibPrefix_,
                        swVersion,
                        somaticScoreDebugInput_));
                shared_ptr<std::ostream> pSom =
                    OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "SomaticOutput.tsv");
                somaticScoreCalc->dump(*pSom);
                return 0;
            }
            else
                somaticScoreCalc.reset(
                    new SomaticScoreCalc(
                        exportRootA_, exportRootB_, crr, eaf_, debugFileName,
                        calibPrefix_,
                        swVersion));
        }
#if 0
        // Maf code no longer supported. Switching to VCF.
        boost::scoped_ptr<MafReport> mafReport;
        if (0 != reports.count("Maf"))
        {
            if (GenomeMetadata(exportRootA_).getFormatVersion() < 1004 ||
                GenomeMetadata(exportRootB_).getFormatVersion() < 1004)
            {
                throw Exception("format version too old: for Maf report, required format version is >= 1.4");
            }
            mafReport.reset(new MafReport(oPrefix_+"Maf.tsv", exportRootA_, exportRootB_,
                                          *somaticScoreCalc, crr));
        }
#endif

        size_t slCount = 0;
        for(slIt.seekFirst(); !slIt.eof(); ++slIt, slCount++)
        {
//             if (slCount >= 10000)
//                 break;
//             if (slIt->getRange().chromosome_ > 1)
//                 break;

            const Superlocus& sl = *slIt;
            Range range = sl.getRange();

            vector< vector<PhasedHypothesis> > hypotheses;
            sl.buildPhasedHypotheses(hypotheses, maxHypothesisCount_, true);
            CGA_ASSERT(hypotheses.size() == 2);

            bool matchedRegion =
                !hypotheses[0].empty() && !hypotheses[1].empty() &&
                hypotheses[0][0].size() == hypotheses[1][0].size();
            bool matchedHaploidRegion = matchedRegion && hypotheses[0][0].size() == 1;

            if (somaticScoreCalc && matchedRegion)
                somaticScoreCalc->moveToSuperlocus(range, slIt, matchedHaploidRegion);

            if (0 != superDebug.get())
                *superDebug << range.chromosome_ << " " << range.begin_ << " " << range.end_ << " "
                            << hypotheses[0].size() << " " << hypotheses[1].size() << endl;
            CallDiffResult dr;
            PhasedHypothesis::findBestDiff(sl, maxHypothesisCount_, hypotheses[0], hypotheses[1], crr, dr);
            string classification = CallDiffResult::getMatchTypeString(dr.matchType_);
            superlocusStats[dr.matchType_]++;
            string refAllele = crr.getSequence(range);
            size_t refBAlleleCount = 0;
            for(size_t ii=0; ii<dr.hyp_[1].size(); ii++)
            {
                if (dr.hyp_[1][ii].allele() == refAllele)
                    refBAlleleCount++;
            }
#if 0
            if (0 != mafReport.get() && matchedRegion)
                mafReport->addMafRecords(dr, sl);
#endif
            if (0 != superOut.get())
            {
                *superOut << sl.getId() << "\t"
//                           << "A" << boost::size(sl.getLoci(0)) << "B" << boost::size(sl.getLoci(1)) << "\t"
                          << crr.listChromosomes()[range.chromosome_].getName() << "\t"
                          << range.begin_ << "\t" << range.end_ << "\t"
                          << classification << "\t" << crr.getSequence(range);
                *superOut << "\t";
                for(size_t ii=0; ii<dr.hyp_[0].size(); ii++)
                {
                    if (ii > 0)
                        *superOut << ";";
                    *superOut << dr.hyp_[0][ii].allele();
                }
                *superOut << "\t";
                for(size_t ii=0; ii<dr.hyp_[1].size(); ii++)
                {
                    if (ii > 0)
                        *superOut << ";";
                    *superOut << dr.hyp_[1][ii].allele();
                }
                *superOut << "\n";
            }
            if (0 != callOut.get())
            {
                for(size_t ii=0; ii<2; ii++)
                {
                    BOOST_FOREACH(const Locus& locus, sl.getLoci(ii))
                    {
                        BOOST_FOREACH(const Call& call, locus.getCalls())
                        {
                            vector<cdmt::MatchType> mtCall;
                            for(size_t jj=0; jj<dr.callClass_[ii].size(); jj++)
                            {
                                for(size_t kk=0; kk<dr.callClass_[ii][jj].size(); kk++)
                                {
                                    if (dr.callClass_[ii][jj][kk].second == &call)
                                        mtCall.push_back(dr.callClass_[ii][jj][kk].first);
                                }
                            }
                            string callClass = CallDiffResult::getMatchTypeString(mtCall);
                            *callOut << sl.getId() << "\t" << char('A'+ii) << "\t" << classification << "\t";
                            *callOut << callClass << "\t";
                            call.write(*callOut, crr);
                            *callOut << "\n";
                        }
                    }
                }
            }

            for(size_t ii=0; ii<2; ii++)
            {
                BOOST_FOREACH(const Locus& locus, sl.getLoci(ii))
                {
                    if (locus.isRefConsistent())
                        continue;

                    vector<cdmt::MatchType> mtLocus;
                    BOOST_FOREACH(const Allele& allele, locus.getAlleles())
                    {
                        cdmt::MatchType mt = cdmt::REF_IDENTICAL;
                        BOOST_FOREACH(size_t offset, allele.getCallOffsets())
                        {
                            const Call& call = locus.getCalls()[offset];
                            for(size_t jj=0; jj<dr.callClass_[ii].size(); jj++)
                            {
                                for(size_t kk=0; kk<dr.callClass_[ii][jj].size(); kk++)
                                {
                                    if (dr.callClass_[ii][jj][kk].second == &call)
                                        mt = CallDiffResult::mergeMatchTypes(
                                            mt, dr.callClass_[ii][jj][kk].first);
                                }
                            }
                        }
                        mtLocus.push_back(mt);
                    }
                    std::sort(mtLocus.begin(), mtLocus.end());
                    if (needVariantOutput)
                        locusDiffClassCache[ii][locus.getId()] = mtLocus;

                    string locusClass = getClassification(locus, crr);
                    locusStats[ii][locusClass][mtLocus]++;
                    if (0 != locusOut.get())
                    {
                        string locusDiffClass = CallDiffResult::getMatchTypeString(mtLocus);
                        BOOST_FOREACH(const Call& call, locus.getCalls())
                        {
                            *locusOut << sl.getId() << "\t" << char('A'+ii) << "\t";
                            *locusOut << locusClass << "\t" << locusDiffClass << "\t";
                            call.write(*locusOut, crr);
                            *locusOut << "\n";
                        }
                    }
                    if (somaticScoreCalc && matchedRegion && 0 == ii)
                        somaticScoreCalc->processLocus(
                            sl.getId(), locus, locItA, locusClass, mtLocus, refBAlleleCount);
                }
            }
        }

        if (0 != reports.count("LocusStats"))
        {
            shared_ptr<std::ostream> pStats =
                OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "LocusStats.tsv");
            dumpLocusStats(*pStats, locusStats, true);
            dumpLocusStats(*pStats, locusStats, false);
        }

        if (0 != reports.count("SuperlocusStats"))
        {
            shared_ptr<std::ostream> pStats =
                OutputStream::openCompressedOutputStreamByExtension(
                    oPrefix_ + "SuperlocusStats.tsv");
            dumpSuperlocusStats(*pStats, superlocusStats);
        }

        if (0 != reports.count("SomaticOutput"))
        {
            shared_ptr<std::ostream> pSom =
                OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "SomaticOutput.tsv");
            somaticScoreCalc->dump(*pSom);
        }

        if (needVariantOutput)
        {
            writeVariantOutput(crr, variantFileNameA_, oPrefix_ + "VariantsA.tsv",
                               locusDiffClassCache[0], somaticScoreCalc.get());
            writeVariantOutput(crr, variantFileNameB_, oPrefix_ + "VariantsB.tsv",
                               locusDiffClassCache[1], 0);
        }

#if 0
        if (0 != reports.count("Maf"))
        {
            mafReport->finish();
        }
#endif

        return 0;
    }

    std::string CallDiff::getClassification(
        const variants::Locus& locus, const reference::CrrFile& crr) const
    {
        string hetModifier = "";
        string varType = "no-call";
        const Call* prevCall = 0;
        const vector<Allele> alleles = locus.getAlleles();
        const string allele0 = alleles[0].getAlleleSequence();

        if (alleles.size() > 1)
            hetModifier = "hom-";
        else
            hetModifier = "hap-";

        for(size_t ii=0; ii<alleles.size(); ii++)
        {
            if (alleles[ii].hasNoCalls())
            {
                hetModifier = "no-call-";
                continue;
            }

            if (allele0 != alleles[ii].getAlleleSequence() && hetModifier != "no-call-")
                hetModifier = "het-";

            BOOST_FOREACH(const Call& call, locus.getCalls())
            {
                if (!call.isRefConsistent(crr))
                {
                    if (0 != prevCall)
                    {
                        if (prevCall->range_ != call.range_ || prevCall->alleleSeq_ != call.alleleSeq_)
                            varType = "other";
                        continue;
                    }
                    varType = call.varType_;
                    prevCall = &call;
                }
            }
        }
        if (1 == alleles.size())
        {
            if (alleles[0].hasNoCalls())
                return hetModifier + "no-call";
        }

        return hetModifier + varType;
    }

    void CallDiff::dumpLocusStats(std::ostream& out,
                                  LocusStats& locusStats,
                                  bool pct) const
    {
        // Determine totals by match classification and locus
        // classification.
        map< vector<cdmt::MatchType>, int > mtMap;
        array< map< vector<cdmt::MatchType>, int >, 2> mtMapByFile;
        array< map< string, int >, 2 > ltMapByFile;
        array< int, 2 > sumTotalByFile;
        sumTotalByFile.assign(0);
        for(size_t ii=0; ii<locusStats.size(); ii++)
        {
            map<string, map<vector<cdmt::MatchType>, int> >::const_iterator first,
                last=locusStats[ii].end();
            for(first=locusStats[ii].begin(); first!=last; ++first)
            {
                map<vector<cdmt::MatchType>, int>::const_iterator first2,
                    last2=first->second.end();
                for(first2=first->second.begin(); first2!=last2; ++first2)
                {
                    mtMap[first2->first] += first2->second;
                    mtMapByFile[ii][first2->first] += first2->second;
                    ltMapByFile[ii][first->first] += first2->second;
                    sumTotalByFile[ii] += first2->second;
                }
            }
        }

        // Determine the ordering for column headers.
        vector< std::pair<int, vector<cdmt::MatchType> > > mtMapByCount;
        {
            map< vector<cdmt::MatchType>, int >::const_iterator first, last=mtMap.end();
            for(first=mtMap.begin(); first!=last; ++first)
                mtMapByCount.push_back(make_pair(first->second, first->first));
        }
        typedef std::greater< std::pair<int, vector<cdmt::MatchType> > > ReverseOrder;
        std::sort(mtMapByCount.begin(), mtMapByCount.end(), ReverseOrder());

        // Write the stats.
        for(size_t ii=0; ii<locusStats.size(); ii++)
        {
            out << "Locus stats for file " << char('A'+ii);
            if (pct)
                out << " by pct of LocusClassification\n";
            else
                out << " by count\n";
            out << "File\tLocusClassification";
            for(size_t jj=0; jj<mtMapByCount.size() && jj<statsColumnCount_; jj++)
                out << "\t" << CallDiffResult::getMatchTypeString(mtMapByCount[jj].second);
            out << "\tother\ttotal\n";

            map<string, map<vector<cdmt::MatchType>, int> >::iterator first,
                last=locusStats[ii].end();
            for(first=locusStats[ii].begin(); first!=last; ++first)
            {
                out << char('A'+ii) << "\t" << first->first;

                int otherSum = 0;
                int totalSum = ltMapByFile[ii][first->first];
                for(size_t jj=0; jj<mtMapByCount.size(); jj++)
                {
                    int count = first->second[mtMapByCount[jj].second];
                    if (jj < statsColumnCount_)
                    {
                        if (pct)
                            out << boost::format("\t%.2f%%") % (100.0 * count / totalSum);
                        else
                            out << boost::format("\t%d") % count;
                    }
                    else
                        otherSum += count;
                }
                if (pct)
                {
                    out << boost::format("\t%.2f%%") % (100.0 * otherSum / totalSum);
                    out << boost::format("\t%.2f%%\n") % 100.0;
                }
                else
                {
                    out << boost::format("\t%d") % otherSum;
                    out << boost::format("\t%d\n") % totalSum;
                }
            }

            // Print totals line.
            out << char('A'+ii) << "\ttotal";
            int otherSum = 0;
            for(size_t jj=0; jj<mtMapByCount.size(); jj++)
            {
                int count = mtMapByFile[ii][mtMapByCount[jj].second];
                if (jj < statsColumnCount_)
                {
                    if (pct)
                        out << boost::format("\t%.2f%%") % (100.0 * count / sumTotalByFile[ii]);
                    else
                        out << boost::format("\t%d") % count;
                }
                else
                    otherSum += mtMapByFile[ii][mtMapByCount[jj].second];
            }
            if (pct)
            {
                out << boost::format("\t%.2f%%") % (100.0 * otherSum / sumTotalByFile[ii]);
                out << boost::format("\t%.2f%%\n") % 100.0;
            }
            else
            {
                out << boost::format("\t%d") % otherSum;
                out << boost::format("\t%d\n") % sumTotalByFile[ii];
            }
            out << "\n";
        }
    }

    void CallDiff::dumpSuperlocusStats(std::ostream& out,
                                       SuperlocusStats& superlocusStats) const
    {
        vector<string> simpleCategories;
        simpleCategories.push_back("identical");
        simpleCategories.push_back("consistent");
        simpleCategories.push_back("onlyA");
        simpleCategories.push_back("onlyB");
        simpleCategories.push_back("mismatch");
        simpleCategories.push_back("phase-mismatch");
        simpleCategories.push_back("ploidy-mismatch");
        map<string, int> simpleSuperlocusStats;
        int sumTotal = 0;
        map< vector<cdmt::MatchType>,int >::const_iterator first, last=superlocusStats.end();
        for(first=superlocusStats.begin(); first!=last; ++first)
        {
            string cat = CallDiffResult::getSimpleMatchTypeString(first->first);
            CGA_ASSERT(std::find(simpleCategories.begin(), simpleCategories.end(), cat) !=
                       simpleCategories.end());
            simpleSuperlocusStats[cat] += first->second;
            sumTotal += first->second;
        }

        out << "SimpleClassification\tCount\tPct\n";
        BOOST_FOREACH(const std::string& cat, simpleCategories)
        {
            out << cat << "\t" << simpleSuperlocusStats[cat]
                << boost::format("\t%.02f%%") % (100.0 * simpleSuperlocusStats[cat] / sumTotal)
                << "\n";
        }
        out << "\n";

        vector< std::pair< int, vector<cdmt::MatchType> > > statsByCount;
        for(first=superlocusStats.begin(); first!=last; ++first)
            statsByCount.push_back(make_pair(first->second, first->first));
        std::sort(statsByCount.begin(), statsByCount.end(),
                  std::greater< pair<int, vector<cdmt::MatchType> > >());

        out << "Classification\tCount\tPct\n";
        for(size_t ii=0; ii<statsByCount.size(); ii++)
            out << CallDiffResult::getMatchTypeString(statsByCount[ii].second)
                << "\t" << statsByCount[ii].first
                << boost::format("\t%.02f%%") % (100.0 * statsByCount[ii].first / sumTotal)
                << "\n";
    }

} } // cgatools::command
