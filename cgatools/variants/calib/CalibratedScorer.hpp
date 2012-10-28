// Copyright 2011 Complete Genomics, Inc.
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

#ifndef CGATOOLS_VARIANTS_CALIB_CALIBRATEDSCORER_HPP_
#define CGATOOLS_VARIANTS_CALIB_CALIBRATEDSCORER_HPP_ 1

//! @file CalibratedScorer.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"

#include <vector>
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace variants { namespace calib {

    class CoverageBinner
    {
    public:
        CoverageBinner()
            : cvgLevels_(1, 0)
        {
        }

        CoverageBinner(const std::vector<int32_t>& cvgLevels);

        size_t getBin(int32_t cvg) const
        {
            size_t result = 0;
            for(size_t ii=0; ii<cvgLevels_.size(); ii++)
            {
                if (cvgLevels_[ii] <= cvg)
                    result = ii;
            }
            return result;
        }

        size_t getBinCount() const
        {
            return cvgLevels_.size();
        }

        int32_t getMinCvg(size_t bin) const
        {
            return cvgLevels_[bin];
        }

    private:
        std::vector< int32_t > cvgLevels_;
    };

    class CalibratedScorer
    {
    public:
        // Constructs CalibratedScorer.
        // @varType The varType from the var file. This is one of snp,
        //          ins, del, or sub.
        // @scoreType The type of score calibration to get. This is one
        //            of:
        //            - fp -- False positive rate of het calls, given
        //              totalReadCount from masterVar file and varScore
        //              of the variant allele from variant file (either
        //              varScoreVAF or varScoreEAF).
        //            - fn -- False negative rate of het loci, given
        //              uniqueSequenceCoverage and refScore.
        //            - uc -- Undercall rate of het calls, given
        //              totalReadCount from masterVar file and varScore
        //              of the reference allele from the variant file.
        //            - oc -- Overcall rate of homozygous calls, given
        //              totalReadCount from masterVar file and the
        //              minimum varScore, out of the two alleles of the
        //              homozygous locus.
        // @eaf True if using EAF scores, false otherwise.
        // @dataPath A path to a directory containing the calibration
        //           data.
        // @softwareVersion The SOFTWARE_VERSION used to assemble this
        //                  genome. E.g. "2.0.0.10".
        // @a20Mixture The fraction of variants expected to be at allele
        //             fraction 20%. All other variants are modeled as
        //             being present at 50% allele fraction.
        // @refBasesPerHetVariant Estimate of number of true reference
        //                        bases for each heterozygous
        //                        variant. Used to shift the estimation
        //                        of fp rate and fn for genomes with
        //                        substantially more or fewer variants
        //                        than the genome used to calibrate
        //                        scores. Generally, fn count varies as
        //                        the count of variants, and fp count
        //                        varies as the count of bases.
        CalibratedScorer(
            const std::string& varType,
            const std::string& scoreType,
            bool eaf,
            const std::string& dataPath,
            const std::string& softwareVersion,
            double a20Mixture = 0.0,
            double refBasesPerHetVariant = 0.0);

        double getCalibratedScore(int32_t cvg, int32_t score) const
        {
            return S_[binner_.getBin(cvg)][scoreToOffset(score)];
        }

        double getCalibratedLikelihoodRatio(int32_t cvg, int32_t score) const
        {
            return L_[binner_.getBin(cvg)][scoreToOffset(score)];
        }

        double getPTrue(int32_t cvg, int32_t score) const
        {
            return PTrue_[binner_.getBin(cvg)][scoreToOffset(score)];
        }

        const CoverageBinner& getBinner() const
        {
            return binner_;
        }

        int32_t getMinScore() const
        {
            return minScore_;
        }

        int32_t getMaxScore() const
        {
            return minScore_ + S_[0].size()-1;
        }

    private:
        void getScoreStreams(
            const std::string& calibId,
            const std::string& dataPath,
            const std::string& softwareVersion,
            boost::shared_ptr<std::istream>& inMetrics,
            boost::shared_ptr<std::istream>& inData,
            std::string& fnMetrics,
            std::string& fnData,
            std::string& fnDataA20) const;

        void readData(
            std::istream& inData,
            const std::string& fnData,
            CoverageBinner& binner,
            std::vector< std::vector<double> >& SS,
            int32_t& minScore);

        size_t scoreToOffset(int32_t score) const
        {
            score = std::max(score, minScore_);
            size_t offset = score-minScore_;
            offset = std::min(offset, S_[0].size()-1);
            return offset;
        }

        int swVersionToInt(const std::string& versionStr) const;

        CoverageBinner binner_;
        std::vector< std::vector<double> > S_, L_, PTrue_;
        int32_t minScore_;
    };

} } } // cgatools::variants::calib

#endif // CGATOOLS_VARIANTS_CALIB_CALIBRATEDSCORER_HPP_
