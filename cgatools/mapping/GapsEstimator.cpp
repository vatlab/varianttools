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

#include "GapsEstimator.hpp"

#include "cgatools/util/Streams.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

namespace cgatools { namespace mapping { namespace GapEst {

    using std::map;
    using std::string;
    using std::vector;

    using boost::lexical_cast;
    using boost::uint16_t;
    using boost::uint32_t;
    using boost::uint64_t;
    using boost::int16_t;
    using boost::int32_t;
    using boost::int64_t;

    namespace ba = boost::algorithm;
    namespace bf = boost::filesystem;

    //////////////////////////////////////////////////////////////////
    // NoCallSequenceRetriever
    //////////////////////////////////////////////////////////////////

    void NoCallSequenceRetriever::getSequence(Sequence& seq, int cloneEndOffset, int sequenceLength) const
    {
        for(int ii=0; ii<sequenceLength; ii++)
            seq[ii] = 'N';
    }

    //////////////////////////////////////////////////////////////////
    // GapsEstimatorData
    //////////////////////////////////////////////////////////////////

    struct GapsEstimatorData
    {
        GapsEstimatorData(const cgdata::HalfDnbStructure::Reads& readBaseCount,
                          boost::uint32_t seed)
            : readBaseCount_(readBaseCount),
              gapCount_(readBaseCount.size()-1),
              rand_(RandomNumberGenerator(seed), UniformRealDistribution())
        {
        }

        // Read base count, in order from end of clone toward mate gap.
        cgdata::HalfDnbStructure::Reads readBaseCount_;
        boost::uint32_t gapCount_;
        typedef boost::mt19937 RandomNumberGenerator;
        typedef boost::uniform_real<double> UniformRealDistribution;
        typedef boost::
            variate_generator<RandomNumberGenerator, UniformRealDistribution> RandomRealUniformGenerator;
        RandomRealUniformGenerator rand_;
    };

    //////////////////////////////////////////////////////////////////
    // DependentGaps
    //////////////////////////////////////////////////////////////////

    class DependentGaps
    {
    public:
        typedef GapsEstimator::GapTuple GapTuple;
        typedef GapsEstimator::Sequence Sequence;

        DependentGaps(size_t firstGap);

        void loadConfig(const std::string& filename);

        //! Returns P(gaps|sequence).
        inline double getPGaps(const Sequence& sequence, const GapTuple& gaps) const;

        //! Returns P(gaps|sequence), for all gap combos.
        inline const std::vector<double>& getPGaps(const Sequence& sequence) const;

        //! Returns P(gaps) where the sequence is unknown.
        double getPGaps(const GapTuple& gaps) const;

        //! Returns P(gaps) where the sequence is unknown, for all
        //! gap combos.
        const std::vector<double>& getPGaps() const;

        //! getGapTuple(ii) is the partial gap tuple associated with the
        //! probability getPGaps(seq)[ii]
        inline const GapTuple& getGapTuple(size_t index) const;

        //! Return the first gap described by this instance.
        inline size_t firstGap() const;

        //! Return the number of gaps described by this instance.
        inline size_t gapCount() const;

        //! Return offset of sequence of interest, from current
        //! anchor site.
        inline int getSequenceOffset() const;

        //! Return length of sequence of interest.
        inline int getSequenceLength() const;

    private:
        void init(size_t firstGap);
        Sequence getSequence(const std::string& seq) const;
        size_t getSeqCount(size_t seqLen) const;
        size_t getSequenceIndex(const Sequence& sequence) const;
        size_t getGapTupleIndex(const GapTuple& gaps) const;

        size_t firstGap_;
        size_t count_;
        int offset_;
        int seqLen_;
        // PGaps_[sequenceIndex][gapsCombo]
        std::vector< std::vector<double> > PGaps_;
        std::vector< GapTuple > indices_;
    };

    DependentGaps::DependentGaps(size_t firstGap)
    {
        init(firstGap);
    }

    void DependentGaps::init(size_t firstGap)
    {
        firstGap_ = firstGap;
        count_ = 0;
        offset_ = 0;
        seqLen_ = 0;
        PGaps_.clear();
        indices_.clear();
    }

    void DependentGaps::loadConfig(const std::string& filename)
    {
        init(firstGap_);

        util::InputStream ifs(filename.c_str());
        if (!ifs)
        {
            throw util::Exception("can't load gaps file: " + filename);
        }

        // Read header.
        string header;
        while(getline(ifs, header)) 
            if ((!header.empty()) && header[0]=='>')
            {
                header.erase(0,1);
                break;
            }

        vector<string> headerFields;
        ba::split(headerFields, header, ba::is_any_of("\t"));

        if (headerFields.size() < 2)
            throw util::Exception("bad gaps header: "+header);

        vector<string> headerFirstField;
        ba::split(headerFirstField, headerFields[0], ba::is_any_of(";"));
        if (!ba::starts_with(headerFirstField[0], "sequence:"))
            throw util::Exception("bad gaps header: "+header);

        for(size_t ii=0; ii<headerFirstField.size(); ii++)
        {
            vector<string> parts;
            ba::split(parts, headerFirstField[ii], ba::is_any_of(":"));
            string key = ba::to_lower_copy(parts[0]);
            if ("firstgap" == key)
            {
                if (parts.size() != 2 || firstGap_ != boost::lexical_cast<size_t>(parts[1]))
                    throw util::Exception("firstGap mismatch: "+headerFirstField[ii]);
            }
            else if ("gapcount" == key)
            {
                // Ignore. Gap count is implicit in the other fields.
            }
        }

        string sequenceSpec = headerFirstField[0].substr(9);
        vector<string> sequenceFields;
        ba::split(sequenceFields, sequenceSpec, ba::is_any_of("-"));
        if (sequenceFields.size() != 2)
            throw util::Exception("bad gaps header: "+header);
        offset_ = lexical_cast<int>(sequenceFields[0]);
        seqLen_ = lexical_cast<int>(sequenceFields[1]) - offset_;
        for(size_t ii=1; ii<headerFields.size(); ii++)
        {
            if (!ba::starts_with(headerFields[ii], "gaps:"))
                throw util::Exception("bad gaps header: "+header);

            string gapSpec = headerFields[ii].substr(5);
            vector<string> parts;
            ba::split(parts, gapSpec, ba::is_any_of(";"));

            if (1 == ii)
                count_ = parts.size();
            else if (parts.size() != count_)
                throw util::Exception("bad gaps header: "+header);

            GapTuple gaps;
            gaps.assign(0);
            for(size_t ii=0; ii<count_; ii++)
                gaps[firstGap_+ii] = lexical_cast<int16_t>(parts[ii]);
            indices_.push_back(gaps);
        }

        // Read data.
        map< int, vector<double> > seqMapping;
        size_t rowCount = 0;
        std::string line;
        while (getline(ifs, line))
        {
            vector<string> fields;
            ba::split(fields, line, ba::is_any_of("\t"));

            if (fields.size() != indices_.size() + 1)
                throw util::Exception("wrong field count in gaps file: "+filename);
            if ((int)fields[0].size() != seqLen_)
                throw util::Exception("sequence length mismatch in gaps file: "+filename);

            vector<double> probabilities(indices_.size());
            double sum = 0.0;
            for(size_t ii=0; ii<probabilities.size(); ii++)
            {
                probabilities[ii] = lexical_cast<double>(fields[1+ii]);
                sum += probabilities[ii];
            }
            for(size_t ii=0; ii<probabilities.size(); ii++)
                probabilities[ii] /= sum;

            Sequence seq = getSequence(fields[0]);
            size_t seqIndex = getSequenceIndex(seq);
            seqMapping[seqIndex] = probabilities;
            ++rowCount;
        }

        // Apply seqMapping to PGaps_;
        size_t seqCount = getSeqCount(seqLen_);

        if (rowCount==1 && seqMapping[seqCount-1].size() == indices_.size())
        {
            seqMapping[0] = seqMapping[seqCount-1];
            seqMapping.erase(seqCount-1);
            seqLen_ = 0;
            seqCount = getSeqCount(seqLen_);
        }

        PGaps_.resize(seqCount);
        for(size_t ii=0; ii<seqCount; ii++)
        {
            vector<double>& probabilities = seqMapping[ii];
            if (probabilities.size() != indices_.size())
                throw util::Exception("gaps file missing sequence: "+filename);

            PGaps_[ii] = probabilities;
        }
    }

    DependentGaps::Sequence DependentGaps::getSequence(const std::string& seq) const
    {
        CGA_ASSERT_EQ((int)seq.size(),seqLen_);
        Sequence result;
        result.assign('N');
        for(int ii=0; ii<seqLen_; ii++)
        {
            result[ii] = seq[ii];
        }
        return result;
    }

    size_t DependentGaps::getSeqCount(size_t seqLen) const
    {
        size_t result = 1;
        for (size_t ii=0; ii<seqLen; ii++)
            result *= 5;
        return result;
    }

    double DependentGaps::getPGaps(const Sequence& sequence, const GapTuple& gaps) const
    {
        return PGaps_[getSequenceIndex(sequence)][getGapTupleIndex(gaps)];
    }

    const std::vector<double>& DependentGaps::getPGaps(const Sequence& sequence) const
    {
        size_t index = getSequenceIndex(sequence);
        CGA_ASSERT_L(index,PGaps_.size());
        return PGaps_[index];
    }

    double DependentGaps::getPGaps(const GapTuple& gaps) const
    {
        return getPGaps()[getGapTupleIndex(gaps)];
    }

    const std::vector<double>& DependentGaps::getPGaps() const
    {
        return PGaps_.back();
    }

    const DependentGaps::GapTuple& DependentGaps::getGapTuple(size_t index) const
    {
        return indices_[index];
    }

    size_t DependentGaps::firstGap() const
    {
        return firstGap_;
    }

    size_t DependentGaps::gapCount() const
    {
        return count_;
    }

    int DependentGaps::getSequenceOffset() const
    {
        return offset_;
    }

    int DependentGaps::getSequenceLength() const
    {
        return seqLen_;
    }

    size_t DependentGaps::getSequenceIndex(const Sequence& sequence) const
    {
        size_t result = 0;
        for(int ii=0; ii<seqLen_; ii++)
            result = result*5 + baseTypeIndex(sequence[ii]);

        return result;
    }

    size_t DependentGaps::getGapTupleIndex(const GapTuple& gaps) const
    {
        for(size_t ii=0; ii<indices_.size(); ii++)
        {
            const GapTuple& indexedGaps = indices_[ii];
            size_t jj;
            for(jj=firstGap_; jj<firstGap_+count_; jj++)
            {
                if (gaps[jj] != indexedGaps[jj])
                    break;
            }
            if (firstGap_+count_ == jj)
                return ii;
        }
        // Failed to find index for gaps.
        CGA_ASSERT_MSG(false,"Failed to find index for gaps");
        return 0;
    }

    //////////////////////////////////////////////////////////////////
    // GapsEstimator
    //////////////////////////////////////////////////////////////////

    GapsEstimator::GapsEstimator(const ArmData& armData, size_t side, boost::uint32_t seed)
        : side_ (side)
    {
        cgdata::HalfDnbStructure::Reads readBaseCount = armData.halfDnbs_[side].reads_;

        if (0 != side_)
        {
            for(size_t ii=0; ii<readBaseCount.size()/2; ii++)
                std::swap(readBaseCount[ii], readBaseCount[readBaseCount.size()-ii-1]);
        }

        data_.reset(new GapsEstimatorData(readBaseCount, seed));
    }

    GapsEstimator::~GapsEstimator()
    {
    }

    void GapsEstimator::loadGaps(const std::string& gapFile0,const std::string& gapFile1)
    {
        dg_.clear();

        dg_.push_back(DependentGaps(0));
        dg_.back().loadConfig(gapFile0);
        CGA_ASSERT_EQ(2,dg_.back().gapCount());
        dg_.push_back(DependentGaps(2));
        dg_.back().loadConfig(gapFile1);
        CGA_ASSERT_EQ(1,dg_.back().gapCount());

        if (3 != data_->gapCount_)
            throw util::Exception("wrong gap count in spec'd gaps: 3");
    }

    double GapsEstimator::getProbability(GapTuple gaps, const SequenceRetriever& sr) const
    {
        if (side_>0)
            std::reverse(gaps.begin(),gaps.end());
        double result = 1.0;
        Sequence seq;
        int currentAnchor = 0;
        BOOST_FOREACH(const DependentGaps& dg, dg_)
        {
            int seqLen = dg.getSequenceLength();
            if (seqLen != 0)
            {
                int seqOffset = dg.getSequenceOffset();
                sr.getSequence(seq, currentAnchor+seqOffset, seqLen);
            }
            for(size_t ii=dg.firstGap(); ii<dg.firstGap()+dg.gapCount(); ii++)
                currentAnchor += data_->readBaseCount_[ii] + gaps[ii];
            result *= dg.getPGaps(seq, gaps);
        }
        return result;
    }

    double GapsEstimator::getProbability(GapTuple gaps) const
    {
        if (side_>0)
            std::reverse(gaps.begin(),gaps.end());
        double result = 1.0;
        BOOST_FOREACH(const DependentGaps& dg, dg_)
        {
            result *= dg.getPGaps(gaps);
        }
        return result;
    }

    class DecreasingByDouble
    {
    public:
        DecreasingByDouble(const std::vector<double>& aray)
            : aray_(aray)
        {
        }

        bool operator()(size_t ii, size_t jj)
        {
            return aray_[ii] > aray_[jj];
        }

    private:
        const std::vector<double>& aray_;
    };

    void GapsEstimator::getGapTuples(double fraction, std::vector<GapTuple>& result) const
    {
        result.clear();

        // Build the full set of gap tuples and their likelihoods.
        vector<GapTuple> gapTuples;
        vector<double> probabilities;
        BOOST_FOREACH(const DependentGaps& dg, dg_)
        {
            if (0 == gapTuples.size())
            {
                probabilities = dg.getPGaps();
                for(size_t ii=0; ii<probabilities.size(); ii++)
                    gapTuples.push_back(dg.getGapTuple(ii));
            }
            else
            {
                vector<GapTuple> newGapTuples;
                vector<double> newProbabilities;
                const vector<double>& gts = dg.getPGaps();
                for(size_t ii=0; ii<gapTuples.size(); ii++)
                {
                    for(size_t jj=0; jj<gts.size(); jj++)
                    {
                        newGapTuples.push_back(gapTuples[ii]);
                        GapTuple gt = dg.getGapTuple(jj);
                        for(size_t kk=0; kk<dg.gapCount(); kk++)
                            newGapTuples.back()[dg.firstGap()+kk] = gt[dg.firstGap()+kk];

                        newProbabilities.push_back(probabilities[ii] * gts[jj]);
                    }
                }
                std::swap(gapTuples, newGapTuples);
                std::swap(probabilities, newProbabilities);
            }
        }

        // Order the gapTuples in order of decreasing likelihood.
        vector<size_t> order;
        for(size_t ii=0; ii<gapTuples.size(); ii++)
            order.push_back(ii);
        std::sort(order.begin(), order.end(), DecreasingByDouble(probabilities));

        double fractionSoFar = 0.0;
        for(size_t ii=0; ii<order.size(); ii++)
        {
            result.push_back(gapTuples[order[ii]]);

            fractionSoFar += probabilities[order[ii]];
            if (fractionSoFar >= fraction)
                break;
        }
    }

    GapsEstimator::GapTuple GapsEstimator::randomGaps(const SequenceRetriever& sr) const
    {
        GapTuple result;
        Sequence seq;
        int currentAnchor = 0;
        BOOST_FOREACH(const DependentGaps& dg, dg_)
        {
            int seqLen = dg.getSequenceLength();
            if (seqLen != 0)
            {
                int seqOffset = dg.getSequenceOffset();
                sr.getSequence(seq, currentAnchor+seqOffset, seqLen);
            }
            const vector<double>& PGaps = dg.getPGaps(seq);
            double rVal = data_->rand_();
            size_t ii;
            for(ii=0; ii<PGaps.size(); ii++)
            {
                if (rVal < PGaps[ii])
                    break;
                rVal -= PGaps[ii];
            }
            if (PGaps.size() == ii)
                ii--;
            const GapTuple& resultThisDg = dg.getGapTuple(ii);
            std::copy(resultThisDg.begin()+dg.firstGap(),
                      resultThisDg.begin()+dg.firstGap()+dg.gapCount(),
                      result.begin()+dg.firstGap());
            for(size_t ii=dg.firstGap(); ii<dg.firstGap()+dg.gapCount(); ii++)
                currentAnchor += data_->readBaseCount_[ii] + result[ii];
        }
        return result;
    }

    size_t GapsEstimator::gapCount() const
    {
        size_t result = 0;
        BOOST_FOREACH(const DependentGaps& dg, dg_)
        {
            result += dg.gapCount();
        }
        return result;
    }

}}} // GapEst
