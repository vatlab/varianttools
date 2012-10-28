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

#ifndef CGA_TOOLS_GAPS_ESTIMATOR_HPP
#define CGA_TOOLS_GAPS_ESTIMATOR_HPP

//! @file GapsEstimator.hpp

#include "cgatools/core.hpp"

#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/cgdata/LibraryReader.hpp"
#include "cgatools/cgdata/Dnb.hpp"
#include "cgatools/util/BaseUtil.hpp"

namespace cgatools { namespace mapping { namespace GapEst {

    inline bool isValidBase(char base) {
        return 
            base=='A' || base=='C' || base=='G' || base=='T' ||
            base=='a' || base=='c' || base=='g' || base=='t';
    }

    inline uint8_t baseTypeIndex(char base) {
        switch (base) {
            case 'A':
            case 'a':
                return 0;
            case 'C':
            case 'c':
                return 1;
            case 'G':
            case 'g':
                return 2;
            case 'T':
            case 't':
                return 3;
            default:
                return 4;
        }
    }
    //! Retrieves sequence at a given offset from the end of the clone.
    class SequenceRetriever
    {
    public:
        static const size_t WOBBLE_SEQUENCE_BASE_COUNT = 10;
        typedef boost::array<char,WOBBLE_SEQUENCE_BASE_COUNT> Sequence;

        virtual ~SequenceRetriever() { }

        //! Store in seq the sequence starting at cloneEndOffset from
        //! the end of the clone, working towards the mate gap
        //! sequenceLength bases.
        virtual void getSequence(Sequence& seq, int cloneEndOffset, int sequenceLength) const = 0;
    };

    //! Retrieves all no-calls.
    class NoCallSequenceRetriever : public SequenceRetriever
    {
    public:
        void getSequence(Sequence& seq, int cloneEndOffset, int sequenceLength) const;
    };

    //! Retrieves sequence at a given offset from the end of the clone.
    template <class TSequence>
    class TemplateSequenceRetriever : public SequenceRetriever
    {
    public:
        TemplateSequenceRetriever(const TSequence& sequence,
                                  bool isCircular,
                                  boost::uint32_t cloneEnd,
                                  int side,
                                  int strand)
            : sequence_(sequence),
              isCircular_(isCircular),
              cloneEnd_(cloneEnd),
              strand_(side == strand ? 0 : 1)
        {
        }

        void getSequence(Sequence& seq, int cloneEndOffset, int sequenceLength) const;

    private:
        const TSequence& sequence_;
        bool isCircular_;
        int cloneEnd_;
        int strand_;
    };

    template <class TSequence>
    void TemplateSequenceRetriever<TSequence>::
    getSequence(Sequence& seq, int cloneEndOffset, int sequenceLength) const
    {
        if (0 == strand_)
        {
            int offset = cloneEnd_ + cloneEndOffset;
            for(int ii=0; ii<sequenceLength; ii++, offset++)
            {
                if (offset >= (int)sequence_.endSequencePosition())
                {
                    if (!isCircular_)
                        seq[ii] = 'N';
                    else
                        seq[ii] = sequence_[offset];
                }
                else if (offset < sequence_.startSequencePosition())
                {
                    if (!isCircular_)
                        seq[ii] = 'N';
                    else
                        seq[ii] = sequence_[offset+sequence_.endSequencePosition()];
                }
                else
                    seq[ii] = sequence_[offset];
            }
        }
        else
        {
            int offset = cloneEnd_ - cloneEndOffset - 1;
            for(int ii=0; ii<sequenceLength; ii++, offset--)
            {
                if (offset >= (int)sequence_.endSequencePosition())
                {
                    if (!isCircular_)
                        seq[ii] = 'N';
                    else
                        seq[ii] = util::baseutil::complement(sequence_[offset]);
                }
                else if (offset < sequence_.startSequencePosition())
                {
                    if (!isCircular_)
                        seq[ii] = 'N';
                    else
                        seq[ii] = util::baseutil::
                                complement(sequence_[offset+sequence_.endSequencePosition()]);
                }
                else
                    seq[ii] = util::baseutil::complement(sequence_[offset]);
            }
        }
    }

    class DependentGaps;
    struct GapsEstimatorData;

    typedef cgdata::DnbStructure ArmData;

    //! Estimates wobble (small) gap probabilities, based on the
    //! sequence at a particular offset from the end of the clone.
    class GapsEstimator
    {
    public:
        typedef cgdata::SmallGapTuple GapTuple;
        typedef SequenceRetriever::Sequence Sequence;

        GapsEstimator(const ArmData& armData, size_t side, boost::uint32_t seed);
        //GapsEstimator(const cgdata::HalfDnbStructure::Reads& readConfig,
        //              boost::uint32_t seed);
        ~GapsEstimator();

        //! Load gaps according to the spec.
        void loadGaps(const std::string& gapFile0,const std::string& gapFile1);

        //! Get P(gaps|sequence), where sequence is retrieved using sr.
        double getProbability(GapTuple gaps, const SequenceRetriever& sr) const;

        //! Get P(gaps) where the sequence is unknown.
        double getProbability(GapTuple gaps) const;

        //! Returns as few gap pairs as possible to cover at least the
        //! specified fraction of hDNB population.
        void getGapTuples(double fraction, std::vector<GapTuple>& result) const;

        //! Get a random GapTuple, according to P(gaps|sequence), where
        //! sequence is retrieved using sr.
        GapTuple randomGaps(const SequenceRetriever& sr) const;

        //! The number of wobble (small) gaps described by this GapsEstimator.
        size_t gapCount() const;

    private:
        typedef GapEst::DependentGaps DependentGaps;
        typedef GapEst::GapsEstimatorData GapsEstimatorData;

        std::vector<DependentGaps> dg_;
        boost::shared_ptr<GapsEstimatorData> data_;
        size_t side_;
    };

}}}

#endif      // end of ifndef CGA_TOOLS_GAPS_ESTIMATOR_HPP
