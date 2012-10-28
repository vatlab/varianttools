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
#include "cgatools/variants/Superlocus.hpp"

#include <boost/foreach.hpp>

namespace cgatools { namespace variants {

    using util::Exception;
    using namespace cgatools::reference;

    using namespace std;

    Superlocus::Superlocus()
        : range_(0,0,0),
          id_(0)
    {
    }

    std::pair<std::deque<Locus>::const_iterator, std::deque<Locus>::const_iterator>
    Superlocus::getLoci(size_t fileOffset) const
    {
        // Skip to first locus.
        deque<Locus>::const_iterator qIter, qLast = queues_[fileOffset].end();
        for(qIter=queues_[fileOffset].begin(); qIter!=qLast; ++qIter)
        {
            const Range& iterRange = qIter->getRange();
            if (iterRange.beginLocation() >= range_.beginLocation() ||
                iterRange.endLocation() > range_.beginLocation())
                break;
        }

        deque<Locus>::const_iterator first = qIter;
        for(; qIter!=qLast; ++qIter)
        {
            if (range_.endLocation() < qIter->getRange().endLocation() &&
                range_.endLocation() <= qIter->getRange().beginLocation())
                break;
        }

        return make_pair(first, qIter);
    }

    uint32_t Superlocus::extendLeftBySuffixMatching(
        const reference::CrrFile& crr, const reference::Location& loc, const std::string& sequence)
    {
        if (0 == sequence.size())
            return loc.offset_;

        const CompactDnaSequence& chromosome = crr.listChromosomes()[loc.chromosome_];
        uint32_t pos = loc.offset_;
        uint32_t alleleOffset = sequence.size();
        while (pos > 0)
        {
            pos--;
            if (0 == alleleOffset)
                alleleOffset = sequence.size();
            alleleOffset--;
            if (sequence[alleleOffset] != chromosome.getBase(pos))
                return pos+1;
        }

        return pos;
    }

    uint32_t Superlocus::extendRightByPrefixMatching(
        const reference::CrrFile& crr, const reference::Location& loc, const std::string& sequence)
    {
        if (0 == sequence.size())
            return loc.offset_;

        const CompactDnaSequence& chromosome = crr.listChromosomes()[loc.chromosome_];
        uint32_t pos = loc.offset_;
        uint32_t alleleOffset = 0;
        while (pos < chromosome.length())
        {
            if (sequence[alleleOffset] != chromosome.getBase(pos))
                break;
            pos++;
            alleleOffset++;
            if (sequence.size() == alleleOffset)
                alleleOffset = 0;
        }

        return pos;
    }

    void Superlocus::buildPhasedHypotheses(std::vector< std::vector<PhasedHypothesis> >& hypotheses,
                                           size_t maxHypothesisCount,
                                           bool useHapLinks) const
    {
        hypotheses.resize(queues_.size());
        for(size_t ii=0; ii<queues_.size(); ii++)
            buildPhasedHypotheses(hypotheses[ii], maxHypothesisCount, useHapLinks, ii);
    }

    void Superlocus::buildPhasedHypotheses(std::vector<PhasedHypothesis>& hypotheses,
                                           size_t maxHypothesisCount,
                                           bool useHapLinks,
                                           size_t fileOffset) const
    {
        hypotheses.clear();

        // Find loci of interest.
        pair<deque<Locus>::const_iterator, deque<Locus>::const_iterator> iters = getLoci(fileOffset);
        deque<Locus>::const_iterator qIter = iters.first, qLast = iters.second;

        if (qIter == qLast)
        {
            // Skipped chromosome.
            return;
        }
        CGA_ASSERT(qIter->getRange().beginLocation() <= range_.beginLocation());
        CGA_ASSERT((qLast-1)->getRange().endLocation() >= range_.endLocation());

        // Determine ploidy.
        uint16_t ploidy = qIter->getPloidy();
        if (Call::UNKNOWN_PLOIDY == ploidy)
        {
            return;
        }
        hypotheses.resize(1, PhasedHypothesis(range_, ploidy));
        vector<size_t> perm(ploidy);

        for(; qIter!=qLast; ++qIter)
        {
            if (0 == hypotheses.size())
                break;

            if (range_.endLocation() < qIter->getRange().endLocation() &&
                range_.endLocation() <= qIter->getRange().beginLocation())
                break;

            const Locus& locus = *qIter;
            if (locus.getPloidy() != ploidy)
            {
                // Ploidy mismatch.
                hypotheses.clear();
                return;
            }

            for(size_t ii=0; ii<perm.size(); ii++)
                perm[ii] = ii;

            if (range_.beginLocation() <= locus.getRange().beginLocation() &&
                range_.endLocation() >= locus.getRange().endLocation())
            {
                size_t hypothesisCount = hypotheses.size();

                vector<size_t> deletedHypotheses;
                for(size_t ii=0; ii<hypothesisCount; ii++)
                {
                    if (hypothesisPermutationsAreEqual(hypotheses[ii]) || allelePermutationsAreEqual(locus))
                        addCalls(hypotheses[ii], locus, perm);
                    else
                    {
                        // Add sequence for all allele permutations that
                        // are consistent with the haplinks.
                        size_t consistentCount = 0;
                        PhasedHypothesis hypothesis = hypotheses[ii]; // intentional copy
                        do
                        {
                            if ( (!useHapLinks) || areHapLinksConsistent(hypothesis, locus, perm) )
                            {
                                if (0 == consistentCount)
                                {
                                    addCalls(hypotheses[ii], locus, perm);
                                }
                                else
                                {
                                    hypotheses.push_back(hypothesis);
                                    addCalls(hypotheses.back(), locus, perm);
                                }
                                consistentCount++;
                            }
                        } while (std::next_permutation(perm.begin(), perm.end()));
                        if (0 == consistentCount)
                        {
                            // No permutations were consistent with
                            // haplinks. Mark this hypothesis as
                            // deleted.
                            deletedHypotheses.push_back(ii);

                            if (0 == ii)
                            {
                                addCalls(hypotheses[ii], locus, perm);
                            }
                        }
                    }
                }

                // Delete hypotheses, starting from end.
                for(int ii=deletedHypotheses.size()-1; ii>=0; ii--)
                {
                    // Don't delete *all* hypotheses. Due to down
                    // sampling, it's possible no remaining hypotheses
                    // are consistent with hapLinks. In that case, just
                    // do a best effort (i.e. possibly return a
                    // hypothesis that is inconsistent with hapLinks).
                    if (hypotheses.size() > 1)
                        hypotheses.erase(hypotheses.begin() + deletedHypotheses[ii]);
                }
            }
            else
            {
                // In the case where the interval doesn't contain the
                // call, this had better either be a big no-call block
                // or ref-call block.
                CGA_ASSERT(locus.isRefCallLocus() || locus.isNoCallLocus());
                Range rg = range_.overlappingRange(locus.getRange());
                string sequence = "?";
                if (locus.isRefCallLocus())
                    sequence = locus.getReference().getSequence(rg);
                BOOST_FOREACH(PhasedHypothesis& hypothesis, hypotheses)
                {
                    addSequence(hypothesis, locus, perm, sequence);
                }
            }

            if (hypotheses.size() > maxHypothesisCount)
                downSample(hypotheses, maxHypothesisCount);
        }

        if (0 == hypotheses.size())
            throw Exception("failed to phase loci: no permutations were consistent with hapLinks");
    }

    // Returns true if permuting the PhasedAlleles of the hypothesis
    // results in an equivalent hypothesis.
    bool Superlocus::hypothesisPermutationsAreEqual(const PhasedHypothesis& hypothesis) const
    {
        for(size_t ii=1; ii<hypothesis.size(); ii++)
        {
            if ( hypothesis[ii].allele() != hypothesis[0].allele() )
                return false;
            if ( hypothesis[ii].calls().size() != hypothesis[0].calls().size() )
                return false;
            for(size_t jj=0; jj<hypothesis[ii].calls().size(); jj++)
            {
                if (!callPermutationsAreEqual(*hypothesis[ii].calls()[jj], *hypothesis[0].calls()[jj]))
                    return false;
            }
        }
        return true;
    }

    // Returns true if permuting the calls results in an equivalent
    // hypothesis.
    bool Superlocus::callPermutationsAreEqual(const Call& lhs, const Call& rhs) const
    {
        if (lhs.range_ != rhs.range_)
            return false;
        if (lhs.alleleSeq_ != rhs.alleleSeq_)
            return false; // conservative check
        if (lhs.varScoreVAF_ != rhs.varScoreVAF_)
            return false;
        if (lhs.varScoreEAF_ != rhs.varScoreEAF_)
            return false;
        if (lhs.varQuality_ != rhs.varQuality_)
            return false;
        if (lhs.hapLink_ != rhs.hapLink_)
            return false;
        return true;
    }

    // Returns true if permuting the alleles results in an equivalent
    // hypothesis.
    bool Superlocus::allelePermutationsAreEqual(const Locus& locus) const
    {
        const vector<Allele>& alleles = locus.getAlleles();
        const vector<Call>& calls = locus.getCalls();
        const vector<size_t>& offsets0 = alleles[0].getCallOffsets();
        for(size_t ii=1; ii<alleles.size(); ii++)
        {
            const vector<size_t>& offsets = alleles[ii].getCallOffsets();
            if (offsets.size() != offsets0.size())
                return false;
            for(size_t jj=0; jj<offsets.size(); jj++)
            {
                if (!callPermutationsAreEqual(calls[offsets[jj]], calls[offsets0[jj]]))
                    return false;
            }
        }
        return true;
    }

    void Superlocus::addCalls(PhasedHypothesis& hypothesis,
                              const Locus& locus,
                              const std::vector<size_t>& perm) const
    {
        const vector<Allele>& alleles = locus.getAlleles();
        const vector<Call>& calls = locus.getCalls();
        for(size_t jj=0; jj<hypothesis.size(); jj++)
        {
            const vector<size_t>& offsets = alleles[perm[jj]].getCallOffsets();
            BOOST_FOREACH(size_t offset, offsets)
            {
                hypothesis[jj].addCall(calls[offset], locus, locus.getReference());
            }
        }
    }

    void Superlocus::addSequence(PhasedHypothesis& hypothesis,
                                 const Locus& locus,
                                 const std::vector<size_t>& perm,
                                 const std::string& sequence) const
    {
        const vector<Allele>& alleles = locus.getAlleles();
        const vector<Call>& calls = locus.getCalls();
        for(size_t jj=0; jj<hypothesis.size(); jj++)
        {
            const vector<size_t>& offsets = alleles[perm[jj]].getCallOffsets();
            CGA_ASSERT(1 == offsets.size());
            hypothesis[jj].addSequence(calls[offsets[0]], locus, sequence);
        }
    }

    bool Superlocus::areHapLinksConsistent(const PhasedHypothesis& hypothesis,
                                           const Locus& locus,
                                           const std::vector<size_t>& perm) const
    {
        const vector<Allele>& alleles = locus.getAlleles();

        for(size_t jj=0; jj<hypothesis.size(); jj++)
        {
            const string& hapLink = alleles[perm[jj]].getHapLink();
            if (hapLink.size() > 0)
            {
                for(size_t ii=0; ii<hypothesis.size(); ii++)
                {
                    if (ii == jj)
                        continue;

                    if (hypothesis[ii].hasHapLink(hapLink))
                        return false;
                }
            }
        }

        return true;
    }

    void Superlocus::downSample(std::vector<PhasedHypothesis>& hypotheses,
                                size_t maxHypothesisCount) const
    {
        vector< pair<uint32_t, size_t> > hypHashes;
        for(size_t ii=0; ii<hypotheses.size(); ii++)
            hypHashes.push_back( make_pair(hashHypothesis(hypotheses[ii]), ii) );
        std::sort(hypHashes.begin(), hypHashes.end());
        vector<PhasedHypothesis> tmp;
        tmp.reserve(2*maxHypothesisCount);
        tmp.resize(maxHypothesisCount, PhasedHypothesis(Range(), 0));
        tmp.swap(hypotheses);
        for(size_t ii=0; ii<maxHypothesisCount; ii++)
            hypotheses[ii].swap(tmp[hypHashes[ii].second]);
    }

    uint32_t Superlocus::hashHypothesis(const PhasedHypothesis& hyp) const
    {
        // based on djb2 hash
        uint32_t hash = 5381;

        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            hash = ((hash << 5) + hash) + 17;
            const string& allele = hyp[ii].allele();
            for(size_t jj=0; jj<allele.size(); jj++)
                hash = ((hash << 5) + hash) + uint8_t(allele[jj]);
        }

        return hash;
    }

} } // cgatools::variants
