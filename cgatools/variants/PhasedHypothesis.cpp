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
#include "cgatools/variants/PhasedHypothesis.hpp"
#include "cgatools/variants/CallDiffResult.hpp"
#include "cgatools/variants/Superlocus.hpp"

#include <boost/foreach.hpp>

namespace cgatools { namespace variants {

    using namespace std;
    namespace bu = cgatools::util::baseutil;

    using reference::Range;

    PhasedHypothesis::PhasedHypothesis()
    {
    }

    PhasedHypothesis::PhasedHypothesis(const reference::Range& range, size_t ploidy)
        : range_(range),
          alleles_(ploidy)
    {
    }

    void PhasedHypothesis::swap(PhasedHypothesis& other)
    {
        std::swap(range_, other.range_);
        alleles_.swap(other.alleles_);
    }

    void PhasedHypothesis::findBestDiff(const Superlocus& sl,
                                        size_t maxHypothesisCount,
                                        const std::vector<PhasedHypothesis>& lhs,
                                        const std::vector<PhasedHypothesis>& rhs,
                                        const reference::CrrFile& crr,
                                        CallDiffResult& result)
    {
        int bestScore = -1;
        vector<CallDiffResult> bestResults;

        if (0 == lhs.size() || 0 == rhs.size() ||
            lhs[0].size() != rhs[0].size())
        {
            // Ploidy mismatch.
            result.matchTypeBySegment_.clear();
            result.matchType_.clear();
            result.matchType_.push_back(cdmt::PLOIDY_MISMATCH);

            for(size_t ii=0; ii<2; ii++)
            {
                result.hyp_[ii] = PhasedHypothesis();
                result.callClass_[ii].clear();
                BOOST_FOREACH(const Locus& locus, sl.getLoci(ii))
                {
                    const vector<Call>& calls = locus.getCalls();
                    const vector<Allele>& alleles = locus.getAlleles();
                    for(size_t jj=0; jj<alleles.size(); jj++)
                    {
                        if (result.callClass_[ii].size() <= jj)
                            result.callClass_[ii].resize(jj+1);

                        const vector<size_t> callOffsets = alleles[jj].getCallOffsets();
                        for(size_t kk=0; kk<callOffsets.size(); kk++)
                        {
                            result.callClass_[ii][jj].push_back(
                                make_pair(cdmt::PLOIDY_MISMATCH, &calls[callOffsets[kk]]));
                        }
                    }
                }
            }

            return;
        }

        CGA_ASSERT(lhs[0].getRange() == rhs[0].getRange());

        vector<size_t> perm(lhs[0].size());
        for(size_t ii=0; ii<perm.size(); ii++)
            perm[ii] = ii;

        string refSequence = crr.getSequence(lhs[0].getRange());

        do
        {
            for(size_t ii=0; ii<lhs.size(); ii++)
            {
                for(size_t jj=0; jj<rhs.size(); jj++)
                {
                    using std::count;
                    CallDiffResult dr;
                    diffHypotheses(lhs[ii], rhs[jj], refSequence, perm, dr.matchType_);
                    int countRefIdentical =
                        count(dr.matchType_.begin(), dr.matchType_.end(), cdmt::REF_IDENTICAL);
                    int countAltIdentical =
                        count(dr.matchType_.begin(), dr.matchType_.end(), cdmt::ALT_IDENTICAL);
                    int countRefConsistent =
                        count(dr.matchType_.begin(), dr.matchType_.end(), cdmt::REF_CONSISTENT);
                    int countAltConsistent =
                        count(dr.matchType_.begin(), dr.matchType_.end(), cdmt::ALT_CONSISTENT);
                    int score = 1000 * countAltConsistent + 1010 * countAltIdentical +
                        1001 * countRefConsistent + 1012 * countRefIdentical;
                    if ( score > bestScore )
                    {
                        dr.setHypotheses(lhs[ii], rhs[jj], perm);
                        bestResults.clear();
                        bestResults.push_back(dr);
                        bestScore = score;
                    }
                    else if (score == bestScore)
                    {
                        dr.setHypotheses(lhs[ii], rhs[jj], perm);
                        bestResults.push_back(dr);
                    }
                }
            }
        } while (std::next_permutation(perm.begin(), perm.end()));

        // Determine call-level classifications, and use them to
        // make the final decision about which phasing is best.
        bestScore = -1;
        bool hasMismatch = false;
        for(size_t ii=0; ii<bestResults.size(); ii++)
        {
            CallDiffResult& dr = bestResults[ii];
            int countRefIdentical = 0;
            int countRefConsistent = 0;
            int countAltIdentical = 0;
            int countAltConsistent = 0;
            int countOther = 0;
            dr.callClass_[0].resize(dr.matchType_.size());
            dr.callClass_[1].resize(dr.matchType_.size());
            dr.matchTypeBySegment_.resize(dr.matchType_.size());
            for(size_t jj=0; jj<dr.matchType_.size(); jj++)
            {
                classifyCalls(dr.hyp_[0][jj], dr.hyp_[1][jj], refSequence, crr,
                              lhs[0].getRange(),
                              dr.matchType_[jj], dr.callClass_[0][jj], dr.callClass_[1][jj],
                              dr.matchTypeBySegment_[jj]);
                for(size_t mm=0; mm<dr.callClass_.size(); mm++)
                {
                    for(size_t kk=0; kk<dr.callClass_[mm][jj].size(); kk++)
                    {
                        if (cdmt::REF_IDENTICAL == dr.callClass_[mm][jj][kk].first)
                            countRefIdentical++;
                        else if (cdmt::REF_CONSISTENT == dr.callClass_[mm][jj][kk].first)
                            countRefConsistent++;
                        else if (cdmt::ALT_IDENTICAL == dr.callClass_[mm][jj][kk].first)
                            countAltIdentical++;
                        else if (cdmt::ALT_CONSISTENT == dr.callClass_[mm][jj][kk].first)
                            countAltConsistent++;
                        else
                            countOther++;
                    }
                }
            }
            int score = 1000 * countAltConsistent + 1010 * countAltIdentical +
                1001 * countRefConsistent + 1012 * countRefIdentical;
            if (score > bestScore)
            {
                bestScore = score;
                result = bestResults[ii];
                hasMismatch = countOther > 0;
            }
        }

        // Look for phase-mismatch.
        if (hasMismatch)
        {
            vector< vector<PhasedHypothesis> > hypotheses2;
            sl.buildPhasedHypotheses(hypotheses2, maxHypothesisCount, false);
            const std::vector<PhasedHypothesis>& lhs2 = hypotheses2[0];
            const std::vector<PhasedHypothesis>& rhs2 = hypotheses2[1];

            CGA_ASSERT ( 0 != lhs2.size() && 0 != rhs2.size() &&
                         lhs2[0].size() == rhs2[0].size() );

            for(size_t ii=0; ii<perm.size(); ii++)
                perm[ii] = ii;

            bool isPhaseMismatch = false;
            vector<cdmt::MatchType> mt2;
            do
            {
                for(size_t ii=0; ii<lhs2.size() && !isPhaseMismatch; ii++)
                {
                    for(size_t jj=0; jj<rhs2.size() && !isPhaseMismatch; jj++)
                    {
                        mt2.clear();
                        diffHypotheses(lhs2[ii], rhs2[jj], refSequence, perm, mt2);

                        bool mm = false;
                        for(size_t kk=0; kk<mt2.size(); kk++)
                        {
                            if (mt2[kk] >= cdmt::ONLY_A)
                            {
                                mm = true;
                                break;
                            }
                        }
                        if (!mm)
                            isPhaseMismatch = true;
                    }
                }
            } while (std::next_permutation(perm.begin(), perm.end()) && !isPhaseMismatch);

            if (isPhaseMismatch)
            {
                for(size_t ii=0; ii<result.matchType_.size(); ii++)
                {
                    if (result.matchType_[ii] >= cdmt::ONLY_A)
                        result.matchType_[ii] = cdmt::PHASE_MISMATCH;
                }
                for(size_t ii=0; ii<result.callClass_.size(); ii++)
                {
                    for(size_t jj=0; jj<result.callClass_[ii].size(); jj++)
                    {
                        for(size_t kk=0; kk<result.callClass_[ii][jj].size(); kk++)
                        {
                            if (result.callClass_[ii][jj][kk].first >= cdmt::ONLY_A)
                                result.callClass_[ii][jj][kk].first = cdmt::PHASE_MISMATCH;
                        }
                    }
                }
            }
        }

        CGA_ASSERT(bestResults.size() > 0);
    }

    void PhasedHypothesis::testVariant(
        const Superlocus& sl,
        const std::vector<PhasedHypothesis>& hyps,
        const reference::Range& varRange,
        const std::string& varAllele,
        const reference::CrrFile& crr,
        std::string& result)
    {
        result.clear();

        if (0 == hyps.size())
        {
            // Ploidy mismatch.
            result = "N";
            return;
        }

        CGA_ASSERT(hyps[0].getRange().beginLocation() <= varRange.beginLocation());
        CGA_ASSERT(hyps[0].getRange().endLocation() >= varRange.endLocation());

        PhasedAllele rhs;

        Locus lRefCall1(crr), lVarCall(crr), lRefCall2(crr);
        Call refCall1, varCall, refCall2;
        if (hyps[0].getRange().beginLocation() < varRange.beginLocation())
        {
            refCall1.locus_ = 1;
            refCall1.ploidy_ = hyps[0].size();
            refCall1.haplotype_ = Call::ALL_HAPLOTYPES;
            refCall1.range_ = Range(hyps[0].getRange().beginLocation(),
                                   varRange.beginLocation());
            refCall1.varType_ = "ref";
            refCall1.reference_ = "=";
            refCall1.alleleSeq_ = "=";
            lRefCall1.addCall(refCall1);
            lRefCall1.initFromCalls();
            rhs.addCall(refCall1, lRefCall1, crr);
        }
        varCall.locus_ = 2;
        varCall.ploidy_ = hyps[0].size();
        varCall.haplotype_ = Call::ALL_HAPLOTYPES;
        varCall.range_ = varRange;
        varCall.reference_ = crr.getSequence(varRange);
        varCall.alleleSeq_ = varAllele;
        lVarCall.addCall(varCall);
        lVarCall.initFromCalls();
        rhs.addCall(varCall, lVarCall, crr);
        if (varRange.endLocation() < hyps[0].getRange().endLocation())
        {
            refCall2.locus_ = 3;
            refCall2.ploidy_ = hyps[0].size();
            refCall2.haplotype_ = Call::ALL_HAPLOTYPES;
            refCall2.range_ = Range(varRange.endLocation(),
                                   hyps[0].getRange().endLocation());
            refCall2.varType_ = "ref";
            refCall2.reference_ = "=";
            refCall2.alleleSeq_ = "=";
            lRefCall2.addCall(refCall2);
            lRefCall2.initFromCalls();
            rhs.addCall(refCall2, lRefCall2, crr);
        }

        Range extendedVariantRange = varRange;
        if ( (0 == varRange.length() && 0 != varAllele.size()) ||
             (0 == varAllele.size() && 0 != varRange.length()) )
        {
            string matchAllele = (0 == varRange.length()) ? varAllele : crr.getSequence(varRange);
            extendedVariantRange.begin_ = Superlocus::extendLeftBySuffixMatching(
                crr, extendedVariantRange.beginLocation(), matchAllele);
            extendedVariantRange.end_ = Superlocus::extendRightByPrefixMatching(
                crr, extendedVariantRange.endLocation(), matchAllele);
            if (extendedVariantRange.begin_ > 0)
                extendedVariantRange.begin_--;
            if (extendedVariantRange.end_ < crr.listChromosomes()[extendedVariantRange.chromosome_].length())
                extendedVariantRange.end_++;
            if (extendedVariantRange.begin_ < hyps[0].getRange().begin_)
                extendedVariantRange.begin_ = hyps[0].getRange().begin_;
            if (extendedVariantRange.end_ > hyps[0].getRange().end_)
                extendedVariantRange.end_ = hyps[0].getRange().end_;
        }
        CGA_ASSERT(hyps[0].getRange().beginLocation() <= extendedVariantRange.beginLocation());
        CGA_ASSERT(hyps[0].getRange().endLocation() >= extendedVariantRange.endLocation());

        int resultScore = -1;
        string refSequence = crr.getSequence(hyps[0].getRange());
        for(size_t ii=0; ii<hyps.size(); ii++)
        {
            string tmp;
            testVariantOneHypothesis(
                hyps[ii], rhs, varRange, extendedVariantRange, refSequence, tmp);
            int score = 0;
            for(size_t jj=0; jj<tmp.size(); jj++)
            {
                switch(tmp[jj])
                {
                case '1':
                    score += 1000;
                    break;
                case '0':
                    score += 100;
                    break;
                case 'N':
                    score += 10;
                    break;
                case 'V':
                    score += 1;
                    break;
                default:
                    break;
                }
            }
            if (score > resultScore)
            {
                result = tmp;
                std::sort(result.begin(), result.end());
                resultScore = score;
            }
        }
    }

    void PhasedHypothesis::testVariantOneHypothesis(
        const PhasedHypothesis& hyp,
        const PhasedAllele& rhs,
        const reference::Range& variantRange,
        const reference::Range& extendedVariantRange,
        const std::string& refSequence,
        std::string& result)
    {
        result.clear();
        for(size_t ii=0; ii<hyp.size(); ii++)
        {
            const PhasedAllele& lhs = hyp[ii];
            vector<AlleleDiffSegment> segs;
            splitIntoSegments(lhs, rhs, refSequence, hyp.getRange(), segs);
            CGA_ASSERT(segs.size() > 0);

            // Find minimal set of segments overlapping extendedVariantRange.
            int lSeg = 0, rSeg = segs.size()-1;
            for(size_t ii=0; ii<segs.size(); ii++)
            {
                if ( segs[ii].ref_.second > extendedVariantRange.begin_ ||
                     (segs[ii].ref_.second == extendedVariantRange.begin_ &&
                      segs[ii].ref_.second == segs[ii].ref_.first) )
                {
                    lSeg = ii;
                    break;
                }
            }
            for(size_t ii=0; ii<segs.size(); ii++)
            {
                if ( segs[segs.size()-ii-1].ref_.first < extendedVariantRange.end_ ||
                     (segs[segs.size()-ii-1].ref_.first == extendedVariantRange.end_ &&
                      segs[segs.size()-ii-1].ref_.first == segs[segs.size()-ii-1].ref_.second) )
                {
                    rSeg = segs.size()-ii-1;
                    break;
                }
            }

            CGA_ASSERT(lSeg <= rSeg);
            CGA_ASSERT(rSeg < int(segs.size()));

            bool called =
                bu::isCalledSequence(lhs.allele(), segs[lSeg].lhs_.first, segs[rSeg].lhs_.second);
            if ( (segs[lSeg].ref_.first < extendedVariantRange.begin_ ||
                  segs[rSeg].ref_.second > extendedVariantRange.end_) && !called)
            {
                result.push_back('N');
                continue;
            }

            char val = '0';
            for(int ii=0; ii<=lSeg; ii++)
            {
                for(int jj=segs.size()-1; jj>=rSeg; jj--)
                {
                    char tmp = testPhasedAllele(lhs, rhs, segs, ii, jj);
                    if ('1' == tmp)
                        val = '1';
                    else if ('N' == tmp && '1' != val)
                        val = 'N';
                }
            }

            result.push_back(val);
        }
    }

    char PhasedHypothesis::testPhasedAllele(const PhasedAllele& lhs,
                                            const PhasedAllele& rhs,
                                            const std::vector<AlleleDiffSegment>& segs,
                                            int lSeg,
                                            int rSeg)
    {
        bool compat = bu::isConsistent(
            lhs.allele(), segs[lSeg].lhs_.first, segs[rSeg].lhs_.second,
            rhs.allele(), segs[lSeg].rhs_.first, segs[rSeg].rhs_.second);
        if (compat)
        {
            bool called =
                bu::isCalledSequence(lhs.allele(), segs[lSeg].lhs_.first, segs[rSeg].lhs_.second);
            if (!called)
                return 'N';

            return '1';
        }

        return '0';
    }

    void PhasedHypothesis::classifyCalls(
        const PhasedAllele& lhs,
        const PhasedAllele& rhs,
        const std::string& refSequence,
        const reference::CrrFile& crr,
        const reference::Range& range,
        cdmt::MatchType mt,
        std::vector< std::pair<cdmt::MatchType, const Call*> >& lhsCallClass,
        std::vector< std::pair<cdmt::MatchType, const Call*> >& rhsCallClass,
        cdmt::MatchType& mtBySeg)
    {
        lhsCallClass.clear();
        rhsCallClass.clear();

        vector<AlleleDiffSegment> segs;
        splitIntoSegments(lhs, rhs, refSequence, range, segs);

        // Compare each segment.
        mtBySeg = cdmt::REF_IDENTICAL;
        BOOST_FOREACH(AlleleDiffSegment& seg, segs)
        {
            cdmt::MatchType mtSeg;
            bool lhsRefMatch = bu::isConsistent(
                lhs.allele(), seg.lhs_.first, seg.lhs_.second,
                refSequence, seg.ref_.first-range.begin_, seg.ref_.second-range.begin_);
            bool rhsRefMatch = bu::isConsistent(
                rhs.allele(), seg.rhs_.first, seg.rhs_.second,
                refSequence, seg.ref_.first-range.begin_, seg.ref_.second-range.begin_);
            if (!bu::isConsistent(lhs.allele(), seg.lhs_.first, seg.lhs_.second,
                                  rhs.allele(), seg.rhs_.first, seg.rhs_.second))
            {
                if (lhsRefMatch && rhsRefMatch)
                {
                    // This case happens, for example, if fileA is A,
                    // fileB is C, and reference is N.
                    mtSeg = cdmt::MISMATCH;
                }
                else if (lhsRefMatch)
                    mtSeg = cdmt::ONLY_B;
                else if (rhsRefMatch)
                    mtSeg = cdmt::ONLY_A;
                else
                    mtSeg = cdmt::MISMATCH;
            }
            else if ( (!bu::isCalledSequence(lhs.allele(), seg.lhs_.first, seg.lhs_.second)) ||
                      (!bu::isCalledSequence(rhs.allele(), seg.rhs_.first, seg.rhs_.second)) )
            {
                if (lhsRefMatch && rhsRefMatch)
                    mtSeg = cdmt::REF_CONSISTENT;
                else
                    mtSeg = cdmt::ALT_CONSISTENT;                
            }
            else
            {
                CGA_ASSERT(lhsRefMatch == rhsRefMatch);
                if (lhsRefMatch && rhsRefMatch)
                    mtSeg = cdmt::REF_IDENTICAL;
                else
                    mtSeg = cdmt::ALT_IDENTICAL;
            }
            mtBySeg = CallDiffResult::mergeMatchTypes(mtBySeg, mtSeg);
            if (CallDiffResult::isConsistent(mt) && mtSeg > mt)
                mtSeg = mt;
            for(size_t ii=seg.lhsCalls_.first; ii<seg.lhsCalls_.second; ii++)
            {
                if (lhsCallClass.size() > 0 && lhsCallClass.back().second == lhs.calls()[ii])
                    lhsCallClass.back().first =
                        CallDiffResult::mergeMatchTypes(matchTypeForCall(mtSeg, crr,
                                                                         lhsCallClass.back().second),
                                                        lhsCallClass.back().first);
                else
                    lhsCallClass.push_back(make_pair(matchTypeForCall(mtSeg, crr, lhs.calls()[ii]),
                                                     lhs.calls()[ii]));
            }
            for(size_t ii=seg.rhsCalls_.first; ii<seg.rhsCalls_.second; ii++)
            {
                if (rhsCallClass.size() > 0 && rhsCallClass.back().second == rhs.calls()[ii])
                    rhsCallClass.back().first =
                        CallDiffResult::mergeMatchTypes(matchTypeForCall(mtSeg, crr,
                                                                         rhsCallClass.back().second),
                                                        rhsCallClass.back().first);
                else
                    rhsCallClass.push_back(make_pair(matchTypeForCall(mtSeg, crr, rhs.calls()[ii]),
                                                     rhs.calls()[ii]));
            }
        }
    }

    cdmt::MatchType PhasedHypothesis::matchTypeForCall(cdmt::MatchType mt,
                                                       const reference::CrrFile& crr,
                                                       const Call* call)
    {
        if (!CallDiffResult::isConsistent(mt))
            return mt;
        if (cdmt::ALT_CONSISTENT == mt && call->isRefConsistent(crr))
            return cdmt::REF_CONSISTENT;
        if (cdmt::ALT_IDENTICAL == mt && call->isRefConsistent(crr))
            return cdmt::REF_IDENTICAL;
        return mt;
    }

    void PhasedHypothesis::splitIntoSegments(const PhasedAllele& lhs,
                                             const PhasedAllele& rhs,
                                             const std::string& refSequence,
                                             const reference::Range& range,
                                             std::vector<AlleleDiffSegment>& segs)
    {
        segs.clear();
        uint32_t ii=0, jj=0, refPos = range.begin_, lhsPos = 0, rhsPos = 0;
        for(; ii<lhs.calls().size() || jj<rhs.calls().size(); )
        {
            uint32_t iiNext = ii, jjNext = jj;
            uint32_t refPosNext = refPos, lhsPosNext = lhsPos, rhsPosNext = rhsPos;
            uint32_t iiEnd = iiNext, jjEnd = jjNext;
            uint32_t lhsPosEnd = lhsPosNext, rhsPosEnd = rhsPosNext;

            // Eat all 0-length calls into same segment.
            while (iiNext+1 < lhs.calls().size() &&
                   0 == lhs.calls()[iiNext]->range_.length() &&
                   0 == lhs.calls()[iiNext+1]->range_.length())
            {
                iiNext++;
                iiEnd = iiNext;
                lhsPosEnd = lhsPosNext = lhs.pos()[iiNext];
            }
            while (jjNext+1 < rhs.calls().size() &&
                   0 == rhs.calls()[jjNext]->range_.length() &&
                   0 == rhs.calls()[jjNext+1]->range_.length())
            {
                jjNext++;
                jjEnd = jjNext;
                rhsPosEnd = rhsPosNext = rhs.pos()[jjNext];
            }

            for(;;)
            {
                if (iiNext >= lhs.calls().size())
                {
                    CGA_ASSERT(jjNext < rhs.calls().size());
                    CGA_ASSERT(0 == rhs.calls()[jjNext]->range_.length());
                    refPosNext = rhs.calls().back()->range_.end_;
                    jjEnd = jjNext = rhs.calls().size();
                    rhsPosEnd = rhsPosNext = rhs.pos()[jjNext];
                    break;
                }
                if (jjNext >= rhs.calls().size())
                {
                    CGA_ASSERT(iiNext < lhs.calls().size());
                    CGA_ASSERT(0 == lhs.calls()[ii]->range_.length());
                    refPosNext = lhs.calls().back()->range_.end_;
                    iiEnd = iiNext = lhs.calls().size();
                    lhsPosEnd = lhsPosNext = lhs.pos()[iiNext];
                    break;
                }
                if (lhs.calls()[iiNext]->range_.end_ >= range.end_ &&
                    rhs.calls()[jjNext]->range_.end_ >= range.end_)
                {
                    refPosNext = range.end_;
                    iiNext++;
                    jjNext++;
                    lhsPosEnd = lhsPosNext = lhs.pos()[iiNext];
                    rhsPosEnd = rhsPosNext = rhs.pos()[jjNext];
                    iiEnd = iiNext;
                    jjEnd = jjNext;
                    break;
                }
                if (lhs.calls()[iiNext]->range_.end_ == rhs.calls()[jjNext]->range_.end_)
                {
                    refPosNext = lhs.calls()[iiNext]->range_.end_;
                    iiNext++;
                    jjNext++;
                    lhsPosEnd = lhsPosNext = lhs.pos()[iiNext];
                    rhsPosEnd = rhsPosNext = rhs.pos()[jjNext];
                    iiEnd = iiNext;
                    jjEnd = jjNext;
                    break;
                }
                if (lhs.calls()[iiNext]->range_.end_ < rhs.calls()[jjNext]->range_.end_)
                {
                    if (segCanSplit(*rhs.calls()[jjNext], lhs.calls()[iiNext]->range_.end_))
                    {
                        if ("?" == rhs.calls()[jjNext]->alleleSeq_)
                        {
                            rhsPosEnd = rhs.pos()[jjNext+1];
                        }
                        else
                        {
                            rhsPosEnd = rhsPosNext = rhs.pos()[jjNext] + lhs.calls()[iiNext]->range_.end_ -
                                std::max(rhs.calls()[jjNext]->range_.begin_, range.begin_);
                        }
                        refPosNext = lhs.calls()[iiNext]->range_.end_;
                        iiNext++;
                        iiEnd = iiNext;
                        jjEnd = jjNext+1;
                        lhsPosEnd = lhsPosNext = lhs.pos()[iiNext];
                        break;
                    }
                    iiNext++;
                    iiEnd = iiNext;
                    continue;
                }
                if (lhs.calls()[iiNext]->range_.end_ > rhs.calls()[jjNext]->range_.end_)
                {
                    if (segCanSplit(*lhs.calls()[iiNext], rhs.calls()[jjNext]->range_.end_))
                    {
                        if ("?" == lhs.calls()[iiNext]->alleleSeq_)
                        {
                            lhsPosEnd = lhs.pos()[iiNext+1];
                        }
                        else
                        {
                            lhsPosEnd = lhsPosNext = lhs.pos()[iiNext] + rhs.calls()[jjNext]->range_.end_ -
                                std::max(lhs.calls()[iiNext]->range_.begin_, range.begin_);
                        }
                        refPosNext = rhs.calls()[jjNext]->range_.end_;
                        jjNext++;
                        jjEnd = jjNext;
                        iiEnd = iiNext+1;
                        rhsPosEnd = rhsPosNext = rhs.pos()[jjNext];
                        break;
                    }
                    jjNext++;
                    jjEnd = jjNext;
                    continue;
                }
                CGA_ASSERT(false);
            }
            segs.push_back(AlleleDiffSegment(make_pair(lhsPos, lhsPosEnd),
                                             make_pair(rhsPos, rhsPosEnd),
                                             make_pair(refPos, refPosNext),
                                             make_pair(ii, iiEnd),
                                             make_pair(jj, jjEnd)));
            lhsPos = lhsPosNext;
            rhsPos = rhsPosNext;
            refPos = refPosNext;
            ii = iiNext;
            jj = jjNext;
            if (refPos == range.end_)
            {
                // We may need to end in the middle of a call or calls.
                bool lhsEnd = ii >= lhs.calls().size() || lhs.calls()[ii]->range_.length() > 0;
                bool rhsEnd = jj >= rhs.calls().size() || rhs.calls()[jj]->range_.length() > 0;
                if (lhsEnd && rhsEnd)
                    break;
            }
        }
    }

    bool PhasedHypothesis::segCanSplit(const Call& call, uint32_t pos)
    {
        CGA_ASSERT(call.range_.begin_ <= pos && pos <= call.range_.end_);
        if (call.range_.begin_ == pos)
            return true;
        if (call.range_.end_ == pos)
            return true;
        if ("=" == call.alleleSeq_ || "?" == call.alleleSeq_)
            return true;
        if (string::npos != call.alleleSeq_.find('?'))
            return false;
        return call.alleleSeq_.size() == call.range_.length();
    }

    void PhasedHypothesis::diffHypotheses(const PhasedHypothesis& lhs,
                                          const PhasedHypothesis& rhs,
                                          const std::string& refSequence,
                                          const std::vector<size_t>& perm,
                                          std::vector<cdmt::MatchType>& result)
    {
        result.resize(perm.size());
        for(size_t ii=0; ii<perm.size(); ii++)
        {
            bool lhsRefMatch = bu::isConsistent(lhs[perm[ii]].allele(), refSequence);
            bool rhsRefMatch = bu::isConsistent(rhs[ii].allele(), refSequence);
            if (!bu::isConsistent(lhs[perm[ii]].allele(), rhs[ii].allele()))
            {
                if (lhsRefMatch && rhsRefMatch)
                {
                    // This case happens, for example, if fileA is A,
                    // fileB is C, and reference is N.
                    result[ii] = cdmt::MISMATCH;
                }
                else if (lhsRefMatch)
                    result[ii] = cdmt::ONLY_B;
                else if (rhsRefMatch)
                    result[ii] = cdmt::ONLY_A;
                else
                    result[ii] = cdmt::MISMATCH;
                continue;
            }

            if (bu::isCalledSequence(lhs[perm[ii]].allele()) &&
                bu::isCalledSequence(rhs[ii].allele()))
            {
                CGA_ASSERT(lhsRefMatch == rhsRefMatch);
                if (lhsRefMatch && rhsRefMatch)
                    result[ii] = cdmt::REF_IDENTICAL;
                else
                    result[ii] = cdmt::ALT_IDENTICAL;
                continue;
            }

            if (lhsRefMatch && rhsRefMatch)
                result[ii] = cdmt::REF_CONSISTENT;
            else
                result[ii] = cdmt::ALT_CONSISTENT;                
        }
    }

} } // cgatools::variants
