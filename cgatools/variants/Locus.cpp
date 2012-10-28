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
#include "cgatools/variants/Locus.hpp"
#include "cgatools/variants/Call.hpp"
#include "cgatools/variants/Allele.hpp"

#include <sstream>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <set>

namespace cgatools { namespace variants {

    using namespace std;

    using boost::lexical_cast;

    using util::Exception;

    namespace bu = util::baseutil;
    using namespace cgatools::cgdata;

    LocusAnnotations::LocusAnnotations(
        const std::vector<std::string>& annotationColumnHeaders,
        bool xRefIsAlleleSpecificInit)
        : names_(annotationColumnHeaders),
          xRefIsAlleleSpecific_(xRefIsAlleleSpecificInit)
    {
        // Find allele-specific annotations.
        map< string, pair<size_t,size_t> > asAnn;
        for(size_t ii=0; ii<names_.size(); ii++)
        {
            swapOrder_.push_back(ii);
            if ( !boost::starts_with(names_[ii], "allele1") &&
                 !boost::starts_with(names_[ii], "allele2") )
                continue;

            string nm = names_[ii].substr(7);
            if (asAnn.find(nm) == asAnn.end())
                asAnn[nm] = make_pair(names_.size(), names_.size());
            if (boost::starts_with(names_[ii], "allele1"))
                asAnn[nm].first = ii;
            else
                asAnn[nm].second = ii;
        }

        // If alleles are swapped on output, swapOrder_ tells us which
        // annotations must also be swapped.
        map< string, pair<size_t,size_t> >::const_iterator first, last=asAnn.end();
        for(first=asAnn.begin(); first!=last; ++first)
        {
            if ( first->second.first  != names_.size() &&
                 first->second.second != names_.size() )
            {
                std::swap(swapOrder_[first->second.first], swapOrder_[first->second.second]);
            }
        }
    }

    Locus::Locus(const Locus& other)
        : range_(other.range_),
          crr_(other.crr_),
          calls_(other.calls_),
          alleles_(other.alleles_),
          olplType_(other.olplType_),
          ann_(other.ann_),
          extras_(other.extras_)
    {
        for(size_t ii=0; ii<alleles_.size(); ii++)
            alleles_[ii].setLocus(this);
    }

    Locus& Locus::operator=(const Locus& other)
    {
        if (this != &other)
        {
            range_ = other.range_;
            crr_ = other.crr_;
            calls_ = other.calls_;
            alleles_ = other.alleles_;
            olplType_ = other.olplType_;
            ann_ = other.ann_;
            extras_ = other.extras_;
            for(size_t ii=0; ii<alleles_.size(); ii++)
                alleles_[ii].setLocus(this);
        }
        return *this;
    }

    void Locus::setRange(const Range& range)
    {
        CGA_ASSERT(1 == calls_.size());
        calls_[0].range_ = range;
        range_ = range;
    }

    boost::uint32_t Locus::getId() const
    {
        CGA_ASSERT(calls_.size() > 0);
        return calls_[0].locus_;
    }

    void Locus::setId(boost::uint32_t locusId)
    {
        for(size_t ii=0; ii<calls_.size(); ii++)
            calls_[ii].locus_ = locusId;
    }

    boost::uint16_t Locus::getPloidy() const
    {
        CGA_ASSERT(calls_.size() > 0);
        return calls_[0].ploidy_;
    }

    void Locus::setHapLink(size_t alleleOffset, const std::string& hapLink)
    {
        const std::vector<size_t>& callOffsets = alleles_[alleleOffset].getCallOffsets();
        for(size_t ii=0; ii<callOffsets.size(); ii++)
            calls_[callOffsets[ii]].hapLink_ = hapLink;
    }

    bool Locus::isRefCallLocus() const
    {
        if (calls_.size() != 1)
            return false;

        if ("=" == calls_[0].alleleSeq_)
            return true;

        if ("?" == calls_[0].alleleSeq_)
            return false;

        if ("=" != calls_[0].reference_)
            return calls_[0].alleleSeq_ == calls_[0].reference_;

        return calls_[0].alleleSeq_ == calls_[0].refSequence(*crr_);
    }

    bool Locus::isNoCallLocus() const
    {
        if (calls_.size() != 1)
            return false;

        return "?" == calls_[0].alleleSeq_;
    }

    bool Locus::isRefConsistent() const
    {
        for(size_t ii=0; ii<calls_.size(); ii++)
        {
            if (!calls_[ii].isRefConsistent(*crr_))
                return false;
        }
        return true;
    }

    bool Locus::hasNoCalls() const
    {
        for(size_t ii=0; ii<calls_.size(); ii++)
        {
            if (calls_[ii].hasNoCalls())
                return true;
        }
        return false;
    }

    void Locus::computeReadCounts(int32_t* allele1ReadCount,
                                  int32_t* allele2ReadCount,
                                  int32_t* referenceAlleleReadCount,
                                  int32_t* totalReadCount,
                                  cgdata::EvidenceReader& evidence) const
    {
        *allele1ReadCount = -1;
        *allele2ReadCount = -1;
        *referenceAlleleReadCount = -1;
        *totalReadCount = -1;
        if (isNoCallLocus() || isRefCallLocus())
            return;
        evidence.seek(getRange());
        if (!evidence.inInterval())
            return;

        const EvidenceReader::IntervalRecord& interval = evidence.getInterval();
        const std::vector<EvidenceReader::DnbRecord>& dnbs = evidence.getDnbs();
        if (!hasFullyCalledNonReferenceAllele())
            return;

        bool hasLengthChange;
        EvidenceMatchMap refMatchMap;
        AlleleMatchMap matchMap;
        matchAlleles(interval, hasLengthChange, matchMap, refMatchMap);

        boost::array<size_t, 2> counts;
        counts.assign(0);
        set<string> seenDnbs;
        *referenceAlleleReadCount = 0;
        *totalReadCount = 0;
        BOOST_FOREACH(const EvidenceReader::DnbRecord& dnb, dnbs)
        {
            string dnbid = dnb.getId();
            if (seenDnbs.count(dnbid) != 0)
                continue;
            seenDnbs.insert(dnbid);

            if (!dnb.hasOverlap(getRange(), hasLengthChange))
                continue;

            ++(*totalReadCount);

            for (size_t locAllele = 0; locAllele < matchMap.size(); ++locAllele)
                counts[locAllele] += countSupport(dnb, matchMap[locAllele]);

            *referenceAlleleReadCount += countSupport(dnb, refMatchMap);
        }
        if (hasCount(matchMap, 0))
            *allele1ReadCount = counts[0];
        if (hasCount(matchMap, 1))
            *allele2ReadCount = counts[1];
    }

    bool Locus::hasCount(const AlleleMatchMap& m, size_t alleleIdx) const
    {
        return alleleIdx < m.size() &&
            std::count(m[alleleIdx].begin(), m[alleleIdx].end(), true) != 0;
    }

    size_t Locus::countSupport(
        const EvidenceReader::DnbRecord& dnb, const EvidenceMatchMap& mm) const
    {
        const int32_t SCORE_THRESHOLD = 3;

        int32_t maxCompat = 0, maxIncompat = 0;
        for (size_t intAllele = 0; intAllele < mm.size(); ++intAllele)
        {
            int dnbScore = dnb.scoreAllele_[intAllele];
            if (mm[intAllele])
                maxCompat = std::max<int32_t>(maxCompat, dnbScore);
            else
                maxIncompat = std::max<int32_t>(maxIncompat, dnbScore);
        }
        if (maxCompat >= maxIncompat + SCORE_THRESHOLD)
            return 1;
        else
            return 0;
    }

    // For each allele in the locus, find the compatible alleles in the
    // evidence file interval
    void Locus::matchAlleles(const EvidenceReader::IntervalRecord& interval,
                             bool& hasLengthChange,
                             AlleleMatchMap& m,
                             EvidenceMatchMap& rm) const
    {
        size_t alleleCount = 0;
        size_t maxEvAllele = 0;
        std::set<size_t> allelesWithMissingAlignment;
        for(size_t ii=0; ii<=interval.alleleIndexes_.size(); ii++)
        {
            size_t iiAlleleCount = std::count(interval.alleleIndexes_.begin(),
                                              interval.alleleIndexes_.end(),
                                              ii);

            // When the call is 0;1 we don't put the alignment of reference
            // into the evidence file. When the call is 1;1, we don't put
            // the second alignment of the alt allele. Old versions of the
            // assembler would sometimes align reference non-trivially, or
            // align the alt allele differently for the two variant file
            // alleles.
            if (0 == ii)
            {
                // There better be at least one alt.
                CGA_ASSERT(iiAlleleCount != interval.alleleIndexes_.size());
                if (iiAlleleCount > 0)
                    allelesWithMissingAlignment.insert(ii);
            }
            else
            {
                CGA_ASSERT(iiAlleleCount > 0);
                if (iiAlleleCount > 1)
                    allelesWithMissingAlignment.insert(ii);
            }

            alleleCount += iiAlleleCount;
            if (iiAlleleCount > 0)
                maxEvAllele = ii;

            if (alleleCount == interval.alleleIndexes_.size())
                break;
        }

        CGA_ASSERT(alleleCount == interval.alleleIndexes_.size());

        hasLengthChange = false;
        m.resize(getAlleles().size());
        for (size_t ii = 0; ii < getAlleles().size(); ++ii)
        {
            const Allele& a = getAlleles()[ii];
            if (a.hasNoCalls())
            {
                m[ii].assign(false);
                continue;
            }
            m[ii].assign(true);
            BOOST_FOREACH(size_t ofs, a.getCallOffsets())
            {
                const Call& c = getCalls()[ofs];
                if (c.alleleSeq_.size() != c.reference_.size())
                    hasLengthChange = true;
                if (c.varType_ != "ref")
                    m[ii][0] = false;
                for(size_t evAllele=1; evAllele<=maxEvAllele; evAllele++)
                {
                    if (!interval.isCompatible(evAllele, c, *crr_))
                        m[ii][evAllele] = false;
                }
            }
            for(size_t evAllele=maxEvAllele+1; evAllele<m[ii].size(); evAllele++)
                m[ii][evAllele] = false;
            if ( 0 == std::count(m[ii].begin(), m[ii].end(), true) )
            {
                // If we could find no allele in the evidence file for this
                // variant file allele, assume it's because we are missing
                // the correct alignment string in the evidence file, and
                // attribute variant allele to the evidence allele with a
                // missing alignment.

                if (0 == allelesWithMissingAlignment.size())
                    throw Exception("misalign var");

                CGA_ASSERT(allelesWithMissingAlignment.size() > 0);
                CGA_ASSERT(allelesWithMissingAlignment.size() == 1);
                BOOST_FOREACH(size_t evAllele, allelesWithMissingAlignment)
                {
                    m[ii][evAllele] = true;
                }
            }
        }

        rm.assign(false);
        rm[0] = true;

        string refSeq = crr_->getSequence(getRange());
        // Note that we assume that interesting loci are never preceded
        // or followed by an insertion call, so we demand a strict
        // comparison (passing ignoreAdjacentInsertion = false).
        for(size_t evAllele=1; evAllele<=maxEvAllele; evAllele++)
        {
            if (interval.isCompatible(evAllele, getRange(), refSeq, false))
                rm[evAllele] = true;
        }
    }

    bool Locus::hasFullyCalledNonReferenceAllele() const
    {
        if (Call::UNKNOWN_PLOIDY == getPloidy())
            return false;
        BOOST_FOREACH(const Allele& a, getAlleles())
        {
            if (!a.hasNoCalls() && a.getAlleleSequence() != a.getRefSequence())
                return true;
        }
        return false;
    }

    void Locus::clearCalls()
    {
        calls_.clear();
    }

    void Locus::addCall(const Call& call)
    {
        calls_.push_back(call);
    }

    void Locus::initFromCalls(bool relaxedReferenceValidation)
    {
        CGA_ASSERT(calls_.size() > 0);

        if (relaxedReferenceValidation)
        {
            for(size_t ii=0; ii<calls_.size(); ii++)
            {
                if (calls_[ii].range_.length() > 0 && 0 == calls_[ii].reference_.size())
                    calls_[ii].reference_ = crr_->getSequence(calls_[ii].range_);
            }
        }

        validateCalls(relaxedReferenceValidation);

        alleles_.resize(calls_[0].ploidy_, Allele(this));
        for(size_t ii=0; ii<alleles_.size(); ii++)
        {
            alleles_[ii].clearCalls();

            for(size_t jj=0; jj<calls_.size(); jj++)
            {
                if (Call::ALL_HAPLOTYPES == calls_[jj].haplotype_ ||
                    static_cast<int>(ii) == calls_[jj].haplotype_-1)
                {
                    alleles_[ii].addCallOffset(jj);
                }
            }
        }

        validateAllelesAndInitRange();
    }

    enum OneBaseWalk
    {
        WALK_EOS = 0,
        WALK_OK,
        WALK_INCOMPATIBLE,
        WALK_LENGTH_NOCALL
    };

    void Locus::locationCalls(
        const Location& loc, std::vector< std::pair<char, const Call*> >& calls) const
    {
        calls.clear();
        size_t ploidy = Call::UNKNOWN_PLOIDY == calls_[0].ploidy_ ? 2 : calls_[0].ploidy_;

        if (1 == calls_.size() &&
            Call::ALL_HAPLOTYPES == calls_[0].haplotype_)
        {
            if ("?" == calls_[0].alleleSeq_)
            {
                calls.resize(ploidy, make_pair('N', &calls_[0]));
                return;
            }
            else if ("=" == calls_[0].alleleSeq_)
            {
                calls.resize(ploidy, make_pair(crr_->getBase(loc), &calls_[0]));
                return;
            }
            else if (Call::UNKNOWN_PLOIDY == calls_[0].ploidy_)
            {
                calls.resize(ploidy, make_pair('N', &calls_[0]));
                return;
            }
        }

        calls.resize(alleles_.size(), make_pair('N', static_cast<const Call*>(0)));
        for(size_t ii=0; ii<alleles_.size(); ii++)
        {
            bool foundCall = false;
            const vector<size_t>& offsets = alleles_[ii].getCallOffsets();
            BOOST_FOREACH(size_t offset, offsets)
            {
                const Call& call = calls_[offset];
                if (loc < call.range_.beginLocation() || call.range_.endLocation() <= loc)
                    continue;

                foundCall = true;
                calls[ii].second = &call;

                string refSequence = call.refSequence(*crr_);
                string varSequence = call.calledSequence(*crr_);

                if (varSequence.length() == 0)
                {
                    calls[ii].first = '-';
                }
                else
                {
                    // The base of interest is somewhere in the middle
                    // of the call. Compare to ref to try to find the
                    // base of interest.
                    char lrCall[2];
                    OneBaseWalk lrCalled[2];
                    lrCall[0] = lrCall[1] = 'N';
                    lrCalled[0] = lrCalled[1] = WALK_EOS;

                    // Walk in from the left, looking for
                    // reference-compatible bases.
                    for(size_t jj=0; jj<refSequence.size() && jj<varSequence.size(); jj++)
                    {
                        if (call.range_.begin_ + jj == loc.offset_)
                        {
                            // This is the base in question.
                            if (varSequence[jj] != '?')
                                lrCall[0] = varSequence[jj];
                            lrCalled[0] = WALK_OK;

                            break;
                        }

                        if (varSequence[jj] == '?')
                        {
                            lrCalled[0] = WALK_LENGTH_NOCALL;
                            break;
                        }
                        if (!bu::isConsistent(varSequence[jj], refSequence[jj]))
                        {
                            lrCalled[0] = WALK_INCOMPATIBLE;
                            break;
                        }
                    }

                    // Walk in from the right, looking for
                    // reference-compatible bases.
                    for(size_t jj=0; jj<refSequence.size() && jj<varSequence.size(); jj++)
                    {
                        size_t refPos = refSequence.size() - jj - 1;
                        size_t varPos = varSequence.size() - jj - 1;
                        if (call.range_.end_ - jj - 1 == loc.offset_)
                        {
                            // This is the base in question.
                            if (varSequence[varPos] != '?')
                                lrCall[1] = varSequence[varPos];
                            lrCalled[1] = WALK_OK;

                            break;
                        }

                        if (varSequence[varPos] == '?')
                        {
                            lrCalled[1] = WALK_LENGTH_NOCALL;
                            break;
                        }
                        if (!bu::isConsistent(varSequence[varPos], refSequence[refPos]))
                        {
                            lrCalled[1] = WALK_INCOMPATIBLE;
                            break;
                        }
                    }

                    // If the call walking from the left and
                    // from the right are compatible, call
                    // this base.
                    if (WALK_OK == lrCalled[0] &&
                        WALK_OK == lrCalled[1])
                    {
                        if (lrCall[0] == lrCall[1])
                            calls[ii].first = lrCall[0];
                        else if (lrCall[0] == 'N')
                            calls[ii].first = lrCall[1];
                        else if (lrCall[1] == 'N')
                            calls[ii].first = lrCall[0];
                        else
                            calls[ii].first = '.';
                    }
                    else if (WALK_OK == lrCalled[0])
                    {
                        calls[ii].first = lrCall[0];
                    }
                    else if (WALK_OK == lrCalled[1])
                    {
                        calls[ii].first = lrCall[1];
                    }
                    else
                    {
                        // Walk failed from both
                        // directions. Is this a no-call or
                        // a complex variation?
                        calls[ii].first = '.';
                        if (WALK_LENGTH_NOCALL == lrCalled[0] ||
                            WALK_LENGTH_NOCALL == lrCalled[1])
                            calls[ii].first = 'N';
                        if (WALK_EOS == lrCalled[0] || WALK_EOS == lrCalled[1])
                            calls[ii].first = '-';
                    }
                }

                break;
            }
            CGA_ASSERT(foundCall);
        }
    }

    void Locus::validateCalls(bool relaxedReferenceValidation) const
    {
        if (0 == calls_.size())
            locusError("No calls in locus");
        if (Call::ALL_HAPLOTYPES == calls_[0].haplotype_ && calls_.size() > 1)
            locusError("Locus has a call for all haplotypes and multiple calls");
        for(size_t ii=0; ii<calls_.size(); ii++)
        {
            if (calls_[ii].range_.begin_ > calls_[ii].range_.end_)
                callError("Call begin is greater than end", calls_[ii]);
            if (calls_[ii].reference_ != "=" &&
                crr_->getSequence(calls_[ii].range_) != calls_[ii].reference_)
            {
                if ( (!relaxedReferenceValidation) ||
                     (!bu::isConsistent(crr_->getSequence(calls_[ii].range_), calls_[ii].reference_)) )
                    callError("Call reference sequence does not match reference", calls_[ii]);
            }
            if (calls_[ii].haplotype_ > calls_[ii].ploidy_)
                callError("Call has haplotype greater than ploidy", calls_[ii]);
            if (ii > 0)
            {
                if (Call::ALL_HAPLOTYPES == calls_[ii].haplotype_)
                    locusError("Locus has a call for all haplotypes and multiple calls");
                if (calls_[ii-1].haplotype_ > calls_[ii].haplotype_)
                    locusError("Badly ordered calls");
                if (calls_[ii-1].haplotype_ == calls_[ii].haplotype_ &&
                    calls_[ii-1].range_ >= calls_[ii].range_)
                    locusError("Badly ordered calls");
                if (calls_[ii-1].haplotype_ == calls_[ii].haplotype_ &&
                    calls_[ii-1].range_.endLocation() != calls_[ii].range_.beginLocation())
                    locusError("Haplotype call ranges do not abut");
            }
        }
    }

    void Locus::validateAllelesAndInitRange()
    {
        for(size_t ii=0; ii<alleles_.size(); ii++)
        {
            if (0 == alleles_[ii].getCallOffsets().size())
                locusError("Not all haplotypes of locus have calls");
        }
        for(size_t ii=0; ii<alleles_.size(); ii++)
        {
            const vector<size_t>& offsets = alleles_[ii].getCallOffsets();
            if (ii > 0)
            {
                if (calls_[offsets[0]].range_.beginLocation() != range_.beginLocation() ||
                    calls_[offsets.back()].range_.endLocation() != range_.endLocation())
                    locusError("Not all haplotypes of locus have same range");
            }
            else
            {
                range_ = Range(calls_[offsets[0]].range_.beginLocation(),
                               calls_[offsets.back()].range_.endLocation());
            }
        }
        if (0 == alleles_.size())
        {
            // Special case: ploidy unknown.
            CGA_ASSERT(1 == calls_.size());
            CGA_ASSERT(Call::UNKNOWN_PLOIDY == calls_[0].ploidy_);
            CGA_ASSERT(Call::ALL_HAPLOTYPES == calls_[0].haplotype_);
            CGA_ASSERT(calls_[0].alleleSeq_ == "?");

            range_ = calls_[0].range_;
        }
    }

    void Locus::callError(const std::string& msg, const Call& call) const
    {
        ostringstream out;
        call.write(out, *crr_);
        throw Exception(msg + ":\n" + out.str());
    }

    void Locus::locusError(const std::string& msg) const
    {
        ostringstream out;
        out << *this;
        throw Exception(msg + ":\n" + out.str());
    }

    const std::string VARTYPE_COMPLEX = "complex";

    std::string Locus::getType() const
    {
        if ("" != olplType_)
            return olplType_;

        size_t calledAlleleCount = 0;
        const Call* call = 0;
        BOOST_FOREACH(const Allele& a, alleles_)
        {
            size_t ricc = 0; // ref inconsistent call count
            BOOST_FOREACH(const size_t ofs, a.getCallOffsets())
            {
                const Call& c = calls_[ofs];

                // When reference itself is not fully called the type is "complex"
                if ("=" != c.reference_ && !bu::isCalledSequence(c.reference_))
                    return VARTYPE_COMPLEX;

                if (!c.isRefConsistent(*crr_))
                {
                    ++ricc;
                    if (0 == call)
                        call = &c;
                    else if (call->varType_ != c.varType_ || call->range_ != c.range_)
                        return VARTYPE_COMPLEX;
                }

            }
            if (0 != ricc && a.hasNoCalls())
                return VARTYPE_COMPLEX;
            if (!a.hasNoCalls())
                ++calledAlleleCount;
        }
        if (0 == calledAlleleCount)
        {
            if ("PAR-called-in-X" == calls_[0].varType_ || "no-ref" == calls_[0].varType_)
                return calls_[0].varType_;
            else
                return VARTYPE_COMPLEX;
        }
        if (0 != call)
            return call->varType_;
        else
            return "ref";
    }

    void Locus::setType(const std::string& olplType)
    {
        olplType_ = olplType;
    }

    std::string Locus::getZygosity() const
    {
        if (getPloidy() == Call::UNKNOWN_PLOIDY)
        {
            return "no-call";
        }
        else if (getPloidy() == 1)
        {
            if (alleles_[0].hasNoCalls())
                return "no-call";
            else
                return "hap";
        }
        else
        {
            size_t calledCount = 0;
            BOOST_FOREACH(const Allele& a, alleles_)
            {
                if (!a.hasNoCalls())
                    ++calledCount;
            }
            if (calledCount == 0)
                return "no-call";
            else if (calledCount == 1)
                return "half";
            if (isRefCallLocus())
                return "hom";
            std::string a1 = alleles_[0].getAlleleSequence(),
                        a2 = alleles_[1].getAlleleSequence();
            if (a1 == a2)
                return "hom";
            std::string ref = crr_->getSequence(range_);
            if (a1 == ref || a2 == ref)
                return "het-ref";
            else
                return "het-alt";
        }
    }

    namespace {
        int order(const Locus& loc, const Allele& a)
        {
            if (a.hasNoCalls())
            {
                bool hasCalls = false;
                BOOST_FOREACH(char c, a.getAlleleSequence())
                {
                    if (bu::isValidBase(c)) {
                        hasCalls = true;
                        break;
                    }
                }
                if (hasCalls)
                    return 3;
                else
                    return 4;
            }
            else
            {
                if (a.getAlleleSequence() != a.getRefSequence())
                    return 1;
                else
                    return 2;
            }
        }
    }

    void Locus::reorderAlleles()
    {
        if (getPloidy() < 2 || calls_.size() < 2)
            return;
        if (order(*this, alleles_[0]) > order(*this, alleles_[1]))
        {
            std::swap(alleles_[0], alleles_[1]);

            // Also need to re-order allele-specific annotations.
            CGA_ASSERT(0 != ann_.get());
            CGA_ASSERT(extras_.size() == ann_->count());
            for(size_t ii=0; ii<extras_.size(); ii++)
            {
                size_t ss = ann_->getSwapOrder(ii);
                if (ss > ii)
                    std::swap(extras_[ss], extras_[ii]);
            }
        }
    }

    void Locus::writeAsOneLine(std::ostream& out, bool writeExtras, char sep) const
    {
        // Reorder alleles
        boost::array<const Allele*, 2> alleles;
        alleles.assign(0);
        for (size_t ii = 0; ii < getPloidy(); ++ii)
            alleles[ii] = &alleles_[ii];
        bool swapped = false;
        if (getPloidy() >= 2 && calls_.size() > 1 &&
                order(*this, alleles_[0]) > order(*this, alleles_[1]))
        {
            std::swap(alleles[0], alleles[1]);
            swapped = true;
        }

        size_t ploidy = std::max(getPloidy(), static_cast<uint16_t>(1));
        out << getId() << sep << ploidy << sep
            << crr_->listChromosomes()[range_.chromosome_].getName() << sep
            << range_.begin_ << sep
            << range_.end_ << sep
            << getZygosity() << sep
            << getType();
        std::string alleleSeq;
        bool printScore = false;
        if (isRefCallLocus()) {
            out << sep << '=';
            alleleSeq = "=";
        }
        else if (isNoCallLocus())
        {
            out << sep << '=';
            alleleSeq = "?";
        }
        else
        {
            out << sep << crr_->getSequence(range_);
            printScore = true;
        }
        // alleles
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << sep;
            if (ii < ploidy)
                out << (alleleSeq.empty() ? alleles[ii]->getAlleleSequence() :
                                            alleleSeq);
        }
        // varScoreVAF
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << sep;
            if (ii < getPloidy() && printScore)
                out << alleles[ii]->getMinVarScoreVAF();
        }
        // varScoreEAF
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << sep;
            if (ii < getPloidy() && printScore)
                out << alleles[ii]->getMinVarScoreEAF();
        }
        // varQuality
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << sep;
            if (ii < getPloidy() && printScore)
                out << alleles[ii]->getMinVarQuality();
        }
        // haplinks
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << sep;
            if (ii < getPloidy())
                out << alleles[ii]->getHapLink();
        }
        // xrefs
        if (xRefIsAlleleSpecific())
        {
            for (size_t ii = 0; ii < 2; ++ii)
            {
                out << sep;
                if (ii < getPloidy())
                    out << alleles[ii]->getXRef();
            }
        }
        else
            out << sep << getAllXRef();

        // extra stuff
        if (writeExtras)
        {
            if (swapped)
            {
                CGA_ASSERT(0 != ann_.get());
                CGA_ASSERT(extras_.size() == ann_->count());
                for(size_t ii=0; ii<extras_.size(); ii++)
                    out << sep << extras_[ann_->getSwapOrder(ii)];
            }
            else
            {
                BOOST_FOREACH(const std::string& ex, extras_)
                {
                    out << sep << ex;
                }
            }
        }
    }

    bool Locus::hasAnnotation(const std::string& name) const
    {
        return ann_->getIndex(name) != ann_->count();
    }

    const std::string& Locus::getAnnotation(const std::string& name) const
    {
        size_t idx = ann_->getIndex(name);
        if (ann_->count() == idx)
            return ann_->emptyString();
        return extras_[idx];
    }

    std::string Locus::getAllXRef() const
    {
        set<string> xrefs;
        BOOST_FOREACH(const Call& c, getCalls())
        {
            if (!c.xRef_.empty())
            {
                string::const_iterator first = c.xRef_.begin();
                size_t pos = 0;
                while (true)
                {
                    size_t next = c.xRef_.find_first_of(';', pos);
                    if (string::npos == next)
                    {
                        xrefs.insert(string(first+pos, c.xRef_.end()));
                        break;
                    }
                    xrefs.insert(string(first+pos, first+next));
                    pos = next+1;
                }
            }
        }
        return boost::join(xrefs, ";");
    }

    bool Locus::xRefIsAlleleSpecific() const
    {
        return ann_->xRefIsAlleleSpecific();
    }

    void Locus::setLocusAnnotations(const LocusAnnotations& ann)
    {
        ann_.reset(new LocusAnnotations(ann));
        extras_.resize(ann_->count());
    }

    /* static */
    void Locus::writeOneLineFileHeader(
        std::ostream& out, bool alleleSpecificXRef, char sep)
    {
        out << ">locus" << sep << "ploidy"
            << sep << "chromosome" << sep << "begin" << sep << "end"
            << sep << "zygosity" << sep << "varType"
            << sep << "reference" << sep << "allele1Seq" << sep << "allele2Seq"
            << sep << "allele1VarScoreVAF" << sep << "allele2VarScoreVAF"
            << sep << "allele1VarScoreEAF" << sep << "allele2VarScoreEAF"
            << sep << "allele1VarQuality" << sep << "allele2VarQuality"
            << sep << "allele1HapLink" << sep << "allele2HapLink";
        if (alleleSpecificXRef)
            out << sep << "allele1XRef" << sep << "allele2XRef";
        else
            out << sep << "XRef";
    }

    std::ostream& operator<<(std::ostream& out, const Locus& locus)
    {
        const vector<Call>& calls = locus.getCalls();
        for(size_t ii=0; ii<calls.size(); ii++)
        {
            calls[ii].write(out, locus.getReference());
            out << "\n";
        }
        return out;
    }


} } // cgatools::variants
