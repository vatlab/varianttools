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
#include "cgatools/variants/SuperlocusIterator.hpp"
#include "cgatools/util/BaseUtil.hpp"

#include <queue>
#include <boost/foreach.hpp>

namespace cgatools { namespace variants {

    using namespace std;

    using namespace cgatools::util;
    namespace bu = cgatools::util::baseutil;
    using reference::CompactDnaSequence;
    using reference::Location;
    using reference::Range;

    class SuperlocusQueueIterator
    {
    public:
        SuperlocusQueueIterator(std::deque<Locus>& queue,
                                size_t& countRetired,
                                VariantFileIterator& iter)
            : queue_(&queue),
              countRetired_(&countRetired),
              iter_(&iter),
              idx_(countRetired)
        {
            if (0 == queue_->size() && !iter_->eof())
            {
                queue_->push_back(*(*iter_));
                ++(*iter_);
            }
        }

        bool eof() const
        {
            CGA_ASSERT( (*countRetired_) <= idx_ );
            return idx_-(*countRetired_) == queue_->size() && iter_->eof();
        }

        const Locus& operator*() const
        {
            CGA_ASSERT( (*countRetired_) <= idx_ );
            return (*queue_)[idx_-(*countRetired_)];
        }

        const Locus* operator->() const
        {
            CGA_ASSERT( (*countRetired_) <= idx_ );
            return &(*queue_)[idx_-(*countRetired_)];
        }

        SuperlocusQueueIterator& operator++()
        {
            ++idx_;
            CGA_ASSERT( (*countRetired_) <= idx_ );
            if (idx_-(*countRetired_) == queue_->size() && !iter_->eof())
            {
                queue_->push_back(*(*iter_));
                ++(*iter_);
            }
            return *this;
        }

        size_t getOffsetInQueue() const
        {
            CGA_ASSERT( (*countRetired_) <= idx_ );
            return idx_ - *countRetired_;
        }

        bool operator<(const SuperlocusQueueIterator& rhs) const
        {
            if (rhs.eof())
                return true;
            if (eof())
                return false;
            const Range& lhsRange = (*this)->getRange();
            const Range& rhsRange = rhs->getRange();
            if (rhsRange.beginLocation() != lhsRange.beginLocation())
                return rhsRange.beginLocation() < lhsRange.beginLocation();
            return rhsRange.endLocation() < lhsRange.endLocation();
        }

    private:
        std::deque<Locus>* queue_;
        size_t* countRetired_;
        VariantFileIterator* iter_;
        size_t idx_;
    };

    SuperlocusIterator::SuperlocusIterator(
        uint32_t extend3Mers,
        uint32_t extendBases,
        int iterationFlags)
        : queues_(sl_.queues_),
          countRetired_(sl_.queues_.size(), 0),
          started_(false),
          eof_(false),
          extend3Mers_(extend3Mers),
          extendBases_(extendBases),
          iterationFlags_(iterationFlags),
          padBases_(50)
    {
    }

    void SuperlocusIterator::setVariantFile(VariantFileIterator& iter)
    {
        vector<VariantFileIterator*> iters;
        iters.push_back(&iter);
        setVariantFiles(iters);
    }

    void SuperlocusIterator::setVariantFiles(VariantFileIterator& iterA, VariantFileIterator& iterB)
    {
        vector<VariantFileIterator*> iters;
        iters.push_back(&iterA);
        iters.push_back(&iterB);
        setVariantFiles(iters);
    }

    void SuperlocusIterator::setVariantFiles(const std::vector<VariantFileIterator*> iters)
    {
        if (0 == iters.size())
            throw Exception("failed to construct SuperlocusIterator: no variant files specified");
        if (0 != iters_.size())
            throw Exception("failed to construct SuperlocusIterator: variant files already specified");
        iters_ = iters;
        queues_.resize(iters_.size());
        countRetired_.resize(iters_.size(), 0);
        crr_ = &(*iters[0])->getReference();
        contigs_ = crr_->listContigs();

        // Make sure all the iterators are for the same reference.
        for(size_t ii=1; ii<iters.size(); ii++)
            CGA_ASSERT(crr_ == &(*iters[ii])->getReference());

        sl_.id_ = 1;
    }

    void SuperlocusIterator::seekFirst()
    {
        if (0 == iters_.size())
            throw Exception("failed to iterate: SuperlocusIterator has no variant files");

        if (started_)
            throw Exception("SuperlocusIterator::seekFirst called twice");

        started_ = true;
        precedingRefStart_ = Location();
        if (!findVariant(Location(), true, sl_.range_, sl_.searchRange_))
        {
            eof_ = true;
            return;
        }
        extendVariant();
    }

    void SuperlocusIterator::next()
    {
        if (0 == iters_.size())
            throw Exception("failed to iterate: SuperlocusIterator has no variant files");

        if (!started_)
            throw Exception("failed to iterate: SuperlocusIterator::seekFirst not called");

        sl_.id_++;
        slRetired_ = sl_.range_.endLocation();
        Location loc(sl_.range_.endLocation());
        loc.offset_++;

        if (!findVariant(loc, true, sl_.range_, sl_.searchRange_))
        {
            eof_ = true;
            return;
        }
        extendVariant();
        CGA_ASSERT(slRetired_ < sl_.range_.beginLocation());
        CGA_ASSERT(slRetired_ <= precedingRefStart_);
    }

    void SuperlocusIterator::skipToVariant(const Range& range, const std::string& alleleSeq)
    {
        if (0 == iters_.size())
            throw Exception("failed to iterate: SuperlocusIterator has no variant files");

        if (range.beginLocation() < sl_.range_.beginLocation())
            throw Exception("failed to seek in variant files: input variants out of order");

        sl_.id_++;

        sl_.range_ = range;
        Call call;
        call.locus_ = 1;
        call.ploidy_ = 1;
        call.haplotype_ = 1;
        call.range_ = range;
        call.reference_ = "=";
        call.alleleSeq_ = alleleSeq;
        Locus locus;
        locus.addCall(call);
        locus.initFromCalls();
        sl_.searchRange_ = extendSearchRange(range, locus);
        extendVariant();
    }

    // Finds the first Locus L that is not consistent with the
    // reference, such that L.getRange().endLocation() >= loc. Returns
    // true for success, in which case range is set to L.getRange() and
    // the searchRange is set to the range that must be searched for
    // no-calls and nearby variants.
    // If trackRefStart is true, also adjusts precedingRefStart; this is
    // done only when we are not on the look-ahead codepath.
    bool SuperlocusIterator::findVariant(const Location& loc, bool trackRefStart,
                                         Range& range, Range& searchRange)
    {
        priority_queue< SuperlocusQueueIterator > pq;
        for(size_t ii=0; ii<queues_.size(); ii++)
            pq.push(SuperlocusQueueIterator(queues_[ii], countRetired_[ii], *iters_[ii]));

        while (!pq.empty())
        {
            SuperlocusQueueIterator qit = pq.top();
            pq.pop();

            if (qit.eof())
                continue;

            const Locus& locus = *qit;
            bool isInteresting = false;
            if (iterationFlags_ & ITER_INCONSISTENT_ONLY)
                isInteresting = !locus.isRefConsistent();
            else
                isInteresting = !locus.isRefCallLocus();
            if (locus.getRange().endLocation() >= loc && isInteresting)
            {
                range = locus.getRange();
                searchRange = extendSearchRange(range, locus);
                return true;
            }
            else if (trackRefStart)
            {
                if (locus.isRefCallLocus())
                    precedingRefStart_ =
                        std::max(precedingRefStart_, locus.getRange().beginLocation());
                else
                    precedingRefStart_ =
                        std::max(precedingRefStart_, locus.getRange().endLocation());
            }

            ++qit;
            pq.push(qit);
        }

        return false;
    }

    void SuperlocusIterator::retireLoci()
    {
        queueCutoff_ = sl_.searchRange_.beginLocation();
        if (queueCutoff_.offset_ > padBases_)
            queueCutoff_.offset_ -= padBases_;
        else
            queueCutoff_.offset_ = 0;
        // Refine queueCutoff_ -- make sure it's on a cuttable boundary.
        for(;;)
        {
            bool newCutoff = false;
            for(size_t ii=0; ii<queues_.size(); ii++)
            {
                deque<Locus>::const_iterator first, last=queues_[ii].end();
                for(first=queues_[ii].begin(); first!=last; ++first)
                {
                    if (first->getRange().endLocation() >= queueCutoff_)
                        break;
                }
                if (first != last)
                {
                    if (first->getRange().endLocation() == queueCutoff_ ||
                        first->isRefCallLocus() || first->isNoCallLocus())
                        continue;
                    if (first->getRange().beginLocation() < queueCutoff_)
                    {
                        queueCutoff_ = first->getRange().beginLocation();
                        newCutoff = true;
                    }
                }
                else if (!iters_[ii]->eof())
                {
                    Location nextLoc = (*iters_[ii])->getRange().beginLocation();
                    if (nextLoc < queueCutoff_)
                    {
                        queueCutoff_ = nextLoc;
                        newCutoff = true;
                    }
                }
            }
            if (!newCutoff)
                break;
        }
        for(size_t ii=0; ii<queues_.size(); ii++)
        {
            while (queues_[ii].size() > 1 && queues_[ii].front().getRange().endLocation() < queueCutoff_)
            {
                queues_[ii].pop_front();
                countRetired_[ii]++;
            }
        }
    }

    // Extends sl_.range_ and sl_.searchRange_ to the right and left
    // according to extend3Mers_, extendBases_, and circular
    // prefix/suffix matching. Checks forward one superlocus for
    // superlocus merging.
    void SuperlocusIterator::extendVariant()
    {
        retireLoci();
        extendVariant(sl_.range_, sl_.searchRange_);

        // Check forward one superlocus for superlocus merging.
        for(;;)
        {
            Location loc(sl_.range_.endLocation());
            loc.offset_++;

            Range nextRange, searchRange;
            bool ok = findVariant(loc, false, nextRange, searchRange);
            if (!ok)
                return;

            extendVariant(nextRange, searchRange);
            if (sl_.range_.endLocation() < nextRange.beginLocation())
                return;

            CGA_ASSERT(sl_.range_.chromosome_ == nextRange.chromosome_);
            sl_.range_.end_ = nextRange.end_;
            sl_.searchRange_.end_ = std::max(sl_.range_.end_, sl_.searchRange_.end_);
        }
    }

    // Extends the given range to the right and left according to
    // extend3Mers_, extendBases_, and circular prefix/suffix matching.
    void SuperlocusIterator::extendVariant(Range& slRange, Range& searchRange)
    {
        // Read in more loci, extending pendingRange_ right as far as is
        // necessary to incorporate no-calls and additional nearby
        // variants. The pendingRange_ must be on a Locus boundary,
        // unless the locus is a ref-call locus or a no-call locus.
        priority_queue< SuperlocusQueueIterator > pq;
        for(size_t ii=0; ii<queues_.size(); ii++)
            pq.push(SuperlocusQueueIterator(queues_[ii], countRetired_[ii], *iters_[ii]));

        size_t iterationCount = 0;
        while (!pq.empty())
        {
            SuperlocusQueueIterator qit = pq.top();
            pq.pop();

            if (qit.eof())
                continue;

            const Locus& locus = *qit;
            if (searchRange.endLocation() <= locus.getRange().beginLocation())
                break;

            if (searchRange.beginLocation() < locus.getRange().endLocation() && !locus.isRefCallLocus())
            {
                // We need to update slRange.
                bool isInteresting = false;
                if (iterationFlags_ & ITER_INCONSISTENT_ONLY)
                    isInteresting = !locus.isRefConsistent();
                else
                    isInteresting = !locus.isRefCallLocus();
                if (isInteresting)
                {
                    slRange.end_ = std::max(slRange.end_, locus.getRange().end_);
                    searchRange = extendSearchRange(searchRange, locus);
                }
                else if (locus.isNoCallLocus())
                {
                    uint32_t end = std::min(locus.getRange().begin_+1, locus.getRange().end_);
                    slRange.end_ = std::max(slRange.end_, end);
                    searchRange.end_ = std::max(slRange.end_, searchRange.end_);
                }
                else
                {
                    slRange.end_ = std::max(slRange.end_, locus.getRange().end_);
                    searchRange.end_ = std::max(slRange.end_, searchRange.end_);
                }
            }

            ++qit;
            pq.push(qit);
            ++iterationCount;
            if (0 == iterationCount % 1000)
            {
                retireLoci();
            }
        }

        // Extend slRange left as far as is necessary to incorporate
        // no-calls. The slRange must be on a Locus boundary, unless the
        // locus is a ref-call locus or a no-call locus.
        for(;;)
        {
            bool newBegin = false;
            for(size_t ii=0; ii<queues_.size(); ii++)
            {
                deque<Locus>::const_iterator first, last = queues_[ii].end();
                for(first=queues_[ii].begin(); first!=last; ++first)
                {
                    if (first->getRange().endLocation() <= searchRange.beginLocation())
                        continue;
                    if (first->getRange().beginLocation() >= slRange.beginLocation())
                        break;

                    const Locus& locus = *first;
                    if (!locus.isRefCallLocus())
                    {
                        // We need to update slRange.
                        if ( locus.isNoCallLocus() &&
                             (iterationFlags_ & ITER_INCONSISTENT_ONLY) )
                        {
                            uint32_t begin = locus.getRange().begin_;
                            if (locus.getRange().end_ > 0)
                                begin = std::max(begin, locus.getRange().end_-1);
                            if (begin < slRange.begin_ &&
                                Location(slRange.chromosome_, begin) >= slRetired_)
                            {
                                slRange.begin_ = begin;
                                searchRange.begin_ = std::min(searchRange.begin_, slRange.begin_);
                                CGA_ASSERT(queueCutoff_ <= searchRange.beginLocation());
                                CGA_ASSERT(slRetired_ <= searchRange.beginLocation());
                                newBegin = true;
                            }
                        }
                        else
                        {
                            if (locus.getRange().begin_ < slRange.begin_ &&
                                Location(slRange.chromosome_, locus.getRange().begin_) >= slRetired_)
                            {
                                slRange.begin_ = locus.getRange().begin_;
                                searchRange.begin_ = std::min(searchRange.begin_, slRange.begin_);
                                CGA_ASSERT(queueCutoff_ <= searchRange.beginLocation());
                                CGA_ASSERT(slRetired_ <= searchRange.beginLocation());
                                newBegin = true;
                            }
                        }
                        break;
                    }
                }
            }
            if (!newBegin)
                break;
        }
    }

    // Extends the search range left and right according to superlocus
    // extension rules for the given locus.
    reference::Range SuperlocusIterator::extendSearchRange(
        const reference::Range& range, const Locus& locus) const
    {
        const Range& lRange = locus.getRange();
        CGA_ASSERT(range.chromosome_ == lRange.chromosome_);
        const CompactDnaSequence& chromosome = crr_->listChromosomes()[range.chromosome_];
        size_t contig = findContig(range.beginLocation());

        // Extend locus left by base count, 3-mer count, and suffix
        // matching.
        uint32_t minLocBases = 0;
        if (lRange.begin_ > extendBases_)
            minLocBases = lRange.begin_ - extendBases_;
        uint32_t minLoc3Mer = chromosome.extendLeftBy3Mers(lRange.begin_, extend3Mers_);
        uint32_t minLoc = std::min(minLocBases, minLoc3Mer);
        if (iterationFlags_ & ITER_CIRCULAR_MATCHING)
        {
            uint32_t minLocSuffix = extendLeftBySuffixMatching(locus);
            minLoc = std::min(minLoc, minLocSuffix);
        }
        if (contig < contigs_.size())
            minLoc = std::max(minLoc, contigs_[contig].begin_);
        if (locus.getRange().begin_ >= padBases_)
            minLoc = std::max(minLoc, locus.getRange().begin_-padBases_);
        if (queueCutoff_.chromosome_ == lRange.chromosome_)
            minLoc = std::max(minLoc, queueCutoff_.offset_);
        if (slRetired_.chromosome_ == lRange.chromosome_)
            minLoc = std::max(minLoc, slRetired_.offset_);

        // Extend locus right by base count, 3-mer count, and prefix
        // matching.
        uint32_t maxLocBases = std::min(lRange.end_+extendBases_, uint32_t(chromosome.length()));
        uint32_t maxLoc3Mer = chromosome.extendRightBy3Mers(lRange.end_, extend3Mers_);
        uint32_t maxLoc = std::max(maxLocBases, maxLoc3Mer);
        if (iterationFlags_ & ITER_CIRCULAR_MATCHING)
        {
            uint32_t maxLocSuffix = extendRightByPrefixMatching(locus);
            maxLoc = std::max(maxLoc, maxLocSuffix);
        }
        if (contig < contigs_.size())
            maxLoc = std::min(maxLoc, contigs_[contig].end_);
        maxLoc = std::min(maxLoc, lRange.end_+padBases_);

        return Range(range.chromosome_,
                     std::min(range.begin_, minLoc),
                     std::max(range.end_, maxLoc));
    }

    uint32_t SuperlocusIterator::extendLeftBySuffixMatching(const Locus& locus) const
    {
        uint32_t result = locus.getRange().begin_;
        BOOST_FOREACH(const Call& call, locus.getCalls())
        {
            if (bu::isCalledSequence(call.alleleSeq_))
            {
                result = std::min(result, Superlocus::extendLeftBySuffixMatching(
                                      *crr_, call.range_.beginLocation(), call.alleleSeq_));
                result = std::min(result, Superlocus::extendLeftBySuffixMatching(
                                      *crr_, call.range_.beginLocation(), call.refSequence(*crr_)));
            }
        }
        if (result > 0)
            result--;
        return result;
    }

    uint32_t SuperlocusIterator::extendRightByPrefixMatching(const Locus& locus) const
    {
        uint32_t result = locus.getRange().end_;
        BOOST_FOREACH(const Call& call, locus.getCalls())
        {
            if (bu::isCalledSequence(call.alleleSeq_))
            {
                result = std::max(result, Superlocus::extendRightByPrefixMatching(
                                      *crr_, call.range_.endLocation(), call.alleleSeq_));
                result = std::max(result, Superlocus::extendRightByPrefixMatching(
                                      *crr_, call.range_.endLocation(), call.refSequence(*crr_)));
            }
        }
        if (result < crr_->listChromosomes()[locus.getRange().chromosome_].length())
            result++;
        return result;
    }

    size_t SuperlocusIterator::findContig(const reference::Location& loc) const
    {
        vector<Range>::const_iterator iter =
            std::lower_bound(contigs_.begin(), contigs_.end(), Range(loc,loc));
        if (iter != contigs_.begin())
            --iter;
        while (iter != contigs_.end() && iter->endLocation() < loc)
            ++iter;
        if (iter->beginLocation() <= loc && loc <= iter->endLocation())
            return iter - contigs_.begin();
        return contigs_.size();
    }

} } // cgatools::variants
