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

#ifndef CGATOOLS_VARIANTS_SUPERLOCUSITERATOR_HPP_
#define CGATOOLS_VARIANTS_SUPERLOCUSITERATOR_HPP_ 1

//! @file SuperlocusIterator.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/Superlocus.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"

#include <deque>

namespace cgatools { namespace variants {

    class SuperlocusQueueIterator;

    //! A class to iterate through one or more variant files at a time,
    //! enabling access to a stream of Superlocus. Example:
    //! @code
    //! VariantFileIterator locItA(crr), locItB(crr);
    //! locItA.open(variantFileNameA_);
    //! locItB.open(variantFileNameB_);
    //! SuperlocusIterator slIt;
    //! for(slIt.setVariantFiles(locItA, locItB); !slIt.eof(); ++slIt)
    //! {
    //!     const Superlocus& sl = *slIt;
    //!     Range range = sl.getRange();
    //! }
    //! @endcode
    class SuperlocusIterator : private boost::noncopyable
    {
    public:
        enum IterationFlags
        {
            //! Only include superloci with at least one inconsistent call.
            ITER_INCONSISTENT_ONLY = 1,

            //! Use circular prefix/suffix matching to extend the
            //! superlocus boundaries to include the possible locations
            //! of equivalent indels.
            ITER_CIRCULAR_MATCHING = 2
        };

        //! Creates an empty SuperlocusIterator, with the given
        //! parameter values.
        //! @param extend3Mers Always extend superloci to the right and
        //!                    left by the given number of distinct
        //!                    reference 3-mers.
        //! @param extendBases Always extend superloci to the right and
        //!                    left by at least the given number of bases.
        SuperlocusIterator(
            uint32_t extend3Mers = 4,
            uint32_t extendBases = 0,
            int iterationFlags = ITER_INCONSISTENT_ONLY | ITER_CIRCULAR_MATCHING);

        //! Initialize the SuperlocusIterator to process one file.
        void setVariantFile(VariantFileIterator& iter);

        //! Initialize the SuperlocusIterator to process two files.
        void setVariantFiles(VariantFileIterator& iterA, VariantFileIterator& iterB);

        //! Initialize the SuperlocusIterator to process any number of files.
        void setVariantFiles(const std::vector<VariantFileIterator*> iters);

        //! Skip to the next variant in a stream. The variants must
        //! stream in reference order.
        void skipToVariant(const reference::Range& range, const std::string& alleleSeq);

        //! Do the initial seek, in the case where you're iterating via ++.
        void seekFirst();

        //! Returns true if this SuperlocusIterator is at end of file.
        bool eof() const
        {
            return eof_;
        }

        //! Returns the Superlocus referred to by this iterator.
        const Superlocus& operator*() const
        {
            return sl_;
        }

        //! Returns the Superlocus referred to by this iterator.
        const Superlocus* operator->() const
        {
            return &sl_;
        }

        //! Returns the beginning of pure reference call upstream of the
        //! current superlocus, or a location at or greater than the start
        //! of the current superlocus if the call upstream of the superlocus
        //! is not a pure reference.
        reference::Location getPrecedingRefStart() const
        {
            return precedingRefStart_;
        }

        //! Moves to the next Superlocus.
        SuperlocusIterator& operator++()
        {
            next();
            return *this;
        }

    private:
        void next();
        bool findVariant(
            const reference::Location& loc, bool trackRefStart,
            reference::Range& range, reference::Range& searchRange);
        void retireLoci();
        void extendVariant();
        void extendVariant(reference::Range& slRange, reference::Range& searchRange);

        reference::Range extendSearchRange(
            const reference::Range& range, const Locus& locus) const;
        uint32_t extendLeftBySuffixMatching(const Locus& locus) const;
        uint32_t extendRightByPrefixMatching(const Locus& locus) const;
        uint32_t extendLeftBySuffixMatching(
            const reference::Location& loc, const std::string& sequence) const;
        uint32_t extendRightByPrefixMatching(
            const reference::Location& loc, const std::string& sequence) const;
        size_t findContig(const reference::Location& loc) const;

        Superlocus sl_;
        std::vector<VariantFileIterator*> iters_;
        std::vector< std::deque<Locus> >& queues_;
        std::vector< size_t > countRetired_;
        std::vector< boost::shared_ptr<SuperlocusQueueIterator> > qIters_;
        bool started_;
        bool eof_;
        uint32_t extend3Mers_;
        uint32_t extendBases_;
        int iterationFlags_;
        uint32_t padBases_;
        reference::Location queueCutoff_;
        reference::Location slRetired_;
        const reference::CrrFile* crr_;
        std::vector<reference::Range> contigs_;
        // If the superlocus is adjacent to a reference call, this
        // location tracks the beginning of that call; if this location
        // is equal or greater than the superlocus range start, then
        // there's a non-reference call to the left of the superlocus.
        reference::Location precedingRefStart_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_SUPERLOCUSITERATOR_HPP_
