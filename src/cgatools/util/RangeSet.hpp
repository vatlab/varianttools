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

#ifndef CGATOOLS_UTIL_RANGE_SET_HPP_
#define CGATOOLS_UTIL_RANGE_SET_HPP_ 1

//! @file RangeSet.hpp
//! File containing definitions of RangeSet class

#include "cgatools/core.hpp"
#include <vector>
#include <string>
#include <set>

namespace cgatools { namespace reference {
    class Range;
    class Location;
    class CrrFile;
}}

namespace cgatools { namespace util {

    typedef std::vector<std::string> StringVector;

    //! the class is designed to support a relatively small set of ranges
    //! and answer the questions: 
    //! - does a given range intersect one of the ranges in the FastRangeSet?
    //! - is a given location contained by one of the ranges in the FastRangeSet?
    class FastRangeSet 
    {
    public:
        typedef std::pair<uint32_t,uint32_t> Range;

        class RangeSet : public std::set<Range> {
        public:
            bool intersects(const Range& r) const;
            bool intersects(const cgatools::reference::Range& r) const;
            bool contains(uint32_t pos) const;
        };

        typedef std::vector<RangeSet>   ChromosomeRanges;

        FastRangeSet(const reference::CrrFile &ref);
        //! Returns true if there is a range in this set that overlaps the specified
        //! range. Zero-length overlap is allowed for empty ranges,
        //! therefore an empty range overlaps with itself.  An empty
        //! range overlap includes begin of other range: [b,b) overlaps
        //! [a,b) and [b,c) and vise versa.
        bool intersects (const reference::Range& inRange) const;

        //! Returns true if any range of this set contains the given location.
        bool contains(const reference::Location& loc) const;

        void add(const reference::Range& r);

        //! add a string range into the set
        void add(const std::string &rangeStr);
        void add(const StringVector &rangeStrSet);

        //! Fills the range set to cover the whole reference
        //@param extendRangeLength - can be used to extend the chr length beyond the reference length
        void addWholeReference(size_t extendRangeLength=0);

        const ChromosomeRanges& getRanges() const {return ranges_;}

        void clear();
        bool empty() const;

        void regressionTest();
    protected:
        const reference::CrrFile&       reference_;
        ChromosomeRanges                ranges_;
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_STREAMS_HPP_
