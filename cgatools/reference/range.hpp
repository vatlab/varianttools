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

#ifndef CGATOOLS_REFERENCE_RANGE_HPP_
#define CGATOOLS_REFERENCE_RANGE_HPP_ 1

//! @file range.hpp
//! File containing definitions of Location and Range classes.

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"

namespace cgatools { namespace reference {

    //! Class representing a location in a CrrFile. Offsets in the
    //! CrrFile are 0-based.
    class Location
    {
    public:
        Location()
            : chromosome_(0),
              offset_(0)
        {
        }

        Location(uint16_t chromosome, uint32_t offset)
            : chromosome_(chromosome),
              offset_(offset)
        {
        }

        //! Returns the position of the other Location minus the position of the current.
        //! Returns std::numeric_limits<int32_t>::max() if the locations are on different chromosomes
        int32_t distanceTo( const Location& other ) const;

        uint16_t chromosome_; //!< Chromosome id in CrrFile.
        uint32_t offset_;     //!< 0-based offset in CrrFile.
    };

    inline int32_t Location::distanceTo( const Location& other ) const
    {
        if (chromosome_!=other.chromosome_)
            return std::numeric_limits<int32_t>::max();
        return other.offset_-offset_;
    }

    inline bool operator==(const Location& lhs, const Location& rhs)
    {
        return lhs.chromosome_ == rhs.chromosome_ &&
            lhs.offset_ == rhs.offset_;
    }

    inline bool operator!=(const Location& lhs, const Location& rhs)
    {
        return !(lhs == rhs);
    }

    inline bool operator<(const Location& lhs, const Location& rhs)
    {
        if (lhs.chromosome_ != rhs.chromosome_)
            return lhs.chromosome_ < rhs.chromosome_;
        return lhs.offset_ < rhs.offset_;
    }

    inline bool operator<=(const Location& lhs, const Location& rhs)
    {
        return !(rhs < lhs);
    }

    inline bool operator>(const Location& lhs, const Location& rhs)
    {
        return rhs < lhs;
    }

    inline bool operator>=(const Location& lhs, const Location& rhs)
    {
        return !(lhs < rhs);
    }

    inline std::ostream& operator<<(std::ostream& ostr, const Location& l)
    {
        return ostr << l.chromosome_<< "," << l.offset_;
    }

    //! Class representing a range of sequence in a CrrFile. Offsets in
    //! the CrrFile are 0-based.
    class Range
    {
    public:
        Range()
            : chromosome_(0),
              begin_(0),
              end_(0)
        {
        }

        Range(uint16_t chromosome, uint32_t begin, uint32_t end)
            : chromosome_(chromosome),
              begin_(begin),
              end_(end)
        {
        }

        Range(const Location& beginLoc, const Location& endLoc)
            : chromosome_(beginLoc.chromosome_),
              begin_(beginLoc.offset_),
              end_(endLoc.offset_)
        {
            CGA_ASSERT_EQ(beginLoc.chromosome_, endLoc.chromosome_);
        }

        //! Return the Location corresponding to the first base of the
        //! Range.
        Location beginLocation() const
        {
            return Location(chromosome_, begin_);
        }

        //! Return the Location corresponding to one past the last base
        //! of the Range.
        Location endLocation() const
        {
            return Location(chromosome_, end_);
        }

        //! Return the length of the Range.
        uint32_t length() const
        {
            return end_-begin_;
        }

        //! Returns true if this range contains the given location.
        bool contains (const Location& loc) const
        {
            return chromosome_ == loc.chromosome_
                && loc.offset_ >= begin_ && loc.offset_ < end_;
        }

        //! Returns true if this range overlaps the specified
        //! range. Zero-length overlap is allowed for empty ranges,
        //! therefore an empty range overlaps with itself.  An empty
        //! range overlap includes begin of other range: [b,b) overlaps
        //! [a,b) and [b,c) and vise versa.
        bool intersects (const Range& r) const
        {
            return chromosome_ == r.chromosome_
                && ((r.begin_ < end_ && begin_ < r.end_)
                    || begin_ == r.begin_ || end_ == r.end_);
        }

        //! Returns overlapping portion of both ranges. Always returns
        //! an empty range if one of ranges is empty.
        Range overlappingRange (const Range& rr) const
        {
            CGA_ASSERT(intersects(rr));
            return Range(chromosome_, std::max(begin_, rr.begin_), std::min(end_, rr.end_));
        }

        uint16_t chromosome_; //!< Chromosome id in CrrFile.
        uint32_t begin_;      //!< 0-based offset in CrrFile of first
                              //!< base of the Range.
        uint32_t end_;        //!< 0-based offset in CrrFile of one past
                              //!< first base of the Range.
    };

    inline bool operator==(const Range& lhs, const Range& rhs)
    {
        return lhs.chromosome_ == rhs.chromosome_ &&
            lhs.begin_ == rhs.begin_ &&
            lhs.end_ == rhs.end_;
    }

    inline bool operator!=(const Range& lhs, const Range& rhs)
    {
        return !(lhs == rhs);
    }

    inline bool operator<(const Range& lhs, const Range& rhs)
    {
        if (lhs.chromosome_ != rhs.chromosome_)
            return lhs.chromosome_ < rhs.chromosome_;
        if (lhs.begin_ != rhs.begin_)
            return lhs.begin_ < rhs.begin_;
        return lhs.end_ < rhs.end_;
    }

    inline bool operator<=(const Range& lhs, const Range& rhs)
    {
        return !(rhs < lhs);
    }

    inline bool operator>(const Range& lhs, const Range& rhs)
    {
        return rhs < lhs;
    }

    inline bool operator>=(const Range& lhs, const Range& rhs)
    {
        return !(lhs < rhs);
    }

    //! Class that can be used as Overlap template parameter for an IntervalTree
    //! built on reference::Range.
    struct RangeOverlap
    {
        bool operator()(const Range& a, const Range& b) const
        {
            return a.intersects(b);
        }
    };

    //! Class that can be used as GetBoundary template parameter for an IntervalTree
    //! built on reference::Range.
    struct GetRangeBoundary
    {
        Location operator()(const Range& r, size_t side) const
        {
            if (0 == side)
                return r.beginLocation();
            else
                return r.endLocation();
        }
    };

} } // cgatools::reference

#endif // CGATOOLS_REFERENCE_RANGE_HPP_
