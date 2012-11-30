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
#include "cgatools/util/RangeSet.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/parse.hpp"

#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/range.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace util {

    bool FastRangeSet::RangeSet::intersects( const Range& r ) const
    {
        if (empty())
            return false;
        RangeSet::const_iterator it = lower_bound(r);
        if (it!=end() && it->first < r.second)
            return true;
        return it!=begin() && r.first < (--it)->second;
    }

    bool FastRangeSet::RangeSet::intersects( const cgatools::reference::Range& r ) const
    {
        return intersects(std::make_pair(r.begin_,r.end_));
    }

    bool FastRangeSet::RangeSet::contains( uint32_t pos ) const
    {
        if (empty())
            return false;
        RangeSet::const_iterator it = lower_bound(std::make_pair(pos,pos));
        if (it!=end() && it->first == pos)
            return true;
        return it!=begin() && pos < (--it)->second;
    }

    bool FastRangeSet::intersects( const reference::Range& inRange ) const
    {
        CGA_ASSERT_L(inRange.chromosome_,ranges_.size());
        return ranges_[inRange.chromosome_].intersects(inRange);
    }

    bool FastRangeSet::contains( const reference::Location& loc ) const
    {
        CGA_ASSERT_L(loc.chromosome_,ranges_.size());
        return ranges_[loc.chromosome_].contains(loc.offset_);
    }

    void FastRangeSet::add( const reference::Range& r )
    {
        ranges_[r.chromosome_].insert(std::make_pair( r.begin_, r.end_));
    }

    void FastRangeSet::add( const std::string &rangeStr )
    {
        StringVector parsedStr;
        boost::split(parsedStr,rangeStr,boost::is_any_of(","));
        size_t chr = reference_.getChromosomeId(parsedStr[0]);
        if (parsedStr.size()==1)
        {
            ranges_[chr].insert(
                std::make_pair( 0, reference_.listChromosomes()[chr].length()));
        } else {
            ranges_[chr].insert(std::make_pair( util::parseValue<uint32_t>(parsedStr[1]),
                                                util::parseValue<uint32_t>(parsedStr[2])));
        }
    }

    void FastRangeSet::add( const StringVector &rangeStrSet )
    {
        BOOST_FOREACH(const std::string &st,rangeStrSet)
            add(st);
    }

    FastRangeSet::FastRangeSet( const reference::CrrFile &ref ) 
        : reference_(ref), ranges_(ref.listChromosomes().size())
    {
    }

    void FastRangeSet::regressionTest()
    {
        add("chr1");
        add("chr2,20,100");
        add("chr3,20,100");
        add("chr3,110,200");
        CGA_ASSERT(intersects(reference::Range(0,500,10)));
        CGA_ASSERT(intersects(reference::Range(1,20,30)));
        CGA_ASSERT(!intersects(reference::Range(1,15,17)));
        CGA_ASSERT(intersects(reference::Range(1,15,30)));
        CGA_ASSERT(intersects(reference::Range(1,15,130)));
        CGA_ASSERT(intersects(reference::Range(1,30,40)));
        CGA_ASSERT(intersects(reference::Range(1,30,400)));
        CGA_ASSERT(intersects(reference::Range(2,30,400)));
        CGA_ASSERT(intersects(reference::Range(2,105,400)));
        CGA_ASSERT(intersects(reference::Range(2,105,140)));
        CGA_ASSERT(intersects(reference::Range(2,130,140)));
        CGA_ASSERT(intersects(reference::Range(2,130,400)));
        CGA_ASSERT(!intersects(reference::Range(2,105,107)));
        CGA_ASSERT(!intersects(reference::Range(2,400,407)));

        CGA_ASSERT(contains(reference::Location(0,400)));
        CGA_ASSERT(!contains(reference::Location(1,400)));
        CGA_ASSERT(!contains(reference::Location(2,400)));
        CGA_ASSERT(contains(reference::Location(2,20)));
        CGA_ASSERT(contains(reference::Location(2,25)));
        CGA_ASSERT(contains(reference::Location(2,99)));
        CGA_ASSERT(!contains(reference::Location(2,15)));
        CGA_ASSERT(!contains(reference::Location(2,100)));
        CGA_ASSERT(!contains(reference::Location(2,105)));
        CGA_ASSERT(contains(reference::Location(2,110)));
        CGA_ASSERT(contains(reference::Location(2,111)));
        CGA_ASSERT(!contains(reference::Location(2,200)));
    }

    void FastRangeSet::clear()
    {
        BOOST_FOREACH(RangeSet &s, ranges_)
            s.clear();
    }

    void FastRangeSet::addWholeReference(size_t extendRangeLength)
    {
        for (size_t i=0; i<reference_.listChromosomes().size(); ++i)
            add(reference::Range(i,0,reference_.listChromosomes()[i].length()+extendRangeLength));
    }

    bool FastRangeSet::empty() const
    {
        BOOST_FOREACH(const RangeSet &s, ranges_)
            if (!s.empty())
                return false;
        return true;
    }


} } // cgatools::util
