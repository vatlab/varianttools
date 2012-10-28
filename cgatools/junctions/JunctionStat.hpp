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

#ifndef CGA_TOOLS_JUNCTION_STAT_HPP_
#define CGA_TOOLS_JUNCTION_STAT_HPP_

//! @file JunctionStat.hpp
//! This file contains code that collects diferent junction file statistics

// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/util/GenericHistogram.hpp"

// Standard library.
#include <string>
#include <vector>

namespace cgatools { namespace junctions {

    #define ADD_STAT_PARAMS() \
        ADD_STAT_PARAM(TOTAL_JUNCTIONS) \
        ADD_STAT_PARAM(LEFT_SIDES_EQUAL) \
        ADD_STAT_PARAM(RIGHT_SIDES_EQUAL) \
        ADD_STAT_PARAM(BOTH_SIDES_EQUAL) \
        ADD_STAT_PARAM(LEFT_UNKNOWN) \
        ADD_STAT_PARAM(RIGHT_UNKNOWN) \
        ADD_STAT_PARAM(INTERCHOMOSOMIAL) \
        ADD_STAT_PARAM(STRAND_BOTH_POSITIVE) \
        ADD_STAT_PARAM(STRAND_LEFT_NEGATIVE) \
        ADD_STAT_PARAM(STRAND_RIGHT_NEGATIVE) \
        ADD_STAT_PARAM(STRAND_BOTH_NEGATIVE)

    class JunctionCounterFields {
    public:
        enum Fields {
            #undef ADD_STAT_PARAM
            #define ADD_STAT_PARAM(p) p,
            ADD_STAT_PARAMS()
        };

        #undef ADD_STAT_PARAM
        #define ADD_STAT_PARAM(p) +1
        static const size_t NUMBER_OF_FIELDS = ADD_STAT_PARAMS();

        typedef boost::array<std::string,NUMBER_OF_FIELDS> Names;
        static const Names names_;
        static const Names& getFieldNames() {
            return names_;
        }
    };

    template <class Fields, class Type>
    class CounterSet
    {
    public:
        CounterSet() {
            values_.assign(Type(0));
        }

        void print(std::ostream &out) const {
            const typename Fields::Names& fieldNames = Fields::getFieldNames();
            for (size_t i=0; i<values_.size(); ++i)
                out << fieldNames[i] << '\t' << values_[i] << std::endl;
        }

        Type &operator[] (size_t i) {CGA_ASSERT_L(i,values_.size()); return values_[i];}
        boost::array<Type,Fields::NUMBER_OF_FIELDS> values_;
    };

    class JuctionStatisticsCollector 
    {
    public:

        JuctionStatisticsCollector(uint32_t distanceTolerance)
            :distanceTolerance_(distanceTolerance),out_(NULL)
        {
            size_t buckets[7][3] = 
            {
                {0,101,10},
                {200, 1001, 100},
                {2000,10001,1000},
                {20000,   100001,     10000},
                {2000000, 10000001,   1000000},
                {20000000,100000001,  10000000},
                {200000000,1000000001,100000000}
            };
            
            scoreHistogram_.addRange(0,20,1);
            scoreHistogram_.addRange(20,200,10);

            for (size_t i=0; i<6; ++i)
                oneSideNeighbourDist_.addRange(buckets[i][0],buckets[i][1],buckets[i][2]);

            searchDistanceTolerance_ = distanceTolerance*5;
            twoSideNeighbourDist_.addRange(0,searchDistanceTolerance_+10,10);

            for (size_t i=0; i<4; ++i)
                for (size_t j=0; j<clusterSizeInBases_.size(); ++j)
                    clusterSizeInBases_[j].addRange(buckets[i][0],buckets[i][1],buckets[i][2]);

            for (size_t j=0; j<clusterSizeInJunctions_.size(); ++j)
                clusterSizeInJunctions_[j].addRange(0,15,1);
        }

        void findOneSideNeighbours(const Junctions& junctions, int firstSide);
        void findNeighbours(const Junctions& junctions, int sideToCompare=JUNCTION_BOTH_SIDES);
        void run(const Junctions& junctions, std::ostream &out);
    protected:
        uint32_t    distanceTolerance_;
        uint32_t    searchDistanceTolerance_;
        CounterSet<JunctionCounterFields,int64_t> counters_;
        util::GenericHistogram<size_t,size_t> scoreHistogram_;
        util::GenericHistogram<size_t,size_t> oneSideNeighbourDist_;
        util::GenericHistogram<size_t,size_t> twoSideNeighbourDist_;

        boost::array<util::GenericHistogram<size_t,size_t>,3> clusterSizeInBases_;
        boost::array<util::GenericHistogram<size_t,size_t>,3> clusterSizeInJunctions_;
        std::ostream *out_;
    };

}}


#endif //CGA_TOOLS_JUNCTION_STAT_HPP_

