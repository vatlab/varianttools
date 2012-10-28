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


// This file contains implementations of the classes defined in
// JunctionStat.hpp.


// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/JunctionStat.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"

namespace cgatools { namespace junctions {

    #undef ADD_STAT_PARAM
    #define ADD_STAT_PARAM(p) #p,
    const JunctionCounterFields::Names JunctionCounterFields::names_ = {
        { ADD_STAT_PARAMS() }
    };


    void JuctionStatisticsCollector::findOneSideNeighbours(const Junctions& junctions, int side)
    {
        if (junctions.size()<2)
            return;
        JunctionRefs junctionIts(junctions);
        std::sort(junctionIts.begin(),junctionIts.end(),CompareBySide(side));
        for(size_t prevI=0, i=prevI; ++i<junctionIts.size();++prevI) {
            const reference::Location& prevLoc = CompareBySide::getLocation(*junctionIts[prevI],side);
            const reference::Location& loc = CompareBySide::getLocation(*junctionIts[i],side);
            uint32_t dist = std::abs(prevLoc.distanceTo(loc));
            if (dist<=distanceTolerance_) 
            {
                ++counters_[side==JUNCTION_LEFT_SIDE ? 
                            JunctionCounterFields::LEFT_SIDES_EQUAL:
                            JunctionCounterFields::RIGHT_SIDES_EQUAL
                            ];
            }

            ++oneSideNeighbourDist_[dist];
        }
    }

    void JuctionStatisticsCollector::findNeighbours(const Junctions& junctions, int sideToCompare)
    {
        if (junctions.size()<2)
            return;
        int firstSide = (sideToCompare==JUNCTION_BOTH_SIDES) ? JUNCTION_LEFT_SIDE : sideToCompare;
        int secondSide = JUNCTION_RIGHT_SIDE-firstSide;
        typedef std::vector<Junctions::const_iterator> JunctionIts;
        JunctionIts junctionIts;
        junctionIts.reserve(junctions.size());
        for (Junctions::const_iterator it=junctions.begin(),endIt=junctions.end(); it!=endIt;++it)
            junctionIts.push_back(it);
        std::sort(junctionIts.begin(),junctionIts.end(), CompareBySide(firstSide));

        reference::Location clusterFrom(999,0);
        reference::Location clusterTo(999,0);
        size_t   clusterSize = 1;
        size_t last =0;
        size_t numberOfJunctions = junctionIts.size();
        while (last<numberOfJunctions) {

            int32_t bestDist = -1;
            reference::Location lastLocationSecondSide = 
                CompareBySide::getLocation(*junctionIts[last],secondSide);
            reference::Location lastLocationFirstSide = 
                CompareBySide::getLocation(*junctionIts[last],firstSide);
            for (size_t i=last+1; i<numberOfJunctions; ++i) {
                uint32_t distFirst = std::abs(lastLocationFirstSide.distanceTo(
                    CompareBySide::getLocation(*junctionIts[i],firstSide)));
                if (distFirst>searchDistanceTolerance_)
                    break;
                uint32_t distSecond = std::abs(lastLocationSecondSide.distanceTo(
                    CompareBySide::getLocation(*junctionIts[i],secondSide)));
                uint32_t maxDist = (sideToCompare==JUNCTION_BOTH_SIDES) ? 
                                        std::max(distFirst,distSecond) : distFirst;
                if (bestDist<0)
                    bestDist = maxDist;
                else
                    bestDist = std::min<int32_t>(bestDist,maxDist);
            }

            if (bestDist>=0 && (uint32_t)bestDist<=distanceTolerance_) {
                if (sideToCompare==JUNCTION_BOTH_SIDES)
                    ++counters_[JunctionCounterFields::BOTH_SIDES_EQUAL];
                if (clusterFrom.chromosome_==999) {
                    clusterFrom = lastLocationFirstSide;
                }
                clusterTo = CompareBySide::getLocation(*junctionIts[last+1],firstSide);
                ++clusterSize;
            } else {
                if (clusterFrom.chromosome_!=999) {
                    CGA_ASSERT_EQ(clusterFrom.chromosome_,clusterTo.chromosome_);
                    ++clusterSizeInBases_[sideToCompare][clusterTo.offset_-clusterFrom.offset_];
                    ++clusterSizeInJunctions_[sideToCompare][clusterSize];
                    clusterFrom.chromosome_ = 999;
                    clusterSize=1;
                }
            }

            if (bestDist>=0 && (uint32_t)bestDist<=searchDistanceTolerance_)
                ++twoSideNeighbourDist_[bestDist];
            ++last;
        }
    }

    void JuctionStatisticsCollector::run( const Junctions& junctions, std::ostream &out )
    {
        out_ = &out;
        counters_[JunctionCounterFields::TOTAL_JUNCTIONS] = junctions.size();
        for(size_t i=0; i<junctions.size(); i++) {
            const Junction& j = junctions[i];
            const JunctionSideSection& left = j.sideSections_[JUNCTION_LEFT_SIDE];
            const JunctionSideSection& right = j.sideSections_[JUNCTION_RIGHT_SIDE];
            ++scoreHistogram_[size_t(j.score_)];
            if(left.strand_==JUNCTION_MINUS_STRAND) {
                if(right.strand_==JUNCTION_MINUS_STRAND)
                    ++counters_[JunctionCounterFields::STRAND_BOTH_NEGATIVE];
                else
                    ++counters_[JunctionCounterFields::STRAND_LEFT_NEGATIVE];
            } else if (right.strand_==JUNCTION_MINUS_STRAND) {
                ++counters_[JunctionCounterFields::STRAND_RIGHT_NEGATIVE];
            } else {
                ++counters_[JunctionCounterFields::STRAND_BOTH_POSITIVE];
            }

            if (left.position_.chromosome_ != right.position_.chromosome_)
                ++counters_[JunctionCounterFields::INTERCHOMOSOMIAL];
        }

        findOneSideNeighbours(junctions, JUNCTION_LEFT_SIDE);
        findOneSideNeighbours(junctions, JUNCTION_RIGHT_SIDE);
        for (int i=JUNCTION_LEFT_SIDE; i<=JUNCTION_BOTH_SIDES; ++i)
            findNeighbours(junctions,i);
        counters_.print(out);
        out << std::endl << "Scores:" << std::endl << scoreHistogram_;
        out << std::endl << "One side distance:" << std::endl << oneSideNeighbourDist_;
        out << std::endl << "Two sides distance:" << std::endl << twoSideNeighbourDist_;
        out << std::endl << "Cluster size in bases(left side):" 
            << std::endl << clusterSizeInBases_[JUNCTION_LEFT_SIDE];
        out << std::endl << "Cluster size in junctions(left side):" 
            << std::endl << clusterSizeInJunctions_[JUNCTION_LEFT_SIDE];
        out << std::endl << "Cluster size in bases(right side):" 
            << std::endl << clusterSizeInBases_[JUNCTION_RIGHT_SIDE];
        out << std::endl << "Cluster size in junctions(right side):" 
            << std::endl << clusterSizeInJunctions_[JUNCTION_RIGHT_SIDE];
        out << std::endl << "Cluster size in bases(both sides):" 
            << std::endl << clusterSizeInBases_[JUNCTION_BOTH_SIDES];
        out << std::endl << "Cluster size in junctions(both sides):" 
            << std::endl << clusterSizeInJunctions_[JUNCTION_BOTH_SIDES];
    }

}}
