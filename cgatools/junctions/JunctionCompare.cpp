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
// JunctionCompare.hpp.


// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/parse.hpp"

// Boost.
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <set>

//#define CGA_JUNCTIONDIFF_CHECK_CONSISTENCY
//#define CGA_JUNCTIONDIFF_FIND_1SIDE_COMPATIBILITY

namespace cgatools { namespace junctions {

    Junctions JunctionRef::nullValue_(1);

    void FullGenomeJunctionComparator::compareBySide( const JunctionsArray& junctions, int firstSide, 
        JunctionRefs &oneSideCompatibleCombined, 
        JunctionRefs &twoSidesCompatibleCombined,
        JunctionRefs &twoSidesCompatible,
        JunctionRefs &twoSidesIncompatible
        ) const
    {
        size_t numberOfFiles = junctions.size();
        CGA_ASSERT_EQ(numberOfFiles,2);

        CompareBySide cmp(firstSide);

        JunctionRefArray jIts(numberOfFiles);
        for (size_t i=0; i<numberOfFiles; ++i) {
            if (junctions[i].empty())
                return;
            jIts[i].addJunctions(junctions[i], i);
            std::sort(jIts[i].begin(),jIts[i].end(),cmp);
        }

        std::set<std::string> matchedJunctionIds;

        //iterators to move along the junction lists
        std::vector<JunctionRefs::const_iterator> its;
        for (size_t i=0; i<numberOfFiles; ++i)
            its.push_back(jIts[i].begin());

        bool iterate = true;
        while (iterate) 
        {

            iterate = false;
            for (size_t fileNo = 0; fileNo<numberOfFiles; ++fileNo) 
            {
                size_t otherFileNo = numberOfFiles-1-fileNo;

                //skip the empty file
                if (its[fileNo]==jIts[fileNo].end())
                    continue;
                iterate = true;

                const JunctionRef &minIt = *its[fileNo];
                std::vector<JunctionRef> minItCompatibleRefs;

                //check all the possible pairs built by the iterator its[fileNo]
                //and all the elements from a corresponding list within a window defined by
                //distanceTolerance_
                if (its[otherFileNo]!=jIts[otherFileNo].end())
                {
                    //move only the iterator that is behind the others
                    if (cmp(**its[otherFileNo],**its[fileNo]))
                        continue;
                    for (JunctionRefs::const_iterator tempIt = its[otherFileNo];
                                            tempIt<jIts[otherFileNo].end(); ++tempIt) 
                    {
                        uint32_t distanceByFirstSide = std::abs(
                            CompareBySide::getLocation(*minIt,firstSide).
                            distanceTo(CompareBySide::getLocation(**tempIt,firstSide)));
                        if (distanceByFirstSide>distanceTolerance_)
                            break;

                        //the current implementation is simplified to the two-file case
                        //please extend this part for generalization to N
                        const JunctionRef &firstIt = fileNo==0 ? minIt : *tempIt;
                        const JunctionRef &secondIt = fileNo==1 ? minIt : *tempIt;


                        double compatibility = getCompatibility(*firstIt,*secondIt);
                        switch (int(compatibility)) {
                            case 0 : 
                                break;
                            case 1 :
                                {
                                    oneSideCompatibleCombined.push_back(firstIt);
                                    oneSideCompatibleCombined.back().matchedJunctionIts_.push_back(secondIt);
                                }
                                break;
                            case 2 : 
                                {
                                    twoSidesCompatibleCombined.push_back(firstIt);
                                    twoSidesCompatibleCombined.back().matchedJunctionIts_.push_back(secondIt);

                                    minItCompatibleRefs.push_back(*tempIt);

                                    //add ids of the matched junction into the temporary buffer
                                    matchedJunctionIds.insert(firstIt.getId());
                                    matchedJunctionIds.insert(secondIt.getId());
                                }
                                break;
                            default:
                                CGA_ASSERT_MSG(false, "Unexpected compatibility result: "<<compatibility);
                        }
                    }
                }

                //mark the junction that doesn't have a match as incompatible
                if (matchedJunctionIds.count(minIt.getId())==0)
                {
                    twoSidesIncompatible.push_back(minIt);
                }
                else
                {
                    twoSidesCompatible.push_back(minIt);
                    twoSidesCompatible.back().matchedJunctionIts_.swap(minItCompatibleRefs);
                }

                //move the processed iterator
                matchedJunctionIds.erase(its[fileNo]->getId());
                ++its[fileNo];
            }
        }
    }

    void FullGenomeJunctionComparator::mergeResults(JunctionRefArray &input, JunctionRefs &output) const 
    {
        //Sort right side selected junctions by the left side in order to merge later
        CompareBySide cmp(JUNCTION_LEFT_SIDE);

        std::sort(input[0].begin(),input[0].end(), cmp);
        std::sort(input[1].begin(),input[1].end(), cmp);

        std::set_union(input[0].begin(),input[0].end(),
            input[1].begin(),input[1].end(),
            std::back_insert_iterator<JunctionRefs>(output),
            cmp);
    }

    class BelowTheScoreThreshold
    {
    public:
        BelowTheScoreThreshold(size_t threshold):threshold_(threshold) {}
        bool operator() (const Junction& j) {return j.score_<threshold_;}
        size_t threshold_;
    };

    void FullGenomeJunctionComparator::prefilterByScore(JunctionsArray& junctionArray) const
    {
        for (size_t i=0; i<junctionArray.size(); ++i) 
        {
            Junctions &junctions = junctionArray[i];
            junctions.erase(
                std::remove_if(junctions.begin(),junctions.end(),BelowTheScoreThreshold(scoreThresholds_[i])),
                junctions.end());
        }
    }

    void FullGenomeJunctionComparator::compare( const JunctionFiles& junctionFiles )
    {
        prefilteredJunctions_.resize(junctionFiles.size());
        for (size_t i=0; i<junctionFiles.size();++i)
            prefilteredJunctions_[i] = junctionFiles[i].junctions_;
        prefilterByScore(prefilteredJunctions_);

        JunctionRefArray oneSideCompatibleCombined(2);
        JunctionRefArray twoSidesCompatibleCombined(2);
        JunctionRefArray twoSidesCompatible(2);
        JunctionRefArray twoSidesIncompatible(2);

        // do compare for 2 (can be extended for N) individuals finding common junctions and differences
        compareBySide(prefilteredJunctions_, JUNCTION_LEFT_SIDE, // JUNCTION_LEFT_SIDE == 0
                        oneSideCompatibleCombined[0], 
                        twoSidesCompatibleCombined[0],
                        twoSidesCompatible[0],
                        twoSidesIncompatible[0]);

        #if defined(CGA_JUNCTIONDIFF_CHECK_CONSISTENCY) || defined(CGA_JUNCTIONDIFF_FIND_1SIDE_COMPATIBILITY)
            compareBySide(prefilteredJunctions_, JUNCTION_RIGHT_SIDE, // JUNCTION_RIGHT_SIDE == 1
                oneSideCompatibleCombined[1], 
                twoSidesCompatibleCombined[1],
                twoSidesIncompatible[1]);

            #ifdef CGA_JUNCTIONDIFF_FIND_1SIDE_COMPATIBILITY
                JunctionRefs oneSideCompatibleRes;
                printList("compatible1Side0", "compatible (one side,left)", oneSideCompatibleCombined[0]);
                printList("compatible1Side1", "compatible (one side,right)", oneSideCompatibleCombined[1]);
                mergeResults(oneSideCompatibleCombined, oneSideCompatibleRes);
                printList("compatibleOneSide", "compatible (one side,merged)", oneSideCompatibleRes);
            #endif

            #ifdef CGA_JUNCTIONDIFF_CHECK_CONSISTENCY
                if (twoSidesCompatibleCombined[0].size()!=twoSidesCompatibleCombined[1].size()) {
                    report_  << "Warning! Sets of compatible don't match." << std::endl;
                    printList("compatibleL", "compatible (Two sides) Left", twoSidesCompatibleCombined[0]);
                    printList("compatibleR", "compatible (Two sides) Right", twoSidesCompatibleCombined[1]);
                }

                if (twoSidesIncompatible[0].size()!=twoSidesIncompatible[1].size()) {
                    report_  << "Warning! Sets of incompatible don't match." << std::endl;
                    printList("incompatibleL", "incompatible (left)", twoSidesIncompatible[0]);
                    printList("incompatibleR", "incompatible (right)", twoSidesIncompatible[1]);
                }
            #endif
        #endif

        // The usage of splitList commands below look overcomplicated and could be
        // replaced by simply sorting out the junctions in compareBySide. It will be probably
        // done when the code settles.

        // print compatible junctions for each file
        compatiblePerFile_.resize(2);
        splitList(twoSidesCompatible[0],compatiblePerFile_);

        // print compatible junctions
        comatibleCombined_.swap(twoSidesCompatibleCombined[0]);

        // print incompatible junctions for each file
        incompatiblePerFile_.resize(2);
        splitList(twoSidesIncompatible[0],incompatiblePerFile_);
    }

    void FullGenomeJunctionComparator::splitList(const JunctionRefs& list, JunctionRefArray& result) const
    {
        BOOST_FOREACH(const JunctionRef &it, list) {
            uint32_t fileNo = it.sourceId_;
            if (result.size()<=fileNo)
                result.resize(fileNo+1);
            result[fileNo].push_back(it);
        }
    }

    void FullGenomeJunctionComparator::printList(
        const std::string& fileName, 
        const JunctionRefs& list,
        const util::DelimitedFile::Metadata& sourceMetadata,
        const std::vector<std::string> &extraHeaderColumns
        ) const
    {
        util::DelimitedFileMetadata metadata(sourceMetadata);
        metadata.initDefaults();

        util::OutputStream out(fileName);

        out << metadata;
        out << JunctionFile::header_;
        BOOST_FOREACH(const std::string& s, extraHeaderColumns)
            out << JunctionFile::SEP << s;
        out << std::endl;

        BOOST_FOREACH(const JunctionRef& it, list) {
            CGA_ASSERT_EQ((*it).annotations_.size(),extraHeaderColumns.size());
            (*it).write(out, reference_, 0);
            BOOST_FOREACH(const std::string& annotation, (*it).annotations_)
                out << JunctionFile::SEP << annotation;
            out << std::endl;
        }
        out.close();
    }

    void FullGenomeJunctionComparator::printMultiList(
        const std::string& fileName, 
        const JunctionRefs& list,
        const util::DelimitedFile::Metadata& sourceMetadata,
        size_t multiListSize
        ) const
    {
        util::DelimitedFileMetadata metadata(sourceMetadata);
        metadata.initDefaults();

        util::OutputStream out(fileName);

        out << metadata;

        // print matched junctions from "multiListSize" files as several lists next to each other
        // the column names in list headers are postfixed by the list index
        if (multiListSize>1) {
            std::vector<std::string> listHeader(multiListSize);
            for (std::string::const_iterator it=JunctionFile::header_.begin()+1,
                itEnd=JunctionFile::header_.end(); ; ++it)
            {
                if (it==itEnd || *it=='\t') {
                    for (size_t i=0; i<multiListSize; ++i)
                        listHeader[i].push_back('0'+i);
                    if (it==itEnd) break;
                }
                for (size_t i=0; i<multiListSize; ++i)
                    listHeader[i].push_back(*it);
            }
            out << '>';
            for (size_t i=0; i<multiListSize; ++i)
            {
                out << listHeader[i];
                if (i<multiListSize-1)
                    out << "\t|\t";
            }
            out << std::endl;
        } else {
            out << JunctionFile::header_ << std::endl;
        }

        BOOST_FOREACH(const JunctionRef& it, list) {
            (*it).write(out, reference_, 0);
            if (multiListSize>1) {
                CGA_ASSERT_EQ(multiListSize-1,it.matchedJunctionIts_.size());
                for (size_t i=0; i<multiListSize-1; ++i) {
                    out << "\t|\t";
                    (*it.matchedJunctionIts_[i]).write(out,reference_, 0);
                }
            }
            out << std::endl;
        }
        out.close();
    }

    //returns:
    // [0,1) - incompatible junctions
    // [1,2) - one side is compatible
    // [2,3) - both sides are compatible
    double FullGenomeJunctionComparator::getCompatibility( const Junction& j0, const Junction& j1 ) const
    {
        double result = 0;
        uint64_t totDistance = 0;
        if (size_t(j0.score_+1E-5) >= scoreThresholds_[0] && 
            size_t(j1.score_+1E-5) >= scoreThresholds_[1]) 
        {
            for (size_t i=size_t(JUNCTION_LEFT_SIDE); i<=size_t(JUNCTION_RIGHT_SIDE); ++i) {
                uint32_t dist = std::abs(CompareBySide::getLocation(j0,i).distanceTo(
                    CompareBySide::getLocation(j1,i)));
                if (dist <= distanceTolerance_
                    && j0.sideSections_[i].strand_ == j1.sideSections_[i].strand_)
                    ++result;
                totDistance+=dist*dist;
            }
            result+= 10/(sqrt((double)totDistance)+20); //put the weight value into [0,1/2]
        }
        return result;
    }

    //-------------------------------- stream output -----------------------------------

    std::ostream& operator<<(std::ostream &ostr, const JunctionRef::It& j) 
    {
        const JunctionSideSection &left = (*j).sideSections_[JUNCTION_LEFT_SIDE];
        const JunctionSideSection &right = (*j).sideSections_[JUNCTION_RIGHT_SIDE];
        return ostr << (*j).id_
            << '|' << (*j).score_
            << '|' << CompareBySide::getLocation((*j),JUNCTION_LEFT_SIDE)
            << "|s:" << left.strand_
            << '|' << CompareBySide::getLocation((*j),JUNCTION_RIGHT_SIDE)
            << "|s:" << right.strand_
            ;
    }

    std::ostream& operator<<(std::ostream &ostr, const JunctionRef& j) 
    {
        ostr << '(';
        ostr << j.it_;
        ostr << ")->{";
        BOOST_FOREACH(const JunctionRef& jit, j.matchedJunctionIts_)
            ostr << jit;
        ostr << '}';
        return ostr;
    }

}}
