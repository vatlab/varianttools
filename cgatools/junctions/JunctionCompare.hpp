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

#ifndef CGA_TOOLS_JUNCTION_COMPARE_HPP_
#define CGA_TOOLS_JUNCTION_COMPARE_HPP_

//! @file JunctionCompare.hpp
//! This file contains definitions of classes used to compare sets of junctions.


// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/util/GenericHistogram.hpp"

// Standard library.
#include <string>
#include <vector>
#include <map>

namespace cgatools { namespace reference {
    class CrrFile;
}}

namespace cgatools { namespace junctions {

    //! A utility class used as a lightweight reference to a junction record and
    //! might contain a list of references to associated junctions
    class JunctionRef 
    {
    public:
        typedef Junctions::const_iterator It;
        JunctionRef() :it_(nullValue_.begin()) {}
        JunctionRef(It it) :it_(it) {}

        const Junction& operator* () const {return *it_;}
        bool operator< (const JunctionRef& other) const
        {
            if (sourceId_ != other.sourceId_)
                return sourceId_ < other.sourceId_;
            else
                return it_->id_ < other.it_->id_;
        }
        bool operator==(const JunctionRef& other) const
        {
            return sourceId_ == other.sourceId_ &&
                   it_ == other.it_;
        }

        It                      it_;
        std::vector<JunctionRef> matchedJunctionIts_;

        //! id of a source (for example: a file)
        uint32_t                sourceId_;

        bool isNull() const {return it_==nullValue_.begin() && matchedJunctionIts_.empty();}
        //! generates a unique id based on a junction id and a source id (file number)
        std::string getId() const {return it_->id_ + "-" + boost::lexical_cast<std::string>(sourceId_);}

    protected:
        static Junctions nullValue_;

        friend std::ostream& operator<<(std::ostream &ostr, const JunctionRef::It& j);
        friend std::ostream& operator<<(std::ostream &ostr, const JunctionRef& j);
    };


    //The input junctions have to be canonicalized
    class CompareBySide : public CompareJunctionsBySide {
    public:
        CompareBySide(int side=JUNCTION_LEFT_SIDE) 
            : CompareJunctionsBySide(side)
        {}
        using CompareJunctionsBySide::operator();

        bool operator() (JunctionRef::It j1,JunctionRef::It j2) const {
            return CompareJunctionsBySide::operator() (*j1,*j2);
        }

        bool operator() (const JunctionRef &j1,const JunctionRef &j2) const {
            return CompareJunctionsBySide::operator() (*j1,*j2);
        }
    };

    class JunctionRefs : public std::vector<JunctionRef> 
    {
    public:
        JunctionRefs() : std::vector<JunctionRef>() {}

        JunctionRefs(const Junctions& junctions, uint32_t id = 0) 
            : std::vector<JunctionRef>()
        {
            addJunctions(junctions, id);
        }

        void addJunctions(const Junctions& junctions, uint32_t id) {
            reserve(size()+junctions.size());
            for (Junctions::const_iterator it=junctions.begin(),endIt=junctions.end(); it!=endIt;++it) {
                push_back(it);
                back().sourceId_ = id;
            }
        }
    };

    class FullGenomeJunctionComparator
    {
    public:
        typedef std::vector<Junctions> JunctionsArray;
        typedef std::vector<JunctionRefs> JunctionRefArray;


        FullGenomeJunctionComparator(
            uint32_t distanceTolerance,
            const std::vector<size_t> &scoreThresholds,
            const reference::CrrFile& reference,
            const std::string &outPrefix
        )
        :   distanceTolerance_(distanceTolerance),
            scoreThresholds_(scoreThresholds),
            outPrefix_(outPrefix),
            reference_(reference)
        {}

        static const char SEPARATOR = '\t';


        void compare(const JunctionFiles& junctionFiles);

        //! prints a junction reference list into a file
        void printList(
            const std::string& fileName, 
            const JunctionRefs& list,
            const util::DelimitedFile::Metadata& sourceMetadata,
            const std::vector<std::string> &extraHeaderColumns
            ) const;

        //! prints a combined junction reference list into a file
        //! @multiListSize - the number of synchronized lists printed. 
        //!                  The additional lists are the matched junctions inside the "list" members
        void printMultiList(const std::string& fileName, 
            const JunctionRefs& list,
            const util::DelimitedFile::Metadata& sourceMetadata,
            size_t multiListSize = 1
            ) const;

        //! returns references on incompatible junctions one list per file
        const JunctionRefArray& getCompatible() const {return compatiblePerFile_;}
        //! returns references on incompatible junctions one list per file
        const JunctionRefArray& getIncompatible() const {return incompatiblePerFile_;}

        //! returns references on compatible junctions
        const JunctionRefs&     getCompatibleCombined()  const  {return comatibleCombined_;}
        //! returns prefiltered junctions
        const JunctionsArray&     getPrefiltered()  const  {return prefilteredJunctions_;}

    protected:
        uint32_t            distanceTolerance_;
        std::vector<size_t> scoreThresholds_;
        std::string         outPrefix_;

        const reference::CrrFile& reference_;

        JunctionsArray      prefilteredJunctions_;

        JunctionRefArray    compatiblePerFile_;
        JunctionRefArray    incompatiblePerFile_;
        JunctionRefs        comatibleCombined_;

        //!Compares two junctions
        //!returns:
        //! [0,1) - incompatible junctions
        //! [1,2) - one side is compatible
        //! [2,3) - both sides are compatible
        double getCompatibility(const Junction& j0, const Junction& j1) const;

        //! compares several junction files based on the given side
        //! returns:
        //! oneSideCompatible - the junctions that have one side compatible 
        //! twoSideCompatible - the junctions that have both sides compatible 
        //! twoSidesIncompatible - the junctions that didn't go into twoSideCompatible
        void compareBySide(const JunctionsArray& junctions, int firstSide,
            JunctionRefs &oneSideCompatibleCombined,
            JunctionRefs &twoSidesCompatibleCombined,
            JunctionRefs &twoSidesCompatible,
            JunctionRefs &twoSidesIncompatible) const;

        //! split junctions by the sourceId
        void splitList(const JunctionRefs& list, JunctionRefArray& result) const;
        void prefilterByScore(JunctionsArray& junctionArray) const; 

        //! used to merge several oneSideCompatible lists
        void mergeResults(JunctionRefArray &input, JunctionRefs &output) const;

    };

}}


#endif //CGA_TOOLS_JUNCTION_COMPARE_HPP_

