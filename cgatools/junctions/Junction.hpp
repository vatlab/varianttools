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

#ifndef CGA_TOOLS_JUNCTION_HPP_
#define CGA_TOOLS_JUNCTION_HPP_

//! @file Junction.hpp
//! This file contains definitions of classes used to manipulate junctions
//! and files containing them.
//! A Junction describes a section of assembled sequence consisting of:
//! - A left section which maps, exactly or approximately,
//!   to a section of the reference genome (either strand).
//! - A transition section..
//! - A right section which maps, exactly or approximately,
//!   to a section of the reference genome (either strand).

// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/DelimitedFile.hpp"

// Boost.
#include <boost/array.hpp>

// Standard library.
#include <string>
#include <vector>
#include <map>


namespace cgatools { namespace junctions {

    // Constants for the two possible sides of a Junction.
    const int JUNCTION_LEFT_SIDE = 0;
    const int JUNCTION_RIGHT_SIDE = 1;
    const int JUNCTION_BOTH_SIDES = 2;

    // The two possible strands of a Junction.
    enum JunctionStrand {
        JUNCTION_MINUS_STRAND = -1,
        JUNCTION_UNKNOWN_STRAND = 0,
        JUNCTION_PLUS_STRAND = +1
    };

    //! Class describing the left or right side of a junction.
    class JunctionSideSection {
    public:

        // Strand of the reference that this junction side maps to.
        JunctionStrand strand_;

        // Position of the base closest to the transition section.
        reference::Location position_;

        // Number of bases in this side of the junction, or -1 if unknown.
        int length_;

        //list of the intersecting repeats
        std::string repeatClassification_;

        //list of the intersecting genes
        std::string genes_;

        //! The default constructor constructs an unknown junction side.
        JunctionSideSection();

        //! Constructor from strand/position/length
        JunctionSideSection(JunctionStrand, const reference::Location& position, int length, 
                    const std::string& repeatClassification,
                    const std::string& genes
                    );

        //! serializes the junction side into a stream
        void write(std::ostream&, const reference::CrrFile&) const;

        //! Returns 1 if the junction side points to increasing reference
        //! offsets and -1 otherwise.
        int getDir ( size_t side ) const;

        //! Returns the reference position closest to the breakpoint from
        //! the given junction side.
        reference::Location getBasePos( size_t side ) const;
    };



    //! Main class describing a single Junction.
    class Junction {
    public:
        std::string id_;

        // The three sections that make up this junction.
        boost::array<JunctionSideSection,2> sideSections_;

        //! The base sequence of the transition region which "bridges"
        //! between the left and right side. Can be empty if there
        //! is no transition section.
        std::string transitionSequence_;

        //! The length of the transition sequence, or the estimated
        //! transition length from paired end analysis.
        size_t transitionLength_;

        //! Flag that indicates whether this junction transition region is known.
        //! If false, the sequence is invalid.
        bool            transitionIsKnown_;

        //! Junction quality score. Corresponds to "discordantMatePairAlignments" column in junctions file
        uint32_t        score_;

        //! Various annotation fields
        std::string     xRef_;
        std::string     deletedTransposableElement_;
        std::string     knownUnderrepresentedRepeat_;
        double          frequencyInBaselineGenomeSet_;

        //! The best estimation of the reference sequence of the junction 
        std::string     assembledSequence_;

        //! Extra columns in the junction file are treated as opaque annotations
        std::vector<std::string> annotations_;

        //! The default constructor constructs an invalid Junction
        //! missing both side sections and the transition section.
        Junction();

        //! constructs a complete junction
        Junction(
            const std::string& junctionId,
            const JunctionSideSection& leftSection,
            const JunctionSideSection& rightSection,
            const std::string& transitionSequence,
            size_t transitionLength,
            bool transitionIsKnown,
            uint32_t score,
            const std::string& xRef,
            const std::string& deletedTransposableElement,
            const std::string& knownUnderrepresentedRepeat,
            double frequencyInBaselineGenomeSet,
            const std::string& assembledSequence
            );

        //! returns location of the given side
        const reference::Location &getLocation(int side) const {
            return sideSections_[side].position_;
        }

        //! returns distance for a canonicalized junction
        size_t getDistance() const;

        //! serializes the junction into a stream
        void write(std::ostream&, const reference::CrrFile&,
                   size_t expectedAnnotationCount) const;
    };

    //! A default container for multiple junctions
    typedef std::vector<Junction> Junctions;

    //! Compares two junctions by a given side
    //! Uses the other side and then junction id as secondary keys
    class CompareJunctionsBySide {
    public:
        CompareJunctionsBySide(int side=JUNCTION_LEFT_SIDE) 
            :side_(side), otherSide_(JUNCTION_RIGHT_SIDE-side)
        { 
            CGA_ASSERT_MSG(side==JUNCTION_LEFT_SIDE || side>=JUNCTION_RIGHT_SIDE, 
                "wrong junction side: " << side);
        }

        static const reference::Location &getLocation(const Junction& j, int side) {
            return j.sideSections_[side].position_;
        }

        bool operator() (const Junction& j1,const Junction& j2) const {
            const reference::Location &f1 = j1.getLocation(side_);
            const reference::Location &f2 = j2.getLocation(side_);
            if (f1!=f2)
                return f1<f2;
            const reference::Location &s1 = j1.getLocation(otherSide_);
            const reference::Location &s2 = j2.getLocation(otherSide_);
            if (s1!=s2)
                return s1<s2;
            return j1.id_ < j2.id_; //compare by id to guarantee the sort order
        }

    protected:
        int side_;
        int otherSide_;
    };

    //! Class describing a Junction file.
    class JunctionFile {
    public:
        static const std::string SEP;

        JunctionFile() {}

        //! Create a populated JunctionFile object in memory by
        //! reading an existing junction file.
        JunctionFile(const std::string& name, const reference::CrrFile&);

        //! Read junctions from a file.
        void read(const std::string& name, const reference::CrrFile& reference);

        //! Write a junction file representing this JunctionFile object.
        void write(const std::string& name, const reference::CrrFile&);

        //! Add a new junction.
        void add(const Junction&);

        //! The junctions in this file.
        Junctions junctions_;

        //! Header line.
        static const std::string header_;

        //! Metadata.
        util::DelimitedFile::Metadata metadata_;

        //! The name of the file the junctions were read from
        std::string fileName_;

        //! List of the names of the extra "annotation" columns
        std::vector<std::string> annotationHeaders_;
    };

    //! A default container for multiple junction files
    typedef std::vector<JunctionFile> JunctionFiles;
}}


#endif //CGA_TOOLS_JUNCTION_HPP_

