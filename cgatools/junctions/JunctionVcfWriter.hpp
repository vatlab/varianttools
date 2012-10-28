// Copyright 2012 Complete Genomics, Inc.
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

#ifndef CGA_TOOLS_JUNCTION_VCF_WRITER_HPP_
#define CGA_TOOLS_JUNCTION_VCF_WRITER_HPP_

//! @file JunctionVcfWriter.hpp
//! This file contains procedure for writing junctions to vcf file


// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"

// Standard library.
#include <string>
#include <vector>
#include <map>
#include <set>

// Boost
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace junctions {

//! Map from a junction ID in a single file to all compatible
//! junctions in all other files.
    typedef std::map<std::string, junctions::JunctionRefs> JunctionCompatMap;

//! For each file, the per-ID map to compatible junctions.
typedef std::vector<JunctionCompatMap> JunctionCompatMapPerFile;

//! Data type to keep track of junctions
typedef std::set<JunctionRef> JunctionRefSet;

//! One junction side, sortable in the reference order.
struct JunctionRefSide
{
    JunctionRef jr_;
    size_t side_;

    JunctionRefSide(const JunctionRef& jr, size_t side)
        :   jr_(jr), side_(side)
    {}

    bool operator<(const JunctionRefSide& b) const
    {
        reference::Location
            apos = jr_.it_->sideSections_  [  side_].getBasePos(side_),
            bpos = b.jr_.it_->sideSections_[b.side_].getBasePos(b.side_);
        if (apos != bpos)
            return apos < bpos;
        else
            return side_ < b.side_;
    }
};


class JunctionVcfWriter
{
    public:

        JunctionVcfWriter ( boost::shared_ptr<reference::CrrFile> ref, 
                            boost::shared_ptr<junctions::JunctionFiles> junctionFiles );

        std::string  fileFieldSeparator_;
        size_t  filterScoreThreshold_;
        size_t  filterSideLength_;


        //! Writes VCF headers.
        void writeJunctionVcfHeaders(std::ostream& out) const;

        //! Returns VCF id of the junction side. VCF has a separate line
        //! for each side of the "adjacency", and we generate the ID out of
        //! the junction file ID by appending "L" or "R" depending on the side.
        std::string createJunctionVcfId( const JunctionRef& jref, size_t side) const;


        //! Return genomic position in the VCF format: chromosome name
        //! without "chr" and 1-based offset. Separator between the
        //! chromosome subfield and the position can be specified; if the
        //! separator is empty, the chromosome subfield is not printed at all.
        std::string formatPositionForVcf( reference::Location pos, const std::string& sep ) const;

        void writeJunctionPositionToVcf( const junctions::JunctionSideSection& jss, 
                                         size_t side, std::ostream& out) const;


        //! Writes the contents of the ALT column for a VCF "adjacency" line
        //! that corresponds to the given junction side. This column describes
        //! the orientation of the current side, the transition sequence, and
        //! the orientation of the other side. It consist of the sequence
        //! field SF, transition field TF and the oriented adjacent position
        //! field (OAP). The SF field is the closest base to the junction from
        //! the current side. The TF sequence is in strand orientation required
        //! to convert the current side of the junction to the reference strand.
        //! The OAP is the 1-based position of the other-side's base that's
        //! closest to the junction, surrounded by either [ or ] characters,
        //! depending on whether the adjacent sequence extends right- or
        //! leftward from the specified position. The field order is
        //! SF TF OAP if the current side points right, and OAP TF SF otherwise.
        void writeJunctionAltFieldToVcf(const JunctionRef& jref,
            size_t side, std::ostream& out, bool suppressChrom=false) const;

        //! Parse mobile element deletion information in our format, and
        //! return a string that contains the same information in more VCF-ish
        //! format: comma-separated, no "chr" in chromosome name, 1-based
        //! closed-interval coordinates.
        std::string convertMobileElementToVcf(const std::string& med) const;

        //! Write semicolon-separated list of all INFO subfields. Currently
        //! we write the type (always BND for "breakend"), ID of the other
        //! side of the junction, frequency in the baseline, and xref and
        //! deleted mobile element fields if present.
        void writeJunctionInfoFieldToVcf(const JunctionRef& jref,
            size_t side, std::ostream& out) const;

        //! Write semicolon-separated list of all INFO subfields. Currently
        //! we write the type (always BND for "breakend"), ID of the other
        //! side of the junction, frequency in the baseline, and xref and
        //! deleted mobile element fields if present.
        //! Filter field helper.
        void addFilterFlag( std::ostream& out, const std::string& flag, bool& filtered) const;

        //! Write semicolon-separated list of all filters that this junction
        //! failed or PASS if it passes all filters. The set of filters is the
        //! same as the high-confidence file filters, with the exception of
        //! not removing baseline cross-chr junctions.
        void writeJunctionFilterFieldToVcf( const JunctionRef& jref, std::ostream& out) const;

        //! most of the junction information squeezed into a single field
        void writeJunctionComparisonField( const JunctionRef& jref, size_t side, 
                                           std::ostream& out) const;

        //! Writes one side of a VCF "adjacency" to the given stream.
        //! Note that, regardless of our junction side strand, VCF record
        //! is always written relative to the primary strand.
        void writeJunctionToVcf(  const JunctionRef& jref, size_t side,
                                  const JunctionCompatMapPerFile& compat,
                                  std::ostream& out) const;

protected:

        JunctionVcfWriter();

        //! Initializes internal data such as sample IDs
        void init();

        boost::shared_ptr<reference::CrrFile>       reference_;
        boost::shared_ptr<junctions::JunctionFiles> junctionFiles_;
        std::vector<std::string>                    sampleIds_;
};

}}


#endif //CGA_TOOLS_JUNCTION_VCF_WRITER_HPP_

