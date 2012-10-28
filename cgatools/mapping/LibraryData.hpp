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

#ifndef CGA_TOOLS_COMMAND_LIBRARYDATA_HPP_
#define CGA_TOOLS_COMMAND_LIBRARYDATA_HPP_ 1

//! @file LibraryData.hpp

#include "cgatools/core.hpp"
#include "cgatools/cgdata/Dnb.hpp"
#include "cgatools/cgdata/LibraryReader.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/ptr_container/ptr_map.hpp>
#include <string>

namespace cgatools { namespace mapping { namespace GapEst {
    class GapsEstimator;
}}}

namespace cgatools { namespace cgdata { 
    class GenomeMetadata;
    class LibraryMetadata;
}}

namespace cgatools { namespace mapping {

    class LibraryData
    {
    public:
        LibraryData(const cgdata::LibraryMetadata& lm, bool loadDnbStructureOnly);

        boost::array<boost::shared_ptr<GapEst::GapsEstimator>,2>    gapsEstimators_;
        boost::shared_ptr<cgdata::MateGapTable>                     mateGapTable_;
        cgdata::DnbStructure                                        dnbStructure_;
    };

    class LibraryMetadataContainer 
    {
    public:
        //for some reason ptr_map.insert(..) doesn't accept const std::string& as a key
        //so, I added the const into the declaration
        typedef boost::ptr_map<const std::string,LibraryData> Libraries; 
        typedef std::map<std::string, std::string>      Lane2Library;

        LibraryMetadataContainer(const cgdata::GenomeMetadata &gm, bool loadDnbStructureOnly);
        
        const LibraryData& getLibraryData(const std::string& laneName) const;

        void loadLane2LibraryData(const cgdata::GenomeMetadata &gm);

        Libraries       libraryData_;
        Lane2Library    lane2libraryData_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_COMMAND_LIBRARYDATA_HPP_
