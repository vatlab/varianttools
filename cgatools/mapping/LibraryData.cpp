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
#include "LibraryData.hpp"
#include "GapsEstimator.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/util/DelimitedFile.hpp"

#include <boost/foreach.hpp>

namespace cgatools { namespace mapping {

    LibraryMetadataContainer::LibraryMetadataContainer( 
        const cgdata::GenomeMetadata &gm, bool loadDnbStructureOnly )
    {
        std::vector<std::string> libraryNames = gm.getLibraryNames();
        BOOST_FOREACH(const std::string& lib, libraryNames)
            libraryData_.insert(lib,new LibraryData(gm.getLibraryMetadata(lib),loadDnbStructureOnly));
        loadLane2LibraryData(gm);
    }

    const LibraryData& 
        LibraryMetadataContainer::getLibraryData( const std::string& laneName ) const
    {
        Lane2Library::const_iterator it = lane2libraryData_.find(laneName);
        if (it==lane2libraryData_.end())
            throw util::Exception("Library for the lane not found. Lane: "+laneName);
        Libraries::const_iterator itLib = libraryData_.find(it->second);
        if (itLib==libraryData_.end())
            throw util::Exception("Library data not found: "+it->second);
        return *itLib->second;
    }

    void LibraryMetadataContainer::loadLane2LibraryData( const cgdata::GenomeMetadata &gm )
    {
        std::vector<std::string> laneNames = gm.getLaneNames();
        BOOST_FOREACH(const std::string& lane, laneNames)
            lane2libraryData_[lane] = gm.getLibraryName(lane);
    }


    LibraryData::LibraryData( const cgdata::LibraryMetadata& lm, bool loadDnbStructureOnly )
    {
        boost::shared_ptr<std::istream> libFileStream = util::InputStream::
            openCompressedInputStreamByExtension(lm.getDnbStructureFileName());
        util::DelimitedFile libFile(*libFileStream, lm.getDnbStructureFileName());
        dnbStructure_.init(libFile);
        if (!loadDnbStructureOnly)
        {
            for (uint8_t i=0; i<gapsEstimators_.size(); ++i)
            {
                gapsEstimators_[i].reset(new mapping::GapEst::GapsEstimator(dnbStructure_,i,0));
                cgdata::LibraryMetadata::SequenceDependentGapFileNames gapFiles = 
                    lm.getSequenceGapDistributionFileName(i);
                gapsEstimators_[i]->loadGaps(gapFiles[0],gapFiles[1]);
            }

            mateGapTable_.reset(
                new cgdata::MateGapTable(lm.getMainGapDistributionFileName()/*,0.001*/));
        }
    }



} } // cgatools::mapping
