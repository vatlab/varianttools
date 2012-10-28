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

#ifndef CGA_TOOLS_COMMAND_JUNCTIONRATIONALIZATION_HPP_
#define CGA_TOOLS_COMMAND_JUNCTIONRATIONALIZATION_HPP_ 1

//! @file JunctionType.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"

namespace cgatools { namespace command {

    class JunctionType : public Command
    {
    public:
        JunctionType(const std::string& name);

    protected:
        int run(po::variables_map& vm);

    private:
        std::string referenceFileName_, junctionFileName_,
                    fullJunctionFileName_, outputPrefix_,
                    repmaskFileName_, mobileElementNames_,
                    geneFileName_;

        int maxPairingDistance_, maxCopyTargetLength_, maxDistanceToTransposon_,
            maxSimpleEventDistance_, maxRelatedJunctionDistance_,
            maxContainedGenesRangeLength_, maxRelatedJunctionOutput_;

        uint32_t regulatoryRegionLength_;



    };
    
} } // cgatools::command

#endif // CGA_TOOLS_JUNCTIONRATIONALIZATION_HPP_
