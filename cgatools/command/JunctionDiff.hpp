// Copyright 2010-2012 Complete Genomics, Inc.
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

#ifndef CGA_TOOLS_COMMAND_JUNCTIONDIFF_HPP_
#define CGA_TOOLS_COMMAND_JUNCTIONDIFF_HPP_ 1

//! @file JunctionDiff.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"

namespace cgatools { namespace command {

    class JunctionDiffConfig {
    public:
        JunctionDiffConfig() 
            :   calcFileStat_(false),
                junctionFiles_(2),
                scoreThresholds_(2),
                distanceTolerance_(0)
        {}

        bool                        calcFileStat_;
        std::string                 referenceFileName_;
        std::string                 outPrefix_;
        std::vector<std::string>    junctionFiles_;
        std::vector<size_t>         scoreThresholds_;
        uint32_t                    distanceTolerance_;
        uint32_t                    minJunctionLength_;
        bool                        outputVcf_;
        size_t                      vcfFilterScoreThreshold_;
        size_t                      vcfFilterSideLength_;
        bool                        vcfNormalPriorityOutput_;
    };

    class JunctionDiff : public Command
    {
    public:
        JunctionDiff(const std::string& name);

    protected:
        int run(po::variables_map& vm);

    private:
        JunctionDiffConfig config_;
    };
    
} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_JUNCTIONDIFF_HPP_
