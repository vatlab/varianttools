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

#ifndef CGA_TOOLS_COMMAND_EVIDENCE2CACHE_HPP_
#define CGA_TOOLS_COMMAND_EVIDENCE2CACHE_HPP_

//! @file Evidence2Cache.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/util/RangeSet.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "boost/scoped_ptr.hpp"

namespace cgatools { namespace command {

    class Evidence2CacheConfig {
    public:
        Evidence2CacheConfig() 
            :verboseOutput_(false)
        {}

        void init();

        bool   addAlleleId_;
        bool   verboseOutput_;

        boost::scoped_ptr<util::FastRangeSet>  
                            exportRegions_;
        util::StringVector  exportRegionList_;
        reference::CrrFile  reference_;

        std::string referenceFileName_;
        std::string outputDirName_;
        std::string exportRootDirName_;
    };

    class Evidence2Cache : public Command
    {
    public:
        Evidence2Cache(const std::string& name);

    protected:
        int run(po::variables_map& vm);

    private:
        Evidence2CacheConfig config_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_EVIDENCE2CACHE_HPP_
