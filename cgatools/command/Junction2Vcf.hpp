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

#ifndef CGA_TOOLS_COMMAND_JUNCTION2VCF_HPP_
#define CGA_TOOLS_COMMAND_JUNCTION2VCF_HPP_ 1

//! @file Junction2Vcf.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/junctions/Junction.hpp"
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace command {

    class Junction2Vcf : public Command
    {
        public:

            Junction2Vcf(const std::string& name);

        protected:

            int run(po::variables_map& vm);
            void saveVcf ( std::ostream& out );

            std::string referenceFileName_;
            std::string outputFileName_;
            std::string junctionsFileName_;
            size_t scoreThreshold_;
            size_t minJunctionLength_;

            boost::shared_ptr<reference::CrrFile>       reference_;
            boost::shared_ptr<junctions::JunctionFiles> junctionFiles_;
    };
}} // cgatools::command

#endif // CGA_TOOLS_COMMAND_JUNCTION2VCF_HPP_
