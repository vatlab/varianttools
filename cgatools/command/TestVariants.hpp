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

#ifndef CGA_TOOLS_COMMAND_TESTVARIANTS_HPP_
#define CGA_TOOLS_COMMAND_TESTVARIANTS_HPP_ 1

//! @file TestVariants.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"

namespace cgatools { namespace command {

    class TestVariants : public Command
    {
    public:
        TestVariants(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        void writeHeader(std::ostream& out,
                         const std::string& line,
                         const std::vector< boost::shared_ptr<variants::VariantFileIterator> >& vfi);

        std::string referenceFileName_;
        std::string inputFileName_;
        std::string outputFileName_;
        size_t maxHypothesisCount_;

        std::vector<std::string> variantFileNames_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_TESTVARIANTS_HPP_
