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

#ifndef CGA_TOOLS_COMMAND_DECODECRR_HPP_
#define CGA_TOOLS_COMMAND_DECODECRR_HPP_ 1

//! @file DecodeCrr.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/reference/CrrFile.hpp"

namespace cgatools { namespace command {

    class DecodeCrr : public Command
    {
    public:
        DecodeCrr(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        reference::Range parseRangeParam(
            const reference::CrrFile& crr, const std::string& param) const;

        std::string referenceFileName_;
        std::string outputFileName_;
        std::string range_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_DECODECRR_HPP_
