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

#ifndef CGA_TOOLS_COMMAND_LISTVARIANTS_HPP_
#define CGA_TOOLS_COMMAND_LISTVARIANTS_HPP_ 1

//! @file ListVariants.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/variants/Call.hpp"
#include "cgatools/reference/CrrFile.hpp"

#include <set>

namespace cgatools { namespace command {

    class ListVariants : public Command
    {
    public:
        class CallOrder
        {
        public:
            bool operator()(const cgatools::variants::Call& lhs,
                            const cgatools::variants::Call& rhs) const
            {
                if (lhs.range_ != rhs.range_)
                    return lhs.range_ < rhs.range_;
                return lhs.alleleSeq_ < rhs.alleleSeq_;
            }
        };

        ListVariants(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        void queueCall(
            std::set<cgatools::variants::Call, CallOrder>& callSet,
            cgatools::variants::Call call);

        void retireQueuedCalls(
            std::ostream& out,
            std::set<cgatools::variants::Call, CallOrder>& callSet);

        void purgeQueuedCalls(
            std::ostream& out,
            std::set<cgatools::variants::Call, CallOrder>& callSet);

        void retireCall(
            std::ostream& out,
            const cgatools::variants::Call& call);

        void canonicalizeCall(
            cgatools::variants::Call& call) const;

        std::string rotateLeft(
            const std::string& str,
            int32_t count) const;

        std::string referenceFileName_;
        std::string outputFileName_;

        std::vector<std::string> variantFileNames_;
        std::vector<std::string> variantListingFileNames_;
        bool listLongVariants_;
        uint32_t queueBaseCount_;
        uint32_t variantId_;
        cgatools::variants::Call prevCall_;
        cgatools::reference::CrrFile crr_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_LISTVARIANTS_HPP_
