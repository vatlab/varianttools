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

#ifndef CGA_TOOLS_COMMAND_JOIN_HPP_
#define CGA_TOOLS_COMMAND_JOIN_HPP_ 1

//! @file Join.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/util/DelimitedFile.hpp"

namespace cgatools { namespace command {

    class Join : public Command
    {
    public:
        Join(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        struct QueryPlan
        {
            // Offsets of fields to be matched.
            std::vector<int> matchIdx_;
            // Offsets of fields for use in overlap comparison.
            std::pair<int,int> overlapIdx_;
        };

        void dumpRecord(std::ostream& out,
                        const std::vector<std::string>& aFields,
                        const std::vector<std::string>& bFields);
        bool overlap(const std::pair<int64_t, int64_t>& lhs,
                     const std::pair<int64_t, int64_t>& rhs);
        void parseJoinFields(const std::vector<std::string>& fields,
                             const QueryPlan& qp,
                             std::string& matchKey,
                             std::pair<int64_t, int64_t>& range);
        void initQueryPlan(util::DelimitedFile& aa, util::DelimitedFile& bb);
        void parseMatchFields(const util::DelimitedFile& df,
                              const std::string& fieldNameList,
                              QueryPlan& qp);
        void parseOverlapFields(const util::DelimitedFile& df,
                                const std::string& fieldNameList,
                                QueryPlan& qp);

        std::vector<std::string> inputFileNames_;
        std::string outputFileName_;

        std::vector<std::string> match_;
        std::string overlapSpec_;
        std::string outputMode_;
        std::string overlapMode_;
        std::string select_;
        bool alwaysDump_;
        double overlapFractionA_, overlapFractionB_;
        int64_t boundaryUncertaintyA_, boundaryUncertaintyB_;

        std::vector<QueryPlan> qp_;
        std::vector<std::string> aFields_, bFields_;
        std::vector< std::pair<size_t, size_t> > transform_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_JOIN_HPP_
