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

#ifndef CGA_TOOLS_COMMAND_SNPDIFF_HPP_
#define CGA_TOOLS_COMMAND_SNPDIFF_HPP_ 1

//! @file SnpDiff.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"

namespace cgatools { namespace command {

    class SnpDiff : public Command
    {
    public:
        SnpDiff(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        std::string formatCount(uint32_t count, uint32_t total, bool fraction) const;
        void printMatchRow(std::ostream& out,
                           const char* rowName,
                           std::map<std::string,int>& matchTable,
                           const std::vector< std::pair<std::string,std::string> >& matchClassNames,
                           bool fraction) const;
        void printMatchTable(std::ostream& out,
                             std::vector< std::map<std::string,int> >& stats,
                             const std::vector< std::pair<std::string,std::string> >& matchClassNames,
                             bool fraction) const;
        void printConcordanceRow(std::ostream& out,
                                 const char* rowName,
                                 const std::map<std::string,int>& matchTable,
                                 bool fraction,
                                 bool alleleConcordance) const;
        void printConcordanceTable(std::ostream& out,
                                   const std::vector< std::map<std::string,int> >& stats,
                                   bool fraction,
                                   bool alleleConcordance) const;

        std::string referenceFileName_;
        std::string variantFileName_;
        std::string genotypesFileName_;
        std::string oPrefix_;
        std::string reports_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_SNPDIFF_HPP_
