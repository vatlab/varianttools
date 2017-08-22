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

#ifndef CGATOOLS_REFERENCE_CHROMOSOMEIDFIELD_HPP_
#define CGATOOLS_REFERENCE_CHROMOSOMEIDFIELD_HPP_ 1

//! @file ChromosomeIdField.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/DelimitedLineParser.hpp"

namespace cgatools { namespace reference {

    class CrrFile;

    class ChromosomeIdField : public util::DelimitedFieldParser
    {
    public:
        ChromosomeIdField(const std::string& name, uint16_t* id, const CrrFile& ref)
            : DelimitedFieldParser(name),
              ref_(ref),
              id_(id)
        {
        }

        void parse(const char* first, const char* last);

    private:
        const CrrFile& ref_;
        uint16_t* id_;
        std::string buf_;
    };

} } // cgatools::reference

#endif // CGATOOLS_REFERENCE_CHROMOSOMEIDFIELD_HPP_
