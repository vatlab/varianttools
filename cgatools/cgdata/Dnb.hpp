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

#ifndef CGA_TOOLS_DNB_HPP_
#define CGA_TOOLS_DNB_HPP_ 1

//! @file Dnb.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"

#include <vector>
#include <boost/array.hpp>

namespace cgatools { namespace util {
    class DelimitedFile;
}};

namespace cgatools { namespace cgdata {

    class LibRecord {
        friend std::ostream& operator<< (std::ostream& out, const LibRecord& r);
    public:
        LibRecord():id_(0),arm_(0),indInArm_(0),objInArm_(0),minSize_(0),maxSize_(0) {}

        void initParser(util::DelimitedFile &delimitedFile);

        size_t          id_;
        std::string     type_;
        size_t          arm_;
        size_t          indInArm_;
        size_t          objInArm_;
        int             minSize_;
        int             maxSize_;
    };

    class HalfDnbStructure {
    public:
        class Gap {
        public:
            Gap():minSize_(0),maxSize_(0) {}
            int minSize_;
            int maxSize_;
        };
        typedef std::vector<Gap>    Gaps;
        typedef std::vector<size_t> Reads;

        HalfDnbStructure():totReadLength_(0) {}

        Reads getReadsForStrand(bool negativeStrand) const {
            if (negativeStrand) {
                Reads result(reads_);
                std::reverse(result.begin(),result.end());
                return result;
            } else 
                return reads_;
        }

        Gaps    gaps_;
        Reads   reads_;
        size_t  totReadLength_;
    };

    class DnbStructure {
    public:
        typedef std::vector<HalfDnbStructure>   HalfDnbs;

        void init(const std::string &libFile);
        void init(util::DelimitedFile &libFile);

        HalfDnbs                halfDnbs_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_DNB_HPP_
