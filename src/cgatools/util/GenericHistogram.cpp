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

#include "cgatools/core.hpp"
#include "GenericHistogram.hpp"
#include <fstream>

namespace cgatools { namespace util {

    std::ostream& operator<<(std::ostream& out, const SimpleHistogram& src) {
        src.write(out);
        return out;
    }
    void SimpleHistogram::write( std::ostream& out ) const
    {
        out << "#sum," << sum_ << ",overall," << number_ << std::endl;
        out << std::endl;
        out << ">bucket,count" << std::endl;
        for (size_t ii = 0; ii < count_.size()-1; ++ii) {
            out << ii << ',' << count_[ii] << std::endl;
        }
        out << "over," << count_.back() << std::endl;
    }
} }
