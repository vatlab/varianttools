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
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/command/DecodeCrr.hpp"
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace command {

    using std::endl;

    using std::string;
    using std::vector;

    using reference::CrrFile;
    using reference::Range;

    using util::Exception;

    namespace ba = boost::algorithm;

    DecodeCrr::DecodeCrr(const std::string& name)
        : Command(name,
                  "Prints the reference sequence for a given reference range.",
                  "",
                  "")
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file (may be passed in as argument at the end of the command).")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ("range", po::value<string>(&range_),
             "The range of bases to print (chr,begin,end or chr:begin-end).")
            ;

        positionalOptions_.add("reference", 1);
    }

    reference::Range DecodeCrr::parseRangeParam(
        const reference::CrrFile& crr, const std::string& param) const
    {
        vector<string> tokens;
        if (param.find(':') == string::npos)
        {
            ba::split(tokens, param, ba::is_any_of(","));
        }
        else
        {
            ba::split(tokens, param, ba::is_any_of(":-"));
        }
        if (tokens.size() != 3)
            throw Exception("failed to parse range: "+param);
        return Range(crr.getChromosomeId(tokens[0]),
                     util::parseValue<uint32_t>(tokens[1]),
                     util::parseValue<uint32_t>(tokens[2]));
    }

    int DecodeCrr::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "range");

        CrrFile crr(referenceFileName_);
        std::ostream& out = openStdout(outputFileName_);

        Range range = parseRangeParam(crr, range_);
        out << crr.getSequence(range) << endl;

        return 0;
    }

} } // cgatools::command
