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
#include "cgatools/mapping/SamOptions.hpp"
#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace cgatools { namespace mapping {


void addSamConfigOptions(SamGeneratorConfig& config, boost::program_options::options_description &options)
{
    options.add_options()
        ("skip-not-mapped", po::bool_switch(&config.skipNotMapped_),
        "Skip not mapped records")

        ("add-mate-sequence", po::bool_switch(
        &config.addMateSequenceAndScore_)->default_value(false),
        "Generate mate sequence and score tags.")

        ("mate-sv-candidates", po::bool_switch(
        &config.mateSvCandidates_)->default_value(false),
        "Inconsistent mappings are normally converted as single arm mappings with no mate "
        "information provided. If the option is used map2sam will mate unique single arm "
        "mappings in SAM including those on different stands and chromosomes. "
        "To distinguish these \"artificially\" mated records a tag \"XS:i:1\" is used. "
        "The MAPQ provided for these records is a single arm mapping weight.")

        ("add-unmapped-mate-info",po::bool_switch(
        &config.addMateSequenceAndScoreIfNotLocallyAvailable_)->default_value(false),
        "works like add-mate-sequence, but is applied to inconsistent mappings only")

        ("primary-mappings-only",po::bool_switch(&config.primaryMappingsOnly_)->default_value(false),
        "report only the best mappings")

        ("consistent-mapping-range",po::value<int>(&config.maxConsistentRange_)->default_value(1300),
        "limit the maximum distance between consistent mates")

        ;
}

} } // cgatools::mapping
