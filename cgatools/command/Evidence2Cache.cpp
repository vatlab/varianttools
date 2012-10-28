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
#include "cgatools/command/Evidence2Cache.hpp"
#include "cgatools/mapping/EvidenceSamUtil.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/CrrFile.hpp"

#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/mapping/EvidenceCache.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time.hpp>
#include <boost/filesystem.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

namespace cgatools { namespace command {

    Evidence2Cache::Evidence2Cache(const std::string& name)
        : Command(name,
                  "Converts CGI variant evidence data into internal cache format. "
                  "The output of this command is used to merge evidence mapping data with the base mappings",
                  "0.3 or later",
        "The evidence2cache converter takes as input evidence mapping files as a part as CGI export package "
        "and generates a directory containing cache data. "
        )
    {
        options_.add_options()
            ("reference,s", po::value<std::string>(&config_.referenceFileName_),
             "The reference crr file.")
            ("genome-root,r", po::value<std::string>(&config_.exportRootDirName_),
             "The genome directory, for example /data/GS00118-DNA_A01; "
             "this directory is expected to contain ASM subdirectory.")
            ("output,o", po::value<std::string>(&config_.outputDirName_),"Output directory.")
            ("extract-genomic-region", po::value<util::StringVector>(&config_.exportRegionList_),
             "Defines a region as a half-open interval 'chr,from,to' or 'chr'. "
             "Multiple intervals can be defined byusing the option multiple times.")
            ("debug-output,v", po::bool_switch(&config_.verboseOutput_)->default_value(false),
             "Generate verbose debug output. Please don't rely on this option in production.")
            ;
    }

    int Evidence2Cache::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "genome-root");
        requireParam(vm, "output");

        requireParam(vm, "output");
        config_.init();

        cgdata::GenomeMetadata  genomeMetadata(config_.exportRootDirName_);

        mapping::EvidenceCacheBuilder   cacheBuilder(genomeMetadata,config_.reference_);
        cacheBuilder.exportRanges(config_.outputDirName_, *config_.exportRegions_);

        return 0;
    }

    void Evidence2CacheConfig::init()
    {
        CGA_ASSERT(!referenceFileName_.empty());
        reference_.open(referenceFileName_);
        exportRegions_.reset(new util::FastRangeSet(reference_));

        if (exportRegionList_.empty())
            exportRegions_->addWholeReference();
        else
            exportRegions_->add(exportRegionList_);
    }

} } // cgatools::command
