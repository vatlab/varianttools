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
#include "cgatools/command/Map2Sam.hpp"
#include "cgatools/mapping/Map2SamConverter.hpp"

#include <boost/filesystem/path.hpp>

namespace cgatools { namespace command {

    Map2Sam::Map2Sam(const std::string& name)
        : Command(name,
                  "Converts CGI initial reference mappings into SAM format.",
                  "0.3 or later",
        "The Map2Sam converter takes as input Reads and Mappings files from "
        "a Complete Genomics data package, the library files (found automatically "
        "inside the package) and a crr reference file and generates one SAM file as an output. "
        "The output is sent into stdout by default. All the mapping records from the input are "
        "converted into corresponding SAM records one to one. In addition, the unmapped DNB records "
        "are reported as SAM records having appropriate indication. Map2Sam converter tries "
        "to identify primary mappings and highlight them using the appropriate flag. "
        "The negative gaps in CGI mappings are represented using GS/GQ/GC tags."
        )
        ,config_(new mapping::Map2SamConfig())
    {
        options_.add_options()
            ("reads,r", po::value<std::string>(&config_->inputReadsFileName_),
             "Input reads file.")

            ("mappings,m", po::value<std::string>(&config_->inputMappingsFileName_),
             "Input mappings file.")

            ("genome-root", po::value<std::string>(&config_->exportRootDirName_),
             "The genome directory, for example /data/GS00118-DNA_A01; this directory is "
             "expected to contain a LIB directory.")

            ("reference,s", po::value<std::string>(&config_->referenceFileName_),
             "Reference file.")

            ("output,o", po::value<std::string>(&config_->outputFileName_)->default_value("STDOUT"),
             "The output SAM file (may be omitted for stdout).")

            ("from,f", po::value<size_t>(&config_->recordsFrom_)->default_value(0),
             "Defines start read record.")

            ("to,t", po::value<size_t>(&config_->recordsTo_)->
             default_value(std::numeric_limits<size_t>::max()),
             "Defines end read record (the end record is not included in the results).")

            ("extract-genomic-region,e", po::value<util::StringVector>(&config_->exportRegions_),
             "Defines a region as a half-open interval 'chr,from,to'")

            ("cgi-bam-caveat-acknowledgement", po::bool_switch(&config_->caveat_)->default_value(false),
             "Supply this option to enable map2sam. Before enabling this option, please "
             "contact support@completegenomics.com for how "
             "to run map2sam and training on the caveats to interpreting "
             "the resulting BAM files.")
            ;

        addSamConfigOptions(config_->samGeneratorConfig_, options_);

        hiddenOptions_.add_options()
            ("debug-info,d", po::bool_switch(&config_->dumpDebugInfo_)->default_value(false),
            "Dump debug information together with output. ")
            ;

    }

    int Map2Sam::run(po::variables_map& vm)
    {
        requireParam(vm, "reads");
        requireParam(vm, "mappings");
        requireParam(vm, "genome-root");
        if (!config_->caveat_)
            throw util::Exception("Please contact support@completegenomics.com for how "
                                  "to run map2sam and training on the caveats to interpreting "
                                  "the resulting BAM files.");

        config_->commandLine_ = getCommandLine();
        mapping::BaseMap2SamConverter converter(*config_, openStdout(config_->outputFileName_));
        converter.init();
        converter.run();

        return 0;
    }

} } // cgatools::command
