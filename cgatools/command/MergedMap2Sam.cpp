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
#include "cgatools/command/MergedMap2Sam.hpp"
#include "cgatools/mapping/MergeMap2SamConverter.hpp"

namespace cgatools { namespace command {

    MergedMap2Sam::MergedMap2Sam(const std::string& name)
        : Command(name,
                  "Converts CGI initial reference mappings and variant evidence mappings into SAM format.",
                  "0.3 or later",
        "The MergedMap2Sam converter takes as input the root directory of an export package, "
        "a lane batch id, "
        "a CRR reference file, the root directory of the evidence cache created by evidence2cache command "
        "and generates one SAM file as an output. The output is sent into stdout by default. "
        "The MAPQ in the mapping records from the both sources are recomputed, "
        "duplicated mappings are removed. "
        "The result subset can be optionally filtered (min-mapq, primary-mappings-only,...). "
        "The result mappings are converted into corresponding SAM records: "
        "one half DNB mapping to one SAM record. "
        "The negative gaps in CGI mappings are represented using GS/GQ/GC tags."
        )
        ,config_(new mapping::MergedMap2SamConfig())
    {
        options_.add_options()
            ("lane-batch-id,l", po::value<std::string>(&config_->inputReadsBatchId_),
            "Input lane batch in the format SLIDE-LANE_BATCH. Example: GS10364-FS3-L01_001")

            ("genome-root,r", po::value<std::string>(&config_->exportRootDirName_),
             "The export package root, for example /data/GS00118-DNA_A01; "
             "this directory is expected to contain MAP, ASM and LIB subdirectories.")

            ("reference,s", po::value<std::string>(&config_->referenceFileName_),
             "Reference file in CRR format.")

            ("evidence-cache-root,c", 
             po::value<std::string>(&config_->evidenceCacheRoot_)->default_value(""),
             "Root directory of the evidence cache created by evidence2cache. "
             "Several cache directories can be combined by providing the common root.")

            ("output,o", po::value<std::string>(&config_->outputFileName_)->default_value("STDOUT"),
             "The output SAM file (may be omitted for stdout).")

            ("outputStream,O",po::value<util::StringVector>(&config_->outputStreamNames_),
             "The output SAM stream name for a chromosome. Format: 'chrName,streamName'. "
             "Example: '--outputStream chr1,./chr1file'. "
             "Use 'all' for all the chromosomes that are not listed (--output value by default).")

            ("extract-genomic-region,e", po::value<util::StringVector>(&config_->exportRegions_),
             "Defines a region as a half-open interval 'chr,from,to'. "
             "The option can be used multiple times to define multiple regions")

            ("add-allele-id", po::bool_switch(&config_->samGeneratorConfig_.printAlleleInfo_),
             "generate interval ID and allele ID tags for SAM records representing evidence mappings")
        ;

        addSamConfigOptions(config_->samGeneratorConfig_, options_);

        options_.add_options()
            ("min-mapq", po::value<uint16_t>(&tmpMapq_)->default_value(0),
             "filter out mappings with MAPQ < min-mapq.")

            ("alpha", po::value<double>(&config_->mapqAlpha_)->default_value(1E-12),
             "affects the mapping score: decreasing the value of alpha "
             "makes it the score more sensitive to bad mappings")

            ("from,f", po::value<size_t>(&config_->recordsFrom_)->default_value(0),
             "Defines start read record.")

            ("to,t", po::value<size_t>(&config_->recordsTo_)->
             default_value(std::numeric_limits<size_t>::max()),
             "Defines end read record (the end record is not converted).")
            ;

        hiddenOptions_.add_options()
            ("debug-info,d", po::bool_switch(&config_->dumpDebugInfo_)->default_value(false),
            "Dump debug information together with output.")
            ;
    }

    int MergedMap2Sam::run(po::variables_map& vm)
    {
        requireParam(vm, "lane-batch-id");
        requireParam(vm, "genome-root");
        requireParam(vm, "reference");

        if (tmpMapq_>255)
            CGA_ERROR_EX("Invalid MAPQ:"<< tmpMapq_ <<". MAPQ should be in the range 0-255");
        config_->minMapQ_ = boost::uint8_t(tmpMapq_);

        {
            cgdata::GenomeMetadata gm(config_->exportRootDirName_);
            config_->inputReadsFileName_ = gm.getReadsFileName(config_->inputReadsBatchId_);
            config_->inputMappingsFileName_ = gm.getMappingsFileName(config_->inputReadsBatchId_);
        }

        config_->commandLine_ = getCommandLine();
        mapping::MergedMap2SamConverter converter(*config_, openStdout(config_->outputFileName_));
        converter.init();
        converter.run();

        return 0;
    }

} } // cgatools::command
