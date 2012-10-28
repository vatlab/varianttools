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

#ifndef CGA_TOOLS_SAM_OPTIONS_HPP_
#define CGA_TOOLS_SAM_OPTIONS_HPP_ 1

//! @file SamOptions.hpp

#include "cgatools/core.hpp"

namespace boost { namespace program_options {
    class options_description;
}}

namespace cgatools { namespace mapping {

    //! Sam Generator configuration parameters
    class SamGeneratorConfig
    {
    public:
        SamGeneratorConfig()
            :   
            mateSvCandidates_(false)
            ,addMateSequenceAndScore_(false)
            ,addMateSequenceAndScoreIfNotLocallyAvailable_(false)
            ,printAlleleInfo_(true)
            ,packCigar_(true)
            ,removePaddingAtCigarEnds_(true)
            ,addAlternativeMates_(false)
            ,addAlternativeMappings_(false)
            ,primaryMappingsOnly_(true)
            ,maxConsistentRange_(1300)
        {}

        bool    mateSvCandidates_;
        bool    addMateSequenceAndScore_;
        bool    addMateSequenceAndScoreIfNotLocallyAvailable_;
        bool    printAlleleInfo_;
        //! if packed the cigar looses information about 0-sized gaps 
        //! but does not cause format errors in GATK
        bool    packCigar_;
        bool    removePaddingAtCigarEnds_;
        bool    addAlternativeMates_;
        bool    addAlternativeMappings_;
        bool    skipNotMapped_;
        //generate only one mapping pair per read
        bool    primaryMappingsOnly_;
        int     maxConsistentRange_;
    };

    void addSamConfigOptions(SamGeneratorConfig& config, 
        boost::program_options::options_description &options);

} } // cgatools::mapping

#endif // CGA_TOOLS_SAM_OPTIONS_HPP_
