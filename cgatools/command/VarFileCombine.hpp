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

#ifndef CGA_TOOLS_COMMAND_VARFILECOMBINE_HPP_
#define CGA_TOOLS_COMMAND_VARFILECOMBINE_HPP_ 1

//! @file VarFileCombine.hpp

#include "cgatools/core.hpp"

#include <boost/shared_ptr.hpp>

#include "cgatools/command/Command.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/variants/Locus.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/util/DelimitedFile.hpp"

namespace cgatools { namespace command {

    class AnnotationSource;

    class VarFileCombine : public Command
    {
    public:
        VarFileCombine(const std::string& name);

    protected:
        int run(po::variables_map& vm);

    private:
        std::string referenceFileName_, variantFileName_, exportRoot_;
        std::string outputFileName_, annotations_;
        std::string repmaskFileName_, segdupFileName_;

        reference::CrrFile crr_;
        boost::shared_ptr<cgdata::GenomeMetadata> exp_;
        std::vector< boost::shared_ptr<AnnotationSource> > annotators_;

        void addAnnotators(const variants::VariantFileIterator& srcFile);
        void fillMetadata(util::DelimitedFileMetadata& meta,
                          const util::DelimitedFileMetadata& vfMeta);
        void printLocus(std::ostream& out, const variants::Locus& loc);
    };

} } // cgatools::command

#endif
