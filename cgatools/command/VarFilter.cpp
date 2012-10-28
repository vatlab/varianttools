// Copyright 2011 Complete Genomics, Inc.
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
#include "cgatools/command/VarFilter.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"
#include "cgatools/util/DelimitedFile.hpp"

using namespace std;
using namespace cgatools::util;
using namespace cgatools::reference;
using namespace cgatools::variants;

namespace cgatools { namespace command {

    VarFilter::VarFilter(const std::string& name)
        : Command(name,
                  "Copies input var file or masterVar file to output, "
                  "applying specified filters.",
                  "0.3 or later",

                  "Copies input var file or masterVar file to output, applying specified "
                  "filters (which are available to all cgatools commands that read a var "
                  "file or masterVar file as input). "
                  "Filters are specified by appending the filter "
                  "specification to the var file name on the command line. For example:\n\n"
                  "/path/to/var.tsv.bz2#varQuality!=VQHIGH\n\n"
                  "The preceding example filters out any calls marked as VQLOW. The filter "
                  "specification follows the \"#\" sign, and consists of a list of filters "
                  "to apply, separated by a comma. Each filter is a colon-separated list of "
                  "call selectors. Any scored call that passes all the colon-separated call "
                  "selectors for one or more of the comma-separated filters is turned into a "
                  "no-call. The following call selectors are available:\n\n"
                  "    hom            \tSelects only calls in homozygous loci.\n"
                  "    het            \tSelects any scored call not selected by the hom selector.\n"
                  "    varType=XX     \tSelects calls whose varType is XX.\n"
                  "    varScoreVAF&lt;XX \tSelects calls whose varScoreVAF&lt;XX.\n"
                  "    varScoreEAF&lt;XX \tSelects calls whose varScoreEAF&lt;XX.\n"
                  "    varQuality!=XX \tSelects calls whose varQuality is not XX.\n"
                  "\nHere is an example that filters homozygous SNPs with varScoreVAF "
                  "&lt; 25 and heterozygous insertions with varScoreEAF &lt; 50:\n\n"
                  "&apos;/path/to/var.tsv.bz2#hom:varType=snp:varScoreVAF&lt;25,"
                  "het:varType=ins:varScoreEAF&lt;50&apos;\n"
            )
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file.")
            ("input", po::value<string>(&inputFileName_),
             "The input var file or masterVar file (typically with filters specified).")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ;

        positionalOptions_.add("input", -1);
    }

    int VarFilter::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "input");

        CrrFile crr(referenceFileName_);
        std::ostream& out = openStdout(outputFileName_);

        VariantFileIterator vf(crr);
        vf.setReferenceCoverValidation(false);
        vf.open(inputFileName_);

        DelimitedFileMetadata metaOut;
        metaOut.initDefaults(vf.getMetadata());
        if ("" != vf.getFilterString())
            metaOut.add("VAR_FILTER", vf.getFilterString());
        out << metaOut;

        if (vf.isOlpl())
        {
            vf.writeOlplFileHeader(out);
            out << endl;
            for(; !vf.eof(); ++vf)
            {
                vf->writeAsOneLine(out);
                out << "\n";
            }
        }
        else
        {
            out << ">" << Call::getHeader() << endl;
            for(; !vf.eof(); ++vf)
                out << (*vf);
        }

        return 0;
    }

} } // cgatools::command
