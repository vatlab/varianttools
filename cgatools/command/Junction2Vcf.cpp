// Copyright 2012 Complete Genomics, Inc.
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
#include "cgatools/command/Junction2Vcf.hpp"
#include "cgatools/junctions/JunctionVcfWriter.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

using namespace std;
using namespace cgatools::util;
using namespace cgatools::reference;

using boost::array;
using boost::lexical_cast;
using boost::shared_ptr;
using cgatools::junctions::JunctionVcfWriter;
using cgatools::junctions::JunctionRef;
using cgatools::junctions::JunctionRefs;
using cgatools::junctions::JunctionRefSide;
using cgatools::junctions::JunctionCompatMapPerFile;

namespace cgatools { namespace command {

    Junction2Vcf::Junction2Vcf(const std::string& name)
        : Command(name,
                  "Converts junction file to VCF.",
                  "1.5 or later",
                  "Converts junction file to VCF."
            )
    {
        junctionsFileName_ = "";
        outputFileName_    = "STDOUT";
        scoreThreshold_    = 10;
        minJunctionLength_ = 70;

        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file.")
            ("input,i", po::value<string>(&junctionsFileName_),
                "junctions file to be converted")
            ("score-threshold,s", po::value<size_t>(&scoreThreshold_)->default_value(scoreThreshold_),
                "score threshold value for the input junctions")
            ("min-length,l", po::value<size_t>(&minJunctionLength_)->default_value(minJunctionLength_),
                "length threshold value for the input junctions")
            ("output", po::value<string>(&outputFileName_)->default_value(outputFileName_),
             "The output file (may be omitted for stdout).")
            ;

        positionalOptions_.add("variants", -1);
    }

    int Junction2Vcf::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "input");

        reference_.reset(new reference::CrrFile());
        reference_->open(referenceFileName_);

        junctionFiles_.reset( new junctions::JunctionFiles() );
        junctionFiles_->resize(1);
        (*junctionFiles_)[0].read(junctionsFileName_,*reference_);

        std::ostream& out = openStdout(outputFileName_);
        saveVcf(out);

        return 0;
    }

    void Junction2Vcf::saveVcf ( std::ostream& out )
    {
        JunctionVcfWriter vcfWriter( reference_, junctionFiles_ );
        vcfWriter.fileFieldSeparator_ = "\t";
        vcfWriter.filterSideLength_      = minJunctionLength_;
        vcfWriter.filterScoreThreshold_  = scoreThreshold_;

        vcfWriter.writeJunctionVcfHeaders(out);

        JunctionRefs junctionRefs( (*junctionFiles_)[0].junctions_, 0 );
        vector<JunctionRefSide> all;
        all.reserve(junctionRefs.size()*2);

        BOOST_FOREACH(const JunctionRef& jr, junctionRefs)
        {
            for (size_t side = 0; side < 2; ++side)
            {
                all.push_back(JunctionRefSide(jr, side));
            }
        }
        std::sort(all.begin(), all.end());

        JunctionCompatMapPerFile compat(1);

        BOOST_FOREACH(const JunctionRefSide& jsr, all)
        {
            vcfWriter.writeJunctionToVcf(jsr.jr_, jsr.side_, compat, out);
        }

    }

} } // cgatools::command
