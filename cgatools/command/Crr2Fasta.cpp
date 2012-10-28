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
#include "cgatools/command/Crr2Fasta.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/Exception.hpp"

#include <iostream>
#include <vector>

#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace command {

    using std::cerr;
    using std::cin;
    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;

    using util::Exception;

    namespace ba = boost::algorithm;

    Crr2Fasta::Crr2Fasta(const std::string& name)
        : Command(name,
                  "Converts a crr reference file to the fasta format.",
                  "",
                  "")
    {
        options_.add_options()
            ("input", po::value<string>(&inputFileName_),
             "The input crr file (may be passed in as argument at the end of the command).")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output fasta file (may be omitted for stdout).")
            ("line-width", po::value<size_t>(&maxLineWidth_)->default_value(50),
             "The maximum width of a line of sequence.")
            ;

        positionalOptions_.add("input", 1);
    }

    int Crr2Fasta::run(po::variables_map& vm)
    {
        requireParam(vm, "input");

        reference::CrrFile crr(inputFileName_);
        crr.validate();

        std::ostream& out = openStdout(outputFileName_);
        const vector<reference::CompactDnaSequence>& chromosomes = crr.listChromosomes();
        BOOST_FOREACH(const reference::CompactDnaSequence& chromosome, chromosomes)
        {
            out << ">" << chromosome.getName() << endl;
            string sequence;
            for(size_t pos=0; pos<chromosome.length(); pos+=maxLineWidth_)
            {
                size_t posEnd = std::min(pos+50, chromosome.length());
                sequence.clear();
                chromosome.appendSequence(sequence, pos, posEnd-pos);
                out << sequence << "\n";
            }
        }

        return 0;
    }

} } // cgatools::command
