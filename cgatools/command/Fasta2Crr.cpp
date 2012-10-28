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
#include "cgatools/command/Fasta2Crr.hpp"
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
    using util::InputStream;

    namespace ba = boost::algorithm;

    Fasta2Crr::Fasta2Crr(const std::string& name)
        : Command(name,
                  "Converts fasta reference files to the crr format.",
                  "",
                  "")
    {
        options_.add_options()
            ("input", po::value< vector<string> >(&inputFileNames_),
             "The input fasta files (may be passed in as arguments at the end of the command, "
             "or omitted for stdin). "
             "Take care to specify the fasta files in chromosome order; ordering is important. "
             "To work with human Complete Genomics data, the chromosome order should be "
             "chr1...chr22, chrX, chrY, chrM.")
            ("output", po::value<string>(&outputFileName_),
             "The output crr file.")
            ("circular", po::value<string>(&circularChromosomes_),
             "A comma-separated list of circular chromosome names. If ommitted, defaults to chrM.")
            ;

        positionalOptions_.add("input", -1);
    }

    int Fasta2Crr::run(po::variables_map& vm)
    {

        if (0 == vm.count("input"))
            inputFileNames_.push_back("STDIN");

        requireParam(vm, "output");
        util::OutputStream out(outputFileName_);

        vector<string> circulars;
        circulars.push_back("chrM");

        if (vm.count("circular") != 0)
        {
            circulars.clear();
            ba::split(circulars, circularChromosomes_, ba::is_any_of(","));
        }

        reference::CrrFileWriter writer(&out);
        BOOST_FOREACH(const std::string& inputFileName, inputFileNames_)
        {
            std::istream& in = openStdin(inputFileName);

            string line;
            for(int lineNumber=1; InputStream::getline(in, line); lineNumber++)
            {
                try
                {
                    if (1 == lineNumber || '>' == line[0])
                    {
                        string name = parseFastaHeader(line);
                        bool circular =
                            std::find(circulars.begin(), circulars.end(), name) != circulars.end();
                        writer.newChromosome(name, circular);
                        continue;
                    }

                    writer.addSequence(line);
                }
                catch(std::exception& ee)
                {
                    throw Exception((boost::format("%s:%d: %s") %
                                     inputFileName % lineNumber % ee.what()).str());
                }
            }
        }

        return 0;
    }

    std::string Fasta2Crr::parseFastaHeader(const std::string& line) const
    {
        if (line.length() < 2 || line[0] != '>') {
            throw Exception("expected FASTA header, found: " + line);
        }

        if (line.find('|') != string::npos)
            throw Exception("NCBI-style fasta headers are not supported: " + line);

        string name = line.substr(1);
        boost::trim(name);
        return name;
    }

} } // cgatools::command
