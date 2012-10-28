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
#include "cgatools/command/ListCrr.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/Exception.hpp"

#include <boost/format.hpp>
#include <boost/foreach.hpp>

namespace cgatools { namespace command {

    using std::cerr;
    using std::cin;
    using std::cout;
    using std::endl;
    using std::string;
    using std::vector;

    using util::Exception;
    using reference::Location;
    using reference::Range;
    using namespace cgatools::reference;

    ListCrr::ListCrr(const std::string& name)
        : Command(name,
                  "Lists chromosomes, contigs, or ambiguous sequences of a crr file.",
                  "",

                  "For mode=chromosome, prints a space-separated table describing each "
                  "chromosome within the reference. The columns are defined as follows:\n\n"
                  "    ChromosomeId \tA numeric identifier for the chromosome.\n"
                  "    Chromosome   \tThe name of the chromosome.\n"
                  "    Length       \tThe length in bases of the chromosome.\n"
                  "    Circular     \tBoolean indicating if the chromosome is circular.\n"
                  "    Md5          \tMd5 of the string containing the upper case IUPAC "
                  "code for each base in the chromosome (spaces and dashes are omitted).\n"
                  "\n"
                  "For mode=contig, prints a space-separated table describing each "
                  "gap and each contig within the reference. Here, a gap between "
                  "contigs is defined as any stretch of min-contig-gap-length or "
                  "more no-called reference bases (N character). The columns are defined "
                  "as follows:\n\n"
                  "    ChromosomeId \tA numeric identifier for the chromosome.\n"
                  "    Chromosome   \tThe name of the chromosome.\n"
                  "    Type         \tEither CONTIG or GAP.\n"
                  "    Offset       \tThe 0-based offset of the start of the contig or "
                  "gap within the chromosome.\n"
                  "    Length       \tThe length in bases of the contig or gap.\n"
                  "\n"
                  "For mode=ambiguity, prints a space-separated table describing each "
                  "run of ambiguity codes within the reference. The columns are defined "
                  "as follows:\n\n"
                  "    ChromosomeId \tA numeric identifier for the chromosome.\n"
                  "    Chromosome   \tThe name of the chromosome.\n"
                  "    Code         \tThe IUPAC code for the region.\n"
                  "    Offset       \tThe 0-based offset of the run of ambiguity codes "
                  "in the chromosome.\n"
                  "    Length       \tThe length in bases of the run of ambiguity codes."
            )
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file (may be passed in as argument at the end of the command).")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ("mode", po::value<string>(&mode_)->default_value("chromosome"),
             "One of chromosome, contig, or ambiguity.")
            ("min-contig-gap-length", po::value<uint32_t>(&minContigGapLength_)->default_value(50),
             "Minimum length of gap between reference contigs, for mode=contig.")
            ;

        positionalOptions_.add("reference", 1);
    }

    int ListCrr::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");

        CrrFile crr(referenceFileName_);
        std::ostream& out = openStdout(outputFileName_);

        if ("chromosome" == mode_)
        {
            out << "ChromosomeId Chromosome    Length Circular Md5" << endl;
            const vector<CompactDnaSequence>& chromosomes = crr.listChromosomes();
            for(size_t ii=0; ii<chromosomes.size(); ii++)
            {
                const CompactDnaSequence& chromosome = chromosomes[ii];
                out << boost::format("%12d %10s %9d %8s %s") %
                    ii % chromosome.getName() % chromosome.length() %
                    ( chromosome.isCircular() ? "true" : "false" ) %
                    chromosome.getMd5Digest().hex() << endl;
            }
        }
        else if ("contig" == mode_)
        {
            out << "ChromosomeId Chromosome   Type    Offset    Length" << endl;
            const vector<CompactDnaSequence>& chromosomes = crr.listChromosomes();
            vector<Range> contigs = crr.listContigs(minContigGapLength_);
            Location loc(0,0);
            BOOST_FOREACH(const Range& contig, contigs)
            {
                while (loc < contig.beginLocation())
                {
                    if (loc.chromosome_ < contig.chromosome_)
                    {
                        if (loc.offset_ != chromosomes[loc.chromosome_].length())
                            out << boost::format("%12d %10s    GAP %9d %9d") %
                                loc.chromosome_ % chromosomes[loc.chromosome_].getName() %
                                loc.offset_ % (chromosomes[loc.chromosome_].length()-loc.offset_)
                                 << endl;
                        loc = Location(loc.chromosome_+1, 0);
                    }
                    else
                    {
                        out << boost::format("%12d %10s    GAP %9d %9d") %
                            loc.chromosome_ % chromosomes[loc.chromosome_].getName() %
                            loc.offset_ % (contig.begin_-loc.offset_)
                             << endl;
                        loc.offset_ = contig.begin_;
                    }
                }

                out << boost::format("%12d %10s CONTIG %9d %9d") %
                    contig.chromosome_ % chromosomes[contig.chromosome_].getName() %
                    contig.begin_ % contig.length()
                     << endl;
                loc.offset_ = contig.end_;
            }
        }
        else if ("ambiguity" == mode_)
        {
            out << "ChromosomeId Chromosome Code    Offset    Length" << endl;
            const vector<CompactDnaSequence>& chromosomes = crr.listChromosomes();
            for(size_t ii=0; ii<chromosomes.size(); ii++)
            {
                BOOST_FOREACH(const AmbiguousRegion& amb, chromosomes[ii].getAmbiguousRegions())
                {
                    out << boost::format("%12d %10s    %c %9d %9d") %
                        ii % chromosomes[ii].getName() % amb.code_ % amb.offset_ % amb.length_
                         << endl;
                }
            }
        }
        else
            throw Exception("unrecognized mode: "+mode_);

        return 0;
    }

} } // cgatools::command
