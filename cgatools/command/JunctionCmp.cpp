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
#include "cgatools/command/JunctionCmp.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"
#include "cgatools/junctions/JunctionStat.hpp"

#include <boost/filesystem.hpp>

namespace cgatools { namespace command {


    class CompareJunctions 
    {
    public:
        typedef boost::shared_ptr<std::istream> InStream;

        CompareJunctions(const JunctionCmpConfig& config)
            :config_(config), fileFieldSeparator_('\t')
        {
            init();
        }

        void run() 
        {
            if (config_.calcFileStat_)
                for (size_t i=0; i<junctFiles_.size();++i) {
                    util::OutputStream statStr(config_.outPrefix_+"stat-"+
                        boost::filesystem::path(config_.junctionFiles_[i]).leaf());
                    statStr <<"Statistics for " << config_.junctionFiles_[i] << ":" 
                        << std::endl 
                        << std::endl;
                    junctions::JuctionStatisticsCollector(config_.distanceTolerance_).
                                                            run(junctFiles_[i].junctions_, statStr);
                }

            junctions::FullGenomeJunctionComparator 
                cmp(config_.distanceTolerance_,config_.scoreThresholds_,reference_,config_.outPrefix_);
            cmp.compare(junctFiles_);

            // print compatible junctions
            cmp.printMultiList(config_.outPrefix_+"compatible_pairs.tsv", 
                cmp.getCompatibleCombined(), junctFiles_[0].metadata_, 2);

            // print compatible junctions
            for (size_t i=0; i<cmp.getCompatible().size(); ++i)
                cmp.printList(config_.outPrefix_+"compatible_"
                    +boost::filesystem::path(config_.junctionFiles_[i]).leaf(),
                cmp.getCompatible()[i], junctFiles_[i].metadata_, junctFiles_[i].annotationHeaders_);

            // print incompatible junctions
            for (size_t i=0; i<cmp.getIncompatible().size(); ++i)
                cmp.printList(config_.outPrefix_+"incompatible_"
                    +boost::filesystem::path(config_.junctionFiles_[i]).leaf(),
                cmp.getIncompatible()[i], junctFiles_[i].metadata_, junctFiles_[i].annotationHeaders_);

        }

    protected:

        const JunctionCmpConfig&   config_;

        junctions::JunctionFiles    junctFiles_;

        reference::CrrFile          reference_;

        char fileFieldSeparator_;

        void init() 
        {
            reference_.open(config_.referenceFileName_);
            junctFiles_.resize(config_.junctionFiles_.size());
            for (size_t i=0; i<config_.junctionFiles_.size(); ++i) 
            {
                junctFiles_[i].read(config_.junctionFiles_[i],reference_);
            }
        }
    };

    JunctionCmp::JunctionCmp(const std::string& name)
        : Command(name,
                  "Compares junction calls of Complete Genomics junctions files ",
                  "1.5 or later",
                  "junctioncmp takes several junction files as input "
                  "(limited to 2 in the current implementation) and produces the following output:\n"
                  "  - compatible-j.tsv - the records of the file are combinations of compatible junctions "
                  "(one junction from each file in a record). "
                  " A junction from one file is included into several records if there are several "
                  "compatible junctions in the other files.\n"
                  "  - incompatible[0..N]-j.tsv - the junctions from an input file "
                                    "that don't have a match in any of the other input files."
                                    "The junctions from incompatible*-j\n"
                  "  - report.txt - a brief summary report\n"
                  "\nThe input junction files must be in the canonical form: the left side start position "
                  "is less then right side start position.\n"
                  "\nTwo junctions are compatible if:\n"
                  "  - they come from different files\n"
                  "  - left and right positions of one junction are not more then "
                  "\"--distance\" bases apart from the corresponding positions of another junction\n"
                  "  - the junction scores are equal or above the <scoreThreshold>\n"
                  "  - they are on the same strands"
        )
    {
        options_.add_options()
            ("reference,s", po::value<std::string>(&config_.referenceFileName_),
             "Reference file.")
            ("junctions,j", po::value<std::vector<std::string> >(&config_.junctionFiles_),
             "input junction files. The order is important.")
            ("scoreThreshold,t", po::value<std::vector<size_t> >(&config_.scoreThresholds_),
             "score threshold value for each input file. The order is important. "
             "(optional. The default threshold is 0)")
            ("distance,d", po::value<uint32_t>(&config_.distanceTolerance_)->default_value(200),
             "Distance between coordinates of potentially compatible junctions.")
            ("output-prefix,o", po::value<std::string>(&config_.outPrefix_)->default_value(""),
             "The path prefix for all the output reports.")
            ("statout,S", po::bool_switch(&config_.calcFileStat_),
             "Report various input file statistics. Experimental feature.")
            ;
    }

    int JunctionCmp::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");

        if (config_.junctionFiles_.size()!=2)
            CGA_ERROR_EX("Wrong number of input junction files: "
                    << config_.junctionFiles_.size()
                    << " The number of input junction files is currently limited to 2."
            );

        if (config_.scoreThresholds_.size() > config_.junctionFiles_.size())
            CGA_ERROR_EX("The number of input score thresholds is bigger "
                << "than the number of input junction files: "
                << config_.scoreThresholds_.size());

        //set score threshold to 0 by default
        config_.scoreThresholds_.resize(config_.junctionFiles_.size(),0);

        CompareJunctions compare(config_);
        compare.run();

        return 0;
    }


} } // cgatools::command
