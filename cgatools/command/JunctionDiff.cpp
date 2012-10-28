// Copyright 2010-2012 Complete Genomics, Inc.
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
#include "cgatools/command/JunctionDiff.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"
#include "cgatools/junctions/JunctionStat.hpp"
#include "cgatools/junctions/JunctionVcfWriter.hpp"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/regex.hpp>
#include <boost/shared_ptr.hpp>

using std::ostream;
using std::string;
using std::vector;
using cgatools::junctions::JunctionRef;
using cgatools::junctions::JunctionRefSide;
using cgatools::junctions::JunctionRefSet;
using cgatools::junctions::JunctionVcfWriter;
using cgatools::junctions::JunctionCompatMapPerFile;

namespace cgatools { namespace command {

    class DiffJunctions
    {
    public:
        typedef boost::shared_ptr<std::istream> InStream;

        DiffJunctions(const JunctionDiffConfig& config)
            :config_(config), fileFieldSeparator_('\t')
        {
            init();
        }

        void run() 
        {
            #ifdef ENABLE_JUNCTIONDIFF_STAT
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
            #endif // ENABLE_JUNCTIONDIFF_STAT

            junctions::FullGenomeJunctionComparator 
                cmp(config_.distanceTolerance_,config_.scoreThresholds_,*reference_,config_.outPrefix_);
            cmp.compare(*junctFiles_);

            junctions::JunctionRefs incompatible;
            BOOST_FOREACH(const junctions::JunctionRef& ref, cmp.getIncompatible()[0])
            {
                const junctions::JunctionSideSection& s0 = 
                    (*ref).sideSections_[junctions::JUNCTION_LEFT_SIDE];
                const junctions::JunctionSideSection& s1 = 
                    (*ref).sideSections_[junctions::JUNCTION_RIGHT_SIDE];

                if ( s0.strand_ != junctions::JUNCTION_PLUS_STRAND
                     || s1.strand_ != junctions::JUNCTION_PLUS_STRAND
                     || s0.position_.chromosome_ != s1.position_.chromosome_
                     || (*ref).getDistance()>=config_.minJunctionLength_)
                {
                    incompatible.push_back(ref);
                }
            }

            // print incompatible junctions
            cmp.printList(config_.outPrefix_+"diff-"+
                boost::filesystem::path(config_.junctionFiles_[0]).leaf(), 
                incompatible, (*junctFiles_)[0].metadata_, (*junctFiles_)[0].annotationHeaders_);

            if (config_.outputVcf_)
            {
                writeAllJunctionsToVcf(cmp, config_.outPrefix_+"junctiondiff.vcf");
            }

            if (config_.calcFileStat_)
            { // print report
                util::OutputStream reportStr_(config_.outPrefix_+"report.tsv");

                util::DelimitedFileMetadata reportMetadata;
                reportMetadata.initDefaults();
                reportMetadata.set("TYPE","JUNCTIONDIFF_REPORT");
                reportStr_ << reportMetadata;

                // write short summary report
                reportStr_ << ">fileId"
                    << '\t' << "inputJunctions"
                    << '\t' << "filteredJunctions"
                    << '\t' << "incompatible"
                    << '\t' << "filteredIncompatible"
                    << '\t' << "compatible" 
                    << '\t' << "scoreThreshold" 
                    << '\t' << "maxDistance" 
                    << '\t' << "minDelLength"
                    << std::endl;
                for (size_t i=0; i<cmp.getPrefiltered().size(); ++i)
                    reportStr_ << 'A'+i
                    << '\t' << (*junctFiles_)[i].junctions_.size()
                    << '\t' << cmp.getPrefiltered()[i].size()
                    << '\t' << cmp.getIncompatible()[i].size()
                    << '\t' << (i==0 ? incompatible.size() : 0)
                    << '\t' << cmp.getPrefiltered()[i].size()-cmp.getIncompatible()[i].size()
                    << '\t' << config_.scoreThresholds_[i]
                    << '\t' << config_.distanceTolerance_
                    << '\t' << config_.minJunctionLength_
                    << std::endl;
            }
                
        }

    protected:

        const JunctionDiffConfig&   config_;

        boost::shared_ptr<junctions::JunctionFiles> junctFiles_;

        boost::shared_ptr<reference::CrrFile>  reference_;

        char fileFieldSeparator_;

        void init() 
        {
            reference_.reset(new reference::CrrFile());
            reference_->open(config_.referenceFileName_);

            junctFiles_.reset( new junctions::JunctionFiles() );
            junctFiles_->resize(config_.junctionFiles_.size());

            for (size_t i=0; i<config_.junctionFiles_.size(); ++i) 
            {
                (*junctFiles_)[i].read(config_.junctionFiles_[i],*reference_);
            }

        }

    private:

        //! Writes one side of a VCF "adjacency" to the given stream.
        //! Note that, regardless of our junction side strand, VCF record
        //! is always written relative to the primary strand.
        //! Forks each junction into two single-side structures so that we
        //! can sort the sides in the reference order. Also populates the
        //! compatibility maps for each file, which are necessary because
        //! the comparison algorithm populates the matchedJunctionId_ field
        //! only in one direction, and we want symmetric references in the
        //! VCF file.
        void copyJunctionListForVcf(
                const junctions::JunctionRefs& jrl,
                vector<JunctionRefSide>& out,
                JunctionCompatMapPerFile& compat) const
        {
            BOOST_FOREACH(const JunctionRef& jr, jrl)
            {
                BOOST_FOREACH(const JunctionRef& cjr, jr.matchedJunctionIts_)
                {
                    if (jr.sourceId_ == cjr.sourceId_)
                        continue; // shouldn't really happen
                    compat[jr.sourceId_][jr.it_->id_].
                        push_back(cjr);
                    compat[cjr.sourceId_][cjr.it_->id_].
                        push_back(jr);
                }

                for (size_t side = 0; side < 2; ++side)
                    out.push_back(JunctionRefSide(jr, side));
            }
        }

        //////////////////////////////////////
        void writeAllJunctionsToVcf(
            const junctions::FullGenomeJunctionComparator& cmp,
            const std::string& fn) const
        {
            JunctionVcfWriter vcfWriter( reference_, junctFiles_ );
            vcfWriter.fileFieldSeparator_ = "\t";
            vcfWriter.filterSideLength_   = config_.vcfFilterSideLength_;
            vcfWriter.filterScoreThreshold_  = config_.vcfFilterScoreThreshold_;

            boost::shared_ptr<std::ostream> out =
                util::OutputStream::openCompressedOutputStreamByExtension(fn);
            
            vcfWriter.writeJunctionVcfHeaders(*out);
            
            vector<JunctionRefSide> all;
            JunctionCompatMapPerFile compat(junctFiles_->size());
            for (uint32_t srcId = 0; srcId < junctFiles_->size(); ++srcId)
            {
                copyJunctionListForVcf(cmp.getCompatible()[srcId], all, compat);
                copyJunctionListForVcf(cmp.getIncompatible()[srcId], all, compat);
            }

            std::sort(all.begin(), all.end());

            // Out of each cluster of matched junctions, we output only
            // the first one. The rest are collected in the following set.
            JunctionRefSet junctionsToSuppress;
            BOOST_FOREACH(const JunctionRefSide& jsr, all)
            {
                const JunctionRef& jr = jsr.jr_;
                if (junctionsToSuppress.count(jr) != 0)
                    continue;

                if ( config_.vcfNormalPriorityOutput_ )
                {
                    pickNormalPriorityMatch( jr, junctionsToSuppress );
                }
                else
                {
                    pickDefaultMatch( jr, junctionsToSuppress );
                }
            }

            BOOST_FOREACH(const JunctionRefSide& jsr, all)
            {
                if (0 == junctionsToSuppress.count(jsr.jr_))
                    vcfWriter.writeJunctionToVcf(jsr.jr_, jsr.side_, compat, *out);
            }
        }

        void pickNormalPriorityMatch ( const JunctionRef& jr, JunctionRefSet& junctionsToSuppress ) const
        {
            const JunctionRef* keep = &jr;

            BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
            {
                if (junctionsToSuppress.count(match) != 0)
                    continue;

                if ( keep->sourceId_ < match.sourceId_ ) // normal has higher source ID
                {
                    keep = &match;
                }
                else if ( keep->sourceId_ == match.sourceId_ ) 
                { // if source ID is the same then take junction with lower position
                    if ( match < *keep )
                    {
                        keep = &match;
                    }
                }
            }

            if ( !(*keep == jr) )
                junctionsToSuppress.insert(jr);

            BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
            {
                if ( *keep == match )
                    continue;
                junctionsToSuppress.insert(match);
            }
        }

        void pickDefaultMatch ( const JunctionRef& jr, JunctionRefSet& junctionsToSuppress ) const
        {
            BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
            {
                if (jr == match)
                    continue;
                junctionsToSuppress.insert(match);
            }
        }
        //////////////////////////////////////
    };

    JunctionDiff::JunctionDiff(const std::string& name)
        : Command(name,
                  "Reports difference between junction calls of Complete Genomics junctions files.",
                  "1.5 or later",
                  "junctiondiff takes two junction files A and B as input "
                  "and produces the following output:"
                  "\n  - \t\"diff-inputFileName\" - the junctions from an input file A "
                  "that are not present in input file B."
                  "\n  - \t\"report.txt\" - a brief summary report (if --statout is used)\n"
                  "\nTwo junctions are considered equivalent if:"
                  "\n  - \tthey come from different files"
                  "\n  - \tleft and right positions of one junction are not more than "
                  "\"--distance\" bases apart from the corresponding positions of another junction"
                  "\n  - \tthe junction scores are equal or above the scoreThreshold"
                  "\n  - \tthey are on the same strands"
        )
    {
        options_.add_options()
            ("reference,s", po::value<std::string>(&config_.referenceFileName_),
             "Reference file.")
            ("junctionsA,a", po::value<std::string>(&config_.junctionFiles_[0]),"input junction file A.")
            ("junctionsB,b", po::value<std::string>(&config_.junctionFiles_[1]),"input junction file B.")
            ("scoreThresholdA,A", po::value<size_t>(&config_.scoreThresholds_[0])->default_value(10),
             "score threshold value for the input file A.")
            ("scoreThresholdB,B", po::value<size_t>(&config_.scoreThresholds_[1])->default_value(0),
             "score threshold value for the input file B.")
            ("distance,d", po::value<uint32_t>(&config_.distanceTolerance_)->default_value(200),
             "Max distance between coordinates of potentially compatible junctions.")
            ("minlength,l", po::value<uint32_t>(&config_.minJunctionLength_)->default_value(500),
             "Minimum deletion junction length to be included into the difference file.")
            ("output-prefix,o", po::value<std::string>(&config_.outPrefix_)->default_value(""),
             "The path prefix for all the output reports.")
            ("statout,S", po::bool_switch(&config_.calcFileStat_),
             "(Debug) Report various input file statistics. Experimental feature.")
            ;

        hiddenOptions_.add_options()
            ("output-vcf", po::bool_switch(&config_.outputVcf_)->default_value(false),
                    "Flag to write VCF output.")
            ("vcf-score-threshold",
                    po::value<size_t>(&config_.vcfFilterScoreThreshold_)->default_value(10),
                    "Cutoff value for the VCF quality filter flag, e.g. 10")
            ("vcf-side-len-threshold",
                    po::value<size_t>(&config_.vcfFilterSideLength_)->default_value(70),
                    "Cutoff value for the VCF side length filter flag, e.g. 70")
            ("vcf-normal", po::bool_switch(&config_.vcfNormalPriorityOutput_)->default_value(false),
                    "Normal junction priority for VCF output.")
             ;
    }

    int JunctionDiff::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");

        if (config_.junctionFiles_.size()!=2)
            CGA_ERROR_EX("Wrong number of input junction files: "
                    << config_.junctionFiles_.size()
                    << " The number of input junction files is currently limited to 2."
            );

        if (config_.scoreThresholds_.size()>config_.junctionFiles_.size())
            CGA_ERROR_EX("The number of input score thresholds is bigger "
                << "than the number of input junction files: "
                << CGA_VOUT(config_.scoreThresholds_.size())<<CGA_VOUT(config_.junctionFiles_.size()));

        //set score threshold to 0 by default
        config_.scoreThresholds_.resize(config_.junctionFiles_.size(),0);

        DiffJunctions diff(config_);
        diff.run();

        return 0;
    }


} } // cgatools::command
