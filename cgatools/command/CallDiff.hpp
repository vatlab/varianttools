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

#ifndef CGA_TOOLS_COMMAND_CALLDIFF_HPP_
#define CGA_TOOLS_COMMAND_CALLDIFF_HPP_ 1

//! @file CallDiff.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/variants/Locus.hpp"
#include "cgatools/variants/CallDiffResult.hpp"

namespace cgatools { namespace command {

    class SomaticScoreCalc;

    class CallDiff : public Command
    {
    public:
        CallDiff(const std::string& name);

    protected:
        void addOptions(po::options_description& options,
                        po::positional_options_description& positionalOptions);

        int run(po::variables_map& vm);

    private:
        // locusStats[0 or 1][locusClass][locusMatchClass] = count
        typedef boost::array< std::map<std::string, std::map<
            std::vector<variants::cdmt::MatchType>, int> >, 2 > LocusStats;

        // superlocusStats[classification][matchClass] = count
        typedef std::map< std::vector<variants::cdmt::MatchType>,int> SuperlocusStats;

        // locusId -> locus-match-class
        typedef std::map< uint32_t, std::vector<variants::cdmt::MatchType> > LocusDiffClassStore;

        std::string getClassification(
            const variants::Locus& locus, const reference::CrrFile& crr) const;
        void dumpLocusStats(std::ostream& out,
                            LocusStats& locusStats,
                            bool pct) const;
        void dumpSuperlocusStats(std::ostream& out,
                                 SuperlocusStats& superlocusStats) const;
        void writeVariantOutput(const reference::CrrFile& crr,
                                const std::string& inpFn, const std::string& outFn,
                                const LocusDiffClassStore& ldc,
                                SomaticScoreCalc* somatic) const;

        std::string referenceFileName_;
        std::string variantFileNameA_;
        std::string variantFileNameB_;
        std::string oPrefix_;
        std::string reports_;
        size_t extend3Mers_;
        size_t extendBases_;
        size_t statsColumnCount_;
        size_t maxHypothesisCount_;
        bool noReferenceCoverValidation_;
        bool eaf_;
        std::string exportRootA_;
        std::string exportRootB_;
        std::string calibPrefix_;
        std::string somaticScoreDebugInput_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_CALLDIFF_HPP_
