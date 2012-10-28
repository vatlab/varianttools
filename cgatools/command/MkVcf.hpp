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

#ifndef CGA_TOOLS_COMMAND_MKVCF_HPP_
#define CGA_TOOLS_COMMAND_MKVCF_HPP_ 1

//! @file MkVcf.hpp

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"

#include <queue>

namespace cgatools { namespace command {

    class MkVcf : public Command
    {
    public:
        MkVcf(const std::string& name);

    protected:
        int run(po::variables_map& vm);
        
    private:
        void initSources();
        void validateSources();
        void writeHeaders(std::ostream& out);
        std::string getAssemblyId(size_t idxGenome) const;
        void transferSource(
            std::ostream& out) const;
        void transferMetadataHeaders(
            std::ostream& out) const;
        std::set<std::string> getValues(
            const std::vector< std::vector< std::vector<cgatools::conv::VcfKvHeaderRecord> > >& recs,
            const std::string& key,
            bool mustExist,
            bool mustBeUnique) const;
        void transferInfoFilterHeaders(
            std::ostream& out) const;
        void transferInfoFilterHeadersForKey(
            std::ostream& out,
            const std::vector< std::vector<cgatools::conv::VcfSubFieldHeaderRecord> >& recs,
            cgatools::conv::VcfSubFieldHeaderRecord::Key key) const;
        std::string getPrintKey(cgatools::conv::VcfSubFieldHeaderRecord::Key key) const;
        void checkConsistency(
            const std::vector< std::vector<cgatools::conv::VcfSubFieldHeaderRecord> >& recs,
            size_t ii0,
            size_t jj0) const;

        std::string referenceFileName_;
        std::string outputFileName_;
        std::vector<std::string> genomeRoots_;
        std::vector<std::string> varNames_;
        std::string fieldNames_;
        std::string sourceNames_;
        bool includeNoCalls_;
        std::string calibPrefix_;
        size_t genomeCount_;

        std::vector<std::string> junctionFileNames_;
        size_t junctionScoreThreshold_;
        size_t junctionSideLengthThreshold_;
        size_t junctionDistanceTolerance_;
        size_t junctionLengthThreshold_;
        bool   junctionNormalPriority_;
        bool   junctionTumorHC_;

        cgatools::reference::CrrFile crr_;

        typedef boost::shared_ptr<cgatools::conv::VcfRecordSource> PQItem;
        typedef std::vector<PQItem> PQContainer;
        class PQItemComparator
        {
        public:
            bool operator()(const PQItem&, const PQItem&) const;
        };
        std::priority_queue< PQItem, PQContainer, PQItemComparator > sources_;
        std::vector<PQItem> vSources_;
        std::vector<std::string> vSourceNames_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_MKVCF_HPP_
