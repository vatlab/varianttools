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

#ifndef CGATOOLS_COPYNUMBER_CNVFILEVCFSOURCE_HPP_
#define CGATOOLS_COPYNUMBER_CNVFILEVCFSOURCE_HPP_ 1

//! @file CnvFileVcfSource.hpp

#include "cgatools/core.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/DelimitedFile.hpp"

#include <deque>
#include <map>
#include <set>
#include <vector>

// Forward declarations
namespace cgatools { namespace cgdata {
    class GenomeMetadata;
} } // cgatools::cgdata

namespace cgatools { namespace copynumber {


    // structures for temporary data
    struct diploidData {
        std::string chr;
        int begin;
        int end;
        cgatools::reference::Range range;
        float avgNormCvg;
        float gcCvg;
        float normedGcCvg;
        float relCvg;
        int cP;
        std::string cT;
        std::string pS;
        std::string tS;
    };

    struct nondiploidData {
        std::string chr;
        int begin;
        int end;
        cgatools::reference::Range range;
        float cL;
        std::string lS;
    };

    struct somaticData {
        std::string chr;
        int begin;
        int end;
        cgatools::reference::Range range;
        float cL;
        std::string lS;
        std::string laf;
        std::string llaf;
        std::string ulaf;
    };


    class CnvFileVcfRecordWriter : public cgatools::conv::VcfRecordWriter
    {
    public:
        CnvFileVcfRecordWriter(
            const std::vector< std::string> cnvFieldNames,
            const cgatools::reference::CrrFile& crr,
            int numGenomes,
            bool someSomatic);

        // Overridden methods of VcfRecordWriter.
        cgatools::reference::Location getLocation() const;
        void writeRef(std::ostream& out) const;
        void writeAlt(std::ostream& out) const;
        void writeInfo(std::ostream& out) const;
        void writeFormat(std::ostream& out) const;
        void writeSample(std::ostream& out, size_t gIdx) const;

        // methods for updating current writeable record
        void setDiploidData(diploidData &data, int gIdx);
        void setNondiploidData(nondiploidData &data, int gIdx);
        void setSomaticData(somaticData &data, int gIdx);

        // get ranges for each type
        cgatools::reference::Range getDiploidRange(int gIdx) const { return dipData_[gIdx].range; }
        cgatools::reference::Range getNondiploidRange(int gIdx) const {return nondipData_[gIdx].range; }
        cgatools::reference::Range getSomaticRange(int gIdx) const {return somData_[gIdx].range; }

        // determine whether we have a real record to write for a given genomes
        bool isDeferred(int gIdx) const { return deferred_[gIdx]; }
        // set whether we have a real record to write for a given genomes
        void setDeferred(int gIdx,bool status) { deferred_[gIdx] = status; }

    private:
        const cgatools::reference::CrrFile& crr_;
        bool someSomatic_;
        std::vector< diploidData > dipData_;
        std::vector< nondiploidData > nondipData_;
        std::vector< somaticData > somData_;
        std::vector< std::string > formatIds_;
        std::vector< bool > deferred_;
    };

    class CnvFileVcfRecordSource : public cgatools::conv::VcfRecordSource
    {
    public:
        CnvFileVcfRecordSource(
            const std::vector< std::string >& dipdet,
            const std::vector< std::string >& nondipdet,
            const std::vector< std::string >& somnondipdet,
            const std::vector<std::string> fieldNames,
            const cgatools::reference::CrrFile& crr);

        // Overridden methods of VcfRecordSource.

        std::vector<cgatools::conv::VcfSubFieldHeaderRecord> getSubFieldHeaderRecords() const;
        std::string getSource(size_t idxGenome) const;
        std::vector<cgatools::conv::VcfKvHeaderRecord> getKeyValueHeaderRecords(size_t idxGenome) const;
        std::string getAssemblyId(size_t idxGenome) const;

        bool eof() const;
        cgatools::conv::VcfRecordSource& operator++();
        const cgatools::conv::VcfRecordWriter& operator*() const;
        const cgatools::conv::VcfRecordWriter* operator->() const;

    private:

        void limitFieldNames(const std::vector< std::string > & fieldNames);
        void setUpDelimitedFiles();
        void setupDiploidDetails(boost::shared_ptr<cgatools::util::DelimitedFile> &df);
        void setupNondiploidDetails(boost::shared_ptr<cgatools::util::DelimitedFile> &df);
        void setupSomaticDetails(boost::shared_ptr<cgatools::util::DelimitedFile> &df);
        bool testForSomaticData();
        void computeGcCorrectedMeans();

    private:

        const cgatools::reference::CrrFile& crr_;
        cgatools::reference::Range genomeEnd_;
        std::vector< std::string > dipDetFn_;
        std::vector< std::string > nondipDetFn_;
        std::vector< std::string > somnondipDetFn_;
        std::vector< boost::shared_ptr<std::istream> > dipDetIStr_;
        std::vector< boost::shared_ptr<std::istream> > nondipDetIStr_;
        std::vector< boost::shared_ptr<std::istream> > somnondipDetIStr_;
        std::vector< boost::shared_ptr<cgatools::util::DelimitedFile> > dipDet_;
        std::vector< boost::shared_ptr<cgatools::util::DelimitedFile> > nondipDet_;
        std::vector< boost::shared_ptr<cgatools::util::DelimitedFile> > somnondipDet_;
        std::vector< float > gcCorrectedMean_;
        diploidData tmpDip_;
        nondiploidData tmpNondip_;
        somaticData tmpSom_;
        std::vector< std::string > cnvFieldNames_;
        std::set< std::string > cnvFieldNamesSet_;
        boost::shared_ptr<CnvFileVcfRecordWriter> writer_;
        bool someSomatic_;
        bool eof_;
    };

} } // cgatools::copynumber

#endif // CGATOOLS_COPYNUMBER_CNVFILEVCFSOURCE_HPP_
