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

#ifndef CGATOOLS_MOBILEELEMENT_MEIFILEVCFSOURCE_HPP_
#define CGATOOLS_MOBILEELEMENT_MEIFILEVCFSOURCE_HPP_ 1

//! @file MeiFileVcfSource.hpp

#include "cgatools/core.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/DelimitedFile.hpp"

#include <set>
#include <vector>

// Forward declarations
namespace cgatools { namespace cgdata {
    class GenomeMetadata;
} } // cgatools::cgdata

namespace cgatools { namespace mobileelement {

    // structures for temporary data
    struct meiData {
        std::string chr;
        int begin;
        int end;
        int pos; // middle
        uint16_t chrId;
        cgatools::reference::Location location;
        std::string strand;
        std::string ET;
        std::string IS;
        std::string IDC;
        std::string IDCL;
        std::string IDCR;
        int RDC;
        std::string NBET;
        std::string ETS;
        std::string KES;
        int elemBegin;
        int elemEnd;
    };

    class MeiFileVcfRecordWriter : public cgatools::conv::VcfRecordWriter
    {
    public:
        MeiFileVcfRecordWriter(
            const std::vector< std::string> meiFieldNames,
            const cgatools::reference::CrrFile& crr,
            int numGenomes);

        // Overridden methods of VcfRecordWriter.
        cgatools::reference::Location getLocation() const;
        void writeRef(std::ostream& out) const;
        void writeAlt(std::ostream& out) const;
        void writeInfo(std::ostream& out) const;
        void writeFormat(std::ostream& out) const;
        void writeSample(std::ostream& out, size_t gIdx) const;

        // methods for updating current writeable record
        void setMeiData(meiData &data, int gIdx);

    private:
        std::string superType(meiData) const;

        const cgatools::reference::CrrFile& crr_;
        std::vector< meiData > meiData_;
        std::vector< std::string > formatIds_;
    };

    class MeiFileVcfRecordSource : public cgatools::conv::VcfRecordSource
    {
    public:
        MeiFileVcfRecordSource(
            const std::vector< std::string >& meiInput,
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

    private:

        const cgatools::reference::CrrFile& crr_;
        std::vector< std::string > meiFn_;
        std::vector< boost::shared_ptr<std::istream> > meiIStr_;
        std::vector< boost::shared_ptr<cgatools::util::DelimitedFile> > meiInput_;
        meiData tmpMei_;
        std::vector< std::string > meiFieldNames_;
        std::set< std::string > meiFieldNamesSet_;
        boost::shared_ptr<MeiFileVcfRecordWriter> writer_;
        bool eof_;
    };

} } // cgatools::mobileelement

#endif // CGATOOLS_MOBILEELEMENT_MEIFILEVCFSOURCE_HPP_
