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

#ifndef CGATOOLS_CGDATA_GENOMEMETADATA_HPP_
#define CGATOOLS_CGDATA_GENOMEMETADATA_HPP_ 1

//! @file GenomeMetadata.hpp

#include "cgatools/core.hpp"
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/array.hpp>

namespace cgatools { namespace cgdata {

    class LibraryMetadata
    {
    public:
        typedef boost::array<std::string,2> SequenceDependentGapFileNames;

        LibraryMetadata(const boost::filesystem::path& libraryPath,const std::string libraryName);
        //! left side - 0, right side - 1
        std::string             getGapDistributionFileName(uint8_t side) const;
        SequenceDependentGapFileNames   getSequenceGapDistributionFileName(uint8_t side) const;
        std::string             getMainGapDistributionFileName() const;

        std::string             getDnbStructureFileName() const;
        const std::string &     getName() {return libraryName_;}
    protected:
        boost::filesystem::path libraryPath_;
        std::string             libraryName_;
    };

    class GenomeMetadata : boost::noncopyable
    {
    public:
        GenomeMetadata(const std::string& root);

        std::string getEvidenceIntervalsFileName(const std::string& chromName) const;
        std::string getEvidenceDnbsFileName(const std::string& chromName) const;
        std::string getReferenceSupportFileName(const std::string& chromName) const;
        std::string getGeneFileName() const;
        std::string getNcRNAFileName() const;
        std::string getCnvSegmentsDiploidFileName(bool exnOnFail = true) const;
        std::string getCnvSegmentsNondiploidFileName(bool exnOnFail = true) const;
        std::string getCnvDetailsDiploidFileName(bool exnOnFail = true) const;
        std::string getCnvDetailsNondiploidFileName(bool exnOnFail = true) const;
        std::string getCnvDetailsSomaticDiploidFileName(bool exnOnFail = true) const;
        std::string getCnvDetailsSomaticNondiploidFileName(bool exnOnFail = true) const;
        std::string getMobileElementInsertionFileName() const;
        std::string getMasterVarFileName(bool exnOnFail = true) const;
        std::string getGenomeReference() const;
        std::string getJunctionsFileName( bool highConfidence = false ) const;
        //! the lane name is in the format "SLIDE-LANE": GS14901-FS3-L03
        std::string getLaneDir(const std::string& laneName) const;
        std::string getLibraryName(const std::string& laneName) const;
        std::vector<std::string> getLaneNames() const;

        //
        //@param laneBatchId - an id in a form "A_B" where "A" is a lane id and "B" is a 3 digit batch id 
        //       example: "GS14901-FS3-L03_001"
        std::string getReadsFileName(const std::string& laneBatchId) const;
        std::string getMappingsFileName(const std::string& laneBatchId) const;

        // returns a metadata class for the given library
        LibraryMetadata getLibraryMetadata(const std::string& libraryName) const;
        std::vector<std::string> getLibraryNames() const;

        int         getFormatVersion() const;
        const std::string& getAsmId() const;

    private:
        boost::filesystem::path root_, asmDir_, asmRefDir_, asmEvidenceDir_, mapDir_, libDir_, svDir_;
        mutable std::string asmIdCache_;
    };

} } // cgatools::cgdata

#endif // CGATOOLS_CGDATA_GENOMEMETADATA_HPP_
