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
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/Files.hpp"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

namespace cgatools { namespace cgdata {

    using namespace cgatools::util;

    namespace fs = boost::filesystem;
    using std::string;
    using std::vector;

    LibraryMetadata::LibraryMetadata( 
        const boost::filesystem::path& libraryPath,
        const std::string libraryName )
        :libraryPath_(libraryPath)
        ,libraryName_(libraryName)
    {
    }

    std::string LibraryMetadata::getGapDistributionFileName( uint8_t side ) const
    {
        std::string fn = std::string("lib_gaps_rollup_")+((side==0)?"L":"R")+"_"+libraryName_+".tsv";
        return util::files::findDataFile(libraryPath_, fn);
    }

    std::string LibraryMetadata::getMainGapDistributionFileName() const
    {
        std::string fn = std::string("lib_gaps_M_")+libraryName_+".tsv";
        return util::files::findDataFile(libraryPath_, fn);
    }

    std::string LibraryMetadata::getDnbStructureFileName() const
    {
        std::string fn = "lib_DNB_"+libraryName_+".tsv";
        return util::files::findDataFile(libraryPath_, fn);
    }

    LibraryMetadata::SequenceDependentGapFileNames 
        LibraryMetadata::getSequenceGapDistributionFileName( uint8_t side ) const
    {
        SequenceDependentGapFileNames result;
        std::string prefix = std::string("lib_gaps_")+((side==0)?"L":"R");
        std::string suffix = std::string("_")+libraryName_+".tsv";
        result[0] = util::files::findDataFile(libraryPath_, prefix+'0'+suffix);
        result[1] = util::files::findDataFile(libraryPath_, prefix+'1'+suffix);
        return result;
    }

    GenomeMetadata::GenomeMetadata(const std::string& root)
        : root_(root)
    {
        asmDir_ = root_ / "ASM";
        asmRefDir_ = asmDir_ / "REF";
        asmEvidenceDir_ = asmDir_ / "EVIDENCE";
        mapDir_ = root_ / "MAP";
        libDir_ = root_ / "LIB";
        svDir_  = asmDir_ / "SV";

        util::files::check_dir(root_);
    }

    std::string GenomeMetadata::getEvidenceIntervalsFileName(
        const std::string& chromName) const
    {
        std::string fn = "evidenceIntervals-" + chromName + "-" + getAsmId() + ".tsv";
        return util::files::findDataFile(asmEvidenceDir_, fn);
    }

    std::string GenomeMetadata::getEvidenceDnbsFileName(
        const std::string& chromName) const
    {
        std::string fn = "evidenceDnbs-" + chromName + "-" + getAsmId() + ".tsv";
        return util::files::findDataFile(asmEvidenceDir_, fn);
    }

    std::string GenomeMetadata::getJunctionsFileName( bool highConfidence ) const
    {
        std::string fn = std::string(highConfidence ? "highConfidence" : "all") 
                       + "JunctionsBeta-" + getAsmId() + ".tsv";
        return util::files::findDataFile(svDir_, fn);
    }

    std::string GenomeMetadata::getReferenceSupportFileName(
        const std::string& chromName) const
    {
        std::string fn = "coverageRefScore-" + chromName + "-" + getAsmId() + ".tsv";
        return util::files::findDataFile(asmRefDir_, fn);
    }

    std::string GenomeMetadata::getGeneFileName() const
    {
        std::string fn = "gene-" + getAsmId() + ".tsv";
        return util::files::findDataFile(asmDir_, fn);
    }

    std::string GenomeMetadata::getNcRNAFileName() const
    {
        std::string fn = "ncRNA-" + getAsmId() + ".tsv";
        return util::files::findDataFile(asmDir_, fn);
    }

    std::string GenomeMetadata::getCnvSegmentsDiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("cnvSegmentsDiploidBeta-" + getAsmId() + ".tsv");
        fileNames.push_back("cnvSegmentsBeta-" + getAsmId() + ".tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getCnvSegmentsNondiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("cnvSegmentsNondiploidBeta-" + getAsmId() + ".tsv");
        fileNames.push_back("cnvTumorSegmentsBeta-" + getAsmId() + ".tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getCnvDetailsDiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("cnvDetailsDiploidBeta-" + getAsmId() + ".tsv");
        fileNames.push_back("cnvDetailsBeta-" + getAsmId() + ".tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getCnvDetailsSomaticDiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("somaticCnvDetailsDiploidBeta-" + getAsmId() + "(-[NT][1-9])*.tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getCnvDetailsNondiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("cnvDetailsNondiploidBeta-" + getAsmId() + ".tsv");
        fileNames.push_back("cnvTumorDetailsBeta-" + getAsmId() + ".tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getCnvDetailsSomaticNondiploidFileName(bool exnOnFail) const
    {
        fs::path cnvDir = asmDir_ / "CNV";
        vector<string> fileNames;
        fileNames.push_back("somaticCnvDetailsNondiploidBeta-" + getAsmId() + "(-[NT][1-9])*.tsv");
        return util::files::findDataFileRegexMulti(cnvDir, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getMobileElementInsertionFileName() const
    {
        fs::path meiDir = asmDir_ / "MEI";

        return util::files::findDataFile(meiDir, "mobileElementInsertionsBeta-" + getAsmId() + ".tsv");
    }

    std::string GenomeMetadata::getMasterVarFileName(bool exnOnFail) const
    {
        vector<string> fileNames;
        fileNames.push_back("masterVarBeta-" + getAsmId() + ".tsv");
        fileNames.push_back("masterVarBeta-" + getAsmId() + "(-[NT][1-9])*.tsv");
        return util::files::findDataFileRegexMulti(asmDir_, fileNames, exnOnFail);
    }

    std::string GenomeMetadata::getGenomeReference() const
    {
        std::string fn = getGeneFileName();
        boost::shared_ptr<std::istream> in(InputStream::openCompressedInputStreamByExtension(fn));
        DelimitedFile df(*in, fn);
        return df.getMetadata().get("GENOME_REFERENCE");
    }

    int GenomeMetadata::getFormatVersion() const
    {
        std::string fn = getGeneFileName();
        boost::shared_ptr<std::istream> in(InputStream::openCompressedInputStreamByExtension(fn));
        DelimitedFile df(*in, fn);
        return df.getMetadata().getFormatVersion();
    }

    const std::string& GenomeMetadata::getAsmId() const
    {
        if (!asmIdCache_.empty())
            return asmIdCache_;

        boost::regex re("^var-(.+)\\.tsv(\\.bz2|\\.gz)?$");
        fs::directory_iterator itEnd;
        vector<string> varFiles;
        vector<string> asmIds;
        for (fs::directory_iterator it(asmDir_); it != itEnd; ++it)
        {
            boost::smatch what;
            std::string fname(it->leaf());
            if (!boost::regex_match(fname, what, re)) {
                continue;
            }
            asmIdCache_ = what[1];
            varFiles.push_back(it->path().external_file_string());
            asmIds.push_back(asmIdCache_);
        }
        if (asmIdCache_.empty())
        {
            throw util::Exception("failed to determine ASSEMBLY_ID: var file not found in: " +
                                  asmDir_.string());
        }
        for(size_t ii=1; ii<varFiles.size(); ii++)
        {
            if (asmIds[ii] != asmIds[0])
            {
                string msg = "failed to determine ASSEMBLY_ID for genome "+root_.string()+
                    " -- multiple var files:";
                for(size_t ii=0; ii<varFiles.size(); ii++)
                    msg += "\n" + varFiles[ii];
                throw util::Exception(msg);
            }
        }
        return asmIdCache_;
    }

    std::string GenomeMetadata::getLaneDir( const std::string& laneName ) const
    {
        fs::path laneDir = mapDir_ / laneName;
        util::files::check_dir(laneDir);
        return laneDir.string();
    }

    std::string GenomeMetadata::getLibraryName( const std::string& laneName ) const
    {
        fs::path laneDir = mapDir_ / laneName;
        util::files::check_dir(laneDir);
        std::string fileNamePattern = "reads_"+laneName+".*";
        std::vector<std::string> readsFiles = util::files::listDir(laneDir, fileNamePattern, 1);
        if (readsFiles.empty())
            throw util::Exception("Cannot determine the library, no reads files found at: "
                +(laneDir/fileNamePattern).string());
        std::string fname = (laneDir/readsFiles.front()).string();
        boost::shared_ptr<std::istream> in(
            InputStream::openCompressedInputStreamByExtension(
             fname));
        DelimitedFile df(*in, fname);
        return df.getMetadata().get("LIBRARY");
    }

    LibraryMetadata GenomeMetadata::getLibraryMetadata(const std::string& libraryName) const
    {
        fs::path libraryDir = libDir_ / libraryName;
        util::files::check_dir(libraryDir);
        return LibraryMetadata(libraryDir,libraryName);
    }

    std::vector<std::string> GenomeMetadata::getLibraryNames() const
    {
        return util::files::listDir(libDir_,"[GS|SI].*");
    }

    std::vector<std::string> GenomeMetadata::getLaneNames() const
    {
        return util::files::listDir(mapDir_,".*");
    }

    std::string GenomeMetadata::getReadsFileName( const std::string& laneBatchId ) const
    {
        std::string::size_type strEnd = laneBatchId.rfind('_');

        if (strEnd==std::string::npos || strEnd==0)
            throw util::Exception("invalid lane batch id: " + laneBatchId + 
            ". Correct format: GS14901-FS3-L02_003");

        return util::files::findDataFile(mapDir_/laneBatchId.substr(0,strEnd),
                                "reads_" + laneBatchId + ".tsv");
    }

    std::string GenomeMetadata::getMappingsFileName( const std::string& laneBatchId ) const
    {
        std::string::size_type strEnd = laneBatchId.rfind('_');

        if (strEnd==std::string::npos || strEnd==0)
            throw util::Exception("invalid lane batch id: " + laneBatchId + 
            ". Correct format: GS14901-FS3-L02_003");

        return util::files::findDataFile(mapDir_/laneBatchId.substr(0,strEnd),
                                "mapping_" + laneBatchId + ".tsv");
    }

} } // cgatools::cgdata
