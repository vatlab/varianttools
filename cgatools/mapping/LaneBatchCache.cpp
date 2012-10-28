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
#include "cgatools/util/parse.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Files.hpp"
#include "cgatools/mapping/LaneBatchCache.hpp"

#include <boost/foreach.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/regex.hpp>

using namespace std;

namespace cgatools { namespace mapping {

std::string BaseLaneBatchStreams::getBatchStreamKey( 
    const std::string& slide, const std::string& lane, size_t batchNo)
{
    return slide + "-" + lane + "_" + boost::lexical_cast<string>(batchNo);
}

std::ostream& OutLaneBatchStreams::getBatchStream( 
    const std::string& slide, const std::string& lane, size_t batchNo)
{
    string key = BaseLaneBatchStreams::getBatchStreamKey(slide,lane,batchNo);

    OutputLaneBatches::iterator it = outputLanes_.find(key);
    if (it==outputLanes_.end()) {
        std::string fileNum = boost::lexical_cast<std::string>(batchNo);
        boost::filesystem::create_directories(outputPrefix_);
        boost::filesystem::path outputFile = outputPrefix_ / (slide+"-"+lane+"_"+fileNum+".tsv.gz");
        it = outputLanes_.insert(it,key,
            new util::CompressedOutputStream(outputFile.string()));
        writeMetadata(*it->second, slide, lane, batchNo);
    }
    return *it->second;
}

void OutLaneBatchStreams::writeMetadata( 
    std::ostream &output, 
    const std::string& slide, 
    const std::string& lane, 
    size_t batchNo 
    ) const
{
    std::string fileNum = boost::lexical_cast<std::string>(batchNo);

    util::DelimitedFileMetadata meta;
    meta.initDefaults();
    meta.add("TYPE",fileType_);
    meta.add("SLIDE",slide);
    meta.add("LANE",lane);
    meta.add("BATCH_FILE_NUMBER",fileNum);

    addMetadata(meta);

    output << meta;

    if (!headerLine_.empty())
        output << headerLine_ << std::endl;
}

void InLaneBatchStreams::collectFiles()
{
    boost::regex re("^(.+)-(L[0-9]+)_([0-9]+)\\.tsv(\\.bz2|\\.gz)?$");
    std::vector<boost::filesystem::path> filesFound;
    util::files::findFiles(inputRoot_,"GS.*",filesFound,true,false);
    BOOST_FOREACH(const boost::filesystem::path& fname,filesFound)
    {
        boost::smatch what;
        std::string leaf = fname.leaf();

        if (!boost::regex_match(leaf, what, re)) 
            continue;

        string key = 
            BaseLaneBatchStreams::getBatchStreamKey(what[1],what[2],util::parseValue<size_t>(what[3]));

        InputBatchFiles& files = inputBatches_[key];
        files.push_back(fname);
    }
}

const InLaneBatchStreams::InputBatchFiles& InLaneBatchStreams::getBatchFiles( 
    const std::string& slide, const std::string& lane, size_t batchNo) const
{
    string key = BaseLaneBatchStreams::getBatchStreamKey(slide,lane,batchNo);

    InputLaneBatches::const_iterator it = inputBatches_.find(key);

    if (it == inputBatches_.end())
        return emptyBatch_;
    else
        return *it->second;
}

} } // cgatools::mapping
