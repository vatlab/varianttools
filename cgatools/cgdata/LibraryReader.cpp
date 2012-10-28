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
#include "cgatools/cgdata/LibraryReader.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Streams.hpp"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

namespace cgatools { namespace cgdata {


template <>
void GapTable<SmallGapTuple>::load( const std::string& fname )
{
    boost::shared_ptr<std::istream> gapFileStream = 
        util::InputStream::openCompressedInputStreamByExtension(fname);
    util::DelimitedFile parser(*gapFileStream, fname);
    GapTable::GapsRecord record;
    util::DelimitedLineParser lineParser;
    lineParser.addField(util::ValueField<int8_t>("gap0", &record.gaps_[0]));
    lineParser.addField(util::ValueField<int8_t>("gap1", &record.gaps_[1]));
    lineParser.addField(util::ValueField<int8_t>("gap2", &record.gaps_[2]));
    parser.addField(util::DelimitedField("gaps", ';', lineParser));
    parser.addField(util::ValueField<double>("frequency", &record.frequency_));

    while (parser.next())
        records_.push_back(record);

    std::sort(records_.begin(),records_.end());
}

template <>
void GapTable<MateGapSize>::load( const std::string& fname )
{
    boost::shared_ptr<std::istream> gapFileStream = 
        util::InputStream::openCompressedInputStreamByExtension(fname);
    util::DelimitedFile parser(*gapFileStream, fname);
    GapTable::GapsRecord record;
    parser.addField(util::ValueField<MateGapSize>("mateGap", &record.gaps_));
    parser.addField(util::ValueField<double>("frequency", &record.frequency_));

    while (parser.next())
        records_.push_back(record);

    std::sort(records_.begin(),records_.end());
}

} } // cgatools::cgdata
