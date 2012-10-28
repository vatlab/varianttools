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
#include "cgatools/reference/GeneDataStore.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"

using namespace cgatools::util;

namespace cgatools { namespace reference {

void GeneDataStore::load(const reference::CrrFile& crr, const std::string& fn)
{
    boost::shared_ptr<std::istream> in =
            InputStream::openCompressedInputStreamByExtension(fn);
    DelimitedFile df(*in, fn);
    metadata_ = df.getMetadata();

    Range range;
    GeneDescriptionData data;
    std::vector<uint32_t> exonStarts, exonEnds;
    uint32_t genomeCdsStart, genomeCdsEnd;

    df.addField(StringField("mrnaAcc", &data.mrnaAcc_));
    df.addField(StringField("proteinAcc", &data.proteinAcc_));
    df.addField(StringField("geneId", &data.geneId_));
    df.addField(StringField("symbol", &data.symbol_));
    df.addField(StrandField("orientation", &data.strand_));
    df.addField(ChromosomeIdField("chromosome", &range.chromosome_, crr));
    df.addField(ValueField<uint32_t>("genomeStart", &range.begin_));
    df.addField(ValueField<uint32_t>("genomeEnd", &range.end_));
    df.addField(ValueVectorField<uint32_t>("exonStarts", ';', &exonStarts));
    df.addField(ValueVectorField<uint32_t>("exonEnds", ';', &exonEnds));
    df.addField(ValueField<uint32_t>("genomeCdsStart", &genomeCdsStart));
    df.addField(ValueField<uint32_t>("genomeCdsEnd", &genomeCdsEnd));

    while (df.next())
    {
        if (exonStarts.size() != exonEnds.size())
            throw Exception(fn + ": malformed exon coordinate list: " + df.getLine());
        data.exonRanges_.clear();
        for (size_t ii = 0; ii < exonStarts.size(); ++ii)
        {
            data.exonRanges_.push_back(Range(range.chromosome_,
                                             exonStarts[ii],
                                             exonEnds[ii]));
        }
        data.genomeCodingRange_ = Range(range.chromosome_, genomeCdsStart, genomeCdsEnd);
        tree_.put(range, data);
    }
}

}}
