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
#include "cgatools/command/Sam2Reads.hpp"
#include "cgatools/mapping/SamReader.hpp"
#include "cgatools/mapping/LaneBatchCache.hpp"
#include "cgatools/mapping/Cigar.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/BaseUtil.hpp"

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

namespace cgatools { namespace command {

class SamOutLaneBatchStreams : public mapping::OutLaneBatchStreams 
{
public:
    SamOutLaneBatchStreams (
        const boost::filesystem::path& outputPrefix, 
        const std::string & fileType,
        const std::string & headerLine,
        const mapping::SamFileHeader& samHeader
        )
        : mapping::OutLaneBatchStreams(outputPrefix, fileType,headerLine), samHeader_(samHeader)
    {
    }

protected:
    virtual void addMetadata(util::DelimitedFileMetadata &meta) const
    {
        std::string rgId = meta.get("SLIDE")+'-'+meta.get("LANE");
        const mapping::SamFileHeader::Record& rg = samHeader_.getReadGroup(rgId);
        boost::array<std::string,2> tags = {{"SM","LB"}};

        BOOST_FOREACH(const std::string& s, tags)
        {
            if (!rg.tags_.hasTag(s))
                CGA_ERROR_EX("Missing tag " << s << ", read group: " << rg);
        }

        meta.add("LANE",rgId);
        meta.add("LIBRARY",rg.tags_.getTag("LB").value_);
        meta.add("SAMPLE",rg.tags_.getTag("SM").value_);
        if (rg.tags_.hasTag("DS"))
            meta.add("FIELD_SIZE",rg.tags_.getTag("DS").value_);
    }

    const mapping::SamFileHeader&   samHeader_;
};

Sam2Reads::Sam2Reads(const std::string& name)
    : Command(name,
              "Converts SAM files into lane/batch separated read files.",
              "CGI SAM 1.4 or later",
    "The sam2reads converter takes as input SAM files... "
    )
    ,config_(new Sam2ReadsConfig())
{
    options_.add_options()
        ("input-sam,i", po::value<std::string>(&config_->inputSamFileName_)->default_value("STDIN"),
         "Input SAM file.")

        ("output,o", po::value<std::string>(&config_->outputDirName_),"Output directory.")

        ("from,f", po::value<std::string>(&config_->recordsFrom_),
         "Defines start read record id.")

        ("to,t", po::value<std::string>(&config_->recordsTo_),
         "Defines end read record id (the end record is not included in the results).")

        ("extract-genomic-region,e", po::value<util::StringVector>(&config_->exportRegionList_),
         "Defines a region as a half-open interval 'chr,from,to'")
        ;

    hiddenOptions_.add_options()
        ("debug-info,d", po::bool_switch(&config_->dumpDebugInfo_)->default_value(false),
        "Dump debug information together with output. ")
        ;

}

int Sam2Reads::run(po::variables_map& vm)
{
    requireParam(vm, "input-sam");
    requireParam(vm, "output");

    std::istream& in = openStdin(config_->inputSamFileName_);

    std::string slide;
    std::string lane;
    size_t batchNo;
    size_t offset;

    mapping::SamFileParser parser(in, config_->inputSamFileName_);

    SamOutLaneBatchStreams outStreams(config_->outputDirName_,
        "SAM_READS",
        ">readOffset\tside\tbases\tscores",
        parser.header_);

    while (parser.next())
    {
        const mapping::SamFileRecord & r = parser.getRecord();
        r.parseQName(slide, lane, batchNo, offset);
        if (r.isSecondary() || (r.flag_==0 && r.qname_=="empty"))
            continue;

        std::string updatedSeq;
        std::string updatedQual;
        if (r.samTags_.hasTag("GC"))
        {
            size_t readPos = 0;
            const mapping::SamTag& gapCigar = r.samTags_.getTag("GC");
            mapping::Cigar c(gapCigar.value_);
            for (size_t i=0; i < c.size(); ++i)
            {
                switch(c[i].type_)
                {
                case 'S':
                    updatedSeq.insert(updatedSeq.end(),r.seq_.begin()+readPos,
                        r.seq_.begin()+readPos+c[i].length_);
                    updatedQual.insert(updatedQual.end(),r.qual_.begin()+readPos,
                        r.qual_.begin()+readPos+c[i].length_);
                    readPos += c[i].length_;
                    break;
                case 'G':
                    readPos += c[i].length_;
                    if (!r.samTags_.hasTag("GS"))
                        CGA_ERROR_EX("the record with GC tag but is missing GS tag: " << r);
                    if (!r.samTags_.hasTag("GQ"))
                        CGA_ERROR_EX("the record with GC tag but is missing GQ tag: " << r);
                    updatedSeq.insert(updatedSeq.length(),r.samTags_.getTag("GS").value_);
                    updatedQual.insert(updatedQual.length(),r.samTags_.getTag("GQ").value_);
                    break;
                }
            }
        }

        if (r.isReverseComplemented()) 
        {
            if (updatedSeq.empty())
            {
                updatedSeq = r.seq_;
                updatedQual = r.qual_;
            }

            std::reverse(updatedSeq.begin(),updatedSeq.end());
            BOOST_FOREACH(char &ch, updatedSeq)
                ch = util::baseutil::complement(ch);

            std::reverse(updatedQual.begin(),updatedQual.end());
        }

        std::ostream & ostr = outStreams.getBatchStream(slide, lane, batchNo);
        ostr 
            << offset << '\t'
            << r.getSide() << '\t'
            << (updatedSeq.empty() ? r.seq_ : updatedSeq) << '\t'
            << (updatedQual.empty() ? r.qual_ : updatedQual) << '\t'
            << std::endl;

        if (config_->dumpDebugInfo_)
            std::cout << parser.getRecord() << std::endl;
    }

    return 0;
}

} } // cgatools::command
