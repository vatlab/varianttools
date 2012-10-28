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
#include "cgatools/mapping/SamReader.hpp"

#include <boost/foreach.hpp>

namespace cgatools { namespace mapping {

void SamTag::parse(const char* first, const char* last)
{
    const char* current = first;
    const char* lastStart = first;
    for (size_t ind=0; current!=last; ++current)
    {
        if (*current == ':')
        {
            ++ind;
            switch(ind)
            {
            case 1:
                name_ = std::string(lastStart,current);
                break;
            case 2:
                type_ = *lastStart;
                break;
            }
            lastStart = current+1;
            if (2 == ind)
                break;
        }
    }
    value_.assign(lastStart,last);
}

void SamHeaderTag::parse(const char* first, const char* last)
{
    for (const char* current = first; ; ++current)
        if (*current == ':')
        {
            if (current==first)
                CGA_ERROR_EX("Error tag format: " << std::string(first,last));
            name_.assign(first,current);
            value_.assign(current+1,last);
            return;
        }
    CGA_ERROR_EX("Error tag format: " << std::string(first,last));
}

void SamFileRecord::bindToParser( util::DelimitedLineParser &parser )
{
    parser.addField(util::StringField("QNAME", &qname_));
    parser.addField(util::ValueField<uint32_t>("FLAG", &flag_));
    parser.addField(util::StringField("RNAME", &rname_));
    parser.addField(util::ValueField<uint32_t>("POS", &pos_));
    parser.addField(util::ValueField<uint32_t>("MAPQ", &mapq_));
    parser.addField(util::StringField("CIGAR", &cigar_));
    parser.addField(util::StringField("RNEXT", &rnext_));
    parser.addField(util::ValueField<uint32_t>("PNEXT", &pnext_));
    parser.addField(util::ValueField<int32_t>("TLEN", &tlen_));
    parser.addField(util::StringField("SEQ", &seq_));
    parser.addField(util::StringField("QUAL", &qual_));

    for (char i='0'; i<='9'; ++i)
        parser.addField(SamTagField<SamTag>("Tag"+i, samTags_, i==0));
}

void SamFileRecord::parseQName( std::string& slide, std::string& lane, size_t& batchNo, size_t& offset) const
{
    size_t sep = 0;
    for (std::string::const_iterator itPrev=qname_.begin(), 
            it=qname_.begin(), itEnd=qname_.end(); it!=itEnd; ++it)
    {
        if (*it=='-') 
        {
            ++sep;
            switch (sep)
            {
            case 2:
                slide.assign(itPrev,it);
                itPrev = it+1;
                break;
            case 3:
                lane.assign(itPrev,it);
                itPrev = it+1;
                break;
            }
        } else if (*it==':') {
            CGA_ASSERT_MSG(sep == 3, "Wrong QNAME format: " << qname_ << " sep:" << sep);
            batchNo = util::parseValue<size_t>(std::string(itPrev,it));
            offset = util::parseValue<size_t>(std::string(++it,itEnd));
            break;
        }
    }
}


bool SamFileParser::next()
{
    while(in_) 
    {
        getline(in_,line_);
        ++lineNo_;

        if (line_.empty())
            continue;

        if (line_[0]=='@')
        {
            header_.addLine(line_);
            continue;
        }

        if (!isInitialised_)
        {
            header_.initMetadata();
        }

        samFileRecord_.samTags_.clear();
        dl_.parseLine(line_);

        return true;
    }

    return false;
}

void SamFileHeader::initMetadata()
{
    if (hasRecords("RG"))
    {
        const SamFileHeader::Records& records = getRecords("RG");
        BOOST_FOREACH(const SamFileHeader::Record& r, records)
        {
            if (!r.tags_.hasTag("ID"))
                CGA_ERROR_EX("The header RG record is missin tag ID: " << r);
            const SamHeaderTag& t = r.tags_.getTag("ID");
            readGroups_[t.value_] = &r;
        }
    }
}

void SamFileParser::initMetadata()
{
    std::cout << samFileRecord_;
    header_.initMetadata();
    isInitialised_ = true;
}

bool SamFileHeader::hasReadGroup( const std::string& rgId ) const
{
    return readGroups_.find(rgId) != readGroups_.end();
}

const SamFileHeader::Record& SamFileHeader::getReadGroup( const std::string& rgId ) const
{
    HeaderRecordIndex::const_iterator it = readGroups_.find(rgId);
    if (it == readGroups_.end())
        CGA_ERROR_EX("Read group not found: " << rgId);
    return *it->second;
}

std::ostream& operator<< (std::ostream& ost, const SamFileRecord& r) 
{
    ost << "QNAME: " << r.qname_ << '\t';
    ost << "F:" 
        << (r.isConsistent()?"dbl":"sgl") << '-'
        << (r.isUnmapped()?"U":"M") << '-'
        << (r.isMateUnmapped()?"mU":"mM") << '-'
        << (r.isReverseComplemented()?"Rev":"Nrm") << '-'
        << (r.isMateReverseComplemented()?"mRev":"mNrm") << '-'
        << '\t';

    ost << "RNAME: " << r.rname_ << '\t';
    ost << "POS: " << r.pos_ << '\t';
    ost << "MAPQ: " << r.mapq_ << '\t';
    ost << "CIGAR: " << r.cigar_ << '\t';
    ost << "RNEXT: " << r.rnext_ << '\t';
    ost << "PNEXT: " << r.pnext_ << '\t';
    ost << "TLEN: " << r.tlen_ << '\t';
    ost << r.seq_ << '\t';
    ost << r.qual_;

    
    BOOST_FOREACH(const SamTags<SamTag>::SamTagMap::value_type & t, r.samTags_.samTags_)
    {
        ost << '\t' << t.second.name_ 
            << ':' << t.second.type_ 
            << ':' << t.second.value_ 
            ;
    }

    return ost;
}

void SamFileHeader::Record::bindToParser(util::DelimitedLineParser &parser)
{
    parser.addField(util::StringField("tagId", &id_));

    for (char i='0'; i <= ('0'+30); ++i)
        parser.addField(SamTagField<SamHeaderTag>("Tag"+i, tags_, i==0));
}

void SamFileHeader::addLine( const std::string& headerLine )
{
    samHeaderRecord_.tags_.clear();
    dl_.parseLine(headerLine);

    if (samHeaderRecord_.id_.size() < 2 || samHeaderRecord_.id_[0]!='@')
        CGA_ERROR_EX("Wrong header record format: " << headerLine);

    samHeaderRecord_.id_ = samHeaderRecord_.id_.substr(1);
    recordMap_[samHeaderRecord_.id_].push_back(samHeaderRecord_);
}

std::ostream& operator<< (std::ostream& ost, const SamTag& r)
{
    ost
        << r.name_ << ':'
        << r.type_ << ':'
        << r.value_;
    return ost;
}

std::ostream& operator<< (std::ostream& ost, const SamHeaderTag& r)
{
    ost
        << r.name_ << ':'
        << r.value_;
    return ost;
}

std::ostream& operator<< (std::ostream& ost, const SamFileHeader::Record& r)
{
    ost
        << r.id_ << '\t';
    r.tags_.print(ost);

    return ost;
}

} // namespace mapping
} // cgatools

