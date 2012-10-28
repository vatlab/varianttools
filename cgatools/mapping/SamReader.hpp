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

#ifndef CGA_TOOLS_SAM_READER_HPP_
#define CGA_TOOLS_SAM_READER_HPP_ 1

//! @file SamReader.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/DelimitedLineParser.hpp"

#include <map>

namespace cgatools { namespace mapping {

class SamTag
{
    friend std::ostream& operator<< (std::ostream& ost, const SamTag& r);
public:
    void parse(const char* first, const char* last);

    std::string name_;
    char        type_;
    std::string value_;
};

class SamHeaderTag
{
    friend std::ostream& operator<< (std::ostream& ost, const SamHeaderTag& r);
public:
    void parse(const char* first, const char* last);

    std::string name_;
    std::string value_;
};

template <class SamTag>
class SamTags
{
public:
    typedef typename std::map<std::string, SamTag> SamTagMap;

    bool hasTag(const std::string& t) const
    {
        return samTags_.find(t) != samTags_.end();
    }

    const SamTag& getTag(const std::string& t) const
    {
        typename SamTagMap::const_iterator it = samTags_.find(t);
        CGA_ASSERT_MSG(it!=samTags_.end(), "The tag '" << t << "' not found.");
        return it->second;
    }

    void addTag(const SamTag& t)
    {
        samTags_[t.name_] = t;
    }

    void print(std::ostream& ost) const
    {
        for (typename SamTagMap::const_iterator it=samTags_.begin(),
                itEnd = samTags_.end(); it!=itEnd; ++it)
        {
            if (it != samTags_.begin())
                ost << '\t';
            ost << it->second;
        }
    }

    void clear() {samTags_.clear();}

    SamTagMap samTags_;
};

template <class SamTag>
class SamTagField : public util::DelimitedFieldParser
{
public:
    // samTags should be emptied before each iteration of parsing of tags
    // it is checked by an assertion in a first tag in a list
    // This is the reason to require that the first tag object in a list is flagged
    SamTagField(const std::string& name, SamTags<SamTag>& samTags, bool isFirst)
        : DelimitedFieldParser(name),
        isFirst_(isFirst),
        samTags_(samTags)
    {
    }

    void parse(const char* first, const char* last) {
        CGA_ASSERT_MSG((!isFirst_) || samTags_.samTags_.empty(), "samTags map was not cleaned");
        SamTag tag;
        tag.parse(first, last);
        samTags_.addTag(tag);
    }

private:
    bool                isFirst_;
    SamTags<SamTag>&    samTags_;
};

class SamFileRecord 
{
    friend std::ostream& operator<< (std::ostream& ost, const SamFileRecord& r);

public:
    SamFileRecord() : flag_(-1), pos_(-1), mapq_(-1), pnext_(-1), tlen_(-1) {}

    void bindToParser(util::DelimitedLineParser &parser);

    //each segment properly aligned according to the aligner
    bool isConsistent() const {return (flag_ & 0x2) > 0;}

    //segment unmapped
    bool isUnmapped() const {return (flag_ & 0x4) > 0;}

    //next segment in the template unmapped
    bool isMateUnmapped() const {return (flag_ & 0x8) > 0;}

    //SEQ being reverse complemented
    bool isReverseComplemented() const {return (flag_ & 0x10) > 0;}

    //SEQ of the next segment in the template being reversed
    bool isMateReverseComplemented() const {return (flag_ & 0x20) > 0;}

    //the first segment in the template
    bool isFirstMate() const {return (flag_ & 0x40) > 0;}

    //the last segment in the template
    bool isLastMate() const {return (flag_ & 0x80) > 0;}

    //0 - first mate, 1 - secons mate
    int  getSide() const {return (flag_&0x40)>0 ? 0:((flag_&0x80)>0 ? 1:-1);}

    //secondary alignment
    bool isSecondary() const {return (flag_ & 0x100) > 0;}

    void parseQName(std::string& slide, std::string& lane, size_t& batchNo, size_t& offset) const;

    std::string qname_;
    uint32_t    flag_;
    std::string rname_;
    uint32_t    pos_;
    uint32_t    mapq_;
    std::string cigar_;
    std::string rnext_;
    uint32_t    pnext_;
    int32_t     tlen_;
    std::string seq_;
    std::string qual_;

    SamTags<SamTag> samTags_;
};

class SamFileHeader
{
public:
    class Record 
    {
        friend std::ostream& operator<< (std::ostream& ost, const Record& r);

    public:
        void bindToParser(util::DelimitedLineParser &parser);

        std::string             id_;
        SamTags<SamHeaderTag>   tags_;
    };

    typedef std::vector<Record> Records;
    typedef std::map<std::string, Records> RecordMap;

    typedef std::map<std::string, const Record *> HeaderRecordIndex;

    SamFileHeader()
    : dl_(true) 
    {
        samHeaderRecord_.bindToParser(dl_);
    }

    bool hasRecords(const std::string& section) const
    {
        return recordMap_.find(section)!=recordMap_.end();
    }

    const Records& getRecords(const std::string& section) const
    {
        RecordMap::const_iterator it = recordMap_.find(section);
        CGA_ASSERT_MSG(it!=recordMap_.end(), "Section is missing: " << section);
        return it->second;
    }

    Records& getRecords(const std::string& section)
    {
        return recordMap_[section];
    }

    void addLine(const std::string& headerLine);

    RecordMap   recordMap_;

    bool hasReadGroup(const std::string& rgId) const;
    const Record& getReadGroup(const std::string& rgId) const;

    void initMetadata();

protected:
    Record      samHeaderRecord_;
    util::DelimitedLineParser   dl_;

    HeaderRecordIndex   readGroups_;
};

class SamFileParser {
public:
    SamFileParser(std::istream& in, const std::string& fileName)
        : isInitialised_(false), in_(in), fileName_(fileName),lineNo_(0), dl_(true) 
    {
        samFileRecord_.bindToParser(dl_);
    }

    bool next();

    const SamFileRecord & getRecord() const {return samFileRecord_;}

    void initMetadata();

    SamFileHeader   header_;

protected:
    bool            isInitialised_;

    std::istream&   in_;
    std::string     fileName_;

    std::string     line_;
    size_t          lineNo_;
    SamFileRecord               samFileRecord_;
    util::DelimitedLineParser   dl_;
};

} } // cgatools::mapping

#endif // CGA_TOOLS_SAM_READER_HPP_

