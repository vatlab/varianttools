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
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/util/Exception.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace cgatools { namespace reference {

    using std::string;
    using std::vector;
    using namespace cgatools::util;

    CrrFile::CrrFile()
        : file_(new boost::iostreams::mapped_file_source())
    {
    }

    CrrFile::CrrFile(const std::string& path)
        : file_(new boost::iostreams::mapped_file_source())
    {
        open(path);
    }

    void CrrFile::open(const std::string& path)
    {
        // Read metadata.
        uint64_t chrTableOffset;
        InputStream in(path);
        readHeader(path, in, &chrTableOffset);

        vector<Md5Digest> digests;
        vector<CrrFileWriter::ChromosomeInfo> chromosomeTable;
        in.seekg(chrTableOffset, std::ios_base::beg);
        readChromosomeTable(in, &chromosomeTable, &digests);

        in.close();

        file_->open(path);
        chromosomes_.clear();

        const char* memory = file_->data();
        for(size_t ii=0; ii<chromosomeTable.size(); ii++)
        {
            const CrrFileWriter::ChromosomeInfo& info = chromosomeTable[ii];
            chromosomes_.push_back(CompactDnaSequence(info.name_,
                                                      info.circular_,
                                                      memory + info.fileOffset_,
                                                      digests[ii],
                                                      info.length_,
                                                      info.amb_));
        }
    }

    const std::vector<CompactDnaSequence>& CrrFile::listChromosomes() const
    {
        return chromosomes_;
    }

    std::vector<Range> CrrFile::listContigs(uint32_t minGapLength) const
    {
        vector<Range> result;
        for(size_t ii=0; ii<chromosomes_.size(); ii++)
        {
            uint32_t pos=0;
            BOOST_FOREACH(const AmbiguousRegion& amb, chromosomes_[ii].getAmbiguousRegions())
            {
                if ('N' == amb.code_ && amb.length_ >= minGapLength)
                {
                    if (pos != amb.offset_)
                        result.push_back(Range(ii, pos, amb.offset_));
                    pos = amb.offset_ + amb.length_;
                }
            }
            if (pos != chromosomes_[ii].length())
                result.push_back(Range(ii, pos, chromosomes_[ii].length()));
        }
        return result;
    }

    std::string CrrFile::getSequence(const Range& range) const
    {
        if (range.chromosome_ >= chromosomes_.size())
            throw Exception( (boost::format("unrecognized chromosome id: %d") % range.chromosome_)
                             .str());
        return chromosomes_[range.chromosome_].getSequence(range.begin_, range.length());
    }

    char CrrFile::getBase(const Location& loc) const
    {
        if (loc.chromosome_ >= chromosomes_.size())
            throw Exception( (boost::format("unrecognized chromosome id: %d") % loc.chromosome_).str());
        return chromosomes_[loc.chromosome_].getBase(loc.offset_);
    }

    uint16_t CrrFile::getChromosomeId(const std::string& chromosomeName) const
    {
        for(size_t ii=0; ii<chromosomes_.size(); ii++)
        {
            const std::string & name = chromosomes_[ii].getName();
            if (chromosomeName == name || 
                // add by Bo Peng to handle chrX and X name difference
                (name[0] == 'c' && name[1] == 'h' && name[2] == 'r' && name.compare(3, std::string::npos, chromosomeName)==0) ||
                (chromosomeName[0] == 'c' && chromosomeName[1] == 'h' && chromosomeName[2] == 'r' && chromosomeName.compare(3, std::string::npos, name)==0))
                return ii;
        }
        // added by Bo Peng for better error message.
        std::string names = "";
        for(size_t ii=0; ii<chromosomes_.size(); ii++)
        {
            names += (ii == 0 ? " " : ", ") + chromosomes_[ii].getName();
        }

        throw Exception("unrecognized chromosome name: " + chromosomeName + " (available: " + names + ")");
    }

    void CrrFile::validate() const
    {
        for(size_t ii=0; ii<chromosomes_.size(); ii++)
        {
            chromosomes_[ii].validate();
        }
    }

    void CrrFile::readHeader(const std::string& path, std::istream& in, uint64_t* chrTableOffset)
    {
        char buf[4];
        in.read(buf, 4);
        if (0 != memcmp(buf, "CRR\n", 4))
            throw Exception("failed to open reference "+path+": not a crr file");

        uint32_t version;
        readBinaryInt(in, &version);
        if (version != currentVersion())
            throw Exception("failed to open reference "+path+": version mismatch");

        uint64_t filler;
        readBinaryInt(in, &filler);
        readBinaryInt(in, &filler);

        readBinaryInt(in, chrTableOffset);
    }

    void CrrFile::readChromosomeTable(std::istream& in,
                                      std::vector<CrrFileWriter::ChromosomeInfo>* pInfo,
                                      std::vector<Md5Digest>* pDigests)
    {
        vector<CrrFileWriter::ChromosomeInfo>& info = *pInfo;
        vector<Md5Digest>& digests = *pDigests;

        size_t count;
        readBinaryUIntZC(in, &count);
        info.resize(count);
        digests.resize(count);

        for(size_t ii=0; ii<count; ii++)
        {
            readBinaryString(in, &info[ii].name_);
            readBinaryBool(in, &info[ii].circular_);
            readBinaryUIntZC(in, &info[ii].fileOffset_);

            char buf[16];
            in.read(buf, 16);
            digests[ii].set(buf);

            readBinaryUIntZC(in, &info[ii].length_);

            size_t ambCount;
            readBinaryUIntZC(in, &ambCount);
            info[ii].amb_.resize(ambCount);

            for(size_t jj=0; jj<ambCount; jj++)
            {
                info[ii].amb_[jj].code_ = in.get();
                if (!in.good())
                    throw Exception("failed to open reference: unexpected eof");
                readBinaryUIntZC(in, &info[ii].amb_[jj].offset_);
                readBinaryUIntZC(in, &info[ii].amb_[jj].length_);
            }
        }
    }

    uint32_t CrrFile::currentVersion()
    {
        return 1;
    }

} } // cgatools::reference
