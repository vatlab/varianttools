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
#include "cgatools/util/BaseUtil.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace cgatools { namespace reference {

    using std::string;
    using std::vector;
    using namespace cgatools::util;
    namespace bu = baseutil;

    CrrFileWriter::ChromosomeInfo::ChromosomeInfo()
        : circular_(false),
          fileOffset_(0),
          length_(0)
    {
    }

    CrrFileWriter::ChromosomeInfo::ChromosomeInfo(
        const std::string& name, bool circular, uint64_t fileOffset)
        : name_(name),
          circular_(circular),
          fileOffset_(fileOffset),
          length_(0)
    {
    }

    CrrFileWriter::CrrFileWriter(std::ostream* out)
        : out_(out),
          closed_(false)
    {
        writeHeader(0);
    }

    CrrFileWriter::~CrrFileWriter()
    {
        close();
    }

    void CrrFileWriter::close()
    {
        if (!closed_)
        {
            endChromosome();
            uint64_t chromosomeTable = out_->tellp();
            writeChromosomeTable();
            writeHeader(chromosomeTable);
        }
    }

    void CrrFileWriter::writeHeader(uint64_t chrTableOffset)
    {
        out_->seekp(0, std::ios_base::beg);
        out_->write("CRR\n", 4);
        writeBinaryInt(*out_, CrrFile::currentVersion());
        writeBinaryInt(*out_, uint64_t(0)); // filler for future use
        writeBinaryInt(*out_, uint64_t(0)); // filler for future use
        writeBinaryInt(*out_, chrTableOffset);
    }

    void CrrFileWriter::writeGuard()
    {
        writeBinaryInt(*out_, uint32_t(0xdeadbeef));
        writeBinaryInt(*out_, uint32_t(0));
        writeBinaryInt(*out_, uint32_t(0));
        writeBinaryInt(*out_, uint32_t(0xdeadbeef));
    }

    void CrrFileWriter::writeChromosomeTable()
    {
        writeBinaryUIntZC(*out_, chromosomes_.size());
        BOOST_FOREACH(const ChromosomeInfo& info, chromosomes_)
        {
            writeBinaryString(*out_, info.name_);
            writeBinaryBool(*out_, info.circular_);
            writeBinaryUIntZC(*out_, info.fileOffset_);
            util::Md5Digest digest = info.md5_.getDigest();
            out_->write(static_cast<const char*>(digest.data()), digest.size());
            writeBinaryUIntZC(*out_, info.length_);
            writeBinaryUIntZC(*out_, info.amb_.size());
            BOOST_FOREACH(const AmbiguousRegion& amb, info.amb_)
            {
                out_->put(amb.code_);
                writeBinaryUIntZC(*out_, amb.offset_);
                writeBinaryUIntZC(*out_, amb.length_);
            }
        }
    }

    void CrrFileWriter::newChromosome(const std::string& name, bool circular)
    {
        CGA_ASSERT(!closed_);

        if (chromosomes_.size() != 0)
            endChromosome();
        writeGuard();

        BOOST_FOREACH(const ChromosomeInfo info, chromosomes_)
        {
            if (info.name_ == name)
                throw Exception("repeated chromosome name: "+name);
        }

        chromosomes_.push_back(ChromosomeInfo(name, circular, out_->tellp()));
        packedBases_ = 0;
        packedBaseCount_ = 0;
    }

    void CrrFileWriter::endChromosome()
    {
        if (packedBaseCount_ > 0)
        {
            packedBases_ <<= 2 * (4-packedBaseCount_); // pad A's
            out_->put(packedBases_);
        }
        writeGuard();
    }

    void CrrFileWriter::addSequence(const std::string& sequence)
    {
        CGA_ASSERT(!closed_);
        CGA_ASSERT(0 != chromosomes_.size());

        for(size_t ii=0; ii<sequence.size(); ii++)
        {
            char base = sequence[ii];

            if (isspace(base))
                continue;

            if ('-' == base)
                continue;

            if (0 == bu::disambiguate(base))
                throw Exception("Unrecognized IUPAC code: "+string(1, sequence[ii]));

            addBase(sequence[ii]);
        }
    }

    void CrrFileWriter::addBase(char base)
    {
        ChromosomeInfo& info = chromosomes_.back();

        base = toupper(base);
        char ch = toupper(bu::disambiguate(base));

        // Add base to buffer.
        uint32_t val = bu::pack(ch);
        packedBases_ <<= 2;
        packedBases_ |= val;
        packedBaseCount_++;

        if (packedBaseCount_ >= 4)
        {
            // Flush buffer.
            out_->put(packedBases_);
            packedBases_ = 0;
            packedBaseCount_ = 0;
        }

        if (ch != base)
        {
            // Update ambiguity info.
            if (info.amb_.size() == 0 ||
                info.amb_.back().code_ != base ||
                info.amb_.back().offset_ + info.amb_.back().length_ != info.length_)
                info.amb_.push_back(AmbiguousRegion(base, info.length_, 1));
            else
                info.amb_.back().length_++;
        }

        // Update ChromosomeInfo.
        info.md5_.update(&base, 1);
        info.length_++;
    }

} } // cgatools::reference
