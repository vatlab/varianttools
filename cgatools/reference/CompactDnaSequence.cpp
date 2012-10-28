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
#include "cgatools/reference/CompactDnaSequence.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace cgatools { namespace reference {

    using std::string;
    using std::vector;
    using namespace cgatools::util;
    namespace bu = baseutil;

    class AmbiguousRegionLessByOffset
    {
    public:
        bool operator()(const AmbiguousRegion& lhs, const AmbiguousRegion& rhs) const
        {
            return lhs.offset_ < rhs.offset_;
        }
    };

    CompactDnaSequence::CompactDnaSequence(const std::string& name,
                                           bool circular,
                                           const void* packedData,
                                           const util::Md5Digest& md5,
                                           size_t length,
                                           const std::vector<AmbiguousRegion> amb)
        : name_(name),
          circular_(circular),
          packedData_(static_cast<const uint8_t*>(packedData)),
          md5_(md5),
          length_(length),
          amb_(amb)
    {
        CGA_ASSERT(length_ >= 0);
    }

    int64_t CompactDnaSequence::fixCircularPos(int64_t pos) const
    {
        if (pos < 0 || pos >= length_)
        {
            if (circular_)
            {
                if (pos < 0)
                    pos += length_;
                if (pos >= length_)
                    pos -= length_;
                if (pos >= 0 && pos < length_)
                    return pos;
            }
            throw Exception("failed to get reference sequence: position out of range");
        }
        return pos;
    }

    std::string CompactDnaSequence::getSequence(int64_t pos, int64_t length) const
    {
        string result;
        appendSequence(result, pos, length);
        return result;
    }

    std::string CompactDnaSequence::getUnambiguousSequence(int64_t pos, int64_t length) const
    {
        string result;
        appendUnambiguousSequence(result, pos, length);
        return result;
    }

    void CompactDnaSequence::appendSequence(std::string& seq, int64_t pos, int64_t length) const
    {
        size_t startLength = seq.size();
        appendUnambiguousSequence(seq, pos, length);

        // Apply ambiguity codes.
        applyAmbiguity(&seq[startLength], pos, length);
    }

    void CompactDnaSequence::appendUnambiguousSequence(std::string& seq, int64_t pos, int64_t length) const
    {
        if (length < 0 || length > length_)
            throw Exception("failed to get reference sequence: position out of range");
        seq.reserve(seq.size()+length);
        for(int64_t ii=0; ii<length; ii++)
            seq.push_back(getUnambiguousBase(pos+ii));
    }

    char CompactDnaSequence::getBase(int64_t pos) const
    {
        char base = getUnambiguousBase(pos);
        applyAmbiguity(&base, pos, 1);
        return base;
    }

    char CompactDnaSequence::getUnambiguousBase(int64_t pos) const
    {
        pos = fixCircularPos(pos);
        size_t offset = pos / 4;
        size_t shift = 6 - ( (pos & 0x3) << 1 );
        uint32_t val = (packedData_[offset] >> shift) & 0x3;
        return bu::unpack(val);
    }

    size_t CompactDnaSequence::extendLeftBy3Mers(size_t pos, size_t count) const
    {
        if (0 == count)
            return pos;

        if (pos < 3)
            return 0;

        // kmers is a bitmap
        boost::uint64_t kmers = 0;
        boost::uint32_t kmer = 0;
        for(size_t ii=0; ii<2; ii++)
        {
            pos--;
            char base = getUnambiguousBase(pos);
            kmer = (kmer << 2) | bu::pack(base);
        }
        while (pos > 0)
        {
            pos--;
            char base = getUnambiguousBase(pos);
            kmer = ( (kmer << 2) & 0x3f ) | bu::pack(base);
            uint64_t mask = 1;
            mask <<= kmer;
            if ( 0 == (kmers & mask) )
            {
                count--;
                if (0 == count)
                    break;
            }
            kmers |= mask;
        }
        return pos;
    }

    size_t CompactDnaSequence::extendRightBy3Mers(size_t pos, size_t count) const
    {
        if (0 == count)
            return pos;

        if (size_t(length_) < pos+3)
            return length_;

        // kmers is a bitmap
        boost::uint64_t kmers = 0;
        boost::uint32_t kmer = 0;
        for(size_t ii=0; ii<2; ii++)
        {
            char base = getUnambiguousBase(pos);
            kmer = (kmer << 2) | bu::pack(base);
            pos++;
        }
        while (pos < size_t(length_))
        {
            char base = getUnambiguousBase(pos);
            kmer = ( (kmer << 2) & 0x3f ) | bu::pack(base);
            pos++;
            uint64_t mask = 1;
            mask <<= kmer;
            if ( 0 == (kmers & mask) )
            {
                count--;
                if (0 == count)
                    break;
            }
            kmers |= mask;
        }
        return pos;
    }

    void CompactDnaSequence::validate() const
    {
        const int64_t batchSize = 1024;
        Md5Context md5;
        string sequence;
        for(int64_t pos=0; pos<length_; pos+=batchSize)
        {
            int64_t posEnd = std::min(pos+batchSize, length_);
            sequence.clear();
            appendSequence(sequence, pos, posEnd-pos);
            md5.update(sequence.c_str(), sequence.length());
        }
        if (md5.getDigest() != md5_)
            throw Exception("reference validation failed: md5 mismatch for chromosome: "+name_);
    }

    void CompactDnaSequence::applyAmbiguity(char* seq, int64_t pos, int64_t length) const
    {
        if (0 == length)
            return;

        pos = fixCircularPos(pos);
        while (pos+length > length_)
        {
            applyAmbiguity(seq, pos, length_-pos);
            length -= length_-pos;
            seq += length_-pos;
            pos = 0;
        }

        AmbiguousRegion dummy('N', pos, 0);
        vector<AmbiguousRegion>::const_iterator iter =
            std::lower_bound(amb_.begin(), amb_.end(), dummy, AmbiguousRegionLessByOffset());
        if (iter != amb_.begin())
            --iter;
        while (iter != amb_.end() && iter->offset_ < pos+length)
        {
            if (iter->offset_+iter->length_ > pos)
            {
                int64_t maxStartPos = std::max(pos, int64_t(iter->offset_));
                int64_t minEndPos   = std::min(pos+length, int64_t(iter->offset_+iter->length_));
                std::fill(seq+maxStartPos-pos, seq+minEndPos-pos, iter->code_);
            }

            ++iter;
        }
    }

} } // cgatools::reference
