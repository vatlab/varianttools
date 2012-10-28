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

#ifndef CGATOOLS_REFERENCE_COMPACTDNASEQUENCE_HPP_
#define CGATOOLS_REFERENCE_COMPACTDNASEQUENCE_HPP_ 1

//! @file CompactDnaSequence.hpp
//! File containing definitions CompactDnaSequence class.

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/Md5.hpp"

#include <vector>
#include <string>

namespace cgatools { namespace reference {

    //! Struct to describe an ambiguous region of reference, within a
    //! chromosome.
    struct AmbiguousRegion
    {
        AmbiguousRegion()
            : code_('N'),
              offset_(0),
              length_(0)
        {
        }

        AmbiguousRegion(char code, uint32_t pos, uint32_t length)
            : code_(code),
              offset_(pos),
              length_(length)
        {
        }

        char code_;       //!< The IUPAC code for the region.
        uint32_t offset_; //!< The offset of the ambiguous region.
        uint32_t length_; //!< The length in bases of the ambiguous region.
    };

    //! Class to describe the DNA sequence of a chromosome, in a compact
    //! manner. Used internally by CrrFile class.
    class CompactDnaSequence
    {
    public:
        CompactDnaSequence(const std::string& name,
                           bool circular,
                           const void* packedData,
                           const util::Md5Digest& md5,
                           size_t length,
                           const std::vector<AmbiguousRegion> amb);

        //! Return the sequence of IUPAC codes for the chromosome as if
        //! by repeatedly calling CompactDnaSequence::getBase().
        std::string getSequence            (int64_t pos, int64_t length) const;

        //! Return an unambiguous sequence of base calls as if by
        //! repeatedly calling
        //! CompactDnaSequence::getUnambiguousBase().
        std::string getUnambiguousSequence(int64_t pos, int64_t length) const;

        //! Append the sequence of IUPAC codes as if by repeatedly
        //! calling CompactDnaSequence::getBase().
        void appendSequence            (std::string& seq, int64_t pos, int64_t length) const;

        //! Append an unambiguous sequence of base calls as if by
        //! repeatedly calling
        //! CompactDnaSequence::getUnambiguousBase().
        void appendUnambiguousSequence(std::string& seq, int64_t pos, int64_t length) const;

        //! Get the IUPAC code for this chromosome at position pos. For
        //! circular chromosomes, pos is allowed to range from -length
        //! to 2*length-1.
        char getBase(int64_t pos) const;

        //! Return an unambiguous base call for this chromosome at
        //! position pos (as by util::BaseUtil::disambiguate(char)) that
        //! is consistent with the IUPAC code for the chromosome at this
        //! position. For circular chromosomes, pos is allowed to range
        //! from -length to 2*length-1.
        char getUnambiguousBase(int64_t pos) const;

        //! Return pos, extended to the left until it has passed by
        //! count distinct 3-mers of unambiguous reference
        //! sequence. This function stops at the chromosome end, even
        //! for circular chromosomes.
        size_t extendLeftBy3Mers (size_t pos, size_t count) const;

        //! Return pos, extended to the right until it has passed by
        //! count distinct 3-mers of unambiguous reference
        //! sequence. This function stops at the chromosome end, even
        //! for circular chromosomes.
        size_t extendRightBy3Mers(size_t pos, size_t count) const;

        //! Verify that the md5s recorded in the crr file metadata are
        //! the same as the md5s produced by re-computing them on the
        //! data.
        void validate() const;

        //! Return the name of this chromosome.
        const std::string& getName() const
        {
            return name_;
        }

        //! Return whether this chromosome is circular.
        bool isCircular() const
        {
            return circular_;
        }

        //! Return the md5 digest of the chromosome's sequence. In
        //! particular, this is the md5 of the IUPAC codes of the
        //! chromosome, converted to upper case.
        const util::Md5Digest& getMd5Digest() const
        {
            return md5_;
        }

        //! Return the length in bases of the chromosome.
        size_t length() const
        {
            return size_t(length_);
        }

        //! Return the list of AmbiguousRegion for this chromosome, in
        //! order by position.
        const std::vector<AmbiguousRegion>& getAmbiguousRegions() const
        {
            return amb_;
        }

    private:
        inline int64_t fixCircularPos(int64_t pos) const;
        void applyAmbiguity(char* seq, int64_t pos, int64_t length) const;

        std::string name_;
        bool circular_;
        const uint8_t* packedData_;
        util::Md5Digest md5_;
        int64_t length_;
        std::vector<AmbiguousRegion> amb_;
    };

} } // cgatools::reference

#endif // CGATOOLS_REFERENCE_COMPACTDNASEQUENCE_HPP_
