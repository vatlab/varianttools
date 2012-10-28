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

#ifndef CGATOOLS_UTIL_BASEUTIL_HPP_
#define CGATOOLS_UTIL_BASEUTIL_HPP_ 1

//! @file BaseUtil.hpp
//! File containing definition of the BaseUtil class.

#include "cgatools/core.hpp"
#include <string>

namespace cgatools { namespace util { namespace baseutil {

    //! @cond PRIVATE
    extern char BASE_COMPATIBILITY[256];
    extern char BASE_COMPLEMENT[256];
    //! @endcond

    //! Returns true for A, a, C, c, G, g, T, or t.
    bool isValidBase(char base);

    //! Returns true if all the bases in the sequence are valid.
    bool isCalledSequence(const std::string& sequence);

    //! Returns true if all the bases in the subsequence are valid.
    bool isCalledSequence(const std::string& sequence, size_t start, size_t end);

    //! Returns true if iupacCode is a valid IUPAC (base ambiguity)
    //! code.
    bool isValidIupacCode(char iupacCode);

    //! Returns an integer corresponding to the given base, or
    //! throws an exception if the base is invalid. Valid base
    //! calls:
    //! - A or a -> 0
    //! - C or c -> 1
    //! - G or g -> 2
    //! - T or t -> 3
    uint32_t pack(char base);

    //! Returns the unpacked base call. This can be used with
    //! BaseUtil::pack(char), such that
    //! BaseUtil::unpack(BaseUtil::pack(base)) is equivalent to
    //! toupper(base).
    char unpack(uint32_t packedBase);

    //! Returns an unambiguous base call that is consistent with the
    //! given iupacCode.
    char disambiguate(char iupacCode);

    //! Returns true if the given IUPAC codes are consistent. Two
    //! IUPAC codes are considered consistent if there is some base
    //! A, C, G, or T such that both codes are consistent with that
    //! base.
    inline bool isConsistent(char lhs, char rhs)
    {
        return 0 != ( BASE_COMPATIBILITY[uint8_t(lhs)] &
                      BASE_COMPATIBILITY[uint8_t(rhs)] );
    }

    //! Returns true if lhs and rhs are consistent. Here, lhs and
    //! rhs are sequences of IUPAC codes, and in addition they may
    //! have zero or more '?' characters to indicate an unknown
    //! sequence of zero or more bases.
    bool isConsistent(const std::string& lhs, const std::string& rhs);

    //! Returns true if lhs and rhs are consistent for the given range
    //! of posistions. Here, lhs and rhs are sequences of IUPAC codes,
    //! and in addition they may have zero or more '?' characters to
    //! indicate an unknown sequence of zero or more bases.
    bool isConsistent(const std::string& lhs, size_t lhsStart, size_t lhsEnd,
                      const std::string& rhs, size_t rhsStart, size_t rhsEnd);

    //! Returns the complement of the given IUPAC code. For bases,
    //! the complements are as follows:
    //! - A -> T
    //! - C -> G
    //! - G -> C
    //! - T -> A
    //! An ambiguous IUPAC code's complement is compatible with the
    //! complements of all the bases the original IUPAC code is
    //! compatible with.
    char complement(char iupacCode);

    //! Returns the reverse complement of the given sequence of
    //! IUPAC codes.
    std::string reverseComplement(const std::string& sequence);

} } } // cgatools::util::baseutil

#endif // CGATOOLS_UTIL_BASEUTIL_
