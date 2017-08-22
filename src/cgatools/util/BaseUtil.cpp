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
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/Exception.hpp"
#include <algorithm>

namespace cgatools { namespace util { namespace baseutil {

    using std::string;

    //! Returns true if allele is consistent with fixed. Here,
    //! allele and fixed are sequences of IUPAC codes, and allele
    //! may have zero or more '?' characters.
    bool isConsistentToFixedLength(const std::string& fixed,  size_t fixedPos,  size_t fixedEnd,
                                   const std::string& allele, size_t allelePos, size_t alleleEnd);

    //! Returns true if the given subsequence within lhs and rhs are
    //! consistent. The subsequence must consist of IUPAC codes.
    bool matchFixedLength(const std::string& lhs,
                          const std::string& rhs,
                          size_t lhsPos,
                          size_t rhsPos,
                          size_t len);

    bool isValidBase(char iupacCode)
    {
        switch(iupacCode)
        {
        case 'a': case 'A':
        case 'c': case 'C':
        case 'g': case 'G':
        case 't': case 'T':
            return true;
        default:
            return false;
        }
    }

    bool isCalledSequence(const std::string& sequence)
    {
        return isCalledSequence(sequence, 0, sequence.size());
    }

    bool isCalledSequence(const std::string& sequence, size_t start, size_t end)
    {
        CGA_ASSERT(start <= end && end <= sequence.size());
        for(size_t ii=start; ii<end; ii++)
        {
            if (!isValidBase(sequence[ii]))
                return false;
        }
        return true;
    }

    bool isValidIupacCode(char iupacCode)
    {
        return 0 != BASE_COMPATIBILITY[uint8_t(iupacCode)];
    }

    uint32_t pack(char base)
    {
        switch(base)
        {
        case 'a': case 'A':
            return 0;
        case 'c': case 'C':
            return 1;
        case 'g': case 'G':
            return 2;
        case 't': case 'T':
            return 3;
        default:
            throw Exception("invalid base call: "+string(1, base));
        }
    }

    char unpack(uint32_t packedBase)
    {
        CGA_ASSERT(packedBase < 4);
        char packedBase2Base[] = { 'A', 'C', 'G', 'T' };
        return packedBase2Base[packedBase];
    }

    char disambiguate(char iupacCode)
    {
        switch (iupacCode)
        {
        case 'a': case 'A':
        case 'n': case 'N':
        case 'r': case 'R':
        case 'm': case 'M':
        case 'w': case 'W':
            return 'A';
        case 'c': case 'C':
        case 'y': case 'Y':
        case 'h': case 'H':
        case 'v': case 'V':
            return 'C';
        case 'g': case 'G':
        case 's': case 'S':
        case 'b': case 'B':
            return 'G';
        case 't': case 'T':
        case 'k': case 'K':
        case 'd': case 'D':
            return 'T';
        default:
            throw Exception("invalid IUPAC code: "+string(1, iupacCode));
        }
    }

    bool isConsistent(const std::string& lhs, const std::string& rhs)
    {
        return isConsistent(lhs, 0, lhs.size(), rhs, 0, rhs.size());
    }

    bool isConsistent(const std::string& lhs, size_t lhsStart, size_t lhsEnd,
                      const std::string& rhs, size_t rhsStart, size_t rhsEnd)
    {
        size_t lQ = lhs.find('?', lhsStart);
        if (lhsEnd <= lQ)
            return isConsistentToFixedLength(lhs, lhsStart, lhsEnd, rhs, rhsStart, rhsEnd);
        size_t rQ = rhs.find('?', rhsStart);
        if (rhsEnd <= rQ)
            return isConsistentToFixedLength(rhs, rhsStart, rhsEnd, lhs, lhsStart, lhsEnd);

        // Both have length no-calls. They are consistent if the
        // prefixes and suffixes are consistent.
        if (!matchFixedLength(lhs, rhs, lhsStart, rhsStart, std::min(lQ-lhsStart, rQ-rhsStart)))
            return false;

        size_t lQ2 = lhs.rfind('?', lhsEnd-1);
        CGA_ASSERT(string::npos != lQ2);
        size_t rQ2 = rhs.rfind('?', rhsEnd-1);
        CGA_ASSERT(string::npos != rQ2);
        size_t len = std::min(lhsEnd-lQ2-1, rhsEnd-rQ2-1);
        return matchFixedLength(lhs, rhs, lhsEnd-len, rhsEnd-len, len);
    }

    char complement(char iupacCode)
    {
        return BASE_COMPLEMENT[uint8_t(iupacCode)];
    }


    std::string reverseComplement(const std::string& sequence)
    {
        string result(sequence);
        for(size_t ii=0; ii<result.size(); ii++)
            result[ii] = complement(result[ii]);
        std::reverse(result.begin(), result.end());
        return result;
    }

    bool isConsistentToFixedLength(const std::string& fixed,  size_t fixedStart,  size_t fixedEnd,
                                   const std::string& allele, size_t alleleStart, size_t alleleEnd)
    {
        size_t fixedPos = fixedStart;
        size_t allelePos = alleleStart;
        for (;;)
        {
            size_t nextQ = allele.find('?', allelePos);
            if (alleleEnd <= nextQ)
            {
                // Better have a match from here to end of string.
                if ( (alleleStart == allelePos && fixedEnd-fixedStart != alleleEnd-alleleStart) ||
                     fixedEnd-fixedPos < alleleEnd-allelePos )
                    return false;

                fixedPos = fixedEnd - (alleleEnd-allelePos);
                if (!matchFixedLength(allele, fixed, allelePos, fixedPos, alleleEnd-allelePos))
                    return false;

                return true;
            }

            // There is a question mark. Find the first match.
            size_t len = nextQ - allelePos;
            size_t ii;
            for(ii=fixedPos; ii+len<=fixedEnd; ii++)
            {
                if (matchFixedLength(allele, fixed, allelePos, ii, len))
                    break;

                if (alleleStart == allelePos)
                    return false;
            }
            if (ii+len > fixedEnd)
                return false;
            fixedPos = ii+len;
            allelePos = nextQ;
            while (allelePos < alleleEnd && allele[allelePos] == '?')
                allelePos++;
        }
    }

    bool matchFixedLength(const std::string& lhs,
                                    const std::string& rhs,
                                    size_t lhsPos,
                                    size_t rhsPos,
                                    size_t len)
    {
        if (lhs.size() < lhsPos+len ||
            rhs.size() < rhsPos+len)
            return false;

        for(size_t ii=0; ii<len; ii++)
        {
            char lhsCh = lhs[lhsPos+ii];
            char rhsCh = rhs[rhsPos+ii];
            if (!isConsistent(lhsCh, rhsCh))
                return false;
        }
        return true;
    }


    char BASE_COMPATIBILITY[256] =
    {
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 00..07
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 08..0f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 10..17
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 18..1f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 20..27
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 28..2f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 30..37
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 38..3f
        0x0 /*.*/, 0x1 /*A*/, 0xe /*B*/, 0x2 /*C*/, 0xd /*D*/, 0x0 /*.*/, 0x0 /*.*/, 0x4 /*G*/, // 40..47
        0xb /*H*/, 0x0 /*.*/, 0x0 /*.*/, 0xc /*K*/, 0x0 /*.*/, 0x3 /*M*/, 0xf /*N*/, 0x0 /*.*/, // 48..4f
        0x0 /*.*/, 0x0 /*.*/, 0x5 /*R*/, 0x6 /*S*/, 0x8 /*T*/, 0x0 /*.*/, 0x7 /*V*/, 0x9 /*W*/, // 50..57
        0x0 /*.*/, 0xa /*Y*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 58..5f
        0x0 /*.*/, 0x1 /*a*/, 0xe /*b*/, 0x2 /*c*/, 0xd /*d*/, 0x0 /*.*/, 0x0 /*.*/, 0x4 /*g*/, // 60..67
        0xb /*h*/, 0x0 /*.*/, 0x0 /*.*/, 0xc /*k*/, 0x0 /*.*/, 0x3 /*m*/, 0xf /*n*/, 0x0 /*.*/, // 68..6f
        0x0 /*.*/, 0x0 /*.*/, 0x5 /*r*/, 0x6 /*s*/, 0x8 /*t*/, 0x0 /*.*/, 0x7 /*v*/, 0x9 /*w*/, // 70..77
        0x0 /*.*/, 0xa /*y*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 78..7f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 80..87
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 88..8f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 90..97
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 98..9f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // a0..a7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // a8..af
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // b0..b7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // b8..bf
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // c0..c7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // c8..cf
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // d0..d7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // d8..df
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // e0..e7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // e8..ef
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // f0..f7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // f8..ff
    };

    char BASE_COMPLEMENT[256] =
    {
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 00..07
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 08..0f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 10..17
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 18..1f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 20..27
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 28..2f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 30..37
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 38..3f
        0x0 /*.*/, 'T' /*A*/, 'V' /*B*/, 'G' /*C*/, 'H' /*D*/, 0x0 /*.*/, 0x0 /*.*/, 'C' /*G*/, // 40..47
        'D' /*H*/, 0x0 /*.*/, 0x0 /*.*/, 'M' /*K*/, 0x0 /*.*/, 'K' /*M*/, 'N' /*N*/, 0x0 /*.*/, // 48..4f
        0x0 /*.*/, 0x0 /*.*/, 'Y' /*R*/, 'S' /*S*/, 'A' /*T*/, 0x0 /*.*/, 'B' /*V*/, 'W' /*W*/, // 50..57
        0x0 /*.*/, 'R' /*Y*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 58..5f
        0x0 /*.*/, 'T' /*a*/, 'V' /*b*/, 'G' /*c*/, 'H' /*d*/, 0x0 /*.*/, 0x0 /*.*/, 'C' /*g*/, // 60..67
        'D' /*h*/, 0x0 /*.*/, 0x0 /*.*/, 'M' /*k*/, 0x0 /*.*/, 'K' /*m*/, 'N' /*n*/, 0x0 /*.*/, // 68..6f
        0x0 /*.*/, 0x0 /*.*/, 'Y' /*r*/, 'S' /*s*/, 'A' /*t*/, 0x0 /*.*/, 'B' /*v*/, 'W' /*w*/, // 70..77
        0x0 /*.*/, 'R' /*y*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 78..7f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 80..87
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 88..8f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 90..97
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // 98..9f
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // a0..a7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // a8..af
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // b0..b7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // b8..bf
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // c0..c7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // c8..cf
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // d0..d7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // d8..df
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // e0..e7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // e8..ef
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // f0..f7
        0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, 0x0 /*.*/, // f8..ff
    };

} } } // cgatools::util::baseutil
