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

#ifndef CGATOOLS_UTIL_MD5_HPP_
#define CGATOOLS_UTIL_MD5_HPP_ 1

//! @file Md5.hpp
//! File containing Md5 class for computing md5 checksums.

#include "cgatools/core.hpp"
#include <string>
#include <boost/array.hpp>

namespace cgatools { namespace util {

    //! The MD5 result type.
    class Md5Digest
    {
    public:
        //! Create an Md5Digest with all 0's.
        Md5Digest();

        //! Create an Md5Digest with the given value. Here, val points
        //! to 16 bytes that are read by this function.
        Md5Digest(const void* val);

        //! Set the Md5Digest value. Here, val points to 16 bytes that
        //! are read by this function.
        void set(const void* val);

        //! Get the Md5Digest binary value. It is 16 bytes.
        const void* data() const;

        //! Get the length of the Md5Digest binary value. It is 16 bytes.
        size_t size() const
        {
            return 16;
        }

        //! Get the hexadecimal string representation of the Md5Digest,
        //! as is commonly used by programs like md5sum.
        std::string hex() const;

    private:
        boost::array<uint8_t, 16> digest_;
    };

    bool operator==(const Md5Digest& lhs, const Md5Digest& rhs);
    inline bool operator!=(const Md5Digest& lhs, const Md5Digest& rhs)
    {
        return ! (lhs == rhs);
    }

    //! A class to compute the MD5 hash of a sequence of bytes. The idea
    //! is you stream bytes to Md5Context::update(const void*, size_t),
    //! then get the digest.
    class Md5Context
    {
    public:
        //! Initialize a new MD5 context.
        Md5Context();

        //! Re-initialize the MD5 context.
        void init();

        //! Stream more bytes.
        void update(const void *buf, size_t length);

        //! Get the MD5 hash of the streamed bytes.
        Md5Digest getDigest() const;

        //! Get the MD5 hash of the streamed bytes, in hex format, as is
        //! commonly used by programs like md5sum.
        std::string hexDigest() const;

    private:
        void final();
        void transform();

        boost::array<uint32_t, 4> buf_;
        boost::array<uint32_t, 2> bits_;
        boost::array<uint8_t, 64> in_;
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_MD5_HPP_
