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

#ifndef CGATOOLS_UTIL_PARSE_HPP_
#define CGATOOLS_UTIL_PARSE_HPP_ 1

//! @file parse.hpp
//! File containing general purpose functions parsing strings.

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"
#include <string>
#include <vector>
#include <string.h>
#include <boost/static_assert.hpp>

namespace cgatools { namespace util {

    template <typename Value, bool IsSigned>
    struct ValueParser
    {
        static Value parse(const char* first, const char* last)
        {
            // You need to use one of the specializations of this
            // class. This assertion should fail whenever this class is
            // instantiated.
            BOOST_STATIC_ASSERT(0 == sizeof(Value));
        }
    };

    template <typename Value>
    struct ValueParser<Value, true>
    {
        static Value parse(const char* first, const char* last)
        {
            if (first == last)
                throw Exception("failed to parse int: empty string");

            bool neg = '-' == *first;

            const char* current = first;
            if (neg)
                ++current;

            if (current == last)
                throw Exception("failed to parse int: empty string");

            Value val(0);
            Value maxVal = std::numeric_limits<Value>::max() / 10;
            for(; current<last; current++)
            {
                if (*current < '0' || *current > '9')
                    throw Exception("failed to parse int: "+std::string(first, last));
                Value digit = *current - '0';
                if (val >= maxVal)
                {
                    Value maxDigit = std::numeric_limits<Value>::max() % 10;
                    if (neg && current+1 == last)
                        maxDigit++;
                    if (val > maxVal || digit > maxDigit)
                        throw Exception("failed to parse int: overflow: "+std::string(first, last));
                }
                val = val*10 + digit;
            }
            if (neg)
                val = -val;
            return val;
        }
    };

    template <typename Value>
    struct ValueParser<Value, false>
    {
        static Value parse(const char* first, const char* last)
        {
            if (first == last)
                throw Exception("failed to parse int: empty string");

            Value val(0);
            Value maxVal = std::numeric_limits<Value>::max() / 10;
            for(const char* current=first; current<last; current++)
            {
                if (*current < '0' || *current > '9')
                    throw Exception("failed to parse int: "+std::string(first, last));
                Value digit = *current - '0';
                if (val >= maxVal)
                {
                    Value maxDigit = std::numeric_limits<Value>::max() % 10;
                    if (val > maxVal || digit > maxDigit)
                        throw Exception("failed to parse int: overflow: "+std::string(first, last));
                }
                val = val*10 + digit;
            }
            return val;
        }
    };

    float  parseFloat (const char* first, const char* last);
    double parseDouble(const char* first, const char* last);

    template <>
    struct ValueParser<float, true>
    {
        static float parse(const char* first, const char* last)
        {
            return parseFloat(first, last);
        }
    };

    template <>
    struct ValueParser<double, true>
    {
        static double parse(const char* first, const char* last)
        {
            return parseDouble(first, last);
        }
    };

    //! Parse the string to the given Value type. Calling
    //! parseValue<Value>(str) is much like calling
    //! boost::lexical_cast<Value>(str), with the following exceptions:
    //! - If Value is an 8-bit integer, boost may treat this as a string
    //!   to character conversion, while parseValue still treats it as a
    //!   string to integer conversion.
    //! - Integer conversion by parseValue is very strict about only
    //!   finding digits, possibly starting with a '-' character.
    //! - This is much faster than boost::lexical_cast.
    //! - This is only implemented for the following types:
    //!   - cgatools::int8_t
    //!   - cgatools::int16_t
    //!   - cgatools::int32_t
    //!   - cgatools::int64_t
    //!   - cgatools::uint8_t
    //!   - cgatools::uint16_t
    //!   - cgatools::uint32_t
    //!   - cgatools::uint64_t
    //!   - float
    //!   - double
    template <typename Value>
    Value parseValue(const char* first, const char* last)
    {
        return ValueParser<Value, std::numeric_limits<Value>::is_signed>::parse(first, last);
    }

    //! Parse the string to the given Value type.
    //! @see parseValue(const char*, const char*)
    template <typename Value> Value parseValue(const char* str)
    {
        const char* strEnd = str + strlen(str);
        return parseValue<Value>(str, strEnd);
    }

    //! Parse the string to the given Value type.
    //! @see parseValue(const char*, const char*)
    template <typename Value> Value parseValue(const std::string& str)
    {
        return parseValue<Value>(str.c_str(), str.c_str()+str.length());
    }

    //! Parse a user-token-specified separated string
    //! The function can be used if the input string format is not strictly defined
    void parseTokenString(std::string const & inputString, char const * sep, 
        std::vector<std::string> & tokens, 
        bool trimWhitespace = false, 
        bool compressEmpty = false);

} } // cgatools::util

#endif // CGATOOLS_UTIL_PARSE_HPP_
