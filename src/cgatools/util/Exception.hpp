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

#ifndef CGATOOLS_UTIL_EXCEPTION_HPP_
#define CGATOOLS_UTIL_EXCEPTION_HPP_ 1

//! @file Exception.hpp
//! File containing Exception class and assertion macros.

#include "cgatools/core.hpp"
#include <cstdlib>
#include <exception>
#include <string>
#include <iostream>
#include <sstream>

namespace cgatools { namespace util {

    //! Assertion that always runs. On failure, crash.
    #define CGA_ASSERT_MSG(expr, message)                               \
        do                                                              \
        {                                                               \
            if (!(expr)) {                                              \
                std::cerr << "assert failed: " <<                       \
                    __FILE__ << ":" <<                                  \
                    __LINE__ << ": " <<                                 \
                    #expr   << std::endl <<                             \
                    message << std::endl;                               \
                volatile char *p = 0, ch = *p;                          \
                if (ch == '\0') *p = ch; p++; ch = *p;                  \
                std::exit(1);                                           \
            }                                                           \
        } while (0)

    #define CGA_WARN_MSG(expr, message)                                 \
        do                                                              \
        {                                                               \
            if (!(expr)) {                                              \
                std::cerr << "Warning. Assertion failed: " <<           \
                    __FILE__ << ":" <<                                  \
                    __LINE__ << ": " <<                                 \
                    #expr   << std::endl <<                             \
                    message << std::endl;                               \
            }                                                           \
        } while (0)

    #define CGA_ASSERT(expr) CGA_ASSERT_MSG(expr, "")

    #ifdef NDEBUG

        //! Debug-only assert.
        #define CGA_DBG_ASSERT(expr) (void(0))

    #else // NDEBUG

        //! Debug-only assert.
        #define CGA_DBG_ASSERT(expr) CGA_ASSERT(expr)

    #endif // NDEBUG

    //! A basic exception type.
    class Exception: public std::exception
    {
    public:
        Exception(const std::string& message);
        virtual ~Exception() throw () {}

        const char* what() const throw ();

    private:
        std::string message_;
    };

    //prints an expression and its value into an output stream
    #define CGA_VOUT(var) "{"<<#var<<"="<<var<<"}"
    #define CGA_ASSERT_EQ(var1,var2) CGA_ASSERT_MSG((var1)==(var2), CGA_VOUT(var1)<<CGA_VOUT(var2))
    #define CGA_ASSERT_NEQ(var1,var2) CGA_ASSERT_MSG((var1)!=(var2), CGA_VOUT(var1)<<CGA_VOUT(var2))
    #define CGA_ASSERT_L(var1,var2) CGA_ASSERT_MSG((var1)<(var2), CGA_VOUT(var1)<<CGA_VOUT(var2))
    #define CGA_ASSERT_LE(var1,var2) CGA_ASSERT_MSG((var1)<=(var2), CGA_VOUT(var1)<<CGA_VOUT(var2))

    //throws an exception. 
    #define CGA_ERROR_EX(m_message) \
    { \
        std::stringstream stx; \
        stx << m_message; \
        throw cgatools::util::Exception(stx.str()); \
    }

} } // cgatools::util

#endif // CGATOOLS_UTIL_EXCEPTION_HPP_

