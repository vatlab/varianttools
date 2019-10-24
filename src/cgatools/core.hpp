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

#ifndef CGATOOLS_CORE_HPP_
#define CGATOOLS_CORE_HPP_ 1

//! @file core.hpp
//! First header included by all source files. Intended for generally
//! useful preprocessor definitions and typedefs.

#if defined(_WIN32) || defined(_WIN64)
    //! Macro that specifies whether to use Win32/64 API vs. something
    //! else.  This can be tested, for example, to switch between
    //! CloseHandle() and close() calls when closing a system-level file
    //! handle.
    #ifndef CGA_USE_WIN_API
        #define CGA_USE_WIN_API 1
    #endif

    //! On Windows, min and max are macros by default, unless you use this
    //! #define:
    #ifndef NOMINMAX
        #define NOMINMAX
    #endif
#endif

//! Make sure boost defines UINT64_C, etc.
#ifndef __STDC_CONSTANT_MACROS
    #define __STDC_CONSTANT_MACROS 1
#endif

// Use version 2 of the Boost Filesystem library.
// This stopped being the default as of Boost 1.44 and
// will no longer be supported in Boost 1.48.

// Changed by Bo Peng from 2 to 3 to support newer version of boost
#define BOOST_FILESYSTEM_VERSION 3

#include <boost/cstdint.hpp>
#include <limits>

//! Namespace containing the libcgatools API.
namespace cgatools {

    typedef boost::int8_t   int8_t;  //!< 8-bit int.
    typedef boost::int16_t  int16_t; //!< 16-bit int.
    typedef boost::int32_t  int32_t; //!< 32-bit int.
    typedef boost::int64_t  int64_t; //!< 64-bit int.

    typedef boost::uint8_t  uint8_t;  //!< 8-bit unsigned int.
    typedef boost::uint16_t uint16_t; //!< 16-bit unsigned int.
    typedef boost::uint32_t uint32_t; //!< 32-bit unsigned int.
    typedef boost::uint64_t uint64_t; //!< 64-bit unsigned int.

} // cgatools

#endif // CGATOOLS_CORE_HPP_
