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

#ifndef CGATOOLS_MAINPAGE_HPP_
#define CGATOOLS_MAINPAGE_HPP_ 1

//! @file mainpage.hpp
//! File containing documentation for the doxygen main page.

#include "cgatools/core.hpp"

//! Namespace containing the libcgatools API.
namespace cgatools {

    //! Namespace for general-purpose utilities.
    namespace util
    {
        //! DNA nucleotide (base) utility functions.
        namespace baseutil { }
    }

    //! Namespace for reference-specific utilities.
    namespace reference { }

    //! Namespace for variant file processing utilities.
    namespace variants { }

    //! Namespace for command line utility classes.
    namespace command { }

} // cgatools

//! @mainpage cgatools API documentation
//! The cgatools library provides APIs for analysis of genomes assembled
//! by Complete Genomics, which are used by the Complete Genomics
//! Analysis Tools (cgatools).

#endif // CGATOOLS_MAINPAGE_HPP_
