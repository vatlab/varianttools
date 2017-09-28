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

#ifndef CGATOOLS_UTIL_FILES_HPP_
#define CGATOOLS_UTIL_FILES_HPP_ 1

//! @file Files.hpp
//! File containing file system service routines

#include "cgatools/core.hpp"

#include <vector>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

namespace cgatools { namespace util { namespace files {

    //! returns a list of the file names matching the given regular expression file name pattern
    bool findFiles( boost::filesystem::path const & directory,
        std::string const & filename,
        std::vector<boost::filesystem::path> & filesFound,
        bool recursive = false,
        bool includeDir = false,
        boost::regex_constants::syntax_option_type regExOptions = boost::regex_constants::basic);

    //! returns a list of the file names in the dir, doesn't support recursive scanning
    std::vector<std::string> listDir(const boost::filesystem::path& dir,
        const std::string &pattern, size_t maxFiles = std::numeric_limits<size_t>::max());

    //! checks if the specifide directory exists
    void check_dir(const boost::filesystem::path& p);

    //! check existence and returns path to the specified data file. 
    //! The function tries first the compressed versions (baseName+'.bz2) or (baseName+'.gz')
    std::string findDataFile(const boost::filesystem::path& dir, const std::string& baseName);

    //! check existence and returns path to the specified data file. 
    //! The function tries first the compressed versions (baseName+'.bz2') or (baseName+'.gz')
    //! The first matching file is returned, or if multiple files match
    //! a single regex, throws an exception. On failure, returns "".
    std::string findDataFileRegex(const boost::filesystem::path& dir, const std::string& baseName);

    //! Tries to find existing files matching the given base name
    //! regex. If no match is found, returns empty string or throws an
    //! exception, depending on value of exnOnFail.
    std::string findDataFileRegexMulti(
        const boost::filesystem::path& dir,
        const std::vector<std::string>& possibleBaseNames,
        bool exnOnFail);

} } } // cgatools::util::files

#endif // CGATOOLS_UTIL_FILES_HPP_
