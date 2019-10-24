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
#include "Files.hpp"
#include "Exception.hpp"
#include "Streams.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

using namespace std;

namespace cgatools { namespace util { namespace files {

    //-----------------------------------------------------------------------
    // File search
    //-----------------------------------------------------------------------

    bool findFiles( boost::filesystem::path const & directory,
        std::string const & filename,
        std::vector<boost::filesystem::path> & filesFound,
        bool recursive,
        bool includeDir,
        boost::regex_constants::syntax_option_type regExOptions)
    {
        bool completePaths(directory.is_complete());
        boost::regex re(filename, regExOptions);

        if (!exists(directory)) 
        {
            return(false);
        }
        try {
            boost::filesystem::directory_iterator end_itr; // default construction yields past-the-end
            for ( boost::filesystem::directory_iterator itr( directory ); itr != end_itr; ++itr )
            {
                if (!is_directory(itr->status())) {
                    // changed by Bo Peng from leaf() to filename() to support newer version of boost
                    if (boost::regex_search(itr->path().filename().c_str(), re))
                    {
                        if (completePaths) {
                            filesFound.push_back(boost::filesystem::system_complete(itr->path()));
                        }
                        else {
                            filesFound.push_back(itr->path());
                        }
                    }
                }
                else {
                    if (includeDir) {
                        // changed by Bo Peng from leaf() to filename() to support newer version of boost
                        if (boost::regex_search(itr->path().filename().c_str(), re)) {
                            if (completePaths) {
                                filesFound.push_back(boost::filesystem::system_complete(itr->path()));
                            }
                            else {
                                filesFound.push_back(itr->path());
                            }
                        }
                    }
                    if (recursive)
                    {
                        (void)findFiles(itr->path(), filename, filesFound, 
                            recursive, includeDir, regExOptions);
                    }
                }
            }
        }
        catch (...)
        {
            // issues such as premissions etc can result in exceptions but however can
            // still be ignored
        }
        return (0 != filesFound.size());
    }

    // returns a list of the file names in the dir
    std::vector<std::string> listDir(const boost::filesystem::path& dir,
        const std::string &pattern, size_t maxFiles)
    {
        CGA_ASSERT_L(0,maxFiles);
        std::vector<std::string> result;
        boost::regex re(pattern);

        for (boost::filesystem::directory_iterator it(dir), itEnd; it != itEnd; ++it)
        {
            boost::smatch what;
            // changed by Bo Peng from leaf() to path().filename() to support newer version of boost
            std::string fname(it->path().filename().string());
            if (boost::regex_match(fname, re)) {
                result.push_back(fname);
                if (--maxFiles==0)
                    break;
            }
        }

        return result;
    }

    void check_dir(const boost::filesystem::path& p)
    {
        if (!boost::filesystem::exists(p) || !boost::filesystem::is_directory(p))
        {
            throw util::Exception("not an existing directory: " + p.string());
        }
    }

    std::string findDataFile(const boost::filesystem::path& dir, const std::string& baseName)
    {
        static const char* exts[] = {".bz2", ".gz"};
        static const size_t extCount = sizeof(exts) / sizeof(exts[0]);

        for (size_t ii = 0; ii < extCount; ++ii)
        {
            boost::filesystem::path p = dir / (baseName + exts[ii]);
            // changed by Bo Peng from external_file_string() to native() to support newer version of boost
            if (InputStream::isReadable(p.native()))
                // changed by Bo Peng from external_file_string() to native() to support newer version of boost
                return p.native();
        }
        // changed by Bo Peng from external_file_string() to native() to support newer version of boost
        return (dir / baseName).native();
    }

    std::string findDataFileRegex(const boost::filesystem::path& dir, const std::string& baseName)
    {
        static const char* exts[] = {"\\.bz2", "\\.gz", ""};
        static const size_t extCount = sizeof(exts) / sizeof(exts[0]);

        for (size_t ii = 0; ii < extCount; ++ii)
        {
            vector<string> fns = listDir(dir, baseName+exts[ii], 2);
            if (fns.size() > 1)
                // changed by Bo Peng from external_file_string() to native() to support newer version of boost
                throw Exception("multiple files matching regex "+baseName+exts[ii]+
                                " in directory "+dir.native());
            if (1 == fns.size())
                // changed by Bo Peng from external_file_string() to native() to support newer version of boost
                return (dir / fns[0]).native();
        }
        return "";
    }

    std::string findDataFileRegexMulti(
        const boost::filesystem::path& dir,
        const std::vector<std::string>& possibleBaseNames,
        bool exnOnFail)
    {
        BOOST_FOREACH(const string& baseName, possibleBaseNames)
        {
            string fn = findDataFileRegex(dir, baseName);
            if (""  != fn)
                return fn;
        }

        CGA_ASSERT(possibleBaseNames.size() > 0);
        if (exnOnFail)
        {
            string nameList = boost::join(possibleBaseNames, " or ");
            // changed by Bo Peng from external_file_string() to native() to support newer version of boost
            throw Exception("failed to find "+nameList+" in "+dir.native());
        }
        return "";
    }

} } } // cgatools::util::files
