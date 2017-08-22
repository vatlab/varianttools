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
#include "cgatools/util/parse.hpp"
#include "cgatools/util/Exception.hpp"

#include <string>
#include <errno.h>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

namespace cgatools { namespace util {

    using std::string;

#if defined(CGA_USE_WIN_API)
    inline float strtof(_In_z_ const char * _Str, _Out_opt_ _Deref_post_z_ char ** _EndPtr) 
    {
        return (float)strtod(_Str, _EndPtr);
    }
#endif

    float  parseFloat (const char* first, const char* last)
    {
        if (first==last || *first==' ')
            throw Exception("failed to parse float: " + string(first, last));

        std::string stringCopy(first,last);

        char* lastParsed;
        errno = 0;
        float result = strtof(stringCopy.c_str(), &lastParsed);
        if (ERANGE == errno || lastParsed-stringCopy.c_str() != int(stringCopy.size()))
            throw Exception("failed to parse float: " + string(first, last));
        return result;
    }

    double parseDouble(const char* first, const char* last)
    {
        if (first==last || *first==' ')
            throw Exception("failed to parse double: " + string(first, last));

        std::string stringCopy(first,last);

        char* lastParsed;
        errno = 0;
        double result = strtod(stringCopy.c_str(), &lastParsed);
        if (ERANGE == errno || lastParsed-stringCopy.c_str() != int(stringCopy.size()))
            throw Exception("failed to parse double: " + string(first, last));
        return result;
    }

    // Parse a user-token-specified separated string
    //
    void parseTokenString(std::string const & inputString, char const * sep, 
        std::vector<std::string> & tokens, bool trimWhitespace, 
        bool compressEmpty) 
    {
        const boost::char_separator<char> separator(sep, "", 
            compressEmpty ? boost::drop_empty_tokens : boost::keep_empty_tokens);
        typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
        Tokenizer tokenizer(inputString, separator);
        Tokenizer::iterator it;
        tokens.clear();

        for(it=tokenizer.begin(); it!=tokenizer.end(); ++it) {
            std::string s = *it;
            if (trimWhitespace)
                boost::trim(s);
            tokens.push_back(s);
        }
    }

} } // cgatools::util
