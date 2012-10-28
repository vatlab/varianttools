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
#include "cgatools/util/StringSet.hpp"
#include "cgatools/util/Exception.hpp"
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

namespace cgatools { namespace util {

    using namespace std;
    namespace ba = boost::algorithm;

    StringSet::StringSet(const std::string& csValues,
                         const std::string& csUniverse,
                         const std::string& notFoundMsg)
        : notFoundMsg_(notFoundMsg)
    {
        if ("" != csUniverse)
        {
            vector<string> parts;
            ba::split(parts, csUniverse, boost::is_any_of(","));
            BOOST_FOREACH(const string& part, parts)
            {
                universe_.insert(part);
            }
        }

        if ("" == csValues)
            return;

        if ("all" == csValues)
        {
            BOOST_FOREACH(const string& value, universe_)
            {
                this->insert(value);
            }
            return;
        }

        vector<string> parts;
        ba::split(parts, csValues, boost::is_any_of(","));
        BOOST_FOREACH(const string& part, parts)
        {
            this->insert(part);
        }
    }

    std::pair<StringSet::iterator, bool> StringSet::insert(const std::string& value)
    {
        if (!universe_.empty() && universe_.find(value) == universe_.end())
            throw Exception(notFoundMsg_+": "+value);
        return set<string>::insert(value);
    }

} } // cgatools::util
