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

#ifndef CGATOOLS_UTIL_STRINGSET_HPP_
#define CGATOOLS_UTIL_STRINGSET_HPP_ 1

//! @file StringSet.hpp

#include "cgatools/core.hpp"
#include <string>
#include <set>

namespace cgatools { namespace util {

    //! StringSet convenience class.
    class StringSet : public std::set<std::string>
    {
    public:
        //! Constructs a StringSet.
        //! @param csValues A comma-separated list of values to
        //!                 insert. The speciall value "all" means to
        //!                 include all strings from csUniverse.
        //! @param csUniverse A comma-separated list of values allowed
        //!                   in this StringSet. If this is empty, then
        //!                   any value is allowed. If an item is
        //!                   inserted that is not allowed, an exception
        //!                   is thrown.
        //! @param notFoundMsg A message to include in the exception
        //!                    when an invalid item is inserted into
        //!                    this set.
        StringSet(const std::string& csValues,
                  const std::string& csUniverse,
                  const std::string& notFoundMsg);

        //! Inserts the string, checking against the allowed values
        //! first.
        std::pair<iterator, bool> insert(const std::string& value);

        bool contains(const std::string& value) const
        {
            return 0 != count(value);
        }

    private:
        std::set<std::string> universe_;
        std::string notFoundMsg_;
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_STRINGSET_HPP_
