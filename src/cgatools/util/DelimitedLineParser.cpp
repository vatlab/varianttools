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
#include "cgatools/util/DelimitedLineParser.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/util/Exception.hpp"

namespace cgatools { namespace util {

    using boost::array;
    using std::string;

    void StringField::parse(const char* first, const char* last)
    {
        val_->assign(first, last);
    }


    void CharField::parse( const char* first, const char* last )
    {
        if (last != first+1)
            throw Exception("unable to parse char: "+string(first,last));
        *val_ = *first;
    }

    void StrandField::parse(const char* first, const char* last)
    {
        if (last == first+1)
        {
            if ('+' == *first)
            {
                *val_ = false;
                return;
            }
            else if ('-' == *first)
            {
                *val_ = true;
                return;
            }
        }
        throw Exception("unable to parse strand: "+string(first,last));
    }

    void SideField::parse( const char* first, const char* last )
    {
        if (last == first+1)
        {
            if ('L' == *first)
            {
                *val_ = 0;
                return;
            }
            else if ('R' == *first)
            {
                *val_ = 1;
                return;
            }
        }
        throw Exception("unable to parse side: "+string(first,last));
    }

    DelimitedLineParser::DelimitedLineParser(bool allowOptionalFieldsInLine)
        : allowOptionalFieldsInLine_(allowOptionalFieldsInLine)
    {
    }

    void DelimitedLineParser::parseLine(
        const std::string& line,
        char delimiter,
        EmptyFieldHandling emptyHandling,
        StrictnessChecking strictnessChecking)
    {
        parseLine(&line[0], &line[0]+line.size(), delimiter, emptyHandling, strictnessChecking);
    }

    void DelimitedLineParser::parseLine(
        const char* first, const char* last,
        char delimiter,
        EmptyFieldHandling emptyHandling,
        StrictnessChecking strictnessChecking)
    {
        size_t fieldCount = 0;
        while (first <= last)
        {
            if (SKIP_EMPTY_FIELDS == emptyHandling)
            {
                while (first < last && *first == delimiter)
                    ++first;

                if (first == last)
                    break;
            }

            const char* next = first;
            while (next < last && *next != delimiter)
                ++next;

            if (fields_.size() == fieldCount)
            {
                if (STRICT_CHECKING == strictnessChecking || next != first)
                    throw Exception("failed to parse line: too many fields");

                // Ignore extra trailing empty fields.
                first = next+1;
                continue;
            }

            try
            {
                fields_[fieldCount]->parse(first, next);
            }
            catch(const std::exception&)
            {
                throw;
            }

            fieldCount++;
            first = next+1;
        }
        
        if (fields_.size() != fieldCount && !allowOptionalFieldsInLine_)
            throw Exception("failed to parse line: not enough fields");
    }


    void DelimitedField::parse( const char* first, const char* last )
    {
        lineParser_.parseLine(first, last, delimiter_);
    }

} } // cgatools::util
