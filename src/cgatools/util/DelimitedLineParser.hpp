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

#ifndef CGATOOLS_UTIL_DELIMITEDLINEPARSER_HPP_
#define CGATOOLS_UTIL_DELIMITEDLINEPARSER_HPP_ 1

//! @file DelimitedLineParser.hpp
//! File containing definitions of DelimitedLineParser and general
//! purpose DelimitedFieldParsers.

#include "cgatools/core.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/util/Exception.hpp"

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>

#include <cstring>

namespace cgatools { namespace util {

    //! Class that parses one field of a delimited line.
    class DelimitedFieldParser
    {
    public:
        DelimitedFieldParser(const std::string& name)
            : name_(name)
        {
        }

        virtual ~DelimitedFieldParser()
        {
        }

        const std::string& getName() const
        {
            return name_;
        }

        //! Called once per field by DelimitedLineParser.
        //! @param first The first char of the field.
        //! @param last  One past the end of the field. When parse is
        //!              called, *last is the null character.
        virtual void parse(const char* first, const char* last) = 0;

    protected:
        //! The name of the field.
        std::string name_;
    };

    //! A no-op DelimitedFieldParser.
    class IgnoreField : public DelimitedFieldParser
    {
    public:
        IgnoreField(const std::string& name)
            : DelimitedFieldParser(name)
        {
        }

        void parse(const char* first, const char* last)
        {
        }
    };

    //! A DelimitedFieldParser that uses cgatools::util::parseValue to
    //! parse its value.
    template <typename Value>
    class ValueField : public DelimitedFieldParser
    {
    public:
        ValueField(const std::string& name, Value* val)
            : DelimitedFieldParser(name),
              val_(val)
        {
        }

        void parse(const char* first, const char* last)
        {
            for(size_t ii=0; ii<exceptions_.size(); ii++)
            {
                size_t strSize = exceptions_[ii].first.length();
                if (size_t(last-first) == strSize && 
                    0 == std::strncmp(first, exceptions_[ii].first.c_str(), strSize))
                {
                    *val_ = exceptions_[ii].second;
                    return;
                }
            }

            *val_ = parseValue<Value>(first, last);
        }

        //! Add an exception, such that if the field equals key, set the
        //! value to val.
        ValueField<Value>& exception(const std::string& key, const Value& val)
        {
            exceptions_.push_back(make_pair(key, val));
            return *this;
        }

    private:
        Value* val_;
        std::vector< std::pair<std::string, Value> > exceptions_;
    };

    //! DelimitedFieldParser that parses a variable length list of values
    //! separated by its own delimiter (distinct from the field delimiter).
    //! Each token is assumed to be of the same type and is parsed by
    //! cgatools::util::parseValue. Results are stored in an std::vector.
    template <typename Value>
    class ValueVectorField : public DelimitedFieldParser
    {
    public:
        ValueVectorField(const std::string& name, char delimiter, std::vector<Value>* val)
            : DelimitedFieldParser(name),
              delimiter_(delimiter),
              val_(val)
        {
        }

        void parse(const char* first, const char* last)
        {
            val_->clear();

            const char* ii = first;
            while (ii < last)
            {
                const char* jj;
                for (jj = ii; jj < last && *jj != delimiter_; ++jj) ;

                Value v = parseValue<Value>(ii, jj);
                val_->push_back(v);

                ii = jj + 1;
            }
        }
    private:
        char delimiter_;
        std::vector<Value>* val_;
    };

    //! A DelimitedFieldParser that records a copy of the field.
    class StringField : public DelimitedFieldParser
    {
    public:
        StringField(const std::string& name, std::string* val)
            : DelimitedFieldParser(name),
              val_(val)
        {
        }

        void parse(const char* first, const char* last);

    private:
        std::string* val_;
    };

    //! A DelimitedFieldParser that records a character (not an unsigned integer).
    class CharField : public DelimitedFieldParser
    {
    public:
        CharField(const std::string& name, char* val)
            : DelimitedFieldParser(name),
            val_(val)
        {
        }

        void parse(const char* first, const char* last);

    private:
        char* val_;
    };

    class StrandField : public DelimitedFieldParser
    {
    public:
        StrandField(const std::string& name, bool* val)
            : DelimitedFieldParser(name),
              val_(val)
        {
        }

        void parse(const char* first, const char* last);

    private:
        bool* val_;
    };

    class SideField : public DelimitedFieldParser
    {
    public:
        SideField(const std::string& name, uint8_t* val)
            : DelimitedFieldParser(name),
            val_(val)
        {
        }

        void parse(const char* first, const char* last);

    private:
        uint8_t* val_;
    };

    //! A class that parses delimited lines.
    class DelimitedLineParser
    {
        friend class DelimitedFile;
    public:
        //! Enumeration to describe how empty fields are handled in
        //! DelimitedLineParser::parseLine().
        enum EmptyFieldHandling
        {
            //! Treat empty fields like any other field. Fields with
            //! empty names at the end of the line are still ignored.
            PROCESS_EMPTY_FIELDS = 0,
            //! Skip empty fields, such that they do not contribute to
            //! the field count for the line.
            SKIP_EMPTY_FIELDS = 1
        };

        //! Enumeration to describe how strictly we should check for
        //! badly formed input.
        enum StrictnessChecking
        {
            //! Strict syntax checking -- no fields at the end of the
            //! line, for example.
            STRICT_CHECKING = 0,
            //! Relaxed syntax checking -- allows empty header fields at the
            //! end of the line, for example.
            RELAXED_CHECKING = 1,
            //! Very relaxed syntax checking -- also allows empty header
            //! fields elsewhere.
            VERY_RELAXED_CHECKING = 2
        };

        typedef std::vector< boost::shared_ptr<DelimitedFieldParser> > Fields;

        //! Construct a DelimitedLineParser with no fields.
        //! @param allowOptionalFieldsInLine allow optional fields at the end of a line that can be missing
        DelimitedLineParser(bool allowOptionalFieldsInLine = false);

        //! Adds another DelimitedFieldParser to this
        //! DelimitedLineParser. The DelimitedFieldParsers must be added
        //! in order from left to right.
        template <class Field>
        DelimitedLineParser& addField(const Field& field)
        {
            boost::shared_ptr<Field> ptr(new Field(field));
            fields_.push_back(ptr);
            return *this;
        }

        //! Overrides the DelimitedFieldParser at a given field offset.
        template <class Field>
        DelimitedLineParser& setField(size_t offset, const Field& field)
        {
            CGA_ASSERT(offset < fields_.size());
            boost::shared_ptr<Field> ptr(new Field(field));
            fields_[offset] = ptr;
            return *this;
        }

        //! Parses the given line.
        void parseLine( const char* first, const char* last,
                        char delimiter = '\t',
                        EmptyFieldHandling emptyHandling = PROCESS_EMPTY_FIELDS,
                        StrictnessChecking strictnessChecking = RELAXED_CHECKING);

        void parseLine( const std::string& line,
                        char delimiter = '\t',
                        EmptyFieldHandling emptyHandling = PROCESS_EMPTY_FIELDS,
                        StrictnessChecking strictnessChecking = RELAXED_CHECKING);

        //! Provides read-only access to the list of currently defined fields
        const Fields& getFields() const {return fields_;}

    private:
        Fields  fields_;
        bool    allowOptionalFieldsInLine_;
    };

    class DelimitedField : public DelimitedFieldParser
    {
    public:
        DelimitedField(const std::string& name, char delimiter, const DelimitedLineParser& lineParser)
            : DelimitedFieldParser(name),
            delimiter_(delimiter),
            lineParser_(lineParser)
        {
        }

        void parse(const char* first, const char* last);

    private:
        char delimiter_;
        DelimitedLineParser lineParser_;
    };

} } // cgatools::util

#endif // CGATOOLS_UTIL_DELIMITEDLINEPARSER_HPP_
