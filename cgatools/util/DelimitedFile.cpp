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
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/Streams.hpp"

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <sstream>
#include "cgatools/util/StringSet.hpp"

namespace cgatools { namespace util {

    const std::string DelimitedFileMetadata::OUTPUT_FORMAT_VERSION = "2.2";
    std::string DelimitedFileMetadata::PIPELINE_VERSION = "";

    using std::string;

    namespace ba = boost::algorithm;

    bool DelimitedFileMetadata::hasKey(const std::string& key) const
    {
        for(size_t ii=0; ii<kv_.size(); ii++)
        {
            if (kv_[ii].first == key)
                return true;
        }
        return false;
    }

    const std::string& DelimitedFileMetadata::get(const std::string& key) const
    {
        for(size_t ii=0; ii<kv_.size(); ii++)
        {
            if (kv_[ii].first == key)
                return kv_[ii].second;
        }
        reportError("metadata key not found: "+key);
        CGA_ASSERT(false); // never reached
    }

    void DelimitedFileMetadata::set(const std::string& key, const std::string& value)
    {
        for(size_t ii=0; ii<kv_.size(); ii++)
        {
            if (kv_[ii].first == key)
            {
                kv_[ii].second = value;
                return;
            }
        }
        add(key, value);
    }

    void DelimitedFileMetadata::add(const std::string& key, const std::string& value)
    {
        kv_.push_back(make_pair(key, value));
    }

    void DelimitedFileMetadata::removeAll(const std::string& key)
    {
        for(size_t ii=0; ii<kv_.size(); )
        {
            if (kv_[ii].first == key)
            {
                kv_.erase(kv_.begin()+ii);
                continue;
            }
            ii++;
        }
    }

    int DelimitedFileMetadata::getFormatVersion() const
    {
        std::string version;
        if (hasKey("FORMAT_VERSION"))
            version = get("FORMAT_VERSION");
        else if (hasKey("VERSION"))
            version = get("VERSION");
        else
            reportError("file has no version header");
        size_t dot = version.find('.');
        if (dot == std::string::npos)
            reportError("file contains bad version header value: " + version);;
        const char* p = version.c_str();
        int major = util::parseValue<int>(p, p + dot);
        int minor = util::parseValue<int>(p + dot + 1, p + version.size());
        if (minor >= 1000 || minor < 0 || major < 0)
            reportError("file contains bad version header value: " + version);;
        return major * 1000 + minor;
    }

    std::string DelimitedFileMetadata::getSoftwareVersionString() const
    {
        if (hasKey("SRC_SOFTWARE_VERSION"))
            return get("SRC_SOFTWARE_VERSION");
        if (hasKey("SOFTWARE_VERSION"))
            return get("SOFTWARE_VERSION");
        if (hasKey("BUILD"))
            return get("BUILD");
        return "";
    }

    void DelimitedFileMetadata::initDefaults(const DelimitedFileMetadata& meta)
    {
        set("GENERATED_BY","cgatools");
        std::ostringstream sstr;
        sstr << boost::posix_time::microsec_clock::local_time();
        set("GENERATED_AT", sstr.str());
        if ("" == PIPELINE_VERSION)
        {
            // If this version of cgatools is running in the pipeline,
            // you need to set DelimitedFileMetadata::PIPELINE_VERSION
            // to the correct version string before calling this function.
            CGA_ASSERT(!CGA_TOOLS_IS_PIPELINE);
            set("CGATOOLS_VERSION", CGA_TOOLS_VERSION);
        }
        else
            set("SOFTWARE_VERSION", PIPELINE_VERSION);
        set("FORMAT_VERSION", OUTPUT_FORMAT_VERSION);

        typedef std::pair<std::string, std::string> KV;
        BOOST_FOREACH(const KV& kv, meta.getMap())
        {
            if ("GENERATED_BY" == kv.first)
                continue;
            if ("GENERATED_AT" == kv.first)
                continue;
            if ("CGATOOLS_VERSION" == kv.first)
                continue;
            if ("FORMAT_VERSION" == kv.first)
                continue;
            if ("" != PIPELINE_VERSION &&
                "SOFTWARE_VERSION" == kv.first)
                continue;

            add(kv.first, kv.second);
        }
    }

    void DelimitedFileMetadata::sort()
    {
        std::sort(kv_.begin(), kv_.end());
    }

    void DelimitedFileMetadata::transfer(const DelimitedFileMetadata& src,
                                         const std::string& keys,
                                         const std::string& prefix)
    {
        typedef std::pair<std::string, std::string> KV;

        util::StringSet xferKeys(keys, "", "");
        BOOST_FOREACH(const KV& kv, src.getMap())
        {
            if (xferKeys.contains(kv.first))
                set(prefix + kv.first, kv.second);
        }
    }

    void DelimitedFileMetadata::reportError(const std::string& error) const
    {
        if (fileName_.empty())
            throw Exception(error);
        else
            throw Exception(error + " (in file " + fileName_ + ")");
    }

    std::ostream& operator<< (std::ostream& out, const DelimitedFileMetadata& meta)
    {
        const std::vector< std::pair<std::string, std::string> >& mm = meta.getMap();
        if (mm.size() > 0)
        {
            for(size_t ii=0; ii<mm.size(); ii++)
                out << "#" << mm[ii].first << "\t" << mm[ii].second << "\n";
            out << "\n";
        }

        return out;
    }

    DelimitedFile::DelimitedFile(std::istream& in,
                                 const std::string& fileName,
                                 char delimiter,
                                 EmptyFieldHandling emptyHandling,
                                 StrictnessChecking strictnessChecking)
        : in_(in),
          fileName_(fileName),
          delimiter_(delimiter),
          emptyHandling_(emptyHandling),
          strictnessChecking_(strictnessChecking),
          lineNo_(0),
          activeLineSetId_(),
          withinActiveLineSet_(true)
    {
        readHeader();
    }

    void DelimitedFile::reportError(const std::string& error) const
    {
        std::ostringstream errs;
        errs << error << " (in file " << fileName_ << ", line " << lineNo_ << ")";
        throw Exception(errs.str());
    }

    void DelimitedFile::addAllFields(std::vector<std::string>& fields)
    {
        fields.resize(columnHeaders_.size());
        for(size_t ii=0; ii<columnHeaders_.size(); ii++)
            lp_.setField(ii, StringField(columnHeaders_[ii], &fields[ii]));
    }

    bool DelimitedFile::next()
    {
        for (;;)
        {
            bool ok = InputStream::getline(in_, line_);
            if (!ok)
                return false;
            ++lineNo_;

            if (DelimitedLineParser::STRICT_CHECKING == strictnessChecking_)
                break;

            if (!line_.empty())
            {
                if (activeLineSetId_.empty())
                    break;
                else if (line_ == "#ON " + activeLineSetId_)
                    withinActiveLineSet_ = true;
                else if (line_ == "#OFF " + activeLineSetId_)
                    withinActiveLineSet_ = false;
                else if (withinActiveLineSet_)
                    break;
            }
        }

        if (withinActiveLineSet_)
        {
            try
            {
                lp_.parseLine(line_, delimiter_, emptyHandling_, strictnessChecking_);
            }
            catch(std::exception& ee)
            {
                reportError(std::string(ee.what()) + ": " + line_);
            }
        }
        return true;
    }

    size_t DelimitedFile::getFieldOffset(const std::string& fieldName) const
    {
        size_t result = columnHeaders_.size();
        for(size_t ii=0; ii<columnHeaders_.size(); ii++)
        {
            if (columnHeadersEqual(columnHeaders_[ii], fieldName))
            {
                if (columnHeaders_.size() != result)
                    reportError("multiple fields with same name: "+fieldName);
                result = ii;
            }
        }
        if (columnHeaders_.size() == result)
            reportError("missing required field: "+fieldName);
        return result;
    }

    bool DelimitedFile::hasField(const std::string& fieldName) const
    {
        for(size_t ii=0; ii<columnHeaders_.size(); ii++)
        {
            if (columnHeadersEqual(fieldName, columnHeaders_[ii]))
                return true;
        }
        return false;
    }

    bool DelimitedFile::columnHeadersEqual(const std::string& h1, const std::string& h2) const
    {
        if (h1 == h2)
            return true;
        if (ba::to_lower_copy(h1) == boost::to_lower_copy(h2))
            return true;
        return false;
    }

    void DelimitedFile::readHeader()
    {
        metadata_.setFileName(fileName_);
        bool ok;
        for (ok = InputStream::getline(in_, line_), ba::trim_right(line_), ++lineNo_;
             ok && (0 == line_.size() || '#' == line_[0]);
             ok = InputStream::getline(in_, line_), ba::trim_right(line_), ++lineNo_)
        {
            if (0 == line_.size())
                continue;

            line_ = line_.substr(1);
            if (0 != line_.size())
            {
                size_t idxTab = line_.find('\t');
                if (string::npos == idxTab)
                {
                    metadata_.add(line_, "");
                }
                else
                {
                    const string k = line_.substr(0, idxTab);
                    const string v = line_.substr(idxTab+1);
                    metadata_.add(k, v);
                    if ("LINE_SET" == k)
                    {
                        activeLineSetId_ = v;
                        withinActiveLineSet_ = false;
                    }
                }
            }
        }
        if (!ok)
            reportError("header missing");

        if (DelimitedLineParser::STRICT_CHECKING == strictnessChecking_ && '>' != line_[0])
            reportError("header missing \">\" start character");

        if ('>' == line_[0])
            line_ = line_.substr(1);
        string dString(1, delimiter_);
        ba::split(columnHeaders_, line_, ba::is_any_of(dString.c_str()));
        for(size_t ii=0; ii<columnHeaders_.size(); ii++)
        {
            if (DelimitedLineParser::SKIP_EMPTY_FIELDS == emptyHandling_ &&
                columnHeaders_[ii] == "")
            {
                columnHeaders_.erase(columnHeaders_.begin()+ii);
                ii--;
                continue;
            }
            lp_.addField(IgnoreField(columnHeaders_[ii]));
        }

        if (DelimitedLineParser::STRICT_CHECKING == strictnessChecking_)
        {
            for(size_t ii=0; ii<lp_.fields_.size(); ii++)
            {
                if (lp_.fields_[ii]->getName().empty())
                    reportError("empty field names in header");
            }
        }

        // Ignore trailing fields with empty field names.
        while (!lp_.fields_.empty() && lp_.fields_.back()->getName().empty())
        {
            lp_.fields_.pop_back();
            columnHeaders_.pop_back();
        }

        // Fields with empty field names are errors, except for any
        // trailing fields with empty names, which are ignored.
        if (DelimitedLineParser::VERY_RELAXED_CHECKING != strictnessChecking_)
        {
            for (size_t i=0; i<lp_.getFields().size(); ++i)
            {
                if (lp_.getFields()[i]->getName().empty())
                {
                    std::ostringstream nonameMsg;
                    nonameMsg << "column header contains fields with empty field names: ";
                    for (size_t j=0; j<lp_.getFields().size(); ++j)
                        nonameMsg << "'" <<lp_.getFields()[j]->getName() << "' ";
                    reportError(nonameMsg.str());
                }
            }
        }
    }

} } // cgatools::util
