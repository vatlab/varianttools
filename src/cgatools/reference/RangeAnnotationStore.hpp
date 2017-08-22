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

#ifndef CGATOOLS_REFERENCE_RANGEANNOTATIONSTORE_HPP_
#define CGATOOLS_REFERENCE_RANGEANNOTATIONSTORE_HPP_ 1

//! @file RangeAnnotationStore.hpp

#include "cgatools/core.hpp"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

#include "cgatools/util/Streams.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/RangeIntersector.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"

using cgatools::util::InputStream;

namespace cgatools { namespace reference {

//! Front end to IntervalTree based on reference ranges that allows
//! to load annotations from a file. To build your own annotation store,
//! derive from this class and provide an implementation for bindColumns
//! method to associate column names to fields of the TValue structure.
//! See class RepMaskStore in VarFileCombine.cpp for an example.
template <typename Derived, typename TValue>
class RangeAnnotationStore
{
private:
    typedef util::IntervalTree<reference::Range,
                               reference::Location,
                               TValue,
                               reference::RangeOverlap,
                               reference::GetRangeBoundary > DataStore;
public:
    typedef typename DataStore::QueryResultType QueryResultType;

    RangeAnnotationStore(const reference::CrrFile& crr)
        : crr_(crr)
    {}

    //! Returns metadata of the source delimited file.
    const util::DelimitedFileMetadata& getMetadata() const
    {
        return metadata_;
    }

    //! Fills the vector with iterators that point to the ranges intersecting
    //! a given range.
    void intersect(const reference::Range& range,
                   std::vector<QueryResultType>& result) const
    {
        data_.intersect(range, result);
    }

#if 0
    //! Derived class must implement this function to add column parsers for
    //! the range and association columns.
    void bindColumns(util::DelimitedFile& df, reference::Range& range, TValue& data)
#endif

protected:
    //! Typedef for this class to simplify constructor invocation
    //! in descendants
    typedef RangeAnnotationStore<Derived, TValue> Base;

    //! Load data from the file with a given name.
    void load(const std::string& fn, char delimiter = '\t')
    {
        boost::shared_ptr<std::istream> in =
                InputStream::openCompressedInputStreamByExtension(fn);
        util::DelimitedFile df(*in, fn, delimiter);
        metadata_ = df.getMetadata();

        reference::Range range;
        TValue payload;
        static_cast<Derived*>(this)->bindColumns(df, range, payload);

        while (df.next())
            data_.put(range, payload);
    }

    //! Utility function that the implementation of the derived class
    //! can call from bindColumns to add the parsers for the columns that
    //! correspond to range data.
    void bindRangeColumns(util::DelimitedFile& df,
                          reference::Range& range,
                          const std::string& chrColName = "chromosome",
                          const std::string& begColName = "begin",
                          const std::string& endColName = "end")
    {
        df.addField(ChromosomeIdField(chrColName, &range.chromosome_, crr_));
        df.addField(util::ValueField<uint32_t>(begColName, &range.begin_));
        df.addField(util::ValueField<uint32_t>(endColName, &range.end_));
    }

private:
    const reference::CrrFile& crr_;
    util::DelimitedFileMetadata metadata_;
    DataStore data_;
};

}}

#endif
