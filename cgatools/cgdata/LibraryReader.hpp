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

#ifndef CGATOOLS_CGDATA_LIBRARY_READER_HPP_
#define CGATOOLS_CGDATA_LIBRARY_READER_HPP_ 1

//! @file LibraryReader.hpp

#include "cgatools/core.hpp"
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/array.hpp>

namespace cgatools { namespace cgdata {

//! The class contains simple gap estimation table template. 
template <class T>
class GapTable
{
public:
    class GapsRecord
    {
        friend class GapTable;
    public:
        bool operator== (const T& gaps) const {
            return gaps_ == gaps;
        }

        bool operator< (const GapsRecord& r) const {
            return gaps_ < r.gaps_;
        }

        GapsRecord(const T &gaps, double frequency_ = 0)
            : gaps_(gaps), frequency_(frequency_) {}

        T       gaps_;
        double  frequency_;

    protected:
        GapsRecord()
            : frequency_(0) {}

    };
    typedef std::vector<GapsRecord> Table;

    GapTable(const std::string& fname, double outliersFraction = 0) 
        :fileName_(fname),
        outliersFraction_(outliersFraction)
    {
        load(fileName_);

        if (outliersFraction_!=0)
            cutOutliers();
    }

    double getFrequency( const T& gaps, double outsideTheRange=0 )
    {
        GapsRecord record(gaps);
        typename Table::const_iterator it=std::lower_bound(records_.begin(),records_.end(),GapsRecord(gaps));

        if (it!=records_.end())
        {
            if (it->gaps_==gaps)
                return it->frequency_;
            if (records_.front()<it->gaps_)
                return 0;
        }
        return outsideTheRange;
    }


protected:
    void load(const std::string& fname);

    void cutOutliers()
    {
        double sideCutFraction = outliersFraction_/2;
        size_t  left=0,
                right=0;

        double sideFractionLeft = 0;
        for (left = 0; left < records_.size(); ++left)
        {
            sideFractionLeft += records_[left].frequency_;
            if (sideFractionLeft > sideCutFraction) {
                break;
            }
        }

        double sideFractionRight = 0;
        for (right = records_.size(); right > 0; )
        {
            --right;
            sideFractionRight += records_[right].frequency_;
            if (sideFractionRight > sideCutFraction) {
                ++right;
                break;
            }
        }
        records_.erase(records_.begin()+right,records_.end());
        records_.erase(records_.begin(),records_.begin()+left);
    }

    Table       records_;

    std::string fileName_;
    double      outliersFraction_;
};

//! The following types are used to parametrize simple gap estimation table
typedef boost::array<int8_t,3>  SmallGapTuple;
typedef int                     MateGapSize;

typedef GapTable<SmallGapTuple> WobbleGapTable;
typedef GapTable<MateGapSize> MateGapTable;

} } // cgatools::cgdata

#endif // CGATOOLS_CGDATA_LIBRARY_READER_HPP_
