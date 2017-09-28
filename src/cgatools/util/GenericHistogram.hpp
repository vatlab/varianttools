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

#ifndef CGA_TOOLS_GENERIC_HISTOGRAM_HPP_
#define CGA_TOOLS_GENERIC_HISTOGRAM_HPP_

//! @file GenericHistogram.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/Exception.hpp"
#include <boost/cstdint.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace cgatools { namespace util {

    //the index and value must support stream output operators
    //merge function also requires support of the "Value& operator+=(const Value &rhs)" for the value type

    template <class Value=size_t, class Index=size_t>
    class GenericHistogram {
    public:
        typedef Value value_type;
        typedef Index index_type;
        typedef GenericHistogram<value_type,index_type> my_type;

        typedef std::pair<Index, Value> Bucket;
        typedef std::vector<Bucket>     Buckets;

        typedef typename Buckets::iterator          iterator;
        typedef typename Buckets::const_iterator    const_iterator;

        GenericHistogram() {
            initBuckets();
        }

        GenericHistogram(const Index& min, const Index& max, const Index& bucketSize) {
            initBuckets();
            addRange(min, max, bucketSize);
        }

        //0,7,2 will produce the buckets: 0,2,4,6
        void addRange(const Index& min, const Index& max, const Index& bucketSize) {
            CGA_ASSERT_MSG(bucketSize<bucketSize+bucketSize,"Invalid bucket size:"<<bucketSize);
            for(Index i=min; i<max; i+=bucketSize)
                addBucket(i);
        }

        iterator begin() {return buckets_.begin();}
        const_iterator begin() const {return buckets_.begin();}
        iterator end() {return buckets_.end();}
        const_iterator end() const {return buckets_.end();}
        iterator last() {return buckets_.end()-1;}
        const_iterator last() const {return buckets_.end()-1;}

        bool hasBucket(index_type index) {
            iterator b = findBucket(index);
            return b!=buckets_.end() && index==b->first;
        }

        void addBucket(index_type index) {
            iterator b = findBucket(index);
            if (b==buckets_.end() || !(index==b->first)
                //handle the case if the default bucket index equals the given index
                || (buckets_.size()==1 && index==index_type()) 
                )
                buckets_.insert(b,std::make_pair(index,value_type()));
        }

        value_type& operator[] (const index_type& index) {
            return findBucket(index)->second;
        }

        const value_type& operator[] (const index_type& index) const {
            return findBucket(index)->second;
        }

        iterator findBucket(Index index) {
            return std::lower_bound(buckets_.begin(), last(), index, compare_);
        }

        const_iterator findBucket(Index index) const {
            return std::lower_bound(buckets_.begin(), last(), index, compare_);
        }

        void merge(const my_type& src) {
            CGA_ASSERT_MSG(src.buckets_.size()==buckets_.size(),
                "Merge for different size histograms is not implemented. src:"
                <<src.buckets_.size() << " dst:"<<buckets_.size());
            for(int i=0, iEnd=src.buckets_.size(); i<iEnd; ++i) {
                const Bucket& srcObj = src.buckets_[i];
                Bucket& myObj = buckets_[i];
                CGA_ASSERT_MSG(srcObj.first == myObj.first,
                    "Different bucket ranges. src:"<<srcObj.first<<" dst:"<<myObj.first);
                myObj.second+=srcObj.second;
            }

        }

        //accumulate values in the buckets
        value_type totalInBuckets(value_type initialValue) const {
            for (const_iterator it=buckets_.begin(), itEnd=buckets_.end(); it!=itEnd; ++it)
                initialValue+=it->second;
            return initialValue;
        }

        void write( std::ostream& out, char separator = ',' ) const {
            out << ">bucket" << separator << "count" << separator << "cumulative count" << std::endl;
            if (buckets_.empty())
                return;
            Value cumulativeValue=begin()->second;
            for(typename GenericHistogram<Value,Index>::const_iterator it=begin(), 
                itEnd=last(); it!=itEnd; ++it, cumulativeValue+=it->second) 
            {
                out << it->first << separator << it->second << separator << cumulativeValue << std::endl;
            }
            out << "out of range" << separator << last()->second << separator << cumulativeValue <<std::endl;
        }

    protected:
        struct CompareBuckets {
            bool operator() (const Bucket& b0, const Bucket& b1) {return b0.first < b1.first;}
            bool operator() (const index_type& i, const Bucket& b) {return i < b.first;}
            bool operator() (const Bucket& b, const index_type& i) {return b.first < i;}
        };

        void initBuckets() {
            buckets_.push_back(std::make_pair(index_type(),value_type()));
        }

        Buckets buckets_;
        CompareBuckets compare_;
    };

    template <class Value, class Index>
    std::ostream& operator<<(std::ostream& out, const GenericHistogram<Value,Index>& r) {
        r.write(out);
        return out;
    }

    class SimpleHistogram
    {
        friend std::ostream& operator<<(std::ostream& out, const SimpleHistogram& src);
    public:
        SimpleHistogram(size_t maxValue)
            : count_(maxValue+1, 0), sum_(0), number_(0)
        {}

        void operator()(size_t val)
        {
            if (val+1 >= count_.size()) {
                ++( count_.back() );
            } else {
                ++( count_[val] );
            }
            ++number_;
            sum_ += val;
        }

        void operator+=(const SimpleHistogram& other)
        {
            CGA_ASSERT(other.count_.size() == count_.size());
            for (size_t ii = 0; ii < count_.size(); ++ii) {
                count_[ii] += other.count_[ii];
            }
            sum_ += other.sum_;
            number_ += other.number_;
        }

        void write(std::ostream& out) const;
    private:
        std::vector<size_t> count_;
        size_t sum_;
        size_t number_;
    };

} }
#endif //CGA_TOOLS_GENERIC_HISTOGRAM_HPP_
