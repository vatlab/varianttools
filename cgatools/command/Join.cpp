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
#include "cgatools/command/Join.hpp"
#include "cgatools/util/RangeIntersector.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/IndirectComparator.hpp"

#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <math.h>

namespace cgatools { namespace command {

    using namespace cgatools::util;
    using namespace std;
    using boost::array;
    using boost::shared_ptr;
    namespace ba = boost::algorithm;

    template <class T>
    class PairOverlap
    {
    public:
        bool operator()(const std::pair<T,T>& lhs, const std::pair<T,T>& rhs) const
        {
            if (lhs.second < rhs.first)
                return false;
            if (rhs.second < lhs.first)
                return false;
            if (lhs.first != lhs.second && rhs.first != rhs.second)
            {
                if (lhs.second == rhs.first || rhs.second == lhs.first)
                    return false;
            }
            return true;
        }
    };

    template <class T>
    struct GetPairBoundary
    {
        T operator()(const std::pair<T,T>& r, size_t side) const
        {
            if (0 == side)
                return r.first;
            else
                return r.second;
        }
    };

    class JoinOutputStream
    {
    public:
        JoinOutputStream(std::ostream& out,
                         const std::vector< std::pair<size_t, size_t> >& transform,
                         bool alwaysDump)
            : out_(out),
              transform_(transform),
              hitCount_(0),
              alwaysDump_(alwaysDump),
              delimiter_('\t')
        {
            for(size_t ii=0; ii<transform.size(); ii++)
            {
                if (1 == transform[ii].first)
                {
                    bMask_.resize(transform[ii].second+1);
                    bMask_[transform[ii].second] = true;
                }
            }
            bFields_.resize(bMask_.size());
        }

        virtual ~JoinOutputStream()
        {
        }

        void writeHeaders(const DelimitedFile& aa, const DelimitedFile& bb)
        {
            DelimitedFile::Metadata combined(aa.getMetadata());
            const vector< pair<string,string> >& kv = bb.getMetadata().getMap();
            for(size_t ii=0; ii<kv.size(); ii++)
            {
                if (!combined.hasKey(kv[ii].first))
                    combined.set(kv[ii].first, kv[ii].second);
            }
            out_ << combined;
            out_ << ">";
            writeData(aa.getColumnHeaders(), bb.getColumnHeaders());
        }

        virtual void nextB(const std::vector<std::string>& aFields,
                           const std::vector<std::string>& bFields,
                           const std::pair<int64_t,int64_t>& rangeA,
                           const std::pair<int64_t,int64_t>& rangeB) = 0;
        virtual void finishRecord(const std::vector<std::string>& aFields,
                                  const std::pair<int64_t,int64_t>& rangeA) = 0;

    protected:
        void writeData(const std::vector<std::string>& aFields,
                       const std::vector<std::string>& bFields)
        {
            array< const vector<std::string>*, 2 > inFields;
            inFields[0] = &aFields;
            inFields[1] = &bFields;
            for(size_t ii=0; ii<transform_.size(); ii++)
            {
                if (ii > 0)
                    out_ << delimiter_;

                pair<size_t, size_t>& transform = transform_[ii];
                out_ << (*inFields[transform.first])[transform.second];
            }
            out_ << "\n";
        }

        std::ostream& out_;
        std::vector< std::pair<size_t, size_t> > transform_;

        // Fields of B for use by derived classes. Initially contains
        // enough empty values for use with writeData().
        std::vector<std::string> bFields_;

        // Fields used by the transform_.
        std::vector<bool> bMask_;

        size_t hitCount_; //!< # B records this A record
        bool alwaysDump_;
        char delimiter_;
    };

    class FullJoinOutputStream : public JoinOutputStream
    {
    public:
        FullJoinOutputStream(std::ostream& out,
                             const std::vector< std::pair<size_t, size_t> >& transform,
                             bool alwaysDump)
            : JoinOutputStream(out, transform, alwaysDump)
        {
        }

        void nextB(const std::vector<std::string>& aFields,
                   const std::vector<std::string>& bFields,
                   const std::pair<int64_t,int64_t>& rangeA,
                   const std::pair<int64_t,int64_t>& rangeB)
        {
            writeData(aFields, bFields);
            hitCount_++;
        }

        void finishRecord(const std::vector<std::string>& aFields,
                          const std::pair<int64_t,int64_t>& rangeA)
        {
            if (0 == hitCount_ && alwaysDump_)
                writeData(aFields, bFields_);
            hitCount_ = 0;
        }
    };

    class CompactJoinOutputStream : public JoinOutputStream
    {
    public:
        CompactJoinOutputStream(std::ostream& out,
                                const std::vector< std::pair<size_t, size_t> >& transform,
                                bool alwaysDump)
            : JoinOutputStream(out, transform, alwaysDump)
        {
            bAllFields_.resize(bFields_.size());
        }

        void nextB(const std::vector<std::string>& aFields,
                   const std::vector<std::string>& bFields,
                   const std::pair<int64_t,int64_t>& rangeA,
                   const std::pair<int64_t,int64_t>& rangeB)
        {
            bMask_.resize(bFields.size());
            bFields_.resize(bFields.size());
            bAllFields_.resize(bFields.size());

            for(size_t ii=0; ii<bFields.size(); ii++)
            {
                if (0 == bFields[ii].size())
                    continue;
                if (!bMask_[ii])
                    continue;

                bool skip = false;
                for(size_t jj=0; jj<bAllFields_[ii].size(); jj++)
                {
                    if (bFields[ii] == bAllFields_[ii][jj])
                    {
                        skip = true;
                        continue;
                    }
                }
                if (skip)
                    continue;

                bAllFields_[ii].push_back(bFields[ii]);
            }
            hitCount_++;
        }

        void finishRecord(const std::vector<std::string>& aFields,
                          const std::pair<int64_t,int64_t>& rangeA)
        {
            if (hitCount_ > 0 || alwaysDump_)
            {
                for(size_t ii=0; ii<bAllFields_.size(); ii++)
                {
                    vector<size_t> order(bAllFields_[ii].size());
                    for(size_t jj=0; jj<order.size(); jj++)
                        order[jj] = jj;
                    std::sort(order.begin(), order.end(),
                              IndirectComparator< vector<string> >(bAllFields_[ii]));

                    bFields_[ii].clear();
                    for(size_t jj=0; jj<bAllFields_[ii].size(); jj++)
                    {
                        if (jj > 0)
                            bFields_[ii].push_back(';');
                        bFields_[ii] += bAllFields_[ii][order[jj]];
                    }
                    bAllFields_[ii].clear();
                }
                writeData(aFields, bFields_);
            }
            hitCount_ = 0;
        }

    private:
        std::vector< std::vector<std::string> > bAllFields_;
    };

    class CompactPctJoinOutputStream : public JoinOutputStream
    {
        struct FieldRec
        {
            FieldRec(const std::string& name)
                : name_(name)
            {
            }

            bool operator<(const FieldRec& rhs) const
            {
                return name_ < rhs.name_;
            }

            std::string name_;
            std::vector< std::pair<int64_t,int64_t> > ranges_;
        };

    public:
        CompactPctJoinOutputStream(std::ostream& out,
                                   const std::vector< std::pair<size_t, size_t> >& transform,
                                   bool alwaysDump)
            : JoinOutputStream(out, transform, alwaysDump)
        {
        }

        void nextB(const std::vector<std::string>& aFields,
                   const std::vector<std::string>& bFields,
                   const std::pair<int64_t,int64_t>& rangeA,
                   const std::pair<int64_t,int64_t>& rangeB)
        {
            bMask_.resize(bFields.size());
            bFields_.resize(bFields.size());
            bAllFields_.resize(bFields.size());

            for(size_t ii=0; ii<bFields.size(); ii++)
            {
                if (0 == bFields[ii].size())
                    continue;
                if (!bMask_[ii])
                    continue;

                bool skip = false;
                for(size_t jj=0; jj<bAllFields_[ii].size(); jj++)
                {
                    if (bFields[ii] == bAllFields_[ii][jj].name_)
                    {
                        bAllFields_[ii][jj].ranges_.push_back(rangeB);
                        skip = true;
                        continue;
                    }
                }
                if (skip)
                    continue;

                bAllFields_[ii].push_back(FieldRec(bFields[ii]));
                bAllFields_[ii].back().ranges_.push_back(rangeB);
            }
            hitCount_++;
        }

        void finishRecord(const std::vector<std::string>& aFields,
                          const std::pair<int64_t,int64_t>& rangeA)
        {
            if (hitCount_ > 0 || alwaysDump_)
            {
                for(size_t ii=0; ii<bAllFields_.size(); ii++)
                {
                    vector<size_t> order(bAllFields_[ii].size());
                    for(size_t jj=0; jj<order.size(); jj++)
                        order[jj] = jj;
                    std::sort(order.begin(), order.end(),
                              IndirectComparator< vector< FieldRec> >(bAllFields_[ii]));

                    bFields_[ii].clear();
                    for(size_t jj=0; jj<bAllFields_[ii].size(); jj++)
                    {
                        if (jj > 0)
                            bFields_[ii].push_back(';');
                        bFields_[ii] += bAllFields_[ii][order[jj]].name_;
                        bFields_[ii].push_back(':');
                        bFields_[ii] += getPctOverlapString(rangeA, bAllFields_[ii][order[jj]].ranges_);
                    }
                    bAllFields_[ii].clear();
                }
                writeData(aFields, bFields_);
            }
            hitCount_ = 0;
        }

        std::string getPctOverlapString(const std::pair<int64_t,int64_t>& rangeA,
                                        std::vector< std::pair<int64_t,int64_t> >& ranges)
        {
            std::sort(ranges.begin(), ranges.end());
            int64_t overlap = 0;
            size_t iiNext=0;
            for(size_t ii=iiNext; ii<ranges.size(); ii=iiNext)
            {
                pair<int64_t, int64_t> rangeB = ranges[ii];
                for(iiNext=ii+1; iiNext<ranges.size(); iiNext++)
                {
                    if (rangeB.second < ranges[iiNext].first)
                        break;

                    rangeB.first = std::min(rangeB.first, ranges[iiNext].first);
                    rangeB.second = std::max(rangeB.second, ranges[iiNext].second);
                }

                CGA_ASSERT(rangeB.first <= rangeA.second);
                CGA_ASSERT(rangeA.first <= rangeB.second);
                overlap += std::min(rangeA.second, rangeB.second) - std::max(rangeA.first, rangeB.first);
            }

            if (rangeA.first == rangeA.second)
                return "100";

            double overlapFraction = double(overlap) / (rangeA.second-rangeA.first);
            int overlapPct = int(ceil(100.0 * overlapFraction));
            overlapPct = std::min(overlapPct, 100);
            overlapPct = std::max(overlapPct, 1);
            return boost::lexical_cast<string>(overlapPct);
        }

    private:
        std::vector< std::vector< FieldRec > > bAllFields_;
    };

    Join::Join(const std::string& name)
        : Command(name,
                  "Joins two tab-delimited files based on equal fields or overlapping regions.",
                  "Any",
                  "Joins two tab-delimited files based on equal fields or overlapping regions. "
                  "By default, an output record is produced for each match found between "
                  "file A and file B, but output format can be controlled by the "
                  "--output-mode parameter."
            )
    {
        options_.add_options()
            ("input", po::value< vector<string> >(&inputFileNames_),
             "File name to use as input (may be passed in as arguments at the end "
             "of the command), or omitted "
             "for stdin). There must be exactly two input files to join. If only "
             "one file is specified by name, file A is taken to be stdin and "
             "file B is the named file. File B is read fully into memory, and "
             "file A is streamed. File A's columns appear first in the output.")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file name (may be omitted for stdout).")
            ("match", po::value< vector<string> >(&match_),
             "A match specification, which is a column from A and a column from B "
             "separated by a colon.")
            ("overlap", po::value<string>(&overlapSpec_)->default_value(""),
             "Overlap specification. An overlap specification consists of a range "
             "definition for files A and B, separated by a colon. A range definition "
             "may be two columns, in which case they are interpreted as the beginning "
             "and end of the range. Or it may be one column, in which case the range "
             "is defined as the 1-base range starting at the given value. The records "
             "from the two files must overlap in order to be considered for output. "
             "Two ranges are considered to overlap if the overlap is at least one base "
             "long, or if one of the ranges is length 0 and the ranges overlap or abut. "
             "For example, \"begin,end:offset\" will match wherever end-begin &gt; 0, "
             "begin&lt;offset+1, and end&gt;offset, or wherever end-begin = 0, begin&lt;=offset+1, "
             "and end&gt;=offset.")
            ("output-mode,m", po::value<string>(&outputMode_)->default_value("full"),
             "Output mode, one of the following:\n"
             "    full        \tPrint an output record for each match found between "
             "file A and file B.\n"
             "    compact     \tPrint at most one record for each record of file A, "
             "joining the file B values by a semicolon and suppressing repeated B "
             "values and empty B values.\n"
             "    compact-pct \tSame as compact, but for each distinct B value, "
             "annotate with the percentage of the A record that is overlapped by "
             "B records with that B value. Percentage is rounded up to nearest integer.")
            ("overlap-mode", po::value<string>(&overlapMode_)->default_value("strict"),
             "Overlap mode, one of the following:\n"
             "    strict                \tRange A and B overlap if A.begin &lt; B.end "
             "and B.begin &lt; A.end.\n"
             "    allow-abutting-points \tRange A and B overlap they meet the strict requirements, or "
             "if A.begin &lt;= B.end and B.begin &lt;= A.end and either A or B has zero length.")
            ("select", po::value<string>(&select_)->default_value("A.*,B.*"),
             "Set of fields to select for output.")
            ("always-dump,a", po::bool_switch(&alwaysDump_)->default_value(false),
             "Dump every record of A, even if there are no matches with file B.")
            ("overlap-fraction-A", po::value<double>(&overlapFractionA_)->default_value(0.0),
             "Minimum fraction of A region overlap for filtering output.")
            ("boundary-uncertainty-A", po::value<int64_t>(&boundaryUncertaintyA_)->default_value(0),
             "Boundary uncertainty for overlap filtering. Specifically, records "
             "failing the following predicate are filtered away: "
             "overlap &gt;= overlap-fraction-A * ( A-range-length - boundary-uncertainty-A )")
            ("overlap-fraction-B", po::value<double>(&overlapFractionB_)->default_value(0.0),
             "Minimum fraction of B region overlap for filtering output.")
            ("boundary-uncertainty-B", po::value<int64_t>(&boundaryUncertaintyB_)->default_value(0),
             "Boundary uncertainty for overlap filtering. Specifically, records "
             "failing the following predicate are filtered away: "
             "overlap &gt;= overlap-fraction-B * ( B-range-length - boundary-uncertainty-B )")
            ;

        positionalOptions_.add("input", -1);
    }

    int Join::run(po::variables_map& vm)
    {
        requireParam(vm, "input");
        qp_.resize(2);

        if (inputFileNames_.size() > 2)
            throw Exception("too many file names: expected exactly two input files");

        std::string fnA;
        shared_ptr<istream> fB;
        if (1 == inputFileNames_.size())
            fnA = "STDIN";
        else
            fnA = inputFileNames_[0];
        fB = InputStream::openCompressedInputStreamByExtension(inputFileNames_.back());

        DelimitedFile aa(openStdin(fnA), fnA);
        DelimitedFile bb(*fB, inputFileNames_.back());

        // Check that the metadata for both files is as expected up front.
        aa.addAllFields(aFields_);
        bb.addAllFields(bFields_);
        initQueryPlan(aa, bb);

        std::ostream& of = openStdout(outputFileName_);
        shared_ptr<JoinOutputStream> out;
        if ("full" == outputMode_)
            out.reset(new FullJoinOutputStream(of, transform_, alwaysDump_));
        else if ("compact" == outputMode_)
            out.reset(new CompactJoinOutputStream(of, transform_, alwaysDump_));
        else if ("compact-pct" == outputMode_)
            out.reset(new CompactPctJoinOutputStream(of, transform_, alwaysDump_));
        else
            throw Exception("unrecognized output-mode: "+outputMode_);
        out->writeHeaders(aa, bb);

        // Read file B into memory.
        string matchKey;
        pair<int64_t, int64_t> range;
        vector< string > bLinesAll;
        typedef IntervalTree<pair<int64_t, int64_t>, int64_t, size_t,
                             PairOverlap<int64_t>,
                             GetPairBoundary<int64_t> > RI;
        map<string, RI> ri;
        while (bb.next())
        {
            parseJoinFields(bFields_, qp_[1], matchKey, range);
            if (range.second < range.first)
                throw Exception("range begins after it ends: "+bb.getLine());
            ri[matchKey].put(range, bLinesAll.size());
            bLinesAll.push_back(bb.getLine());
        }

        // Stream file A, filtering and annotating as necessary.
        std::vector<RI::QueryResultType> matches;
        while (aa.next())
        {
            parseJoinFields(aFields_, qp_[0], matchKey, range);
            if (range.second < range.first)
                throw Exception("range begins after it ends: "+aa.getLine());

            if (ri.find(matchKey) == ri.end())
            {
                out->finishRecord(aFields_, range);
                continue;
            }

            ri[matchKey].intersect(range, matches);
            vector< pair< pair<int64_t,int64_t>, size_t > > ordered;
            for(size_t ii=0; ii<matches.size(); ii++)
                ordered.push_back(make_pair(matches[ii]->first, matches[ii]->second));
            std::sort(ordered.begin(), ordered.end());

            for(size_t ii=0; ii<ordered.size(); ii++)
            {
                const pair<int64_t, int64_t>& mrange = ordered[ii].first;
                const string& bLine = bLinesAll[ordered[ii].second];
                vector<string> bFields;
                boost::split(bFields, bLine, boost::is_any_of("\t"));

                if ("strict" == overlapMode_)
                {
                    if ( ! (range.first < mrange.second && mrange.first < range.second) )
                        continue;
                }

                // Check overlap constraints.
                int64_t overlap = std::min(mrange.second, range.second) - std::max(mrange.first, range.first);
                CGA_ASSERT(overlap >= 0);
                if (overlapFractionA_ > 0.0)
                {
                    int64_t len = std::max(int64_t(0), range.second-range.first - boundaryUncertaintyA_);
                    if (overlap < len * overlapFractionA_)
                        continue;
                }
                if (overlapFractionB_ > 0.0)
                {
                    int64_t len = std::max(int64_t(0), mrange.second-mrange.first - boundaryUncertaintyB_);
                    if (overlap < len * overlapFractionB_)
                        continue;
                }

                // Add to bFields.
                out->nextB(aFields_, bFields, range, mrange);
            }

            out->finishRecord(aFields_, range);
        }

        return 0;
    }

    void Join::dumpRecord(std::ostream& out,
                          const std::vector<std::string>& aFields,
                          const std::vector<std::string>& bFields)
    {
        for(size_t ii=0; ii<aFields.size(); ii++)
        {
            if (ii > 0)
                out << "\t";
            out << aFields[ii];
        }

        for(size_t ii=0; ii<bFields.size(); ii++)
        {
            if (ii > 0 || aFields.size() > 0)
                out << "\t";
            out << bFields[ii];
        }

        out << endl;
    }

    bool Join::overlap(const std::pair<int64_t, int64_t>& lhs,
                       const std::pair<int64_t, int64_t>& rhs)
    {
        if (lhs.second < rhs.first || rhs.second < lhs.first)
            return false;
        if (lhs.first < lhs.second && rhs.first < rhs.second)
        {
            if (lhs.second == rhs.first || rhs.second == lhs.second)
                return false;
        }
        return true;
    }                       

    void Join::parseJoinFields(const std::vector<std::string>& fields,
                               const QueryPlan& qp,
                               std::string& matchKey,
                               std::pair<int64_t, int64_t>& range)
    {
        matchKey.clear();
        for(size_t ii=0; ii<qp.matchIdx_.size(); ii++)
        {
            if (ii > 0)
                matchKey.push_back(' ');
            matchKey += fields[qp.matchIdx_[ii]];
        }

        if (qp.overlapIdx_.first < 0)
            range = pair<int64_t,int64_t>(0,1);
        else
        {
            range.first = parseValue<int64_t>(fields[qp.overlapIdx_.first]);
            if (qp.overlapIdx_.second < 0)
                range.second = range.first+1;
            else
                range.second = parseValue<int64_t>(fields[qp.overlapIdx_.second]);
        }

        if (range.second < range.first)
        {
            string msg;
            for(size_t ii=0; ii<fields.size(); ii++)
            {
                if (ii > 0)
                    msg.append("\t");
                msg += fields[ii];
            }
            throw Exception("range end is before range begin: "+msg);
        }
    }

    void Join::initQueryPlan(util::DelimitedFile& aa, util::DelimitedFile& bb)
    {
        BOOST_FOREACH(const string& match, match_)
        {
            try
            {
                vector<string> parts;
                ba::split(parts, match, boost::is_any_of(":"));
                if (parts.size() != 2)
                    throw Exception("expected two colon-separated values");
                qp_[0].matchIdx_.push_back(aa.getFieldOffset(parts[0]));
                qp_[1].matchIdx_.push_back(bb.getFieldOffset(parts[1]));
                if (qp_[0].matchIdx_.size() != qp_[1].matchIdx_.size())
                    throw Exception("match field count mismatch");
            }
            catch(std::exception& ee)
            {
                throw Exception("failed parsing match spec \""+match+"\": "+ee.what());
            }
        }

        qp_[0].overlapIdx_ = make_pair(-1, -1);
        qp_[1].overlapIdx_ = make_pair(-1, -1);
        if (overlapSpec_ != "")
        {
            try
            {
                vector<string> parts;
                ba::split(parts, overlapSpec_, boost::is_any_of(":"));
                if (parts.size() != 2)
                    throw Exception("expected two colon-separated values");
                parseOverlapFields(aa, parts[0], qp_[0]);
                parseOverlapFields(bb, parts[1], qp_[1]);
            }
            catch(std::exception& ee)
            {
                throw Exception("failed parsing overlap spec \""+overlapSpec_+"\": "+ee.what());
            }
        }

        // Initialize transform_.
        vector<string> parts;
        ba::split(parts, select_, boost::is_any_of(","));
        BOOST_FOREACH(const string& part, parts)
        {
            if (boost::to_lower_copy(part) == "a.*")
            {
                for(size_t ii=0; ii<aFields_.size(); ii++)
                    transform_.push_back(pair<size_t,size_t>(0, ii));
                continue;
            }
            if (boost::to_lower_copy(part) == "b.*")
            {
                for(size_t ii=0; ii<bFields_.size(); ii++)
                    transform_.push_back(pair<size_t,size_t>(1, ii));
                continue;
            }
            if (part.size() > 2 && boost::to_lower_copy(part.substr(0, 2)) == "a.")
            {
                string name = part.substr(2);
                transform_.push_back(pair<size_t, size_t>(0, aa.getFieldOffset(name)));
                continue;
            }
            if (part.size() > 2 && boost::to_lower_copy(part.substr(0, 2)) == "b.")
            {
                string name = part.substr(2);
                transform_.push_back(pair<size_t, size_t>(1, bb.getFieldOffset(name)));
                continue;
            }
            size_t aId=aFields_.size(), bId=bFields_.size();
            try
            {
                aId = aa.getFieldOffset(part);
            }
            catch(const std::exception&)
            {
            }
            try
            {
                bId = bb.getFieldOffset(part);
            }
            catch(const std::exception&)
            {
            }
            if (aFields_.size() == aId && bFields_.size() == bId)
                throw Exception("field not found: "+part);
            if (aFields_.size() != aId && bFields_.size() != bId)
                throw Exception("ambiguous field selection: "+part);
            if (aFields_.size() == aId)
                transform_.push_back(pair<size_t, size_t>(1, bId));
            else
                transform_.push_back(pair<size_t, size_t>(0, aId));
        }
    }

    void Join::parseMatchFields(const util::DelimitedFile& df,
                                const std::string& fieldNameList,
                                QueryPlan& qp)
    {
        vector<string> fieldNames;
        ba::split(fieldNames, fieldNameList, boost::is_any_of(","));
        BOOST_FOREACH(const string& fieldName, fieldNames)
        {
            qp.matchIdx_.push_back(df.getFieldOffset(fieldName));
        }
    }

    void Join::parseOverlapFields(const util::DelimitedFile& df,
                                  const std::string& fieldNameList,
                                  QueryPlan& qp)
    {
        vector<string> fieldNames;
        ba::split(fieldNames, fieldNameList, boost::is_any_of(","));
        if (fieldNames.size() > 2)
            throw Exception("expected one or two field names per file");
        qp.overlapIdx_.first = df.getFieldOffset(fieldNames[0]);
        if (fieldNames.size() > 1)
            qp.overlapIdx_.second = df.getFieldOffset(fieldNames[1]);
        else
            qp.overlapIdx_.second = -1;
    }

} } // cgatools::command
