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
#include "cgatools/command/ListVariants.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"
#include "cgatools/variants/SuperlocusIterator.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"

#include <queue>
#include <set>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

using namespace std;

using boost::shared_ptr;

namespace bu = cgatools::util::baseutil;

using namespace cgatools::util;
using namespace cgatools::reference;
using namespace cgatools::variants;

namespace cgatools { namespace command {

    class VariantLister
    {
    public:
        virtual ~VariantLister()
        {
        }

        virtual bool eof() const = 0;

        virtual VariantLister& operator++() = 0;

        virtual const Call& operator*() const = 0;

        virtual const Call* operator->() const
        {
            return & ( * (*this) );
        }

    };

    class VariantFileVariantLister : public VariantLister
    {
    public:
        VariantFileVariantLister(const CrrFile& crr,
                                 const std::string& variantFileName)
            : locIt_(crr),
              idx_(0)
        {
            locIt_.open(variantFileName);
            buildQueue();
        }

        bool eof() const
        {
            return idx_ == queue_.size();
        }

        VariantLister& operator++()
        {
            idx_++;
            if (idx_ < queue_.size())
                return *this;

            idx_ = 0;
            queue_.clear();
            if (!locIt_.eof())
            {
                ++locIt_;
                buildQueue();
            }
            return *this;
        }

        const Call& operator*() const
        {
            return queue_[idx_];
        }

    private:
        void buildQueue();

        VariantFileIterator locIt_;
        std::vector<Call> queue_;
        size_t idx_;
    };

    void VariantFileVariantLister::buildQueue()
    {
        const CrrFile& crr = locIt_->getReference();
        for(; !locIt_.eof(); ++locIt_)
        {
            const vector<Call>& calls = locIt_->getCalls();
            for(size_t ii=0; ii<calls.size(); ii++)
            {
                if (bu::isCalledSequence(calls[ii].alleleSeq_) &&
                    calls[ii].refSequence(crr) != calls[ii].calledSequence(crr))
                {
                    queue_.push_back(calls[ii]);
                }
            }

            if (0 != queue_.size())
                break;
        }

        std::sort(queue_.begin(), queue_.end(), ListVariants::CallOrder());
    }

    class SuperlocusVariantLister : public VariantLister
    {
    public:
        SuperlocusVariantLister(const CrrFile& crr,
                                const std::string& variantFileName)
            : locIt_(crr),
              slIt_(4, 20),
              idx_(0)
        {
            locIt_.open(variantFileName);
            slIt_.setVariantFile(locIt_);
            slIt_.seekFirst();
            buildQueue();
        }

        bool eof() const
        {
            return idx_ == queue_.size();
        }

        VariantLister& operator++()
        {
            idx_++;
            if (idx_ < queue_.size())
                return *this;

            idx_ = 0;
            queue_.clear();
            if (!slIt_.eof())
            {
                ++slIt_;
                buildQueue();
            }
            return *this;
        }

        const Call& operator*() const
        {
            return queue_[idx_];
        }

    private:
        void buildQueue();
        Call alleleToCall(Range range, string allele);

        VariantFileIterator locIt_;
        SuperlocusIterator slIt_;
        std::vector<Call> queue_;
        size_t idx_;
    };

    void SuperlocusVariantLister::buildQueue()
    {
        const CrrFile& crr = locIt_->getReference();
        for(; !slIt_.eof(); ++slIt_)
        {
            const Superlocus& sl = *slIt_;
            Range range = sl.getRange();

            vector< vector<PhasedHypothesis> > hypotheses;
            sl.buildPhasedHypotheses(hypotheses, 4, true);

            if (!hypotheses[0].empty())
            {
                BOOST_FOREACH(const PhasedHypothesis& hypothesis, hypotheses[0])
                {
                    for(size_t ii=0; ii<hypothesis.size(); ii++)
                    {
                        const string& allele = hypothesis[ii].allele();
                        if (!bu::isCalledSequence(allele))
                            continue;

                        Call call = alleleToCall(range, allele);
                        if (call.reference_.size() > 0 || call.alleleSeq_.size() > 0)
                        {
                            queue_.push_back(call);
                        }
                    }
                }
            }

            BOOST_FOREACH(const Locus& locus, sl.getLoci(0))
            {
                BOOST_FOREACH(const Call& call, locus.getCalls())
                {
                    if (bu::isCalledSequence(call.alleleSeq_) &&
                        call.refSequence(crr) != call.calledSequence(crr))
                    {
                        queue_.push_back(call);
                    }
                }
            }

            if (0 != queue_.size())
                break;
        }

        std::sort(queue_.begin(), queue_.end(), ListVariants::CallOrder());
    }

    Call SuperlocusVariantLister::alleleToCall(Range range, string allele)
    {
        const CrrFile& crr = locIt_->getReference();
        string refAllele = crr.getSequence(range);

        // Reduce allele by prefix matching.
        for(size_t ii=0; ; ii++)
        {
            if (ii >= allele.size() || ii >= refAllele.size() || allele[ii] != refAllele[ii])
            {
                if (ii > 0)
                {
                    allele = allele.substr(ii);
                    refAllele = refAllele.substr(ii);
                    range.begin_ += ii;
                }
                break;
            }
        }

        // Reduce allele by suffix matching.
        while (allele.size() > 0 && refAllele.size() > 0 &&
               allele[allele.size()-1] == refAllele[refAllele.size()-1])
        {
            allele.resize(allele.size()-1);
            refAllele.resize(refAllele.size()-1);
            range.end_--;
        }

        Call call;
        call.range_ = range;
        call.reference_ = refAllele;
        call.alleleSeq_ = allele;
        if (call.reference_.size() == call.alleleSeq_.size())
        {
            if (1 == call.reference_.size())
                call.varType_ = "snp";
            else
                call.varType_ = "sub";
        }
        else if (0 == call.reference_.size())
            call.varType_ = "ins";
        else if (0 == call.alleleSeq_.size())
            call.varType_ = "del";
        else
            call.varType_ = "sub";
        return call;
    }

    class VariantListingFileVariantLister : public VariantLister
    {
    public:
        VariantListingFileVariantLister(const CrrFile& crr,
                                        const std::string& variantListingFileName)
            : crr_(crr),
              eofCall_(false),
              idx_(0)
        {
            in_ = InputStream::openCompressedInputStreamByExtension(variantListingFileName);
            df_.reset(new DelimitedFile(*in_, variantListingFileName));
            df_->addField(ChromosomeIdField("chromosome", &nextCall_.range_.chromosome_, crr_));
            df_->addField(ValueField<uint32_t>("begin", &nextCall_.range_.begin_));
            df_->addField(ValueField<uint32_t>("end", &nextCall_.range_.end_));
            df_->addField(StringField("varType", &nextCall_.varType_), DelimitedFile::FPT_OPTIONAL);
            df_->addField(StringField("reference", &nextCall_.reference_));
            df_->addField(StringField("alleleSeq", &nextCall_.alleleSeq_));
            df_->addField(StringField("xRef", &nextCall_.xRef_), DelimitedFile::FPT_OPTIONAL);
            readNext();
            buildQueue();
        }

        bool eof() const
        {
            return idx_ == queue_.size();
        }

        VariantLister& operator++()
        {
            idx_++;
            if (idx_ < queue_.size())
                return *this;

            idx_ = 0;
            queue_.clear();
            if (!eofCall_)
            {
                buildQueue();
            }
            return *this;
        }

        const Call& operator*() const
        {
            return queue_[idx_];
        }

    private:
        void buildQueue();
        void readNext();

        const CrrFile& crr_;
        boost::shared_ptr<std::istream> in_;
        boost::shared_ptr<DelimitedFile> df_;
        Call nextCall_;
        bool eofCall_;
        std::vector<Call> queue_;
        size_t idx_;
    };

    void VariantListingFileVariantLister::buildQueue()
    {
        for(; !eofCall_; )
        {
            CGA_ASSERT(bu::isCalledSequence(nextCall_.alleleSeq_) &&
                       nextCall_.refSequence(crr_) != nextCall_.calledSequence(crr_));
            queue_.push_back(nextCall_);

            readNext();
            if (0 != queue_.size())
                break;
        }

        std::sort(queue_.begin(), queue_.end(), ListVariants::CallOrder());
    }

    void VariantListingFileVariantLister::readNext()
    {
        if (!df_->next())
        {
            eofCall_ = true;
            return;
        }

        string nextReference = crr_.getSequence(nextCall_.range_);
        if (0 == nextCall_.range_.length() || nextCall_.reference_.size() > 0)
        {
            if (nextReference != nextCall_.reference_)
                throw Exception( (boost::format("variant listing file reference sequence mismatch at %s,%d") %
                                  crr_.listChromosomes()[nextCall_.range_.chromosome_].getName() %
                                  nextCall_.range_.begin_).str() );
        }
        nextCall_.reference_ = nextReference;
    }

    class VariantListerDescByPosition
    {
    public:
        bool operator()(const boost::shared_ptr<VariantLister>& lhs,
                        const boost::shared_ptr<VariantLister>& rhs) const
        {
            const Call& lc = *(*lhs);
            const Call& rc = *(*rhs);
            return ListVariants::CallOrder()(rc, lc);
        }
    };

    ListVariants::ListVariants(const std::string& name)
        : Command(name,
                  "Lists the variants present in a variant file.",
                  "0.3 or later",

                  "Lists all called variants present in the specified variant files, in "
                  "a format suitable for processing by the testvariants command. The output is "
                  "a tab-delimited file consisting of the following columns:\n\n"
                  "    variantId  \tSequential id assigned to each variant.\n"
                  "    chromosome \tThe chromosome of the variant.\n"
                  "    begin      \t0-based reference offset of the beginning of the variant.\n"
                  "    end        \t0-based reference offset of the end of the variant.\n"
                  "    varType    \tThe varType as extracted from the variant file.\n"
                  "    reference  \tThe reference sequence.\n"
                  "    alleleSeq  \tThe variant allele sequence as extracted from the variant file.\n"
                  "    xRef       \tThe xRef as extrated from the variant file."
            ),
          queueBaseCount_(1000),
          variantId_(1)
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file.")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ("variants", po::value< vector<string> >(&variantFileNames_),
             "The input variant files (may be passed in as argument at the end of the command).")
            ("variant-listing", po::value< vector<string> >(&variantListingFileNames_),
             "The output of another listvariants run, to be merged in to produce the output of this run.")
            ("list-long-variants", po::bool_switch(&listLongVariants_)->default_value(false),
             "In addition to listing short variants, list longer variants as well "
             "(10's of bases) by concatenating nearby calls.\n")
            ;

        positionalOptions_.add("variants", -1);
    }

    int ListVariants::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");

        crr_.open(referenceFileName_);
        std::ostream& out = openStdout(outputFileName_);

        priority_queue< shared_ptr<VariantLister>,
            vector< shared_ptr<VariantLister> >, VariantListerDescByPosition > pq;
        BOOST_FOREACH(const string& variantFileName, variantFileNames_)
        {
            if (!listLongVariants_)
            {
                shared_ptr<VariantLister> vl(new VariantFileVariantLister(crr_, variantFileName));
                if (!vl->eof())
                    pq.push(vl);
            }
            else
            {
                shared_ptr<VariantLister> vl(new SuperlocusVariantLister(crr_, variantFileName));
                if (!vl->eof())
                    pq.push(vl);
            }
        }
        BOOST_FOREACH(const string& variantListingFileName, variantListingFileNames_)
        {
            shared_ptr<VariantLister> vl(new VariantListingFileVariantLister(crr_, variantListingFileName));
            if (!vl->eof())
                pq.push(vl);
        }

        out << "variantId\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\txRef" << endl;
        set<Call, CallOrder> callSet;
        while (!pq.empty())
        {
            shared_ptr<VariantLister> vl = pq.top();
            pq.pop();

            queueCall(callSet, *(*vl));
            retireQueuedCalls(out, callSet);

            ++(*vl);
            if (!vl->eof())
                pq.push(vl);
        }

        purgeQueuedCalls(out, callSet);

        return 0;
    }

    void ListVariants::queueCall(std::set<Call, CallOrder>& callSet, Call call)
    {
        canonicalizeCall(call);
        callSet.insert(call);
    }

    void ListVariants::retireQueuedCalls(std::ostream& out, std::set<Call, CallOrder>& callSet)
    {
        if (callSet.size() > 0)
        {
            Location highWater = callSet.rbegin()->range_.beginLocation();
            if (highWater.offset_ > queueBaseCount_)
                highWater.offset_ -= queueBaseCount_;
            else
                highWater.offset_ = 0;
            while (!callSet.empty() && callSet.begin()->range_.endLocation() < highWater)
            {
                const Call& call = *callSet.begin();
                retireCall(out, call);
                callSet.erase(callSet.begin());
            }
        }
    }

    void ListVariants::purgeQueuedCalls(std::ostream& out, std::set<Call, CallOrder>& callSet)
    {
        while (!callSet.empty())
        {
            const Call& call = *callSet.begin();
            retireCall(out, call);
            callSet.erase(callSet.begin());
        }
    }

    void ListVariants::retireCall(std::ostream& out, const variants::Call& call)
    {
        if (!CallOrder()(prevCall_, call))
        {
            cerr << "variant dropped because canonical position shifted too far: ";
            call.write(cerr, crr_);
            return;
        }

        out << variantId_ << "\t"
            << crr_.listChromosomes()[call.range_.chromosome_].getName() << "\t"
            << call.range_.begin_ << "\t"
            << call.range_.end_ << "\t"
            << call.varType_ << "\t"
            << call.refSequence(crr_) << "\t"
            << call.calledSequence(crr_) << "\t"
            << call.xRef_ << "\n";
        prevCall_ = call;
        variantId_++;
    }

    void ListVariants::canonicalizeCall(Call& call) const
    {
        const CompactDnaSequence& chrom = crr_.listChromosomes()[call.range_.chromosome_];

        // Remove bases that match reference from the left and right ends
        uint32_t origBegin = call.range_.begin_;
        uint32_t begin = call.range_.begin_;
        uint32_t end = call.range_.end_;
        size_t alleleLen = call.alleleSeq_.size();
        while (alleleLen > 0 && begin < end &&
               chrom.getBase(begin) == call.alleleSeq_[begin - call.range_.begin_])
        {
            ++begin;
            --alleleLen;
        }
        while (alleleLen > 0 && begin < end &&
               chrom.getBase(end - 1) ==
               call.alleleSeq_[alleleLen + (begin - call.range_.begin_) - 1])
        {
            --end;
            --alleleLen;
        }
        if (alleleLen != call.alleleSeq_.size())
        {
            call.alleleSeq_ = call.alleleSeq_.substr(begin - call.range_.begin_, alleleLen);
            call.range_.begin_ = begin;
            call.range_.end_ = end;
            call.reference_ = crr_.getSequence(call.range_);
            if ( 0 == call.range_.length() && 0 == call.alleleSeq_.size() )
            {
                throw Exception( (boost::format("call sequence equals reference for call at %s,%d") %
                                  chrom.getName() % origBegin).str() );
            }
        }

        if (0 == call.range_.length())
            call.varType_ = "ins";
        else if (0 == call.alleleSeq_.size())
            call.varType_ = "del";
        else if (1 == call.range_.length() && 1 == call.alleleSeq_.size())
            call.varType_ = "snp";
        else
            call.varType_ = "sub";

        if (0 == call.range_.length() && 0 != call.alleleSeq_.size())
        {
            // Circular prefix match to the left.
            uint32_t len = call.alleleSeq_.size();
            uint32_t offset = call.range_.begin_;
            size_t sOffset = len-1;
            for(; offset>0; offset--)
            {
                if (chrom.getBase(offset-1) != call.alleleSeq_[sOffset])
                    break;
                if (0 == sOffset)
                    sOffset = len;
                sOffset--;
            }
            call.alleleSeq_ =
                rotateLeft(call.alleleSeq_, int32_t(offset)-int32_t(call.range_.begin_));
            call.range_.begin_ = call.range_.end_ = offset;
        }
        else if (0 == call.alleleSeq_.size() && 0 != call.range_.length())
        {
            // Circular prefix match to the left.
            uint32_t len = call.range_.length();
            uint32_t offset = call.range_.begin_;
            for(; offset>0; offset--)
            {
                if (chrom.getBase(offset-1) != chrom.getBase(offset+len-1))
                    break;
            }
            call.range_.begin_ = offset;
            call.range_.end_ = offset+len;
            call.reference_ = crr_.getSequence(call.range_);
        }
    }

    std::string ListVariants::rotateLeft(const std::string& str, int32_t count) const
    {
        if (count < 0)
        {
            int32_t rightShift = (-count) % int32_t(str.size());
            count = (str.size() - rightShift) % str.size();
        }
        string result(str.size(), '\0');
        for(size_t ii=0; ii<str.size(); ii++)
            result[ii] = str[(ii+count)%str.size()];
        return result;
    }

} } // cgatools::command
