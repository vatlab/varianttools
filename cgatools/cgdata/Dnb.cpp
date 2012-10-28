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
#include "cgatools/cgdata/Dnb.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Streams.hpp"

#include <boost/foreach.hpp>

#include <numeric>

namespace cgatools { namespace cgdata {

    std::ostream& operator <<(std::ostream& out, const LibRecord& r) {
        return  out << r.id_ << '\t' << r.type_ << '\t' << r.arm_ << '\t' 
            << r.indInArm_ << '\t' << r.objInArm_ 
            << '\t' << r.minSize_ << '\t' << r.maxSize_;
    }

    void LibRecord::initParser( util::DelimitedFile &delimitedFile )
    {
        CGA_ASSERT(delimitedFile.getDelimitedLineParser().getFields().size()==7);
        delimitedFile.addField(util::ValueField<size_t>("id",&id_));
        delimitedFile.addField(util::StringField("type",&type_));
        delimitedFile.addField(util::ValueField<size_t>("armId",&arm_));
        delimitedFile.addField(util::ValueField<size_t>("indArm",&indInArm_));
        delimitedFile.addField(util::ValueField<size_t>("objArm",&objInArm_));
        delimitedFile.addField(util::ValueField<int>("min",&minSize_));
        delimitedFile.addField(util::ValueField<int>("max",&maxSize_).exception("", 0));
        CGA_ASSERT(delimitedFile.getDelimitedLineParser().getFields().size()==7);
    }

    void DnbStructure::init(const std::string &libFile)
    {
        boost::shared_ptr<std::istream> libFileStream = 
            util::InputStream::openCompressedInputStreamByExtension(libFile);
        util::DelimitedFile parser(*libFileStream, libFile);
        init(parser);
    }

    void DnbStructure::init( util::DelimitedFile &libFile )
    {
        LibRecord record;
        record.initParser(libFile);
        while (libFile.next()) {
            if (record.type_=="mategap") {
                continue;
            }
            CGA_ASSERT_MSG(record.arm_<2, CGA_VOUT(record.arm_));
            while (halfDnbs_.size()<=record.arm_)
                halfDnbs_.push_back(HalfDnbStructure());
            HalfDnbStructure &h = halfDnbs_[record.arm_];
            if (record.type_=="read") {
                CGA_ASSERT_MSG(record.objInArm_<4, CGA_VOUT(record.objInArm_));
                while (h.reads_.size()<=record.objInArm_)
                    h.reads_.push_back(0);
                CGA_ASSERT_MSG(record.maxSize_==record.minSize_,
                    CGA_VOUT(record.maxSize_)<<CGA_VOUT(record.minSize_));
                h.reads_[record.objInArm_] = record.maxSize_;
            }
            if (record.type_=="gap") {
                CGA_ASSERT_MSG(record.objInArm_<4, CGA_VOUT(record.objInArm_));
                while (h.gaps_.size()<=record.objInArm_)
                    h.gaps_.push_back(HalfDnbStructure::Gap());
                h.gaps_[record.objInArm_].minSize_ = record.minSize_;
                h.gaps_[record.objInArm_].maxSize_ = record.maxSize_;
                CGA_ASSERT_MSG(record.minSize_<=record.maxSize_,
                    CGA_VOUT(record.minSize_)<<CGA_VOUT(record.maxSize_));
            }
        }
        BOOST_FOREACH(HalfDnbStructure &h, halfDnbs_) {
            h.totReadLength_ = std::accumulate(h.reads_.begin(),h.reads_.end(),0);
            CGA_ASSERT_EQ(h.gaps_.size()+1,h.reads_.size());
        }
    }

} } // cgatools::mapping
