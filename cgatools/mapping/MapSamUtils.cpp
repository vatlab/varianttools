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
#include "cgatools/mapping/MapSamUtils.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/mapping/EvidenceSamUtil.hpp"

#include <boost/foreach.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/lexical_cast.hpp>

#include <sstream>
#include <numeric>

namespace cgatools { namespace mapping {

    std::ostream& operator <<(std::ostream& out, const ReadsRecord& r) {
        out << r.flags_.flags_ << '\t' << r.reads_ << '\t' << r.scores_ << '\t' << r.recordIndex_;
        return out;
    }

    std::ostream& operator <<(std::ostream& out, const MappingsRecord& r) {
        out << r.flags_.flags_
            << '\t' << r.chr_
            << '\t' << r.offsetInChr_
            << '\t' << r.gaps_[0]
            << '\t' << r.gaps_[1]
            << '\t' << r.gaps_[2]
            << '\t' << r.weightChar_
            << '\t' << r.bestMate_
            << '\t' << "p:" << r.isPrimary_
            << '\t' << r.recordIndex_ ;
        return out;
    }


    void ReadsRecord::initParser( util::DelimitedFile &delimitedFile )
    {
        delimitedFile.addField(util::ValueField<size_t>("flags",&flags_.flags_));
        delimitedFile.addField(util::StringField("reads",&reads_));
        delimitedFile.addField(util::StringField("scores",&scores_));
    }

    void MappingsRecord::initParser( util::DelimitedFile &delimitedFile )
    {
        delimitedFile.addField(util::ValueField<size_t>("flags",&flags_.flags_));
        delimitedFile.addField(util::StringField("chromosome",&chr_));
        delimitedFile.addField(util::ValueField<int>("offsetInChr",&offsetInChr_));
        delimitedFile.addField(util::ValueField<int>("gap1",&gaps_[0]));
        delimitedFile.addField(util::ValueField<int>("gap2",&gaps_[1]));
        delimitedFile.addField(util::ValueField<int>("gap3",&gaps_[2]));
        delimitedFile.addField(util::CharField("weight",(char *)&weightChar_));
        delimitedFile.addField(util::ValueField<int>("mateRec",&bestMate_));
    }

    std::string MappingsRecord::createCigar( const cgdata::HalfDnbStructure::Reads &readLengths ) const
    {
        std::stringstream cigar;
        for (size_t i=0; i<readLengths.size(); ++i)
        {
            if (i>0)
            {
                int gap = gaps_[i-1];
                if (gap>0)
                    cigar << gap << 'N';
                if (gap<0)
                    cigar << -gap << 'B';
            }
            cigar << readLengths[i] << 'M';
        }
        return cigar.str();
    }


    BaseMappingSamRecord::BaseMappingSamRecord( 
        const std::string& readName, 
        const MappingsRecord &mapping, 
        const cgdata::DnbStructure& dnbStructure, 
        const std::string& fullReadSequence, 
        const std::string& fullReadScores, 
        const reference::CrrFile& reference
        ) 
    :   
        SamRecord(
            readName, 
            mapping.recordIndex_>=0, 
            mapping.flags_.getStrand()==1, 
            mapping.isPrimary_,
            mapping.flags_.getSide(),
            mapping.chr_.empty()? 0 : reference.getChromosomeId(mapping.chr_), 
            mapping.offsetInChr_,
            mapping.createCigar(dnbStructure.halfDnbs_[mapping.flags_.getSide()].
                        getReadsForStrand(mapping.flags_.getStrand()==1)),
            mapping.weightChar_-33,
            mapping.bestMate_!=mapping.recordIndex_,
            fullReadSequence,
            fullReadScores,
            //fix if the number of hlfDNBs > 2
            UInt16Pair(mapping.flags_.getSide()==0 ? 0 : dnbStructure.halfDnbs_[0].totReadLength_,
                        dnbStructure.halfDnbs_[mapping.flags_.getSide()].totReadLength_)
        )
    {
        CGA_ASSERT_L(side_,2); //fix definition of sequenceStartAndLength
        typeId_ = BASE_MAPPING;
    }

} } // cgatools::mapping
