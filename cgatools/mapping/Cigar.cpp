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
#include "cgatools/util/Exception.hpp"
#include "cgatools/mapping/Cigar.hpp"
#include "cgatools/util/parse.hpp"

#include <boost/foreach.hpp>
#include <cmath>
#include <locale>

namespace cgatools { namespace mapping {

    int Cigar::CigarElement::getSequenceLength() const
    {
        switch (type_) {
            case 'N':
            case 'M':
            case 'I':
            case 'S': // skip (in gap tag)
            case 'G': // gap (in gap tag)
                return length_;
            case 'B': // back (in extended cigar)
            case 'D': // delete
            case 'P': // padded sequence
                return 0;
            default:
                return 0;
        }
    }

    void Cigar::parse(const std::string &cigar)
    {
        CGA_ASSERT(!cigar.empty());
        CGA_ASSERT(cigar[0]>='0' && cigar[0]<='9');
        CGA_ASSERT(cigar[cigar.size()-1]>'9');
        const char* cigarArray = cigar.c_str();
        size_t numStart = 0;
        for (size_t i=0; i<cigar.size(); ++i) {
            CGA_ASSERT(cigar[i]>='0');
            if (cigar[i]>'9') {
                char cmdType = std::toupper(cigar[i], std::locale());
                size_t length = util::parseValue<size_t>(cigarArray+numStart,cigarArray+i);
                add(CigarElement(cmdType,length));
                numStart = i+1;
            }
        }
        CGA_ASSERT_EQ(numStart,cigar.size());
    }

    size_t Cigar::getSequenceLength() const
    {
        int result = 0;
        BOOST_FOREACH(const CigarElement& e, getParsedCigar()) {
            result+=e.getSequenceLength();
            CGA_ASSERT_MSG(result>=0,*this);
        }
        return (size_t)result;
    }

    size_t Cigar::getReferenceLength() const
    {
        ParsedCigar cigar(getParsedCigar());
        size_t result = 0;
        for (ParsedCigar::iterator it=cigar.begin(), itEnd=cigar.end(); it!=itEnd; ++it)
        {
            CigarElement& e = *it;
            switch (e.type_) {
            case 'N':
            case 'M':
            case 'D':
            case 'G': // gap (in gap tag)
            case 'S': // skip (in gap tag)
                e.type_='+';
                break;
            case 'B': // back (in extended cigar)
                for (ParsedCigar::iterator it1=it-1; e.length_>0; --it1)
                {
                    int minLen = std::min(it1->length_,e.length_);
                    e.length_-=minLen;
                    it1->length_-=minLen;
                    if (it1==cigar.begin())
                    {
                        CGA_ASSERT_MSG(e.length_==0,CGA_VOUT(*this)<<CGA_VOUT(e.length_));
                        break;
                    }
                }
                break;
            case 'I': // insert
            case 'P': // padded sequence
                e.type_='0';
                break;
            }
        }
        BOOST_FOREACH(const CigarElement& e, cigar)
            if (e.type_=='+')
                result+=e.length_;
        return result;
    }

    void Cigar::merge_back( const CigarElement &e )
    {
        if ((!parsedCigar_.empty()) && parsedCigar_.back().type_==e.type_)
            parsedCigar_.back().length_+=e.length_;
        else
            parsedCigar_.push_back(e);
    }

    void Cigar::trancatePaddings()
    {
        while (!parsedCigar_.empty()) {
            if (parsedCigar_.back().type_ == 'P')
                parsedCigar_.pop_back();
            else
                break;
        }
        while (!parsedCigar_.empty()) {
            if (parsedCigar_.front().type_ == 'P')
                parsedCigar_.erase(parsedCigar_.begin());
            else
                break;
        }
    }

    namespace {
        inline char collapsibleType(char ch) 
        {
            switch (ch)
            {
                case 'I':
                case 'P':
                    return 'I';
                case 'D':
                    return 'D';
                default:
                    return 0;
            }
        }
    }

    cgatools::mapping::Cigar Cigar::pack() const 
    {
        Cigar result(true,true);
        BOOST_FOREACH(const CigarElement &e, parsedCigar_)
        {
            result.add(e);
            size_t prevIndex = result.getParsedCigar().size()-1;
            while (prevIndex>0)
            {
                CigarElement &lastE = result.parsedCigar_[prevIndex];
                char lastType = collapsibleType(lastE.type_);
                if (lastType>0)
                {
                    --prevIndex;
                    CigarElement &prevE = result.parsedCigar_[prevIndex];
                    char prevType = collapsibleType(prevE.type_);
                    if (prevType==0 || lastType==prevType)
                        break;

                    //collapse insertion and deletion and generate 'M'/'N' to replace the overlapping part
                    bool prevEShorter = prevE.length_<lastE.length_;
                    CigarElement &shortE = prevEShorter ? prevE : lastE;
                    CigarElement &longE = !prevEShorter ? prevE : lastE;
                    longE.length_ -= shortE.length_;
                    CGA_ASSERT_MSG(longE.type_=='D' || shortE.type_=='D',
                        CGA_VOUT(longE.type_)<<CGA_VOUT(shortE.type_));
                    if (longE.type_=='I' || shortE.type_=='I')
                    {
                        CGA_ASSERT_MSG(longE.type_!='P' && shortE.type_!='P',
                            CGA_VOUT(longE.type_)<<CGA_VOUT(shortE.type_));
                        shortE.type_ = 'M';
                    } else 
                    {
                        CGA_ASSERT_MSG(longE.type_!='I' && shortE.type_!='I',
                            CGA_VOUT(longE.type_)<<CGA_VOUT(shortE.type_));
                        shortE.type_ = 'N';
                    }
                    if (prevE.length_==0) 
                    {
                        result.parsedCigar_.erase(result.parsedCigar_.begin()+prevIndex);
                        break;
                    } else if (prevIndex>0) 
                    {
                        CigarElement &checkE = result.parsedCigar_[prevIndex-1];
                        if (checkE.type_==prevE.type_) {
                            checkE.length_ += prevE.length_;
                            result.parsedCigar_.erase(result.parsedCigar_.begin()+prevIndex);
                        } else {
                            std::swap(result[prevIndex],result[prevIndex+1]);
                        }
                    }
                } else
                    break;
            }
        }
        return result;
    }

    void Cigar::add( const CigarElement &e )
    {
        if (e.length_!=0 || !removeZeros_) 
        {
            if (mergeNeighbours_)
                merge_back(e);
            else 
                push_back(e);
        }
    }

    std::ostream& operator<< (std::ostream& out, const Cigar& cigar) {
        BOOST_FOREACH(const Cigar::CigarElement& e, cigar.getParsedCigar()) {
            out << e.length_ << e.type_;
        }
        return out;
    }

} } // cgatools::mapping
