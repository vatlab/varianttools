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

#ifndef CGA_TOOLS_CIGAR_HPP_
#define CGA_TOOLS_CIGAR_HPP_ 1

//! @file Cigar.hpp

#include "cgatools/core.hpp"

#include <vector>
#include <string>
#include <map>
#include <boost/shared_ptr.hpp>

namespace cgatools { namespace mapping {

    //! Provides support for manipulations on SAM/BAM and CGI Evidence CIGAR
    class Cigar {
        friend std::ostream& operator<< (std::ostream& out, const Cigar& cigar);
    public:
        class CigarElement {
        public:
            CigarElement(char type, size_t length)
                : type_(type), length_(length) {}
            int getSequenceLength() const;
            bool operator== (const CigarElement& e) const {return type_==e.type_ && length_==e.length_;}

            char    type_;
            size_t  length_;
        };
        typedef std::vector<CigarElement> ParsedCigar;

        Cigar(bool mergeNeighbours = false, bool removeZeros = false) 
            :mergeNeighbours_(mergeNeighbours), removeZeros_(removeZeros)
        {}

        Cigar(const std::string &cigar, bool mergeNeighbours = false, bool removeZeros = false) 
            :mergeNeighbours_(mergeNeighbours), removeZeros_(removeZeros)
        {
            parse(cigar);
        }

        void parse(const std::string &cigar);

        //! computes length of the sequence the cigar describes
        size_t getSequenceLength() const;
        //! computes length of the reference segment the cigar describes
        size_t getReferenceLength() const;

        const ParsedCigar& getParsedCigar() const {return parsedCigar_;}

        //adds or merges based on the parameter given in constructor
        void add(const CigarElement &e);

        //! adds a new cigar element
        void push_back(const CigarElement &e) {parsedCigar_.push_back(e);}
        //! merges the last and new elements if they have equal types
        void merge_back(const CigarElement &e);

        CigarElement &operator[](size_t i) {return parsedCigar_[i];}
        CigarElement &back() {return parsedCigar_.back();}

        size_t size() const {return parsedCigar_.size();}
        //removes paddings in the beginning and at the end of the CIGAR to satisfy SAMValidate from Pickard
        void trancatePaddings();

        //! Removes zero size section merges neighbor sections of the same type and returns a copy
        Cigar pack() const;

        bool operator== (const Cigar& c) const {return parsedCigar_==c.parsedCigar_;}
    protected:

        bool mergeNeighbours_;
        bool removeZeros_;
        ParsedCigar parsedCigar_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_CIGAR_HPP_
