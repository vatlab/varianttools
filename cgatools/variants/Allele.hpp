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

#ifndef CGATOOLS_VARIANTS_ALLELE_HPP_
#define CGATOOLS_VARIANTS_ALLELE_HPP_ 1

//! @file Allele.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/Call.hpp"

#include <string>
#include <vector>

namespace cgatools { namespace variants {

    class Locus;

    //! A class that corresponds to a single allele within a locus in
    //! a Complete Genomics variant file.
    class Allele
    {
    public:
        //! Construct the Allele. The Allele cannot be used until its
        //! locus is set and its call offsets are set.
        Allele()
            : locus_(0)
        {
        }

        //! Construct the Allele. The Allele cannot be used until its
        //! call offsets are set.
        Allele(const Locus* locus)
            : locus_(locus)
        {
        }

        //! Returns the offsets of the calls for this Allele within
        //! locus_->getCalls().
        const std::vector<size_t>& getCallOffsets() const
        {
            return callOffsets_;
        }

        const Call& getCall(size_t index) const;

        //! Returns the sequence of this Allele of this Locus, replacing
        //! "=" with the corresponding reference sequence.
        std::string getAlleleSequence() const;

        //! Returns the reference sequence corresponding to this Locus,
        //! replacing "=" with the corresponding reference sequence.
        std::string getRefSequence() const;

        //! Returns the first hapLink specified for this Allele, in call
        //! order.
        const std::string& getHapLink() const;

        //! Returns the union of xRef annotations for all Calls in this
        //! Allele, joined by ";".
        std::string getXRef() const;

        //! Returns the minimum varScoreVAF of all calls in this Allele.
        int32_t getMinVarScoreVAF() const;

        //! Returns the minimum varScoreEAF of all calls in this Allele.
        int32_t getMinVarScoreEAF() const;

        //! Returns the lowest varQuality of all calls in this Allele.
        Call::VarQuality getMinVarQuality() const;

        //! Returns true iff this Allele has one or more calls whose
        //! sequence contains N or ? characters.
        bool hasNoCalls() const;

        //! Clears out the set of calls associated with this Allele.
        void clearCalls();

        //! Associates a call in locus_->getCalls() with this Allele.
        void addCallOffset(size_t offset);

        //! Sets the locus for this Allele.
        void setLocus(const Locus* locus);

    private:
        const Locus* locus_;
        std::vector<size_t> callOffsets_;
    };

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_ALLELE_HPP_
