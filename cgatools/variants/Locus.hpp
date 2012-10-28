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

#ifndef CGATOOLS_VARIANTS_LOCUS_HPP_
#define CGATOOLS_VARIANTS_LOCUS_HPP_ 1

//! @file Locus.hpp

#include "cgatools/core.hpp"
#include "cgatools/variants/Call.hpp"
#include "cgatools/variants/Allele.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"

#include <map>

namespace cgatools { namespace variants {

    class LocusAnnotations
    {
    public:
        LocusAnnotations()
        {
        }

        LocusAnnotations(
            const std::vector<std::string>& annotationColumnHeaders,
            bool xRefIsAlleleSpecificInit);

        size_t count() const
        {
            return names_.size();
        }

        size_t getIndex(const std::string& name) const
        {
            for(size_t ii=0; ii<names_.size(); ii++)
            {
                if (names_[ii] == name)
                    return ii;
            }
            return names_.size();
        }

        size_t getSwapOrder(size_t idx) const
        {
            return swapOrder_[idx];
        }

        const std::string& emptyString() const
        {
            return emptyString_;
        }

        bool xRefIsAlleleSpecific() const
        {
            return xRefIsAlleleSpecific_;
        }

    private:
        std::vector< std::string > names_;
        std::vector< size_t > swapOrder_;
        std::string emptyString_;
        bool xRefIsAlleleSpecific_;
    };

    //! A class that corresponds to a single locus in a Complete
    //! Genomics variant file.
    class Locus
    {
        typedef reference::CrrFile CrrFile;
        typedef reference::Location Location;
        typedef reference::Range Range;
    public:
        //! Construct a dummy locus (for use with std::vector, for example).
        Locus()
            : crr_(0),
              ann_(new LocusAnnotations())
        {
        }

        //! Construct an empty locus. The Locus cannot be used until its
        //! calls have been added by Locus::addCall() and it has been
        //! initialized by Locus::initFromCalls().
        Locus(const CrrFile& crr)
            : crr_(&crr),
              ann_(new LocusAnnotations())
        {
        }

        //! Copy and assignment must be overridden because Allele has a
        //! pointer to this.
        Locus(const Locus& other);

        //! Copy and assignment must be overridden because Allele has a
        //! pointer to this.
        Locus& operator=(const Locus& other);

        //! Retrieve the calls for this locus, in order as found in the
        //! variant file.
        const std::vector<Call>& getCalls() const
        {
            return calls_;
        }

        //! Retrieve the calls for this locus, in order as found in the
        //! variant file.
        std::vector<Call>& getCalls()
        {
            return calls_;
        }

        //! Retrieve the alleles for this locus, in order by haplotype.
        const std::vector<Allele>& getAlleles() const
        {
            return alleles_;
        }

        //! Get the reference range for this locus.
        const Range& getRange() const
        {
            return range_;
        }

        //! Set range, only for loci with one call.
        void setRange(const Range& range);

        //! Get the reference crr file associated with this locus.
        const CrrFile& getReference() const
        {
            return *crr_;
        }

        //! Get the locus id (locus column in the variant file).
        boost::uint32_t getId() const;

        //! Set the locus id (locus column in the variant file). This
        //! updates all calls in the locus.
        void setId(boost::uint32_t locusId);

        //! Get the ploidy of this locus.
        boost::uint16_t getPloidy() const;

        //! Update all the hapLinks for the given allele of this locus.
        void setHapLink(size_t alleleOffset, const std::string& hapLink);

        //! Returns true iff this locus is a one-call reference-called
        //! locus.
        bool isRefCallLocus() const;

        //! Returns true iff this locus is a one-call no-called locus.
        bool isNoCallLocus() const;

        //! Returns true iff each call in this locus is consistent with
        //! the reference.
        bool isRefConsistent() const;

        //! Returns true iff this locus has one or more calls whose
        //! sequence contains N or ? characters.
        bool hasNoCalls() const;

        //! Returns the requested read counts, or -1 if this read count
        //! is empty as annotated in a masterVar file. This function
        //! computes read counts from the evidence, regardless of
        //! whether or not this locus is annotated with read counts. For
        //! a function that can use the masterVar annotation if
        //! available, see VariantFileIterator::getReadCounts().
        void computeReadCounts(int32_t* allele1ReadCount,
                               int32_t* allele2ReadCount,
                               int32_t* referenceAlleleReadCount,
                               int32_t* totalReadCount,
                               cgdata::EvidenceReader& evidence) const;

        //! Clears the set of calls for this locus.
        void clearCalls();

        //! Adds a call to the locus.
        void addCall(const Call& call);

        //! After adding calls by Locus::addCall(), this function must
        //! be called to prepare the Locus for use.
        void initFromCalls(bool relaxedReferenceValidation = true);

        //! Returns in the calls parameter the bases called by this
        //! locus at the specified location (or - or ., in accordance
        //! with snpdiff interpretation of those values), as well as the
        //! Call used to determine each base call.
        void locationCalls(const Location& loc, std::vector< std::pair<char, const Call*> >& calls) const;

        //! Returns locus type, as specified for simplified variation files.
        std::string getType() const;

        //! Sets the masterVar varType for this locus.
        void setType(const std::string& olplType);

        //! Returns locus zygosity, as specified for simplified variation files.
        std::string getZygosity() const;

        //! Changes the allele order to the one specified for simplified
        //! variation files (i.e. "more interesting allele first").
        void reorderAlleles();

        //! Writes this locus to the stream in the format specified for
        //! simplified one-line-per-locus variation files.
        void writeAsOneLine(std::ostream& out, bool writeExtras = true, char sep = '\t') const;

        //! Returns true iff this locus has the given masterVar annotation.
        bool hasAnnotation(const std::string& name) const;

        //! Returns the masterVar annotation with the given name. If the
        //! annotation is not available, returns empty string.
        const std::string& getAnnotation(const std::string& name) const;

        //! Returns the union of xRef annotations for all Calls in this
        //! Locus, joined by ";".
        std::string getAllXRef() const;

        //! Returns false if this is a masterVar with xRef annotation
        //! but not allele1XRef or allele2XRef. (True in older masterVar
        //! files.) Otherwise, returns true.
        bool xRefIsAlleleSpecific() const;

        //! Set the locus annotations for this locus.
        void setLocusAnnotations(const LocusAnnotations& ann);

        static void writeOneLineFileHeader(
            std::ostream& out, bool alleleSpecificXRef, char sep = '\t');
    private:
        void validateCalls(bool relaxedReferenceValidation) const;
        void validateAllelesAndInitRange();
        void callError (const std::string&, const Call&) const;
        void locusError(const std::string&) const;

        Range range_;
        const CrrFile* crr_;
        std::vector<Call> calls_;
        std::vector<Allele> alleles_;
        std::string olplType_; // varType, from olpl file
        boost::shared_ptr< LocusAnnotations > ann_;

        // For each interval allele, true if it's compatible to the object
        // of interest (variant allele or reference, depending on the context).
        typedef boost::array<bool, 4> EvidenceMatchMap;

        // For each locus allele, true if the allele is compatible to the
        // interval allele 0, 1 and 2. For alleles with no-calls, false in
        // all three slots.
        typedef std::vector< EvidenceMatchMap > AlleleMatchMap;

        // Returns true if the read count should be printed for the specified allele
        bool hasCount(const AlleleMatchMap& m, size_t alleleIdx) const;

        // Returns 1 or 0 depending on whether dnb strongly supports the interval alleles
        // marked as "true" in the match map over the alleles marked as "false".
        size_t countSupport(
            const cgdata::EvidenceReader::DnbRecord& dnb, const EvidenceMatchMap& mm) const;

        // For each allele in the locus, find the compatible alleles in the
        // evidence file interval
        void matchAlleles(const cgdata::EvidenceReader::IntervalRecord& interval,
                          bool& hasLengthChange,
                          AlleleMatchMap& m,
                          EvidenceMatchMap& rm) const;

        bool hasFullyCalledNonReferenceAllele() const;

    public:
        //! Additional annotation columns
        std::vector<std::string> extras_;
    };

    std::ostream& operator<<(std::ostream& out, const Locus& locus);

} } // cgatools::variants

#endif // CGATOOLS_VARIANTS_LOCUS_HPP_
