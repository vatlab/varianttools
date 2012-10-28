// Copyright 2011 Complete Genomics, Inc.
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

#ifndef CGA_TOOLS_CONV_VCFRECORDSOURCE_HPP_
#define CGA_TOOLS_CONV_VCFRECORDSOURCE_HPP_ 1

//! @file VcfRecordSource.hpp

#include "cgatools/core.hpp"
#include "cgatools/reference/CrrFile.hpp"

#include <iostream>

namespace cgatools { namespace conv {

    //! Abstract base class that writes a VCF record. Classes derived
    //! from this include SmallVariantsVcfRecordWriter, etc.
    class VcfRecordWriter
    {
    public:
        //! Write the VCF record to out. Not overridable in sub-classes,
        //! but this just calls all the per-column write() functions,
        //! which are overridable in sub-classes.
        void writeRecord(
            std::ostream& out,
            const cgatools::reference::CrrFile& crr,
            size_t sampleCount) const;

        //! Get 0-based location of first base of VCF record. The value
        //! output to the VCF file is at Location.offset+1.
        virtual cgatools::reference::Location getLocation() const = 0;

        //! Write ID field of VCF record.
        virtual void writeId(std::ostream& out) const
        {
            out << ".";
        }

        //! Write REF field of VCF record.
        virtual void writeRef(std::ostream& out) const = 0;

        //! Write ALT field of VCF record.
        virtual void writeAlt(std::ostream& out) const = 0;

        //! Write QUAL field of VCF record.
        virtual void writeQual(std::ostream& out) const
        {
            out << ".";
        }

        //! Write FILTER field of VCF record.
        virtual void writeFilter(std::ostream& out) const
        {
            out << ".";
        }

        //! Write INFO field of VCF record.
        virtual void writeInfo(std::ostream& out) const = 0;

        //! Write FORMAT field of VCF record.
        virtual void writeFormat(std::ostream& out) const = 0;

        //! Write per-sample field of VCF record, for the given
        //! genome. Genomes are numbered starting with 0 for the first
        //! genome.
        virtual void writeSample(std::ostream& out, size_t idxGenome) const = 0;

        //! Get the chromosome name for the given chromosome ID, as
        //! written in a VCF record. I.e. no "chr" prefix.
        static std::string getVcfChromosomeName(
            size_t chromosomeId,
            const cgatools::reference::CrrFile& crr);
    };

    struct VcfSubFieldHeaderRecord
    {
        enum Key
        {
            VCF_ALT    = 0,
            VCF_FILTER = 1,
            VCF_INFO   = 2,
            VCF_FORMAT = 3,
            VCF_END    = 4
        };

        VcfSubFieldHeaderRecord(
            Key key,
            const std::string& id,
            const std::string& number,
            const std::string& type,
            const std::string& description)
            : key_(key),
              id_(id),
              number_(number),
              type_(type),
              description_(description)
        {
        };

        Key         key_;
        std::string id_;
        std::string number_;
        std::string type_;
        std::string description_;
    };

    //! Struct to represent metadata records like:
    //! ##source_KEY=VALUE
    struct VcfKvHeaderRecord
    {
        VcfKvHeaderRecord(const std::string& key,
                          const std::string& value)
            : key_(key),
              value_(value)
        {
        }

        std::string key_;
        std::string value_;
    };

    //! Abstract base class that provides a stream of
    //! VcfRecordWriter. Classes derived from this include
    //! SmallVariantsVcfRecordSource, etc.
    class VcfRecordSource
    {
    public:
        //! Returns the VcfSubFieldHeaderRecords for all
        //! genomes. Includes, INFO, FORMAT, FILTER, and ALT records.
        virtual std::vector<VcfSubFieldHeaderRecord> getSubFieldHeaderRecords() const = 0;

        //! Returns the ##source=VALUE header value for the given
        //! genome, which is the same as the input file's
        //! getMetadata().getSoftwareVersionString().
        virtual std::string getSource(size_t idxGenome) const = 0;

        //! Returns ##KEY=VALUE records for the given genome. Must
        //! include at least the source_GENOME_REFERENCE key.
        virtual std::vector<VcfKvHeaderRecord> getKeyValueHeaderRecords(size_t idxGenome) const = 0;

        //! Returns the input file's getMetadata().get("ASSEMBLY_ID").
        virtual std::string getAssemblyId(size_t idxGenome) const = 0;

        //! Record iteration

        virtual bool eof() const = 0;
        virtual VcfRecordSource& operator++() = 0;
        virtual const VcfRecordWriter& operator*() const = 0;
        virtual const VcfRecordWriter* operator->() const = 0;
    };

} } // cgatools::conv

#endif // CGA_TOOLS_CONV_VCFRECORDSOURCE_HPP_
