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

#ifndef CGATOOLS_REFERENCE_CRRFILEWRITER_HPP_
#define CGATOOLS_REFERENCE_CRRFILEWRITER_HPP_ 1

//! @file CrrFileWriter.hpp
//! File containing definition of the CrrFileWriter class.

#include "cgatools/core.hpp"
#include "cgatools/util/Md5.hpp"
#include "cgatools/reference/CompactDnaSequence.hpp"

namespace cgatools { namespace reference {

    class CrrFile;

    //! Writes a CrrFile to an output stream.
    class CrrFileWriter
    {
    public:
        //! Construct this CrrFileWriter -- out must be randomly
        //! accessible (not stdout).
        CrrFileWriter(std::ostream* out);
        ~CrrFileWriter();

        //! Start adding sequence for a new chromosome.
        void newChromosome(const std::string& name, bool circular);

        //! Function to stream additional sequence into the CrrFileWriter.
        void addSequence(const std::string& sequence);

        //! Finish writing the CrrFile (automatically called by
        //! destructor).
        void close();

    private:
        struct ChromosomeInfo
        {
            ChromosomeInfo();
            ChromosomeInfo(const std::string& name, bool circular, uint64_t fileOffset);

            std::string name_;
            bool circular_;
            uint64_t fileOffset_;
            util::Md5Context md5_;
            uint64_t length_;
            std::vector<AmbiguousRegion> amb_;
        };

        void writeHeader(uint64_t chrTableOffset);
        void writeGuard();
        void writeChromosomeTable();
        void endChromosome();
        void addBase(char base);

        std::ostream* out_;
        std::vector<ChromosomeInfo> chromosomes_;
        uint8_t packedBases_;
        uint32_t packedBaseCount_;
        bool closed_;

        friend class CrrFile;
    };

} } // cgatools::reference

#endif // CGATOOLS_REFERENCE_CRRFILEWRITER_HPP_
