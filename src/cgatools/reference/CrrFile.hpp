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

#ifndef CGATOOLS_REFERENCE_CRRFILE_HPP_
#define CGATOOLS_REFERENCE_CRRFILE_HPP_ 1

//! @file CrrFile.hpp
//! File containing definitions CrrFile class.

#include "cgatools/core.hpp"
#include "cgatools/reference/range.hpp"
#include "cgatools/reference/CrrFileWriter.hpp"
#include "cgatools/reference/CompactDnaSequence.hpp"

#include <vector>
#include <string>
#include <boost/iostreams/device/mapped_file.hpp>

namespace cgatools { namespace reference {

    //! Compact Randomly-accessible Reference file. It's a binary
    //! format.
    class CrrFile
    {
    public:
        CrrFile();
        CrrFile(const std::string& path);
        void open(const std::string& path);

        //! Return the list of chromosomes (and their metadata).
        const std::vector<CompactDnaSequence>& listChromosomes() const;

        //! Return the list of ranges that may be considered contigs, if
        //! you split contigs by stretches of reference N's that are at
        //! list minGapLength long.
        std::vector<Range> listContigs(uint32_t minGapLength = 50) const;

        //! Get the reference sequence for the given range.
        std::string getSequence(const Range& range) const;

        //! Get the reference sequence for the given location.
        char getBase(const Location& loc) const;

        //! Get the chromosome id for the given chromosome name.
        uint16_t getChromosomeId(const std::string& chromosomeName) const;

        //! Validate the md5 checksums in the chromosome metadata.
        void validate() const;

        //! The expected .crr file version.
        static uint32_t currentVersion();

        //! Returns the internal location by chromosome name
        //! and the zero-based offset into that chromosome.
        Location getLocation(const std::string& chromosomeName,
            boost::uint32_t offset) const
        {
            return Location(getChromosomeId(chromosomeName), offset);
        }

    private:
        void readHeader(const std::string& path, std::istream& in, uint64_t* chrTableOffset);
        void readChromosomeTable(std::istream& in,
                                 std::vector<CrrFileWriter::ChromosomeInfo>* pInfo,
                                 std::vector<util::Md5Digest>* pDigests);

        boost::shared_ptr<boost::iostreams::mapped_file_source> file_;
        std::vector<CompactDnaSequence> chromosomes_;
    };

} } // cgatools::reference

#endif // CGATOOLS_REFERENCE_CRRFILE_HPP_
