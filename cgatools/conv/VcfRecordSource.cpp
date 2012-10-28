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

#include "cgatools/core.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"

#include <boost/algorithm/string.hpp>

using namespace std;

using cgatools::reference::Location;

namespace cgatools { namespace conv {

    void VcfRecordWriter::writeRecord(
        std::ostream& out,
        const cgatools::reference::CrrFile& crr,
        size_t sampleCount) const
    {
        Location loc = getLocation();

        // CHROM
        string chrom = getVcfChromosomeName(loc.chromosome_, crr);
        out << chrom;

        // POS
        out << "\t" << (loc.offset_+1);

        // ID
        out << "\t";
        writeId(out);

        // REF
        out << "\t";
        writeRef(out);

        // ALT
        out << "\t";
        writeAlt(out);

        // QUAL
        out << "\t";
        writeQual(out);

        // FILTER
        out << "\t";
        writeFilter(out);

        // INFO
        out << "\t";
        writeInfo(out);

        // FORMAT
        out << "\t";
        writeFormat(out);

        // Samples
        for(size_t ii=0; ii<sampleCount; ii++)
        {
            out << "\t";
            writeSample(out, ii);
        }

        out << endl;
    }

    std::string VcfRecordWriter::getVcfChromosomeName(
        size_t chromosomeId,
        const cgatools::reference::CrrFile& crr)
    {
        CGA_ASSERT(chromosomeId < crr.listChromosomes().size());
        string chrom = crr.listChromosomes()[chromosomeId].getName();
        if (boost::starts_with(chrom, "chr"))
            chrom = chrom.substr(3);
        return chrom;
    }

} } // cgatools::conv
