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
#include "cgatools/command/TestVariants.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"
#include "cgatools/variants/SuperlocusIterator.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/util/BaseUtil.hpp"

#include <queue>
#include <boost/foreach.hpp>

namespace cgatools { namespace command {

    using namespace std;
    using namespace cgatools::util;
    using namespace cgatools::reference;
    using namespace cgatools::variants;

    using boost::shared_ptr;

    TestVariants::TestVariants(const std::string& name)
        : Command(name,
                  "Tests variant files for presence of variants.",
                  "0.3 or later",

                  "Tests variant files for presence of variants. The output is a tab-delimited "
                  "file consisting of the columns of the input variants file, plus a column for "
                  "each assembly results file that contains a character code for each allele. "
                  "The character codes have meaning as follows:\n\n"
                  "    0 \tThis allele of this genome is consistent with the reference at this locus "
                  "but inconsistent with the variant.\n"
                  "    1 \tThis allele of this genome has the input variant at this locus.\n"
                  "    N \tThis allele of this genome has no-calls but is consistent with "
                  "the input variant.\n"
            ),
          maxHypothesisCount_(32)
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file.")
            ("input", po::value<string>(&inputFileName_)->default_value("STDIN"),
             "The input variants to test for.")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ("variants", po::value< vector<string> >(&variantFileNames_),
             "The input variant files (may be passed in as arguments at the end of the command).")
            ;

        positionalOptions_.add("variants", -1);
    }

    int TestVariants::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "variants");

        CrrFile crr(referenceFileName_);
        std::istream& in = openStdin(inputFileName_);
        std::ostream& out = openStdout(outputFileName_);

        vector< shared_ptr<VariantFileIterator> > vfi;
        vector< shared_ptr<SuperlocusIterator> > sli;
        BOOST_FOREACH(const string& variantFileName, variantFileNames_)
        {
            shared_ptr<VariantFileIterator> it(new VariantFileIterator(crr));
            it->open(variantFileName);
            vfi.push_back(it);

            shared_ptr<SuperlocusIterator> slit(new SuperlocusIterator(6, 0));
            slit->setVariantFile(*it);
            sli.push_back(slit);
        }

        Range range;
        string reference;
        string alleleSeq;
        DelimitedFile df(in, inputFileName_);
        df.addField(ChromosomeIdField("chromosome", &range.chromosome_, crr));
        df.addField(ValueField<uint32_t>("begin", &range.begin_));
        df.addField(ValueField<uint32_t>("end", &range.end_));
        df.addField(StringField("reference", &reference), DelimitedFile::FPT_OPTIONAL);
        df.addField(StringField("alleleSeq", &alleleSeq));
        writeHeader(out, df.getLine(), vfi);
        while (df.next())
        {
            if (reference != "" && reference != "?" && reference != crr.getSequence(range))
                throw Exception("variant list reference sequence mismatch: "+df.getLine());
            out << df.getLine();
            for(size_t ii=0; ii<sli.size(); ii++)
            {
                SuperlocusIterator& slit = *(sli[ii]);
                slit.skipToVariant(range, alleleSeq);
                const Superlocus& sl = *slit;
                vector< vector<PhasedHypothesis> > hypotheses;
                sl.buildPhasedHypotheses(hypotheses, maxHypothesisCount_, true);

                string result;
                PhasedHypothesis::testVariant(
                    sl, hypotheses[0],
                    range,
                    alleleSeq,
                    crr,
                    result);
                out << "\t" << result;
            }
            out << "\n";
        }

        return 0;
    }

    void TestVariants::writeHeader(std::ostream& out,
                                   const std::string& line,
                                   const std::vector< boost::shared_ptr<variants::VariantFileIterator> >& vfi)
    {
        out << line;
        for(size_t ii=0; ii<vfi.size(); ii++)
        {
            const DelimitedFile::Metadata& meta = vfi[ii]->getMetadata();
            if (!meta.hasKey("ASSEMBLY_ID"))
                out << "\tASM" << (ii+1);
            else
                out << "\t" << meta.get("ASSEMBLY_ID");
        }
        out << "\n";
    }

} } // cgatools::command
