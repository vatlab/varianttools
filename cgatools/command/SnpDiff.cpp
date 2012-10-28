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
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/command/SnpDiff.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"

#include <boost/foreach.hpp>
#include <boost/format.hpp>

namespace cgatools { namespace command {

    using namespace cgatools::util;
    namespace bu = baseutil;

    using reference::Location;
    using reference::ChromosomeIdField;
    using variants::Call;
    using variants::VariantFileIterator;

    using namespace std;

    using boost::shared_ptr;

    enum GenotypeClass
    {
        GENOTYPE_CLASS_REF_REF     = 0,
        GENOTYPE_CLASS_HET_REF_ALT = 1,
        GENOTYPE_CLASS_HET_ALT_ALT = 2,
        GENOTYPE_CLASS_HOM_ALT_ALT = 3,
        GENOTYPE_CLASS_OTHER       = 4,
        GENOTYPE_CLASS_COUNT
    };

    static const char* genotypeClassName[] =
    {
        "ref-ref",
        "het-ref-alt",
        "het-alt-alt",
        "hom-alt-alt",
        "other"
    };

    struct GenotypeCall
    {
        reference::Location loc_;
        std::string genotypes_;
        std::string line_;

        GenotypeCall(const reference::Location& loc,
                     const std::string& genotypes,
                     const std::string& line)
            : loc_(loc),
              genotypes_(genotypes),
              line_(line)
        {
        }

        GenotypeClass classifyGenotypes(char refBase) const
        {
            if (2 != genotypes_.size())
                return GENOTYPE_CLASS_OTHER;
            if ( (!bu::isValidBase(genotypes_[0])) || (!bu::isValidBase(genotypes_[1])) )
                return GENOTYPE_CLASS_OTHER;

            if ( bu::isConsistent(genotypes_[0], refBase) &&
                 bu::isConsistent(genotypes_[1], refBase) )
                return GENOTYPE_CLASS_REF_REF;
            if ( bu::isConsistent(genotypes_[0], refBase) ||
                 bu::isConsistent(genotypes_[1], refBase) )
                return GENOTYPE_CLASS_HET_REF_ALT;
            if ( genotypes_[0] == genotypes_[1] )
                return GENOTYPE_CLASS_HOM_ALT_ALT;
            return GENOTYPE_CLASS_HET_ALT_ALT;
        }

        //! Returns "match class" for this genotypes line, given the
        //! variant calls. Returns one match class code per allele in
        //! the genotypes. The codes are as follows:
        //! - C for concordant
        //! - D for discordant
        //! - N for no-call
        //! - . for other
        std::string getMatchClass(const std::string& variantCalls) const
        {
            string vc = variantCalls;
            if (vc.size() != genotypes_.size())
            {
                if (0 == genotypes_.size() || 0 != genotypes_.size() % vc.size())
                    return string(genotypes_.size(), '.');

                size_t count = genotypes_.size() / vc.size() - 1;
                for(size_t ii=0; ii<count; ii++)
                    vc += variantCalls;
            }

            CGA_ASSERT(vc.size() == genotypes_.size());
            vector<size_t> gIdx;
            gIdx.reserve(vc.size());
            for(size_t ii=0; ii<vc.size(); ii++)
                gIdx.push_back(ii);
            size_t minDiscordant = gIdx.size()+1;
            string bestMatch(gIdx.size(), '.');
            string currentMatch(gIdx.size(), '.');
            do
            {
                size_t discordant = 0;
                for(size_t ii=0; ii<gIdx.size(); ii++)
                {
                    if ('N' == vc[ii])
                        currentMatch[ii] = 'N';
                    else if (vc[ii] == genotypes_[gIdx[ii]])
                        currentMatch[ii] = 'C';
                    else
                    {
                        if ('-' == vc[ii])
                            currentMatch[ii] = '.';
                        else if ('.' == vc[ii])
                            currentMatch[ii] = '.';
                        else if (bu::isConsistent(genotypes_[gIdx[ii]], vc[ii]))
                            currentMatch[ii] = 'C';
                        else
                            currentMatch[ii] = 'D';
                        if (currentMatch[ii] != 'C')
                            discordant++;
                    }
                }
                if (discordant < minDiscordant)
                {
                    minDiscordant = discordant;
                    bestMatch = currentMatch;
                }
            } while(std::next_permutation(gIdx.begin(), gIdx.end()));

            std::sort(bestMatch.begin(), bestMatch.end());
            return bestMatch;
        }

        bool operator<(const GenotypeCall& a) const
        {
            return loc_ < a.loc_;
        }
    };

    class GenotypesField : public StringField
    {
    public:
        GenotypesField(const std::string& name, std::string* val)
            : StringField(name, val),
              vv_(val),
              allowed_("ACGTacgtNn-.")
        {
        }

        void parse(const char* first, const char* last)
        {
            StringField::parse(first, last);
            if (vv_->size() > 3)
                throw Exception("invalid genotypes field: too many genotypes");
            for(size_t ii=0; ii<vv_->size(); ii++)
            {
                if ( std::find(allowed_.begin(), allowed_.end(), (*vv_)[ii]) ==
                     allowed_.end() )
                    throw Exception("invalid genotypes field: invalid character");
            }
        }

    private:
        std::string* vv_;
        std::string allowed_;
    };

    class GenotypeFile
    {
    public:
        GenotypeFile(const std::string& fn, const reference::CrrFile& crr)
            : fn_(fn),
              strand_(false)
        {
            load(crr);
        }

        std::vector<GenotypeCall>& getData()
        {
            return data_;
        }

        const std::vector<std::string>& getColumnHeaders() const
        {
            return columnHeaders_;
        }

        bool hasGenotypes() const
        {
            return hasGenotypes_;
        }

    private:
        std::string fn_;
        std::vector<GenotypeCall> data_;
        reference::Location loc_;
        std::string genotypes_;
        bool strand_;
        std::vector<string> columnHeaders_;
        bool hasGenotypes_;

        void load(const reference::CrrFile& crr)
        {
            shared_ptr<istream> in(InputStream::openCompressedInputStreamByExtension(fn_));

            DelimitedFile df(*in, fn_);
            parseHeader(df, crr);
            columnHeaders_ = df.getColumnHeaders();
            while (df.next())
            {
                string genotypes = genotypes_;
                if (strand_)
                    genotypes = bu::reverseComplement(genotypes_);
                std::sort(genotypes.begin(), genotypes.end());
                data_.push_back(GenotypeCall(loc_, genotypes, df.getLine()));
            }

            std::sort(data_.begin(), data_.end());
        }

        void parseHeader(util::DelimitedFile& df, const reference::CrrFile& crr)
        {
            df.addField(ChromosomeIdField("Chromosome", &loc_.chromosome_, crr));
            df.addField(ValueField<uint32_t>("Offset0Based", &loc_.offset_));
            df.addField(GenotypesField("Genotypes", &genotypes_), DelimitedFile::FPT_OPTIONAL);
            df.addField(StrandField("GenotypesStrand", &strand_), DelimitedFile::FPT_OPTIONAL);
            hasGenotypes_ = false;
            const vector<string>& ch = df.getColumnHeaders();
            for(size_t ii=0; ii<ch.size(); ii++)
            {
                if (ch[ii] == "Genotypes")
                {
                    hasGenotypes_ = true;
                    break;
                }
            }
        }
    };

    SnpDiff::SnpDiff(const std::string& name)
        : Command(name,
                  "Compares snp calls to a Complete Genomics variant file.",
                  "0.3 or later",
                  "Compares the snp calls in the \"genotypes\" file to the "
                  "calls in a Complete Genomics variant file. The genotypes file is "
                  "a tab-delimited file with at least the following columns "
                  "(additional columns may be given):\n\n"
                  "    Chromosome      \t(Required) The name of the chromosome.\n"
                  "    Offset0Based    \t(Required) The 0-based offset in the chromosome.\n"
                  "    GenotypesStrand \t(Optional) The strand of the calls in the "
                  "Genotypes column (+ or -, defaults to +).\n"
                  "    Genotypes       \t(Optional) The calls, one per allele. The following "
                  "calls are recognized:\n"
                  "                    A,C,G,T \tA called base.\n"
                  "                    N       \tA no-call.\n"
                  "                    -       \tA deleted base.\n"
                  "                    .       \tA non-snp variation.\n"
                  "\n"
                  "The output is a tab-delimited file consisting of the columns of the "
                  "original genotypes file, plus the following additional columns:\n\n"
                  "    Reference         \tThe reference base at the given position.\n"
                  "    VariantFile       \tThe calls made by the variant file, one per "
                  "allele. The character codes are the same as is described for the "
                  "Genotypes column.\n"
                  "    DiscordantAlleles \t(Only if Genotypes is present) The number of "
                  "Genotypes alleles that are discordant with calls in the VariantFile. If "
                  "the VariantFile is described as haploid at the given position but the "
                  "Genotypes is diploid, then each genotype allele is compared against the "
                  "haploid call of the VariantFile.\n"
                  "    NoCallAlleles     \t(Only if Genotypes is present) The number of "
                  "Genotypes alleles that were no-called by the VariantFile. If "
                  "the VariantFile is described as haploid at the given position but the "
                  "Genotypes is diploid, then a VariantFile no-call is counted twice.\n"
                  "\n"
                  "The verbose output is a tab-delimited file consisting of the columns "
                  "of the original genotypes file, plus the following additional columns:\n\n"
                  "    Reference   \tThe reference base at the given position.\n"
                  "    VariantFile \tThe call made by the variant file for one allele "
                  "(there is a line in this file for each allele). The character codes "
                  "are the same as is described for the Genotypes column.\n"
                  "    [CALLS]     \tThe rest of the columns are pasted in from the "
                  "VariantFile, describing the variant file line used to make the call.\n"
                  "\n"
                  "The stats output is a comma-separated file with several tables "
                  "describing the results of the snp comparison, for each diploid "
                  "genotype. The tables all describe the comparison result (column "
                  "headers) versus the genotype classification (row labels) in "
                  "different ways. The \"Locus classification\" tables have the most "
                  "detailed match classifications, while the \"Locus concordance\" "
                  "tables roll these match classifications up into \"discordance\" "
                  "and \"no-call\". A locus is considered discordant if it is "
                  "discordant for either allele. A locus is considered no-call if "
                  "it is concordant for both alleles but has a no-call on either allele. "
                  "The \"Allele concordance\" describes the comparison result on a "
                  "per-allele basis."
            )
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The input crr file.")
            ("variants", po::value<string>(&variantFileName_),
             "The input variant file.")
            ("genotypes", po::value<string>(&genotypesFileName_),
             "The input genotypes file.")
            ("output-prefix", po::value<string>(&oPrefix_),
             "The path prefix for all output reports.")
            ("reports", po::value<string>(&reports_)->default_value("Output,Verbose,Stats"),
             "Comma-separated list of reports to generate. A report is one of:\n"
             "    Output  The output genotypes file.\t\n"
             "    Verbose The verbose output file.\t\n"
             "    Stats   The stats output file.\t\n")
             ;
    }

    int SnpDiff::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "variants");
        requireParam(vm, "genotypes");

        StringSet reports(reports_,
                          "Output,Verbose,Stats",
                          "unknown report type");

        reference::CrrFile crr(referenceFileName_);

        GenotypeFile genotypes(genotypesFileName_, crr);
        vector< map<string,int> > stats(GENOTYPE_CLASS_COUNT);

        VariantFileIterator locIt(crr);
        locIt.open(variantFileName_);

        string asmId;
        if (locIt.getMetadata().hasKey("ASSEMBLY_ID"))
            asmId = locIt.getMetadata().get("ASSEMBLY_ID");
        else
            asmId = variantFileName_;

        shared_ptr<std::ostream> out, vOut;

        if (0 != reports.count("Output"))
        {
            out =
                OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "Output.tsv");

            const vector<string>& gch = genotypes.getColumnHeaders();
            for(size_t ii=0; ii<gch.size(); ii++)
                *out << gch[ii] << "\t";
            *out << "Reference";
            *out << "\tVariantFile";
            if (genotypes.hasGenotypes())
                *out << "\tDiscordantAlleles\tNoCallAlleles";
            *out << "\n";
        }

        if (0 != reports.count("Verbose"))
        {
            vOut = OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "Verbose.tsv");

            const vector<string>& gch = genotypes.getColumnHeaders();
            for(size_t ii=0; ii<gch.size(); ii++)
                *vOut << gch[ii] << "\t";
            *vOut << "Reference\tVariantFile";
            *vOut << "\t" << Call::getHeader();
            *vOut << "\n";
        }
        else
            vOut.reset(static_cast<std::ostream*>(0));

        BOOST_FOREACH(GenotypeCall& gc, genotypes.getData())
        {
            char refBase = crr.getBase(gc.loc_);

            while ( (!locIt.eof()) &&
                    locIt->getRange().endLocation() <= gc.loc_ )
                ++locIt;

            vector< std::pair<char, const variants::Call*> > calls;
            if ( (!locIt.eof()) && locIt->getRange().contains(gc.loc_) )
                locIt->locationCalls(gc.loc_, calls);
            else
                calls.push_back(pair<char, const variants::Call*>('N', 0));
            std::sort(calls.begin(), calls.end());

            string callString(calls.size(), '.');
            for(size_t ii=0; ii<calls.size(); ii++)
                callString[ii] = calls[ii].first;
            string matchClass = gc.getMatchClass(callString);

            if (0 != vOut.get())
            {
                for(size_t ii=0; ii<calls.size(); ii++)
                {
                    *vOut << gc.line_ << "\t";
                    *vOut << refBase << "\t" << calls[ii].first << "\t";
                    if (0 != calls[ii].second)
                        calls[ii].second->write(*vOut, crr);
                    *vOut << "\n";
                }
            }

            if (0 != out.get())
            {
                *out << gc.line_ << "\t";
                *out << refBase << "\t";
                for(size_t ii=0; ii<calls.size(); ii++)
                    *out << calls[ii].first;
                if (genotypes.hasGenotypes())
                {
                    int discordant = 0;
                    int noCall = 0;
                    for(size_t ii=0; ii<matchClass.size(); ii++)
                    {
                        switch(matchClass[ii])
                        {
                        case 'D': case '.':
                            discordant++;
                            break;
                        case 'N':
                            noCall++;
                            break;
                        default:
                            break;
                        }
                    }
                    *out << "\t" << discordant << "\t" << noCall;
                }
                *out << "\n";
            }

            if (genotypes.hasGenotypes() && 2 == gc.genotypes_.size())
            {
                // Collect stats.
                stats[gc.classifyGenotypes(refBase)][matchClass]++;
            }
        }

        // Finish processing the entire file, as a sanity check that the
        // file is well-formatted.
        while (!locIt.eof())
            ++locIt;

        // Write stats.
        if (0 != reports.count("Stats"))
        {
            shared_ptr<std::ostream> sOut =
                OutputStream::openCompressedOutputStreamByExtension(oPrefix_ + "Stats.tsv");
            vector< std::pair<string,string> > matchClassNames;
            matchClassNames.push_back(pair<string,string>("CC", "match-match"));
            matchClassNames.push_back(pair<string,string>("CN", "nocall-match"));
            matchClassNames.push_back(pair<string,string>("NN", "nocall-nocall"));
            matchClassNames.push_back(pair<string,string>("CD", "match-mismatch"));
            matchClassNames.push_back(pair<string,string>("DD", "mismatch-mismatch"));
            matchClassNames.push_back(pair<string,string>("DN", "nocall-mismatch"));
            matchClassNames.push_back(pair<string,string>(".C", "match-nonsnp"));
            matchClassNames.push_back(pair<string,string>(".D", "mismatch-nonsnp"));
            matchClassNames.push_back(pair<string,string>("..", "nonsnp-nonsnp"));
            matchClassNames.push_back(pair<string,string>(".N", "nocall-nonsnp"));

            printMatchTable(*sOut, stats, matchClassNames, true);
            *sOut << "\n";
            printMatchTable(*sOut, stats, matchClassNames, false);
            *sOut << "\n";

            printConcordanceTable(*sOut, stats, true, false);
            *sOut << "\n";
            printConcordanceTable(*sOut, stats, false, false);
            *sOut << "\n";

            printConcordanceTable(*sOut, stats, true, true);
            *sOut << "\n";
            printConcordanceTable(*sOut, stats, false, true);
        }

        return 0;
    }

    std::string SnpDiff::formatCount(uint32_t count, uint32_t total, bool fraction) const
    {
        if (fraction)
        {
            double val = 0.0;
            if (0 == total)
            {
                CGA_ASSERT(0 == count);
            }
            else
                val = double(count) / double(total);
            return (boost::format("%.5f") % val).str();
        }
        return boost::lexical_cast<std::string>(count);
    }

    void SnpDiff::printMatchRow(std::ostream& out,
                                const char* rowName,
                                std::map<std::string,int>& matchTable,
                                const std::vector< std::pair<std::string,std::string> >& matchClassNames,
                                bool fraction) const
    {
        out << rowName;
        uint32_t total = 0;
        for(size_t jj=0; jj<matchClassNames.size(); jj++)
            total += matchTable[matchClassNames[jj].first];
        for(size_t jj=0; jj<matchClassNames.size(); jj++)
            out << "\t" << formatCount(matchTable[matchClassNames[jj].first], total, fraction);
        out << "\t" << formatCount(total, total, fraction);
        out << "\n";
    }

    void SnpDiff::printMatchTable(std::ostream& out,
                                  std::vector< std::map<std::string,int> >& stats,
                                  const std::vector< std::pair<std::string,std::string> >& matchClassNames,
                                  bool fraction) const
    {
        // Print header.
        out << "Locus classification ";
        if (fraction)
            out << "by fraction of total";
        else
            out << "by count";
        out << "\nGenotype";
        for(size_t jj=0; jj<matchClassNames.size(); jj++)
            out << "\t" << matchClassNames[jj].second;
        out << "\ttotal\n";

        // Print table.
        map<string,int> totals;
        for(int ii=0; ii<GENOTYPE_CLASS_OTHER; ii++)
        {
            for(size_t jj=0; jj<matchClassNames.size(); jj++)
                totals[matchClassNames[jj].first] += stats[ii][matchClassNames[jj].first];

            printMatchRow(out, genotypeClassName[ii], stats[ii], matchClassNames, fraction);
        }
        printMatchRow(out, "total", totals, matchClassNames, fraction);
    }

    void SnpDiff::printConcordanceRow(std::ostream& out,
                                      const char* rowName,
                                      const std::map<std::string,int>& matchTable,
                                      bool fraction,
                                      bool alleleConcordance) const
    {
        out << rowName;
        uint32_t total = 0;
        uint32_t concordant = 0;
        uint32_t nocall = 0;
        uint32_t discordant = 0;
        typedef pair<string,int> McT;
        BOOST_FOREACH(const McT& matchCount, matchTable)
        {
            if (alleleConcordance)
            {
                // Allele concordance table.
                BOOST_FOREACH(char ch, matchCount.first)
                {
                    total += matchCount.second;
                    if ('.' == ch || 'D' == ch)
                        discordant += matchCount.second;
                    else if ('N' == ch)
                        nocall += matchCount.second;
                    else
                        concordant += matchCount.second;
                }
            }
            else
            {
                // Locus concordance table.
                total += matchCount.second;
                if (string::npos != matchCount.first.find('.') ||
                    string::npos != matchCount.first.find('D'))
                    discordant += matchCount.second;
                else if (string::npos != matchCount.first.find('N'))
                    nocall += matchCount.second;
                else
                    concordant += matchCount.second;
            }
        }
        if (fraction)
        {
            out << "\t" << formatCount(discordant, discordant+concordant, fraction);
            out << "\t" << formatCount(nocall, total, fraction);
        }
        else
        {
            out << "\t" << formatCount(concordant, total, fraction);
            out << "\t" << formatCount(discordant, total, fraction);
            out << "\t" << formatCount(nocall, total, fraction);
            out << "\t" << formatCount(total, total, fraction);
        }
        out << "\n";
    }

    void SnpDiff::printConcordanceTable(std::ostream& out,
                                        const std::vector< std::map<std::string,int> >& stats,
                                        bool fraction,
                                        bool alleleConcordance) const
    {
        // Print header.
        if (alleleConcordance)
            out << "Allele concordance ";
        else
            out << "Locus concordance ";
        if (fraction)
            out << "by fraction";
        else
            out << "by count";
        out << "\nGenotype";
        if (fraction)
            out << "\tdiscordance\tnocall\n";
        else
            out << "\tconcordance\tdiscordance\tnocall\ttotal\n";

        // Print table.
        map<string,int> totals;
        for(int ii=0; ii<GENOTYPE_CLASS_OTHER; ii++)
        {
            typedef pair<string,int> McT;
            BOOST_FOREACH(const McT& matchCount, stats[ii])
            {
                totals[matchCount.first] += matchCount.second;
            }

            printConcordanceRow(out, genotypeClassName[ii], stats[ii], fraction, alleleConcordance);
        }
        printConcordanceRow(out, "total", totals, fraction, alleleConcordance);
    }

} } // cgatools::command
