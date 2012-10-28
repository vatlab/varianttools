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

#include <boost/array.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>

#include "cgatools/command/VarFileCombine.hpp"

#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"
#include "cgatools/reference/RepeatMaskerStore.hpp"
#include "cgatools/variants/VariantFileIterator.hpp"
#include "cgatools/cgdata/CnvSegmentStore.hpp"
#include "cgatools/cgdata/CnvDetailStore.hpp"
#include "cgatools/cgdata/EvidenceReader.hpp"
#include "cgatools/cgdata/ReferenceSupportReader.hpp"

using std::string;
using std::vector;
using std::set;
using std::map;
using boost::shared_ptr;
using namespace cgatools::variants;
using namespace cgatools::reference;
using namespace cgatools::cgdata;
using namespace cgatools::util;

namespace ba = boost::algorithm;
namespace fs = boost::filesystem;

using namespace std;

namespace cgatools { namespace command {

namespace {

const char SEP = '\t';
const string ALL_ANNOTATIONS = "copy,evidence,ref,gene,ncrna,repeat,segdup,"
                               "cnv,cnvDiploid,cnvNondiploid,cnvSomNondiploid";
const string DEFAULT_ANNOTATIONS = "copy,evidence,gene,ncrna,repeat,segdup,cnv";

void checkReferenceMatch(const string& src,
                         const DelimitedFileMetadata& srcMeta,
                         const DelimitedFileMetadata& varMeta)
{
    const string REFKEY = "GENOME_REFERENCE";
    const string OLD_BUILD = "NCBI build 36";

    const string& storeref = srcMeta.get(REFKEY);
    string varref;
    if (varMeta.hasKey(REFKEY))
        varref = varMeta.get(REFKEY);
    else
        varref = OLD_BUILD;

    if (varref != storeref)
    {
        throw Exception(src + " reference build '" + storeref +
                        "' doesn't match variation file reference '" +
                        varref + "'");
    }
}

void checkAsmId(const std::string& data,
                 const DelimitedFileMetadata& dataMeta,
                 const DelimitedFileMetadata& varMeta)
{
    const string KEY = "ASSEMBLY_ID";

    // Old files didn't have assembly ID, so there's nothing I can check
    if (!varMeta.hasKey(KEY) && !dataMeta.hasKey(KEY))
        return;

    if (varMeta.get(KEY) != dataMeta.get(KEY))
    {
        throw Exception("variation file assembly ID '" + varMeta.get(KEY) +
                        "' doesn't match " + data + " assembly ID '" +
                        dataMeta.get(KEY) + "'");
    }
}

} // local namespace

class AnnotationSource
{
public:
    AnnotationSource(const CrrFile& crr)
        :   crr_(crr)
    {}

    virtual void addMeta(DelimitedFileMetadata& meta)
    {}

    virtual void addHeaders(std::ostream& out) = 0;
    virtual void addColumns(std::ostream& out, const Locus& loc) = 0;
    virtual ~AnnotationSource() {}

protected:
    const CrrFile& crr_;

    void addEmpty(std::ostream& out, size_t columnCount)
    {
        for (size_t ii = 0; ii < columnCount; ++ii)
            out << SEP;
    }
};

class EvidenceAnnotation : public AnnotationSource
{
private:
    static const size_t COLUMN_COUNT = 5;
public:
    EvidenceAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr), reader_(crr, exp)
    {}

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "evidenceIntervalId"
            << SEP << "allele1ReadCount"
            << SEP << "allele2ReadCount"
            << SEP << "referenceAlleleReadCount"
            << SEP << "totalReadCount";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        if (loc.isNoCallLocus() || loc.isRefCallLocus())
        {
            addEmpty(out, COLUMN_COUNT);
            return;
        }
        reader_.seek(loc.getRange());
        if (!reader_.inInterval())
        {
            addEmpty(out, COLUMN_COUNT);
            return;
        }
        const EvidenceReader::IntervalRecord& interval = reader_.getInterval();
        boost::array<int32_t, 4> readCounts;
        readCounts.assign(-1);
        loc.computeReadCounts(&readCounts[0],
                              &readCounts[1],
                              &readCounts[2],
                              &readCounts[3],
                              reader_);
        out << SEP << interval.intervalId_;
        for(size_t ii=0; ii<readCounts.size(); ii++)
        {
            out << SEP;
            if (-1 != readCounts[ii])
                out << readCounts[ii];
        }
    }

private:
    EvidenceReader reader_;
};

class RefSupportAnnotation : public AnnotationSource
{
public:
    RefSupportAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr), reader_(crr, exp)
    {}

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "minReferenceScore";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        static const int32_t MISSING = std::numeric_limits<int32_t>::min();
        Range xrg = loc.getRange();
        if (xrg.begin_ > 0)
            --xrg.begin_;
        ++xrg.end_;
        reader_.seek(xrg);
        int32_t score = reader_.getMinScore(xrg, MISSING);
        out << SEP;
        if (score != MISSING)
            out << score;
    }
private:
    ReferenceSupportReader reader_;
};

class CopyAnnotation : public AnnotationSource
{
public:
    CopyAnnotation(const CrrFile& crr, const VariantFileIterator& srcFile)
        : AnnotationSource(crr)
    {
        headers_ = srcFile.getAnnotationColumnHeaders();
        srcMeta_ = srcFile.getMetadata();
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        meta.transfer(srcMeta_,
                      "MIRBASE_VERSION,GENE_ANNOTATIONS,PFAM_DATE,"
                      "REPMASK_GENERATED_AT,SEGDUP_GENERATED_AT,"
                      "DGV_VERSION,CNV_WINDOW_WIDTH,"
                      "CNV_DIPLOID_WINDOW_WIDTH,CNV_NONDIPLOID_WINDOW_WIDTH,"
                      "CNV_SOMATIC_NONDIPLOID_WINDOW_WIDTH");
    }

    virtual void addHeaders(std::ostream& out)
    {
        if (headers_.size() > 0)
            out << SEP << ba::join(headers_, string(1, SEP));
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        size_t cols = 0;
        BOOST_FOREACH(const string& ex, loc.extras_)
        {
            out << SEP << ex;
            ++cols;
        }
        CGA_ASSERT(cols == headers_.size());
    }
private:
    vector<string> headers_;
    DelimitedFileMetadata srcMeta_;
};

class GeneAnnotation : public AnnotationSource
{
public:
    GeneAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr)
    {
        loadAnnotations(exp.getGeneFileName());
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkAsmId("gene", srcMeta_, meta);
        meta.transfer(srcMeta_, "GENE_ANNOTATIONS,PFAM_DATE");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "allele1Gene"
            << SEP << "allele2Gene"
            << SEP << "pfam";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        for (size_t ii = 0; ii < 2; ++ii)
        {
            out << SEP;
            if (ii >= loc.getPloidy())
                continue;
            uint32_t alleleIdx = loc.getAlleles()[ii].getCall(0).haplotype_;
            if (alleleIdx == Call::ALL_HAPLOTYPES)
                alleleIdx = ii + 1;
            AlleleId aid(loc.getId(), alleleIdx);
            Annotations::const_iterator it = annotations_.find(aid);
            if (it != annotations_.end())
                out << ba::join(it->second, ";");
        }
        out << SEP;
        Pfams::const_iterator it = pfams_.find(loc.getId());
        if (it != pfams_.end())
            out << ba::join(it->second, ";");
    }

private:
    typedef std::pair< uint32_t, uint32_t > AlleleId;
    typedef map< AlleleId, set<string> > Annotations;
    Annotations annotations_;
    typedef map< uint32_t, set<string> > Pfams;
    Pfams pfams_;
    DelimitedFileMetadata srcMeta_;

    void loadAnnotations(const string& fn)
    {
        boost::shared_ptr<std::istream> in = InputStream::openCompressedInputStreamByExtension(fn);
        DelimitedFile df(*in, fn);
        srcMeta_ = df.getMetadata();

        AlleleId aid;
        vector<string> tokens;
        tokens.resize(5);
        string pfam;

        df.addField(ValueField<uint32_t>("locus", &aid.first));
        if (df.hasField("allele"))
            df.addField(ValueField<uint32_t>("allele", &aid.second));
        else
            df.addField(ValueField<uint32_t>("haplotype", &aid.second));
        df.addField(StringField("geneId", &tokens[0]));
        df.addField(StringField("mrnaAcc", &tokens[1]));
        df.addField(StringField("symbol", &tokens[2]), DelimitedFile::FPT_OPTIONAL);
        df.addField(StringField("pfam", &pfam), DelimitedFile::FPT_OPTIONAL);
        if (df.hasField("component"))
            df.addField(StringField("component", &tokens[3]));
        else
            df.addField(StringField("exonCategory", &tokens[3]));
        if (df.hasField("impact"))
            df.addField(StringField("impact", &tokens[4]));
        else
            df.addField(StringField("aaCategory", &tokens[4]));

        while (df.next())
        {
            if (!pfam.empty())
            {
                ba::split(pfams_[aid.first], pfam, ba::is_any_of(";"), ba::token_compress_on);
            }

            bool empty = true;
            BOOST_FOREACH(const string& tok, tokens)
            {
                if (!tok.empty())
                {
                    empty = false;
                    break;
                }
            }
            if (empty)
                continue;

            string record = ba::join(tokens, ":");
            annotations_[aid].insert(record);
        }
    }
};

class MirnaAnnotation : public AnnotationSource
{
public:
    MirnaAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr)
    {
        loadAnnotations(exp.getNcRNAFileName());
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkAsmId("miRNA", srcMeta_, meta);
        meta.transfer(srcMeta_, "MIRBASE_VERSION");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "miRBaseId";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        out << SEP;
        Annotations::const_iterator it = annotations_.find(loc.getId());
        if (it != annotations_.end())
            out << ba::join(it->second, ";");
    }

protected:
    typedef map< uint32_t, set<string> > Annotations;
    Annotations annotations_;
    DelimitedFileMetadata srcMeta_;

    void loadAnnotations(const string& fn)
    {
        boost::shared_ptr<std::istream> in = InputStream::openCompressedInputStreamByExtension(fn);
        DelimitedFile df(*in, fn);
        srcMeta_ = df.getMetadata();

        uint32_t locusId;
        string id;

        df.addField(ValueField<uint32_t>("locus", &locusId));
        df.addField(StringField("miRBaseId", &id));

        while (df.next())
        {
            ba::split(annotations_[locusId], id, ba::is_any_of(";"), ba::token_compress_on);
        }
    }
};

class CnvDiploidAnnotation : public AnnotationSource
{
public:
    CnvDiploidAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr), store_(crr, exp, true)
    {
        windowWidth_ = parseValue<uint32_t>(store_.getMetadata().get("WINDOW_WIDTH"));
        CGA_ASSERT(store_.hasCalledPloidy());
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkAsmId("CNV", store_.getMetadata(), meta);
        meta.transfer(store_.getMetadata(), "DGV_VERSION");
        meta.transfer(store_.getMetadata(), "WINDOW_WIDTH", "CNV_DIPLOID_");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "relativeCoverageDiploid";
        out << SEP << "calledPloidy";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        if (loc.getRange().length() >= windowWidth_ * 2)
        {
            addEmpty(out, 2);
            return;
        }

        const CnvSegmentData* seg =
                store_.getBestOverlappingSegment(loc.getRange());
        if (0 == seg)
        {
            addEmpty(out, 2);
            return;
        }

        if (seg->relativeCoverage_ >= .0)
            out << SEP << boost::format("%.2f") % seg->relativeCoverage_;
        else
            out << SEP << "N";

        out << SEP << seg->calledPloidy_;
    }
private:
    CnvSegmentStore store_;
    uint32_t windowWidth_;
};

class CnvNondiploidAnnotation : public AnnotationSource
{
public:
    CnvNondiploidAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr), store_(crr, exp, false)
    {
        windowWidth_ = parseValue<uint32_t>(store_.getMetadata().get("WINDOW_WIDTH"));
        CGA_ASSERT(store_.hasCalledLevel());
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkAsmId("CNV", store_.getMetadata(), meta);
        meta.transfer(store_.getMetadata(), "DGV_VERSION");
        meta.transfer(store_.getMetadata(), "WINDOW_WIDTH", "CNV_NONDIPLOID_");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "relativeCoverageNondiploid";
        out << SEP << "calledLevel";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        if (loc.getRange().length() >= windowWidth_ * 2)
        {
            addEmpty(out, 2);
            return;
        }

        const CnvSegmentData* seg =
                store_.getBestOverlappingSegment(loc.getRange());
        if (0 == seg)
        {
            addEmpty(out, 2);
            return;
        }

        if (seg->relativeCoverage_ >= .0)
            out << SEP << boost::format("%.2f") % seg->relativeCoverage_;
        else
            out << SEP << "N";

        out << SEP << seg->calledLevel_;
    }
private:
    CnvSegmentStore store_;
    uint32_t windowWidth_;
};

class CnvSomaticNondiploidAnnotation : public AnnotationSource
{
public:
    CnvSomaticNondiploidAnnotation(const CrrFile& crr, const GenomeMetadata& exp)
        : AnnotationSource(crr), store_(crr, exp, false,true)
    {
        windowWidth_ = parseValue<uint32_t>(store_.getMetadata().get("WINDOW_WIDTH"));
        CGA_ASSERT(store_.hasCalledLevel());
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkAsmId("CNV", store_.getMetadata(), meta);
        meta.transfer(store_.getMetadata(), "DGV_VERSION");
        meta.transfer(store_.getMetadata(), "WINDOW_WIDTH", "CNV_SOMATIC_NONDIPLOID_");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "relativeCoverageSomaticNondiploid";
        out << SEP << "somaticCalledLevel";
        out << SEP << "bestLAF";
        out << SEP << "lowLAF";
        out << SEP << "highLAF";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        if (loc.getRange().length() >= windowWidth_ * 2)
        {
            addEmpty(out, 5);
            return;
        }

        const CnvDetailData* det =
                store_.getBestOverlappingDetail(loc.getRange());

        if (0 == det)
        {
            addEmpty(out, 5);
            return;
        }

        if (det->relativeCoverage_ >= .0)
            out << SEP << boost::format("%.2f") % det->relativeCoverage_;
        else
            out << SEP << "N";

        out << SEP << det->calledLevel_;
        out << SEP << det->bestLAF_;
        out << SEP << det->lowLAF_;
        out << SEP << det->highLAF_;
    }
private:
    CnvDetailStore store_;
    uint32_t windowWidth_;
};

class RepMaskAnnotation : public AnnotationSource
{
public:
    RepMaskAnnotation(const CrrFile& crr, const string& fn)
        : AnnotationSource(crr), store_(crr, fn)
    {
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkReferenceMatch("repeat masker", store_.getMetadata(), meta);
        meta.transfer(store_.getMetadata(), "GENERATED_AT", "REPMASK_");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "repeatMasker";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        out << SEP;

        if (loc.isRefCallLocus() || loc.isNoCallLocus())
            return;

        std::vector<RepeatMaskerStore::QueryResultType> buffer;
        store_.intersect(loc.getRange(), buffer);

        set<string> annotations;
        BOOST_FOREACH(RepeatMaskerStore::QueryResultType r, buffer)
        {
            const RepeatMaskerAnnotation& s = r->second;
            string ann = (boost::format("%s:%s:%.1f") %
                             s.name_ % s.family_ % s.divergence_).str();
            annotations.insert(ann);
        }

        out << ba::join(annotations, ";");
    }
private:
    RepeatMaskerStore store_;
};

struct SegDupData
{
    uint32_t overlap_;
};

class SegDupStore : public RangeAnnotationStore<SegDupStore, SegDupData>
{
public:
    SegDupStore(const reference::CrrFile& crr, const std::string& fn)
        : Base(crr)
    {
        load(fn);
    }

    void bindColumns(util::DelimitedFile& df,
                     reference::Range& range,
                     SegDupData& data)
    {
        bindRangeColumns(df, range);
        df.addField(ValueField<uint32_t>("count", &data.overlap_));
    }
};

class SegDupAnnotation : public AnnotationSource
{
public:
    SegDupAnnotation(const CrrFile& crr, const string& fn)
        : AnnotationSource(crr), store_(crr, fn)
    {
    }

    virtual void addMeta(DelimitedFileMetadata& meta)
    {
        checkReferenceMatch("segdup", store_.getMetadata(), meta);
        meta.transfer(store_.getMetadata(), "GENERATED_AT", "SEGDUP_");
    }

    virtual void addHeaders(std::ostream& out)
    {
        out << SEP << "segDupOverlap";
    }

    virtual void addColumns(std::ostream& out, const Locus& loc)
    {
        out << SEP;

        if (loc.isRefCallLocus() || loc.isNoCallLocus())
            return;

        std::vector<SegDupStore::QueryResultType> buffer;
        store_.intersect(loc.getRange(), buffer);

        uint32_t overlap = 0;
        BOOST_FOREACH(SegDupStore::QueryResultType r, buffer)
        {
            overlap = std::max(overlap, r->second.overlap_);
        }

        if (0 != overlap)
            out << overlap;
    }
private:
    SegDupStore store_;
};

VarFileCombine::VarFileCombine(const std::string& name)
    : Command(name,
              "Converts a variation file to a one-line-per-locus format.",
              "0.3 or later",

              "The output file contains one line for each locus in the input variation "
              "file. The following columns are always present:\n\n"
              "    locus           \tLocus ID, as in the input file.\n"
              "    chromosome      \tThe name of the chromosome.\n"
              "    begin           \tThe first base of the locus interval, 0-based.\n"
              "    end             \tThe first base past the locus interval, 0-based.\n"
              "    zygosity        \tOne of the following values:\n"
              "        no-call         \tboth alleles contain no-calls\n"
              "        half            \tone allele fully called\n"
              "        hap             \thaploid region\n"
              "        hom             \thomozygous region\n"
              "        het-ref         \theterozygous region, one allele is reference\n"
              "        het-alt         \theterozygous region, neither allele is reference\n"
              "    varType         \tFor simple loci, one of \"ref\", \"snp\", "
              "\"del\", \"ins\" or \"sub\". For more complex regions, \"complex\".\n"
              "    reference       \tReference sequence, or \"=\" for pure reference or "
              "pure no call regions.\n"
              "    allele1Seq      \tSequence of the first allele, may contain \"?\" or "
              "\"N\" characters for unknown-length and known-length no-calls, respectively.\n"
              "    allele2Seq      \tSequence of the second allele.\n"
              "    allele1VarScoreVAF  \tThe varScoreVAF of the first allele. For pre-2.0 "
              "var files, which have totalScore instead of varScoreVAF, this column is "
              "filled in with totalScore. For the loci that contain "
              "multiple calls, this is the minimum score across all calls.\n"
              "    allele2VarScoreVAF  \tThe varScoreVAF of the first allele. For pre-2.0 "
              "var files, which have totalScore instead of varScoreVAF, this column is "
              "filled in with totalScore. For the loci that contain "
              "multiple calls, this is the minimum score across all calls.\n"
              "    allele1VarScoreEAF  \tThe varScoreEAF of the first allele. For pre-2.0 "
              "var files, which have totalScore instead of varScoreEAF, this column is "
              "filled in with totalScore. For the loci that contain "
              "multiple calls, this is the minimum score across all calls.\n"
              "    allele2VarScoreEAF  \tThe varScoreEAF of the first allele. For pre-2.0 "
              "var files, which have totalScore instead of varScoreEAF, this column is "
              "filled in with totalScore. For the loci that contain "
              "multiple calls, this is the minimum score across all calls.\n"
              "    allele1VarQuality   \tThe varQuality of the first allele. For pre-2.0 "
              "var files, which do not have a varQuality column, this field is empty. For "
              "multiple calls, this is the lowest quality across all calls (and empty is "
              "considered lower quality than VQLOW.\n"
              "    allele2VarQuality   \tThe varQuality of the first allele. For pre-2.0 "
              "var files, which do not have a varQuality column, this field is empty. For "
              "multiple calls, this is the lowest quality across all calls (and empty is "
              "considered lower quality than VQLOW.\n"
              "    allele1HapLink  \tHaplink ID of the first allele. Alleles with the same "
              "ID are known to reside on the same haplotype.\n"
              "    allele2HapLink  \tHaplink ID of the second allele.\n"
              "\n"
              "The allele to be placed first is chosen according to the following "
              "priority list: \n"
              "\n"
              "    fully called variant allele;\n"
              "    fully called reference allele;\n"
              "    partially called allele;\n"
              "    completely no-called allele.\n"
              "\n"
              "In addition to the mandatory columns above, various annotation columns "
              "can be added to the file using \"annotations\" parameter. The supported "
              "annotation sources and the corresponding additional columns are listed "
              "below.\n"
              "\n"
              "    copy             \tAdds column \"xRef\" that contains a concatenation "
              "of all dbSNP annotations for this locus from the input variant file. If "
              "the source file is already in one-line-per-locus format, also copies over "
              "all other annotations already present in the source.\n"
              "    evidence         \tAdds columns:\n"
              "        evidenceIntervalId   \tID of the corresponding evidence interval.\n"
              "        allele1ReadCount     \tNumber of evidence reads that support the "
              "first allele\n"
              "        allele2ReadCount     \tNumber of evidence reads that support the "
              "second allele\n"
              "        referenceAlleleReadCount \tNumber of evidence reads that support the "
              "reference\n"
              "        totalReadCount       \tTotal number of evidence reads that overlap "
              "the locus. This includes reads that don't show strong support for either "
              "of the called alleles.\n"
              "    ref              \tAdds column \"minReferenceScore\" that contains the "
              "minimum value of the reference score over the locus interval extended by one "
              "base in either direction. Off by default.\n"
              "    gene             \tAdds columns \"allele1Gene\" and \"allele2Gene\" that "
              "contain summarized information about the overlap and impact on known genes. "
              "Derived from the gene annotation in the CGI data package.\n"
              "    ncrna            \tAdds column \"miRBaseId\" that contains summarized "
              "information about overlap with non-coding RNA. Derived from the ncRNA file in "
              "the CGI data package.\n"
              "    repeat           \tAdds column \"repeatMasker\" that contains information "
              "about the repeats overlapping the locus. Requires a data file available from "
              "the Complete Genomics site: ftp://ftp.completegenomics.com/AnnotationFiles/\n"
              "    segdup           \tAdds column \"segDupOverlap\" that contains the number "
              "segmental duplications overlapping the locus. Requires a data file available "
              "from Complete Genomics site: ftp://ftp.completegenomics.com/AnnotationFiles/\n"
              "    cnv              \tAdds the cnvDiploid, cnvNondiploid and (if present) "
              "cnvSomNondiploid annotations (described below), as available in the CGI data package.\n"
              "    cnvDiploid       \tAdds columns \"relativeCoverageDiploid\" and \"calledPloidy\" "
              "derived from the diploid CNV calls in the CGI data package.\n"
              "    cnvNondiploid    \tAdds columns \"relativeCoverageNondiploid\" and \"calledLevel\" "
              "derived from the nondiploid CNV calls in the CGI data package.\n"
              "    cnvSomNondiploid \tAdds columns \"relativeCoverageSomaticNondiploid\" "
              "and \"somaticCalledLevel\" derived from the somatic nondiploid CNV call details "
              "in the CGI data package.\n"
              "\n"
              "The column groups are added in the order they are listed in the "
              "\"annotations\" command line parameter. By default the tool will attempt "
              "to add all annotations. For older data packages that do not contain some "
              "of the necessary files remove the corresponding annotation source from "
              "the list."
        )
{
    options_.add_options()
        ("reference", po::value<string>(&referenceFileName_),
         "The reference crr file.")
        ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
         "The output file (may be omitted for stdout).")
        ("variants", po::value<string>(&variantFileName_),
         "The input variant file.")
        ("annotations", po::value<string>(&annotations_)->default_value(DEFAULT_ANNOTATIONS),
         "Comma-separated list of annotations to add to each line.")
        ("genome-root", po::value<string>(&exportRoot_),
         "The genome directory, for example /data/GS00118-DNA_A01; "
         "this directory is expected to contain an intact ASM subdirectory.")
        ("repmask-data", po::value<string>(&repmaskFileName_),
         "The file that contains repeat masker data.")
        ("segdup-data", po::value<string>(&segdupFileName_),
         "The file that contains segdup data.")
        ;
}

void VarFileCombine::addAnnotators(const VariantFileIterator& srcFile)
{
    StringSet expset("evidence,gene,ncrna,cnv,cnvDiploid,cnvNondiploid,cnvSomNondiploid,ref", 
                     ALL_ANNOTATIONS, "");

    // Create a set of annotations to check the names
    StringSet annset(annotations_, ALL_ANNOTATIONS, "not a valid annotation name");

    // But also need a list of annotations since the order matters
    vector<string> annlist;
    ba::split(annlist, annotations_, ba::is_any_of(","));

    if (annlist.size() != annset.size())
        throw util::Exception("annotation list contains duplicates: " + annotations_);

    // Do all the checks first to error out early, loading annotations
    // may take a while
    BOOST_FOREACH(const string& s, annlist)
    {
        if (0 != expset.count(s) && !exp_)
            throw util::Exception("annotation requires genome-root to be set: " + s);

        if (s == "repeat")
        {
            if (repmaskFileName_.empty())
                throw util::Exception("annotation requires repmask-data to be set: " + s);
        }
        else if (s == "segdup")
        {
            if (segdupFileName_.empty())
                throw util::Exception("annotation requires segdup-data to be set: " + s);
        }
        else if (s == "ncrna" && !InputStream::isReadable(exp_->getNcRNAFileName()))
        {
            throw util::Exception("file required for 'ncrna' annotation not found: " +
                                        exp_->getNcRNAFileName());
        }
        else if (s == "cnvDiploid" && "" == exp_->getCnvSegmentsDiploidFileName(false))
        {
            cerr << "file required for 'cnvDiploid' annotation not found" << endl;
            exp_->getCnvSegmentsDiploidFileName();
        }
        else if (s == "cnvNondiploid" && "" == exp_->getCnvSegmentsNondiploidFileName(false))
        {
            cerr << "file required for 'cnvNondiploid' annotation not found" << endl;
            exp_->getCnvSegmentsNondiploidFileName();
        }
        else if (s == "cnvSomNondiploid" && "" == exp_->getCnvDetailsSomaticNondiploidFileName(false))
        {
            cerr << "file required for 'cnvSomNondiploid' annotation not found" << endl;
            exp_->getCnvDetailsSomaticNondiploidFileName();
        }
        else if (s == "cnv" &&
                 "" == exp_->getCnvSegmentsDiploidFileName(false) &&
                 "" == exp_->getCnvSegmentsNondiploidFileName(false))
        {
            cerr << "one or more files required for 'cnv' annotation not found" << endl;
            exp_->getCnvSegmentsDiploidFileName();
            exp_->getCnvSegmentsNondiploidFileName();
        }
    }

    BOOST_FOREACH(const string& s, annlist)
    {
        shared_ptr<AnnotationSource> p;

        if (s == "evidence")
            p.reset(new EvidenceAnnotation(crr_, *exp_));
        else if (s == "copy")
            p.reset(new CopyAnnotation(crr_, srcFile));
        else if (s == "ref")
            p.reset(new RefSupportAnnotation(crr_, *exp_));
        else if (s == "gene")
            p.reset(new GeneAnnotation(crr_, *exp_));
        else if (s == "ncrna")
            p.reset(new MirnaAnnotation(crr_, *exp_));
        else if (s == "cnv" || s == "cnvDiploid" || s == "cnvNondiploid" || s == "cnvSomNondiploid" )
        {
            if (s == "cnvDiploid" ||
                (s == "cnv" && "" != exp_->getCnvSegmentsDiploidFileName(false)))
            {
                p.reset(new CnvDiploidAnnotation(crr_, *exp_));
                annotators_.push_back(p);
            }
            if (s == "cnvNondiploid" ||
                (s == "cnv" && "" != exp_->getCnvSegmentsNondiploidFileName(false)))
            {
                p.reset(new CnvNondiploidAnnotation(crr_, *exp_));
                annotators_.push_back(p);
            }
            if (s == "cnvSomNondiploid" ||
                (s == "cnv" && "" != exp_->getCnvDetailsSomaticNondiploidFileName(false)))
            {
                p.reset(new CnvSomaticNondiploidAnnotation(crr_, *exp_));
                annotators_.push_back(p);
            }
            continue;
        }
        else if (s == "repeat")
            p.reset(new RepMaskAnnotation(crr_, repmaskFileName_));
        else if (s == "segdup")
            p.reset(new SegDupAnnotation(crr_, segdupFileName_));
        else
            CGA_ASSERT(false);
        annotators_.push_back(p);
    }
}

void VarFileCombine::fillMetadata(DelimitedFileMetadata& meta, const DelimitedFileMetadata& vfMeta)
{
    meta.initDefaults();

    meta.transfer(vfMeta, "ASSEMBLY_ID,COSMIC,DBSNP_BUILD,GENOME_REFERENCE,"
                          "SAMPLE,SOFTWARE_VERSION,SRC_SOFTWARE_VERSION");

    meta.set("TYPE", "VAR-OLPL");

    BOOST_FOREACH(const shared_ptr<AnnotationSource>& asrc, annotators_)
    {
        asrc->addMeta(meta);
    }
    meta.sort();
}

int VarFileCombine::run(po::variables_map& vm)
{
    requireParam(vm, "reference");
    requireParam(vm, "variants");

    crr_.open(referenceFileName_);
    if (!exportRoot_.empty())
        exp_.reset(new GenomeMetadata(exportRoot_));

    VariantFileIterator locIt(crr_);
    locIt.setReferenceCoverValidation(false);

    locIt.open(variantFileName_);

    if (exp_ && locIt.getMetadata().hasKey("ASSEMBLY_ID"))
    {
        const std::string& varAsmId = locIt.getMetadata().get("ASSEMBLY_ID");
        if (exp_->getAsmId() != varAsmId)
        {
            throw Exception("variation file assembly ID '" + varAsmId +
                            "' doesn't match data package assembly ID '" +
                            exp_->getAsmId() + "'");
        }
    }

    addAnnotators(locIt);

    DelimitedFileMetadata meta;
    fillMetadata(meta, locIt.getMetadata());

    std::ostream& out = openStdout(outputFileName_);

    out << meta;
    Locus::writeOneLineFileHeader(out, locIt.xRefIsAlleleSpecific(), SEP);
    BOOST_FOREACH(const shared_ptr<AnnotationSource>& asrc, annotators_)
    {
        asrc->addHeaders(out);
    }
    out << '\n';

    Locus loc;
    for (; !locIt.eof(); ++locIt)
    {
        loc = *locIt;
        loc.reorderAlleles();
        loc.writeAsOneLine(out, false, SEP);
        BOOST_FOREACH(const shared_ptr<AnnotationSource>& asrc, annotators_)
        {
            asrc->addColumns(out, loc);
        }
        out << '\n';
    }

    return 0;
}

}}
