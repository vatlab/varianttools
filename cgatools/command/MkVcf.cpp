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
#include "cgatools/command/MkVcf.hpp"
#include "cgatools/variants/VariantFileVcfSource.hpp"
#include "cgatools/junctions/JunctionVcfRecord.hpp"
#include "cgatools/copynumber/CnvFileVcfSource.hpp"
#include "cgatools/mobileelement/MeiFileVcfSource.hpp"
#include "cgatools/util/parse.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

using namespace std;
using namespace cgatools::util;
using namespace cgatools::variants;
using namespace cgatools::copynumber;
using namespace cgatools::mobileelement;
using cgatools::conv::VcfRecordWriter;
using cgatools::conv::VcfRecordSource;
using cgatools::conv::VcfKvHeaderRecord;
using cgatools::conv::VcfSubFieldHeaderRecord;
using cgatools::cgdata::GenomeMetadata;
using cgatools::variants::calib::CalibratedScorer;
using namespace cgatools::reference;
using boost::array;
using boost::lexical_cast;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

namespace cgatools { namespace command {

    bool MkVcf::PQItemComparator::operator()(const PQItem& i1, const PQItem& i2) const
    {
        bool eof1 = i1->eof();
        bool eof2 = i2->eof();

        if (eof1 || eof2)
            return eof1;

        Location loc1 = (*i1)->getLocation();
        Location loc2 = (*i2)->getLocation();
        return loc1 > loc2;
    }

    MkVcf::MkVcf(const std::string& name)
        : Command(name,
                  "Converts var file(s) or masterVar file(s) to VCF.",
                  "0.3 or later",

                  "Converts var file(s) or masterVar file(s) to VCF."
            )
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
             "The reference crr file.")
            ("output", po::value<string>(&outputFileName_)->default_value("STDOUT"),
             "The output file (may be omitted for stdout).")
            ("field-names", po::value<string>(&fieldNames_)->default_value(
                 // variant fields
                 "GT,PS,NS,AN,AC,SS,FT,CGA_XR,CGA_FI,GQ,HQ,EHQ,CGA_CEHQ,GL,CGA_CEGL,"
                 "DP,AD,CGA_RDP,CGA_ODP,CGA_OAD,CGA_ORDP,"
                 "CGA_PFAM,CGA_MIRB,CGA_RPT,CGA_SDO,CGA_SOMC,CGA_SOMR,CGA_SOMS,"

                 // single-genome CNV fields
                 "GT,CGA_GP,CGA_NP,CGA_CP,CGA_PS,CGA_CT,CGA_TS,CGA_CL,CGA_LS,"

                 // somatic CNV fields
                 "CGA_SCL,CGA_SLS,CGA_LAF,CGA_LLAF,CGA_ULAF,"

                 // MEI 
                 "GT,FT,CGA_IS,CGA_IDC,CGA_IDCL,CGA_IDCR,CGA_RDC,CGA_NBET,CGA_ETS,CGA_KES,"

                 // SV
                 "GT,FT,CGA_BF,CGA_MEDEL,MATEID,SVTYPE,CGA_BNDG,CGA_BNDGO,CGA_BNDMPC,"
                 "CGA_BNDPOS,CGA_BNDDEF,CGA_BNDP"
               ),
             "Comma-separated list of field names. By default, all fields are included, "
             "but you may override this option to ensure only a subset of the fields is "
             "included in the VCF output. For a description of each field, see the "
             "cgatools user guide.")
            ("source-names", po::value<string>(&sourceNames_)->default_value("masterVar,CNV,SV,MEI"),
             "Comma-separated list of source names. The following source names are available:\n"
             "  masterVar - \tIncludes records extracted from the masterVar file.\n"
             "  CNV       - \tIncludes CNV-related records.\n"
             "  SV        - \tIncludes records derived from junctions files.\n"
             "  MEI       - \tIncludes records describing mobile element insertions.\n"
             "Some of these source types are only available for more recent pipeline "
             "versions, and some of these source types do not support multi-genome VCFs. "
             "For more information about which source types are available for which "
             "versions of the Complete Genomics pipeline software, see the cgatools "
             "user guide.")
            ("genome-root", po::value< vector<string> >(&genomeRoots_),
             "For each genome to include in the VCF, the genome root directory, for "
             "example /data/GS00118-DNA_A01; this directory is expected to contain "
             "the ASM and LIB subdirectories, for example. You must supply "
             "this option for each genome in the VCF, unless you are using "
             "--source-names=masterVar and you have specified the --master-var option "
             "for each genome in the VCF.")
            ("master-var", po::value< vector<string> >(&varNames_),
             "For each genome to include in the VCF, the masterVar file. If "
             "genome-roots parameter is given, this parameter defaults to the masterVar "
             "in the given genome-root.")
            ("include-no-calls", po::bool_switch(&includeNoCalls_)->default_value(false),
             "Small variants VCF records include loci that have no reference-inconsistent calls.\n")
            ("calibration-root", po::value<string>(&calibPrefix_),
             "The directory containing calibration data. For example, there should "
             "exist a file calibration-root/0.0.0/metrics.tsv. This option is only required "
             "if CGA_CEHQ or CGA_CEGL are included in the --field-names parameter.")
            ("junction-file", po::value< vector<string> >(&junctionFileNames_),
             "For each genome to include in the VCF, the junctions file. If genome-roots "
             "parameter is given, this parameter defaults to the respective junctions file "
             "in the export directory.")
            ("junction-score-threshold", 
             po::value<size_t>(&junctionScoreThreshold_)->default_value((size_t)10),
             "Junction score thresholds (discordant mate pair count).")
            ("junction-side-length-threshold", 
             po::value<size_t>(&junctionSideLengthThreshold_)->default_value((size_t)70),
             "Junction side length threshold.")
            ("junction-distance-tolerance", 
             po::value<size_t>(&junctionDistanceTolerance_)->default_value((size_t)200),
             "Distance tolerance for junction compatibility.")
            ("junction-length-threshold", 
             po::value<size_t>(&junctionLengthThreshold_)->default_value((size_t)500),
             "Length threshold for compatible junctions.")
            ("junction-normal-priority", 
             po::bool_switch(&junctionNormalPriority_)->default_value(false),
             "Normal junction priority for vcf output.")
             ("junction-tumor-hc", 
             po::bool_switch(&junctionTumorHC_)->default_value(false),
             "use high confidence junctions for tumors.")
            ;

        positionalOptions_.add("genome-root", -1);
    }

    void MkVcf::initSources()
    {
        vector<string> fieldNames;
        boost::split(fieldNames, fieldNames_, boost::is_any_of(","));
        boost::split(vSourceNames_, sourceNames_, boost::is_any_of(","));

        BOOST_FOREACH(const string& sourceName, vSourceNames_)
        {
            try
            {
                if ("masterVar" == sourceName)
                {
                    std::vector< boost::shared_ptr<VariantFileIterator> > var(genomeCount_);
                    if (varNames_.size() > 0)
                    {
                        for(size_t ii=0; ii<var.size(); ii++)
                        {
                            var[ii].reset(new VariantFileIterator(crr_));
                            var[ii]->open(varNames_[ii]);
                        }
                    }
                    else
                    {
                        for(size_t ii=0; ii<var.size(); ii++)
                        {
                            GenomeMetadata meta(genomeRoots_[ii]);
                            var[ii].reset(new VariantFileIterator(crr_));
                            var[ii]->open(GenomeMetadata(genomeRoots_[ii]).getMasterVarFileName());
                        }
                    }

                    if (0 == var.size())
                        throw Exception("no genome-root or master-var parameters specified");
                    shared_ptr<VcfRecordSource> source(
                        new VariantFileVcfRecordSource(
                            var, fieldNames, crr_, calibPrefix_, includeNoCalls_));
                    sources_.push(source);
                    vSources_.push_back(source);
                }
                else if ("SV" == sourceName)
                {
                    std::vector< boost::shared_ptr<cgatools::cgdata::GenomeMetadata> > genomes;
                    BOOST_FOREACH(const std::string& groot,genomeRoots_)
                    {
                        genomes.push_back
                            ( boost::shared_ptr<cgatools::cgdata::GenomeMetadata>
                              (new GenomeMetadata(groot)) 
                                );
                    }

                    shared_ptr<VcfRecordSource> source 
                        (
                            new cgatools::junctions::JunctionVcfRecordSource 
                            ( 
                                genomes, junctionFileNames_, fieldNames, crr_, 
                                junctionScoreThreshold_, junctionSideLengthThreshold_,
                                junctionDistanceTolerance_, junctionLengthThreshold_,
                                junctionNormalPriority_, junctionTumorHC_
                                )
                            );
                    
                    sources_.push(source);
                    vSources_.push_back(source);
                }
                else if ("CNV" == sourceName)
                {
                    std::vector< std::string > dipdetFn(0);
                    std::vector< std::string > nondipdetFn(0);
                    std::vector< std::string > somnondipdetFn(0);

                    for(size_t ii=0; ii<genomeRoots_.size(); ii++)
                    {
                        GenomeMetadata gmd(genomeRoots_[ii]);
                        dipdetFn.push_back(gmd.getCnvDetailsDiploidFileName());
                        nondipdetFn.push_back(gmd.getCnvDetailsNondiploidFileName());
                        somnondipdetFn.push_back(gmd.getCnvDetailsSomaticNondiploidFileName(false));
                    }

                    shared_ptr<VcfRecordSource> source(
                        new CnvFileVcfRecordSource(
                            dipdetFn,nondipdetFn,somnondipdetFn,
                            fieldNames, crr_) );
                    sources_.push(source);
                    vSources_.push_back(source);
                
                }
                else if ("MEI" == sourceName)
                {
                    std::vector< std::string > meiFn(0);

                    for(size_t ii=0; ii<genomeRoots_.size(); ii++)
                    {
                        GenomeMetadata gmd(genomeRoots_[ii]);
                        meiFn.push_back(gmd.getMobileElementInsertionFileName());
                    }
                    shared_ptr<VcfRecordSource> source
                        (
                            new MeiFileVcfRecordSource
                            ( 
                                meiFn,fieldNames, crr_
                                )
                            );
                    sources_.push(source);
                    vSources_.push_back(source);
                }
                else
                {
                    throw Exception("source name not recognized");
                }
            }
            catch(std::exception& ee)
            {
                throw Exception(sourceName+" source could not be loaded: "+ee.what());
            }
        }
    }

    void MkVcf::validateSources()
    {
        if (genomeCount_ > 1)
        {
            for(size_t ii=0; ii<vSources_.size(); ii++)
            {
                for(size_t jj=0; jj<genomeCount_; jj++)
                {
                    string softwareVersion = vSources_[ii]->getSource(jj);
                    boost::regex re("[0-9]+(\\.[0-9]+).*");
                    if (boost::regex_match(softwareVersion, re)) {
                        int majorVersion = parseValue<int>(
                            softwareVersion.substr(0, softwareVersion.find('.')));
                        if (majorVersion < 2)
                            throw Exception("SOFTWARE_VERSION="+softwareVersion+
                                            " for genome "+lexical_cast<string>(jj+1)+
                                            ", multi-genome VCF not supported for SOFTWARE_VERSION "
                                            "before 2.0.0.0");
                    }
                }
            }
        }
    }

    void MkVcf::writeHeaders(std::ostream& out)
    {
        // VCF version
        out << "##fileformat=VCFv4.1" << endl;

        // Date of creation.
        boost::posix_time::ptime tt = boost::posix_time::microsec_clock::local_time();
        string ttStr = boost::posix_time::to_iso_string(tt);
        out << "##fileDate=" << ttStr.substr(0, 8) << endl;

        // Center, required by TCGA VCF.
        out << "##center=Complete Genomics" << endl;

        transferSource(out);
        transferMetadataHeaders(out);
        transferInfoFilterHeaders(out);

        // Column headers.
        out << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";

        for(size_t ii=0; ii<genomeCount_; ii++)
            out << "\t" << getAssemblyId(ii);
        out << endl;
    }

    std::string MkVcf::getAssemblyId(size_t idxGenome) const
    {
        for(size_t ii=1; ii<vSources_.size(); ii++)
        {
            if (vSources_[0]->getAssemblyId(idxGenome) != vSources_[ii]->getAssemblyId(idxGenome))
                throw Exception("genome "+lexical_cast<string>(idxGenome+1)+
                                " ASSEMBLY_ID mismatch for sources "+
                                vSourceNames_[0]+" and "+vSourceNames_[ii]);
        }
        return vSources_[0]->getAssemblyId(idxGenome);
    }

    void MkVcf::transferSource(
        std::ostream& out) const
    {
        set<string> pipelineVer;
        for(size_t ii=0; ii<vSources_.size(); ii++)
        {
            for(size_t jj=0; jj<genomeCount_; jj++)
                pipelineVer.insert("CGAPipeline_"+vSources_[ii]->getSource(jj));
        }

        out << "##source=" << boost::join(pipelineVer, ";");

        if (!CGA_TOOLS_IS_PIPELINE)
            out << ";cgatools_" << CGA_TOOLS_VERSION;

        out << endl;
    }

    int MkVcf::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");

        if ( genomeRoots_.size() != varNames_.size() &&
             0 != genomeRoots_.size() &&
             0 != varNames_.size() )
            throw Exception("Count of master-var parameters differs from count of genome-root parameters.");
        genomeCount_ = std::max(genomeRoots_.size(), varNames_.size());

        if ( genomeCount_ != junctionFileNames_.size() &&
            0 != genomeCount_ &&
            0 != junctionFileNames_.size() )
            throw Exception("Junction file count is inconsistent with other parameters.");
        genomeCount_ = std::max(genomeCount_, junctionFileNames_.size());

        crr_.open(referenceFileName_);
        initSources();
        validateSources();
        std::ostream& out = openStdout(outputFileName_);
        writeHeaders(out);

        while (!sources_.empty())
        {
            shared_ptr<VcfRecordSource> recs = sources_.top();
            sources_.pop();

            if (recs->eof())
                continue;

            (*recs)->writeRecord(out, crr_, genomeCount_);

            ++(*recs);
            sources_.push(recs);
        }

        return 0;
    }

    void MkVcf::transferMetadataHeaders(
        std::ostream& out) const
    {
        // recs[sourceId][genomeId][ii] -> key,value
        vector< vector< vector<VcfKvHeaderRecord> > > recs(vSources_.size());
        for(size_t ii=0; ii<vSources_.size(); ii++)
        {
            recs[ii].resize(genomeCount_);
            for(size_t jj=0; jj<genomeCount_; jj++)
            {
                recs[ii][jj] = vSources_[ii]->getKeyValueHeaderRecords(jj);
            }
        }

        set<string> genomeReference = getValues(recs, "source_GENOME_REFERENCE", true, true);
        CGA_ASSERT(1 == genomeReference.size());

        set<string> keysSeen;
        for(size_t ii=0; ii<vSources_.size(); ii++)
        {
            for(size_t jj=0; jj<genomeCount_; jj++)
            {
                for(size_t kk=0; kk<recs[ii][jj].size(); kk++)
                {
                    const VcfKvHeaderRecord& rec = recs[ii][jj][kk];
                    if (keysSeen.find(rec.key_) == keysSeen.end())
                    {
                        set<string> vals = getValues(recs, rec.key_, false, false);
                        out << "##" << rec.key_ << "=" << boost::join(vals, ";") << endl;
                        keysSeen.insert(rec.key_);
                    }
                }
            }
        }

        string gref = *genomeReference.begin();
        string assembly = gref;
        if ("NCBI build 36" == gref)
        {
            assembly = "B36";
            out << "##reference=ftp://ftp.completegenomics.com/ReferenceFiles/build36.fa.bz2" << endl;
        }
        else if ("NCBI build 37" == gref)
        {
            assembly = "B37";
            out << "##reference=ftp://ftp.completegenomics.com/ReferenceFiles/build37.fa.bz2" << endl;
        }
        else
            out << "##reference=" << gref << endl;

        for(size_t ii=0; ii<crr_.listChromosomes().size(); ii++)
        {
            out << "##contig=<ID=" << VcfRecordWriter::getVcfChromosomeName(ii, crr_)
                << ",length=" << crr_.listChromosomes()[ii].length()
                << ",assembly=" << assembly
                << ",md5=" << crr_.listChromosomes()[ii].getMd5Digest().hex()
                << ",species=\"Homo sapiens\">"
                << endl;
        }
    }

    std::set<std::string> MkVcf::getValues(
        const vector< vector< vector<VcfKvHeaderRecord> > >& recs,
        const std::string& key,
        bool mustExist,
        bool mustBeUnique) const
    {
        set<string> result;
        for(size_t ii=0; ii<recs.size(); ii++)
        {
            for(size_t jj=0; jj<recs[ii].size(); jj++)
            {
                bool exists = false;
                for(size_t kk=0; kk<recs[ii][jj].size(); kk++)
                {
                    const VcfKvHeaderRecord& rec = recs[ii][jj][kk];
                    if (rec.key_ == key)
                    {
                        result.insert(rec.value_);
                        exists = true;
                    }
                }
                if (mustExist && !exists)
                    throw Exception(key+" not supplied by "+vSourceNames_[ii]+
                                    " source for genome "+lexical_cast<string>(jj+1));
            }
        }
        if (mustBeUnique && result.size() > 1)
            throw Exception(key+" mismatch");
        return result;
    }

    void MkVcf::transferInfoFilterHeaders(
        std::ostream& out) const
    {
        // recs[sourceId][genomeId][ii] -> key,value
        vector< vector<VcfSubFieldHeaderRecord> > recs(vSources_.size());
        for(size_t ii=0; ii<vSources_.size(); ii++)
        {
            recs[ii] = vSources_[ii]->getSubFieldHeaderRecords();
        }

        for(int ii=0; ii<VcfSubFieldHeaderRecord::VCF_END; ii++)
        {
            VcfSubFieldHeaderRecord::Key key = (VcfSubFieldHeaderRecord::Key)ii;
            transferInfoFilterHeadersForKey(out, recs, key);
        }
    }

    void MkVcf::transferInfoFilterHeadersForKey(
        std::ostream& out,
        const std::vector< std::vector<cgatools::conv::VcfSubFieldHeaderRecord> >& recs,
        cgatools::conv::VcfSubFieldHeaderRecord::Key key) const
    {
        string printKey = getPrintKey(key);
        set<string> seen;
        for(size_t ii=0; ii<recs.size(); ii++)
        {
            for(size_t jj=0; jj<recs[ii].size(); jj++)
            {
                const VcfSubFieldHeaderRecord& rec = recs[ii][jj];
                if (rec.key_ != key)
                    continue;
                if (seen.find(rec.id_) == seen.end())
                {
                    checkConsistency(recs, ii, jj);
                    out << "##" << printKey << "=<ID=" << rec.id_;
                    if (rec.number_ != "")
                        out << ",Number=" << rec.number_;
                    if (rec.type_ != "")
                        out << ",Type=" << rec.type_;
                    out << ",Description=\"" << rec.description_ << "\"";
                    out << ">" << endl;
                    seen.insert(rec.id_);
                }
            }
        }
    }

    std::string MkVcf::getPrintKey(cgatools::conv::VcfSubFieldHeaderRecord::Key key) const
    {
        if (VcfSubFieldHeaderRecord::VCF_INFO == key)
            return "INFO";
        else if (VcfSubFieldHeaderRecord::VCF_FORMAT == key)
            return "FORMAT";
        else if (VcfSubFieldHeaderRecord::VCF_FILTER == key)
            return "FILTER";
        else if (VcfSubFieldHeaderRecord::VCF_ALT == key)
            return "ALT";
        CGA_ASSERT(false);
        return "UNKNOWN";
    }

    void MkVcf::checkConsistency(
        const std::vector< std::vector<cgatools::conv::VcfSubFieldHeaderRecord> >& recs,
        size_t ii0,
        size_t jj0) const
    {
        const VcfSubFieldHeaderRecord& rec = recs[ii0][jj0];
        for(size_t ii=0; ii<recs.size(); ii++)
        {
            for(size_t jj=0; jj<recs[ii].size(); jj++)
            {
                const VcfSubFieldHeaderRecord& rec2 = recs[ii][jj];
                if (rec.key_ == rec2.key_ && rec.id_ == rec2.id_)
                {
                    if (rec.number_ != rec2.number_)
                        throw Exception("mismatching Number for "+rec.id_+
                                        " for sources "+vSourceNames_[ii0]+" and "+vSourceNames_[ii]);
                    if (rec.type_ != rec2.type_)
                        throw Exception("mismatching Type for "+rec.id_+
                                        " for sources "+vSourceNames_[ii0]+" and "+vSourceNames_[ii]);
                    if (rec.description_ != rec2.description_)
                        throw Exception("mismatching Description for "+rec.id_+
                                        " for sources "+vSourceNames_[ii0]+" and "+vSourceNames_[ii]);
                }
            }
        }
    }

} } // cgatools::command
