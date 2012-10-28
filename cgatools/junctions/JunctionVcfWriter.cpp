// Copyright 2012 Complete Genomics, Inc.
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
#include "cgatools/junctions//JunctionVcfWriter.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/regex.hpp>

#include "cgatools/util/Streams.hpp"
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/RepeatMaskerStore.hpp"
#include "cgatools/reference/GeneDataStore.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/util/BaseUtil.hpp"

using std::string;
using std::vector;
using std::set;
using std::map;
using boost::shared_ptr;
using namespace cgatools::junctions;
using namespace cgatools::reference;
namespace ba = boost::algorithm;

namespace cgatools { namespace junctions {


    JunctionVcfWriter::JunctionVcfWriter()
{
    CGA_ERROR_EX("The default constructor or JunctionVcfWriter should never be invoked");
}


JunctionVcfWriter::JunctionVcfWriter ( boost::shared_ptr<CrrFile> ref, 
                                       boost::shared_ptr<junctions::JunctionFiles> junctionFiles )
: fileFieldSeparator_   ("\t"), 
  filterScoreThreshold_ (0),
  filterSideLength_     (0)
{
    reference_     = ref;
    junctionFiles_ = junctionFiles;

    init();
}

void JunctionVcfWriter::init()
{
    sampleIds_.reserve( junctionFiles_->size() );
    BOOST_FOREACH( const JunctionFile& jf, *junctionFiles_)
    {
        sampleIds_.push_back( jf.metadata_.get("ASSEMBLY_ID") );
    }
}

void JunctionVcfWriter::writeJunctionVcfHeaders(std::ostream& out) const
{
    out <<
        "##fileformat=VCFv4.1\n"
        "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n"
        "##INFO=<ID=CGA_BF,Number=1,Type=Float,Description=\"Frequency in baseline\">\n"
        "##INFO=<ID=CGA_MEDEL,Number=4,Type=String,Description=\"Consistent with deletion "
        "of mobile element; type,chromosome,start,end\">\n"
        "##INFO=<ID=CGA_XR,Number=A,Type=String,"
        "Description=\"Per-ALT external database reference (dbSNP, COSMIC, etc)\">\n"
        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakend\">\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n"
        "##INFO=<ID=CGA_BNDG,Number=A,Type=String,"
        "Description=\"Transcript name and strand of genes containing breakend\">\n"
        "##INFO=<ID=CGA_BNDGO,Number=A,Type=String,"
        "Description=\"Transcript name and strand of genes containing mate breakend\">\n"
        "##FILTER=<ID=URR,Description=\"Too close to an underrepresented repeat\">\n"
        "##FILTER=<ID=MPCBT,Description=\"Mate pair count below "
        << filterScoreThreshold_ << "\">\n" <<
        "##FILTER=<ID=SHORT,Description=\"Junction side length below "
        << filterSideLength_  << "\">\n" <<
        "##FILTER=<ID=TSNR,Description=\"Transition sequence not resolved\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Genotype filters\">\n"
        "##FORMAT=<ID=CGA_BNDMPC,Number=1,Type=Integer,"
        "Description=\"Mate pair count supporting breakend\">\n"
        "##FORMAT=<ID=CGA_BNDPOS,Number=1,Type=Integer,Description=\"Breakend position\">\n"
        "##FORMAT=<ID=CGA_BNDDEF,Number=1,Type=String,Description=\"Breakend definition\">\n"
        "##FORMAT=<ID=CGA_BNDP,Number=1,Type=String,Description=\"Precision of breakend\">\n"
        ;
    out << "#CHROM" << fileFieldSeparator_
        << "POS" << fileFieldSeparator_
        << "ID" << fileFieldSeparator_
        << "REF" << fileFieldSeparator_
        << "ALT" << fileFieldSeparator_
        << "QUAL" << fileFieldSeparator_
        << "FILTER" << fileFieldSeparator_
        << "INFO" << fileFieldSeparator_
        << "FORMAT";

    for (uint32_t ii=0; ii<junctionFiles_->size(); ii++)
    {
        // Expectation is for command line arguments for
        // junctiondiff to give tumor first, normal second. But
        // the VCF has normal first, tumor second.
        uint32_t srcId = junctionFiles_->size()-1-ii;
        out << fileFieldSeparator_ << sampleIds_[srcId];
    }
    out << "\n";
}


std::string JunctionVcfWriter::createJunctionVcfId( const JunctionRef& jref, size_t side) const
{
    std::ostringstream ss;
    ss << sampleIds_[jref.sourceId_]
    << '_' << jref.it_->id_
        << (side == 0 ? "_L" : "_R");
    return ss.str();
}

std::string JunctionVcfWriter::formatPositionForVcf( reference::Location pos, 
                                                     const std::string& sep ) const
{
    std::ostringstream ss;
    if ( !sep.empty()) {
        std::string chrom = reference_->listChromosomes()[pos.chromosome_].getName();
        if (chrom.find("chr") == 0)
            chrom = chrom.substr(3);
        ss << chrom << sep;
    }
    ss << (pos.offset_ + 1);
    return ss.str();
}

void JunctionVcfWriter::writeJunctionPositionToVcf(const junctions::JunctionSideSection& jss,
                                size_t side, std::ostream& out) const
{
    reference::Location pos = jss.getBasePos(side);
    out << formatPositionForVcf(pos,fileFieldSeparator_);
}

void JunctionVcfWriter::writeJunctionAltFieldToVcf(const JunctionRef& jref,
                                size_t side, std::ostream& out, bool suppressChrom ) const
{
    size_t otherSide = 1 - side;
    const junctions::JunctionSideSection& jss = jref.it_->sideSections_[side];
    const junctions::JunctionSideSection& otherJss =
        jref.it_->sideSections_[otherSide];
    int dir = jss.getDir(side),
        otherDir = otherJss.getDir(otherSide);
    reference::Location basePos = jss.getBasePos(side),
        otherBasePos = otherJss.getBasePos(otherSide);

    string sf(1, reference_->getBase(basePos));
    string tf;
    if (jref.it_->transitionIsKnown_)
    {
        tf = jref.it_->transitionSequence_;
        if (jss.strand_ != junctions::JUNCTION_PLUS_STRAND)
            tf = util::baseutil::reverseComplement(tf);
    }
    else
    {
        tf = "N";
    }
    string oap;
    string sep = suppressChrom ? "" : ":";
    if (otherDir == -1)
        oap = "[" + formatPositionForVcf(otherBasePos, sep) + "[";
    else
        oap = "]" + formatPositionForVcf(otherBasePos, sep) + "]";
    if (dir == 1)
        out << sf << tf << oap;
    else
        out << oap << tf << sf;
}


string JunctionVcfWriter::convertMobileElementToVcf(const std::string& med) const
{
    boost::regex pattern("([^ ]+) [^ ]+ \\(chr([^:]+):([\\d]+)-([\\d]+)\\)");
    boost::smatch match;
    if (!boost::regex_match(med, match, pattern))
        return "";
    std::ostringstream ss;
    ss << match[1]
    << "," << match[2]
    << "," << (1 + boost::lexical_cast<size_t>(match[3]))
        << "," << match[4];
    return ss.str();
}

void JunctionVcfWriter::writeJunctionInfoFieldToVcf(const JunctionRef& jref,
                                 size_t side, std::ostream& out) const
{
    namespace ba = boost::algorithm;
    const junctions::Junction& j = *jref.it_;
    size_t otherSide = 1 - side;
    out << "NS=2;SVTYPE=BND"
        << ";MATEID=" << createJunctionVcfId(jref, otherSide)
        << ";CGA_BF=" << (boost::format("%.2f") % j.frequencyInBaselineGenomeSet_);
    if (!j.xRef_.empty())
    {
        vector<string> dbtokens;
        ba::split(dbtokens, j.xRef_, ba::is_any_of(" "));
        if (!dbtokens.empty())
            out << ";CGA_XR=" << dbtokens[0];
    }
    if (!j.deletedTransposableElement_.empty())
    {
        string meinfo =
            convertMobileElementToVcf(j.deletedTransposableElement_);
        if (!meinfo.empty())
            out << ";CGA_MEDEL=" << meinfo;
    }
    for(size_t ii=0; ii<2; ii++)
    {
        size_t ss = 0 == ii ? side : otherSide;
        const junctions::JunctionSideSection& jss = jref.it_->sideSections_[ss];
        if ("" == jss.genes_)
            continue;
        string genes = jss.genes_;
        boost::replace_all(genes, ":", "|");
        boost::replace_all(genes, ";", "&");
        out << (ss == side ? ";CGA_BNDG=" : ";CGA_BNDGO=")
            << genes;
    }
}

void JunctionVcfWriter::addFilterFlag( std::ostream& out, 
                                       const std::string& flag, bool& filtered) const
{
    if (filtered)
        out << ";";
    out << flag;
    filtered = true;
}

void JunctionVcfWriter::writeJunctionFilterFieldToVcf( const JunctionRef& jref, 
                                                       std::ostream& out) const
{
    const junctions::Junction& j = *jref.it_;
    bool filtered = false;
    if (j.score_ < filterScoreThreshold_)
        addFilterFlag(out, "MPCBT", filtered);
    if (j.knownUnderrepresentedRepeat_ == "Y")
        addFilterFlag(out, "URR", filtered);
    if (!j.transitionIsKnown_)
        addFilterFlag(out, "TSNR", filtered);
    int filterSideLength = filterSideLength_;
    if (j.sideSections_[0].length_ < filterSideLength ||
        j.sideSections_[1].length_ < filterSideLength)
    {
        addFilterFlag(out, "SHORT", filtered);
    }
    if (!filtered)
        out << "PASS";
}

void JunctionVcfWriter::writeJunctionComparisonField(const JunctionRef& jref,
                                  size_t side, std::ostream& out) const
{
    // GT:FT:CGA_BNDMPC:CGA_BNDPOS:CGA_BNDDEF:CGA_BNDP
    out << "1:"; // GT subfield
    writeJunctionFilterFieldToVcf(jref, out); // FT - filters
    out << ":";
    out << jref.it_->score_ << ":"; // CGA_BNDMPC - count of DNBs

    // CGA_BNDPOS - position of the junction, without chromosome
    const junctions::JunctionSideSection& jss = jref.it_->sideSections_[side];
    reference::Location basePos = jss.getBasePos(side);
    out << formatPositionForVcf(basePos, "") << ":";
    // CGA_BNDDEF - same as the ALT field, but no chromosome name
    writeJunctionAltFieldToVcf(jref, side, out, true);
    if (jref.it_->transitionIsKnown_)
        out << ":PRECISE";
    else
        out << ":IMPRECISE";
}

void JunctionVcfWriter::writeJunctionToVcf(const JunctionRef& jref,
                        size_t side,
                        const JunctionCompatMapPerFile& compat,
                        std::ostream& out) const
{
    const junctions::JunctionSideSection& jss = jref.it_->sideSections_[side];
    writeJunctionPositionToVcf(jss, side, out);
    out << fileFieldSeparator_ << createJunctionVcfId(jref, side);
    reference::Location basePos = jss.getBasePos(side);
    out << fileFieldSeparator_ << reference_->getBase(basePos);
    out << fileFieldSeparator_;
    writeJunctionAltFieldToVcf(jref, side, out);
    out << fileFieldSeparator_ << "."; // QUAL is mandatory, but not known
    out << fileFieldSeparator_ << "."; // FILTER is not defined (looks
    // like it was designed for
    // single-genome VCFs)
    out << fileFieldSeparator_;
    writeJunctionInfoFieldToVcf(jref, side, out);
    out << fileFieldSeparator_ << "GT:FT:CGA_BNDMPC:CGA_BNDPOS:CGA_BNDDEF:CGA_BNDP";

    // Write comparison columns
    CGA_ASSERT( junctionFiles_->size() <= 2 );

    for (uint32_t ii=0; ii<junctionFiles_->size(); ii++)
    {
        // Expectation is for command line arguments for
        // junctiondiff to give tumor first, normal second. But
        // the VCF has normal first, tumor second.
        uint32_t srcId = junctionFiles_->size()-1-ii;
        out << fileFieldSeparator_;
        // Skip the column that corresponds to the file that produced
        // this junction
        if (jref.sourceId_ == srcId) {
            writeJunctionComparisonField(jref, side, out);
            continue;
        }
        const JunctionRef* match = 0;
        JunctionCompatMap::const_iterator matchMapIt =
            compat[jref.sourceId_].find(jref.it_->id_);
        if (matchMapIt != compat[jref.sourceId_].end())
        {
            const junctions::JunctionRefs& matchedRefs = matchMapIt->second;
            BOOST_FOREACH(const JunctionRef& mjr, matchedRefs)
            {
                if (mjr.sourceId_ == srcId)
                {
                    match = &mjr;
                    break;
                }
            }
        }
        if (match != 0)
            writeJunctionComparisonField(*match, side, out);
        else
            out << ".:.:.:.:.:.";;
    }
    out << "\n";


}

}} 
