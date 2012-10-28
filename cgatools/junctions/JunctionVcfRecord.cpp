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

//! @file JunctionVcfRecord.cpp

#include "cgatools/core.hpp"

// boost
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

// Complete Genomics.
#include "JunctionVcfRecord.hpp"
#include "cgatools/cgdata/GenomeMetadata.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"

namespace cgatools { namespace junctions {

using namespace std;
using namespace cgatools::conv;
using namespace cgatools::util;

/// implementation of class JunctionVcfRecordWriter ///////////////////

JunctionVcfRecordWriter::JunctionVcfRecordWriter()
{
}

cgatools::reference::Location JunctionVcfRecordWriter::getLocation() const
{
    return pos_;
}

void JunctionVcfRecordWriter::writeId     ( std::ostream& out ) const
{
    out << id_;
}

void JunctionVcfRecordWriter::writeRef    ( std::ostream& out ) const
{
    out << ref_;
}

void JunctionVcfRecordWriter::writeAlt    ( std::ostream& out ) const
{
    out << alt_;
}

void JunctionVcfRecordWriter::writeQual   ( std::ostream& out ) const
{
    out << qual_;
}

void JunctionVcfRecordWriter::writeInfo   ( std::ostream& out ) const
{
    out << info_;
}

void JunctionVcfRecordWriter::writeFormat ( std::ostream& out ) const
{
    out << format_;
}

void JunctionVcfRecordWriter::writeFilter ( std::ostream& out ) const
{
    out << filter_;
}

void JunctionVcfRecordWriter::writeSample ( std::ostream& out,  size_t gIdx ) const
{
    out << sample_[gIdx];
}

/// implementation of class JunctionVcfRecordSource ///////////////////

JunctionVcfRecordSource::JunctionVcfRecordSource 
( 
     const std::vector< boost::shared_ptr<cgatools::cgdata::GenomeMetadata> >& genomes,
     const std::vector<std::string>& junctionFileNames,
     const std::vector<std::string> fieldNames,
     const cgatools::reference::CrrFile& crr,
     size_t scoreThreshold,
     size_t sideLengthThreshold,
     size_t distanceTolerance,
     size_t junctionLengthThreshold,
     bool   normalPriorityOutput,
     bool   useHighConfidenceJunctionsForTumor
) 
  : crr_(crr),
    scoreThreshold_(scoreThreshold),
    sideLengthThreshold_(sideLengthThreshold),
    junctionLengthThreshold_(junctionLengthThreshold),
    distanceTolerance_(distanceTolerance),
    normalPriorityOutput_(normalPriorityOutput),
    useHighConfidenceJunctionsForTumor_(useHighConfidenceJunctionsForTumor)
{
    currentRecord_        = 0;
    fieldNames_           = fieldNames;
    writeAllFields_       = false;
    
    BOOST_FOREACH(const string& str, fieldNames_)
    {
        fieldNameSet_.insert(str);
    }

    if ( fieldNameSet_.count("SV-ALL") || fieldNameSet_.size() == 0 )
    {
        writeAllFields(true);
    }

    if ( need("FT") )
    {
        const char* filterFields[] = { "URR", "MPCBT", "SHORT", "TSNR", "INTERBL", NULL };
        for ( int i = 0; filterFields[i] != NULL; ++i )
        {
            fieldNameSet_.insert(filterFields[i]);
        }
    }

    if ( junctionFileNames.size() > 0 )
    {
        if( genomes.size() != 0 &&  genomes.size() != junctionFileNames.size() )
            throw Exception("number of input genomes is inconsistent with the number of junction files");

        if ( junctionFileNames.size() < 1 || 2 < junctionFileNames.size() )
            throw Exception("SV conversion supports only 1 or 2 genome case");

        junctionFileNames_ = junctionFileNames;
    }
    else
    {
        if ( genomes.size() < 1 || 2 < genomes.size() )
            throw Exception("SV conversion supports only 1 or 2 genome case");

        typedef boost::shared_ptr<cgatools::cgdata::GenomeMetadata> GenomeMetadataPtr;
        BOOST_FOREACH(const GenomeMetadataPtr& g, genomes)
        {
            // In case of 2 genomes, we take highConfidenceJunctions from the second genome (#1),
            // and allJunctions from the normal one (#0)
            bool highConfidence = ( useHighConfidenceJunctionsForTumor_ && 
                                    genomes.size()==2 && 
                                    junctionFileNames_.size() == 1); 
            junctionFileNames_.push_back( g->getJunctionsFileName(highConfidence) );
            //sampleIds_.push_back( g->getAsmId() );
        }
    }

    run();
}

void JunctionVcfRecordSource::writeAllFields ( bool v )
{
    writeAllFields_ = v;
}

int JunctionVcfRecordSource::run()
{
    junctionFiles_.resize(junctionFileNames_.size());
    sampleIds_.resize(junctionFileNames_.size());

    for ( int i = 0, n = junctionFileNames_.size(); i < n; ++i )
    {
        junctionFiles_[i].read(junctionFileNames_[i],crr_);
        sampleIds_[i] =  junctionFiles_[i].metadata_.get("ASSEMBLY_ID");
    }

    records_.resize(0);
    currentRecord_ = 0;

    if ( junctionFiles_.size() == 1 )
    {
        run1Genome();
    }
    else
    {
        run2Genomes();
    }

    return 0;
}

int JunctionVcfRecordSource::run1Genome()
{
    JunctionRefs junctionRefs( junctionFiles_[0].junctions_, 0 );
    vector<JunctionRefSide> all;
    all.reserve(junctionRefs.size()*2);

    BOOST_FOREACH(const JunctionRef& jr, junctionRefs)
    {
        for (size_t side = 0; side < 2; ++side)
        {
            all.push_back(JunctionRefSide(jr, side));
        }
    }
    std::sort(all.begin(), all.end());

    JunctionCompatMapPerFile compat(1);

    records_.reserve(all.size());
    BOOST_FOREACH(const JunctionRefSide& jsr, all)
    {
        JunctionVcfRecordWriter newRecord = createRecord(jsr.jr_, jsr.side_, compat);
        records_.push_back(newRecord);
    }

    return 0;
}

int JunctionVcfRecordSource::run2Genomes()
{
    vector<size_t> scoreThresholds(2,0);
    junctions::FullGenomeJunctionComparator 
        cmp( distanceTolerance_, scoreThresholds, crr_ , "");
    cmp.compare(junctionFiles_);

    vector<JunctionRefSide> all;
    JunctionCompatMapPerFile compat(junctionFiles_.size());
    for (uint32_t srcId = 0; srcId < junctionFiles_.size(); ++srcId)
    {
        copyJunctionListForVcf(cmp.getCompatible()[srcId], all, compat);
        copyJunctionListForVcf(cmp.getIncompatible()[srcId], all, compat);
    }


    std::sort(all.begin(), all.end());

    // Out of each cluster of matched junctions, we output only
    // the first one. The rest are collected in the following set.
    JunctionRefSet junctionsToSuppress;
    BOOST_FOREACH(const JunctionRefSide& jsr, all)
    {
        const JunctionRef& jr = jsr.jr_;
        if (junctionsToSuppress.count(jr) != 0)
            continue;

        if ( normalPriorityOutput_ )
        {
            pickNormalPriorityMatch( jr, junctionsToSuppress );
        }
        else
        {
            pickDefaultMatch( jr, junctionsToSuppress );
        }
    }

    records_.reserve(all.size());
    BOOST_FOREACH(const JunctionRefSide& jsr, all)
    {
        if (0 == junctionsToSuppress.count(jsr.jr_))
        {
            JunctionVcfRecordWriter newRecord = createRecord(jsr.jr_, jsr.side_, compat);
            records_.push_back(newRecord);
        }
    }

    return 0;
}

//! Writes one side of a VCF "adjacency" to the given stream.
//! Note that, regardless of our junction side strand, VCF record
//! is always written relative to the primary strand.
//! Forks each junction into two single-side structures so that we
//! can sort the sides in the reference order. Also populates the
//! compatibility maps for each file, which are necessary because
//! the comparison algorithm populates the matchedJunctionId_ field
//! only in one direction, and we want symmetric references in the
//! VCF file.
void JunctionVcfRecordSource::copyJunctionListForVcf
(
    const junctions::JunctionRefs& jrl,
    vector<JunctionRefSide>& out,
    JunctionCompatMapPerFile& compat
)   const
{
    BOOST_FOREACH(const JunctionRef& jr, jrl)
    {
        BOOST_FOREACH(const JunctionRef& cjr, jr.matchedJunctionIts_)
        {
            if (jr.sourceId_ == cjr.sourceId_)
                continue; // shouldn't really happen
            compat[jr.sourceId_][jr.it_->id_].
                push_back(cjr);
            compat[cjr.sourceId_][cjr.it_->id_].
                push_back(jr);
        }

        for (size_t side = 0; side < 2; ++side)
            out.push_back(JunctionRefSide(jr, side));
    }
}

void JunctionVcfRecordSource::pickNormalPriorityMatch
( 
    const JunctionRef& jr, 
    JunctionRefSet& junctionsToSuppress 
)   const
{
    const JunctionRef* keep = &jr;

    BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
    {
        if (junctionsToSuppress.count(match) != 0)
            continue;

        if ( keep->sourceId_ > match.sourceId_ ) // normal has lower source ID
        {
            keep = &match;
        }
        else if ( keep->sourceId_ == match.sourceId_ ) 
        { // if source ID is the same then take junction with lower position
            if ( match < *keep )
            {
                keep = &match;
            }
        }
    }

    if ( !(*keep == jr) )
        junctionsToSuppress.insert(jr);

    BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
    {
        if ( *keep == match )
            continue;
        junctionsToSuppress.insert(match);
    }
}

void JunctionVcfRecordSource::pickDefaultMatch
( 
    const JunctionRef& jr, 
    JunctionRefSet& junctionsToSuppress 
)   const
{
    BOOST_FOREACH(const JunctionRef& match, jr.matchedJunctionIts_)
    {
        if (jr == match)
            continue;
        junctionsToSuppress.insert(match);
    }
}

JunctionVcfRecordWriter JunctionVcfRecordSource::createRecord
( 
  const JunctionRef& jref, size_t side,   
  const JunctionCompatMapPerFile& compat ) const
{
    JunctionVcfRecordWriter record;

    const junctions::JunctionSideSection& jss = jref.it_->sideSections_[side];

    record.pos_ = jss.getBasePos(side);
    record.id_  = getId(jref,side);

    record.ref_ = crr_.getBase(jss.getBasePos(side));
    record.alt_ = getAltField(jref,side,false);

    record.qual_ = ".";
    record.filter_ = ".";

    record.info_ = getInfo(jref,side);
    record.format_ = getFormat(jref,side);

    record.sample_.resize(sampleIds_.size());
    for ( size_t idx = 0, n = sampleIds_.size(); idx < n; ++idx )
    {
        record.sample_[idx] = getSample(jref,side,idx, compat);
    }

    return record;
}

std::string JunctionVcfRecordSource::addFilterFlag( const std::string& flag, bool& filtered) const
{
    string out;
    if ( filtered ) out +=  ";";
    out += flag;
    filtered = true;
    return out;
}

std::string JunctionVcfRecordSource::getSampleFilter( const JunctionRef& jref ) const
{
    string out;
    const junctions::Junction& j = *jref.it_;
    bool filtered = false;

    if (  j.score_ < scoreThreshold_ ) out += addFilterFlag( "MPCBT", filtered);
    if (  j.knownUnderrepresentedRepeat_ == "Y" ) out += addFilterFlag(  "URR" , filtered);
    if ( !j.transitionIsKnown_                  ) out += addFilterFlag( "TSNR" , filtered);
    if (  j.sideSections_[0].length_ < (int)sideLengthThreshold_ ||
          j.sideSections_[1].length_ < (int)sideLengthThreshold_ )
                                                  out += addFilterFlag( "SHORT", filtered);
#if 1
    if ( j.frequencyInBaselineGenomeSet_ > 0.0 &&
         ( j.sideSections_[0].position_.chromosome_ !=
           j.sideSections_[1].position_.chromosome_ ))
                                                  out += addFilterFlag( "INTERBL", filtered);
#endif
    if (!filtered)
        out += "PASS";
    
    return out;
}

bool JunctionVcfRecordSource::need ( const std::string& fieldName ) const
{
    return writeAllFields_ 
           || (fieldNameSet_.count(fieldName) > 0);
}

std::string JunctionVcfRecordSource::getSample ( const JunctionRef& jref, size_t side ) const
{
    // GT:FT:CGA_BNDMPC:CGA_BNDPOS:CGA_BNDDEF:CGA_BNDP
    string out;
    
    if ( need("GT") ) out += "1:";

    if ( need("FT") ) out += getSampleFilter(jref) + ":";

    if ( need("CGA_BNDMPC") ) 
        out += boost::lexical_cast<string>(jref.it_->score_) + ":";

    // CGA_BNDPOS - position of the junction, without chromosome
    const junctions::JunctionSideSection& jss = jref.it_->sideSections_[side];
    reference::Location basePos = jss.getBasePos(side);
    if ( need("CGA_BNDPOS") ) out += getPosition(basePos, "") + ":";

    // CGA_BNDDEF - same as the ALT field, but no chromosome name
    if ( need("CGA_BNDDEF") ) out += getAltField(jref, side, true) + ":";

    if ( need("CGA_BNDP") ) 
        out += (jref.it_->transitionIsKnown_ ? "PRECISE" : "IMPRECISE");

    if (!out.empty() && out[out.length()-1] == ':') out = out.substr(0,out.length()-1);
    return out;
}

std::string JunctionVcfRecordSource::getSample
( 
    const JunctionRef& jref, 
    size_t side, 
    size_t idx, 
    const JunctionCompatMapPerFile& compat 
)    const
{
    uint32_t srcId = idx; //junctionFiles_.size()-1-idx;
    
    if (jref.sourceId_ == srcId) 
    {
        return getSample(jref, side);
    }
    const JunctionRef* match = 0;

    JunctionCompatMap::const_iterator matchMapIt = compat[jref.sourceId_].find(jref.it_->id_);

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
    {
        return getSample(*match, side);
    }    

    //return string(".:.:.:.:.:.");
    string out;

    for ( int i = 0; sampleFieldIDs_[i] != NULL; ++i )
    {
        if ( need(sampleFieldIDs_[i]) )
        {
            if ( !out.empty() ) out += ":";
            out += ".";
        }
    }

    return out;
}


std::string JunctionVcfRecordSource::getId( const JunctionRef& jref, size_t side) const
{
    return sampleIds_[jref.sourceId_] + "_" + jref.it_->id_ + (side == 0 ? "_L" : "_R");
}


std::string JunctionVcfRecordSource::getMEI(const std::string& med) const
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

const char* JunctionVcfRecordSource::sampleFieldIDs_[] =
{
    "GT", "FT", "CGA_BNDMPC", "CGA_BNDPOS", "CGA_BNDDEF", "CGA_BNDP", NULL
};

std::string JunctionVcfRecordSource::getFormat( const JunctionRef& jref, size_t side) const
{

    string out;

    for ( int i = 0; sampleFieldIDs_[i] != NULL; ++i )
    {
        if ( need(sampleFieldIDs_[i]) )
        {
            if ( !out.empty() ) out += ":";
            out += sampleFieldIDs_[i];
        }
    }

    return out;
}


std::string JunctionVcfRecordSource::getInfo( const JunctionRef& jref, size_t side) const
{
    namespace ba = boost::algorithm;
    const junctions::Junction& j = *jref.it_;
    size_t otherSide = 1 - side;

    string out;

    if ( need("NS"       ) ) out  = out + ";NS=" + boost::lexical_cast<string>(sampleIds_.size());
    if ( need("SVTYPE"   ) ) out += ";SVTYPE=BND";
    if ( need("MATEID"   ) ) out  = out + ";MATEID=" + getId(jref, otherSide);
    if ( need("CGA_BF"   ) ) 
    {
        out  = out + ";CGA_BF=" 
             + (boost::format("%.2f") % j.frequencyInBaselineGenomeSet_).str();
    }

    if ( need("CGA_XR"   ) ) 
    {
        if (!j.xRef_.empty())
        {
            vector<string> dbtokens;
            ba::split(dbtokens, j.xRef_, ba::is_any_of(" "));
            if (!dbtokens.empty())
                out = out +  ";CGA_XR=" + dbtokens[0];
        }
    }

    if ( need("CGA_MEDEL")) 
    {
        if (!j.deletedTransposableElement_.empty())
        {
            string meinfo = getMEI(j.deletedTransposableElement_);
            if (!meinfo.empty())
                out = out + ";CGA_MEDEL=" + meinfo;
        }
    }


    if ( need("CGA_BNDG") || need("CGA_BNDGO") ) 
    {
        for(size_t ii=0; ii<2; ii++)
        {
            size_t ss = 0 == ii ? side : otherSide;
            const junctions::JunctionSideSection& jss = jref.it_->sideSections_[ss];
            if ("" == jss.genes_)
                continue;
            string genes = jss.genes_;
            boost::replace_all(genes, ":", "|");
            boost::replace_all(genes, ";", "&");

            if ( ss == side ) 
            {
                if ( need("CGA_BNDG") )
                    out = out + ";CGA_BNDG=" + genes;
            }
            else
            {
                if ( need("CGA_BNDGO") )
                    out = out + ";CGA_BNDGO=" + genes;
            }
        }
    }

    if (out.length() > 0 && out[0] == ';') out = out.substr(1);

    return out;

}

std::string JunctionVcfRecordSource::getAltField
( 
    const JunctionRef& jref, 
    size_t side, 
    bool suppressChrom 
)   const
{
    size_t otherSide = 1 - side;

    const junctions::JunctionSideSection& jss      = jref.it_->sideSections_[side];
    const junctions::JunctionSideSection& otherJss = jref.it_->sideSections_[otherSide];

    int dir      = jss.getDir(side),
        otherDir = otherJss.getDir(otherSide);

    reference::Location basePos      = jss.getBasePos(side),
                        otherBasePos = otherJss.getBasePos(otherSide);

    string sf(1,crr_.getBase(basePos));
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
        oap = "[" + getPosition(otherBasePos, sep) + "[";
    else
        oap = "]" + getPosition(otherBasePos, sep) + "]";

    string out = (dir == 1)
               ? (sf + tf + oap)
               : (oap + tf + sf);
    return out;
}


std::string JunctionVcfRecordSource::getPosition
    ( 
        reference::Location pos, 
        const std::string& sep 
    )   const
{
    std::ostringstream ss;
    if ( !sep.empty()) 
    {
        std::string chrom = crr_.listChromosomes()[pos.chromosome_].getName();
        if (chrom.find("chr") == 0)
            chrom = chrom.substr(3);
            ss << chrom << sep;
    }
    ss << (pos.offset_ + 1);
    return ss.str();
}

bool JunctionVcfRecordSource::add 
( 
    vector<VcfSubFieldHeaderRecord>& result, 
    VcfSubFieldHeaderRecord::Key key,
    const string& id, const string& number,
    const string& type, 
    const string& description
)   const
{
    if ( !need(id))
    {
        return false;
    }

    result.push_back(VcfSubFieldHeaderRecord( key, id, number, type, description) );

    return true;
}

std::vector<cgatools::conv::VcfSubFieldHeaderRecord> 
    JunctionVcfRecordSource::getSubFieldHeaderRecords() const
{
    vector<VcfSubFieldHeaderRecord> r;

    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "NS"          , "1", "Integer", 
        "Number of Samples With Data");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "CGA_BF"      , "1", "Float"  , 
        "Frequency in baseline");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "CGA_MEDEL"   , "4", "String" , 
        "Consistent with deletion of mobile element; type,chromosome,start,end");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "CGA_XR"      , "A", "String" , 
        "Per-ALT external database reference (dbSNP, COSMIC, etc)");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "MATEID"      , "1", "String" , 
        "ID of mate breakend");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "SVTYPE"      , "1", "String" , 
        "Type of structural variant");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "CGA_BNDG"    , "A", "String" , 
        "Transcript name and strand of genes containing breakend");
    add( r,  VcfSubFieldHeaderRecord::VCF_INFO,   "CGA_BNDGO"   , "A", "String" , 
        "Transcript name and strand of genes containing mate breakend");
                                                             
    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "GT"          , "1", "String", 
        "Genotype" );
    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "FT"          , "1", "String" , 
        "Genotype filters" );

    if ( need("FT") )
    {
        add( r,  VcfSubFieldHeaderRecord::VCF_FILTER, "URR"         , "",  ""       , 
            "Too close to an underrepresented repeat" );
        add( r,  VcfSubFieldHeaderRecord::VCF_FILTER, "MPCBT"       , "",  ""       , 
            "Mate pair count below "+boost::lexical_cast<string>(scoreThreshold_));
        add( r,  VcfSubFieldHeaderRecord::VCF_FILTER, "SHORT"       , "",  ""       , 
            "Junction side length below "+boost::lexical_cast<string>(sideLengthThreshold_) );
        add( r,  VcfSubFieldHeaderRecord::VCF_FILTER, "TSNR"        , "",  ""       , 
            "Transition sequence not resolved" );
#if 1
        add( r,  VcfSubFieldHeaderRecord::VCF_FILTER, "INTERBL" , "",  ""       , 
            "Interchromosomal junction in baseline" );
#endif
    }

    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "CGA_BNDMPC"  , "1", "Integer", 
        "Mate pair count supporting breakend" );
    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "CGA_BNDPOS"  , "1", "Integer",
        "Breakend position" );
    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "CGA_BNDDEF"  , "1", "String" , 
        "Breakend definition" );
    add( r,  VcfSubFieldHeaderRecord::VCF_FORMAT, "CGA_BNDP"    , "1", "String" , 
        "Precision of breakend" );
  
    return r;
}


std::string JunctionVcfRecordSource::getSource(size_t idxGenome) const
{
    return junctionFiles_[idxGenome].metadata_.getSoftwareVersionString();
}


std::vector<cgatools::conv::VcfKvHeaderRecord> 
    JunctionVcfRecordSource::getKeyValueHeaderRecords(size_t idxGenome) const
{
    vector<VcfKvHeaderRecord> result;

    const char* transferKeys[] =
    {
        "GENOME_REFERENCE",
        "GENE_ANNOTATIONS",
        "DBSNP_BUILD",
        "COSMIC",
        "DGV_VERSION",
        "MIRBASE_VERSION",
        "PFAM_DATE",
        "REPMASK_GENERATED_AT",
        "SEGDUP_GENERATED_AT"
    };

    const util::DelimitedFile::Metadata& md = junctionFiles_[idxGenome].metadata_;

    BOOST_FOREACH(const char* key, transferKeys)
    {
        if ( md.hasKey(key) )
        {
            result.push_back(VcfKvHeaderRecord(string("source_")+key,md.get(key)));
        }
    }

    return result;
}


std::string JunctionVcfRecordSource::getAssemblyId(size_t idxGenome) const
{
    return sampleIds_[idxGenome];
}

bool JunctionVcfRecordSource::eof() const
{
    return (currentRecord_ >= records_.size() );
}

cgatools::conv::VcfRecordSource& JunctionVcfRecordSource::operator++()
{
    if ( !eof() ) ++currentRecord_;
    return (*this);
}

const cgatools::conv::VcfRecordWriter& JunctionVcfRecordSource::operator*() const
{
    return records_[currentRecord_];
}

const cgatools::conv::VcfRecordWriter* JunctionVcfRecordSource::operator->() const
{
    return &records_[currentRecord_];
}

}} // end namespace
