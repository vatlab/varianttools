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
#include "cgatools/mobileelement/MeiFileVcfSource.hpp"
#include "cgatools/util/parse.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using namespace cgatools::util;
using namespace cgatools::reference;
using namespace cgatools::conv;

using boost::lexical_cast;
using boost::shared_ptr;

namespace cgatools { namespace mobileelement {

    //////////////////////////////////////////////////////////////////
    // MeiFileVcfRecordWriter
    //////////////////////////////////////////////////////////////////

    MeiFileVcfRecordWriter::MeiFileVcfRecordWriter(
        const std::vector< std::string> meiFieldNames,
        const cgatools::reference::CrrFile& crr,
        int numGenomes)
        : crr_(crr),
          meiData_(numGenomes),
          formatIds_(meiFieldNames)
    {
    }

    Location MeiFileVcfRecordWriter::getLocation() const
    {
        CGA_ASSERT(meiData_.size()==1);

        return meiData_[0].location;
    }

    void MeiFileVcfRecordWriter::writeRef(std::ostream& out) const
    {
        CGA_ASSERT(meiData_.size()==1);

        if (meiData_[0].pos == 0){
            out << ".";
        } else {
            out << crr_.getBase(meiData_[0].location);
        }
    }

    void MeiFileVcfRecordWriter::writeAlt(std::ostream& out) const
    {
        CGA_ASSERT(meiData_.size()==1);
        out << "<INS:ME:"+superType(meiData_[0])+">";
    }

    void MeiFileVcfRecordWriter::writeInfo(std::ostream& out) const
    {

        CGA_ASSERT(meiData_.size()==1);

        int numSamples = 1; // need to change this for multisample

        int elemLen = meiData_[0].elemEnd - meiData_[0].elemBegin;
        out << boost::format("IMPRECISE;SVTYPE=INS;END=%d;SVLEN=%d;CIPOS=%d,%d;MEINFO=%s,%d,%d,%s;NS=%d")
            % (meiData_[0].pos+1)
            % elemLen
            % (meiData_[0].begin-meiData_[0].pos)
            % (meiData_[0].end-meiData_[0].pos)
            % meiData_[0].ET
            % (meiData_[0].elemBegin+1)
            % meiData_[0].elemEnd
            % meiData_[0].strand 
            % numSamples;

    }

    void MeiFileVcfRecordWriter::writeFormat(std::ostream& out) const
    {
        bool wroteOne = false;
        for(size_t ii=0; ii<formatIds_.size(); ii++)
        {
            if (wroteOne)
                out << ":";
            out << formatIds_[ii];
            wroteOne = true;
        }
    }

    void MeiFileVcfRecordWriter::writeSample(std::ostream& out, size_t gIdx) const
    {
        CGA_ASSERT(meiData_.size()==1);

        CGA_ASSERT(gIdx < meiData_.size());
        bool wroteOne = false;
        BOOST_FOREACH(const std::string &f, formatIds_ ){
            if (wroteOne)
                out << ":";
            if ( f == "GT" ){
                out << ".";
            } else if ( f == "FT" ) {
                float nes = boost::lexical_cast<float>(meiData_[gIdx].KES);
                if (nes <= .75){
                    out << ".";
                } else if (nes <= .95){
                    out << "sns75";
                } else {
                    out << "sns95";
                }
            } else if ( f == "CGA_IS" ) {
                out << meiData_[gIdx].IS;
            } else if ( f == "CGA_IDC" ) {
                out << meiData_[gIdx].IDC;
            } else if ( f == "CGA_IDCL" ) {
                out << meiData_[gIdx].IDCL;
            } else if ( f == "CGA_IDCR" ) {
                out << meiData_[gIdx].IDCR;
            } else if ( f == "CGA_RDC" ) {
                if ( meiData_[gIdx].RDC < 0 ){
                    out << ".";
                } else {
                    out << meiData_[gIdx].RDC;
                }
            } else if ( f == "CGA_NBET" ) {
                if (meiData_[gIdx].NBET != ""){
                    out << meiData_[gIdx].NBET;
                } else {
                    out << ".";
                }
            } else if ( f == "CGA_ETS" ) {
                out << meiData_[gIdx].ETS;
            } else if ( f == "CGA_KES" ) {
                out << meiData_[gIdx].KES;
            } else {
                throw Exception("Could not write MEI field name " + f);
            }
            wroteOne = true;
        }
    }

    void MeiFileVcfRecordWriter::setMeiData(meiData &data, int gIdx)
    { 
        meiData_[gIdx] = data; 
        meiData_[gIdx].chrId = crr_.getChromosomeId(data.chr);
        meiData_[gIdx].pos = (int) ( (meiData_[0].begin + meiData_[0].end) / 2 );
        meiData_[gIdx].location = Location(meiData_[gIdx].chrId,meiData_[gIdx].pos);
    }

    std::string MeiFileVcfRecordWriter::superType(meiData) const
    {
        CGA_ASSERT(meiData_.size()==1);

        if( boost::to_upper_copy(meiData_[0].ET).compare(0,3,"ALU") == 0 ){
            return( std::string("ALU"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,2,"L1") == 0 ){
            return( std::string("L1"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,3,"SVA") == 0 ){
            return( std::string("SVA"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,3,"MER") == 0 ){
            return( std::string("MER"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,3,"LTR") == 0 ){
            return( std::string("LTR"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,5,"POLYA") == 0 ){
            return( std::string("POLYA"));
        }
        if( boost::to_upper_copy(meiData_[0].ET).compare(0,4,"HERV") == 0 ){
            return( std::string("HERV"));
        }
        throw Exception("Could not identify mobile element type " + meiData_[0].ET);

    }


    //////////////////////////////////////////////////////////////////
    // MeiFileVcfRecordSource
    //////////////////////////////////////////////////////////////////

    MeiFileVcfRecordSource::MeiFileVcfRecordSource(
        const std::vector< std::string >& meiFn,
        const std::vector<std::string> fieldNames,
        const CrrFile& crr)
        : crr_(crr),
          meiFn_(meiFn),
          meiIStr_(meiFn.size()),
          meiInput_(meiFn.size()),
          meiFieldNames_(0),
          eof_(false)
    {
        if (meiFn_.size()!=1){
            throw Exception("MEI conversion to VCF can only handle one genome in this version of cgatools.");
        }
        limitFieldNames(fieldNames);
        setUpDelimitedFiles();
        writer_.reset(new MeiFileVcfRecordWriter(meiFieldNames_,crr_,meiInput_.size()));
        ++(*this);
    }

    void MeiFileVcfRecordSource::limitFieldNames(const std::vector< std::string > & fieldNames){
        BOOST_FOREACH(const std::string &f,fieldNames){
            if ( f == "GT" ||
                 f == "FT" ||
                 f == "CGA_IS" ||
                 f == "CGA_IDC" ||
                 f == "CGA_IDCL" ||
                 f == "CGA_IDCR" ||
                 f == "CGA_RDC" ||
                 f == "CGA_NBET" ||
                 f == "CGA_ETS" ||
                 f == "CGA_KES"
               )
            {
                if ( meiFieldNamesSet_.find(f) == meiFieldNamesSet_.end() ){
                    meiFieldNamesSet_.insert(f);
                    meiFieldNames_.push_back(f);
                }
            }
        }
    }

    void MeiFileVcfRecordSource::setUpDelimitedFiles(){
        for(size_t ii=0;ii<meiFn_.size();++ii){
            if ( InputStream::isReadable( meiFn_[ii] ) ){
                meiIStr_[ii] = InputStream::openCompressedInputStreamByExtension(meiFn_[ii]);
                meiInput_[ii].reset( new DelimitedFile(*(meiIStr_[ii]),meiFn_[ii]) );
                meiInput_[ii]->addField(StringField("Chromosome",&(tmpMei_.chr)));
                meiInput_[ii]->addField(ValueField<int>("InsertRangeBegin",&(tmpMei_.begin)));
                meiInput_[ii]->addField(ValueField<int>("InsertRangeEnd",&(tmpMei_.end)));
                meiInput_[ii]->addField(StringField("InsertionScore",&(tmpMei_.IS)));
                meiInput_[ii]->addField(StringField("InsertionDnbCount",&(tmpMei_.IDC)));
                meiInput_[ii]->addField(StringField("InsertionLeftDnbCount",&(tmpMei_.IDCL)));
                meiInput_[ii]->addField(StringField("InsertionRightDnbCount",&(tmpMei_.IDCR)));
                meiInput_[ii]->addField(ValueField<int>("ReferenceDnbCount",&(tmpMei_.RDC))
                                        .exception("N",-1));
                meiInput_[ii]->addField(StringField("NextBestElementType",&(tmpMei_.NBET)));
                meiInput_[ii]->addField(StringField("ElementTypeScore",&(tmpMei_.ETS)));
                meiInput_[ii]->addField(StringField("KnownEventSensitivityForInsertionScore",
                                                    &(tmpMei_.KES)));
                meiInput_[ii]->addField(ValueField<int>("ElementSequenceBegin",&(tmpMei_.elemBegin)));
                meiInput_[ii]->addField(ValueField<int>("ElementSequenceEnd",&(tmpMei_.elemEnd)));
                meiInput_[ii]->addField(StringField("ElementType",&(tmpMei_.ET)));
                meiInput_[ii]->addField(StringField("Strand",&(tmpMei_.strand)));
                
            } else {
                throw Exception("Could not read file " + meiFn_[ii] );
            }
        }
    }

    cgatools::conv::VcfRecordSource& MeiFileVcfRecordSource::operator++()
    {

        // for now, we only know how to handle one genome; 
        // assert here to draw attention to code that needs to change in
        // order to support multiple genomes
        CGA_ASSERT(meiInput_.size()==1);

        if( ! meiInput_[0]->next() ){
            eof_=true;
            return *this;
        }
        writer_->setMeiData(tmpMei_,0);

        return *this;

    }

    std::vector<cgatools::conv::VcfSubFieldHeaderRecord>
    MeiFileVcfRecordSource::getSubFieldHeaderRecords() const
    {

        vector<VcfSubFieldHeaderRecord> result;

        // INFO records ...
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "NS", "1", "Integer",
                             "Number of Samples With Data"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "CIPOS", "2", "Integer",
                             "Confidence interval around POS for imprecise variants"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "END", "1", "Integer",
                             "End position of the variant described in this record"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "IMPRECISE", "0", "Flag",
                             "Imprecise structural variation"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "MEINFO", "4", "String",
                             "Mobile element info of the form NAME,START,END,POLARITY"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "SVLEN", ".", "Integer",
                             "Difference in length between REF and ALT alleles"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "SVTYPE", "1", "String",
                             "Type of structural variant"));


        // ALT records ...
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:ALU", "", "",
                             "Insertion of ALU element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:L1", "", "",
                             "Insertion of L1 element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:SVA", "", "",
                             "Insertion of SVA element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:MER", "", "",
                             "Insertion of MER element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:LTR", "", "",
                             "Insertion of LTR element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:PolyA", "", "",
                             "Insertion of PolyA element"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "INS:ME:HERV", "", "",
                             "Insertion of HERV element"));

        // FILTER records ...
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FILTER, "sns75", "", "",
                             "Sensitivity to known MEI calls in range (.75,.95] i.e. medium FDR"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FILTER, "sns95", "", "",
                             "Sensitivity to known MEI calls in range (.95,1.00] "
                             "i.e. high to very high FDR"));


        // FORMAT records ...
        if ( meiFieldNamesSet_.find("GT") != meiFieldNamesSet_.end() )
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"GT","1","String",
                             "Genotype"));
        if (meiFieldNamesSet_.find("FT") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"FT","1","String",
                             "Genotype filters"));
        if (meiFieldNamesSet_.find("CGA_IS") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_IS","1","Float",
                             "MEI InsertionScore: confidence in occurrence of an insertion"));
        if (meiFieldNamesSet_.find("CGA_IDC") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_IDC","1","Float",
                             "MEI InsertionDnbCount: count of paired ends supporting insertion"));
        if (meiFieldNamesSet_.find("CGA_IDCL") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_IDCL","1","Float",
                             "MEI InsertionLeftDnbCount: count of paired ends supporting "
                             "insertion on 5' end of insertion point"));
        if (meiFieldNamesSet_.find("CGA_IDCR") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_IDCR","1","Float",
                             "MEI InsertionRightDnbCount: count of paired ends supporting insertion "
                             "on 3' end of insertion point"));
        if (meiFieldNamesSet_.find("CGA_RDC") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_RDC","1","Integer",
                             "MEI ReferenceDnbCount: count of paired ends supporting "
                             "reference allele"));
        if (meiFieldNamesSet_.find("CGA_NBET") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_NBET","1","String",
                             "MEI NextBestElementType: (sub)type of second-most-likely inserted "
                             "mobile element"));
        if (meiFieldNamesSet_.find("CGA_ETS") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_ETS","1","Float",
                             "MEI ElementTypeScore: confidence that insertion is of type "
                             "indicated by CGA_ET/ElementType"));
        if (meiFieldNamesSet_.find("CGA_KES") != meiFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_KES","1","Float",
                             "MEI KnownEventSensitivityForInsertionScore: fraction of known "
                             "MEI insertion polymorphisms called for this sample with CGA_IS "
                             "at least as high as for the current call"));

        return result;
    }

    std::string
    MeiFileVcfRecordSource::getSource(size_t idxGenome) const
    {
        return meiInput_[idxGenome]->getMetadata().getSoftwareVersionString();
    }

    std::vector<cgatools::conv::VcfKvHeaderRecord>
    MeiFileVcfRecordSource::getKeyValueHeaderRecords(size_t idxGenome) const
    {
        vector<VcfKvHeaderRecord> result;

        const char* transferKeys[] =
        {
            "GENOME_REFERENCE",
            "MEI_1000G_ANNOTATIONS"
        };

        BOOST_FOREACH(const char* key, transferKeys)
        {
            if (meiInput_[idxGenome]->getMetadata().hasKey(key))
                result.push_back(VcfKvHeaderRecord(string("source_")+key,
                                                   meiInput_[idxGenome]->getMetadata().get(key)));
            else 
                throw Exception("Could not find metadata key " + std::string(key) + 
                                " in " + meiFn_[idxGenome]);
        }

        return result;
    }

    std::string
    MeiFileVcfRecordSource::getAssemblyId(size_t idxGenome) const
    {
        return meiInput_[idxGenome]->getMetadata().get("ASSEMBLY_ID");
    }

    bool MeiFileVcfRecordSource::eof() const
    {
        return eof_;
    }

    const cgatools::conv::VcfRecordWriter& MeiFileVcfRecordSource::operator*() const
    {
        return *writer_;
    }

    const cgatools::conv::VcfRecordWriter* MeiFileVcfRecordSource::operator->() const
    {
        return writer_.get();
    }

} } // cgatools::mobileelement
