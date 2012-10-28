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
#include "cgatools/copynumber/CnvFileVcfSource.hpp"
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

namespace cgatools { namespace copynumber {

    /////////////////////////////////////////////////////////////////
    // CnvFileVcfRecordWriter
    /////////////////////////////////////////////////////////////////

    CnvFileVcfRecordWriter::CnvFileVcfRecordWriter(
        const std::vector< std::string> cnvFieldNames,
        const cgatools::reference::CrrFile& crr,
        int numGenomes,
        bool someSomatic)
        : crr_(crr),
          someSomatic_(someSomatic),
          dipData_(numGenomes),
          nondipData_(numGenomes),
          somData_(numGenomes),
          formatIds_(cnvFieldNames),
          deferred_(numGenomes,false)
    {
    }

    Location CnvFileVcfRecordWriter::getLocation() const
    {
        size_t ii=0;
        for(;ii<dipData_.size();++ii){
            if( ! deferred_[ii] ) break;
        }
        CGA_ASSERT( ii < dipData_.size() );

        return dipData_[ii].range.beginLocation();
    }

    void CnvFileVcfRecordWriter::writeRef(std::ostream& out) const
    {
        out << crr_.getBase(getLocation());
    }

    void CnvFileVcfRecordWriter::writeAlt(std::ostream& out) const
    {
        out << "<CGA_CNVWIN>";
    }

    void CnvFileVcfRecordWriter::writeInfo(std::ostream& out) const
    {
        size_t ii=0;
        for(;ii<dipData_.size();++ii){
            if( ! deferred_[ii] ) break;
        }
        CGA_ASSERT( ii < dipData_.size() );

        out << boost::format("NS=%d;CGA_WINEND=%d") % dipData_.size() % dipData_[ii].end ;
    }

    void CnvFileVcfRecordWriter::writeFormat(std::ostream& out) const
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

    void CnvFileVcfRecordWriter::writeSample(std::ostream& out, size_t gIdx) const
    {
        CGA_ASSERT(gIdx < dipData_.size());
        bool wroteOne = false;
        BOOST_FOREACH(const std::string &f, formatIds_ ){
            if (wroteOne)
                out << ":";
            if ( f == "GT" ){
                out << ".";
            } else if ( f == "CGA_GP" ) {
                if ( ! deferred_[gIdx] && dipData_[gIdx].normedGcCvg >= 0 ){
                    out << boost::format("%.2f") % dipData_[gIdx].normedGcCvg;
                } else {
                    out << ".";
                }
            } else if ( f == "CGA_NP" ) {
                if ( ! deferred_[gIdx] && dipData_[gIdx].relCvg >= 0 ){
                    out << boost::format("%.2f") % dipData_[gIdx].relCvg;
                } else {
                    out << ".";
                }
            } else if ( f == "CGA_CP" ) {
                if ( ! deferred_[gIdx] && dipData_[gIdx].cP >= 0 ){
                    out << boost::format("%d") % dipData_[gIdx].cP;
                } else {
                    out << ".";
                }
            } else if ( f == "CGA_PS" ) {
                if ( ! deferred_[gIdx] ){
                    out << dipData_[gIdx].pS;
                } else {
                    out << "0";
                }
            } else if ( f == "CGA_CT" ) {
                if ( deferred_[gIdx] ||
                     dipData_[gIdx].cT == "hypervariable" ||
                     dipData_[gIdx].cT == "invariant" ||
                     dipData_[gIdx].cT == "N" ){
                    out << ".";
                } else {
                    out << dipData_[gIdx].cT;
                }
            } else if ( f == "CGA_TS" ) {
                if ( ! deferred_[gIdx] ){
                    out << dipData_[gIdx].tS;
                } else {
                    out << "0";
                }
            } else if ( f == "CGA_CL" ) {
                if ( ! deferred_[gIdx] && nondipData_[gIdx].cL >= 0 ){
                    out << boost::format("%.3f") % nondipData_[gIdx].cL;
                } else {
                    out << ".";
                }
             } else if ( f == "CGA_LS" ) {
                if ( ! deferred_[gIdx] ){
                    out << nondipData_[gIdx].lS;
                } else {
                    out << "0";
                }
            } else if ( f == "CGA_SCL" ) {
                if ( ! deferred_[gIdx] && somData_[gIdx].range != Range() && somData_[gIdx].cL >= 0 ){
                    out << boost::format("%.3f") % somData_[gIdx].cL;
                } else {
                    out << ".";
                }
            } else if ( f == "CGA_SLS" ) {
                if ( ! deferred_[gIdx] && somData_[gIdx].range != Range() )  out << somData_[gIdx].lS;
                else if (someSomatic_) out << ".";
            } else if ( f == "CGA_LAF" ) {
                if ( ! deferred_[gIdx] && somData_[gIdx].range != Range() )  out << somData_[gIdx].laf;
                else if (someSomatic_) out << ".";
            } else if ( f == "CGA_LLAF" ) {
                if ( ! deferred_[gIdx] && somData_[gIdx].range != Range() )  out << somData_[gIdx].llaf;
                else if (someSomatic_) out << ".";
            } else if ( f == "CGA_ULAF" ) {
                if ( ! deferred_[gIdx] && somData_[gIdx].range != Range() )  out << somData_[gIdx].ulaf;
                else if (someSomatic_) out << ".";
            } else {
                throw Exception("Could not write CNV field name " + f);
            }
            wroteOne = true;
        }
    }

    void CnvFileVcfRecordWriter::setDiploidData(diploidData &data, int gIdx)
    { 
        dipData_[gIdx] = data; 
        uint16_t chr = crr_.getChromosomeId(data.chr);
        dipData_[gIdx].range = Range(chr,data.begin,data.end);
    }

    void CnvFileVcfRecordWriter::setNondiploidData(nondiploidData &data, int gIdx)
    { 
        nondipData_[gIdx] = data; 
        uint16_t chr = crr_.getChromosomeId(data.chr);
        nondipData_[gIdx].range = Range(chr,data.begin,data.end);
    }

    void CnvFileVcfRecordWriter::setSomaticData(somaticData &data, int gIdx)
    { 
        somData_[gIdx] = data; 
        uint16_t chr = crr_.getChromosomeId(data.chr);
        somData_[gIdx].range = Range(chr,data.begin,data.end);
    }


    //////////////////////////////////////////////////////////////////
    // CnvFileVcfRecordSource
    //////////////////////////////////////////////////////////////////

    CnvFileVcfRecordSource::CnvFileVcfRecordSource(
        const std::vector< std::string >& dipdetFn,
        const std::vector< std::string >& nondipdetFn,
        const std::vector< std::string >& somnondipdetFn,
        const std::vector<std::string> fieldNames,
        const CrrFile& crr)
        : crr_(crr),
          genomeEnd_(cgatools::reference::Range(crr_.listChromosomes().size(),0,0)),
          dipDetFn_(dipdetFn),
          nondipDetFn_(nondipdetFn),
          somnondipDetFn_(somnondipdetFn),
          dipDetIStr_(dipdetFn.size()),
          nondipDetIStr_(nondipdetFn.size()),
          somnondipDetIStr_(dipdetFn.size()),
          dipDet_(dipdetFn.size()),
          nondipDet_(nondipdetFn.size()),
          somnondipDet_(dipdetFn.size()),
          gcCorrectedMean_(dipdetFn.size()),
          cnvFieldNames_(0),
          eof_(false)
    {
        CGA_ASSERT(dipDetFn_.size()==nondipDetFn_.size());
        CGA_ASSERT(dipDetFn_.size()==somnondipDetFn_.size());
        computeGcCorrectedMeans();
        setUpDelimitedFiles();
        someSomatic_ = testForSomaticData();
        limitFieldNames(fieldNames);
        writer_.reset(new CnvFileVcfRecordWriter(cnvFieldNames_,crr_,dipDet_.size(),someSomatic_));
        ++(*this);
    }

    void CnvFileVcfRecordSource::limitFieldNames(const std::vector< std::string > & fieldNames){
        BOOST_FOREACH(const std::string &f,fieldNames){
            if ( f == "GT" ||
                 f == "CGA_GP" ||
                 f == "CGA_NP" ||
                 f == "CGA_CP" ||
                 f == "CGA_PS" ||
                 f == "CGA_CT" ||
                 f == "CGA_TS" ||
                 f == "CGA_CL" ||
                 f == "CGA_LS" ||
                 ( someSomatic_ && ( f == "CGA_SCL" ||
                                     f == "CGA_SLS" ||
                                     f == "CGA_LAF" ||
                                     f == "CGA_LLAF" ||
                                     f == "CGA_ULAF" ) 
                 ) 
               )
            {
                if ( cnvFieldNamesSet_.find(f) == cnvFieldNamesSet_.end() ){
                    cnvFieldNamesSet_.insert(f);
                    cnvFieldNames_.push_back(f);
                }
            }
        }
    }

    bool CnvFileVcfRecordSource::testForSomaticData(){
        for(size_t ii=0;ii<dipDet_.size();++ii){
            if ( somnondipDet_[ii] != shared_ptr<DelimitedFile>() ){
                return true;
            }
        }
        return false;
    }

    void CnvFileVcfRecordSource::computeGcCorrectedMeans(){
        for(size_t ii=0;ii<dipDet_.size();++ii){
            double cvg;
            double sumCvg = 0;
            int numCvg = 0;
            boost::shared_ptr<std::istream> tmpStrm = 
                InputStream::openCompressedInputStreamByExtension(dipDetFn_[ii]);
            DelimitedFile detDF(*tmpStrm,dipDetFn_[ii]);
            detDF.addField(ValueField<double>("gcCorrectedCvg",&cvg).exception("N", -1));
            while(detDF.next()){
                if (cvg >= 0){
                    numCvg++;
                    sumCvg+=cvg;
                }
            }
            if ( numCvg == 0 ){
                throw Exception("CNV details file " + dipDetFn_[ii] 
                                + " had no non-N GC-corrected coverage values.");
            }
            gcCorrectedMean_[ii] = sumCvg/(double)numCvg;
        }
    }

    void CnvFileVcfRecordSource::setUpDelimitedFiles(){

        for(size_t ii=0;ii<dipDet_.size();++ii){

            dipDetIStr_[ii] = InputStream::openCompressedInputStreamByExtension(dipDetFn_[ii]);
            dipDet_[ii].reset( new DelimitedFile(*(dipDetIStr_[ii]),dipDetFn_[ii]) );
            setupDiploidDetails(dipDet_[ii]);

            nondipDetIStr_[ii] = InputStream::openCompressedInputStreamByExtension(nondipDetFn_[ii]);
            nondipDet_[ii].reset( new DelimitedFile(*(nondipDetIStr_[ii]),nondipDetFn_[ii]));
            setupNondiploidDetails(nondipDet_[ii]);

            if ( "" != somnondipDetFn_[ii] ){
                somnondipDetIStr_[ii]=InputStream::openCompressedInputStreamByExtension(somnondipDetFn_[ii]);
                somnondipDet_[ii].reset( new DelimitedFile(*(somnondipDetIStr_[ii]),somnondipDetFn_[ii]));
                setupSomaticDetails(somnondipDet_[ii]);
            } else {
                somnondipDet_[ii] = shared_ptr<DelimitedFile>();
            }

        }

    }

    void CnvFileVcfRecordSource::setupDiploidDetails(boost::shared_ptr<DelimitedFile> &df){
        df->addField(StringField("chr",&(tmpDip_.chr)));
        df->addField(ValueField<int>("begin",&(tmpDip_.begin)));
        df->addField(ValueField<int>("end",&(tmpDip_.end)));
        df->addField(ValueField<float>("avgNormalizedCvg",&(tmpDip_.avgNormCvg)).exception("N",-1));
        df->addField(ValueField<float>("gcCorrectedCvg",&(tmpDip_.gcCvg)).exception("N", -1));
        df->addField(ValueField<float>("relativeCvg",&(tmpDip_.relCvg)).exception("N", -1));
        df->addField(ValueField<int>("calledPloidy",&(tmpDip_.cP)).exception("N", -1));
        df->addField(StringField("calledCNVType",&(tmpDip_.cT)));
        df->addField(StringField("ploidyScore",&(tmpDip_.pS)));
        df->addField(StringField("CNVTypeScore",&(tmpDip_.tS)));
    }


    void CnvFileVcfRecordSource::setupNondiploidDetails(boost::shared_ptr<DelimitedFile> &df){
        df->addField(StringField("chr",&(tmpNondip_.chr)));
        df->addField(ValueField<int>("begin",&(tmpNondip_.begin)));
        df->addField(ValueField<int>("end",&(tmpNondip_.end)));
        df->addField(ValueField<float>("calledLevel",&(tmpNondip_.cL)).exception("N",-1));
        df->addField(StringField("levelScore",&(tmpNondip_.lS)));
    }

    void CnvFileVcfRecordSource::setupSomaticDetails(boost::shared_ptr<DelimitedFile> &df){
        df->addField(StringField("chr",&(tmpSom_.chr)));
        df->addField(ValueField<int>("begin",&(tmpSom_.begin)));
        df->addField(ValueField<int>("end",&(tmpSom_.end)));
        df->addField(ValueField<float>("calledLevel",&(tmpSom_.cL)).exception("N",-1));
        df->addField(StringField("levelScore",&(tmpSom_.lS)));
        df->addField(StringField("bestLAF",&(tmpSom_.laf)));
        df->addField(StringField("lowLAF",&(tmpSom_.llaf)));
        df->addField(StringField("highLAF",&(tmpSom_.ulaf)));
    }

    cgatools::conv::VcfRecordSource& CnvFileVcfRecordSource::operator++()
    {

        cgatools::reference::Range minRange = genomeEnd_;
        
        for(size_t ii=0;ii<dipDet_.size();++ii)
        {

            if( ! writer_->isDeferred(ii) ){
                
                // For diploid details files, we always advance one line
                if( ! dipDet_[ii]->next() ){
                    // all samples should be at EOF in same time, 
                    // so we should hit this first for the first sample
                    CGA_ASSERT(ii==0);
                    eof_=true;
                    return *this;
                }
                tmpDip_.normedGcCvg = tmpDip_.gcCvg / gcCorrectedMean_[ii];
                writer_->setDiploidData(tmpDip_,ii);
            }
            if ( writer_->getDiploidRange(ii) < minRange ){
                minRange = writer_->getDiploidRange(ii);
            }
        }

        for(size_t ii=0;ii<dipDet_.size();++ii)
        {
            writer_->setDeferred(ii, writer_->getDiploidRange(ii) > minRange );

            if (writer_->isDeferred(ii)){
                continue;
            }

            // For nondiploid details files, where the windows are larger but begin and end 
            // should correspond to the begin and end of some windows, we advance a line 
            // only if the diploid range has gone past the end of the current nondiploid range
            Location dipBeg = writer_->getDiploidRange(ii).beginLocation();
            Location dipEnd = writer_->getDiploidRange(ii).endLocation();
            Location nondipEnd = writer_->getNondiploidRange(ii).endLocation();
            Location nondipBeg = writer_->getNondiploidRange(ii).beginLocation();
            if(dipBeg >= nondipEnd){
                bool isGood = nondipDet_[ii]->next();
                // make sure we did not run out of data .. this should not happen
                CGA_ASSERT(isGood);
                writer_->setNondiploidData(tmpNondip_,ii);
                // test that window beginnings match
                nondipEnd = writer_->getNondiploidRange(ii).endLocation();
                nondipBeg = writer_->getNondiploidRange(ii).beginLocation();
                CGA_ASSERT(dipBeg==nondipBeg);
            }
            CGA_ASSERT( dipBeg >= nondipBeg && dipEnd <= nondipEnd );

            // same rules apply for somatic nondiploid details files
            Location somEnd = writer_->getSomaticRange(ii).endLocation();
            Location somBeg = writer_->getSomaticRange(ii).beginLocation();
            if(someSomatic_ && somnondipDet_[ii] != shared_ptr<DelimitedFile>() ){
                if(dipBeg >= somEnd){
                    bool isGood = somnondipDet_[ii]->next();
                    CGA_ASSERT(isGood);
                    writer_->setSomaticData(tmpSom_,ii);
                    somEnd = writer_->getSomaticRange(ii).endLocation();
                    somBeg = writer_->getSomaticRange(ii).beginLocation();
                    CGA_ASSERT(dipBeg==somBeg);
                }
                CGA_ASSERT( dipBeg >= somBeg && dipEnd <= somEnd );
            }

        }

        return *this;

    }

    std::vector<cgatools::conv::VcfSubFieldHeaderRecord>
    CnvFileVcfRecordSource::getSubFieldHeaderRecords() const
    {

        vector<VcfSubFieldHeaderRecord> result;

        //      '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">\n' +
        //      '##INFO=<ID=CGA_WINEND,Number=1,Type=Integer,Description="End of coverage window">\n' +
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "NS", "1", "Integer",
                             "Number of Samples With Data"));
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_INFO, "CGA_WINEND", "1", "Integer",
                             "End of coverage window"));

        //      '##ALT=<ID=CGA_CNVWIN,Description="Copy number analysis window">\n' +
        result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_ALT, "CGA_CNVWIN", "", "",
                             "Copy number analysis window"));

        // No FILTER records ... but they would go here ...



        // FORMAT records ...
        if ( cnvFieldNamesSet_.find("GT") != cnvFieldNamesSet_.end() )
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"GT","1","String",
                             "Genotype"));
        if (cnvFieldNamesSet_.find("CGA_GP") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_GP","1","Float",
                             "Depth of coverage for 2k window GC normalized to mean"));
        if (cnvFieldNamesSet_.find("CGA_NP") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_NP","1","Float",
                             "Coverage for 2k window, GC-corrected and normalized "
                             "relative to copy-number-corrected multi-sample baseline"));
        if (cnvFieldNamesSet_.find("CGA_CL") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_CL","1","Float",
                             "Nondiploid-model called level"));
        if (cnvFieldNamesSet_.find("CGA_LS") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_LS","1","Integer",
                             "Nondiploid-model called level score"));
        if (cnvFieldNamesSet_.find("CGA_CP") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_CP","1","Integer",
                             "Diploid-model called ploidy"));
        if (cnvFieldNamesSet_.find("CGA_PS") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_PS","1","Integer",
                             "Diploid-model called ploidy score"));
        if (cnvFieldNamesSet_.find("CGA_CT") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_CT","1","String",
                             "Diploid-model CNV type"));
        if (cnvFieldNamesSet_.find("CGA_TS") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_TS","1","Integer",
                             "Diploid-model CNV type score"));

        if (cnvFieldNamesSet_.find("CGA_SCL") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_SCL","1","Float",
                             "Nondiploid-model somatic called level"));
        if (cnvFieldNamesSet_.find("CGA_SLS") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_SLS","1","Integer",
                             "Nondiploid-model somatic called level score"));
        std::string winwidth = nondipDet_[0]->getMetadata().get("WINDOW_WIDTH");
        if (cnvFieldNamesSet_.find("CGA_LAF") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_LAF","1","Float",
                             "Lesser allele fraction estimate, "+winwidth+"bp window"));
        if (cnvFieldNamesSet_.find("CGA_LLAF") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_LLAF","1","Float",
                             "Lesser allele fraction lower bound, "+winwidth+"bp window"));
        if (cnvFieldNamesSet_.find("CGA_ULAF") != cnvFieldNamesSet_.end())
            result.push_back(VcfSubFieldHeaderRecord(
                             VcfSubFieldHeaderRecord::VCF_FORMAT,"CGA_ULAF","1","Float",
                             "Lesser allele fraction upper bound, "+winwidth+"bp window"));

        return result;
    }

    std::string
    CnvFileVcfRecordSource::getSource(size_t idxGenome) const
    {
        return dipDet_[idxGenome]->getMetadata().getSoftwareVersionString();
    }

    std::vector<cgatools::conv::VcfKvHeaderRecord>
    CnvFileVcfRecordSource::getKeyValueHeaderRecords(size_t idxGenome) const
    {
        vector<VcfKvHeaderRecord> result;

        const char* transferKeys[] =
        {
            "GENOME_REFERENCE",
            "MAX_PLOIDY",
            "NUMBER_LEVELS"
        };

        BOOST_FOREACH(const char* key, transferKeys)
        {
            std::string keyStr(key);
            std::string sampleId("");

            if(keyStr == "NUMBER_LEVELS"){
                sampleId = dipDet_[idxGenome]->getMetadata().get("SAMPLE") + ":";
            }

            std::string val("");
            if (dipDet_[idxGenome]->getMetadata().hasKey(key)){
                val = sampleId+dipDet_[idxGenome]->getMetadata().get(key);

                if(keyStr == "MAXIMUM_PLOIDY" && key != dipDet_[0]->getMetadata().get(key)){
                    throw Exception( "CNV VCF conversion requires all "
                                     "genomes to have the same MAXIMUM_PLOIDY");
                }

            } else if (nondipDet_[idxGenome]->getMetadata().hasKey(key)){

                val = nondipDet_[idxGenome]->getMetadata().get(key);

            } else if (somnondipDet_[idxGenome] != shared_ptr<DelimitedFile>() && 
                       somnondipDet_[idxGenome]->getMetadata().hasKey(key)) {

                val = somnondipDet_[idxGenome]->getMetadata().get(key);

            }

            if ( val != "" ){
                result.push_back(VcfKvHeaderRecord(string("source_")+key, sampleId+val));
            }
                    
        }


        std::string extrakey("WINDOW_WIDTH");
        std::string extrakeyrename("NONDIPLOID_WINDOW_WIDTH");
        result.push_back(VcfKvHeaderRecord(string("source_")+extrakeyrename,
                                           nondipDet_[idxGenome]->getMetadata().get(extrakey)));
        // make sure all nondiploid windows are same size ... 
        // otherwise the FORMAT definitions of LAF, ULAF, LLAF will be problematic
        if(nondipDet_[idxGenome]->getMetadata().get(extrakey) != 
           nondipDet_[0]->getMetadata().get(extrakey))
        {
            throw Exception(
                "CNV VCF conversion requires all genomes to have the same nondiploid WINDOW_WIDTH");
        }
        if(somnondipDet_[idxGenome] != shared_ptr<DelimitedFile>() ){
            if(nondipDet_[idxGenome]->getMetadata().get(extrakey) != 
                       somnondipDet_[idxGenome]->getMetadata().get(extrakey))
            {
                throw Exception( "CNV VCF conversion requires normal and somatic calls "
                                 "to have the same nondiploid WINDOW_WIDTH");
            }
        }


        std::string gcmeankey("MEAN_GC_CORRECTED_CVG");
                               //std::string formatter = "%s:%.2f";
        std::string gcmeanval = ((boost::format("%s:%.2f") 
                                  %  dipDet_[idxGenome]->getMetadata().get("SAMPLE") 
                                  % gcCorrectedMean_[idxGenome])).str();
        result.push_back(VcfKvHeaderRecord(string("source_")+gcmeankey,gcmeanval));


        return result;
    }

    std::string
    CnvFileVcfRecordSource::getAssemblyId(size_t idxGenome) const
    {
        CGA_ASSERT(dipDet_[idxGenome]->getMetadata().get("ASSEMBLY_ID") == 
                   nondipDet_[idxGenome]->getMetadata().get("ASSEMBLY_ID"));
        CGA_ASSERT(somnondipDet_[idxGenome] == shared_ptr<DelimitedFile>() ||
                   ( dipDet_[idxGenome]->getMetadata().get("ASSEMBLY_ID") == 
                     somnondipDet_[idxGenome]->getMetadata().get("ASSEMBLY_ID") ) );
        return dipDet_[idxGenome]->getMetadata().get("ASSEMBLY_ID");
    }

    bool CnvFileVcfRecordSource::eof() const
    {
        return eof_;
    }

    const cgatools::conv::VcfRecordWriter& CnvFileVcfRecordSource::operator*() const
    {
        return *writer_;
    }

    const cgatools::conv::VcfRecordWriter* CnvFileVcfRecordSource::operator->() const
    {
        return writer_.get();
    }

} } // cgatools::copynumber
