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
// implied. See the License for the specific language governing0
// permissions and limitations under the License.

#ifndef CGA_TOOLS_JUNCTION_VCF_RECORD_HPP_
#define CGA_TOOLS_JUNCTION_VCF_RECORD_HPP_

//! @file JunctionVcfRecord.hpp

// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/junctions/JunctionCompare.hpp"
#include "cgatools/junctions/JunctionVcfWriter.hpp"
#include "cgatools/conv/VcfRecordSource.hpp"

// Standard library.
#include <string>
#include <vector>
#include <map>
#include <set>

// Boost
#include <boost/shared_ptr.hpp>

// Forward declarations
namespace cgatools { namespace cgdata { class GenomeMetadata; }} 


namespace cgatools { namespace junctions {

///////////////////////////////////////////////////////////////////////

class JunctionVcfRecordWriter : public cgatools::conv::VcfRecordWriter
{
    public: 

        JunctionVcfRecordWriter();

        // Overridden methods of VcfRecordWriter.
        cgatools::reference::Location getLocation() const;

        void writeId     ( std::ostream& out ) const;
        void writeRef    ( std::ostream& out ) const;
        void writeAlt    ( std::ostream& out ) const;
        void writeQual   ( std::ostream& out ) const;
        void writeFilter ( std::ostream& out ) const;
        void writeInfo   ( std::ostream& out ) const;
        void writeFormat ( std::ostream& out ) const;
        void writeSample ( std::ostream& out,  size_t gIdx ) const;

        friend class JunctionVcfRecordSource;

    protected:

        cgatools::reference::Location   pos_;
        std::vector<std::string>        sample_;
        std::string                     id_, 
                                        ref_, 
                                        alt_,
                                        qual_,
                                        filter_,
                                        info_,
                                        format_;
};

///////////////////////////////////////////////////////////////////////

class JunctionVcfRecordSource : public cgatools::conv::VcfRecordSource
{
    public:

        JunctionVcfRecordSource 
        ( 
            const std::vector< boost::shared_ptr<cgatools::cgdata::GenomeMetadata> >& genomes,
            const std::vector<std::string>& junctionFileNames,
            const std::vector<std::string> fieldNames,
            const cgatools::reference::CrrFile& crr,
            size_t scoreThreshold = 10,
            size_t sideLengthThreshold_ = 70,
            size_t distanceTolerance = 200,
            size_t junctionLengthThreshold = 500,
            bool   normalPriorityOutput = false,
            bool   useHighConfidenceJunctionsForTumor = true
        );

        // Overridden methods of VcfRecordSource.

        std::vector<cgatools::conv::VcfSubFieldHeaderRecord> getSubFieldHeaderRecords() const;
        std::string getSource(size_t idxGenome) const;
        std::vector<cgatools::conv::VcfKvHeaderRecord> getKeyValueHeaderRecords
            ( size_t idxGenome) const;
        std::string getAssemblyId(size_t idxGenome) const;

        bool eof() const;
        cgatools::conv::VcfRecordSource& operator++();
        const cgatools::conv::VcfRecordWriter& operator*() const;
        const cgatools::conv::VcfRecordWriter* operator->() const;

        // JunctionVcf-specific methods

        /// tells the writer whether all eligible fields should be written 
        /// even if they are not in the list of field names
        void writeAllFields ( bool v );

        int run();

    protected:

        int run1Genome();
        int run2Genomes();

        bool add 
            (   
                std::vector<cgatools::conv::VcfSubFieldHeaderRecord>& result, 
                cgatools::conv::VcfSubFieldHeaderRecord::Key key,
                const std::string& id, 
                const std::string& number, 
                const std::string& type, 
                const std::string& description
            )   const;

        JunctionVcfRecordWriter createRecord
            (   
                const JunctionRef& jref,
                size_t side,         
                const JunctionCompatMapPerFile& compat 
            )   const;

        std::string getAltField
            ( 
                const JunctionRef& jref, 
                size_t side, 
                bool suppressChrom = false
            )  const;

        std::string getPosition ( reference::Location pos, const std::string& sep )   const;

        std::string getInfo  ( const JunctionRef& jref, size_t side )   const;
        std::string getMEI   ( const std::string& med ) const;
        std::string getId    ( const JunctionRef& jref, size_t side) const;
        std::string getFormat( const JunctionRef& jref, size_t side) const;
        
        std::string getSample      
            ( 
                const JunctionRef& jref, 
                size_t side, 
                size_t idx, 
                const JunctionCompatMapPerFile& compat  
            )   const;

        std::string getSample      ( const JunctionRef& jref, size_t side ) const;
        std::string getSampleFilter( const JunctionRef& jref ) const;
        std::string addFilterFlag  ( const std::string& flag, bool& filtered) const;

        bool need ( const std::string& fieldName ) const;

        void copyJunctionListForVcf
            (
                const junctions::JunctionRefs& jrl,
                std::vector<JunctionRefSide>& out,
                JunctionCompatMapPerFile& compat
            ) const;

        void pickNormalPriorityMatch 
            ( 
                const JunctionRef& jr, 
                JunctionRefSet& junctionsToSuppress 
            )   const;

        void pickDefaultMatch 
            ( 
            const JunctionRef& jr, 
            JunctionRefSet& junctionsToSuppress 
            )   const;

        /// reference genome
        const cgatools::reference::CrrFile& crr_;

        /// Input genomes' sample IDs
        std::vector<std::string> sampleIds_;
        
        /// Junction file mames
        std::vector<std::string> junctionFileNames_;

        /// Field names
        std::vector<std::string> fieldNames_;
        std::set<std::string>    fieldNameSet_;

        /// Flag that enables writing of all eligible fields even if they are not in the list of field names
        bool writeAllFields_;

        /// PASS thresholds for discordant DNBs count
        size_t scoreThreshold_;

        /// PASS threshold for minimum junction side length, in base pairs
        size_t sideLengthThreshold_;

        /// Length threshold for junction compatibility
        size_t junctionLengthThreshold_;

        /// Distance tolerance for junction compatibility
        size_t distanceTolerance_;

        /// Normal (non-tumor) junction priority for the output
        bool normalPriorityOutput_;

        bool useHighConfidenceJunctionsForTumor_;

        /// Loaded junction files
        JunctionFiles junctionFiles_;

        /// VCF records
        std::vector<JunctionVcfRecordWriter> records_;

        /// Current VCF record
        size_t currentRecord_;

        static const char* sampleFieldIDs_[];
};

}} // end namespace
#endif //CGA_TOOLS_JUNCTION_VCF_SOURCE_HPP_
