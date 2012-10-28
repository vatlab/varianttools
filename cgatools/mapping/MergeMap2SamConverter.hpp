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

#ifndef CGA_TOOLS_COMMAND_MERGEMAP2SAM_CONVERTER_HPP_
#define CGA_TOOLS_COMMAND_MERGEMAP2SAM_CONVERTER_HPP_ 1

//! @file MergeMap2SamConverter.hpp

#include "cgatools/core.hpp"
#include "cgatools/util/RangeSet.hpp"
#include "Map2SamConverter.hpp"
#include "EvidenceCache.hpp"

#include <boost/array.hpp>

namespace cgatools { namespace mapping {

    class DoubleMatrix;

    class MergedMap2SamConfig : public Map2SamConfig 
    {
    public:
        MergedMap2SamConfig()
            :   Map2SamConfig()
            ,mapqAlpha_(1E-12)
            ,minMateGapFrequency_(1E-7)
            ,minMapQ_(0)
            ,debugOutStream_(&std::cerr)
        {}

        double  mapqAlpha_;
        double  minMateGapFrequency_;

        uint8_t     minMapQ_;
        std::string evidenceCacheRoot_;
        std::ostream* debugOutStream_;
    };

    class MergedMap2SamConverter : public Map2SamConverter 
    {
    public:
        typedef boost::array<std::vector<size_t>,2> MappingIndicesBySide;
        typedef boost::array<std::vector<double>,2> MappingWeightsBySide;

        MergedMap2SamConverter(const MergedMap2SamConfig &config, std::ostream &outSamFile)
            : Map2SamConverter(config, outSamFile), config_(config)
        {}

        virtual void init();

    protected:
        //! override the method to merge the evidence mappings
        virtual bool processMappings(const mapping::ReadsRecord& readsRecord, 
            SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const;

        //! Export mapping record in SAM format
        virtual void writeMappingRecord(const SamRecord &m) const;

        //! Merge evidence records into the base mappings
        void addEvidenceRecords(const mapping::ReadsRecord& readsRecord, SamRecordArray& records) const;
        //! Process merged mappings: recompute MAPQ, deduplicate, select best mapping pairs
        void regroupMappings(const mapping::ReadsRecord& readsRecord, SamRecordArray& records) const;

        //! computes P(DNBs) * P(gaps)
        void computeMappingWeights(const SamRecordArray& records, 
            const MappingIndicesBySide& mappings, MappingWeightsBySide& weights) const;
        //! filter the mappings that have equal position and alignment, 
        //! leave Evidence mappings if equal to base
        void filterDuplicatedMappings(SamRecordArray& records) const;

        void computeWeightMatrix( 
            const SamRecordArray& records, const MappingIndicesBySide& mappings,
            const MappingWeightsBySide& weights, DoubleMatrix& weightMatrix, DoubleMatrix& mateGaps) const;

        double computeTotalWeight(const DoubleMatrix& weightMatrix, MappingWeightsBySide& sumAllPairs) const;

        void findBestCandidates(const DoubleMatrix& m, const DoubleMatrix& mateGaps,
            const MappingIndicesBySide& mappingsBySide, const MappingWeightsBySide& sumAllPairs, 
            MappingIndicesBySide& bestGroupMappings, MappingIndicesBySide& groups) const;

        //! compute weights of double arm mappings and detect primary mappings for each group
        void reweightSingleMappings(SamRecordArray& records, const MappingIndicesBySide& mappingsBySide, 
            const MappingWeightsBySide& dnbArmWeights,const MappingWeightsBySide& sumAllPairs) const;

        //! compute weights of single arm mappings
        void reweightDoubleArmMappings(SamRecordArray& records, const MappingIndicesBySide& mappingsBySide, 
            const DoubleMatrix& weightMatrix, const DoubleMatrix& mateGapFrequency,
            const MappingWeightsBySide& sumAllPairs, const MappingIndicesBySide& bestFullMappedArms,
            const MappingIndicesBySide& groups, const double &totalWeight ) const;

        //! the function selects the best mapping record based on MAPQ set by reweightDoubleArmMappings
        void detectTheBestMappingAndMakeItPrimary(SamRecordArray& records) const;

        //! output debug information to stdout
        void dumpBestLocalMatches(const mapping::ReadsRecord& readsRecord, const SamRecordArray& records, 
            const MappingIndicesBySide& mappingsBySide, const MappingWeightsBySide& dnbArmWeights,
            const DoubleMatrix& weights, const DoubleMatrix& mateGaps,
            const MappingWeightsBySide& sumAllPairs, const MappingIndicesBySide& bestFullMappedArms,
            const MappingIndicesBySide& groups,
            const double& totalWeight) const;

        void printMatrix(const DoubleMatrix& m, std::ostream& ostr) const;

        const MergedMap2SamConfig & config_;

        mapping::BatchRecords   evidenceBatchRecords_;
    };

} } // cgatools::mapping

#endif // CGA_TOOLS_COMMAND_MERGEMAP2SAM_CONVERTER_HPP_
