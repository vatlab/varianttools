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
#include "MergeMap2SamConverter.hpp"
#include "AlleleAlignment.hpp"
#include "LibraryData.hpp"
#include "EvidenceSamUtil.hpp"

#include <boost/date_time.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <iomanip>
#include <cmath>
#include <numeric>


namespace cgatools { namespace mapping {

    void MergedMap2SamConverter::addEvidenceRecords( 
        const mapping::ReadsRecord& readsRecord, SamRecordArray& records ) const
    {
        SamRecordArray evidenceRecords;
        CGA_ASSERT(!config_.evidenceCacheRoot_.empty());

        //int totalEvidenceMappings = evidenceBatchRecords_.count(readsRecord.recordIndex_);
        std::pair<mapping::BatchRecords::const_iterator,mapping::BatchRecords::const_iterator> pii
            = evidenceBatchRecords_.equal_range(readsRecord.recordIndex_);
        for (mapping::BatchRecords::const_iterator it=pii.first; it!=pii.second; ++it)
        {
            const EvidenceCacheDnbRecord &e = it->second;

            const cgdata::HalfDnbStructure &hDnbStruct = 
                library_->dnbStructure_.halfDnbs_[e.side_];

            std::auto_ptr<EvidenceSamRecord> 
                erOwner(new EvidenceSamRecord(e,readsRecord.reads_,readsRecord.scores_
                                                ,false,0,hDnbStruct.totReadLength_,reference_)
                );
            EvidenceSamRecord *er = erOwner.get();
            if (!er->correctPosition(reference_))
                continue;
            
            //deduplicate exact copies
            for (SamRecordArray::iterator findIt = evidenceRecords.begin();;++findIt)
            {
                if (findIt==evidenceRecords.end())
                {
                    evidenceRecords.push_back(erOwner);
                    break;
                }
                if (findIt->side_==er->side_
                    && findIt->onNegativeStrand_==er->onNegativeStrand_
                    && findIt->chr_==er->chr_ 
                    && findIt->position_==er->position_
                    && findIt->extCigar_==er->extCigar_) 
                {
                    if (static_cast<const EvidenceSamRecord &>(*findIt).armWeight_ < er->armWeight_)
                    {
                        if (config_.dumpDebugInfo_)
                            *config_.debugOutStream_ << "Remove duplicate: " 
                                << static_cast<const EvidenceSamRecord &>(*findIt) << std::endl;
                        findIt = evidenceRecords.erase(findIt);
                        evidenceRecords.insert(findIt,erOwner);
                    } else
                    {
                        if (config_.dumpDebugInfo_)
                            *config_.debugOutStream_ << "Remove duplicate: " << *er << std::endl;
                    }

                    break;
                }
            }

#ifdef TRACE_MAPQ_COMPUTATION
            std::auto_ptr<EvidenceSamRecord> 
                er1(new EvidenceSamRecord(e,false,1,hDnbStruct.totReadLength_,reference_));

            size_t mateGapSize = abs(e.offsetInReference_[0]-e.offsetInReference_[1]);
            double mateGapFrequency = 
                library_->mateGapTable_->getFrequency(mateGapSize);

            double mapWeight0 = er->armWeight_;

            AlleleSequenceAligner mateAlleleAligner(reference_);
            mateAlleleAligner.setInterval(er1->chr_, 0, 0, NULL);

            GapProbabilityAndConcordanceExtractor 
                mateConcordanceExtractor(mateAlleleAligner, 
                library_->gapsEstimators_);
            mateConcordanceExtractor.run(1-e.side_, e.strand_, e.offsetInReference_[1], 
                e.referenceAlignment_[1], e.sequence_, e.scores_);

            double mapWeight1 = mateConcordanceExtractor.concordance_ * 
                                        mateConcordanceExtractor.gapProbability_;


            double fullMapWeight = mapWeight0*mapWeight1*mateGapFrequency;
            double alpha = 1E-9;
            double singleMappingMapQ0 = alpha/(alpha+fullMapWeight);
            double mapq = -10*log10(singleMappingMapQ0);
            int mapqOrig = e.mappingQuality_-33;
            CGA_ASSERT(fabs(mapq-mapqOrig)<100 && totalEvidenceMappings>=0);//DEBUG

            temporary:
            records.push_back(er1);
            er1->mates_.push_back(er);
            er->mates_.push_back(er1);
#endif //TRACE_MAPQ_COMPUTATION
        }
        records.transfer(records.end(),evidenceRecords);
    }


    namespace {

        template <class T>
        void sumDimension(T begin, T end, std::vector<double>& result)
        {
            result.clear();
            result.reserve(end-begin);
            for(;begin!=end;++begin)
                result.push_back(std::accumulate(begin.begin(),begin.end(),0.0));
        }

        class CompareSamRecordsByStrandAndLocation : public std::binary_function<SamRecord, SamRecord, bool>
        {
        public:
            bool operator()(const SamRecord& r0,const SamRecord& r1) const {
                if (r0.onNegativeStrand_!=r1.onNegativeStrand_)
                    return r0.onNegativeStrand_ < r1.onNegativeStrand_;
                if (r0.chr_!=r1.chr_)
                    return r0.chr_ < r1.chr_;
                return r0.position_ < r1.position_;
            }
        };

        struct IsDefaultMapping : public std::unary_function<const SamRecord&, bool> {
            bool operator() (const SamRecord& r) const {return r.typeId_==SamRecord::DEFAULT;}
        };

    }


    class DoubleMatrix : public boost::numeric::ublas::matrix<double> {};


    void MergedMap2SamConverter::computeMappingWeights( 
        const SamRecordArray& records, 
        const MappingIndicesBySide& mappings,
        MappingWeightsBySide& weights
        ) const
    {
        AlleleSequenceAligner mateAlleleAligner(reference_);
        GapProbabilityAndConcordanceExtractor 
            concordanceExtractor(mateAlleleAligner,library_->gapsEstimators_);

        //compute PDNBs * Pgaps
        for(size_t side=0; side<mappings.size(); ++side)
        {
            size_t armDnbCount = mappings[side].size();
            weights[side].resize(armDnbCount);
            for(size_t j=0; j<armDnbCount; ++j)
            {
                const SamRecord &r = records[mappings[side][j]];
                if (r.typeId_ != SamRecord::EVIDENCE_CACHE)
                {
                    mateAlleleAligner.setInterval(r.chr_, 0, 0, NULL);
                    concordanceExtractor.run(r.side_, r.onNegativeStrand_, r.position_, 
                        r.extCigar_, r.fullReadSequence_, r.fullReadScores_);
                    weights[side][j] = 
                        concordanceExtractor.concordance_ * concordanceExtractor.gapProbability_;
                } else
                    weights[side][j] = static_cast<const EvidenceSamRecord &>(r).armWeight_;
            }
        }
    }


    void MergedMap2SamConverter::filterDuplicatedMappings(SamRecordArray& records) const
    {
        MappingIndicesBySide mappings;
        for (size_t i=0; i<records.size(); ++i)
            mappings[records[i].side_].push_back(i);
        //remove base mappings equal to the evidence mappings
        for(size_t side=0; side<mappings.size(); ++side)
        {
            if (mappings[side].size()<2)
                continue;
            for(int j=mappings[side].size()-1; j>0; --j)
            {
                const SamRecord& r0 = records[mappings[side][j-1]];
                const SamRecord& r1 = records[mappings[side][j]];
                if (r0.onNegativeStrand_==r1.onNegativeStrand_
                    && r0.chr_==r1.chr_ 
                    && r0.position_==r1.position_
                    && Cigar(r0.extCigar_,true,true)==Cigar(r1.extCigar_,true,true))
                {
                    int ind = r0.typeId_==SamRecord::BASE_MAPPING ? j-1 : j;
                    records[mappings[side][ind]].typeId_ = SamRecord::DEFAULT;
                }   
            }
        }
        records.erase_if(IsDefaultMapping());
    }


    void MergedMap2SamConverter::computeWeightMatrix( 
        const SamRecordArray& records, 
        const MappingIndicesBySide& mappings,
        const MappingWeightsBySide& weights,
        DoubleMatrix& weightMatrix,
        DoubleMatrix& mateGaps
        ) const
    {
        size_t size0 = mappings[0].size();
        size_t size1 = mappings[1].size();
        weightMatrix.resize(size0,size1);
        mateGaps.resize(size0,size1);
        // compute weight matrix m(X,Y)=weight of the pair 
        // where X-is a mapping of side 0, Y is a mapping of side 1
        for (size_t i=0; i<size0; ++i)
            for (size_t j=0; j<size1; ++j)
            {
                const SamRecord &r0 = records[mappings[0][i]];
                const SamRecord &r1 = records[mappings[1][j]];
                reference::Location l0(r0.chr_, r0.position_);
                reference::Location l1(r1.chr_, r1.position_);

                double mateGapFrequency = -1;
                if (r0.onNegativeStrand_==r1.onNegativeStrand_ && r0.chr_==r1.chr_)
                {
                    int32_t mateGapSize = r0.onNegativeStrand_ ?
                        r0.position_ - r1.position_ - Cigar(r1.extCigar_).getReferenceLength()
                        :   r1.position_ - r0.position_ - Cigar(r0.extCigar_).getReferenceLength();
                    if (reference_.listChromosomes()[r0.chr_].isCircular())
                    {
                        if (mateGapSize < 0)
                            mateGapSize = -mateGapSize;
                        mateGapSize = std::min(mateGapSize,
                            int32_t(reference_.listChromosomes()[r0.chr_].length()) - mateGapSize);
                    }
                    mateGapFrequency = library_->mateGapTable_->getFrequency(mateGapSize,-1);
                }

                if (mateGapFrequency < -0.5)
                    weightMatrix(i,j) = 0;
                else
                    weightMatrix(i,j) = weights[0][i] * weights[1][j] * mateGapFrequency;

                mateGaps(i,j) = mateGapFrequency;
            }
    }


    double MergedMap2SamConverter::computeTotalWeight(
        const DoubleMatrix& weightMatrix, MappingWeightsBySide& sumAllPairs) const 
    {
        // find sum of the all full mappings related to a particular side mapping
        sumDimension(weightMatrix.begin1(),weightMatrix.end1(), sumAllPairs[0]);
        sumDimension(weightMatrix.begin2(),weightMatrix.end2(), sumAllPairs[1]);

        double totalWeight0 = std::accumulate(sumAllPairs[0].begin(),sumAllPairs[0].end(),0.0);
        double totalWeight1 = std::accumulate(sumAllPairs[1].begin(),sumAllPairs[1].end(),0.0);
        CGA_ASSERT_MSG(fabs(totalWeight0-totalWeight1)<1E-8,CGA_VOUT(totalWeight0)<<CGA_VOUT(totalWeight1));

        if (!(sumAllPairs[0].empty()||sumAllPairs[0].empty()))
        {
            CGA_ASSERT_EQ(sumAllPairs[0].size(),weightMatrix.size1());
            CGA_ASSERT_EQ(sumAllPairs[1].size(),weightMatrix.size2());
        }
        return totalWeight0;
    }


    void MergedMap2SamConverter::dumpBestLocalMatches( 
        const mapping::ReadsRecord& readsRecord,
        const SamRecordArray& records, 
        const MappingIndicesBySide& mappingsBySide,
        const MappingWeightsBySide& dnbArmWeights,
        const DoubleMatrix& weights, 
        const DoubleMatrix& mateGaps,
        const MappingWeightsBySide& sumAllPairs,
        const MappingIndicesBySide& bestFullMappedArms,
        const MappingIndicesBySide& groups,
        const double& totalWeight
        ) const
    {
        std::ostream &ostr = *config_.debugOutStream_;
        ostr << std::endl 
            << "#####: " << readsRecord.reads_ << '\t' << readsRecord.scores_ << std::endl;
        for (size_t side=0; side<mappingsBySide.size(); ++side)
        {
            ostr << "side: " << side;
            if (!bestFullMappedArms[side].empty())
            {
                ostr << " Best full mapping arms(" << bestFullMappedArms[side].size() << "): ";
                for (size_t i=0; i<bestFullMappedArms[side].size(); ++i){
                    double weight = sumAllPairs[side][bestFullMappedArms[side][i]];
                    int weightDb = int(-10*log10(weight));
                    double weightPhred = 1-(weight/
                        (config_.mapqAlpha_+totalWeight));
                    double weightPhredDb = int(-10*log10(weightPhred));
                    ostr << bestFullMappedArms[side][i]
                    << ' ' << weightDb << '=' << weight 
                        << "," << weightPhredDb << '=' << weightPhred 
                        << "/ ";
                }
            } else {
                ostr << " no full mapping found:";
            }
            ostr << std::endl;
            for (size_t i=0; i<mappingsBySide[side].size(); ++i)
                ostr 
                << "  " << i 
                << " (" << dnbArmWeights[side][i] << ")"
                << ": " << records[mappingsBySide[side][i]] << std::endl;
        }
        size_t lastGroupNo = groups[0].size();
        if (lastGroupNo>1){
            ostr << "Found multiple groups: " << lastGroupNo << std::endl;
            ostr << "Weights: ";
            printMatrix(weights,ostr);
            ostr << "Mate gaps: ";
            printMatrix(mateGaps,ostr);
            ostr << "Best mappings:";
            for(size_t i=0; i<lastGroupNo; ++i)
                ostr << "(" << bestFullMappedArms[0][i] << "," << bestFullMappedArms[1][i] << ")";
            ostr << std::endl;
            ostr << "Groups:";
            for(size_t i=0; i<lastGroupNo; ++i)
                ostr << "(" << groups[0][i] << "," << groups[1][i] << ")";
            ostr << std::endl;
            ostr << std::endl;
        }
    }


    void MergedMap2SamConverter::printMatrix(const DoubleMatrix& m, std::ostream& ostr) const
    {
        size_t width = ostr.width();
        std::ios_base::fmtflags flags = ostr.flags();
        ostr << "[" << m.size1() << "," << m.size2() << "]" << std::endl;
        ostr << "     ";
        for (size_t j=0; j<m.size2(); ++j)
            ostr << std::setiosflags(std::ios_base::right) << std::setw(15) << j;
        ostr << std::endl;
        for (size_t i=0; i<m.size1(); ++i)
        {
            ostr << std::setiosflags(std::ios_base::right) << std::setw(3) << i << " (";
            for (size_t j=0; j<m.size2(); ++j)
                ostr << std::setiosflags(std::ios_base::right) << std::setw(15) << m(i,j);
            ostr << ")" << std::endl;
        }
        ostr << std::setw(width) << std::setiosflags(flags);
    }


    void MergedMap2SamConverter::findBestCandidates(
        const DoubleMatrix& m,
        const DoubleMatrix& mateGaps,
        const MappingIndicesBySide& mappingsBySide,
        const MappingWeightsBySide& sumAllPairs,
        MappingIndicesBySide& bestGroupMappings,
        MappingIndicesBySide& groups
    ) const
    {
        size_t maxGroups = std::min(m.size1(),m.size2());
        size_t lastGroupNo = 0;
        double maxGroupWeight = -1;

        for (size_t i=0; i<groups.size(); ++i)
        {
            groups[i].resize(maxGroups);
            bestGroupMappings[i].resize(maxGroups);
        }

        //find maximum index halfdnb index for each column to define the group border
        std::vector<size_t> maxColumnValue(m.size2(),0);
        for (size_t j=0; j<m.size2(); ++j)
        {
            if (j>0)
                maxColumnValue[j] = maxColumnValue[j-1];
            for (size_t i=maxColumnValue[j]+1; i<m.size1(); ++i)
                if (mateGaps(i,j)>-0.5 && i>maxColumnValue[j])
                    maxColumnValue[j] = i;
        }

        for (size_t i=0; i<m.size1(); ++i)
            for (size_t j=0; j<m.size2(); ++j)
            {
                if (mateGaps(i,j)>-0.5) //if a pair exists (compensate double imprecision)
                {
                    if (i<1 || j<1) {
                        lastGroupNo = 0;
                    } else 
                        if (i>maxColumnValue[j-1] && j>groups[1][lastGroupNo]) 
                        { //allow only quasi-diagonal grouping in a 2d matrix 
                          //(it is possible while the halfdnbs are sorted)
                            ++lastGroupNo;
                            CGA_ASSERT_L(lastGroupNo,maxGroups);
                            maxGroupWeight = -1;
                        }

                    groups[0][lastGroupNo]=i;
                    groups[1][lastGroupNo]=j;
                    double currentWeight = sumAllPairs[0][i]*sumAllPairs[1][j]*mateGaps(i,j);

                    if (currentWeight>maxGroupWeight)
                    {
                        bestGroupMappings[0][lastGroupNo]=i;
                        bestGroupMappings[1][lastGroupNo]=j;
                        maxGroupWeight = currentWeight;
                    }
                }
            }
        ++lastGroupNo;
        for (size_t i=0; i<groups.size(); ++i)
        {
            groups[i].resize(lastGroupNo);
            bestGroupMappings[i].resize(lastGroupNo);
        }
    }


    void MergedMap2SamConverter::reweightSingleMappings(SamRecordArray& records,
        const MappingIndicesBySide& mappingsBySide, 
        const MappingWeightsBySide& dnbArmWeights,
        const MappingWeightsBySide& sumAllPairs
        ) const
    {
        //Compute total sum of single arm mapping weight
        double singleArmWeightSum = 0;
        for (size_t side=0; side<mappingsBySide.size(); ++side) 
            for (size_t i=0; i<mappingsBySide[side].size(); ++i)
                singleArmWeightSum += std::pow(dnbArmWeights[side][i],2)*config_.minMateGapFrequency_;

        for (size_t side=0; side<mappingsBySide.size(); ++side) //for side = side0, side1
            for (size_t i=0; i<mappingsBySide[side].size(); ++i)
                if (sumAllPairs[side].empty() || sumAllPairs[side][i]==0) // Single arm mapping
                {
                    SamRecord& r = records[mappingsBySide[side][i]];
                    r.isPrimary_ = false;
                    double mappingWeight = std::pow(dnbArmWeights[side][i],2)*config_.minMateGapFrequency_;
                    r.setMappingQuality(uint8_t(
                        -10*log10(1-(mappingWeight/(singleArmWeightSum+config_.mapqAlpha_)))),
                        false);
                }
    }


    //! compute the weight of double arm mappings and detect primary mappings for each group
    void MergedMap2SamConverter::reweightDoubleArmMappings(SamRecordArray& records,
        const MappingIndicesBySide& mappingsBySide, 
        const DoubleMatrix& weightMatrix,
        const DoubleMatrix& mateGapFrequency,
        const MappingWeightsBySide& sumAllPairs,
        const MappingIndicesBySide& bestFullMappedArms,
        const MappingIndicesBySide& groups,
        const double &totalWeight
        ) const
    {
        boost::array<size_t,2> currentInd;
        currentInd.assign(0);

        // for all the found groups recompute MAPQ
        // and select primary and non-primary mappings
        for (size_t group=0; group<groups[0].size(); ++group)
        {
            //for each side => side 0, side 1
            for (size_t side=0; side<currentInd.size(); ++side)
            {
                //for each side mapping find the best mate pair and use it as weight
                //move the best mate to the first position in the mate list 
                for(size_t& ind = currentInd[side]; 
                            ind<=groups[side][group]; ++ind)//within the current group
                {
                    SamRecord& r = records[mappingsBySide[side][ind]];
                    double mappingWeight = 0;
                    size_t best_mate = 0;

                    // set single arm mappings to 0
                    if (sumAllPairs[side].empty() || sumAllPairs[side][ind] == 0) // Single arm mapping
                    {
                        r.setMappingQuality(0, false);
                        continue;
                    }

                    // process double mappings
                    r.mates_.clear();

                    //the current mapping is the best
                    r.isGroupPrimary_ = (ind == bestFullMappedArms[side][group]); 

                    size_t prev_ind = group==0 ? 0 : groups[1-side][group-1]+1;

                    //find all mates
                    for (size_t mateInd=prev_ind; mateInd<=groups[1-side][group]; ++mateInd)
                    {

                        double weight = side==0 ? weightMatrix(ind,mateInd) : weightMatrix(mateInd,ind);
                        double mateGapFrq = side==0 ? mateGapFrequency(ind,mateInd) 
                                                                        : mateGapFrequency(mateInd,ind);

                        if (mateGapFrq<-0.5)
                            continue;

                        r.mates_.push_back(&records[mappingsBySide[1-side][mateInd]]);

                        //for the group primary mappings the best mate is already found
                        //for the non-primary mappings the best mate is the one with the highest score
                        if ((r.isGroupPrimary_ && mateInd == bestFullMappedArms[1-side][group])
                            || (weight>mappingWeight && !r.isGroupPrimary_))
                        {
                            best_mate = r.mates_.size()-1;
                            mappingWeight = weight;
                            //set a mate to group primary if the mapping is a group primary
                            //the consistency is guaranteed by the grouping algorithm
                            r.mates_.back()->isGroupPrimary_ = r.isGroupPrimary_; 
                        }

                    }

                    r.setMappingQuality(uint8_t(
                        -10*log10(1-(mappingWeight/(totalWeight + config_.mapqAlpha_)))),
                        true);

                    if (r.mates_.size()>1)
                        std::swap(r.mates_[0],r.mates_[best_mate]);
                }
            }
        }
    }


    void MergedMap2SamConverter::regroupMappings( 
        const mapping::ReadsRecord& readsRecord, SamRecordArray& records ) const
    {
        if (records.empty())
            return;

        records.sort(CompareSamRecordsByStrandAndLocation());

        filterDuplicatedMappings(records);

        //split mappings by side
        MappingIndicesBySide mappingsBySide;
        for (size_t i=0; i<records.size(); ++i)
            mappingsBySide[records[i].side_].push_back(i);

        MappingWeightsBySide dnbArmWeights;
        computeMappingWeights(records, mappingsBySide, dnbArmWeights);

        DoubleMatrix weightMatrix;
        DoubleMatrix mateGapFrequency;
        computeWeightMatrix(records, mappingsBySide, dnbArmWeights, weightMatrix, mateGapFrequency);

        // find sum of the all full mappings related to a particular side mapping
        MappingWeightsBySide sumAllPairs;
        double totalWeight = computeTotalWeight(weightMatrix, sumAllPairs);

        MappingIndicesBySide bestFullMappedArms;
        MappingIndicesBySide groups;

        if (config_.dumpDebugInfo_)
        {
            dumpBestLocalMatches(readsRecord, records, mappingsBySide, dnbArmWeights, 
                weightMatrix, mateGapFrequency, sumAllPairs, bestFullMappedArms, groups, totalWeight);
        }

        if (mateGapFrequency.size1()>0 && mateGapFrequency.size2()>0 && totalWeight>0)
        {
            //handle a common case: one full consistent mapping
            if (mateGapFrequency.size1()==1 && mateGapFrequency.size2()==1) 
            {
                if (mateGapFrequency(0,0)>-0.5)
                {
                    bestFullMappedArms[0].push_back(0);
                    bestFullMappedArms[1].push_back(0);
                    groups = bestFullMappedArms;
                }
            } else //handle all the other cases
            {
                findBestCandidates(weightMatrix, mateGapFrequency, 
                                mappingsBySide, sumAllPairs, bestFullMappedArms,groups);
            }
        }

        CGA_ASSERT_EQ(groups[0].empty(), (totalWeight==0.0) );

        if (totalWeight > 0) //contains double mappings
        {
            //compute MAPQ for double arm mappings, set MAPQ to 0 for single arm mappings
            reweightDoubleArmMappings(records, mappingsBySide, weightMatrix, mateGapFrequency,
                                        sumAllPairs, bestFullMappedArms, groups, totalWeight);
        } else {
            //compute MAPQ for single arm mappings. Set mapping weight as 
            //"default mate gap weight" * sqr("arm mapping weight") to make it comparable to double arm MAPQs
            reweightSingleMappings(records,mappingsBySide,dnbArmWeights,sumAllPairs);
        }


        if (config_.dumpDebugInfo_)
        {
            dumpBestLocalMatches(readsRecord, records, mappingsBySide, dnbArmWeights, 
                weightMatrix, mateGapFrequency, sumAllPairs, bestFullMappedArms, groups, totalWeight);
        }

    }

    namespace {
        struct BestHalfMappingInfo {
            BestHalfMappingInfo() : bestMapInd_(-1), bestMapQ_(-1), isConsistent_(false) {}
            size_t  bestMapInd_;
            int     bestMapQ_;
            bool    isConsistent_;
        };
    }

    void MergedMap2SamConverter::detectTheBestMappingAndMakeItPrimary(SamRecordArray& records) const
    {
        size_t bestSide = -1;
        boost::array<BestHalfMappingInfo, Map2SamConverter::MAX_SIDES> bestHalfMappings;

        //find the best mappings, reset primary flags
        for (size_t i=0; i<records.size(); ++i)
        {
            SamRecord& r = records[i];
            bool isConsistent = !r.mates_.empty();
            BestHalfMappingInfo & info = bestHalfMappings[r.side_];
            if (int(r.getMappingQuality()) > info.bestMapQ_ && isConsistent >= info.isConsistent_)
            {
                info.bestMapQ_ = r.getMappingQuality();
                info.bestMapInd_ = i;
                info.isConsistent_ = isConsistent;
                bestSide = r.side_;
            }
            r.isPrimary_ = false;
        }

        CGA_ASSERT_NEQ(size_t(-1),bestSide);
        CGA_ASSERT_EQ(Map2SamConverter::MAX_SIDES,2);

        //set the best mapping as primary
        SamRecord& bestR = records[bestHalfMappings[bestSide].bestMapInd_];
        bestR.isPrimary_ = true;

        if (!bestR.mates_.empty()) //if the best mapping has a mate - make it primary
        {
            bestR.mates_[0]->isPrimary_ = true;

            //connect the best mate back
            SamRecord::SamRecords &mateMates = bestR.mates_[0]->mates_;
            for (size_t i=0; i<mateMates.size(); ++i)
                if (mateMates[i]==&bestR)
                {
                    if (i>0)
                        std::swap(mateMates[i],mateMates[0]);
                    break;
                }

            CGA_ASSERT_MSG(bestR.mates_[0]->mates_[0]==&bestR, "Best mate's best mate is not this record");
        } 
        //else  //otherwise take the best mapping of the mate
        //{
        //    size_t  bestOtherMappingInd = bestHalfMappings[1-bestSide].bestMapInd_;
        //    //CGA_ASSERT_L(bestOtherMappingInd,records.size()); //no primary mapping found for the arm
        //    if (bestOtherMappingInd<records.size()) //primary mapping found for the arm
        //        records[bestOtherMappingInd].isPrimary_ = true;
        //}
    }


    bool MergedMap2SamConverter::processMappings(const mapping::ReadsRecord& readsRecord, 
        SamRecordArray& samMappings, const mapping::MappingsRecords& baseMappingRecords) const
    {
        convertBaseMappingsIntoSamMappings(readsRecord, samMappings, baseMappingRecords);

        CGA_ASSERT(!config_.evidenceCacheRoot_.empty());
        addEvidenceRecords(readsRecord, samMappings);
        regroupMappings(readsRecord, samMappings);

        detectTheBestMappingAndMakeItPrimary(samMappings);
        return true;
    }


    void MergedMap2SamConverter::init()
    {
        Map2SamConverter::init();
        cgdata::GenomeMetadata gm(config_.exportRootDirName_);
        library_.reset(new LibraryData(gm.getLibraryMetadata(gm.getLibraryName(laneId_)),false));

        if (!config_.evidenceCacheRoot_.empty())
        {
            mapping::EvidenceCacheReader cacheReader(config_.evidenceCacheRoot_);
            cacheReader.readBatchRecords(slide_,lane_,batchNumber_,
                                        evidenceBatchRecords_, reference_);
        }
    }


    void MergedMap2SamConverter::writeMappingRecord(const SamRecord &m) const
    {
        //CGA_ASSERT(m.mappingQuality_ == 0 || m.isPrimary_);
        if (config_.minMapQ_ > m.getMappingQuality() && !m.isPrimary_)
            return;

        if (!m.isGroupPrimary_ && !m.isPrimary_)
            return;

        Map2SamConverter::writeMappingRecord(m);
    }

} } // cgatools::mapping
