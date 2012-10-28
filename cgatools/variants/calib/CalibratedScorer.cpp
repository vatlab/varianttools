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
#include "cgatools/variants/calib/CalibratedScorer.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include "cgatools/util/Streams.hpp"

#include <cmath>
#include <algorithm>
#include <map>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace std;
using boost::shared_ptr;
using namespace cgatools::util;
namespace bfs = boost::filesystem;

namespace cgatools { namespace variants { namespace calib {

    CoverageBinner::CoverageBinner(const std::vector<int32_t>& cvgLevels)
        : cvgLevels_(cvgLevels)
    {
        std::sort(cvgLevels_.begin(), cvgLevels_.end());
        CGA_ASSERT(cvgLevels_.size() > 0);
        CGA_ASSERT(0 == cvgLevels_[0]);
    }

    CalibratedScorer::CalibratedScorer(
        const std::string& varType,
        const std::string& scoreType,
        bool eaf,
        const std::string& dataPath,
        const std::string& softwareVersion,
        double a20Mixture,
        double refBasesPerHetVariant)
        : minScore_(std::numeric_limits<int32_t>::max())
    {
        try
        {
            string eafStr = eaf ? "eaf" : "vaf";
            string calibId = ( boost::format("%s-%s-%s")
                               % varType % eafStr % scoreType ).str();

            shared_ptr<istream> inMetrics, inData;
            string fnMetrics, fnData, fnDataA20;
            getScoreStreams(calibId, dataPath, softwareVersion,
                            inMetrics, inData, fnMetrics, fnData, fnDataA20);

            // Gather refBasesPerHetVariantCalibration.
            double refBasesPerHetVariantCalibration;
            {
                DelimitedFile df(*inMetrics, fnMetrics);

                string metric,value;
                df.addField(StringField("metric", &metric));
                df.addField(StringField("value", &value));
                map<string,string> mm;
                while (df.next())
                {
                    mm[metric] = value;
                }
                if (mm.find("calib-format-version") == mm.end() ||
                    mm["calib-format-version"] != "1")
                    throw Exception("unexpected calib-format-version: "+fnMetrics);
                CGA_ASSERT(mm.find("ref-bases") != mm.end());
                CGA_ASSERT(mm.find("count-het-"+varType) != mm.end());

                double refBases = boost::lexical_cast<double>(mm["ref-bases"]);
                double hetCount = boost::lexical_cast<double>(mm["count-het-"+varType]);
                CGA_ASSERT(refBases > 0.0);
                CGA_ASSERT(hetCount > 0.0);
                refBasesPerHetVariantCalibration = refBases / hetCount;
            }

            // Gather scores.
            readData(*inData, fnData, binner_, S_, minScore_);

            // Correct for rate of het variant in this sample. The fn rate
            // should increase as the count of variants, not the count of
            // bases. The fp rate should increase as the count of bases, not
            // variants.
            if (0.0 != refBasesPerHetVariant && "fn" == scoreType)
            {
                double scoreCorrection =
                    10.0 * std::log10(refBasesPerHetVariant / refBasesPerHetVariantCalibration);
                for(size_t ii=0; ii<binner_.getBinCount(); ii++)
                {
                    for(size_t jj=0; jj<S_[ii].size(); jj++)
                    {
                        S_[ii][jj] += scoreCorrection;
                    }
                }
            }
            else if (0.0 != refBasesPerHetVariant && "fp" == scoreType)
            {
                double scoreCorrection =
                    -10.0 * std::log10(refBasesPerHetVariant / refBasesPerHetVariantCalibration);
                for(size_t ii=0; ii<binner_.getBinCount(); ii++)
                {
                    for(size_t jj=0; jj<S_[ii].size(); jj++)
                    {
                        S_[ii][jj] += scoreCorrection;
                    }
                }
            }

            // Correct for mixture.
            if (a20Mixture > 0.0)
            {
                CGA_ASSERT(a20Mixture <= 1.0);
                InputStream inDataA20(fnDataA20);
                CoverageBinner a20Binner;
                vector< vector<double> > a20S;
                int32_t a20MinScore;
                readData(inDataA20, fnDataA20, a20Binner, a20S, a20MinScore);

                CGA_ASSERT(a20Binner.getBinCount() == binner_.getBinCount());
                CGA_ASSERT(a20S.size() == S_.size());
                CGA_ASSERT(a20MinScore == minScore_);

                for(size_t ii=0; ii<S_.size(); ii++)
                {
                    CGA_ASSERT(a20S[ii].size() == S_[ii].size());
                    for(size_t jj=0; jj<S_[ii].size(); jj++)
                    {
                        double a50L = std::pow(10.0, -0.1*S_[ii][jj]);
                        double a20L = std::pow(10.0, -0.1*a20S[ii][jj]);
                        double mixL = 1.0 / ( a20Mixture / a20L + (1-a20Mixture) / a50L );
//                         cout << ii << " " << jj
//                              << " " << (int(ii)+int(minScore_))
//                              << " " << S_[ii][jj]
//                              << " " << a20S[ii][jj]
//                              << " " << (-10.0 * std::log10(mixL))
//                              << endl;
                        S_[ii][jj] = -10.0 * std::log10(mixL);
                    }
                }
            }

            // Init L_, PTrue_.
            L_.resize(S_.size());
            PTrue_.resize(S_.size());
            for(size_t ii=0; ii<S_.size(); ii++)
            {
                for(size_t jj=0; jj<S_[ii].size(); jj++)
                {
                    L_[ii].push_back(std::pow(10.0, -0.1*S_[ii][jj]));
                    PTrue_[ii].push_back(1.0 / (1.0 + L_[ii][jj]));
                    CGA_ASSERT(0.0 <= PTrue_[ii][jj] && PTrue_[ii][jj] <= 1.0);
                }
            }
        }
        catch(std::exception& ee)
        {
            throw Exception("failed to load calibration data "+dataPath+": "+ee.what());
        }
    }

    void CalibratedScorer::getScoreStreams(
        const std::string& calibId,
        const std::string& dataPath,
        const std::string& softwareVersion,
        boost::shared_ptr<std::istream>& inMetrics,
        boost::shared_ptr<std::istream>& inData,
        std::string& fnMetrics,
        std::string& fnData,
        std::string& fnDataA20) const
    {
        int swVersion = swVersionToInt(softwareVersion);
        if (0 == swVersion)
            swVersion = 999999;
        int bestVersion = 0;
        string bestPath = dataPath + "/version0.0.0";
        bfs::directory_iterator last;
        for(bfs::directory_iterator first(dataPath); first!=last; ++first)
        {
            if (InputStream::isReadable( (first->path() / "metrics.tsv").string() ) &&
                boost::starts_with(first->path().leaf(), "version"))
            {
                int iiVersion = swVersionToInt(first->path().leaf().substr(7));
                if (iiVersion <= swVersion && iiVersion > bestVersion)
                {
                    bestVersion = iiVersion;
                    bestPath = first->path().string();
                }
            }
        }

        fnMetrics = bestPath + "/metrics.tsv";
        fnData    = bestPath + "/" + calibId + ".tsv";
        fnDataA20 = bestPath + "/" + calibId + "-af20.tsv";

        inMetrics.reset(new InputStream(fnMetrics));
        inData   .reset(new InputStream(fnData));
    }

    void CalibratedScorer::readData(
        std::istream& inData,
        const std::string& fnData,
        CoverageBinner& binner,
        std::vector< std::vector<double> >& SS,
        int32_t& minScore)
    {
        DelimitedFile df(inData, fnData);

        // Determine cvgLevels.
        vector<int32_t> cvgLevels;
        BOOST_FOREACH(const string& colHeader, df.getColumnHeaders())
        {
            if (boost::starts_with(colHeader, "cvg"))
                cvgLevels.push_back(boost::lexical_cast<int32_t>(colHeader.substr(3)));
        }
        binner = CoverageBinner(cvgLevels);
        SS.resize(binner_.getBinCount());

        int32_t score = 0;
        vector<double> oneScoreData(binner_.getBinCount(), 0.0);
        for(size_t ii=0; ii<binner_.getBinCount(); ii++)
            df.addField(ValueField<double>(
                            "cvg" + boost::lexical_cast<string>(binner_.getMinCvg(ii)),
                            &oneScoreData[ii] ));
        df.addField(ValueField<int32_t>("score", &score));

        int32_t prevScore(std::numeric_limits<int32_t>::max());
        while (df.next())
        {
            if (prevScore != std::numeric_limits<int32_t>::max())
            {
                CGA_ASSERT(score == prevScore+1);
            }
            else
            {
                minScore = score;
            }
            prevScore = score;

            for(size_t ii=0; ii<binner_.getBinCount(); ii++)
                SS[ii].push_back(oneScoreData[ii]);
        }
    }

    int CalibratedScorer::swVersionToInt(const std::string& versionStr) const
    {
        if (versionStr.find("$") != string::npos)
            return 999999;
        vector<string> parts;
        boost::split(parts, versionStr, boost::is_any_of("."));
        if (parts.size() != 4 && parts.size() != 3)
            throw Exception("invalid software version: "+versionStr);
        int result = 0;
        for(size_t ii=0; ii<3; ii++)
        {
            int partInt = boost::lexical_cast<int>(parts[ii]);
            if (ii != 0 && partInt > 100)
                throw Exception("invalid software version: "+versionStr);
            result *= 100;
            result += partInt;
        }
        return result;
    }

} } } // cgatools::variants::calib
