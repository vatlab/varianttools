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
#include "cgatools/command/JunctionType.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/pending/disjoint_sets.hpp>

#include "cgatools/util/Streams.hpp"
#include "cgatools/util/StringSet.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/RepeatMaskerStore.hpp"
#include "cgatools/reference/GeneDataStore.hpp"
#include "cgatools/junctions/Junction.hpp"

using std::string;
using std::vector;
using std::set;
using std::pair;
using std::map;
using boost::shared_ptr;
using boost::disjoint_sets_with_storage;
using boost::identity_property_map;
using namespace cgatools::junctions;
using namespace cgatools::reference;
namespace ba = boost::algorithm;

namespace cgatools { namespace command { namespace {

    // Max. length of target site duplication that may cause overlap
    // between arrows that point towards the same insertion point.
    const int MAX_TSD_LENGTH = 30;

    struct JunctionIdCmp
    {
        bool operator()(const string& a, const string& b) const
        {
            try
            {
                return boost::lexical_cast<int>(a) < boost::lexical_cast<int>(b);
            }
            catch (const boost::bad_lexical_cast&)
            {
                return a < b;
            }
        }
    };

    struct JunctionSideRef
    {
        int side_;
        Junctions::const_iterator junction_;

        JunctionSideRef(int side, Junctions::const_iterator junction)
            :   side_(side), junction_(junction)
        {}
    };

    typedef std::multimap<Location, JunctionSideRef> LocationJunctionMap;

    bool isReferenceArtefact(const Junction& jj)
    {
        return jj.knownUnderrepresentedRepeat_ == "Y";
    }

    void fillJunctionMap(const Junctions& src, LocationJunctionMap& dst)
    {
        for (Junctions::const_iterator jj = src.begin(); jj != src.end(); ++jj)
        {
            if (isReferenceArtefact(*jj))
                continue;

            for (int side = 0; side < 2; ++side)
            {
                Location loc = jj->sideSections_[side].position_;
                dst.insert(std::make_pair(loc, JunctionSideRef(side, jj)));
            }
        }
    }

    enum ContigType
    {
        CT_INWARD, // contig with two junction arrows pointing towards it
        CT_OUTWARD, // contig with two junction arrows pointing towards the outside
        CT_UNDEF  // not the same chromosome, or arrows point in the same direction
    };

    Range arrowPairRange(const Junction& a, int aside,
                         const Junction& b, int bside)
    {
        Location aloc = a.sideSections_[aside].position_;
        Location bloc = b.sideSections_[bside].position_;
        Range r(aloc, bloc);
        if (r.begin_ > r.end_)
            std::swap(r.begin_, r.end_);
        return r;
    }

    // Returns the range from one side of an intrachromosomal junction to the other.
    Range junctionRange(const Junction& jj)
    {
        Range r(jj.sideSections_[0].position_, jj.sideSections_[1].position_);
        CGA_ASSERT(r.begin_ <= r.end_);
        return r;
    }

    struct ContigCopyEvent
    {
        const Junction *a_, *b_;
        Range subject_; // contig that got moved
        Range target_; // contig that got replaced (possibly zero-length)
        bool inverted_;
    };

    struct MobileElementEvent
    {
        const Junction* a_;
        Range repeatRange_;
        RepeatMaskerAnnotation repeatData_;
        Range subject_; // contig that was mobilized
        Range target_;
        bool insertedAfterTarget_;
        bool inverted_;
    };

    bool isForwardArrow(const Junction& x, int side)
    {
        int strand =
                x.sideSections_[side].strand_ == JUNCTION_PLUS_STRAND ? 0 : 1;
        return (strand ^ side) == 0;
    }

    ContigType classifyContig(const Junction& a, int aside,
                              const Junction& b, int bside,
                              int maxPairingDistance, int32_t* actualDistance = NULL)
    {
        Location aloc = a.sideSections_[aside].position_;
        Location bloc = b.sideSections_[bside].position_;
        int32_t distance = aloc.distanceTo(bloc);

        if ( actualDistance != NULL ) *actualDistance = distance;

        if (distance < -maxPairingDistance || distance > maxPairingDistance)
            return CT_UNDEF;

        bool af = isForwardArrow(a, aside);
        bool bf = isForwardArrow(b, bside);
        if (af == bf)
            return CT_UNDEF;

        // Target site duplications may cause arrows that point towards
        // the insertion site to overlap. Treat all cases where arrow
        // heads are very close as "inward" contigs.
        if (std::abs(distance) <= MAX_TSD_LENGTH)
            return CT_INWARD;

        if (distance > 0)
            return af ? CT_INWARD : CT_OUTWARD;
        else
            return bf ? CT_INWARD : CT_OUTWARD;
    }

    // Returns 0 for strand-preserving junctions, +1 for +/- junctions and -1
    // for -/+ junctions
    int junctionStrandChange(const Junction& jj)
    {
        CGA_ASSERT(jj.sideSections_[0].strand_ != JUNCTION_UNKNOWN_STRAND);
        CGA_ASSERT(jj.sideSections_[1].strand_ != JUNCTION_UNKNOWN_STRAND);
        if (jj.sideSections_[0].strand_ == jj.sideSections_[1].strand_)
            return 0;
        else
            return static_cast<int>(jj.sideSections_[0].strand_);
    }

    int invert ( int side, int times )
    {
        return ( times%2 ) ? (1-side) : side;
    }

    bool checkPairing(const Junction& a, int aside,
                      const Junction& b, int bside,
                      int maxPairingDistance, int32_t maxCopyTargetLength)
    {
        if (a.id_ == b.id_)
            return false;

        const int    THIS_SIDE   = 0;
        const int    OTHER_SIDE  = 1 - THIS_SIDE;

        ContigType actualSideType[2] = { CT_UNDEF ,  CT_UNDEF };
        ContigType neededSideType[2] = { CT_OUTWARD, CT_INWARD };
        int32_t    distance[2] = { std::numeric_limits<int32_t>::max(), std::numeric_limits<int32_t>::max() };
        
        for ( int side = 0; side < 2; ++side )
        {
            actualSideType[side] = classifyContig( a, invert(aside,side),
                                                   b, invert(bside,side),
                                                   maxPairingDistance,
                                                   &distance[side] );

            if ( actualSideType[side] != neededSideType[side] )
                return false;
        }

        if ( std::abs(distance[OTHER_SIDE]) > maxCopyTargetLength ) // is it a large inversion?
        {
            // need to verify it is an inversion by comparing 
            // the distances between respective sides
            int64_t dist1 = std::abs(a.sideSections_[1-aside].position_.distanceTo
                                    (b.sideSections_[  bside].position_)),
                    dist2 = std::abs(a.sideSections_[  aside].position_.distanceTo
                                    (b.sideSections_[1-bside].position_));
            if ( dist1+dist2 > 2*MAX_TSD_LENGTH )
            {   // doesn't qualify as an inversion
                return false;
            }
        }
        
        uint16_t fromChr = a.sideSections_[aside].position_.chromosome_;
        uint16_t toChr = a.sideSections_[1-aside].position_.chromosome_;

        // Pairs of intrachromosomal events where junctions are strictly
        // ordered are not very interesting (although interpretation of
        // the pair as a contig movement is formally correct, it's not much
        // more plausible than the two independent events implied by the
        // junctions)
        if (fromChr == toChr)
        {
            Range ar = junctionRange(a);
            Range br = junctionRange(b);
            if (!ar.intersects(br))
                return false;
        }

        // Strand preserving pairings are interesting when a +/+ junction is
        // paired with a -/- junction. Strand changing pairings are interesting
        // when they cancel the strand change.
        int sca = junctionStrandChange(a);
        int scb = junctionStrandChange(b);
        if (sca == 0 && scb == 0)
            return a.sideSections_[aside].strand_ != b.sideSections_[bside].strand_;
        else
            return (0 == sca + scb);
    }

    void findPair(const Junction& start, int side,
                  const LocationJunctionMap& others,
                  int32_t maxPairingDistance,
                  int32_t maxCopyTargetLength,
                  vector<ContigCopyEvent>& result)
    {
        const JunctionSideSection& jside = start.sideSections_[side];
        Location startLoc = jside.position_;

        bool searchForward = !isForwardArrow(start, side);

        LocationJunctionMap::const_iterator p = others.lower_bound(startLoc);
        CGA_ASSERT(p != others.end()); // we must find at least "start" itself
        while (true)
        {
            const Junction& other = *p->second.junction_;
            Location otherLoc = p->first;
            int32_t distance = startLoc.distanceTo(otherLoc);

            if (distance < -maxPairingDistance || distance > maxPairingDistance)
                return;

            if (checkPairing(start, side, other, p->second.side_,
                             maxPairingDistance, maxCopyTargetLength))
            {
                ContigCopyEvent evt;
                evt.a_ = &start;
                evt.b_ = &other;
                evt.subject_ = arrowPairRange(start, side, other, p->second.side_);
                evt.target_ = arrowPairRange(start, 1-side, other, 1-p->second.side_);
                evt.inverted_ = start.sideSections_[0].strand_ != start.sideSections_[1].strand_;
                result.push_back(evt);
                return;
            }

            if (searchForward)
            {
                ++p;
                if (p == others.end())
                    break;
            }
            else
            {
                if (p == others.begin())
                    break;
                --p;
            }
        }

    }

    // Returns a range extended forward or backward from given location
    // by "distance" bases
    Range rangeFromLocation(const Location& startLoc, bool forward, int distance)
    {
        Range r(startLoc.chromosome_, startLoc.offset_, startLoc.offset_);
        if (forward)
        {
            r.end_ += distance;
        }
        else
        {
            int begin = r.begin_;
            begin -= distance;
            if (begin >= 0)
                r.begin_ = begin;
            else
                r.begin_ = 0;
        }
        return r;
    }

    Range junctionSideLocationRange(const Junction& jj, int side)
    {
        return rangeFromLocation(jj.sideSections_[side].position_, true, 0);
    }

    void findMobileElement(const Junction& start, int side,
                           const RepeatMaskerStore& repeats,
                           const util::StringSet& repNames,
                           double maxDivergence,
                           int maxDistance,
                           vector<MobileElementEvent>& result)
    {
        const JunctionSideSection& jside = start.sideSections_[side];
        Location startLoc = jside.position_;
        bool searchForward = !isForwardArrow(start, side);
        Range searchRange = rangeFromLocation(startLoc, searchForward, maxDistance);

        vector<RepeatMaskerStore::QueryResultType> r;
        repeats.intersect(searchRange, r);

        if (!searchForward)
            std::reverse(r.begin(), r.end());

        for (size_t ii = 0; ii < r.size(); ++ii)
        {
            const RepeatMaskerStore::QueryResultType& qrt = r[ii];

            const Range& repRange = qrt->first;
            const RepeatMaskerAnnotation& repData = qrt->second;

            // Skip the repeat if it's not strictly after/before the start location
            if (repRange.contains(startLoc))
                continue;

            // Skip the repeat if it's not on the expected strand; a
            // transposon on positive strand can only mobilize something
            // downstream of it
            if ( (searchForward && !repData.strand_) ||
                 (!searchForward && repData.strand_) )
                continue;

            if (!repNames.contains(repData.name_) || repData.divergence_ > maxDivergence)
                continue;

            Location insLoc = start.sideSections_[1-side].position_;

            MobileElementEvent evt;
            evt.a_ = &start;
            evt.repeatRange_ = repRange;
            evt.repeatData_ = repData;
            evt.inverted_ = (start.sideSections_[0].strand_ != start.sideSections_[1].strand_);
            evt.target_ = Range(insLoc.chromosome_, insLoc.offset_, insLoc.offset_);
            evt.insertedAfterTarget_ = isForwardArrow(start, 1-side);
            evt.subject_ = Range(startLoc.chromosome_, startLoc.offset_, startLoc.offset_);

            if (searchForward)
                evt.subject_.end_ = repRange.begin_;
            else
                evt.subject_.begin_ = repRange.end_;
            result.push_back(evt);
            break;
        }
    }

    class RelatedJunctionGraph
    {
    public:
        typedef vector<string> RelatedJunctions;

        RelatedJunctionGraph(const Junctions& src, const LocationJunctionMap& locMap, int maxDistance)
        {
            for (size_t ii = 0; ii < src.size(); ++ii)
            {
                if (isReferenceArtefact(src[ii]))
                    continue;

                RelatedJunctions r;
                findNeighbors(src[ii], locMap, maxDistance, r);
                related_[src[ii].id_] = r;
            }
        }

        const RelatedJunctions& relatedJunctions(const string& id) const
        {
            Storage::const_iterator r = related_.find(id);
            CGA_ASSERT(r != related_.end());
            return r->second;
        }

        static size_t getRelatedJunctions
        (
            const Junction& src,
            const LocationJunctionMap& locMap,
            int maxDistance,
            RelatedJunctions& result,
            size_t maxNeighbors = std::numeric_limits<int>::max()
        )
        {
            result.resize(0);
            if ( isReferenceArtefact(src) ) return 0;
            findNeighbors(src, locMap, maxDistance, result, maxNeighbors);
            return result.size();
        }

    private:
        typedef map<string, RelatedJunctions> Storage;
        Storage related_; // indexed by junction id

        static void findNeighbors
        (
            const Junction& start,
            const LocationJunctionMap& locMap,
            int maxDistance,
            RelatedJunctions& result,
            size_t maxNeighbors = std::numeric_limits<int>::max()
        )
        {
            for (int side = 0; side < 2; ++side)
            {
                const JunctionSideSection& jside = start.sideSections_[side];
                Location startLoc = jside.position_;
                LocationJunctionMap::const_iterator pos = locMap.lower_bound(startLoc);
                CGA_ASSERT(pos != locMap.end());

                LocationJunctionMap::const_iterator ii = pos;
                while (ii != locMap.end() && result.size() < maxNeighbors)
                {
                    if (ii->second.junction_->id_ != start.id_)
                    {
                        Location otherLoc = ii->first;
                        int32_t distance = startLoc.distanceTo(otherLoc);
                        CGA_ASSERT(distance >= 0);
                        if (distance > maxDistance)
                            break;
                        result.push_back(ii->second.junction_->id_);
                    }
                    ++ii;
                }

                ii = pos;
                while ( result.size() < maxNeighbors )
                {
                    if (ii->second.junction_->id_ != start.id_)
                    {
                        Location otherLoc = ii->first;
                        int32_t distance = startLoc.distanceTo(otherLoc);
                        // Even though we are moving backwards, we still need to
                        // check against the upper limit in case we cross chromosomes
                        if (distance > maxDistance || distance < -maxDistance)
                            break;
                        result.push_back(ii->second.junction_->id_);
                    }

                    if (ii == locMap.begin())
                        break;
                    --ii;
                }
            }
            std::sort(result.begin(), result.end());
            result.erase(std::unique(result.begin(), result.end()),
                         result.end());
        }
    };

    string formatRange(const CrrFile& crr, const Range& r)
    {
        return (boost::format("%s:%d-%d")
                % crr.listChromosomes()[r.chromosome_].getName()
                % r.begin_ % r.end_).str();
    }

    string formatLocation(const CrrFile& crr, const Location& loc)
    {
        return (boost::format("%s:%d")
                % crr.listChromosomes()[loc.chromosome_].getName()
                % loc.offset_).str();
    }

    bool allRelatedJunctionsExplained(const RelatedJunctionGraph::RelatedJunctions& related,
                                      const vector<ContigCopyEvent>& evts)
    {
        BOOST_FOREACH(const string& id, related)
        {
            bool found = false;
            BOOST_FOREACH(const ContigCopyEvent& ev, evts)
            {
                if (ev.a_->id_ == id || ev.b_->id_ == id)
                {
                    found = true;
                    break;
                }
            }
            if (!found)
                return false;
        }
        return true;
    }

    bool preferSimpleEventForJunction(const Junction& jj, int maxSimpleEventDistance)
    {
        Location a = jj.sideSections_[0].position_;
        Location b = jj.sideSections_[1].position_;
        int distance = a.distanceTo(b);
        return (distance < maxSimpleEventDistance && distance > -maxSimpleEventDistance);
    }

    bool pairedToSameJunction(const vector<ContigCopyEvent>& evts)
    {
        CGA_ASSERT(!evts.empty());
        const string id = evts.front().b_->id_;
        for (size_t ii = 1; ii < evts.size(); ++ii)
        {
            CGA_ASSERT(evts.front().a_ == evts[ii].a_);
            if (evts[ii].b_->id_ != id)
                return false;
        }
        return true;
    }

    const ContigCopyEvent& selectPairing(const vector<ContigCopyEvent>& evts)
    {
        const ContigCopyEvent* best = &evts[0];
        BOOST_FOREACH(const ContigCopyEvent& evt, evts)
        {
            // Prefer pairing with smaller subject. This is questionable.
            if (best->subject_.length() > evt.subject_.length())
                best = &evt;
        }
        return *best;
    }

    bool isInPlaceInversion(const ContigCopyEvent& evt)
    {
        Range src = evt.subject_;
        Range dst = evt.target_;

        if (!evt.inverted_)
            return false;

        if (!src.intersects(dst))
            return false;

        Range overlap = src.overlappingRange(dst);
        int srcLen = src.length(),
            overlapLen = overlap.length();
        int nonoverlapLen = srcLen - overlapLen;
        CGA_ASSERT(nonoverlapLen >= 0);
        return nonoverlapLen < 2 * MAX_TSD_LENGTH && nonoverlapLen < .1*srcLen;
    }

    enum EventType
    {
        ET_ARTEFACT,
        ET_COMPLEX,
        ET_DELETION,
        ET_TANDEM_DUPLICATION,
        ET_POSSIBLE_INVERSION,
        ET_INVERSION,
        ET_COPY,
        ET_ME_COPY,
        ET_XCHR
    };

    const char* EVENT_TYPE_STR[] = {
            "artifact",
            "complex",
            "deletion",
            "tandem-duplication",
            "probable-inversion",
            "inversion",
            "distal-duplication",
            "distal-duplication-by-mobile-element",
            "interchromosomal"
    };


    EventType classifyIsolatedJunction(const Junction& jj)
    {
        Location a = jj.sideSections_[0].position_;
        Location b = jj.sideSections_[1].position_;

        if (a.chromosome_ != b.chromosome_)
            return ET_XCHR;
        if (jj.sideSections_[0].strand_ != jj.sideSections_[1].strand_)
            return ET_POSSIBLE_INVERSION;
        if (jj.sideSections_[0].strand_ == JUNCTION_PLUS_STRAND)
            return ET_DELETION;
        else
            return ET_TANDEM_DUPLICATION;
    }

    EventType classifyCopyEvent(const ContigCopyEvent& evt)
    {
        if (isInPlaceInversion(evt))
            return ET_INVERSION;
        else
            return ET_COPY;
    }

    struct SvEvent
    {
        EventType type_;
        set<string> relatedJunctions_;
        Range subject_;
        Range target_;
        bool subjectInverted_;
        bool targetInverted_;
        Range mobileElementRange_;
        RepeatMaskerAnnotation mobileElementData_;
        int index_;
        set<string> affectedGenes_;
        set<string> containedGenes_;
        set<string> geneFusions_;

        // For messy (artefact and complex) junctions
        SvEvent(EventType type, const Junction& jj,
                const vector<string>& relatedJunctions = vector<string>())
            :   type_(type),
                relatedJunctions_(relatedJunctions.begin(), relatedJunctions.end()),
                subject_(0, 1, 0),
                target_(0, 1, 0),
                subjectInverted_(false),
                targetInverted_(false),
                mobileElementRange_(0, 1, 0),
                index_(-1)
        {
            relatedJunctions_.insert(jj.id_);
        }

        // For complex events that arise as a combination of other events
        SvEvent(const set<string>& relatedJunctions)
            :   type_(ET_COMPLEX),
                relatedJunctions_(relatedJunctions),
                subject_(0, 1, 0),
                target_(0, 1, 0),
                subjectInverted_(false),
                targetInverted_(false),
                mobileElementRange_(0, 1, 0),
                index_(-1)
        {
        }

        // For isolated junctions
        SvEvent(const Junction& jj)
            : type_(classifyIsolatedJunction(jj)),
              target_(0, 1, 0),
              targetInverted_(false),
              mobileElementRange_(0, 1, 0),
              index_(-1)
        {
            if (type_ != ET_XCHR)
            {
                subject_ = junctionRange(jj);
                subjectInverted_ = (type_ == ET_POSSIBLE_INVERSION);
            }
            else
            {
                subject_ = junctionSideLocationRange(jj, 0);
                target_ = junctionSideLocationRange(jj, 1);
                subjectInverted_ = (jj.sideSections_[0].strand_ != JUNCTION_PLUS_STRAND);
                targetInverted_ = (jj.sideSections_[1].strand_ != JUNCTION_PLUS_STRAND);
            }
            relatedJunctions_.insert(jj.id_);
        }

        // For "copy" events
        SvEvent(const ContigCopyEvent& evt)
            : type_(classifyCopyEvent(evt)),
              subject_(evt.subject_),
              target_(evt.target_),
              subjectInverted_(evt.inverted_),
              targetInverted_(false),
              mobileElementRange_(0, 1, 0),
              index_(-1)
        {
            relatedJunctions_.insert(evt.a_->id_);
            relatedJunctions_.insert(evt.b_->id_);
        }

        // For mobile element events
        SvEvent(const MobileElementEvent& evt)
            : type_(ET_ME_COPY),
              subject_(evt.subject_),
              target_(evt.target_),
              subjectInverted_(evt.inverted_),
              targetInverted_(false),
              mobileElementRange_(evt.repeatRange_),
              mobileElementData_(evt.repeatData_),
              index_(-1)
        {
            relatedJunctions_.insert(evt.a_->id_);
        }

        void printMobileElementColumns(std::ostream& out, const CrrFile& crr) const
        {
            if (type_ == ET_ME_COPY)
            {
                const Range& r = mobileElementRange_;
                const RepeatMaskerAnnotation& med = mobileElementData_;
                out << boost::format("%s:%s:%.1f\t%s\t%d\t%d\t%s")
                           % med.family_
                           % med.name_
                           % med.divergence_
                           % crr.listChromosomes()[r.chromosome_].getName()
                           % r.begin_
                           % r.end_
                           % (med.strand_ ? "-" : "+");
            }
            else
            {
                out << "\t\t\t\t";
            }
        }
    };

    typedef map<string, const Junction*> IdToJunctionMap;
    typedef map<string, int>             IdToIndexMap;

    void printRange(std::ostream& out, const CrrFile& crr, const Range& r, bool inverted)
    {
        if (r.begin_ <= r.end_)
        {
            out << crr.listChromosomes()[r.chromosome_].getName()
                << '\t' << r.begin_
                << '\t' << r.end_
                << '\t' << (r.end_ - r.begin_)
                << '\t' << (inverted ? "-" : "+");
        }
        else
        {
            out << "\t\t\t\t";
        }
    }

    class SvEventStore
    {
    public:

        void insertEvent(const Junction& jj, shared_ptr<SvEvent> evt)
        {
            events_[jj.id_] = evt;
        }

        // Adds the junction/event pair to the map.
        // Maintains the invariant that all junctions related to an event
        // are pointing to that same event in the map. As the new junctions
        // are processed, some events may merge into complex events.
        void mergeEvent(const Junction& jj, shared_ptr<SvEvent> evt)
        {
            set<string> related;
            BOOST_FOREACH(const string& id, evt->relatedJunctions_)
            {
                JunctionToSvEventMap::const_iterator it = events_.find(id);
                if (it != events_.end())
                {
                    const set<string>& otherRelated = it->second->relatedJunctions_;
                    related.insert(otherRelated.begin(), otherRelated.end());
                }
            }
            shared_ptr<SvEvent> final;
            if (!related.empty() && related != evt->relatedJunctions_)
            {
                related.insert(evt->relatedJunctions_.begin(), evt->relatedJunctions_.end());
                final.reset(new SvEvent(related));
            }
            else
            {
                related = evt->relatedJunctions_;
                final = evt;
            }
            BOOST_FOREACH(const string& id, related)
            {
                events_[id] = final;
            }
        }

        void annotateEvents(const CrrFile& crr,
                const IdToJunctionMap& byId,
                const GeneDataStore& geneData,
                int maxRangeLength,
                uint32_t regulatoryRegionLength)
        {
            int index = 0;
            BOOST_FOREACH(JunctionToSvEventMap::value_type& xx, events_)
            {
                if (xx.second->index_ < 0)
                {
                    ++index;
                    xx.second->index_ = index;
                    uniqueEvents_.push_back(xx.second);

                    SvEvent& evt = *xx.second;
                    collectAffectedGenes(evt, crr, byId, geneData, regulatoryRegionLength);
                    collectContainedGenes(evt, crr, geneData, maxRangeLength);
                }
            }
        }

        shared_ptr<SvEvent> getEvent(const Junction& jj)
        {
            return events_[jj.id_];
        }

        void writeReport(std::ostream& out, const CrrFile& crr, const IdToJunctionMap& byId)
        {
            out << ">EventId\tType\tRelatedJunctionIds\tMatePairCounts\t"
                   "FrequenciesInBaselineGenomeSet\t"
                   "OriginRegionChr\tOriginRegionBegin\tOriginRegionEnd\t"
                   "OriginRegionLength\tOriginRegionStrand\t"
                   "DestinationRegionChr\tDestinationRegionBegin\t"
                   "DestinationRegionEnd\tDestinationRegionLength\t"
                   "DestinationRegionStrand\t"
                   "DisruptedGenes\tContainedGenes\tGeneFusions\t"
                   "RelatedMobileElement\tMobileElementChr\t"
                   "MobileElementBegin\tMobileElementEnd\tMobileElementStrand\n";

            for (size_t ii = 0; ii < uniqueEvents_.size(); ++ii)
            {
                const SvEvent& evt = *(uniqueEvents_[ii]);
                set<string, JunctionIdCmp> related(evt.relatedJunctions_.begin(),
                                                   evt.relatedJunctions_.end());

                out << evt.index_ << '\t'                   //EventId
                    << EVENT_TYPE_STR[evt.type_] << '\t'    //Type
                    << ba::join(related, ";") << '\t';      //RelatedJunctionIds

                std::vector<double> freqs;
                BOOST_FOREACH(const string& id, related)
                {
                    IdToJunctionMap::const_iterator jit = byId.find(id);
                    CGA_ASSERT(jit != byId.end());
                    const Junction* jj = jit->second;
                    if (!freqs.empty())
                        out << ';';
                    out << jj->score_;
                    freqs.push_back(jj->frequencyInBaselineGenomeSet_);
                }
                out << '\t';                                //MatePairCounts
                for (size_t ii = 0; ii < freqs.size(); ++ii)
                {
                    if (0 != ii)
                        out << ';';
                    out << boost::format("%.2f") % freqs[ii];
                }
                out << '\t';                //FrequenciesInBaselineGenomeSet
                printRange(out, crr, evt.subject_, evt.subjectInverted_);
                out << '\t';                //OriginRegionChr,OriginRegionBegin,OriginRegionEnd
                printRange(out, crr, evt.target_, evt.targetInverted_);
                out << '\t';                //DestinationRegionChr,DestinationRegionBegin,DestinationRegionEnd
                out << ba::join(evt.affectedGenes_, ";")
                    << '\t'                                 
                    << ba::join(evt.containedGenes_, ";")
                    << '\t'
                    << ba::join(evt.geneFusions_, ";")
                    << '\t';
                evt.printMobileElementColumns(out, crr);
                out << '\n';
            }
        }

    private:
        typedef map< string, shared_ptr<SvEvent>, JunctionIdCmp > JunctionToSvEventMap;
        JunctionToSvEventMap events_;
        vector< shared_ptr<SvEvent> > uniqueEvents_;

        void collectAffectedGenes(SvEvent& evt,
                const CrrFile& crr,
                const IdToJunctionMap& byId,
                const GeneDataStore& geneData,
                uint32_t regulatoryRegionLength)
        {
            BOOST_FOREACH(const string& jid, evt.relatedJunctions_)
            {
                IdToJunctionMap::const_iterator jit = byId.find(jid);
                CGA_ASSERT(jit != byId.end());
                const Junction* jj = jit->second;

                findGenesOverJunctionSide(*jj, 0, geneData, evt.affectedGenes_);
                findGenesOverJunctionSide(*jj, 1, geneData, evt.affectedGenes_);
                findGeneFusions(evt, *jj, crr, geneData, regulatoryRegionLength);
            }
        }

        void findGenesOverJunctionSide(const Junction& jj, int side,
                                       const GeneDataStore& geneData,
                                       set<string>& result)
        {
            Range r = junctionSideLocationRange(jj, side);
            vector<GeneDataStore::QueryResultType> overlap;
            geneData.intersect(r, overlap);
            BOOST_FOREACH(const GeneDataStore::QueryResultType& ii, overlap)
            {
                if (!ii->second.symbol_.empty())
                    result.insert(ii->second.symbol_);
            }
        }

        // For events that duplicate or remove sequence, find the
        // genes fully contained in that sequence
        void collectContainedGenes(SvEvent& evt,
                const CrrFile& crr,
                const GeneDataStore& geneData,
                int maxRangeLength)
        {
            if (evt.type_ != ET_DELETION &&
                    evt.type_ != ET_TANDEM_DUPLICATION &&
                    evt.type_ != ET_COPY)
                return;
            if (maxRangeLength >= 0 &&
                    evt.subject_.length() >= static_cast<uint32_t>(maxRangeLength))
                return;

            Range r = evt.subject_;
            vector<GeneDataStore::QueryResultType> overlap;
            geneData.intersect(r, overlap);

            BOOST_FOREACH(const GeneDataStore::QueryResultType& ii, overlap)
            {
                if (ii->second.symbol_.empty())
                    continue;
                // Skip partially overlapping genes
                Range common = r.overlappingRange(ii->first);
                if (common != ii->first)
                    continue;
                evt.containedGenes_.insert(ii->second.symbol_);
            }
        }

        void findGeneFusions(SvEvent& evt, const Junction& jj,
                const CrrFile& crr,
                const GeneDataStore& geneData,
                uint32_t regulatoryRegionLength)
        {
            boost::array<vector<GeneDataStore::QueryResultType>, 2> genes;
            for (size_t side = 0; side < 2; ++side)
            {
                Range r = geneSearchRange(jj, side, regulatoryRegionLength);
                geneData.intersect(r, genes[side]);
            }
            BOOST_FOREACH(const GeneDataStore::QueryResultType& left, genes[0])
            {
                if (left->second.symbol_.empty())
                    continue;
                bool leftStrand = (jj.sideSections_[0].strand_ != JUNCTION_PLUS_STRAND);
                bool leftWithGene = (leftStrand == left->second.strand_);
                BOOST_FOREACH(const GeneDataStore::QueryResultType& right, genes[1])
                {
                    if (right->second.symbol_.empty())
                        continue;
                    // Skip genes that already overlap on reference, this
                    // also handles junctions within a single gene
                    if (left->first.intersects(right->first))
                        continue;

                    bool rightStrand = (jj.sideSections_[1].strand_ != JUNCTION_PLUS_STRAND);
                    bool rightWithGene = (rightStrand == right->second.strand_);

                    if (leftWithGene == rightWithGene)
                    {
                        RangeRelationship
                            upRel = rangeRelationship(left->first,
                                                      jj.sideSections_[0].position_,
                                                      left->second.strand_),
                            downRel = rangeRelationship(right->first,
                                                        jj.sideSections_[1].position_,
                                                        right->second.strand_);
                        const string* upSym = &left->second.symbol_;
                        const string* downSym = &right->second.symbol_;
                        Range upRange = extendedRange(left->first,
                                                      regulatoryRegionLength,
                                                      left->second.strand_);
                        Range downRange = extendedRange(right->first,
                                                        regulatoryRegionLength,
                                                        right->second.strand_);

                        if (!leftWithGene)
                        {
                            std::swap(upRel, downRel);
                            std::swap(upSym, downSym);
                            std::swap(upRange, downRange);
                        }

                        string fusion;
                        if (upRel == INSIDE && downRel == INSIDE)
                        {
                            fusion = (*upSym) + "/" + (*downSym);
                        }
                        else if (upRel == UPSTREAM && downRel == UPSTREAM)
                        {
                            // Skip regulatory sequence "fusions" that are present
                            // on unmodified reference. This also skips pairs of genes
                            // on opposite strands, that already share regulatory
                            // sequence; it's not clear whether we want such pairs
                            // to be reported as fusions even if the junction inverts
                            // the strand
                            if (downRange.intersects(upRange))
                                continue;
                            fusion = "TSS-UPSTREAM[" + (*upSym) + "]/" + (*downSym);
                        }
                        else
                        {
                            continue;
                        }
                        evt.geneFusions_.insert(fusion);
                    }
                }
            }
        }

        Range geneSearchRange(const Junction& jj, int side, uint32_t len)
        {
            Range r = junctionSideLocationRange(jj, side);
            if (r.begin_ >= len)
                r.begin_ -= len;
            else
                r.begin_ = 0;
            r.end_ += len;
            return r;
        }

        Range extendedRange(Range r, uint32_t len, bool strand)
        {
            if (strand)
            {
                r.end_ += len;
            }
            else
            {
                if (r.begin_ > len)
                    r.begin_ -= len;
                else
                    r.begin_ = 0;
            }
            return r;
        }

        enum RangeRelationship
        {
            UPSTREAM, INSIDE, DOWNSTREAM
        };

        RangeRelationship rangeRelationship(const Range& range, const Location& p, bool strand)
        {
            CGA_ASSERT(range.chromosome_ == p.chromosome_);
            if (range.contains(p))
                return INSIDE;
            else if ((strand && p.offset_ >= range.end_) ||
                         (!strand && p.offset_ < range.begin_))
                return UPSTREAM;
            else
                return DOWNSTREAM;
        }
    };

    void printJunctionAnnotation
    (   std::ostream&   out, 
        const CrrFile&  crr,
        const Junction& jj, 
        const SvEvent&  evt,
        int             maxCount = std::numeric_limits<int>::max() 
    )
    {
        out << evt.index_ << '\t'
            << EVENT_TYPE_STR[evt.type_] << '\t';

        // array of int-string pairs: int = actial number, string = junction id
        vector<pair<int,string> > ids;
        ids.reserve( std::min<int>(maxCount,(int)evt.relatedJunctions_.size()) );

        // fill the array of pairs so that the lexical cast is only done once per junction
        BOOST_FOREACH(const string& id, evt.relatedJunctions_)
        {
            if ( id != jj.id_ )
            {
                ids.push_back(pair<int,string>( boost::lexical_cast<int>(id),id));
            }
        }

        // if we're only printing a subset of maxCount size, then resize before sorting
        if ( maxCount < (int)ids.size() )
        {
            std::nth_element(ids.begin(),ids.begin()+maxCount, ids.end());
            ids.resize(maxCount);
        }
        std::sort(ids.begin(),ids.end());

        // write semicolon-separated junction ids
        vector<pair<int,string> >::const_iterator it  = ids.begin(),
                                                  end = ids.end();
        if ( it != end )
            for(out << (it++)->second; it != end; ++it )
                out << ';' << it->second;

        out << '\n';
    }

}}} // local namespace

namespace cgatools { namespace command {
    JunctionType::JunctionType(const string& name)
        :
    Command(
      name,
      "Groups and annotates junction calls by event type.",
      "1.5 or later",

      "This tool searches for groups of related junctions and for every "
      "group attempts to determine the event that caused the junctions. "
      "For example, isolated strand-consistent intrachromosomal junction "
      "is likely to be caused by a deletion event.\n"
      "Every junction in the file specified by \"junctions\" parameter "
      "will be annotated. Optionally, the tool can search for the related "
      "junctions in a larger list of junctions specified by \"all-junctions\" "
      "parameter. For example, one may use the high confidence junction file "
      "to restrict the list of events to ones that contain at least one "
      "high-confidence junction, while using the complete list of all "
      "junctions to make sure that even low-confidence junctions will be "
      "taken into account when grouping the junctions and determining the "
      "event type.\n"
      "The output consists of two files, [prefix]AnnotatedJunctions.tsv and "
      "[prefix]Events.tsv. The annotated junction file contains the junctions "
      "from the primary input file annotated by the following columns:\n\n"
      "    EventId           \tInteger id that links the junction file "
      "to the event file\n"
      "    Type              \tType of the event that caused the junction\n"
      "    RelatedJunctions  \tSemicolon-separated list of other junctions "
      "that were grouped with this junction\n\n"
      "The event list file contains the following columns:\n\n"
      "    EventId     \tUnique id of the event\n"
      "    Type        \tType of the event. One of the following values:\n"
      "        artifact  \tcaused by a flaw in the reference\n"
      "        complex   \tevent involves multiple junctions and doesn't "
      "fit the pattern of any simple event type\n"
      "        deletion  \tdeletion of the sequence described by the Origin "
      "columns\n"
      "        tandem-duplication \ttandem duplication of the origin sequence\n"
      "        probable-inversion \tinversion of the origin sequence that "
      "is confirmed from one side of the inversion only\n"
      "        inversion \tinversion of the origin sequence replacing the "
      "sequence described by the Destination columns, confirmed from both "
      "sides\n"
      "        distal-duplication \tcopy of the origin sequence into the area "
      "described by the Destination columns\n"
      "        distal-duplication-by-mobile-element \tcopy of the origin sequence "
      "caused by a known active mobile element\n"
      "        interchromosomal \tisolated junction between different "
      "chromosomes; Origin and Destination columns describe the reference "
      "loci that are brought together by this event.\n"
      "    RelatedJunctionIds \tSemicolon-separated list of the junctions "
      "related to this event.\n"
      "    MatePairCounts \tSemicolon-separated list that contains the read "
      "count for every related junction.\n"
      "    FrequenciesInBaselineGenomeSet \tSemicolon-separated list that "
      "contains the frequency in the baseline set of genomes for every "
      "related junction.\n"
      "    OriginRegion[...] \tDescription of the origin "
      "sequence of the event; the exact semantics of \"origin\" depend on the "
      "event type.\n"
      "    DestinationRegion[...] \tDescription of the "
      "destination region for the event.\n"
      "    DisruptedGenes \tList of all genes that contain one or more of the "
      "locations of the junctions grouped to this event.\n"
      "    ContainedGenes \tFor the events that duplicate or remove regions "
      "of sequence, this column contains the list of genes fully contained "
      "within the deleted or copied region.\n"
      "    GeneFusions \tList of possible gene fusions described as GeneA/GeneB, "
      "and fusions of regulatory sequence of one gene to another gene, "
      "described as TSS-UPSTREAM[GeneA]/GeneB.\n"
      "    RelatedMobileElement \tFor the duplication events caused by a mobile "
      "element, this column contain the description of the element in "
      "Family:Name:DivergencePercent format, for example \"L1:L1HS:0.5\".\n"
      "    MobileElement[...] \tLocation of the mobile "
      "element\n\n"
      "All sequence intervals are described using zero-based, half-open "
      "coordinates.\n\n"
      "Repeat and gene data files necessary to run this command can be "
      "downloaded from the Complete Genomics site:\n\n"
      "ftp://ftp.completegenomics.com/AnnotationFiles/"
    )
    {
        options_.add_options()
            ("reference", po::value<string>(&referenceFileName_),
              "Reference file.")
            ("output-prefix", po::value<string>(&outputPrefix_),
              "The path prefix for all output reports.")
            ("junctions", po::value<string>(&junctionFileName_),
              "Primary input junction file.")
            ("all-junctions", po::value<string>(&fullJunctionFileName_),
              "Superset of the input junction file to use when searching for the "
              "related junctions. The default is to use only the junctions in the "
              "primary junction file.")
            ("repmask-data", po::value<string>(&repmaskFileName_),
              "The file that contains repeat masker data.")
            ("gene-data", po::value<string>(&geneFileName_),
              "The file that contains gene location data.")
            ("regulatory-region-length",
              po::value<uint32_t>(&regulatoryRegionLength_)->default_value(7500),
              "Length of the region upstream of the gene that may contain regulatory "
              "sequence for the gene. Junctions that connect this region to another gene "
              "will be annotated as a special kind of gene fusion.")
            ("contained-genes-max-range",
              po::value<int>(&maxContainedGenesRangeLength_)->default_value(-1),
              "Maximum length of a copy or deletion event to annotate with all genes "
              "that overlap the copied or deleted segment. Negative value causes all "
              "events to be annotated regardless of the length.")
            ("max-related-junction-distance",
              po::value<int>(&maxRelatedJunctionDistance_)->default_value(700),
              "Junctions occurring within this distance are presumed to be related.")
            ("max-pairing-distance",
              po::value<int>(&maxPairingDistance_)->default_value(10*1000*1000),
              "When searching for paired junctions caused by the same event, maximum "
              "allowed distance between junction sides.")
            ("max-copy-target-length",
              po::value<int>(&maxCopyTargetLength_)->default_value(1000),
              "Pairs of junctions will be classified as a copy event only if the "
              "length of the implied copy target region is below this threshold.")
            ("max-simple-event-distance",
              po::value<int>(&maxSimpleEventDistance_)->default_value(10*1000*1000),
              "When given a choice of explaining an event as a mobile element copy or "
              "as a simple deletion/duplication, prefer the latter explanation if the "
              "length of the affected sequence if below this threshold.")
            ("mobile-element-names",
              po::value<string>(&mobileElementNames_)->default_value("L1HS,SVA"),
              "Comma-separated list of the names of the mobile elements that are "
              "known to be active and sometimes copy flanking 3' sequence.")
            ("max-distance-to-m-e",
              po::value<int>(&maxDistanceToTransposon_)->default_value(2000),
              "When searching for a mobile element related to a junction, maximum allowed "
              "distance from the junction side to the element.")
            ("max-related-junction-output",
               po::value<int>(&maxRelatedJunctionOutput_)->default_value(100),
               "Maximum number of related junctions included into annotation field")
            ;
    }

    int JunctionType::run(po::variables_map& vm)
    {
        requireParam(vm, "reference");
        requireParam(vm, "junctions");
        requireParam(vm, "repmask-data");
        requireParam(vm, "gene-data");

        if (fullJunctionFileName_.empty())
            fullJunctionFileName_ = junctionFileName_;

        reference::CrrFile crr;
        crr.open(referenceFileName_);
        
        JunctionFile junctions   (junctionFileName_, crr);
        JunctionFile allJunctions(fullJunctionFileName_, crr);

        // Build ID-to-junction and ID-to-index map for all junctions
        IdToJunctionMap allJunctionsById;
        IdToIndexMap    allIndicesById;

        BOOST_FOREACH(const Junction& jj, allJunctions.junctions_)
        {
            if ( allJunctionsById.count(jj.id_) != 0 )
                throw util::Exception("all-junctions file has duplicate junction id " + jj.id_);

            allJunctionsById[jj.id_] = &jj;
            int next = allIndicesById.size();
            allIndicesById  [jj.id_] = next;
        }

        // Make sure junctions are a subset of allJunctions
        BOOST_FOREACH(const Junction& jj, junctions.junctions_)
        {
            if (allJunctionsById.find(jj.id_) == allJunctionsById.end())
                throw util::Exception("all-junctions file is not a superset of junctions "
                                      "file; id not found: " + jj.id_);
        }

        LocationJunctionMap junctionMap;
        fillJunctionMap(allJunctions.junctions_, junctionMap);

        RepeatMaskerStore repeats(crr, repmaskFileName_);
        util::StringSet mobilizingRepeatNames(mobileElementNames_, "", "");

        GeneDataStore geneDataStore(crr, geneFileName_);
        SvEventStore events;

        // Initialize disjoint set of indices for all junctions
        disjoint_sets_with_storage<identity_property_map,identity_property_map>
            djSet(allJunctions.junctions_.size());

        for ( int i = 0, n = allJunctions.junctions_.size(); i < n; ++i )
            djSet.make_set(i);

        // Process the main list of junctions

        BOOST_FOREACH(const Junction& jj, junctions.junctions_)
        {
            int i = allIndicesById[jj.id_];

            shared_ptr<SvEvent> current;

            if (isReferenceArtefact(jj))
            {
                current.reset(new SvEvent(ET_ARTEFACT, jj));
            }
            else
            {
                // find junctions that are related to jj
                RelatedJunctionGraph::RelatedJunctions related;
                RelatedJunctionGraph::getRelatedJunctions
                ( 
                    jj, junctionMap,
                    maxRelatedJunctionDistance_,related
                );

                // join tojether the sets of related junctions
                BOOST_FOREACH(const string& idrj,related)
                {
                    djSet.union_set(i,allIndicesById[idrj]);
                }

                // complex ??
                bool isComplex = (related.size() > 1);

                if ( isComplex ) // handle as comples
                {
                    current.reset(new SvEvent(ET_COMPLEX, jj));
                }
                else // not complex? find possible pairings
                {
                    vector<ContigCopyEvent> evts;
                    findPair(jj, 0, junctionMap,
                             maxPairingDistance_, maxCopyTargetLength_,
                             evts);
                    findPair(jj, 1, junctionMap,
                             maxPairingDistance_, maxCopyTargetLength_,
                             evts);

                    vector<MobileElementEvent> reps;
                    findMobileElement(jj, 0, repeats, mobilizingRepeatNames,
                            50.0, maxDistanceToTransposon_, reps);
                    findMobileElement(jj, 1, repeats, mobilizingRepeatNames,
                            50.0, maxDistanceToTransposon_, reps);

                    // still could be a complex one
                    isComplex = !allRelatedJunctionsExplained(related, evts);

                    if (isComplex)
                    {
                        current.reset(new SvEvent(ET_COMPLEX, jj));
                    }
                    else if (related.empty() && evts.empty())
                    {
                        if (reps.empty() || preferSimpleEventForJunction(jj, maxSimpleEventDistance_))
                        {
                            current.reset(new SvEvent(jj));
                        }
                        else
                        {
                            const MobileElementEvent& evt = reps[0];
                            current.reset(new SvEvent(evt));
                        }
                    }
                    else
                    {
                        CGA_ASSERT(!evts.empty());
                        const ContigCopyEvent& evt = selectPairing(evts);
                        current.reset(new SvEvent(evt));
                    }
                }
            }

            // just add the event to the list
            events.insertEvent(jj, current);
        }

        // Map type: cluster number to set of related junctions' ids
        typedef std::map<int, std::set<std::string> > IntToStringSetMap;
        IntToStringSetMap junctionClusters;

        // Map type: cluster number to event
        typedef std::map<int, shared_ptr<SvEvent>   > IntToSvEventMap;
        IntToSvEventMap   clusterMap;

        // build the clusters' sets of related junctions' ids
        BOOST_FOREACH(const Junction& jj, allJunctions.junctions_)
        {
            int i = allIndicesById[jj.id_];
            int c = djSet.find_set(i);

            junctionClusters[ c ].insert(jj.id_);
        }

        // build the initial list of clusters' events
        BOOST_FOREACH(const Junction& jj, junctions.junctions_)
        {
            int i = allIndicesById[jj.id_];
            int c = djSet.find_set(i);

            if ( clusterMap.count(c) == 0 && (events.getEvent(jj) != NULL) )
                clusterMap[c] = events.getEvent(jj);
        }

        // build the final list of events
 
        BOOST_FOREACH(const Junction& jj, junctions.junctions_)
        {
            int i = allIndicesById[jj.id_];
            int c = djSet.find_set(i);

            if  ( (clusterMap.count(c) == 0) ||
                  (clusterMap[c] == NULL)    || 
                  (clusterMap[c]->relatedJunctions_ != junctionClusters[c]) )
            {
                clusterMap[c].reset(new SvEvent(junctionClusters[c]) );
            }

            events.insertEvent(jj, clusterMap[c]);
        }

        // annotate events

        events.annotateEvents(
                crr, allJunctionsById, geneDataStore,
                maxContainedGenesRangeLength_,
                regulatoryRegionLength_);

        // write tjh annotated junctions file
        shared_ptr<std::ostream> out =
                util::OutputStream::openCompressedOutputStreamByExtension(
                    outputPrefix_ + "AnnotatedJunctions.tsv");
        util::DelimitedFileMetadata metadata(junctions.metadata_);
        metadata.initDefaults();
        *out << metadata;
        *out << JunctionFile::header_
             << "\tEventId\tType\tRelatedJunctions\n";

        BOOST_FOREACH(const Junction& jj, junctions.junctions_)
        {
            jj.write(*out, crr, 0);
            *out << "\t";

            shared_ptr<SvEvent> evt = events.getEvent(jj);
            printJunctionAnnotation(*out, crr, jj, *evt, maxRelatedJunctionOutput_);
        }

        // write the events file
        out = util::OutputStream::openCompressedOutputStreamByExtension(
                    outputPrefix_ + "Events.tsv");
        util::DelimitedFileMetadata evtMeta;
        evtMeta.initDefaults();
        evtMeta.add("TYPE", "SV-EVENTS");
        if ( !evtMeta.hasKey("SOFTWARE_VERSION") )
        {
            evtMeta.add("SOFTWARE_VERSION", "0.0.0.0");
        }
        *out << evtMeta;
        events.writeReport(*out, crr, allJunctionsById);

        return 0;
    }


} } // cgatools::command
