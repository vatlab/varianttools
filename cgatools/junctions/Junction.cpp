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


// This file contains implementations of the classes defined in
// Junction.hpp.


// Complete Genomics.
#include "cgatools/core.hpp"
#include "cgatools/junctions/Junction.hpp"
#include "cgatools/util/parse.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/BaseUtil.hpp"
#include "cgatools/util/DelimitedFile.hpp"

// Boost.
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

//STL
#include <map>

namespace cgatools { namespace junctions {

    const std::string JunctionFile::header_ =
        ">Id\t"
        "LeftChr\tLeftPosition\tLeftStrand\tLeftLength\t"
        "RightChr\tRightPosition\tRightStrand\tRightLength\t"
        "StrandConsistent\tInterchromosomal\tDistance\t"
        "DiscordantMatePairAlignments\tJunctionSequenceResolved\t"
        "TransitionSequence\tTransitionLength\t"
        "LeftRepeatClassification\tRightRepeatClassification\t"
        "LeftGenes\tRightGenes\t"
        "XRef\tDeletedTransposableElement\t"
        "KnownUnderrepresentedRepeat\t"
        "FrequencyInBaselineGenomeSet\t"
        "AssembledSequence"
        ;

    const std::string JunctionFile::SEP = "\t";


    JunctionSideSection::JunctionSideSection()
        :
        strand_(JUNCTION_UNKNOWN_STRAND),
        length_(-1)
    {
    }


    JunctionSideSection::JunctionSideSection(
        JunctionStrand strand,
        const reference::Location& position,
        int length, 
        const std::string& repeatClassification,
        const std::string& genes)
        :
        strand_(strand),
        position_(position),
        length_(length),
        repeatClassification_(repeatClassification),
        genes_(genes)
    {
    }

    void JunctionSideSection::write(std::ostream& out, 
        const reference::CrrFile& reference) const
    {
        const std::string chromosomeName = reference.listChromosomes()[position_.chromosome_].getName();

        out << chromosomeName << JunctionFile::SEP;
        out << position_.offset_ << JunctionFile::SEP;

        CGA_ASSERT(strand_==JUNCTION_MINUS_STRAND || strand_==JUNCTION_PLUS_STRAND);
        out << (int(strand_)==JUNCTION_MINUS_STRAND ? "-":"+") << JunctionFile::SEP;

        out << length_;
    }

    int JunctionSideSection::getDir( size_t side ) const
    {
        if ((side == 0 && strand_ == JUNCTION_PLUS_STRAND) ||
            (side == 1 && strand_ == JUNCTION_MINUS_STRAND))
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    reference::Location JunctionSideSection::getBasePos( size_t side ) const
    {
        int dir = getDir( side );
        reference::Location pos = position_;
        if (dir == 1)
        {
            CGA_ASSERT(pos.offset_ > 0);
            --pos.offset_;
        }
        return pos;
    }


    Junction::Junction()
    :
        transitionIsKnown_(false),
        score_(0)
    {
    }



    Junction::Junction(
        const std::string& junctionId,
        const JunctionSideSection& leftSection,
        const JunctionSideSection& rightSection,
        const std::string& transitionSequence,
        size_t transitionLength,
        bool transitionIsKnown,
        uint32_t score,
        const std::string& xRef,
        const std::string& deletedTransposableElement,
        const std::string& knownUnderrepresentedRepeat,
        double frequencyInBaselineGenomeSet,
        const std::string& assembledSequence
        )
        :
        id_(junctionId),
        transitionSequence_(transitionSequence),
        transitionLength_(transitionLength),
        transitionIsKnown_(transitionIsKnown),
        score_(score),
        xRef_(xRef),
        deletedTransposableElement_(deletedTransposableElement),
        knownUnderrepresentedRepeat_(knownUnderrepresentedRepeat),
        frequencyInBaselineGenomeSet_(frequencyInBaselineGenomeSet),
        assembledSequence_(assembledSequence)
    {
        sideSections_[JUNCTION_LEFT_SIDE] = leftSection;
        sideSections_[JUNCTION_RIGHT_SIDE] = rightSection;
    }



    // Write a Junction.
    void Junction::write(std::ostream& out, 
        const reference::CrrFile& reference,
        size_t expectedAnnotationCount)  const
    {
        out << id_ << JunctionFile::SEP;
        sideSections_[JUNCTION_LEFT_SIDE].write(out, reference);
        out << JunctionFile::SEP;
        sideSections_[JUNCTION_RIGHT_SIDE].write(out, reference);
        out << JunctionFile::SEP;
        out << (sideSections_[JUNCTION_LEFT_SIDE].strand_ == 
                sideSections_[JUNCTION_RIGHT_SIDE].strand_ ? 'Y':'N') << JunctionFile::SEP;
        bool interChrom = sideSections_[JUNCTION_LEFT_SIDE].position_.chromosome_ != 
                            sideSections_[JUNCTION_RIGHT_SIDE].position_.chromosome_;
        out << (interChrom ? 'Y':'N') << JunctionFile::SEP;
        if (!interChrom)
            out << getDistance(); 
        out << JunctionFile::SEP;
        out << score_ << JunctionFile::SEP;

        out << (transitionIsKnown_ ? 'Y' : 'N') << JunctionFile::SEP;
        out << transitionSequence_ << JunctionFile::SEP;
        if ( transitionLength_ != (size_t)(-1) ) out << transitionLength_;
        out << JunctionFile::SEP;

        out << sideSections_[JUNCTION_LEFT_SIDE].repeatClassification_ << JunctionFile::SEP;
        out << sideSections_[JUNCTION_RIGHT_SIDE].repeatClassification_ << JunctionFile::SEP;

        out << sideSections_[JUNCTION_LEFT_SIDE].genes_ << JunctionFile::SEP;
        out << sideSections_[JUNCTION_RIGHT_SIDE].genes_ << JunctionFile::SEP;

        out << xRef_ << JunctionFile::SEP;
        out << deletedTransposableElement_ << JunctionFile::SEP;
        out << knownUnderrepresentedRepeat_ << JunctionFile::SEP;
        out << boost::format("%.2f") % frequencyInBaselineGenomeSet_ << JunctionFile::SEP;
        out << assembledSequence_;

        if (0 != expectedAnnotationCount)
        {
            CGA_ASSERT(annotations_.size() == expectedAnnotationCount);
            BOOST_FOREACH(const std::string& ann, annotations_)
            {
                out << JunctionFile::SEP << ann;
            }
        }
    }

    size_t Junction::getDistance() const 
    {
        const reference::Location& l0 = sideSections_[JUNCTION_LEFT_SIDE].position_;
        const reference::Location& l1 = sideSections_[JUNCTION_RIGHT_SIDE].position_;
        CGA_ASSERT_EQ(l0.chromosome_, l1.chromosome_);
        CGA_ASSERT_MSG(l1.offset_>=l0.offset_, "The junction is not canonicalized: "<<CGA_VOUT(id_));
        return l1.offset_ - l0.offset_;
    }

    // Read a junction file.
    JunctionFile::JunctionFile(const std::string& name, 
        const reference::CrrFile& reference)
    {
        read(name, reference);
    }

    class JunctionFileRecord {
        friend std::ostream& operator<< (std::ostream& out, const JunctionFileRecord& r);
    public:
        JunctionFileRecord():
            leftStrand_('~'),
            leftPosition_(-1),
            leftLength_(-1),
            rightStrand_('~'),
            rightPosition_(-1),
            rightLength_(-1),
            strandConsistent_('x'),
            interchromosomal_('x'),
            transitionIsKnown_('x'),
            transitionLength_(-1),
            distance_(-1),
            score_(0)
        {}

        void initParser(util::DelimitedFile &delimitedFile, std::vector<std::string>& annotationHeaders) {
            delimitedFile.addField(util::StringField("id",&id_));

            delimitedFile.addField(util::StringField("leftChr",&leftChr_));
            delimitedFile.addField(util::ValueField<int32_t>("leftPosition",&leftPosition_));
            delimitedFile.addField(util::CharField("leftStrand",&leftStrand_));
            delimitedFile.addField(util::ValueField<int32_t>("leftLength",&leftLength_));

            delimitedFile.addField(util::StringField("rightChr",&rightChr_));
            delimitedFile.addField(util::ValueField<int32_t>("rightPosition",&rightPosition_));
            delimitedFile.addField(util::CharField("rightStrand",&rightStrand_));
            delimitedFile.addField(util::ValueField<int32_t>("rightLength",&rightLength_));

            delimitedFile.addField(util::CharField("strandConsistent",&strandConsistent_));
            delimitedFile.addField(util::CharField("interchromosomal",&interchromosomal_));

            delimitedFile.addField(util::ValueField<int32_t>("distance",&distance_)
                .exception("", 0));
            delimitedFile.addField(util::ValueField<uint32_t>("discordantMatePairAlignments",&score_));

            delimitedFile.addField(util::CharField("junctionSequenceResolved",&transitionIsKnown_));

            delimitedFile.addField(util::StringField("transitionSequence",&transitionSequence_));
            delimitedFile.addField(
                util::ValueField<int32_t>("transitionLength",&transitionLength_).exception("", -1));

            delimitedFile.addField(
                util::StringField("leftRepeatClassification",&leftRepeatClassification_));
            delimitedFile.addField(
                util::StringField("rightRepeatClassification",&rightRepeatClassification_));

            delimitedFile.addField(util::StringField("leftGenes",&leftGenes_));
            delimitedFile.addField(util::StringField("rightGenes",&rightGenes_));
            delimitedFile.addField(util::StringField("xRef",&xRef_));

            delimitedFile.addField(
                util::StringField("deletedTransposableElement",&deletedTransposableElement_));
            delimitedFile.addField(
                util::StringField("knownUnderrepresentedRepeat",&knownUnderrepresentedRepeat_));
            delimitedFile.addField(
                util::ValueField<double>("frequencyInBaselineGenomeSet",&frequencyInBaselineGenomeSet_));
            delimitedFile.addField(
                util::StringField("assembledSequence",&assembledSequence_));

            // Add string parsers for all annotations
            size_t lastMandatoryField = delimitedFile.getFieldOffset("assembledSequence");
            const std::vector<std::string>& hdrs = delimitedFile.getColumnHeaders();
            annotationHeaders.clear();
            annotationHeaders.reserve(hdrs.size() - lastMandatoryField - 1);
            annotations_.resize(hdrs.size() - lastMandatoryField - 1);
            for (size_t jj = 0, ii = lastMandatoryField + 1; ii < hdrs.size(); ++ii, ++jj)
            {
                annotationHeaders.push_back(hdrs[ii]);
                delimitedFile.addField(
                        util::StringField(hdrs[ii], &annotations_[jj]));
            }
        }
        std::string     id_;

        char            leftStrand_;
        std::string     leftChr_;
        int32_t         leftPosition_;
        int32_t         leftLength_;
        std::string     leftRepeatClassification_;
        std::string     leftGenes_;

        char            rightStrand_;
        std::string     rightChr_;
        int32_t         rightPosition_;
        int32_t         rightLength_;
        std::string     rightRepeatClassification_;
        std::string     rightGenes_;

        char            strandConsistent_;
        char            interchromosomal_;

        char            transitionIsKnown_;
        std::string     transitionSequence_;
        int32_t         transitionLength_;

        int32_t         distance_;
        uint32_t        score_;

        std::string     xRef_;
        std::string     deletedTransposableElement_;
        std::string     knownUnderrepresentedRepeat_;
        double          frequencyInBaselineGenomeSet_;
        std::string     assembledSequence_;

        std::vector<std::string> annotations_;
    };


    std::ostream& operator<< (std::ostream& out, const JunctionFileRecord& r) {
        const std::string& sep = JunctionFile::SEP;
        out << r.id_
            << sep << r.leftStrand_
            << sep << r.leftChr_
            << sep << r.leftPosition_
            << sep << r.leftLength_

            << sep << r.rightStrand_
            << sep << r.rightChr_
            << sep << r.rightPosition_
            << sep << r.rightLength_

            << sep << r.strandConsistent_
            << sep << r.interchromosomal_
            << sep << r.distance_

            << sep << r.transitionIsKnown_
            << sep << r.transitionSequence_
            << sep << r.transitionLength_

            << sep << r.leftRepeatClassification_
            << sep << r.rightRepeatClassification_

            << sep << r.leftGenes_
            << sep << r.rightGenes_

            << sep << r.score_

            << sep << r.xRef_
            << sep << r.deletedTransposableElement_
            << sep << r.knownUnderrepresentedRepeat_
            << sep << r.frequencyInBaselineGenomeSet_
            << sep << r.assembledSequence_;

        return out;
    }


    // Read a junction file.
    // This should no be used to read a file written using writeSortedWithReplication!
    void JunctionFile::read(const std::string& name, const reference::CrrFile& reference)
    {
        fileName_ = name;
        boost::shared_ptr<std::istream> junctionsFileStream = 
                util::InputStream::openCompressedInputStreamByExtension(name);
        util::DelimitedFile  junctionsFile(*junctionsFileStream, fileName_);

        metadata_ = junctionsFile.getMetadata();

        JunctionFileRecord junctionRecord;
        junctionRecord.initParser(junctionsFile, annotationHeaders_);


        for (bool notEof = junctionsFile.next(); notEof; notEof = junctionsFile.next()) 
        {
            CGA_ASSERT_MSG( 
                    (junctionRecord.leftStrand_ == junctionRecord.rightStrand_ 
                        && junctionRecord.strandConsistent_ == 'Y') 
                    || (junctionRecord.leftStrand_ != junctionRecord.rightStrand_ 
                        && junctionRecord.strandConsistent_ == 'N'),
                    "Incorrect value of 'StrandConsistent'. Record: " << junctionRecord
            );
            CGA_ASSERT_MSG(
                    (junctionRecord.leftChr_ == junctionRecord.rightChr_ 
                        && junctionRecord.interchromosomal_ == 'N')
                    || (junctionRecord.leftChr_ != junctionRecord.rightChr_ 
                        && junctionRecord.interchromosomal_ == 'Y'),
                    "Incorrect value of 'Interchromosomal'. Record: " << junctionRecord
            );
            CGA_ASSERT_MSG(
                ((junctionRecord.interchromosomal_ == 'N' &&
                    junctionRecord.leftPosition_ <= junctionRecord.rightPosition_)
                || (junctionRecord.interchromosomal_ == 'Y')),
                "Incorrect (non-canonical) record form. Record: " << junctionRecord
            );
            CGA_ASSERT_MSG(
                ((junctionRecord.interchromosomal_ == 'N' &&
                 junctionRecord.rightPosition_ - junctionRecord.leftPosition_ == junctionRecord.distance_)
                 || junctionRecord.interchromosomal_ == 'Y'),
                "Incorrect value of 'Distance'. Record: " << junctionRecord
            );
            CGA_ASSERT_MSG(
                ((!junctionRecord.transitionSequence_.empty() && junctionRecord.transitionIsKnown_ == 'Y')
                || junctionRecord.transitionSequence_.empty()),
                "Incorrect value of 'TransitionSequence'. Record: " << junctionRecord
            );
            CGA_ASSERT_MSG(
                (junctionRecord.transitionSequence_.length() == size_t(junctionRecord.transitionLength_)
                || junctionRecord.transitionSequence_.empty()),
                "Incorrect value of 'TransitionLength'. Record: " << junctionRecord
            );

            Junction junction(
                junctionRecord.id_,
                JunctionSideSection(junctionRecord.leftStrand_=='+' ? JUNCTION_PLUS_STRAND
                                        : junctionRecord.leftStrand_=='-' ? JUNCTION_MINUS_STRAND
                                        : JUNCTION_UNKNOWN_STRAND,
                                    reference::Location(
                                        reference.getChromosomeId(junctionRecord.leftChr_),
                                                                junctionRecord.leftPosition_),
                                    junctionRecord.leftLength_,
                                    junctionRecord.leftRepeatClassification_,
                                    junctionRecord.leftGenes_
                ),

                JunctionSideSection(junctionRecord.rightStrand_=='+' ? JUNCTION_PLUS_STRAND
                                        : junctionRecord.rightStrand_=='-' ? JUNCTION_MINUS_STRAND
                                        : JUNCTION_UNKNOWN_STRAND,
                                    reference::Location(
                                        reference.getChromosomeId(junctionRecord.rightChr_),
                                                                junctionRecord.rightPosition_),
                                    junctionRecord.rightLength_,
                                    junctionRecord.rightRepeatClassification_,
                                    junctionRecord.rightGenes_
                ),

                junctionRecord.transitionSequence_,
                junctionRecord.transitionLength_,
                junctionRecord.transitionIsKnown_=='Y',

                junctionRecord.score_,
                junctionRecord.xRef_,
                junctionRecord.deletedTransposableElement_,
                junctionRecord.knownUnderrepresentedRepeat_,
                junctionRecord.frequencyInBaselineGenomeSet_,
                junctionRecord.assembledSequence_
            );
            junction.annotations_ = junctionRecord.annotations_;
            add(junction);
        }
    }


    // Write a JunctionFile in memory to a file.
    // The name passed in as an argument excludes the file extension
    // and the period which precedes it.
    void JunctionFile::write(
        const std::string& name, 
        const reference::CrrFile& reference
        )
    {
        // Open the file and write the header.
        const std::string fileName = name;
        cgatools::util::OutputStream out(fileName);

        out << metadata_ << '\n';
        out << header_;
        BOOST_FOREACH(const std::string& hdr, annotationHeaders_)
        {
            out << SEP << hdr;
        }
        out << '\n';

        // Write all the junctions.
        BOOST_FOREACH(const Junction& junction, junctions_) {
            junction.write(out, reference, annotationHeaders_.size());
            out << '\n';
        }
    }

    void JunctionFile::add(const Junction& junction)
    {
        junctions_.push_back(junction);
    }

}}
