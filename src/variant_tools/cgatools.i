%module cgatools
%{
#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/reference/CompactDnaSequence.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/CrrFileWriter.hpp"
#include "cgatools/reference/GeneDataStore.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"
#include "cgatools/reference/RepeatMaskerStore.hpp"
#include "cgatools/reference/range.hpp"
%}

// for uint16_t
%include "stdint.i"
%include "std_string.i"
%include "std_vector.i"

namespace std
{
    %template()         vector<string>; 
}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(cgatools::util::Exception e)
    {
        SWIG_exception(SWIG_ValueError, e.what());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


//%include "../cgatools/core.hpp"
//%include "../cgatools/util/Streams.hpp"
%ignore cgatools::util::operator==;
%ignore cgatools::util::operator!=;
%ignore cgatools::util::Md5Digest::Md5Digest();
%ignore cgatools::reference::Location::Location();
%ignore cgatools::reference::Range::Range();
%ignore cgatools::reference::Range::Range(const Location&, const Location&);
%ignore cgatools::reference::operator==(const Range&, const Range&);
%ignore cgatools::reference::operator!=(const Range&, const Range&);
%ignore cgatools::reference::operator<(const Range&, const Range&);
%ignore cgatools::reference::operator<=(const Range&, const Range&);
%ignore cgatools::reference::operator>(const Range&, const Range&);
%ignore cgatools::reference::operator>=(const Range&, const Range&);
%ignore cgatools::reference::AmbiguousRegion::AmbiguousRegion();
%ignore cgatools::reference::CrrFile::CrrFile();


%implicitconv Location;
%include "../cgatools/util/Md5.hpp"
%include "../cgatools/reference/range.hpp"
%include "../cgatools/reference/CompactDnaSequence.hpp"
%include "../cgatools/reference/CrrFile.hpp"


%inline %{

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/Exception.hpp"

#include <iostream>
#include <vector>

#include <Python.h>
#include <cstdlib>

#include <boost/format.hpp>
#include "boost/algorithm/string/trim.hpp"

std::string parseFastaHeader(const std::string& line)
{
    if (line.length() < 2 || line[0] != '>') {
        throw cgatools::util::Exception("expected FASTA header, found: " + line);
    }

    if (line.find('|') != std::string::npos)
        throw cgatools::util::Exception("NCBI-style fasta headers are not supported: " + line);

    std::string name = line.substr(1);
    boost::trim(name);
    return name;
}

bool fasta2crr(const std::vector<std::string> & fasta_files, const std::string & crr_file)
{
    cgatools::util::OutputStream out(crr_file);
    cgatools::reference::CrrFileWriter writer(&out);

    std::vector<std::string>::const_iterator it = fasta_files.begin();
    std::vector<std::string>::const_iterator it_end = fasta_files.end();
    for (; it != it_end; ++it) {
        std::cerr << "Reading from " << *it << std::endl;
        const boost::shared_ptr<std::istream> & in = cgatools::util::InputStream::openCompressedInputStreamByExtension(*it);

        std::string line;
        for(int lineNumber=1; InputStream::getline(*in, line); lineNumber++)
        {
            try
            {
                if (1 == lineNumber || '>' == line[0])
                {
                    std::string name = parseFastaHeader(line);
                    // here we only assume chrM is circular...
                    std::cerr << "Adding chromosome " << name << std::endl;
                    writer.newChromosome(name, name == "chrM");
                    continue;
                }

                writer.addSequence(line);
            }
            catch(std::exception& ee)
            {
                throw cgatools::util::Exception((boost::format("%s:%d: %s") %
                                 (*it) % lineNumber % ee.what()).str());
            }
        }
    }

    return 0;
}


/**
 * Right trims or left extend a variant.
 * This is adapted from VT variant_manip.cpp
 */
void right_trim_or_left_extend(cgatools::reference::CrrFile * crr, 
    std::vector<std::string>& alleles, uint32_t& pos1, int chrIdx, uint32_t& left_extended, uint32_t& right_trimmed)
{
    bool to_right_trim = true;
    bool to_left_extend = false;

    while (to_right_trim || to_left_extend)
    {
        //checks if right trimmable or left extendable
        to_right_trim = true;
        to_left_extend = false;
        for (size_t i=0; i<alleles.size(); ++i)
        {
            if (!alleles[i].empty())
            {
                if (alleles[0].at(alleles[0].size()-1) != alleles[i].at(alleles[i].size()-1))
                {
                    to_right_trim = false;
                    //do not break here!!! you need to check for empty alleles that might exist!!!
                }

                if (pos1==1 && alleles[i].size()==1)
                {
                    to_right_trim = false;
                    break;
                }
            }
            else
            {
                to_right_trim = false;
                // if crr is NULL (something wrong with reference genome match, do not left_extend)
                to_left_extend = (crr != NULL);
                break;
            }
        }

        if (to_right_trim)
        {
            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].erase(alleles[i].size()-1);
            }

            ++right_trimmed;
        }

        if (to_left_extend)
        {
            --pos1;
            int ref_len = 0;

            char base = crr->getBase(cgatools::reference::Location(chrIdx, pos1-1));

            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].insert(0, 1, base);
            }

            ++left_extended;
        }
    }
}

/**
 * Left trims a variant.
 */
void left_trim(std::vector<std::string>& alleles, uint32_t& pos1, uint32_t& left_trimmed)
{
    bool to_left_trim =  true;

    while (to_left_trim)
    {
        //checks if left trimmable.
        for (size_t i=0; i<alleles.size(); ++i)
        {
            if (alleles[i].size()==0 || alleles[i].at(0) != alleles[0].at(0))
            {
                to_left_trim = false;
                break;
            }
        }

        if (to_left_trim)
        {
            for (size_t i=0; i<alleles.size(); ++i)
            {
                alleles[i].erase(0, 1);
            }

            ++pos1;
            ++left_trimmed;
        }
    }
};


std::string normalize_variant(cgatools::reference::CrrFile * crr, PyObject * rec, size_t chr_idx, size_t pos_idx, size_t ref_idx, size_t alt_idx)
{
    //
    PyObject * chr_obj = PyList_GetItem(rec, chr_idx);
    if (! PyUnicode_Check(chr_obj)) {
        PyObject * chr_str = PyObject_Repr(chr_obj);
        return (boost::format("Incorrect type of passed chromosome: %s")
            % PyString_AsString(chr_str)).str();
    }

    const char * chr = PyUnicode_AsUTF8(chr_obj);
    if (chr == NULL)
         return("Failed to decode chromosome name.");
    PyObject * pos_obj = PyList_GetItem(rec, pos_idx);
    if (pos_obj == Py_None)
        return "Unrecognized position";
    bool int_pos = PyInt_Check(pos_obj);
    uint32_t pos = int_pos ? PyInt_AsLong(pos_obj) : atol(PyUnicode_AsUTF8(pos_obj));
    PyObject * ref_obj = PyList_GetItem(rec, ref_idx);
    PyObject * alt_obj = PyList_GetItem(rec, alt_idx);
    const char * ref = PyUnicode_AsUTF8(ref_obj);
    const char * alt = PyUnicode_AsUTF8(alt_obj);
    if (ref == NULL)
        return("Invalid reference allele");
    if (alt == NULL)
        return("Invalid alternative allele");

    //
    // first, if chr starts with chr, trim it. This should not be necessary because
    // the adjust function must have done that.
    if (strncmp(chr, "chr", 3) == 0) {
        PyList_SetItem(rec, chr_idx, PyUnicode_FromString(chr + 3));
        chr = PyUnicode_AsUTF8(PyList_GetItem(rec, chr_idx));
    }


    int chrIdx = -1;
    std::string msg;
    try {
        chrIdx = crr->getChromosomeId(chr);
    } catch (std::exception & e) {
        msg = (boost::format("Failed to get chromosome %s from reference genome") % chr).str();
    }
    // Check reference genome
    if (strlen(ref) >= 1 && ref[0] != '-') {
        try {
            if (strlen(ref) == 1) {
                char base = crr->getBase(cgatools::reference::Location(chrIdx, pos - 1));
                if (base != ref[0])
                    msg = (boost::format("Inconsistent base allele %s at %d on chromosome %s") % ref[0] % pos % chr).str();
            } else {
                std::string bases = crr->getSequence(cgatools::reference::Range(chrIdx, pos - 1, pos + strlen(ref)));
                if (strncmp(bases.c_str(), ref, strlen(ref)) != 0)
                    msg = (boost::format("Inconsistent base allele %s at %d on chromosome %s") % ref % pos % chr).str();
            }
        } catch (std::exception & e) {
            msg = (boost::format("Failed to get reference sequence: %s") % e.what()).str();
        }
    }

    std::vector<std::string> alleles;
    // check - and other characters
    for (size_t i = 0; i < 2; ++i) {
        std::string allele;
        const char * a = i == 0 ? ref : alt;
        for (; *a != '\0'; ++a) {
            switch (*a) {
            case '-':
                // pass
                break;
            case 'a':
                allele.append(1, 'A');
                break;
            case 'c':
                allele.append(1, 'C');
                break;
            case 't':
                allele.append(1, 'T');
                break;
            case 'g':
                allele.append(1, 'G');
                break;
            case 'A':
                allele.append(1, 'A');
                break;
            case 'C':
                allele.append(1, 'C');
                break;
            case 'T':
                allele.append(1, 'T');
                break;
            case 'G':
                allele.append(1, 'G');
                break;
            default:
                // the variant will be ignored in this case
                msg = (boost::format("Unrecognized allele %s") % (i == 0 ? ref : alt)).str();
                return msg;
            }
        }
        alleles.push_back(allele);
    }
    //
    uint32_t pos1 = pos;
    uint32_t left_extended = 0;
    uint32_t left_trimmed = 0;
    uint32_t right_trimmed = 0;
    // if there is something wrong with reference genome, DO NOT extend to the left
    right_trim_or_left_extend(msg.empty() ? crr : NULL,  alleles, pos1, chrIdx, left_extended, right_trimmed);

    left_trim(alleles, pos1, left_trimmed);

    if (left_extended > 1)
        msg = (boost::format("Variant %s_%d_%s/%s extended %d bp to the left.") % chr % pos % ref % alt % (pos - pos1 + 1)).str();

    // change values
    if (pos != pos1 || ! int_pos)
        PyList_SetItem(rec, pos_idx, PyInt_FromLong(pos1));
    if (alleles[0].empty()) {
        if (strlen(ref) != 1 || ref[0] != '-')
            PyList_SetItem(rec, ref_idx, PyUnicode_FromString("-"));
    } else {
        if (alleles[0].compare(ref) != 0)
            PyList_SetItem(rec, ref_idx, PyUnicode_FromString(alleles[0].c_str()));
    }
    if (alleles[1].empty()) {
        if (strlen(alt) != 1 || alt[0] != '-')
            PyList_SetItem(rec, alt_idx, PyUnicode_FromString("-"));
    } else {
        if (alleles[1].compare(alt) != 0)
            PyList_SetItem(rec, alt_idx, PyUnicode_FromString(alleles[1].c_str()));
    }

    // no error message
    return msg;
}

%}

bool fasta2crr(const std::vector<std::string> & fasta_files, const std::string & crr_file);
