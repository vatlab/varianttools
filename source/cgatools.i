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


//%include "cgatools/core.hpp"
//%include "cgatools/util/Streams.hpp"
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
%include "cgatools/util/Md5.hpp"
%include "cgatools/reference/range.hpp"
%include "cgatools/reference/CompactDnaSequence.hpp"
%include "cgatools/reference/CrrFile.hpp"


%inline %{

#include "cgatools/core.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/Exception.hpp"

#include <iostream>
#include <vector>

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

%}

bool fasta2crr(const std::vector<std::string> & fasta_files, const std::string & crr_file);
