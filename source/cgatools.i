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
