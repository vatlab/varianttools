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

%include "std_string.i"

//%include "cgatools/core.hpp"
//%include "cgatools/util/Streams.hpp"
%ignore cgatools::util::operator==;
%ignore cgatools::util::operator!=;



%include "cgatools/util/Md5.hpp"
%include "cgatools/reference/range.hpp"
%include "cgatools/reference/CompactDnaSequence.hpp"
%include "cgatools/reference/CrrFile.hpp"
