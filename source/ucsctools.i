//!/usr/bin/env python2.7
//
// $File: ucsctools.i $
// $LastChangedDate: 2013-03-14 13:03:32 -0500 (Thu, 14 Mar 2013) $
// $Rev: 1737 $
//
// This file is part of variant_tools, a software application to annotate,
// summarize, and filter variants for next-gen sequencing ananlysis.
// Please visit http://varianttools.sourceforge.net for details.
//
// Copyright (C) 2011 Gao Wang (wangow@gmail.com) and Bo Peng (bpeng@mdanderson.org)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//


%module ucsctools

%{
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "sqlNum.h"
#include "asParse.h"
#include "udc.h"
#include "obscure.h"
#include "localmem.h"
#include "hmmstats.h"
#include "bigWig.h"
#include "bigBed.h"
#include "bgzf.h"
#include "tabix.h"
#include "knetfile.h"
#include "vcf.h"
%}

%init
%{
%}

%include exception.i

%exception
{
    {
        $function
    }
    catch(vtools::IndexError e)
    {
        SWIG_exception(SWIG_IndexError, e.message());
    }
    catch(vtools::ValueError e)
    {
        SWIG_exception(SWIG_ValueError, e.message());
    }
    catch(vtools::SystemError e)
    {
        SWIG_exception(SWIG_SystemError, e.message());
    }
    catch(vtools::RuntimeError e)
    {
        SWIG_exception(SWIG_RuntimeError, e.message());
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%include "std_string.i"

%inline %{
void showTrack(std::string & track_file)
{
	if (endsWith((char *)track_file.c_str(), ".vcf.gz")) {
        fprintf(stdout, "vcf format");
    } else if (isBigWig((char *)track_file.c_str())) {
        fprintf(stdout, "bigWig format");
    } else {
        fprintf(stdout, "bigBed format");
    }
//         info.file_type = VCFTABIX_FILE;
//         info.file = (void *)vcfTabixFileMayOpen((char *)track_file.c_str(),
//             NULL, 0, 0, VCF_IGNORE_ERRS, 1);
//         if (info.file == NULL) {
//             sqlite3_result_error(context, "cannot open file", -1);
//             return;
//         }
//         info.handler = vcfTabixTrack;
//         info.default_col = 8;
//         //
//         info.name_map["chr"] = 1;
//         info.name_map["chrom"] = 1;
//         info.name_map["pos"] = 2;
//         info.name_map["name"] = 3;
//         info.name_map["ref"] = 4;
//         info.name_map["alt"] = 5;
//         info.name_map["qual"] = 6;
//         info.name_map["filter"] = 7;
//         info.name_map["info"] = 8;
//         info.name_map["format"] = 9;
//         // info fields
//         struct vcfFile * vcff = (struct vcfFile *)info.file;
//         struct vcfInfoDef * def = NULL;
//         for (def = vcff->infoDefs; def != NULL; def = def->next)
//             info.name_map[std::string("info.") + def->key] = 8;
//         // sample
//         for (size_t i = 0; i < vcff->genotypeCount; ++i) {
//             info.name_map[vcff->genotypeIds[i]] = 10 + i;
//             // sample genotype field
//             struct vcfInfoDef * def = NULL;
//             for (def = vcff->gtFormatDefs; def != NULL; def = def->next)
//                 info.name_map[std::string(vcff->genotypeIds[i]) + std::string(".") + def->key] = 10 + i;
//         }
//         // get the first record
//         info.with_leading_chr = false;
//         struct vcfRecord * rec = vcff->records;
//         if (rec != NULL) {
//             if (strncmp(rec->chrom, "chr", 3) == 0)
//                 info.with_leading_chr = true;
//         }
//     } else if (isBigWig((char *)track_file.c_str())) {
//         info.file_type = BIGWIG_FILE;
//         info.file = (void *)bigWigFileOpen((char *)track_file.c_str());
//         if (info.file == NULL)
//             sqlite3_result_error(context, "cannot open file", -1);
//         info.handler = bigWigTrack;
//         info.default_col = 4;
//         //
//         info.name_map["chr"] = 1;
//         info.name_map["chrom"] = 1;
//         info.name_map["start"] = 2;
//         info.name_map["chromStart"] = 2;
//         info.name_map["end"] = 3;
//         info.name_map["chromEnd"] = 3;
//         info.name_map["value"] = 4;
//         //
//         // here we assume that we need to add "chr" to all bigWig files
//         info.with_leading_chr = true;
//     } else {
//         info.file_type = BIGBED_FILE;
//         info.file = (void *)bigBedFileOpen((char *)track_file.c_str());
//         if (info.file == NULL) {
//             sqlite3_result_error(context, "cannot open file", -1);
//             return;
//         }
//         info.handler = bigBedTrack;
//         info.default_col = 0;
//         //
//         info.name_map["chr"] = 1;
//         info.name_map["chrom"] = 1;
//         info.name_map["start"] = 2;
//         info.name_map["chromStart"] = 2;
//         info.name_map["end"] = 3;
//         info.name_map["chromEnd"] = 3;
//         info.name_map["name"] = 4;
//         info.name_map["score"] = 5;
//         info.name_map["strand"] = 6;
//         info.name_map["thickStart"] = 7;
//         info.name_map["thickEnd"] = 8;
//         info.name_map["itemRgb"] = 9;
//         info.name_map["blockCount"] = 10;
//         info.name_map["blockSizes"] = 11;
//         info.name_map["blockStarts"] = 12;
//         //
//         // here we assume that we need to add "chr" to all bigWig files
//         info.with_leading_chr = true;
//    }
}
%} 


void showTrack(std::string & track_file);
