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
extern "C" {
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
}
%}

%init
%{
%}

%include exception.i

%exception
{
    try
    {
        $function
    }
    catch(...)
    {
        SWIG_exception(SWIG_UnknownError, "Unknown runtime error happened.");
    }
}


%include "std_string.i"

%inline %{
void showTrack(const std::string & track_file)
{
    char buf[255];
    if (endsWith((char *)track_file.c_str(), ".vcf.gz")) {
        vcfFile * vcff = vcfTabixFileMayOpen((char *)track_file.c_str(),
             NULL, 0, 0, VCF_IGNORE_ERRS, 1);
        if (vcff == NULL)
             return;
        
        printf("%-23s %d\n\n", "Number of columns:", 8 + (vcff->genotypeCount > 0 ? vcff->genotypeCount + 1 : 0));
        printf("Available columns (with default type VARCHAR):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "1", "chromosome");
        printf("%-23s %s\n", "2 (INTEGER)", "position (1-based)");
        printf("%-23s %s\n", "3", "name of variant");
        printf("%-23s %s\n", "4", "reference allele");
        printf("%-23s %s\n", "5", "alternative alleles");
        printf("%-23s %s\n", "6", "qual");
        printf("%-23s %s\n", "7", "filter");
        printf("%-23s %s\n", "8 (default)", "variant info fields");
        printf("%-23s %s\n", "9", "genotype format");
        //
        for (size_t i = 0; i < vcff->genotypeCount; ++i)
            printf("%-23lu %s %s\n", 10 + i, "genotype for sample", vcff->genotypeIds[i]);
        //
        printf("%-23s %s\n", "'chr'", "chromosome");
        printf("%-23s %s\n", "'chrom'", "chromosome");
        printf("%-23s %s\n", "'pos'", "position (1-based)");
        printf("%-23s %s\n", "'name'", "name of variant");
        printf("%-23s %s\n", "'ref'", "reference allele");
        printf("%-23s %s\n", "'alt'", "alternative alleles");
        printf("%-23s %s\n", "'qual'", "qual");
        printf("%-23s %s\n", "'filter'", "filter");
        printf("%-23s %s\n", "'info'", "variant info fields");
        //
        struct vcfInfoDef * def = NULL;
        for (def = vcff->infoDefs; def != NULL; def = def->next) {
            sprintf(buf, "'info.%s'", def->key);
            printf("%-23s %s\n", buf, def->description);
        }
        printf("%-23s %s\n", "'format'", "genotype format");
        for (size_t i = 0; i < vcff->genotypeCount; ++i) {
            sprintf(buf, "'%s'",  vcff->genotypeIds[i]);
            printf("%-23s %s%s\n", buf, "genotype for sample ",
                vcff->genotypeIds[i]);
            // for all format fields
            struct vcfInfoDef * def = NULL;
            for (def = vcff->gtFormatDefs; def != NULL; def = def->next) {
                sprintf(buf, "'%s.%s'", vcff->genotypeIds[i], def->key);
                printf("%-23s %s for sample %s\n", buf, def->description, vcff->genotypeIds[i]);
            }
        }
     } else if (isBigWig((char *)track_file.c_str())) {
        printf("%-23s %d\n\n", "Number of columns:", 4);
        printf("Available columns (with default type VARCHAR):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "1", "chromosome");
        printf("%-23s %s\n", "2 (INTEGER)", "start position (0-based)");
        printf("%-23s %s\n", "3 (INTEGER)", "end position (1-based)");
        printf("%-23s %s\n", "4", "value");
        printf("%-23s %s\n", "'chr' ('chrom')", "chromosome");
        printf("%-23s %s\n", "'start' ('chromStart')", "start position (0-based)");
        printf("%-23s %s\n", "'end' ('chromEnd')", "end position (1-based)");
        printf("%-23s %s\n", "'value'", "value");
     } else {
        bbiFile * bbf = bigBedFileOpen((char *)track_file.c_str());
        if (bbf == NULL)
            return;
        printf("%-23s %d\n\n", "Number of columns:", bbf->fieldCount);
        printf("%-23s %s\n", "0 (INTEGER, default)", "1 if matched");
        printf("%-23s %s\n", "1", "chromosome");
        printf("%-23s %s\n", "2 (INTEGER)", "start position (0-based)");
        printf("%-23s %s\n", "3 (INTEGER)", "end position (1-based)");
        if (bbf->fieldCount > 3)
            printf("%-23s %s\n", "4", "name");
        if (bbf->fieldCount > 4)
            printf("%-23s %s\n", "5", "score");
        if (bbf->fieldCount > 5)
            printf("%-23s %s\n", "6", "strand");
        if (bbf->fieldCount > 6)
            printf("%-23s %s\n", "7", "thickStart");
        if (bbf->fieldCount > 7)
            printf("%-23s %s\n", "8", "thickEnd");
        if (bbf->fieldCount > 8)
            printf("%-23s %s\n", "9", "itemRgb");
        if (bbf->fieldCount > 9)
            printf("%-23s %s\n", "10", "blockCount");
        if (bbf->fieldCount > 10)
            printf("%-23s %s\n", "11", "blockSize");
        if (bbf->fieldCount > 11)
            printf("%-23s %s\n", "12", "blockStarts");
        //
        printf("%-23s %s\n", "'chr' ('chrom')", "chromosome");
        printf("%-23s %s\n", "'start' ('chromStart')", "start position (0-based)");
        printf("%-23s %s\n", "'end' ('chromEnd')", "end position (1-based)");
        if (bbf->fieldCount > 3)
            printf("%-23s %s\n", "'name'", "name");
        if (bbf->fieldCount > 4)
            printf("%-23s %s\n", "'score'", "score");
        if (bbf->fieldCount > 5)
            printf("%-23s %s\n", "'strand'", "strand");
        if (bbf->fieldCount > 6)
            printf("%-23s %s\n", "'thickStart'", "thickStart");
        if (bbf->fieldCount > 7)
            printf("%-23s %s\n", "'thickEnd'", "thickEnd");
        if (bbf->fieldCount > 8)
            printf("%-23s %s\n", "'itemRgb'", "itemRgb");
        if (bbf->fieldCount > 9)
            printf("%-23s %s\n", "'blockCount'", "blockCount");
        if (bbf->fieldCount > 10)
            printf("%-23s %s\n", "'blockSize'", "blockSize");
        if (bbf->fieldCount > 11)
            printf("%-23s %s\n", "'blockStarts'", "blockStarts");
    }
}
%} 


void showTrack(std::string & track_file);
