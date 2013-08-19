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

std::string wrap(char * input, int start=24)
{
    std::string res;
    int count = start;
    if (count < 24)
        count = 24;
    for (char * p = input; *p != '\0'; ++p, ++count) {
        if (count > 78 && *p == ' ') {
            res += "\n                        ";
            count = 24;
        } else 
            res += *p;
    }
    return res;
}

void showTrack(const std::string & track_file)
{
    char buf[255];
    if (endsWith((char *)track_file.c_str(), ".vcf.gz")) {
        vcfFile * vcff = vcfTabixFileMayOpen((char *)track_file.c_str(),
             NULL, 0, 0, VCF_IGNORE_ERRS, 1);
        if (vcff == NULL)
             return;
        
        printf("%-23s VCF v%d.%d\n", "Version", vcff->majorVersion, vcff->minorVersion);
        printf("%-23s %d\n\n", "Number of fields:", 8 + (vcff->genotypeCount > 0 ? vcff->genotypeCount + 1 : 0));
        printf("Header: (exclude INFO and FORMAT lines)\n");
        char *tok = strtok(vcff->headerString, "\n");
        while (tok != NULL) {
            if ( !(startsWith("##INFO", (char *)tok) || startsWith("##fileformat", (char *)tok)
              || startsWith("##FORMAT", (char *)tok) || startsWith("#CHROM", (char *)tok)))
                printf("%-23s %s\n", "", tok);
            tok = strtok(NULL, "\n");
        }
        printf("\n");
        printf("Available columns (with default type VARCHAR):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "chr (1, chrom)", "chromosome");
        printf("%-23s %s\n", "pos (2 for INTEGER)", "position (1-based)");
        printf("%-23s %s\n", "name (3)", "name of variant");
        printf("%-23s %s\n", "ref (4)", "reference allele");
        printf("%-23s %s\n", "alt (5)", "alternative alleles");
        printf("%-23s %s\n", "qual (6)", "qual");
        printf("%-23s %s\n", "filter (7)", "filter");
        printf("%-23s %s\n", "info (8, default)", "variant info fields");
        //
        struct vcfInfoDef * def = NULL;
        for (def = vcff->infoDefs; def != NULL; def = def->next) {
            sprintf(buf, "info.%s", def->key);
            printf("%-23s %s\n", buf, def->description);
        }
        printf("%-23s %s\n", "format (9)", "genotype format");
        for (size_t i = 0; i < vcff->genotypeCount; ++i) {
            sprintf(buf, "%s (%lu)",  vcff->genotypeIds[i], 10 + i);
            printf("%-23s %s%s\n", buf, "genotype for sample ",
                vcff->genotypeIds[i]);
            // for all format fields
            struct vcfInfoDef * def = NULL;
            for (def = vcff->gtFormatDefs; def != NULL; def = def->next) {
                sprintf(buf, "%s.%s", vcff->genotypeIds[i], def->key);
                printf("%-23s %s for sample %s\n", buf, def->description, vcff->genotypeIds[i]);
            }
        }
     } else if (isBigWig((char *)track_file.c_str())) {
        struct bbiFile *bwf = bigWigFileOpen((char *)track_file.c_str());
        if (bwf == NULL)
            return;
        printf("%-23s %d\n", "Version:", bwf->version);
        printf("%-23s %lu\n", "Primary data size",
            bwf->unzoomedIndexOffset - bwf->unzoomedDataOffset);
        printf("%-23s %lu\n", "Zoom levels:", bwf->zoomLevels);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
        printf("%-23s %lu\n", "Chrom count:", slCount(chromList));
        printf("%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	printf("    %-19s %lu\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bwf);
        printf("%-23s %d\n", "Bases covered:", sum.validCount);
        printf("%-23s %f\n", "Mean:", sum.sumData/sum.validCount);
        printf("%-23s %f\n", "Min:", sum.minVal);
        printf("%-23s %f\n", "Max:", sum.maxVal);
        printf("%-23s %f\n", "std:",
            calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        printf("%-23s %d\n\n", "Number of fields:", 4);
        printf("Available columns (with default type VARCHAR):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "chrom (1)", "chromosome");
        printf("%-23s %s\n", "chromStart (2 as INTEGER)", "start position (0-based)");
        printf("%-23s %s\n", "chromEnd (3 as INTEGER)", "end position (1-based)");
        printf("%-23s %s\n", "value (4 as FLOAT)", "value");
        bbiFileClose(&bwf);
     } else {
        bbiFile * bbi = bigBedFileOpen((char *)track_file.c_str());
        if (bbi == NULL)
            return;

        printf("%-23s %d\n", "Version:", bbi->version);
        printf("%-23s %lu\n", "Item count:", bigBedItemCount(bbi));
        printf("%-23s %lu\n", "Primary data size:",
            bbi->unzoomedIndexOffset - bbi->unzoomedDataOffset);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
        printf("%-23s %d\n", "Zoom levels:", bbi->zoomLevels);
        printf("%-23s %d\n", "Chrom count:", slCount(chromList));
        printf("%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	printf("    %-19s %lu\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bbi);
        printf("%-23s %lu\n", "Bases covered", sum.validCount);
        printf("%-23s %f\n", "Mean depth:", sum.sumData/sum.validCount);
        printf("%-23s %f\n", "Min depth:", sum.minVal);
        printf("%-23s %f\n", "Max depth:", sum.maxVal);
        printf("%-23s %f\n", "Std of depth:", calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        printf("%-23s %d\n\n", "Number of fields:", bbi->fieldCount);
        printf("Available columns (with default type VARCHAR):\n");
        struct asObject *as = bigBedAs(bbi);
        if (as != NULL)
        {
            size_t i = 1;
            struct asColumn *col;
            for (col = as->columnList; col != NULL; col = col->next, ++i)
            {
                struct asTypeInfo * ltype = col->lowType;
                char * typestring = "";
                if (asTypesIsInt(ltype->type))
                    typestring = " as INTEGER";
                else if (asTypesIsFloating(ltype->type))
                    typestring = " as FLOAT";
                sprintf(buf, "%s (%lu%s)",  col->name, i, typestring);
                printf("%-23s %s\n", buf, wrap(col->comment, strlen(buf)).c_str());
            }
        }
        bbiFileClose(&bbi);
    }
}
%} 


void showTrack(std::string & track_file);
