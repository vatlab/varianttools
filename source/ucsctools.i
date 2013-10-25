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
// for tabix
#include "bgzf.h"
#include "tabix.h"
#include "knetfile.h"
// for vcf
#include "vcf.h"
// for bam file
#include "bamFile.h"
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
        printf("Header: (excluding INFO and FORMAT lines)\n");
        char *tok = strtok(vcff->headerString, "\n");
        while (tok != NULL) {
            if ( !(startsWith("##INFO", (char *)tok) || startsWith("##fileformat", (char *)tok)
              || startsWith("##FORMAT", (char *)tok) || startsWith("#CHROM", (char *)tok)))
                printf("%-23s %s\n", "", tok);
            tok = strtok(NULL, "\n");
        }
        printf("\n");
        printf("Available fields (with type VARCHAR if unspecified or all=1):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "chr (1, chrom)", "chromosome");
        printf("%-23s %s\n", "pos (2, INTEGER)", "position (1-based)");
        printf("%-23s %s\n", "name (3)", "name of variant");
        printf("%-23s %s\n", "ref (4)", "reference allele");
        printf("%-23s %s\n", "alt (5)", "alternative alleles");
        printf("%-23s %s\n", "qual (6)", "qual");
        printf("%-23s %s\n", "filter (7)", "filter");
        printf("%-23s %s\n", "info (8, default)", "variant info fields");
        //
        struct vcfInfoDef * def = NULL;
        for (def = vcff->infoDefs; def != NULL; def = def->next) {
            std::string typestring = "";
            if (def->type == vcfInfoFlag)
                typestring = " (INTEGER, flag)";
            else if (def->type == vcfInfoInteger)
                typestring = " (INTEGER)";
            else if (def->type == vcfInfoFloat)
                typestring = " (FLOAT)";
            sprintf(buf, "info.%s%s", def->key, typestring.c_str());
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
                std::string typestring = "";
                if (def->type == vcfInfoFlag)
                    typestring = " (INTEGER, flag)";
                else if (def->type == vcfInfoInteger)
                    typestring = " (INTEGER)";
                else if (def->type == vcfInfoFloat)
                typestring = " (FLOAT)";
                sprintf(buf, "%s.%s%s", vcff->genotypeIds[i], def->key, typestring.c_str());
                printf("%-23s %s for sample %s\n", buf, def->description, vcff->genotypeIds[i]);
            }
        }
        vcfFileFree(&vcff);
     } else if (endsWith((char *)track_file.c_str(), ".bam")) {
        if (!bamFileExists((char *)track_file.c_str()))
            return;
        char * filename;
        samfile_t * bamf = bamOpen((char *)track_file.c_str(), &filename);
        if (bamf == NULL)
            return;
        // bam_init_header_hash(((samfile_t *)bamf)->header);
        bam_index_t * idx = bam_index_load(filename);
        if (idx == NULL)
            return;
        struct bamChromInfo * chrom = bamChromList(bamf);
        printf("%-23s\n%s\n", "Header:", bamf->header->text);
        printf("%-23s %d\n", "Chrom size:", slCount(chrom));
        for (; chrom != NULL; chrom = chrom->next) {
        	printf("    %-19s %u\n", chrom->name, chrom->size);
        }
        printf("\n");
        /* the following piece of code can get number of mapped and unmapped reads
         * if we can include proper header file for khint_t etc
        size_t i;
        for (i = 0; i < idx->n; ++i) {
                khint_t k;
                khash_t(i) *h = idx->index[i];
                k = kh_get(i, h, BAM_MAX_BIN);
                if (k != kh_end(h))
                        printf("\t%llu\t%llu", (long long)kh_val(h, k).list[1].u, (long long)kh_val(h, k).list[1].v);
                else printf("\t0\t0");
                putchar('\n');
        }
        */
        /* available fields */
        printf("Available fields (with type VARCHAR if unspecified or all=1):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if depth is over 0, NULL otherwise");
        printf("%-23s %s\n", "coverage (INTEGER)", "Number of reads that cover the starting position of the variant");
        printf("%-23s %s\n", "calls", "nucleotide of the reads at the variant location");
        printf("%-23s %s\n", "reads", "nucleotide sequence around the variant location");
        printf("%-23s %s\n", "READS", "the same as reads but use colors for insertions and variant allele");
        printf("%-23s %s\n", "qual", "A list of phred base quality of reads at the location");
        printf("%-23s %s\n", "avg_qual (FLOAT)", "Average qual score of all reads");
        printf("%-23s %s\n", "mapq", "A list of phred base quality of alignment at the location");
        printf("%-23s %s\n", "avg_mapq (FLOAT)", "Average mapq score of all reads");
        //
        printf("\nTags and flag that can be outputed or used in filters, with values from the 1st record:\n");
        /* grab the first item and show its properties */
        bam1_t data;
        bam1_t * bam = &data;
        ZeroVar(bam);
        if (bam_read1(bamf->x.bam, bam) > 0) {
            // adapted from part of bam.c bam_format1:
            uint8_t *s = bam1_aux(bam);
            while (s < bam->data + bam->data_len)
            {
                uint8_t type, key[2];
                key[0] = s[0];
                key[1] = s[1];
                s += 2;
                type = *s;
                ++s;
                if (type == 'A') {
                    printf("%c%c                      A (char)   : %c\n", key[0], key[1], *s);
                    ++s;
                } else if (type == 'C') {
                    printf("%c%c                      C (int)    : %u\n", key[0], key[1], *(char*)s);
                    ++s;
                } else if (type == 'c') {
                    printf("%c%c                      c (int8)   : %d\n", key[0], key[1], *(int8_t*)s);
                    ++s;
                } else if (type == 'S') { 
                    printf("%c%c                      S (uint16) : %u\n", key[0], key[1], *(uint16_t*)s);
                    s += 2; 
                } else if (type == 's') {
                    printf("%c%c                      s (int16)  : %d\n", key[0], key[1], *(int16_t*)s);
                    s += 2;
                } else if (type == 'I') {
                    printf("%c%c                      I (uint32  : %u\n", key[0], key[1], *(uint32_t*)s);
                    s += 4;
                } else if (type == 'i') { 
                    printf("%c%c                      i (int32)  : %d\n", key[0], key[1], *(int32_t*)s);
                    s += 4;
                } else if (type == 'f') { 
                    printf("%c%c                      f (float)  : %g\n", key[0], key[1], *(float*)s);
                    s += 4; 
                } else if (type == 'd') { 
                    printf("%c%c                      d (double) : %lg\n", key[0], key[1], *(double*)s);
                    s += 8; 
                } else if (type == 'Z') {
                    printf("%c%c                      Z (string) : %s\n", key[0], key[1], (char *)s);
                    s += strlen((char *)s) + 1;
                } else if (type == 'H') {
                    printf("%c%c                      H (string) : %s\n", key[0], key[1], (char *)s);
                    s += strlen((char *)s) + 1;
                } else if (type == 'B') {
                    // get byte size
                    uint8_t subtype = *s;
                    // The letter can be one of `cCsSiIf', corresponding to int8 t (signed 8-bit
                    // integer), uint8 t (unsigned 8-bit integer), int16 t, uint16 t, int32 t, uint32 t and float,
                    int sz = 0;
                    if (subtype == 'c' || subtype == 'C')
                        sz = 1;
                    else if (subtype == 's' || subtype == 'S')
                        sz = 2;
                    else if (subtype == 'i' || subtype == 'I' || subtype == 'f')
                        sz = 4;
                    ++s;
                    int nItem = (int) *s;
                    s += 4 + nItem * sz;
                    printf("%c%c                      B (array of %c, comparison not supported)\n", key[0], key[1], subtype);
                }
            }
            printf("flag                    int flag   : 0x%X (%spaired, %smapped according to bits 1 & 3)\n",
                bam->core.flag, bam->core.flag & 1 == 0 ? "un" : "", bam->core.flag & 8 != 0 ? "un" : "");
            printf("\n");
        }
        printf("Parameters start (default to 0) and width (default to 5) can be used with reads (or "
            "READS) to adjust the window around variant, with syntax reads?start=-5&width=10. min_qual, "
            "min_mapq and TAG=VAL (or >, >=, <, <=, !=) can be used for all fields to limit the reads "
            "to the ones with mapq and qual scores that exceed the specified value, and tag satisfying "
            "specified conditions.\n");
    } else if (isBigWig((char *)track_file.c_str())) {
        struct bbiFile *bwf = bigWigFileOpen((char *)track_file.c_str());
        if (bwf == NULL)
            return;
        printf("%-23s %d\n", "Version:", bwf->version);
        printf("%-23s %llu\n", "Primary data size",
            bwf->unzoomedIndexOffset - bwf->unzoomedDataOffset);
        printf("%-23s %d\n", "Zoom levels:", bwf->zoomLevels);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
        printf("%-23s %d\n", "Chrom count:", slCount(chromList));
        printf("%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	printf("    %-19s %u\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bwf);
        printf("%-23s %llu\n", "Bases covered:", sum.validCount);
        printf("%-23s %f\n", "Mean:", sum.sumData/sum.validCount);
        printf("%-23s %f\n", "Min:", sum.minVal);
        printf("%-23s %f\n", "Max:", sum.maxVal);
        printf("%-23s %f\n", "std:",
            calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        printf("%-23s %d\n\n", "Number of fields:", 4);
        printf("Available fields (with type VARCHAR if unspecified or all=1):\n");
        printf("%-23s %s\n", "0 (INTEGER)", "1 if matched");
        printf("%-23s %s\n", "chrom (1)", "chromosome");
        printf("%-23s %s\n", "chromStart (2, INTEGER)", "start position (0-based)");
        printf("%-23s %s\n", "chromEnd (3, INTEGER)", "end position (1-based)");
        printf("%-23s %s\n", "value (4, FLOAT)", "value");
        bbiFileClose(&bwf);
     } else if (bigBedFileCheckSigs((char *)track_file.c_str())) {
        bbiFile * bbi = bigBedFileOpen((char *)track_file.c_str());
        if (bbi == NULL)
            return;

        printf("%-23s %d\n", "Version:", bbi->version);
        printf("%-23s %llu\n", "Item count:", bigBedItemCount(bbi));
        printf("%-23s %llu\n", "Primary data size:",
            bbi->unzoomedIndexOffset - bbi->unzoomedDataOffset);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
        printf("%-23s %d\n", "Zoom levels:", bbi->zoomLevels);
        printf("%-23s %d\n", "Chrom count:", slCount(chromList));
        printf("%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	printf("    %-19s %u\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bbi);
        printf("%-23s %llu\n", "Bases covered", sum.validCount);
        printf("%-23s %f\n", "Mean depth:", sum.sumData/sum.validCount);
        printf("%-23s %f\n", "Min depth:", sum.minVal);
        printf("%-23s %f\n", "Max depth:", sum.maxVal);
        printf("%-23s %f\n", "Std of depth:", calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        printf("%-23s %d\n\n", "Number of fields:", bbi->fieldCount);
        printf("Available fields (with type VARCHAR if unspecified or all=1):\n");
        struct asObject *as = bigBedAs(bbi);
        if (as != NULL)
        {
            size_t i = 1;
            struct asColumn *col;
            for (col = as->columnList; col != NULL; col = col->next, ++i)
            {
                struct asTypeInfo * ltype = col->lowType;
                std::string typestring = "";
                if (asTypesIsInt(ltype->type))
                    typestring = ", INTEGER";
                else if (asTypesIsFloating(ltype->type))
                    typestring = ", FLOAT";
                sprintf(buf, "%s (%lu%s)",  col->name, i, typestring.c_str());
                printf("%-23s %s\n", buf, wrap(col->comment, strlen(buf)).c_str());
            }
        }
        bbiFileClose(&bbi);
    } else {
		fprintf(stderr, "Unknown track file type. Only local or remote tabix-indexed "
			" vcf files with extension .vcf.gz, indexed BAM files, bigWig and bigBed files are "
			"supported.");
		return;
	}
}
%} 


void showTrack(std::string & track_file);
