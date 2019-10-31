//!/usr/bin/env python2.7
//
// $File: ucsctools.i $
// $LastChangedDate: 2013-03-14 13:03:32 -0500 (Thu, 14 Mar 2013) $
// $Rev: 1737 $
//
// This file is part of variant_tools, a software application to annotate,
// summarize, and filter variants for next-gen sequencing ananlysis.
// Please visit https://github.com/vatlab/varianttools for details.
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

%include "std_vector.i"
%include "std_string.i"

namespace std
{
    %template()         vector<string>;     /* e.g. infoFields */
}


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

void showTrack(const std::string & track_file, const std::string output_file = std::string())
{
    char buf[255];
    FILE * output = NULL;
    int output_style = 0;
    if (output_file.empty()) {
        output = stdout;
    } else {
        output = fopen(output_file.c_str(), "w");
        if (endsWith((char*)output_file.c_str(), ".fmt"))
            output_style = 2;
        else
            output_style = 1;
    }
    if (endsWith((char *)track_file.c_str(), ".vcf.gz") || 
        endsWith((char *)track_file.c_str(), ".vcf")) {
        /* if it is a local file, but does not exist. This is needed
           because vcfFileMayOpen will simply crash if file does not exist */
        if (track_file.find_first_of(':') == std::string::npos) {
            if (FILE * f = fopen(track_file.c_str(), "r"))
                fclose(f);
            else {
                fprintf(stderr, "File %s does not exist\n", track_file.c_str());
                return;
            }
        }
        vcfFile * vcff = NULL;
        if (endsWith((char *)track_file.c_str(), ".vcf"))
            vcff = vcfFileMayOpen((char *)track_file.c_str(),
                VCF_IGNORE_ERRS, 1, true);
        else 
            vcff = vcfTabixFileMayOpen((char *)track_file.c_str(),
                NULL, 0, 0, VCF_IGNORE_ERRS, 1);
        if (vcff == NULL)
             return;
        
        if (output_style == 0) {
            fprintf(output, "%-23s VCF v%d.%d\n", "Version", vcff->majorVersion, vcff->minorVersion);
            fprintf(output, "%-23s %d\n\n", "Number of fields:", 8 + (vcff->genotypeCount > 0 ? vcff->genotypeCount + 1 : 0));
            fprintf(output, "Header: (excluding INFO and FORMAT lines)\n");
            char *tok = strtok(vcff->headerString, "\n");
            while (tok != NULL) {
                if ( !(startsWith("##INFO", (char *)tok) || startsWith("##fileformat", (char *)tok)
                  || startsWith("##FORMAT", (char *)tok) || startsWith("#CHROM", (char *)tok)))
                    fprintf(output, "%-23s %s\n", "", tok);
                tok = strtok(NULL, "\n");
            }
            fprintf(output, "\n");
            fprintf(output, "Available fields (with type VARCHAR if unspecified or all=1):\n");
            fprintf(output, "%-23s %s\n", "0 (INTEGER)", "1 if matched");
            fprintf(output, "%-23s %s\n", "chr (1, chrom)", "chromosome");
            fprintf(output, "%-23s %s\n", "pos (2, INTEGER)", "position (1-based)");
            fprintf(output, "%-23s %s\n", "name (3)", "name of variant");
            fprintf(output, "%-23s %s\n", "ref (4)", "reference allele");
            fprintf(output, "%-23s %s\n", "alt (5)", "alternative alleles");
            fprintf(output, "%-23s %s\n", "qual (6)", "qual");
            fprintf(output, "%-23s %s\n", "filter (7)", "filter");
            fprintf(output, "%-23s %s\n", "info (8, default)", "variant info fields");
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
                fprintf(output, "%-23s %s\n", buf, def->description);
            }
            fprintf(output, "%-23s %s\n", "format (9)", "genotype format");
            for (size_t i = 0; i < vcff->genotypeCount; ++i) {
                sprintf(buf, "%s (%lu)",  vcff->genotypeIds[i], 10 + i);
                fprintf(output, "%-23s %s%s\n", buf, "genotype for sample ",
                    vcff->genotypeIds[i]);
                // for all format fields
                for (def = vcff->gtFormatDefs; def != NULL; def = def->next) {
                    std::string typestring = "";
                    if (def->type == vcfInfoFlag)
                        typestring = " (INTEGER, flag)";
                    else if (def->type == vcfInfoInteger)
                        typestring = " (INTEGER)";
                    else if (def->type == vcfInfoFloat)
                    typestring = " (FLOAT)";
                    sprintf(buf, "%s.%s%s", vcff->genotypeIds[i], def->key, typestring.c_str());
                    fprintf(output, "%-23s %s for sample %s\n", buf, def->description, vcff->genotypeIds[i]);
                }
            }
        } else if (output_style == 1) {
            /* output fields in .ann format */
            fprintf(output, "# %-23s VCF v%d.%d\n", "Version", vcff->majorVersion, vcff->minorVersion);
            fprintf(output, "# %-23s %d\n", "Number of fields:", 8 + (vcff->genotypeCount > 0 ? vcff->genotypeCount + 1 : 0));
            fprintf(output, "# Header: (excluding INFO and FORMAT lines)\n");
            char *tok = strtok(vcff->headerString, "\n");
            while (tok != NULL) {
                if ( !(startsWith("##INFO", (char *)tok) || startsWith("##fileformat", (char *)tok)
                  || startsWith("##FORMAT", (char *)tok) || startsWith("#CHROM", (char *)tok)))
                    fprintf(output, "# %-23s %s\n", "", tok);
                tok = strtok(NULL, "\n");
            }
            fprintf(output, "\n[chr]\nindex=1\ntype=VARCHAR(255)\nadj=RemoveLeading('chr')\ncomment=Chromosome\n\n");
            fprintf(output, "[pos]\nindex=2\ntype=INTEGER\ncomment=Position (1-based)\n\n");
			fprintf(output, "[name]\nindex=3\ntype=VARCHAR(255) NULL\ncomment=DB SNP ID\n\n");
			fprintf(output, "[ref]\nindex=4\ntype=VARCHAR(255) NOT NULL\ncomment=Reference allele (as on the + strand)\n\n");
			fprintf(output, "[alt]\nindex=5\ntype=VARCHAR(255) NOT NULL\nadj=CheckSplit()\ncomment=Alternative allele (as on the + strand)\n\n");
			fprintf(output, "[qual]\nindex=6\ntype=VARCHAR(255) NOT NULL\ncomment=Quality\n\n");
			fprintf(output, "[filter]\nindex=7\ntype=VARCHAR(255) NOT NULL\ncomment=Filter\n\n");
            //
            struct vcfInfoDef * def = NULL;
            for (def = vcff->infoDefs; def != NULL; def = def->next) {
                fprintf(output, "[%s]\nindex=8\n", def->key);
                if (def->type == vcfInfoFlag)
                    fprintf(output, "type=INTEGER\nadj=ExtractFlag('%s', ';')\n", def->key);
                else if (def->type == vcfInfoInteger)
                    fprintf(output, "type=INTEGER\nadj=ExtractValue('%s=', ';')\n", def->key);
                else if (def->type == vcfInfoFloat)
                    fprintf(output, "type=FLOAT\nadj=ExtractValue('%s=', ';')\n", def->key);
                else 
                    fprintf(output, "type=VARCHAR(255)\nadj=ExtractValue('%s=', ';')\n", def->key);
                fprintf(output, "comment=%s\n\n", def->description);
            }
        } else if (output_style == 2) {
            /* output fields in .fmt format */
            fprintf(output, "[format description]\n");
            fprintf(output, "description=Format file to import all variant and genotype info from %s\n", track_file.c_str());
            fprintf(output, "variant=chr,pos,ref,alt\n");
            /* import all variant info fields */
            struct vcfInfoDef * def = NULL;
            fprintf(output, "variant_info=name,qual,filter");
            for (def = vcff->infoDefs; def != NULL; def = def->next)
                fprintf(output, ",%s", def->key);
            fprintf(output, "\n");
            //
            if (vcff->genotypeCount > 0) {
                fprintf(output, "genotype=GT\n");
                fprintf(output, "genotype_info=");
                bool first_field = true;
                // skip GT
                def = vcff->gtFormatDefs;
                if (def != NULL)
                    def = def->next;
                for (; def != NULL; def = def->next) {
                    if (first_field) {
                        fprintf(output, "%s_geno", def->key);
                        first_field = false;
                    } else
                        fprintf(output, ",%s_geno", def->key);
                }
                fprintf(output, "\n\n");
            } else {
                fprintf(output, "genotype=\n");
                fprintf(output, "genotype_info=\n\n");
            }
            fprintf(output, "\n[chr]\nindex=1\ntype=VARCHAR(255)\nadj=RemoveLeading('chr')\ncomment=Chromosome\n\n");
            fprintf(output, "[pos]\nindex=2\ntype=INTEGER\ncomment=Position (1-based)\n\n");
			fprintf(output, "[name]\nindex=3\ntype=VARCHAR(255) NULL\ncomment=DB SNP ID\n\n");
			fprintf(output, "[ref]\nindex=4\ntype=VARCHAR(255) NOT NULL\ncomment=Reference allele (as on the + strand)\n\n");
			fprintf(output, "[alt]\nindex=5\ntype=VARCHAR(255) NOT NULL\nadj=CheckSplit()\ncomment=Alternative allele (as on the + strand)\n\n");
			fprintf(output, "[qual]\nindex=6\ntype=FLOAT\ncomment=phred-scaled quality score\n\n");
			fprintf(output, "[filter]\nindex=7\ntype=VARCHAR(255)\ncomment=Filter\n\n");
            //
            for (def = vcff->infoDefs; def != NULL; def = def->next) {
                fprintf(output, "[%s]\nindex=8\n", def->key);
                if (def->type == vcfInfoFlag)
                    fprintf(output, "type=INTEGER\nadj=ExtractFlag('%s', ';')\n", def->key);
                else if (def->type == vcfInfoInteger)
                    fprintf(output, "type=INTEGER\nadj=ExtractValue('%s=', ';')\n", def->key);
                else if (def->type == vcfInfoFloat)
                    fprintf(output, "type=FLOAT\nadj=ExtractValue('%s=', ';')\n", def->key);
                else 
                    fprintf(output, "type=VARCHAR(255)\nadj=ExtractValue('%s=', ';')\n", def->key);
                fprintf(output, "comment=%s\n\n", def->description);
            }
            /* if there is genotype */
            if (vcff->genotypeCount > 0) {
                fprintf(output, "[GT]\nindex=10:\ntype=INTEGER\nadj=VcfGenotype(default=('0',))\nfmt=GenoFormatter(style='vcf')\n");
                fprintf(output, "comment=Gentoype coded as 0 (ref ref), 1 (ref alt), 2 (alt alt) or -1 (alt1, alt2), assuming GT is the first FORMAT field in the .vcf file. Missing genotypes will be ignored.\n\n");
                def = vcff->gtFormatDefs;
                // skip GT field, which is handled separaely
                if (def != NULL)
                    def = def->next;
                for (; def != NULL; def = def->next) {
                    fprintf(output, "[%s_geno]\nindex=9,10:\n", def->key);
                    fprintf(output, "adj=FieldFromFormat('%s', ':')\n", def->key);
                    if (def->type == vcfInfoFlag)
                        fprintf(output, "type=INTEGER\n");
                    else if (def->type == vcfInfoInteger)
                        fprintf(output, "type=INTEGER\n");
                    else if (def->type == vcfInfoFloat)
                        fprintf(output, "type=FLOAT\n");
                    else
                        fprintf(output, "type=VARCHAR(255)\n");
                    fprintf(output, "comment=genotype field %s: %s\n\n", def->key, def->description);
                }
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
        fprintf(output, "%-23s\n%s\n", "Header:", bamf->header->text);
        fprintf(output, "%-23s %d\n", "Chrom size:", slCount(chrom));
        for (; chrom != NULL; chrom = chrom->next) {
        	fprintf(output, "    %-19s %u\n", chrom->name, chrom->size);
        }
        fprintf(output, "\n");
        /* the following piece of code can get number of mapped and unmapped reads
         * if we can include proper header file for khint_t etc
        size_t i;
        for (i = 0; i < idx->n; ++i) {
                khint_t k;
                khash_t(i) *h = idx->index[i];
                k = kh_get(i, h, BAM_MAX_BIN);
                if (k != kh_end(h))
                        fprintf(output, "\t%llu\t%llu", (long long)kh_val(h, k).list[1].u, (long long)kh_val(h, k).list[1].v);
                else fprintf(output, "\t0\t0");
                putchar('\n');
        }
        */
        /* available fields */
        fprintf(output, "Available fields (with type VARCHAR if unspecified or all=1):\n");
        fprintf(output, "%-23s %s\n", "0 (INTEGER)", "1 if depth is over 0, NULL otherwise");
        fprintf(output, "%-23s %s\n", "coverage (INTEGER)", "Number of reads that cover the starting position of the variant");
        fprintf(output, "%-23s %s\n", "calls", "nucleotide of the reads at the variant location");
        fprintf(output, "%-23s %s\n", "reads", "nucleotide sequence around the variant location");
        fprintf(output, "%-23s %s\n", "qual", "A list of phred base quality of reads at the location");
        fprintf(output, "%-23s %s\n", "avg_qual (FLOAT)", "Average qual score of all reads");
        fprintf(output, "%-23s %s\n", "mapq", "A list of phred base quality of alignment at the location");
        fprintf(output, "%-23s %s\n", "avg_mapq (FLOAT)", "Average mapq score of all reads");
        //
        fprintf(output, "\nTags and flag that can be outputed or used in filters, with values from the 1st record:\n");
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
                    fprintf(output, "%c%c                      A (char)   : %c\n", key[0], key[1], *s);
                    ++s;
                } else if (type == 'C') {
                    fprintf(output, "%c%c                      C (int)    : %u\n", key[0], key[1], *(char*)s);
                    ++s;
                } else if (type == 'c') {
                    fprintf(output, "%c%c                      c (int8)   : %d\n", key[0], key[1], *(int8_t*)s);
                    ++s;
                } else if (type == 'S') { 
                    fprintf(output, "%c%c                      S (uint16) : %u\n", key[0], key[1], *(uint16_t*)s);
                    s += 2; 
                } else if (type == 's') {
                    fprintf(output, "%c%c                      s (int16)  : %d\n", key[0], key[1], *(int16_t*)s);
                    s += 2;
                } else if (type == 'I') {
                    fprintf(output, "%c%c                      I (uint32  : %u\n", key[0], key[1], *(uint32_t*)s);
                    s += 4;
                } else if (type == 'i') { 
                    fprintf(output, "%c%c                      i (int32)  : %d\n", key[0], key[1], *(int32_t*)s);
                    s += 4;
                } else if (type == 'f') { 
                    fprintf(output, "%c%c                      f (float)  : %g\n", key[0], key[1], *(float*)s);
                    s += 4; 
                } else if (type == 'd') { 
                    fprintf(output, "%c%c                      d (double) : %lg\n", key[0], key[1], *(double*)s);
                    s += 8; 
                } else if (type == 'Z') {
                    fprintf(output, "%c%c                      Z (string) : %s\n", key[0], key[1], (char *)s);
                    s += strlen((char *)s) + 1;
                } else if (type == 'H') {
                    fprintf(output, "%c%c                      H (string) : %s\n", key[0], key[1], (char *)s);
                    s += strlen((char *)s) + 1;
                } else if (type == 'B') {
                    // get byte size
                    uint8_t subtype = *s;
                    // The letter can be one of cCsSiIf, corresponding to int8 t (signed 8-bit
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
                    fprintf(output, "%c%c                      B (array of %c, comparison not supported)\n", key[0], key[1], subtype);
                }
            }
            fprintf(output, "flag                    int flag   : 0x%X (%spaired, %smapped according to bits 1 & 3)\n",
                bam->core.flag, bam->core.flag & 1 == 0 ? "un" : "", bam->core.flag & 8 != 0 ? "un" : "");
            fprintf(output, "\n");
        }
        fprintf(output, "Parameters start (default to 0), width (default to 5) and color (default to 0) can be "
            "used with reads to adjust the window around variant, and use colors for insertions and variant "
            "allele, with syntax reads?start=-5&width=10&color=1. min_qual, min_mapq, strand and TAG=VAL (or >, "
            ">=, <, <=, !=) can be used for all fields to limit the reads to the ones with mapq and qual scores that "
	    "exceed the specified value, on specific strand (+ or -) and tag satisfying specified conditions. Parameter "
            "limit limits the number of reads or calls to display if the depth of coverage is high.\n");
    } else if (isBigWig((char *)track_file.c_str())) {
        struct bbiFile *bwf = bigWigFileOpen((char *)track_file.c_str());
        if (bwf == NULL)
            return;
        fprintf(output, "%-23s %d\n", "Version:", bwf->version);
        fprintf(output, "%-23s %llu\n", "Primary data size",
            bwf->unzoomedIndexOffset - bwf->unzoomedDataOffset);
        fprintf(output, "%-23s %d\n", "Zoom levels:", bwf->zoomLevels);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bwf);
        fprintf(output, "%-23s %d\n", "Chrom count:", slCount(chromList));
        fprintf(output, "%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	fprintf(output, "    %-19s %u\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bwf);
        fprintf(output, "%-23s %llu\n", "Bases covered:", sum.validCount);
        fprintf(output, "%-23s %f\n", "Mean:", sum.sumData/sum.validCount);
        fprintf(output, "%-23s %f\n", "Min:", sum.minVal);
        fprintf(output, "%-23s %f\n", "Max:", sum.maxVal);
        fprintf(output, "%-23s %f\n", "std:",
            calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        fprintf(output, "%-23s %d\n\n", "Number of fields:", 4);
        fprintf(output, "Available fields (with type VARCHAR if unspecified or all=1):\n");
        fprintf(output, "%-23s %s\n", "0 (INTEGER)", "1 if matched");
        fprintf(output, "%-23s %s\n", "chrom (1)", "chromosome");
        fprintf(output, "%-23s %s\n", "chromStart (2, INTEGER)", "start position (0-based)");
        fprintf(output, "%-23s %s\n", "chromEnd (3, INTEGER)", "end position (1-based)");
        fprintf(output, "%-23s %s\n", "value (4, FLOAT)", "value");
        bbiFileClose(&bwf);
     } else if (bigBedFileCheckSigs((char *)track_file.c_str())) {
        bbiFile * bbi = bigBedFileOpen((char *)track_file.c_str());
        if (bbi == NULL)
            return;

        fprintf(output, "%-23s %d\n", "Version:", bbi->version);
        fprintf(output, "%-23s %llu\n", "Item count:", bigBedItemCount(bbi));
        fprintf(output, "%-23s %llu\n", "Primary data size:",
            bbi->unzoomedIndexOffset - bbi->unzoomedDataOffset);
        struct bbiChromInfo *chrom, *chromList = bbiChromList(bbi);
        fprintf(output, "%-23s %d\n", "Zoom levels:", bbi->zoomLevels);
        fprintf(output, "%-23s %d\n", "Chrom count:", slCount(chromList));
        fprintf(output, "%-23s\n", "Chrom size:", slCount(chromList));
        for (chrom=chromList; chrom != NULL; chrom = chrom->next)
        	fprintf(output, "    %-19s %u\n", chrom->name, chrom->size);
        struct bbiSummaryElement sum = bbiTotalSummary(bbi);
        fprintf(output, "%-23s %llu\n", "Bases covered", sum.validCount);
        fprintf(output, "%-23s %f\n", "Mean depth:", sum.sumData/sum.validCount);
        fprintf(output, "%-23s %f\n", "Min depth:", sum.minVal);
        fprintf(output, "%-23s %f\n", "Max depth:", sum.maxVal);
        fprintf(output, "%-23s %f\n", "Std of depth:", calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount));
        fprintf(output, "%-23s %d\n\n", "Number of fields:", bbi->fieldCount);
        fprintf(output, "Available fields (with type VARCHAR if unspecified or all=1):\n");
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
                fprintf(output, "%-23s %s\n", buf, wrap(col->comment, strlen(buf)).c_str());
            }
        }
        bbiFileClose(&bbi);
    } else {
		fprintf(stderr, "Unknown track file type. Only local or remote tabix-indexed "
			" vcf files with extension .vcf.gz, indexed BAM files, bigWig and bigBed files are "
			"supported.");
		return;
	}
    if (!output_file.empty())
        fclose(output);
}


bool tabixFetch(const std::string & filenameOrURL, const std::vector<std::string> & regions, 
    const std::string & output_file = std::string(), bool print_header=true)
{
    FILE * output = NULL;
    if (output_file.empty()) {
        output = stdout;
    } else {
        if (print_header)
            output = fopen(output_file.c_str(), "w");
        else
            output = fopen(output_file.c_str(), "a");
        if (!output) {
            fprintf(stderr, "Failed to open output file %s.", output_file.c_str());
            return false;
        }
    }
    //
	ti_conf_t & conf = ti_conf_vcf;

    std::string fnidx = filenameOrURL + ".tbi";

    tabix_t *t;
    int is_remote = (strncmp(fnidx.c_str(), "ftp://", 6) == 0 || strncmp(fnidx.c_str(), "http://", 7) == 0) ? 1 : 0;
    if (!is_remote )
    {
        // Common source of errors: new VCF is used with an old index
        struct stat stat_tbi,stat_vcf;
        stat(fnidx.c_str(), &stat_tbi);
        stat(filenameOrURL.c_str(), &stat_vcf);
        if (stat_vcf.st_mtime > stat_tbi.st_mtime )
        {
            fprintf(stderr, "[tabix] the index file %s either does not exist or is older than the vcf file. Please reindex.\n",
                fnidx.c_str());
            return false;
        }
    }
    if ((t = ti_open(filenameOrURL.c_str(), 0)) == 0) {
        fprintf(stderr, "[main] fail to open the data file.\n");
        return false;
    }

    // retrieve from specified regions
    int len;
    const char *s;
    
    if (ti_lazy_index_load(t) < 0) {
        fprintf(stderr,"[tabix] failed to load the index file.\n");
        return false;
    }

    const ti_conf_t * idxconf = ti_get_conf(t->idx);
    if ( print_header )
    {
        // If requested, print the header lines here
        ti_iter_t iter = ti_query(t, 0, 0, 0);
        while ((s = ti_read(t, iter, &len)) != 0) {
            if ((int)(*s) != idxconf->meta_char)
                break;
            fputs(s, output);
            fputc('\n', output);
        }
        ti_iter_destroy(iter);
    }
    for (int i = 0; i < regions.size(); ++i) {
        int tid, beg, end;
        if (ti_parse_region(t->idx, regions[i].c_str(), &tid, &beg, &end) == 0) {
            ti_iter_t iter = ti_queryi(t, tid, beg, end);
                while ((s = ti_read(t, iter, &len)) != 0) {
                    fputs(s, output);
                    fputc('\n', output);
            }
            ti_iter_destroy(iter);
        } 
    }
    ti_close(t);
    if (!output_file.empty())
        fclose(output);
	return true;
}

%} 


void showTrack(const std::string & track_file, const std::string output_file = std::string());
bool tabixFetch(const std::string & filenameOrURL, const std::vector<std::string> & regions,
    const std::string & output_file = std::string(), bool print_header=true);

