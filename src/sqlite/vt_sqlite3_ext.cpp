/**
 * $File: vtools_sqlite.c $
 * $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
 * $Rev: 4234 $
 *
 * This file is part of variant_tools, a software application to annotate,
 * summarize, and filter variants for next-gen sequencing ananlysis.
 * Please visit http://varianttools.sourceforge.net for details.
 *
 * Copyright (C) 2011 Bo Peng (bpeng@mdanderson.org)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma GCC diagnostic ignored "-Wwrite-strings"

extern "C"  {
#include <sqlite3ext.h>

int sqlite3_extension_init(sqlite3 * db, char ** pzErrMsg, const sqlite3_api_routines * pApi);

}
SQLITE_EXTENSION_INIT1

/*
** The least() returns the minimal number. It is similar to min(X,Y,...) but ignores NULL values
*/
static void least_not_null_func(
                                sqlite3_context * context,
                                int argc,
                                sqlite3_value ** argv
                                )
{
	int i;
	int iBest = -1;
	double least_value = 0;
	double value = 0;

	for (i = 0; i < argc; i++) {
		// if null, ignore
		if (sqlite3_value_type(argv[i]) == SQLITE_NULL)
			continue;
		// if first non-NULL value, record
		else if (iBest == -1) {
			iBest = i;
			least_value = sqlite3_value_double(argv[i]);
			continue;
		}
		value = sqlite3_value_double(argv[i]);
		if (value < least_value) {
			iBest = i;
			least_value = value;
		}
	}
	sqlite3_result_value(context, argv[iBest == -1 ? 0 : iBest]);
}


#include "cgatools/reference/ChromosomeIdField.hpp"
#include "cgatools/reference/CompactDnaSequence.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/reference/CrrFile.hpp"
#include "cgatools/reference/CrrFileWriter.hpp"
#include "cgatools/reference/GeneDataStore.hpp"
#include "cgatools/reference/RangeAnnotationStore.hpp"
#include "cgatools/reference/RepeatMaskerStore.hpp"
#include "cgatools/reference/range.hpp"

#include <vector>
#include <utility>
#include <iomanip>
#include <sstream>
#include <exception>
#include <fstream>

typedef std::map<std::string, cgatools::reference::CrrFile *> RefGenomeFileMap;
typedef std::map<std::string, int> ChrNameMap;

RefGenomeFileMap refFileMap;
ChrNameMap chrNameMap;

static void ref_sequence(
                         sqlite3_context * context,
                         int argc,
                         sqlite3_value ** argv
                         )
{
	if (argc < 3 || argc > 4 ||
	    sqlite3_value_type(argv[0]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[1]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[2]) == SQLITE_NULL) {
		sqlite3_result_error(context, "Wrong number of parameters.", -1);
		return;
	}
	// auto passed, should be string
	std::string ref_file = std::string((char *)sqlite3_value_text(argv[0]));
	// need to check type
	if (sqlite3_value_type(argv[1]) != SQLITE_TEXT) {
		sqlite3_result_error(context, "A chromosome name is expected.", -1);
		return;
	}
	std::string chr = std::string((char *)sqlite3_value_text(argv[1]));
	if (sqlite3_value_type(argv[2]) != SQLITE_INTEGER) {
		sqlite3_result_error(context, "A 1-based position of integer type is expected.", -1);
		return;
	}
	int start = sqlite3_value_int(argv[2]);
	int end = 0;

	if (argc == 4) {
		end = sqlite3_value_int(argv[3]);

		if (start > end) {
			sqlite3_result_error(context, "incorrect chromosomal range", -1);
			return;
		}
	}

	// get reference genome file
	RefGenomeFileMap::const_iterator it = refFileMap.find(ref_file);
	cgatools::reference::CrrFile * cf = NULL;

	if (it == refFileMap.end()) {
		try {
			cf = new cgatools::reference::CrrFile(ref_file);
		} catch (cgatools::util::Exception & e) {
			sqlite3_result_error(context, e.what(), -1);
			return;
		}
		refFileMap[ref_file] = cf;
	} else {
		cf = it->second;
	}
	//
	// get chromosome index
	int chrIdx = 0;
	ChrNameMap::const_iterator cit = chrNameMap.find(chr);
	if (cit == chrNameMap.end()) {
		try {
			// cgatools' chromosome names have leading chr
			// we add chr to 1, 2, 3, ... etc but not to other
			// names. There might be a problem here.
			if (chr.size() <= 2)
				chrIdx = cf->getChromosomeId("chr" + chr);
			else
				chrIdx = cf->getChromosomeId(chr);
			chrNameMap[chr] = chrIdx;
		} catch (cgatools::util::Exception & e) {
			// unrecognized chromosome
			//sqlite3_result_error(context, e.what(), -1);
			sqlite3_result_null(context);
			return;
		}
	} else {
		chrIdx = cit->second;
	}
	//
	if (argc == 4) {
		// sequence
		try {
			std::string res = cf->getSequence(cgatools::reference::Range(chrIdx, start - 1, end));
			sqlite3_result_text(context, (char *)(res.c_str()), -1, SQLITE_TRANSIENT);
		} catch (cgatools::util::Exception & e) {
			// position out of range
			sqlite3_result_null(context);
		}
	} else {
		// single base
		try {
			char res[2];
			res[0] = cf->getBase(cgatools::reference::Location(chrIdx, start - 1));
			res[1] = '\0';
			sqlite3_result_text(context, (char *)res, -1, SQLITE_TRANSIENT);
		} catch (cgatools::util::Exception & e) {
			// position out of range
			sqlite3_result_null(context);
		}
	}
}


static void mut_sequence(
                         sqlite3_context * context,
                         int argc,
                         sqlite3_value ** argv
                         )
{
	if (argc < 7 || argc > 8 ||
		// ref file
	    sqlite3_value_type(argv[0]) == SQLITE_NULL ||
		// allele chr
	    sqlite3_value_type(argv[1]) == SQLITE_NULL ||
		// allele pos
	    sqlite3_value_type(argv[2]) == SQLITE_NULL ||
		// allele ref
	    sqlite3_value_type(argv[3]) == SQLITE_NULL ||
		// allele alt
	    sqlite3_value_type(argv[4]) == SQLITE_NULL ||
		// range chr
	    sqlite3_value_type(argv[5]) == SQLITE_NULL ||
		// range start
	    sqlite3_value_type(argv[6]) == SQLITE_NULL) {
		sqlite3_result_error(context, "Wrong number of parameters.", -1);
		return;
	}
	// auto passed, should be string
	std::string ref_file = std::string((char *)sqlite3_value_text(argv[0]));
	// chr
	if (sqlite3_value_type(argv[1]) != SQLITE_TEXT) {
		sqlite3_result_error(context, "A chromosome name is expected.", -1);
		return;
	}
	std::string a_chr = std::string((char *)sqlite3_value_text(argv[1]));
	// pos
	if (sqlite3_value_type(argv[2]) != SQLITE_INTEGER) {
		sqlite3_result_error(context, "A 1-based position of integer type is expected.", -1);
		return;
	}
	int a_pos = sqlite3_value_int(argv[2]);
	// ref
	std::string a_ref = std::string((char *)sqlite3_value_text(argv[3]));
	// alt
	std::string a_alt = std::string((char *)sqlite3_value_text(argv[4]));
	// need to check type
	if (sqlite3_value_type(argv[5]) != SQLITE_TEXT) {
		sqlite3_result_error(context, "A chromosome name is expected.", -1);
		return;
	}
	std::string chr = std::string((char *)sqlite3_value_text(argv[5]));
	if (sqlite3_value_type(argv[6]) != SQLITE_INTEGER) {
		sqlite3_result_error(context, "A 1-based position of integer type is expected.", -1);
		return;
	}
	int start = sqlite3_value_int(argv[6]);
	int end = 0;

	if (argc == 8) {
		end = sqlite3_value_int(argv[7]);

		if (start > end) {
			sqlite3_result_error(context, "incorrect chromosomal range", -1);
			return;
		}
	}

	// get reference genome file
	RefGenomeFileMap::const_iterator it = refFileMap.find(ref_file);
	cgatools::reference::CrrFile * cf = NULL;

	if (it == refFileMap.end()) {
		try {
			cf = new cgatools::reference::CrrFile(ref_file);
		} catch (cgatools::util::Exception & e) {
			sqlite3_result_error(context, e.what(), -1);
			return;
		}
		refFileMap[ref_file] = cf;
	} else {
		cf = it->second;
	}
	//
	// get chromosome index
	int chrIdx = 0;
	ChrNameMap::const_iterator cit = chrNameMap.find(chr);
	if (cit == chrNameMap.end()) {
		try {
			// cgatools' chromosome names have leading chr
			// we add chr to 1, 2, 3, ... etc but not to other
			// names. There might be a problem here.
			if (chr.size() <= 2)
				chrIdx = cf->getChromosomeId("chr" + chr);
			else
				chrIdx = cf->getChromosomeId(chr);
			chrNameMap[chr] = chrIdx;
		} catch (cgatools::util::Exception & e) {
			// unrecognized chromosome
			//sqlite3_result_error(context, e.what(), -1);
			sqlite3_result_null(context);
			return;
		}
	} else {
		chrIdx = cit->second;
	}
	//
	std::string res("?");
	try {
		if (argc == 8) 
			res = cf->getSequence(cgatools::reference::Range(chrIdx, start - 1, end));
		else
			res[0] = cf->getBase(cgatools::reference::Location(chrIdx, start - 1));
	} catch (cgatools::util::Exception & e) {
		// position out of range
		sqlite3_result_null(context);
	}
	// if on the same chromosome, and position overlap 
	if (chr == a_chr) {
		// insertion
		if (a_ref[0] == '-') {
			// insertion only appears if the location is included in the range
			if (a_pos >= start && a_pos < start + res.size())
				res.insert(a_pos - start, "\033[32m" + a_alt + "\033[0m");
		} else {
			// case 1
			//        [xxxxxxxxxx----xxxxxxxxx]
			// case 2
			//    ----[-----xxxxxxx]
			// case 3
			//        [xxxxxxxxxxxxxxx] -----
			// case 4
			//   -----[xxxxxxxxxxx]
			//
			int d_start = std::max(start, a_pos) - start;
			int d_end = std::min(start + res.size(), a_pos + a_ref.size()) - start;
			if (d_start < res.size()) {
				if (a_alt[0] == '-')
					for (size_t i = d_start; i < d_end; ++i)
						res[i] = '-';				
				else
					for (size_t i = 0; i < d_end - d_start; ++i)
						res[d_start + i] = a_alt[i];
				res.insert(d_end, "\033[0m");
				res.insert(d_start, "\033[32m");
			}
		}
	}
	sqlite3_result_text(context, (char *)(res.c_str()), -1, SQLITE_TRANSIENT);
}



static void vcf_variant(
                        sqlite3_context * context,
                        int argc,
                        sqlite3_value ** argv
                        )
{
	if (argc < 5 ||
	    sqlite3_value_type(argv[0]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[1]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[2]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[3]) == SQLITE_NULL ||
	    sqlite3_value_type(argv[4]) == SQLITE_NULL) {
		sqlite3_result_null(context);
		return;
	}

	std::string ref_file = std::string((char *)sqlite3_value_text(argv[0]));
	std::string chr = std::string((char *)sqlite3_value_text(argv[1]));
	int pos = sqlite3_value_int(argv[2]);
	std::string ref = std::string((char *)sqlite3_value_text(argv[3]));
	std::string alt = std::string((char *)sqlite3_value_text(argv[4]));
	std::string name;
	if (argc == 6) {
		if (sqlite3_value_type(argv[5]) == SQLITE_NULL)
			name = ".";
		else
			name = std::string((char *)sqlite3_value_text(argv[5]));
	}
	//
	std::stringstream res;
	if (ref != "-" && alt != "-") {
		res << chr << '\t' << pos << '\t';
		if (argc == 6)
			res << name << '\t';
		res << ref << '\t' << alt;
	} else {
		// get reference genome file
		cgatools::reference::CrrFile * cf = NULL;
		RefGenomeFileMap::const_iterator it = refFileMap.find(ref_file);
		if (it == refFileMap.end()) {
			try {
				cf = new cgatools::reference::CrrFile(ref_file);
			} catch (cgatools::util::Exception & e) {
				sqlite3_result_error(context, e.what(), -1);
				return;
			}

			refFileMap[ref_file] = cf;
		} else {
			cf = it->second;
		}
		//
		// get chromosome index
		int chrIdx = 0;
		int succ = true;
		ChrNameMap::const_iterator cit = chrNameMap.find(chr);
		if (cit == chrNameMap.end()) {
			try {
				// cgatools' chromosome names have leading chr
				// we add chr to 1, 2, 3, ... etc but not to other
				// names. There might be a problem here.
				if (chr.size() <= 2)
					chrIdx = cf->getChromosomeId("chr" + chr);
				else
					chrIdx = cf->getChromosomeId(chr);
				chrNameMap[chr] = chrIdx;
			} catch (cgatools::util::Exception & e) {
				// unrecognized chromosome
				//sqlite3_result_error(context, e.what(), -1);
				succ = false;
			}
		} else {
			chrIdx = cit->second;
		}

		if (succ) {
			try {
				char pad = cf->getBase(cgatools::reference::Location(chrIdx, pos - 2));

				// 10 - A
				// 9  G GA
				res << chr << '\t' << (pos - 1) << '\t';
				if (argc == 6)
					res << name << '\t';
				if (ref == "-")
					res << pad << '\t' << pad << alt;
				else
					// 10 A -
					// 9 GA G
					res << pad << ref << '\t' << pad;
			} catch (cgatools::util::Exception & e) {
				// position out of range
				succ = false;
			}
		}
		if (!succ) {
			res << chr << '\t' << pos << '\t';
			if (argc == 6)
				res << name << '\t';
			res << ref << '\t' << alt;
		}
	}
	sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
}


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


struct FieldInfo
{
	bool all;
	int column;
	std::string name;
	std::string query;
};

struct TrackInfo;
typedef void (* track_handler)(void *, char *, int, char *, char *, FieldInfo *, TrackInfo *, sqlite3_context *);

struct TrackInfo
{
	//
	int file_type;
	// pointer to opened file
	void * file;
	// pointer to index file, if any
	void * index_file;
	// function to handle stuff
	track_handler handler;
	// default return column for the file type
	int default_col;
	// name map that map "chrom" etc to column index
	std::map<std::string, int> name_map;
	// if the chromosome has leading "chr"
	bool with_leading_chr;
	//
	cgatools::reference::CrrFile * reference;
	int chrIdx;

	TrackInfo() : file_type(0), file(NULL), index_file(NULL),
		handler(NULL), default_col(0), name_map(), with_leading_chr(false),
		reference(NULL), chrIdx(0)
	{
	}


	//
	void print()
	{
		fprintf(stderr, "Default col: %d\n", default_col);
		fprintf(stderr, "Leading chr: %d\n", with_leading_chr);
		std::map<std::string, int>::iterator it = name_map.begin();
		std::map<std::string, int>::iterator it_end = name_map.end();
		for (; it != it_end; ++it)
			fprintf(stderr, "%s : %d \n", it->first.c_str(), it->second);
	}


};

typedef std::map<std::string, struct TrackInfo> TrackFileMap;
TrackFileMap trackFileMap;

#define BAM_FILE 0
#define BIGBED_FILE 1
#define BIGWIG_FILE 2
#define VCFTABIX_FILE 3

void bigBedTrack(void * track_file, char * chr, int pos, char *, char *, FieldInfo * fi,
                 TrackInfo * info, sqlite3_context * context)
{
	bbiFile * cf = (bbiFile *)track_file;
	struct lm * bbLm = lmInit(0);
	//
	char * chrName = NULL;

	if (info->with_leading_chr) {
		// can do better here (allocate len(chr) + 3
		char buf[255];
		strcpy(buf, "chr");
		strcat(buf, chr);
		chrName = buf;
	} else
		chrName = chr;
	//
	//
	struct bigBedInterval * ivList = bigBedIntervalQuery(cf, chrName, pos - 1, pos, 0, bbLm);
	if (ivList == NULL) {
		sqlite3_result_null(context);
	} else {
		// if returnning typed-valued, only the first record will be considered
		if (fi->column > cf->fieldCount)
			sqlite3_result_null(context);
		else if (!fi->all) {
			if (fi->column == 0)
				// return matched or not
				sqlite3_result_int(context, 1);
			else if (fi->column == 1)
				sqlite3_result_text(context, chrName, -1, SQLITE_TRANSIENT);
			else if (fi->column == 2)
				sqlite3_result_int64(context, ivList->start);
			else if (fi->column == 3)
				sqlite3_result_int64(context, ivList->end);
			else if (fi->column > 3 && fi->column <= cf->fieldCount) {
				int res_type = 0; // varchar
				size_t i = 1;
				struct asColumn * col;
				struct asObject * as = bigBedAs(cf);
				for (col = as->columnList; col != NULL; col = col->next, ++i) {
					struct asTypeInfo * ltype = col->lowType;
					if (i == fi->column) {
						if (asTypesIsInt(ltype->type))
							res_type = 1;   // INTEGER
						else if (asTypesIsFloating(ltype->type))
							res_type = 2;   // FLOAT
						break;
					}
				}

				char * rest = ivList->rest;
				// skip a few \t
				i = 3;
				char * pch = strtok(rest, "\t");
				while (++i < fi->column)
					pch = strtok(NULL, "\t");
				// text fields
				if (res_type == 0)
					sqlite3_result_text(context, pch, -1, SQLITE_TRANSIENT);
				else if (res_type == 1)
					sqlite3_result_int64(context, atoi(pch));
				else
					sqlite3_result_double(context, atof(pch));
			} else
				sqlite3_result_null(context);
		} else {
			std::stringstream res;
			bool first = true;
			struct bigBedInterval * iv;
			for (iv = ivList; iv != NULL; iv = iv->next) {
				if (!first)
					res << '|';
				else
					first = false;
				//
				if (fi->column == 1)
					res << chrName;
				else if (fi->column == 2)
					res << iv->start;
				else if (fi->column == 3)
					res << iv->end;
				else if (fi->column > 3 && fi->column <= cf->fieldCount) {
					char * rest = ivList->rest;
					// skip a few \t
					int i = 3;
					char * pch = strtok(rest, "\t");
					while (++i < fi->column)
						pch = strtok(NULL, "\t");
					res << pch;
				}
			}
			sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
		}
	}
	lmCleanup(&bbLm);
}


static void bigWigTrack(void * track_file, char * chr, int pos, char *, char *, FieldInfo * fi,
                        TrackInfo * info, sqlite3_context * context)
{
	bbiFile * cf = (bbiFile *)track_file;
	struct lm * bbLm = lmInit(0);
	//
	char * chrName = NULL;

	if (info->with_leading_chr) {
		// can do better here (allocate len(chr) + 3
		char buf[255];
		strcpy(buf, "chr");
		strcat(buf, chr);
		chrName = buf;
	} else
		chrName = chr;

	struct bbiInterval * ivList = bigWigIntervalQuery(cf, chrName, pos - 1, pos, bbLm);
	if (ivList == NULL) {
		sqlite3_result_null(context);
	} else {
		if (fi->column > 4)
			sqlite3_result_null(context);
		else if (!fi->all) {
			if (fi->column == 1)
				sqlite3_result_text(context, chrName, -1, SQLITE_TRANSIENT);
			else if (fi->column == 2)
				sqlite3_result_int64(context, ivList->start);
			else if (fi->column == 3)
				sqlite3_result_int64(context, ivList->end);
			else if (fi->column == 4)
				sqlite3_result_double(context, ivList->val);
		} else {
			std::stringstream res;
			struct bbiInterval * iv;
			bool first = true;
			for (iv = ivList; iv != NULL; iv = iv->next) {
				if (!first)
					res << "|";
				else
					first = false;
				//
				if (fi->column == 1)
					res << chrName;
				else if (fi->column == 2)
					res << iv->start;
				else if (fi->column == 3)
					res << iv->end;
				else if (fi->column == 4)
					res << iv->val;
			}
			sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
		}
	}
	lmCleanup(&bbLm);
}


bool vcf_match(char * ref, char * alt, int count, char ** alleles)
{
	if (*ref == '-') {
		// insertion
		int leading = strlen(alleles[0]);
		for (size_t i = 1; i < count; ++i)
			if (strlen(alleles[i]) > leading && strcmp(alt, alleles[i] + leading) == 0)
				return true;
		return false;
	} else if (*alt == '-') {
		int reflen = strlen(alleles[0]);
		// deletion
		for (size_t i = 1; i < count; ++i)
			if (strlen(alleles[i]) < reflen && strcmp(ref, alleles[0] + strlen(alleles[i])) == 0)
				return true;
		return false;
	} else if (strlen(ref) == 1 && strlen(alt) == 1) {
		// regular SNVs
		if (strcmp(ref, alleles[0]) != 0)
			return false;
		for (size_t i = 1; i < count; ++i)
			if (strcmp(alt, alleles[i]) == 0)
				return true;
		return false;
	} else {
		// complex variants
		for (size_t i = 1; i < count; ++i) {
			// find common leading
			char * p = alleles[0];
			char * q = alleles[i];
			for (; *p != '\0' && *q != '\0' && *p == *q; ++p, ++q) ;
			if (strcmp(ref, p) == 0 && strcmp(alt, q) == 0)
				return true;
		}
		return false;
	}
}


static void vcfTabixTrack(void * track_file, char * chr, int pos, char * ref, char * alt, FieldInfo * fi,
                          TrackInfo * info, sqlite3_context * context)
{
	struct vcfFile * cf = (struct vcfFile *)track_file;

	char * chrName = NULL;

	if (info->with_leading_chr) {
		// can do better here (allocate len(chr) + 3
		char buf[255];
		strcpy(buf, "chr");
		strcat(buf, chr);
		chrName = buf;
	} else
		chrName = chr;

	// clear record pool
	vcfFileFlushRecords(cf);
	// read records
	//
	// adjust pos for indels
	int nRecord = 0;
	if (*ref == '-' || *alt == '-') {
		// for insertion, pos, -, aaa
		// it should be   pos - 1 , X, Xaaa
		// in vcf file
		nRecord = vcfTabixBatchRead(cf, chrName, pos - 2, pos - 1, VCF_IGNORE_ERRS, -1);
	} else {
		nRecord = vcfTabixBatchRead(cf, chrName, pos - 1, pos, VCF_IGNORE_ERRS, -1);
	}
	// first case, no result
	if (nRecord == 0) {
		sqlite3_result_null(context);
	} else if (!fi->all) {
		struct vcfRecord * rec = NULL;
		bool match = false;
		for (rec = cf->records; rec != NULL; rec = rec->next) {
			if (vcf_match(ref, alt, rec->alleleCount, rec->alleles)) {
				match = true;
				break;
			}
		}
		if (!match) {
			sqlite3_result_null(context);
			return;
		}

		if (fi->column == 1)
			// chrom
			sqlite3_result_text(context, chrName, -1, SQLITE_TRANSIENT);
		else if (fi->column == 2)
			// pos
			sqlite3_result_int(context, rec->chromStart);
		else if (fi->column == 3)
			// name
			sqlite3_result_text(context, rec->name, -1, SQLITE_TRANSIENT);
		else if (fi->column == 4)
			// ref
			sqlite3_result_text(context, rec->alleles[0], -1, SQLITE_TRANSIENT);
		else if (fi->column == 5) {
			// alt
			char alleles[255];
			alleles[0] = '\0';
			for (size_t i = 1 ; i < rec->alleleCount; ++i) {
				if (i > 1)
					strcat(alleles, ",");
				strcat(alleles, rec->alleles[i]);
			}
			sqlite3_result_text(context, alleles, -1, SQLITE_TRANSIENT);
		} else if (fi->column == 6)
			// qual
			sqlite3_result_text(context, rec->qual, -1, SQLITE_TRANSIENT);
		else if (fi->column == 7) {
			// filter
			char filter[255];
			filter[0] = '\0';
			for (size_t i = 0; i < rec->filterCount; ++i) {
				if (i > 0)
					strcat(filter, ",");
				strcat(filter, rec->filters[i]);
			}
			sqlite3_result_text(context, filter, -1, SQLITE_TRANSIENT);
		} else if (fi->column == 8) {
			// info
			const char * p = fi->name.c_str();
			if (p != NULL)
				while (*p != '\0' && *p != '.') ++p;
			if (p == NULL || *p == '\0') // if no dot
				sqlite3_result_text(context, rec->unparsedInfoElements, -1, SQLITE_TRANSIENT);
			else {
				// p points to the name of info
				++p;
				vcfInfoDef * el = vcfInfoDefForKey(cf, p);
				bool isFlag = el->type == vcfInfoFlag;

				int l = strlen(p);
				// scan INFO and find index of p
				char * info = rec->unparsedInfoElements;
				bool found = false;
				char * q = info;
				for (; *q != '\0'; ) {
					// find it
					if (strncmp(q, p, l) == 0 && (*(q + l) == ';' || *(q + l) == '\0' || *(q + l) == '=')) {
						found = true;
						if (!isFlag) {
							q += l + 1;
							char * r = q + 1;
							while (*r != ';' && *r != '\0') ++r;
							*r = '\0';
						}
						break;
					} else {
						while (*q != ';' && *q != '\0') ++q;
						if (*q != '\0')
							++q;
					}
				}
				if (isFlag)
					sqlite3_result_int64(context, found ? 1 : 0);
				else if (el->type == vcfInfoInteger)
					sqlite3_result_int64(context, atoi(q));
				else if (el->type == vcfInfoFloat)
					sqlite3_result_int64(context, atof(q));
				else
					sqlite3_result_text(context, q, -1, SQLITE_TRANSIENT);
			}
		} else if (fi->column == 9)
			sqlite3_result_text(context, rec->format, -1, SQLITE_TRANSIENT);
		else if (fi->column >= 10 && fi->column < 10 + cf->genotypeCount) {
			// if there is no '.' in the name
			const char * p = fi->name.c_str();
			if (p != NULL)
				while (*p != '\0' && *p != '.') ++p;
			if (p == NULL || *p == '\0') // if no dot
				sqlite3_result_text(context, rec->genotypeUnparsedStrings[fi->column - 10], -1, SQLITE_TRANSIENT);
			else {
				// p point to the name of FMT
				++p;
				vcfInfoDef * el = vcfInfoDefForGtKey(cf, p);
				// length of the format string
				int l = strlen(p);
				// scan FORMAT and find index of p
				char * fmt = rec->format;
				int idx = 0;
				bool found = false;
				for (char * q = fmt; *q != '\0'; ) {
					// find it
					if (strncmp(q, p, l) == 0 && (*(q + l) == '\0' || *(q + l) == ':')) {
						found = true;
						break;
					} else {
						idx += 1;
						while (*q != ':' && *q != '\0') ++q;
						if (*q != '\0')
							++q;
					}
				}
				if (found) {
					// find the idx-th piece in genotype string
					char * geno = rec->genotypeUnparsedStrings[fi->column - 10];
					char * p = geno;
					for (size_t i = 0; i < idx; ++i, ++p)
						while (*p != ':' && *p != '\0')
							++p;
					char * q = p;
					while (*q != ':' && *q != '\0') ++q;
					*q = '\0';
					if (el->type == vcfInfoFlag)
						sqlite3_result_int64(context, p != q ? 1 : 0);
					else if (el->type == vcfInfoInteger)
						sqlite3_result_int64(context, atoi(p));
					else if (el->type == vcfInfoFloat)
						sqlite3_result_int64(context, atof(p));
					else
						sqlite3_result_text(context, p, -1, SQLITE_TRANSIENT);
				}
			}
		} else
			sqlite3_result_null(context);
	} else {
		// third case
		std::stringstream res;
		struct vcfRecord * rec = NULL;
		bool first = true;
		for (rec = cf->records; rec != NULL; rec = rec->next) {
			if (!vcf_match(ref, alt, rec->alleleCount, rec->alleles))
				continue;
			if (!first)
				res << '|';
			else
				first = false;
			if (fi->column == 1)
				// chrom
				res << chrName;
			else if (fi->column == 2)
				// pos
				res << rec->chromStart;
			else if (fi->column == 3)
				// name
				res << rec->name;
			else if (fi->column == 4)
				// ref
				res << rec->alleles[0];
			else if (fi->column == 5) {
				// alt
				for (size_t i = 1 ; i < rec->alleleCount; ++i) {
					if (i > 1)
						res << ',';
					res << rec->alleles[i];
				}
			} else if (fi->column == 6)
				// qual
				res << rec->qual;
			else if (fi->column == 7) {
				// filter
				for (size_t i = 0; i < rec->filterCount; ++i) {
					if (i > 0)
						res << ',';
					res << rec->filters[i];
				}
			} else if (fi->column == 8) {
				const char * p = fi->name.c_str();
				if (p != NULL)
					while (*p != '\0' && *p != '.') ++p;
				if (p == NULL || *p == '\0') // if no dot
					res << rec->unparsedInfoElements;
				else {
					// p point to the name of info
					++p;
					vcfInfoDef * el = vcfInfoDefForKey(cf, p);
					// we should have definiton for all keys
					bool isFlag = el->type == vcfInfoFlag;

					int l = strlen(p);
					// scan INFO and find index of p
					char * info = rec->unparsedInfoElements;
					bool found = false;
					for (char * q = info; *q != '\0'; ) {
						// find it
						if (strncmp(q, p, l) == 0 && (*(q + l) == ';' || *(q + l) == '\0' || *(q + l) == '=')) {
							found = true;
							if (!isFlag) {
								q += l + 1;
								char * r = q + 1;
								while (*r != ';' && *r != '\0') ++r;
								*r = '\0';
								res << q;
							}
							break;
						} else {
							while (*q != ';' && *q != '\0') ++q;
							if (*q != '\0')
								++q;
						}
					}
					if (isFlag)
						res << (found ? 1 : 0);
				}
			} else if (fi->column == 9)
				res << rec->format;
			else if (fi->column >= 10 && fi->column < 10 + cf->genotypeCount) {
				// if there is no '.' in the name
				const char * p = fi->name.c_str();
				if (p != NULL)
					while (*p != '\0' && *p != '.') ++p;
				if (p == NULL || *p == '\0') // if no dot
					res << rec->genotypeUnparsedStrings[fi->column - 10];
				else {
					// p point to the name of FMT
					++p;
					// length of the format string
					int l = strlen(p);
					// scan FORMAT and find index of p
					char * fmt = rec->format;
					int idx = 0;
					bool found = false;
					for (char * q = fmt; *q != '\0'; ) {
						// find it
						if (strncmp(q, p, l) == 0 && (*(q + l) == '\0' || *(q + l) == ':')) {
							found = true;
							break;
						} else {
							idx += 1;
							while (*q != ':' && *q != '\0') ++q;
							if (*q != '\0')
								++q;
						}
					}
					if (found) {
						// find the idx-th piece in genotype string
						char * geno = rec->genotypeUnparsedStrings[fi->column - 10];
						char * p = geno;
						for (size_t i = 0; i < idx; ++i, ++p)
							while (*p != ':' && *p != '\0')
								++p;
						char * q = p;
						while (*q != ':' && *q != '\0') ++q;
						*q = '\0';
						res << p;
					}
				}
			}
		}
		sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
	}
}


#define OP_EQ '='
#define OP_NE '+'
#define OP_GT '>'
#define OP_GE '.'
#define OP_LT '<'
#define OP_LE ','

/* There are fancier and perhaps faster methods such as functors, but I am too tired to
 * figure out the details.
 */

// special case to speed up string comparison
bool compare_val(char type, const char * a, const char * b)
{
	if (type == OP_EQ)
		return strcmp(a, b) == 0;
	else if (type == OP_NE)
		return strcmp(a, b) != 0;
	else if (type == OP_LT)
		return strcmp(a, b) < 0;
	else if (type == OP_LE)
		return strcmp(a, b) <= 0;
	else if (type == OP_GT)
		return strcmp(a, b) > 0;
	else if (type == OP_GE)
		return strcmp(a, b) >= 0;
}


template<typename T>
bool compare_val(char type, T a, T b)
{
	if (type == OP_EQ)
		return a == b;
	else if (type == OP_NE)
		return a != b;
	else if (type == OP_LT)
		return a < b;
	else if (type == OP_LE)
		return a <= b;
	else if (type == OP_GT)
		return a > b;
	else if (type == OP_GE)
		return a >= b;
}


struct BAM_stat
{
	//
	size_t start;   // zero based start position
	size_t counter;
	char call_content[3];
	std::stringstream calls;
	std::vector<int> qual;
	std::vector<int> map_qual;
	// condition
	int shift;
	int width;
	int limit;
	int strand;
	int min_qual;
	int min_mapq;
	char delimiter[3];
	bool show_seq;
	int match_type;
	std::string reference;

	// tags used in conditions, at most 12 tags are allowed, this should be more than enough
	char cond_tag_keys[24];
	// =, ==, !=, >, >=, <= etc
	char cond_tag_op[12];
	std::vector<std::string> cond_tag_values;

	// if use color in output
	bool colorize;

	BAM_stat(size_t pos) :
		start(pos), counter(0), calls(), qual(), map_qual(),
		shift(0), width(1), limit(-1), strand(-1), min_qual(0), min_mapq(),
		show_seq(false), match_type(-1), reference(""), cond_tag_values(), colorize(false)
	{
		call_content[0] = '\0';
		call_content[1] = '\0';
		call_content[2] = '\0';
		delimiter[0] = '|';
		delimiter[1] = '\0';
		delimiter[2] = '\0';
	}
};


char conv_table[] = "=ACMGRSVTWYHKDBN";

static int fetch_func(const bam1_t * b, void * data)
{
	BAM_stat * buf = (BAM_stat *)data;

	if (b->core.n_cigar == 0 || b->core.qual < buf->min_mapq)
		return 0;

	// if on not on the same strand as specified, ignore
	if (buf->strand >= 0 && (buf->strand != (b->core.flag & 0x10)))
		return 0;

	// compare tags
	int ntags = buf->cond_tag_values.size();
	if (ntags > 0) {
		uint8_t * s = bam1_aux(b);
		while (s < b->data + b->data_len) {
			// if the key is not found, continue
			int match = -1;
			for (size_t i = 0; i < ntags + ntags; i += 2) {
				if (buf->cond_tag_keys[i] == s[0] && buf->cond_tag_keys[i + 1] == s[1]) {
					match = i;
					break;
				}
			}
			uint8_t type;
			// other wise,
			s += 2;
			type = *s;
			++s;
			if (type == 'A') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *s, uint8_t(buf->cond_tag_values[match][0])))
					return 0;
				++s;
			} else if (type == 'C') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *s, uint8_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				++s;
			} else if (type == 'c') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(int8_t *)s, int8_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				++s;
			} else if (type == 'S') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(uint16_t *)s, uint16_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 2;
			} else if (type == 's') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(int16_t *)s, int16_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 2;
			} else if (type == 'I') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(uint32_t *)s, uint32_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 4;
			} else if (type == 'i') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(int32_t *)s, int32_t(atoi(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 4;
			} else if (type == 'f') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(float *)s, float(atof(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 4;
			} else if (type == 'd') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], *(double *)s, double(atof(buf->cond_tag_values[match].c_str()))))
					return 0;
				s += 8;
			} else if (type == 'Z' || type == 'H') {
				if (match >= 0 && !compare_val(buf->cond_tag_op[match], (const char *)s, buf->cond_tag_values[match].c_str()))
					return 0;
				s += strlen((char *)s) + 1;
			} else if (type == 'B') {
				//
				if (match >= 0) {
					// sqlite3_result_error(context, "Condition involves array tags (B) is not supported.", -1);
					return 0;
				}
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
				int nItem = (int)*s;
				s += 4 + nItem * sz;
			}
		}
	}

	// record reads into the string
	bool record_reads = buf->call_content[0] == '*' || buf->call_content[0] == '&';
	std::string reads;

	uint32_t outputstart = buf->start + buf->shift;
	uint32_t outputend = outputstart + buf->width;

	// pos is the supposed position
	uint32_t pos = b->core.pos;
	// start position is affected by soft clip
	uint32_t k = 0;
	uint32_t * cigar_p = bam1_cigar(b);

	// qpos is the location into the read itself
	uint32_t qpos = 0;
	// k is the index to cigar string
	// k < b->core.n_cigar
	//
	k = 0;
	const char * ref_p = buf->reference.empty() ? NULL : buf->reference.c_str();
	// output pos is the window within which we would like to output
	int qual;
	int type_matched = -1;
	if (pos > outputstart)
		for (int i = outputstart; i < pos; ++i) {
			if (record_reads)
				reads += ' ';
			if (ref_p)
				++ref_p;
		}
	while (pos < outputend) {
		// if |......ACGT|GAAA
		if (qpos >= b->core.l_qseq) {
			if (record_reads)
				for (int i = pos; i < outputend; ++i)
					reads += ' ';
			break;
		}

		int op = cigar_p[k] & BAM_CIGAR_MASK;
		int l = cigar_p[k] >> BAM_CIGAR_SHIFT;
		if (op == BAM_CMATCH) {
			for (int i = 0; i < l && pos < outputend; ++i) {
				// right at the location, get quality score
				if (pos == buf->start) {
					qual = bam1_qual(b)[qpos];
					if (qual < buf->min_qual)
						return 0;
				}
				if (pos >= outputstart) {
					char allele = conv_table[bam1_seqi(bam1_seq(b), qpos)];
					if (pos == buf->start) {
						if (ref_p) {
							// insertion will override SNP
							if (type_matched != 2) {
								if (allele == *ref_p)
									type_matched = 0;
								else
									type_matched = 1;
							}
						} else {
							// somehow no reference genome is found (extra chromsome?)
							type_matched = 1;
						}
					}
					if (!buf->show_seq && ref_p && allele == *ref_p)
						allele = '.';
					if (record_reads) {
						if (buf->width > 1 && pos == buf->start && buf->colorize)
							reads += std::string("\033[94m") + allele + "\033[0m";
						else
							reads += allele;
					}
					if (ref_p)
						++ref_p;
				}
				++qpos;
				++pos;
			}
		} else if (op == BAM_CINS) {
			// right at the location, get quality score
			if (pos == buf->start) {
				// insertion will show matched single nucleotide
				type_matched = 2;
				qual = bam1_qual(b)[qpos];
				if (qual < buf->min_qual)
					return 0;
			}

			// for insertion, qpos will move l. pos does not move.
			if (record_reads && pos >= outputstart) {
				if (buf->colorize)
					reads += "\033[32m";
				for (int i = 0; i < l; ++i)
					reads += conv_table[bam1_seqi(bam1_seq(b), qpos + i)] ;
				if (buf->colorize)
					reads += "\033[0m";
			}
			qpos += l;
		} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
			// for deletion, qpos will not move, pos move
			for (int i = 0; i < l && pos < outputend; ++i, ++pos) {
				if (pos == buf->start) {
					type_matched = 3;
					qual = bam1_qual(b)[qpos];
					if (qual < buf->min_qual)
						return 0;
				}
				if (record_reads && pos >= outputstart) {
					if (pos == buf->start && buf->colorize)
						reads += "\033[94m*\033[0m";
					else
						reads += '*';
				}
			}
		} else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
			for (int i = 0; i < l && pos < outputend; ++i) {
				// right at the location, get quality score
				if (pos == buf->start) {
					qual = bam1_qual(b)[qpos];
					if (qual < buf->min_qual)
						return 0;
				}
				// soft clip are no seq, hard clip are not
				if (op == BAM_CSOFT_CLIP)
					++qpos;
			}
		}
		// go to the next cigar
		++k;
	}
	if (type_matched == -1) {
		fprintf(stderr, "Failed to fund the type of read. \n");
		return 0;
	}
	if (buf->match_type >= 0 && buf->match_type != type_matched) {
		return 0;
	}

	buf->map_qual.push_back(b->core.qual);
	buf->qual.push_back(qual);
	buf->counter += 1;
	if (buf->limit > 0 && buf->counter > buf->limit)
		return 0;
	if (buf->call_content[0] == '*') {
		if (!buf->calls.str().empty())
			buf->calls << buf->delimiter;
		buf->calls << reads;
	} else if (buf->call_content[0] == '&') {
		// no delimieter
		if (reads.size() == 1 || (reads[0] == '\033' && reads.size() == 10))
			buf->calls << reads;
		else if (buf->colorize)
			buf->calls << "\033[32mI\033[0m";
		else
			buf->calls << 'I';
	} else if (buf->call_content[0] == '%') {
		if (!buf->calls.str().empty())
			buf->calls << buf->delimiter;
		buf->calls << std::hex << b->core.flag;
	} else if (buf->call_content[0] != '\0') {
		if (!buf->calls.str().empty())
			buf->calls << buf->delimiter;
		uint8_t * s = bam1_aux(b);
		while (s < b->data + b->data_len) {
			// if the key is not found, continue
			bool match = false;
			if (buf->call_content[0] == s[0] && buf->call_content[1] == s[1])
				match = true;

			// other wise,
			s += 2;
			uint8_t type = *s;
			++s;
			if (type == 'A') {
				if (match)
					buf->calls << *(char *)s;
				++s;
			} else if (type == 'C') {
				if (match)
					buf->calls << (int)(*s);
				++s;
			} else if (type == 'c') {
				if (match)
					buf->calls << *(int8_t *)s;
				++s;
			} else if (type == 'S') {
				if (match)
					buf->calls << *(uint16_t *)s;
				s += 2;
			} else if (type == 's') {
				if (match)
					buf->calls << *(int16_t *)s;
				s += 2;
			} else if (type == 'I') {
				if (match)
					buf->calls << *(uint32_t *)s;
				s += 4;
			} else if (type == 'i') {
				if (match)
					buf->calls << *(int32_t *)s;
				s += 4;
			} else if (type == 'f') {
				if (match)
					buf->calls << *(float *)s;
				s += 4;
			} else if (type == 'd') {
				if (match)
					buf->calls << *(double *)s;
				s += 8;
			} else if (type == 'Z' || type == 'H') {
				if (match)
					buf->calls << std::string((char *)s);
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
				int nItem = (int)*s;
				s += 4;
				for (size_t k = 0; k < nItem; ++k) {
					if (k > 0 && match)
						buf->calls << ",";
					if (subtype == 'C') {
						if (match)
							buf->calls << (int)(*s);
						++s;
					} else if (subtype == 'c') {
						if (match)
							buf->calls << *(int8_t *)s;
						++s;
					} else if (subtype == 'S') {
						if (match)
							buf->calls << *(uint16_t *)s;
						s += 2;
					} else if (subtype == 's') {
						if (match)
							buf->calls << *(int16_t *)s;
						s += 2;
					} else if (subtype == 'I') {
						if (match)
							buf->calls << *(uint32_t *)s;
						s += 4;
					} else if (subtype == 'i') {
						if (match)
							buf->calls << *(int32_t *)s;
						s += 4;
					} else if (subtype == 'f') {
						if (match)
							buf->calls << *(float *)s;
						s += 4;
					}
				}
			}
		}
	}
	return 0;
}


// calculate depth
static void bamTrack(void * track_file, char * chr, int pos, char *, char *, FieldInfo * fi,
                     TrackInfo * info, sqlite3_context * context)
{
	samfile_t * sf = (samfile_t *)track_file;
	bam_index_t * idx = (bam_index_t *)info->index_file;
	//
	int tid = 0;

	//
	// we can cache this piece of code in the info structure (as chrom map)
	if (info->with_leading_chr) {
		char chrName[255];
		strcpy(chrName, "chr");
		strcat(chrName, chr);
		tid = bam_get_tid(sf->header, chrName);
	} else {
		tid = bam_get_tid(sf->header, chr);
	}
	BAM_stat buf(pos - 1);
	// set parameters
	if (fi->name == "reads") {
		buf.width = 5;
		buf.call_content[0] = '*';
	} else if (fi->name == "calls") {
		buf.width = 1;
		buf.call_content[0] = '&';
	} else if (fi->name == "flag")
		buf.call_content[0] = '%';
	else if (fi->name.size() == 2) {  // a tag?
		buf.call_content[0] = fi->name[0];
		buf.call_content[1] = fi->name[1];
	}
	if (!fi->query.empty()) {
		char * query = strdup(fi->query.c_str());
		char * pch = strtok(query, "&");
		while (pch != NULL) {
			if (strncmp(pch, "min_qual=", 9) == 0)
				buf.min_qual = atoi(pch + 9);
			else if (strncmp(pch, "min_mapq=", 9) == 0)
				buf.min_mapq = atoi(pch + 9);
			else if (strncmp(pch, "show_seq=", 9) == 0)
				buf.show_seq = *(pch + 9) == '1' || *(pch + 9) == 'y' || *(pch + 9) == 'Y';
			else if (strncmp(pch, "type=", 5) == 0) {
				buf.match_type = atoi(pch + 5);
				if (buf.match_type != 0 && buf.match_type != 1 && buf.match_type != 2 && buf.match_type != 3) {
					sqlite3_result_error(context, "Incorrect type of match is specified (0 (match), 1 (mismatch), "
						"2 (insertion) and 3 (deletion) are acceptable)", -1);
					return;
				}
					
			} else if (strncmp(pch, "strand=", 7) == 0) {
				if (pch[7] == '+')
					buf.strand = 0;
				else if (pch[7] == '-')
					buf.strand = 0x10;
				else {
					sqlite3_result_error(context, "Incorrect strand is specified (+ or - are acceptable)", -1);
					return;
				}
			} else if (strncmp(pch, "width=", 6) == 0 && fi->name == "reads") {
				buf.width = atoi(pch + 6);
				if (buf.width < 0) {
					sqlite3_result_error(context, "Width of reads must be positive", -1);
					return;
				}
			} else if (strncmp(pch, "delimiter=", 10) == 0) {
				strncpy(buf.delimiter, pch + 10, 2);
				if (strcmp(buf.delimiter, "\\t") == 0)
					strcpy(buf.delimiter, "\t");
				else if (strcmp(buf.delimiter, "\\n") == 0)
					strcpy(buf.delimiter, "\n");
			} else if (strncmp(pch, "start=", 6) == 0 && fi->name == "reads") {
				buf.shift = atoi(pch + 6);
				if (buf.shift > 0) {
					sqlite3_result_error(context, "Output of reads must cover variant location (start < 0)", -1);
					return;
				}
				// at last shit wide
				if (buf.width + buf.shift <= 0)
					buf.width = - buf.shift + 1;
				// this is a tag
			} else if (strncmp(pch, "limit=", 6) == 0) {
				buf.limit = atoi(pch + 6);
			} else if (strncmp(pch, "color=", 6) == 0) {
				buf.colorize = atoi(pch + 6);
			} else if (strlen(pch) > 3 && (pch[2] == '=' || pch[2] == '!' || pch[2] == '>' || pch[2] == '<')) {
				size_t n = buf.cond_tag_values.size();
				if (n < 12) {
					buf.cond_tag_keys[n + n] = pch[0];
					buf.cond_tag_keys[n + n + 1] = pch[1];
					if (pch[2] == '=') {
						buf.cond_tag_op[n] = OP_EQ;
						if (pch[3] == '=')
							buf.cond_tag_values.push_back(std::string(pch + 4));
						else
							buf.cond_tag_values.push_back(std::string(pch + 3));
					} else if (pch[2] == '!' && pch[3] == '=') {
						buf.cond_tag_op[n] = OP_NE;
						buf.cond_tag_values.push_back(std::string(pch + 4));
					} else if (pch[2] == '>' && pch[3] == '=') {
						buf.cond_tag_op[n] = OP_GE;
						buf.cond_tag_values.push_back(std::string(pch + 4));
					} else if (pch[2] == '<' && pch[3] == '=') {
						buf.cond_tag_op[n] = OP_LE;
						buf.cond_tag_values.push_back(std::string(pch + 4));
					} else if (pch[2] == '>' && pch[3] != '=') {
						buf.cond_tag_op[n] = OP_GT;
						buf.cond_tag_values.push_back(std::string(pch + 3));
					} else if (pch[2] == '<' && pch[3] != '=') {
						buf.cond_tag_op[n] = OP_LT;
						buf.cond_tag_values.push_back(std::string(pch + 3));
					}
				}
			} else {
				fprintf(stderr, "\nERROR: Unrecognized or unused parameter %s\n", pch);
				return;
			}
			// process argument
			pch = strtok(NULL, "&");
		}
		free(query);
	}
	if (buf.shift + buf.width < 0) {
		sqlite3_result_error(context, "Output of reads must cover variant location (width + start < 0)", -1);
		return;
	}
	if (info->reference && info->chrIdx >= 0){
		// get the reference genome
		try {
			buf.reference = info->reference->getSequence(
				cgatools::reference::Range(info->chrIdx, pos - 1 + buf.shift, pos - 1 + buf.shift + buf.width));
		} catch (cgatools::util::Exception & e) {
			// position out of range
			buf.reference.clear();
		}
	}
	// both start and end should be zero based.
	bam_fetch(sf->x.bam, idx, tid, pos - 1, pos,
		&buf, fetch_func);
	if (fi->name.empty() || fi->name == "coverage") {
		if (fi->all) {
			char tmp[48];
			sprintf(tmp, "%lu", buf.counter);
			sqlite3_result_text(context, tmp, -1, SQLITE_TRANSIENT);
		} else
			sqlite3_result_int64(context, buf.counter);
	} else if (fi->name == "qual") {
		std::stringstream res;
		std::vector<int>::iterator it = buf.qual.begin();
		std::vector<int>::iterator it_end = buf.qual.end();
		int index = 0;
		for (; it != it_end && (buf.limit < 0 || index < buf.limit); ++it, ++index) {
			if (index > 0)
				res << ",";
			res << *it;
		}
		sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
	} else if (fi->name == "avg_qual") {
		double res = 0;
		std::vector<int>::iterator it = buf.qual.begin();
		std::vector<int>::iterator it_end = buf.qual.end();
		for (; it != it_end; ++it)
			res += *it;
		res = buf.qual.size() == 0 ? 0 : res / buf.qual.size();
		if (fi->all) {
			char tmp[48];
			sprintf(tmp, "%f", res);
			sqlite3_result_text(context, tmp, -1, SQLITE_TRANSIENT);
		} else
			sqlite3_result_double(context, res);
	} else if (fi->name == "mapq") {
		std::stringstream res;
		std::vector<int>::iterator it = buf.map_qual.begin();
		std::vector<int>::iterator it_end = buf.map_qual.end();
		int index = 0;
		for (; it != it_end && (buf.limit < 0 || index < buf.limit); ++it, ++index) {
			if (index > 0)
				res << ",";
			res << *it;
		}
		sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
	} else if (fi->name == "avg_mapq") {
		double res = 0;
		std::vector<int>::iterator it = buf.map_qual.begin();
		std::vector<int>::iterator it_end = buf.map_qual.end();
		for (; it != it_end; ++it)
			res += *it;
		res = buf.map_qual.size() == 0 ? 0 : res / buf.map_qual.size();
		if (fi->all) {
			char tmp[48];
			sprintf(tmp, "%f", res);
			sqlite3_result_text(context, tmp, -1, SQLITE_TRANSIENT);
		} else
			sqlite3_result_double(context, res);
	} else
		sqlite3_result_text(context, (char *)(buf.calls.str().c_str()), -1, SQLITE_TRANSIENT);
}


extern "C" {
extern void bam_init_header_hash(bam_header_t *);

}
static void track(
                  sqlite3_context * context,
                  int argc,
                  sqlite3_value ** argv
                  )
{
	// do not check the first several parameters (variants) because
	// they are passed automatically and are assumed to be valid
	if (argc < 6 || sqlite3_value_type(argv[5]) == SQLITE_NULL) {
		sqlite3_result_error(context, "please specify at least filename", -1);
		return;
	} else if (argc > 8) {
		sqlite3_result_error(context, "track function accept at most 3 parameters", -1);
		return;
	}

	std::string track_file = std::string((char *)sqlite3_value_text(argv[5]));
	char * chr = (char *)sqlite3_value_text(argv[0]);
	int pos = sqlite3_value_int(argv[1]);
	char * ref = (char *)sqlite3_value_text(argv[2]);
	char * alt = (char *)sqlite3_value_text(argv[3]);
	std::string ref_file = std::string((char *)sqlite3_value_text(argv[4]));

	cgatools::reference::CrrFile * cf = NULL;
	int chrIdx = -1;
	if (!ref_file.empty()) {
		// get reference genome file
		RefGenomeFileMap::const_iterator rit = refFileMap.find(ref_file);

		if (rit == refFileMap.end()) {
			try {
				cf = new cgatools::reference::CrrFile(ref_file);
			} catch (cgatools::util::Exception & e) {
				sqlite3_result_error(context, e.what(), -1);
				return;
			}
			refFileMap[ref_file] = cf;
		} else {
			cf = rit->second;
		}
		//
		// get chromosome index
		ChrNameMap::const_iterator cit = chrNameMap.find(std::string(chr));
		if (cit == chrNameMap.end()) {
			try {
				// cgatools' chromosome names have leading chr
				// we add chr to 1, 2, 3, ... etc but not to other
				// names. There might be a problem here.
				if (strlen(chr) <= 2)
					chrIdx = cf->getChromosomeId("chr" + std::string(chr));
				else
					chrIdx = cf->getChromosomeId(std::string(chr));
			} catch (cgatools::util::Exception & e) {
				// unrecognized chromosome
				//sqlite3_result_error(context, e.what(), -1);
				chrIdx = -1;
			}
			chrNameMap[std::string(chr)] = chrIdx;
		} else {
			chrIdx = cit->second;
		}
	}

	TrackFileMap::const_iterator it = trackFileMap.find(track_file);
	TrackInfo info;
	if (it != trackFileMap.end())
		info = it->second;
	else {
		if (endsWith((char *)track_file.c_str(), ".vcf.gz")) {
			info.file_type = VCFTABIX_FILE;
			info.file = (void *)vcfTabixFileMayOpen((char *)track_file.c_str(),
				NULL, 0, 0, VCF_IGNORE_ERRS, 1);
			if (info.file == NULL) {
				sqlite3_result_error(context, "cannot open file", -1);
				return;
			}
			info.handler = vcfTabixTrack;
			info.default_col = 8;
			//
			info.name_map["chr"] = 1;
			info.name_map["chrom"] = 1;
			info.name_map["pos"] = 2;
			info.name_map["name"] = 3;
			info.name_map["ref"] = 4;
			info.name_map["alt"] = 5;
			info.name_map["qual"] = 6;
			info.name_map["filter"] = 7;
			info.name_map["info"] = 8;
			info.name_map["format"] = 9;
			// info fields
			struct vcfFile * vcff = (struct vcfFile *)info.file;
			struct vcfInfoDef * def = NULL;
			for (def = vcff->infoDefs; def != NULL; def = def->next)
				info.name_map[std::string("info.") + def->key] = 8;
			// sample
			for (size_t i = 0; i < vcff->genotypeCount; ++i) {
				info.name_map[vcff->genotypeIds[i]] = 10 + i;
				// sample genotype field
				struct vcfInfoDef * def = NULL;
				for (def = vcff->gtFormatDefs; def != NULL; def = def->next)
					info.name_map[std::string(vcff->genotypeIds[i]) + std::string(".") + def->key] = 10 + i;
			}
			// get the first record
			info.with_leading_chr = false;
			struct vcfRecord * rec = vcff->records;
			if (rec != NULL) {
				if (strncmp(rec->chrom, "chr", 3) == 0)
					info.with_leading_chr = true;
			}
		} else if (endsWith((char *)track_file.c_str(), ".bam")) {
			info.file_type = BAM_FILE;
			if (!bamFileExists((char *)track_file.c_str()))
				sqlite3_result_error(context, "bam file or its index does not exist", -1);
			char * filename;
			info.file = (void *)bamOpen((char *)track_file.c_str(), &filename);
			if (info.file == NULL) {
				sqlite3_result_error(context, "cannot open file", -1);
				return;
			}
			bam_init_header_hash(((samfile_t *)info.file)->header);
			info.index_file = bam_index_load(filename);
			if (info.index_file == NULL) {
				sqlite3_result_error(context, "cannot open file", -1);
				return;
			}
			info.handler = bamTrack;
			info.default_col = 0;
			// info fields
			//info.with_leading_chr = false;
			//bamChromList
			info.name_map["coverage"] = 1;
			info.name_map["calls"] = 1;
			info.name_map["reads"] = 1;
			info.name_map["qual"] = 1;
			info.name_map["mapq"] = 1;
			info.name_map["avg_qual"] = 1;
			info.name_map["avg_mapq"] = 1;
			info.name_map["flag"] = 1;
			// tags
			bam1_t data;
			bam1_t * bam = &data;
			ZeroVar(bam);
			if (bam_read1(((samfile_t *)info.file)->x.bam, bam) > 0) {
				uint8_t * s = bam1_aux(bam);
				while (s < bam->data + bam->data_len) {
					char key[3];
					key[0] = char(s[0]);
					key[1] = char(s[1]);
					key[2] = '\0';
					info.name_map[key] = 1;
					s += 2;
					uint8_t type = *s;
					++s;
					if (type == 'A' || type == 'C' || type == 'c')
						++s;
					else if (type == 'S' || type == 's')
						s += 2;
					else if (type == 'I' || type == 'i' || type == 'f')
						s += 4;
					else if (type == 'd')
						s += 8;
					else if (type == 'Z' || type == 'H')
						s += strlen((char *)s) + 1;
					else if (type == 'B') {
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
						int nItem = (int)*s;
						s += 4 + nItem * sz;
					}
				}
			}
			struct bamChromInfo * cl = bamChromList((samfile_t *)info.file);
			if (cl != NULL) {
				if (strncmp(cl->name, "chr", 3) == 0)
					info.with_leading_chr = true;
			}
		} else if (isBigWig((char *)track_file.c_str())) {
			info.file_type = BIGWIG_FILE;
			info.file = (void *)bigWigFileOpen((char *)track_file.c_str());
			if (info.file == NULL)
				sqlite3_result_error(context, "cannot open file", -1);
			info.handler = bigWigTrack;
			info.default_col = 4;
			//
			info.name_map["chrom"] = 1;
			info.name_map["chromStart"] = 2;
			info.name_map["chromEnd"] = 3;
			info.name_map["value"] = 4;
			//
			struct bbiChromInfo * chrom = bbiChromList((bbiFile *)info.file);
			info.with_leading_chr = startsWith("chr", chrom->name);
		} else if (bigBedFileCheckSigs((char *)track_file.c_str())) {
			info.file_type = BIGBED_FILE;
			info.file = (void *)bigBedFileOpen((char *)track_file.c_str());
			if (info.file == NULL) {
				sqlite3_result_error(context, "cannot open file", -1);
				return;
			}
			info.handler = bigBedTrack;
			info.default_col = 0;

			struct asObject * as = bigBedAs((bbiFile *)info.file);
			if (as != NULL) {
				size_t i = 1;
				struct asColumn * col;
				for (col = as->columnList; col != NULL; col = col->next, ++i)
					info.name_map[std::string(col->name)] = i;
			}
			//
			struct bbiChromInfo * chrom = bbiChromList((bbiFile *)info.file);
			info.with_leading_chr = startsWith("chr", chrom->name);
		} else {
			sqlite3_result_error(context, "Unknown track file type. Only local or remote tabix-indexed "
				                          " vcf files with extension .vcf.gz, indexed BAM files, bigWig and bigBed files are "
				                          "supported.", -1);
			return;
		}
		trackFileMap[track_file] = info;
		// info.print();
	}
	info.reference = cf;
	info.chrIdx = chrIdx;
	FieldInfo fi;
	fi.column = info.default_col;
	fi.name = "";
	fi.all = false;
	if (argc >= 7) {
		if (sqlite3_value_type(argv[6]) == SQLITE_INTEGER) {
			fi.column = sqlite3_value_int(argv[6]);
			fi.name = "";
		} else {
			char * name = (char *)(sqlite3_value_text(argv[6]));
			char * p = name;
			while (*p != '?' && *p != '\0')
				++p;
			if (*p == '?') {
				// terminate
				*p = '\0';
				fi.name = std::string(name);
				fi.query = std::string(p + 1);
				// look for all=1 and set fi.all
				//
				// NOTE:
				// We are being lazy here. The all option only works for VCF, bigWig and bigBed format
				// and these formats currently only accepts the 'all option'. The BAM format accepts
				// a bunch of options but it does not accept 'all', so we do not have to worry about
				// handling multiple options with 'all'
				if (strncmp(p+1, "all=", 4) == 0)
					fi.all = atoi(p + 5);
				// replace ? to avoid disrupting the string
				*p = '?';
			} else
				fi.name = std::string(name);
			std::map<std::string, int>::iterator it = info.name_map.find(fi.name);
			if (it == info.name_map.end()) {
				char buf[400];
				sprintf(buf, "Unrecognized field '%s', please check available fields with command 'vtools show track %s'."
					         " You might also have used double quote (\") instead of single quote (') for string literal.",
					fi.name.c_str(), track_file.c_str());
				sqlite3_result_error(context, buf, -1);
				return;
			}
			fi.column = it->second;
		}
	}
	// this option is now obsolete
	if (argc >= 8) {
		if (sqlite3_value_type(argv[7]) == SQLITE_INTEGER)
			fi.all = sqlite3_value_int(argv[7]);
		else {
			sqlite3_result_error(context, "wrong datatype for the last parameter. 0 or 1 (all records) is expected.", -1);
			return;
		}
	}

	// call the handler
	(*info.handler)(info.file, chr, pos, ref, alt, &fi, &info, context);
}


static void in_table(
                     sqlite3_context * context,
                     int argc,
                     sqlite3_value ** argv
                     )
{
	// parameters passed:
	// name of variant table
	int variant_id = sqlite3_value_int(argv[0]);
	char * var_table = (char *)sqlite3_value_text(argv[1]);

	char sql[255];
	sprintf(sql, "SELECT 1 FROM %s WHERE variant_id = %d LIMIT 0,1",
			var_table, variant_id);
	sqlite3_stmt * stmt;
	int result = sqlite3_prepare_v2(sqlite3_context_db_handle(context), sql, -1, &stmt, NULL);
	if (result != SQLITE_OK) {
		//sqlite3_result_error(context, sqlite3_errmsg(geno_db), -1);
		sqlite3_result_null(context);
		return;
	}
	// there should be zero or one matching record
	result = sqlite3_step(stmt);
	sqlite3_result_int(context, result == SQLITE_ROW);
	return;
}




sqlite3 * geno_db;
typedef std::map<std::string, int> SampleNameIdMap;
SampleNameIdMap nameIdMap;

class genotypeParams 
{
public:
	genotypeParams(const char * params = NULL) :
		m_params(NULL), m_field(NULL), m_delimiter(NULL), m_missing(NULL)
	{
		if (!params)
			return;

		m_params = strdup(params);
		char * pch = strtok(m_params, "&");
		while (pch != NULL) {
			if (strncmp(pch, "field=", 6) == 0)
				m_field = pch + 6;
			else if (strncmp(pch, "d=", 2) == 0) {
				m_delimiter = pch + 2;
				if (strcmp(m_delimiter, "\\t") == 0)
					m_delimiter = "\t";
			} else if (strncmp(pch, "delimiter=", 10) == 0) {
				m_delimiter = pch + 10;
				if (strcmp(m_delimiter, "\\t") == 0)
					m_delimiter = "\t";
				else if (strcmp(m_delimiter, "\\n") == 0)
					m_delimiter = "\n";
			} else if (strncmp(pch, "missing=", 8) == 0)
				m_missing = pch + 8;
			else
				fprintf(stderr, "Incorrect parameter to function genotype: %s", pch);
			// process argument
			pch = strtok(NULL, "&");
		}
	}
	
	~genotypeParams()
	{
		free(m_params);
	}

	char * field()
	{
		return m_field ? m_field : (char*)"GT";
	}

	char * delimiter()
	{
		return m_delimiter ? m_delimiter : (char*)",";
	}

	char * missing()
	{
		return m_missing;
	}

private:
	// this one holds the copied stuff
	char * m_params;
	// the rest are just pointers
	char * m_field;
	char * m_delimiter;
	char * m_missing;
};


static void genotype(
                     sqlite3_context * context,
                     int argc,
                     sqlite3_value ** argv
                     )
{
	// parameters passed:
	// name of geno_db_file, variant_id, sample_name, and then params
	if (argc < 3 || sqlite3_value_type(argv[2]) == SQLITE_NULL) {
		sqlite3_result_error(context, "please name of a sample", -1);
		return;
	} else if (argc > 4) {
		sqlite3_result_error(context, "Function genotype accepts at most 2 parameter sample_name and field", -1);
		return;
	}

	char * geno_db_file = (char *)sqlite3_value_text(argv[0]);
	int variant_id = sqlite3_value_int(argv[1]);
	// there are two cases, if a single integer is passed
	// otherwise a list of IDs are saved in a file
	bool single_ID = sqlite3_value_type(argv[2]) == SQLITE_INTEGER;
	std::vector<int> sample_IDs;
	if (single_ID) {
		sample_IDs.push_back(sqlite3_value_int(argv[2]));
	} else {
		std::ifstream ids((const char *)sqlite3_value_text(argv[2]));
		while (ids.good()) {
			// read each line and get sample IDs
			int id = -1;
			ids >> id;
			if (id != -1)
				sample_IDs.push_back(id);
		}
	}

	genotypeParams params(argc == 4 ? (char *)sqlite3_value_text(argv[3]) : NULL);

	// open databases
	if (!geno_db) {
		int result = sqlite3_open_v2(geno_db_file, &geno_db, SQLITE_OPEN_READONLY, NULL);
		if (result != SQLITE_OK) {
			sqlite3_result_error(context, "Failed to open genotype database", -1);
			return;
		}
	}

	if (single_ID) {
		int result = 0;
		// run some query
		char sql[255];
		sprintf(sql, "SELECT %s FROM genotype_%d WHERE variant_id = %d LIMIT 0,1",
			params.field(), sample_IDs[0], variant_id);
		sqlite3_stmt * stmt;
		result = sqlite3_prepare_v2(geno_db, sql, -1, &stmt, NULL) ;
		if (result != SQLITE_OK) {
			//sqlite3_result_error(context, sqlite3_errmsg(geno_db), -1);
			sqlite3_result_null(context);
			return;
		}
		// there should be only one matching record
		result = sqlite3_step(stmt);
		if (result == SQLITE_ROW) {
			// how to pass whatever type the query gets to the output???
			switch (sqlite3_column_type(stmt, 0)) {
			case SQLITE_INTEGER:
				sqlite3_result_int(context, sqlite3_column_int(stmt, 0));
				break;
			case SQLITE_FLOAT:
				sqlite3_result_double(context, sqlite3_column_double(stmt, 0));
				break;
			case SQLITE_TEXT:
				sqlite3_result_text(context, (const char *)sqlite3_column_text(stmt, 0), -1, SQLITE_TRANSIENT);
				break;
			case SQLITE_BLOB:
				sqlite3_result_blob(context, sqlite3_column_blob(stmt, 0), -1, SQLITE_TRANSIENT);
				break;
			case SQLITE_NULL:
				sqlite3_result_null(context);
				break;
			}
		} else
			sqlite3_result_null(context);
		return;
	}
	// if there are multiple IDs
	// go through all samples (with id)
	std::stringstream res;
	bool first = true;
	std::vector<int>::iterator it = sample_IDs.begin();
	std::vector<int>::iterator it_end = sample_IDs.end();
	for (; it != it_end; ++it) {
		// run some query
		char sql[255];
		sprintf(sql, "SELECT %s FROM genotype_%d WHERE variant_id = %d LIMIT 0,1",
			params.field(), *it, variant_id);
		sqlite3_stmt * stmt;
		int result = sqlite3_prepare_v2(geno_db, sql, -1, &stmt, NULL) ;
		if (result != SQLITE_OK) {
			if (params.missing() != NULL)
				res << params.missing();
			continue;
		}
		// there should be only one matching record
		result = sqlite3_step(stmt);
		if (result == SQLITE_ROW) {
			if (first)
				first = false;
			else
				res << params.delimiter();
			// how to pass whatever type the query gets to the output???
			switch (sqlite3_column_type(stmt, 0)) {
			case SQLITE_INTEGER:
				res << sqlite3_column_int(stmt, 0);
				break;
			case SQLITE_FLOAT:
				res << sqlite3_column_double(stmt, 0);
				break;
			case SQLITE_TEXT:
				res << (const char *)sqlite3_column_text(stmt, 0);
				break;
			case SQLITE_BLOB:
				res << sqlite3_column_blob(stmt, 0);
				break;
			case SQLITE_NULL:
				res << params.missing();
				break;
			}
		} else if (params.missing() != NULL) {
			if (first)
				first = false;
			else
				res << params.delimiter();
			res << params.missing();
		}
	}
	// we do not close the database because we are readonly and we need the database for
	// next use.
	if (res.str().empty())
		sqlite3_result_null(context);
	else
		sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
}


typedef std::map<int, std::string> IdNameMap;
typedef std::map<std::string, IdNameMap> SampleIdNameMap;
SampleIdNameMap idNameMap;


class samplesParams 
{
public:
	samplesParams(const char * params = NULL) :
		m_params(NULL), m_geno_filter(NULL), m_delimiter(NULL)
	{
		if (!params)
			return;

		m_params = strdup(params);
		char * pch = strtok(m_params, "&");
		while (pch != NULL) {
			if (strncmp(pch, "geno_filter=", 12) == 0 && (*(pch + 12) != '\0'))
				m_geno_filter = pch + 12;
			//else if (strncmp(pch, "sample_filter=", 14) == 0)
			//	m_sample_filter = pch + 14;
			else if (strncmp(pch, "d=", 2) == 0) {
				m_delimiter = pch + 2;
				if (strcmp(m_delimiter, "\\t") == 0)
					m_delimiter = "\t";
			} else
				fprintf(stderr, "Incorrect parameter to function samples: %s", pch);
			// process argument
			pch = strtok(NULL, "&");
		}
	}
	
	~samplesParams()
	{
		free(m_params);
	}

	char * geno_filter()
	{
		return m_geno_filter; 
	}

	//char * sample_filter()
	//{
	//	return m_sample_filter;
	//}

	char * delimiter()
	{
		return m_delimiter ? m_delimiter : (char*)",";
	}

private:
	// this one holds the copied stuff
	char * m_params;
	// the rest are just pointers
	char * m_geno_filter;
	// char * m_sample_filter;
	char * m_delimiter;
};



static void samples(
                    sqlite3_context * context,
                    int argc,
                    sqlite3_value ** argv
                    )
{
	// parameters passed:
	// name of project, variant_id, sample_id_file, [genotype_filter]
	//
	if (argc > 4) {
		sqlite3_result_error(context, "samples function accept at most one parameter", -1);
		return;
	}

	char * geno_db_file = (char *)sqlite3_value_text(argv[0]);
	int variant_id = sqlite3_value_int(argv[1]);
	char * sample_id_file = (char *)sqlite3_value_text(argv[2]);
	samplesParams params(argc > 3 ? (char *)sqlite3_value_text(argv[3]) : NULL);

	// see if the sample_id_file has been loaded
	SampleIdNameMap::iterator it = idNameMap.find(std::string(sample_id_file));
	if (it == idNameMap.end()) {
		// read that file
		IdNameMap inm;
		std::ifstream ids(sample_id_file);
		while (ids.good()) {
			// read each line and fill nameIdMap
			int id = -1;
			std::string name;
			ids >> id >> name;
			if (id != -1)
				inm[id] = name;
		}
		idNameMap[std::string(sample_id_file)] = inm;
		// find the item
		it = idNameMap.find(std::string(sample_id_file));
	}
	IdNameMap & idMap = it->second;

	int result = 0;
	// open databases
	if (!geno_db) {
		result = sqlite3_open_v2(geno_db_file, &geno_db, SQLITE_OPEN_READONLY, NULL);
		if (result != SQLITE_OK) {
			sqlite3_result_error(context, "Failed to open genotype database", -1);
			return;
		}
	}
	//
	// go through all samples (with id)
	std::stringstream res;
	IdNameMap::iterator im = idMap.begin();
	IdNameMap::iterator im_end = idMap.end();
	bool first = true;
	for (; im != im_end; ++im) {
		char sql[255];
		sprintf(sql, "SELECT variant_id FROM genotype_%d WHERE variant_id = %d AND (%s) LIMIT 0,1",
			im->first, variant_id, params.geno_filter() == NULL ? "1" : params.geno_filter());
		//
		sqlite3_stmt * stmt;
		result = sqlite3_prepare_v2(geno_db, sql, -1, &stmt, NULL) ;
		if (result != SQLITE_OK) {
			sqlite3_result_error(context, sqlite3_errmsg(geno_db), -1);
			return;
		}
		result = sqlite3_step(stmt);
		if (result == SQLITE_ROW) {
			if (first)
				first = false;
			else
				res << params.delimiter();
			res << im->second;
		}
	}
	if (res.str().empty())
		sqlite3_result_null(context);
	else
		sqlite3_result_text(context, (char *)(res.str().c_str()), -1, SQLITE_TRANSIENT);
	// output a string with all information
	// we do not close the database because we are readonly and we need the database for
	// next use.
}


/*
** HWE exact test for a bi-allelic locus
*/
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
static void hwe_exact(
                      sqlite3_context * context,
                      int argc,
                      sqlite3_value ** argv
                      )
{
	/*
	   consider this 2X2 table for exact HWE test
	       C       A
	    C  n11     n12_1
	    A  n12_2   n22
	    then n11 = #(CC); n12 = #(CA) + #(AC); n22 = #(AA)
	 */
	double n = sqlite3_value_double(argv[0]);
	double n12 = sqlite3_value_double(argv[1]);
	// collapse double homozygotes to simply homozygote
	double n22 = sqlite3_value_double(argv[2]);

	double n11, n1, n2, pn12, pval, x12, x11, x22, x1, x2, px12;

	if (argc == 4) {
		n22 += sqlite3_value_double(argv[3]);
	}
	n11 = n - n12 - n22;
	//http://www.ncbi.nlm.nih.gov/pmc/articles/PMC1199378/
	//implements equation 2 below
	n1 = 2.0 * n11 + n12;
	n2 = 2.0 * n22 + n12;
	pn12 = exp(log(2.0) * (n12) + gsl_sf_lngamma(n + 1) -
		gsl_sf_lngamma(n11 + 1) - gsl_sf_lngamma(n12 + 1) -
		gsl_sf_lngamma(n22 + 1) - gsl_sf_lngamma(2.0 * n + 1) +
		gsl_sf_lngamma(n1 + 1) + gsl_sf_lngamma(n2 + 1));
	//
	pval = 0.0;
	for (x12 = fmod(n1, 2.0); x12 <= (n1 < n2 ? n1 : n2); x12 = x12 + 2.0) {
		x11 = (n1 - x12) / 2.0;
		x22 = (n2 - x12) / 2.0;
		x1 = 2.0 * x11 + x12;
		x2 = 2.0 * x22 + x12;
		px12 = exp(log(2.0) * (x12) + gsl_sf_lngamma(n + 1) -
			gsl_sf_lngamma(x11 + 1) - gsl_sf_lngamma(x12 + 1) -
			gsl_sf_lngamma(x22 + 1) - gsl_sf_lngamma(2.0 * n + 1) +
			gsl_sf_lngamma(x1 + 1) + gsl_sf_lngamma(x2 + 1));
		if (pn12 >= px12) pval += px12;
	}
	sqlite3_result_double(context, pval);
}


#include "fisher2.h"
static void fisher_exact(
                         sqlite3_context * context,
                         int argc,
                         sqlite3_value ** argv
                         )
{

	/* This is specific for 2x2 tables. contingency_table = matrix(twotwoTable, 2, 2, byrow = T)
	   case       ctrl
	   alt m1         m2
	   ref n1-m1      n2-m2
	 */
	double contingency_table[4] = { 0, 0, 0, 0 };
	double pval, expected, percnt, emin, prt;
	int ok, nrow, ncol, workspace;

	contingency_table[0] = (int)(sqlite3_value_double(argv[0]));
	contingency_table[1] = (int)(sqlite3_value_double(argv[1]));
	contingency_table[2] = (int)(sqlite3_value_double(argv[2])) - contingency_table[0];
	contingency_table[3] = (int)(sqlite3_value_double(argv[3])) - contingency_table[1];
	pval = -99.0;
	ok = (
	      contingency_table[0] >= 0 &&
	      contingency_table[1] >= 0 &&
	      contingency_table[2] >= 0 &&
	      contingency_table[3] >= 0 &&
	      (contingency_table[0] + contingency_table[1] +
	       contingency_table[2] + contingency_table[0] > 0)
	      );
	if (ok) {
		nrow = 2;
		ncol = 2;
		expected = -1.0;
		percnt = 100.0;
		emin = 0.0;
		prt = 0.0;
		workspace = 300000;
		fexact(&nrow, &ncol, contingency_table, &nrow, &expected, &percnt, &emin, &prt, &pval, &workspace);
	}
	sqlite3_result_double(context, pval);
}


// this is a fake function to make this .so file a loadable Python module
// so that it can be imported and recogznied by pyinstaller. The library
// is actually imported into sqlite using SELECT load_extension()

// suppress two warnings
#ifdef _GNU_SOURCE
#  undef _GNU_SOURCE
#endif
#ifdef _REENTRANT
#  undef _REENTRANT
#endif
#include <Python.h>
#if PY_VERSION_HEX >= 0x03000000
struct module_state
{
	PyObject * error;
};

static PyMethodDef myextension_methods[] = {
	{ NULL, NULL }
};

/*
   static int myextension_traverse(PyObject *m, visitproc visit, void *arg)
   {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
   }

   static int myextension_clear(PyObject *m)
   {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
   } */
#endif

PyMODINIT_FUNC init_vt_sqlite3_ext()
{
#if PY_VERSION_HEX >= 0x03000000
	static struct PyModuleDef moduledef = {
		PyModuleDef_HEAD_INIT,
		"_vt_sqlite3_ext",
		NULL,
		sizeof(struct module_state),
		myextension_methods,
		NULL,
		NULL,                       //myextension_traverse,
		NULL,                       //myextension_clear,
		NULL
	};

	PyModule_Create(&moduledef);
#else
	Py_InitModule3("_vt_sqlite3_ext", NULL, "Fake sqlite3 extension module");
#endif
	// initialize global variable
	geno_db = NULL;
	//return NULL;
}


/* SQLite invokes sqlite3_extension_init, which is defined later
** and will call this function to register variant tools defined
** functions.
*/
int sqlite3_my_extension_init(
                              sqlite3 * db,
                              char ** pzErrMsg,
                              const sqlite3_api_routines * pApi
                              )
{
	SQLITE_EXTENSION_INIT2(pApi)
	//http://www.sqlite.org/c3ref/create_function.html
	//The first parameter is the database connection to which the SQL function is to be added
	//The second parameter is the name of the SQL function to be created or redefined.
	//The third parameter (nArg) is the number of arguments that the SQL function or aggregate takes.
	//The fourth parameter, eTextRep, specifies what text encoding this SQL function prefers for its parameters.
	//The fifth parameter is an arbitrary pointer.
	//The sixth, seventh and eighth parameters, xFunc, xStep and xFinal, are pointers to C-language functions that implement the SQL function or aggregate.
	sqlite3_create_function(db, "least_not_null", -1, SQLITE_ANY, 0, least_not_null_func, 0, 0);
	sqlite3_create_function(db, "HWE_exact", -1, SQLITE_ANY, 0, hwe_exact, 0, 0);
	sqlite3_create_function(db, "Fisher_exact", 4, SQLITE_ANY, 0, fisher_exact, 0, 0);
	// ref_sequence(file, chr, start, end)
	sqlite3_create_function(db, "ref_sequence", -1, SQLITE_ANY, 0, ref_sequence, 0, 0);
	// mut_sequence(file, chr_allele, pos_allele, ref_allele, alt_allele, chr, start, end)
	sqlite3_create_function(db, "mut_sequence", -1, SQLITE_ANY, 0, mut_sequence, 0, 0);
	// pad_variant(file, chr, pos, ref, alt, name)  ==> 10 - A ==> 9 name G GA
	sqlite3_create_function(db, "vcf_variant", -1, SQLITE_ANY, 0, vcf_variant, 0, 0);
	sqlite3_create_function(db, "track", -1, SQLITE_ANY, 0, track, 0, 0);
	sqlite3_create_function(db, "in_table", 2, SQLITE_ANY, 0, in_table, 0, 0);
	sqlite3_create_function(db, "genotype", -1, SQLITE_ANY, 0, genotype, 0, 0);
	sqlite3_create_function(db, "samples", -1, SQLITE_ANY, 0, samples, 0, 0);
	return 0;
}


/*
 * The following code is from extension-functions.c (http://www.sqlite.org/contrib)
 *
 * The only change is the removal of duplicated lines
 *
 *   #include "sqlite3ext.h"
 *   SQLITE_EXTENSION_INIT1
 *
 * and change sqlite3_extension_init to call sqlite3_my_extension_init
 * which initialize the functions added by variant tools.
 */

/*
   This library will provide common mathematical and string functions in
   SQL queries using the operating system libraries or provided
   definitions.  It includes the following functions:

   Math: acos, asin, atan, atn2, atan2, acosh, asinh, atanh, difference,
   degrees, radians, cos, sin, tan, cot, cosh, sinh, tanh, coth, exp,
   log, log10, power, sign, sqrt, square, ceil, floor, pi.

   String: replicate, charindex, leftstr, rightstr, ltrim, rtrim, trim,
   replace, reverse, proper, padl, padr, padc, strfilter.

   Aggregate: stdev, variance, mode, median, lower_quartile,
   upper_quartile.

   The string functions ltrim, rtrim, trim, replace are included in
   recent versions of SQLite and so by default do not build.

   Compilation instructions:
   Compile this C source file into a dynamic library as follows:
 * Linux:
   gcc -fPIC -lm -shared extension-functions.c -o libsqlitefunctions.so
 * Mac OS X:
   gcc -fno-common -dynamiclib extension-functions.c -o libsqlitefunctions.dylib
   (You may need to add flags
   -I /opt/local/include/ -L/opt/local/lib -lsqlite3
   if your sqlite3 is installed from Mac ports, or
   -I /sw/include/ -L/sw/lib -lsqlite3
   if installed with Fink.)
 * Windows:
   1. Install MinGW (http://www.mingw.org/) and you will get the gcc
   (gnu compiler collection)
   2. add the path to your path variable (isn't done during the
   installation!)
   3. compile:
   gcc -shared -I "path" -o libsqlitefunctions.so extension-functions.c
   (path = path of sqlite3ext.h; i.e. C:\programs\sqlite)

   Usage instructions for applications calling the sqlite3 API functions:
   In your application, call sqlite3_enable_load_extension(db,1) to
   allow loading external libraries.  Then load the library libsqlitefunctions
   using sqlite3_load_extension; the third argument should be 0.
   See http://www.sqlite.org/cvstrac/wiki?p=LoadableExtensions.
   Select statements may now use these functions, as in
   SELECT cos(radians(inclination)) FROM satsum WHERE satnum = 25544;

   Usage instructions for the sqlite3 program:
   If the program is built so that loading extensions is permitted,
   the following will work:
   sqlite> SELECT load_extension('./libsqlitefunctions.so');
   sqlite> select cos(radians(45));
   0.707106781186548
   Note: Loading extensions is by default prohibited as a
   security measure; see "Security Considerations" in
   http://www.sqlite.org/cvstrac/wiki?p=LoadableExtensions.
   If the sqlite3 program and library are built this
   way, you cannot use these functions from the program, you
   must write your own program using the sqlite3 API, and call
   sqlite3_enable_load_extension as described above, or else
   rebuilt the sqlite3 program to allow loadable extensions.

   Alterations:
   The instructions are for Linux, Mac OS X, and Windows; users of other
   OSes may need to modify this procedure.  In particular, if your math
   library lacks one or more of the needed trig or log functions, comment
   out the appropriate HAVE_ #define at the top of file.  If you do not
   wish to make a loadable module, comment out the define for
   COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE.  If you are using a
   version of SQLite without the trim functions and replace, comment out
   the HAVE_TRIM #define.

   Liam Healy

   History:
   2010-01-06 Correct check for argc in squareFunc, and add Windows
   compilation instructions.
   2009-06-24 Correct check for argc in properFunc.
   2008-09-14 Add check that memory was actually allocated after
   sqlite3_malloc or sqlite3StrDup, call sqlite3_result_error_nomem if
   not.  Thanks to Robert Simpson.
   2008-06-13 Change to instructions to indicate use of the math library
   and that program might work.
   2007-10-01 Minor clarification to instructions.
   2007-09-29 Compilation as loadable module is optional with
   COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE.
   2007-09-28 Use sqlite3_extension_init and macros
   SQLITE_EXTENSION_INIT1, SQLITE_EXTENSION_INIT2, so that it works with
   sqlite3_load_extension.  Thanks to Eric Higashino and Joe Wilson.
   New instructions for Mac compilation.
   2007-09-17 With help from Joe Wilson and Nuno Luca, made use of
   external interfaces so that compilation is no longer dependent on
   SQLite source code.  Merged source, header, and README into a single
   file.  Added casts so that Mac will compile without warnings (unsigned
   and signed char).
   2007-09-05 Included some definitions from sqlite 3.3.13 so that this
   will continue to work in newer versions of sqlite.  Completed
   description of functions available.
   2007-03-27 Revised description.
   2007-03-23 Small cleanup and a bug fix on the code.  This was mainly
   letting errno flag errors encountered in the math library and checking
   the result, rather than pre-checking.  This fixes a bug in power that
   would cause an error if any non-positive number was raised to any
   power.
   2007-02-07 posted by Mikey C to sqlite mailing list.
   Original code 2006 June 05 by relicoder.

 */

//#include "config.h"

#define COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE 1
#define HAVE_ACOSH 1
#define HAVE_ASINH 1
#define HAVE_ATANH 1
#define HAVE_SINH 1
#define HAVE_COSH 1
#define HAVE_TANH 1
#define HAVE_LOG10 1
#define HAVE_ISBLANK 1
#define SQLITE_SOUNDEX 1
#define HAVE_TRIM 1     /* LMH 2007-03-25 if sqlite has trim functions */

#ifdef COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE
/*
   #include "sqlite3ext.h"
   SQLITE_EXTENSION_INIT1
 */
#else
#  include "sqlite3.h"
#endif

#include <ctype.h>
/* relicoder */
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>      /* LMH 2007-03-25 */

#include <stdlib.h>
#include <assert.h>

#ifndef _MAP_H_
#  define _MAP_H_

#  include <stdint.h>

/*
** Simple binary tree implementation to use in median, mode and quartile calculations
** Tree is not necessarily balanced. That would require something like red&black trees of AVL
*/

typedef int (* cmp_func)(const void *, const void *);
typedef void (* map_iterator)(void *, int64_t, void *);

typedef struct node
{
	struct node * l;
	struct node * r;
	void * data;
	int64_t count;
} node;

typedef struct map
{
	node * base;
	cmp_func cmp;
	short free;
} map;

/*
** creates a map given a comparison function
*/
map map_make(cmp_func cmp);

/*
** inserts the element e into map m
*/
void map_insert(map * m, void * e);

/*
** executes function iter over all elements in the map, in key increasing order
*/
void map_iterate(map * m, map_iterator iter, void * p);

/*
** frees all memory used by a map
*/
void map_destroy(map * m);

/*
** compares 2 integers
** to use with map_make
*/
int int_cmp(const void * a, const void * b);

/*
** compares 2 doubles
** to use with map_make
*/
int double_cmp(const void * a, const void * b);

#endif /* _MAP_H_ */

typedef uint8_t u8;
typedef uint16_t u16;
typedef int64_t i64;

static char * sqlite3StrDup(const char * z)
{
	char * res = (char *)sqlite3_malloc(strlen(z) + 1);

	return strcpy(res, z);
}


/*
** These are copied verbatim from fun.c so as to not have the names exported
*/

/* LMH from sqlite3 3.3.13 */
/*
** This table maps from the first byte of a UTF-8 character to the number
** of trailing bytes expected. A value '4' indicates that the table key
** is not a legal first byte for a UTF-8 character.
*/
static const u8 xtra_utf8_bytes[256] = {
	/* 0xxxxxxx */
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,

	/* 10wwwwww */
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,

	/* 110yyyyy */
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,

	/* 1110zzzz */
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,

	/* 11110yyy */
	3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
};


/*
** This table maps from the number of trailing bytes in a UTF-8 character
** to an integer constant that is effectively calculated for each character
** read by a naive implementation of a UTF-8 character reader. The code
** in the READ_UTF8 macro explains things best.
*/
static const int xtra_utf8_bits[] = {
	0,
	12416,      /* (0xC0 << 6) + (0x80) */
	925824,     /* (0xE0 << 12) + (0x80 << 6) + (0x80) */
	63447168    /* (0xF0 << 18) + (0x80 << 12) + (0x80 << 6) + 0x80 */
};

/*
** If a UTF-8 character contains N bytes extra bytes (N bytes follow
** the initial byte so that the total character length is N+1) then
** masking the character with utf8_mask[N] must produce a non-zero
** result.  Otherwise, we have an (illegal) overlong encoding.
*/
static const unsigned int utf_mask[] = {
	0x00000000,
	0xffffff80,
	0xfffff800,
	0xffff0000,
};

/* LMH salvaged from sqlite3 3.3.13 source code src/utf.c */
#define READ_UTF8(zIn, c) { \
		int xtra;                                            \
		c = *(zIn)++;                                        \
		xtra = xtra_utf8_bytes[c];                           \
		switch (xtra) {                                      \
		case 4: c = (int)0xFFFD; break;                    \
		case 3: c = (c << 6) + *(zIn)++;                     \
		case 2: c = (c << 6) + *(zIn)++;                     \
		case 1: c = (c << 6) + *(zIn)++;                     \
			c -= xtra_utf8_bits[xtra];                         \
			if ( (utf_mask[xtra] & c) == 0                          \
			    || (c & 0xFFFFF800) == 0xD800                      \
			    || (c & 0xFFFFFFFE) == 0xFFFE) {  c = 0xFFFD; }    \
		}                                                    \
}

static int sqlite3ReadUtf8(const unsigned char * z)
{
	int c;

	READ_UTF8(z, c);
	return c;
}


#define SKIP_UTF8(zIn) {                               \
		zIn += (xtra_utf8_bytes[*(u8 *)zIn] + 1);            \
}

/*
** pZ is a UTF-8 encoded unicode string. If nByte is less than zero,
** return the number of unicode characters in pZ up to (but not including)
** the first 0x00 byte. If nByte is not less than zero, return the
** number of unicode characters in the first nByte of pZ (or up to
** the first 0x00, whichever comes first).
*/
static int sqlite3Utf8CharLen(const char * z, int nByte)
{
	int r = 0;
	const char * zTerm;

	if (nByte >= 0) {
		zTerm = &z[nByte];
	}else{
		zTerm = (const char *)(-1);
	}
	assert(z <= zTerm);
	while (*z != 0 && z < zTerm) {
		SKIP_UTF8(z);
		r++;
	}
	return r;
}


/*
** X is a pointer to the first byte of a UTF-8 character.  Increment
** X so that it points to the next character.  This only works right
** if X points to a well-formed UTF-8 string.
*/
#define sqliteNextChar(X)  while ( (0xc0 & *++(X)) == 0x80) {}
#define sqliteCharVal(X)   sqlite3ReadUtf8(X)

/*
** This is a macro that facilitates writting wrappers for math.h functions
** it creates code for a function to use in SQlite that gets one numeric input
** and returns a floating point value.
**
** Could have been implemented using pointers to functions but this way it's inline
** and thus more efficient. Lower * ranking though...
**
** Parameters:
** name:      function name to de defined (eg: sinFunc)
** function:  function defined in math.h to wrap (eg: sin)
** domain:    boolean condition that CAN'T happen in terms of the input parameter rVal
**            (eg: rval<0 for sqrt)
*/
/* LMH 2007-03-25 Changed to use errno and remove domain; no pre-checking for errors. */
#define GEN_MATH_WRAP_DOUBLE_1(name, function) \
	static void name(sqlite3_context * context, int argc, sqlite3_value * *argv){ \
		double rVal = 0.0, val; \
		assert(argc == 1); \
		switch (sqlite3_value_type(argv[0]) ) { \
		case SQLITE_NULL: { \
			sqlite3_result_null(context); \
			break; \
		} \
		default: { \
			rVal = sqlite3_value_double(argv[0]); \
			errno = 0; \
			val = function(rVal); \
			if (errno == 0) { \
				sqlite3_result_double(context, val); \
			} else { \
				sqlite3_result_error(context, strerror(errno), errno); \
			} \
			break; \
		} \
		} \
	} \


/*
** Example of GEN_MATH_WRAP_DOUBLE_1 usage
** this creates function sqrtFunc to wrap the math.h standard function sqrt(x)=x^0.5
*/
GEN_MATH_WRAP_DOUBLE_1(sqrtFunc, sqrt)

/* trignometric functions */
GEN_MATH_WRAP_DOUBLE_1(acosFunc, acos)
GEN_MATH_WRAP_DOUBLE_1(asinFunc, asin)
GEN_MATH_WRAP_DOUBLE_1(atanFunc, atan)

/*
** Many of systems don't have inverse hyperbolic trig functions so this will emulate
** them on those systems in terms of log and sqrt (formulas are too trivial to demand
** written proof here)
*/

#ifndef HAVE_ACOSH
static double acosh(double x)
{
	return log(x + sqrt(x * x - 1.0));
}


#endif

GEN_MATH_WRAP_DOUBLE_1(acoshFunc, acosh)

#ifndef HAVE_ASINH
static double asinh(double x)
{
	return log(x + sqrt(x * x + 1.0));
}


#endif

GEN_MATH_WRAP_DOUBLE_1(asinhFunc, asinh)

#ifndef HAVE_ATANH
static double atanh(double x)
{
	return (1.0 / 2.0) * log((1 + x) / (1 - x)) ;
}


#endif

GEN_MATH_WRAP_DOUBLE_1(atanhFunc, atanh)

/*
** math.h doesn't require cot (cotangent) so it's defined here
*/
static double cot(double x)
{
	return 1.0 / tan(x);
}


GEN_MATH_WRAP_DOUBLE_1(sinFunc, sin)
GEN_MATH_WRAP_DOUBLE_1(cosFunc, cos)
GEN_MATH_WRAP_DOUBLE_1(tanFunc, tan)
GEN_MATH_WRAP_DOUBLE_1(cotFunc, cot)

static double coth(double x)
{
	return 1.0 / tanh(x);
}


/*
** Many systems don't have hyperbolic trigonometric functions so this will emulate
** them on those systems directly from the definition in terms of exp
*/
#ifndef HAVE_SINH
static double sinh(double x)
{
	return (exp(x) - exp(-x)) / 2.0;
}


#endif

GEN_MATH_WRAP_DOUBLE_1(sinhFunc, sinh)

#ifndef HAVE_COSH
static double cosh(double x)
{
	return (exp(x) + exp(-x)) / 2.0;
}


#endif

GEN_MATH_WRAP_DOUBLE_1(coshFunc, cosh)

#ifndef HAVE_TANH
static double tanh(double x)
{
	return sinh(x) / cosh(x);
}


#endif

GEN_MATH_WRAP_DOUBLE_1(tanhFunc, tanh)

GEN_MATH_WRAP_DOUBLE_1(cothFunc, coth)

/*
** Some systems lack log in base 10. This will emulate it
*/

#ifndef HAVE_LOG10
static double log10(double x)
{
	static double l10 = -1.0;

	if (l10 < 0.0) {
		l10 = log(10.0);
	}
	return log(x) / l10;
}


#endif

GEN_MATH_WRAP_DOUBLE_1(logFunc, log)
GEN_MATH_WRAP_DOUBLE_1(log10Func, log10)
GEN_MATH_WRAP_DOUBLE_1(expFunc, exp)

/*
** Fallback for systems where math.h doesn't define M_PI
*/
#undef M_PI
#ifndef M_PI
/*
** static double PI = acos(-1.0);
** #define M_PI (PI)
*/
#  define M_PI 3.14159265358979323846
#endif

/* Convert Degrees into Radians */
static double deg2rad(double x)
{
	return x * M_PI / 180.0;
}


/* Convert Radians into Degrees */
static double rad2deg(double x)
{
	return 180.0 * x / M_PI;
}


GEN_MATH_WRAP_DOUBLE_1(rad2degFunc, rad2deg)
GEN_MATH_WRAP_DOUBLE_1(deg2radFunc, deg2rad)

/* constant function that returns the value of PI=3.1415... */
static void piFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	sqlite3_result_double(context, M_PI);
}


/*
** Implements the sqrt function, it has the peculiarity of returning an integer when the
** the argument is an integer.
** Since SQLite isn't strongly typed (almost untyped actually) this is a bit pedantic
*/
static void squareFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	i64 iVal = 0;
	double rVal = 0.0;

	assert(argc == 1);
	switch (sqlite3_value_type(argv[0]) ) {
	case SQLITE_INTEGER: {
		iVal = sqlite3_value_int64(argv[0]);
		sqlite3_result_int64(context, iVal * iVal);
		break;
	}
	case SQLITE_NULL: {
		sqlite3_result_null(context);
		break;
	}
	default: {
		rVal = sqlite3_value_double(argv[0]);
		sqlite3_result_double(context, rVal * rVal);
		break;
	}
	}
}


/*
** Wraps the pow math.h function
** When both the base and the exponent are integers the result should be integer
** (see sqrt just before this). Here the result is always double
*/
/* LMH 2007-03-25 Changed to use errno; no pre-checking for errors.  Also removes
   but that was present in the pre-checking that called sqlite3_result_error on
   a non-positive first argument, which is not always an error. */
static void powerFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	double r1 = 0.0;
	double r2 = 0.0;
	double val;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		r1 = sqlite3_value_double(argv[0]);
		r2 = sqlite3_value_double(argv[1]);
		errno = 0;
		val = pow(r1, r2);
		if (errno == 0) {
			sqlite3_result_double(context, val);
		} else {
			sqlite3_result_error(context, strerror(errno), errno);
		}
	}
}


/*
** atan2 wrapper
*/
static void atn2Func(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	double r1 = 0.0;
	double r2 = 0.0;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		r1 = sqlite3_value_double(argv[0]);
		r2 = sqlite3_value_double(argv[1]);
		sqlite3_result_double(context, atan2(r1, r2));
	}
}


/*
** Implementation of the sign() function
** return one of 3 possibilities +1,0 or -1 when the argument is respectively
** positive, 0 or negative.
** When the argument is NULL the result is also NULL (completly conventional)
*/
static void signFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	double rVal = 0.0;
	i64 iVal = 0;

	assert(argc == 1);
	switch (sqlite3_value_type(argv[0]) ) {
	case SQLITE_INTEGER: {
		iVal = sqlite3_value_int64(argv[0]);
		iVal = (iVal > 0) ? 1 : (iVal < 0) ? -1 : 0;
		sqlite3_result_int64(context, iVal);
		break;
	}
	case SQLITE_NULL: {
		sqlite3_result_null(context);
		break;
	}
	default: {
		/* 2nd change below. Line for abs was: if( rVal<0 ) rVal = rVal * -1.0;  */

		rVal = sqlite3_value_double(argv[0]);
		rVal = (rVal > 0) ? 1 : (rVal < 0) ? -1 : 0;
		sqlite3_result_double(context, rVal);
		break;
	}
	}
}


/*
** smallest integer value not less than argument
*/
static void ceilFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	double rVal = 0.0;
	// i64 iVal = 0;

	assert(argc == 1);
	switch (sqlite3_value_type(argv[0]) ) {
	case SQLITE_INTEGER: {
		i64 iVal = sqlite3_value_int64(argv[0]);
		sqlite3_result_int64(context, iVal);
		break;
	}
	case SQLITE_NULL: {
		sqlite3_result_null(context);
		break;
	}
	default: {
		rVal = sqlite3_value_double(argv[0]);
		sqlite3_result_int64(context, (i64)ceil(rVal));
		break;
	}
	}
}


/*
** largest integer value not greater than argument
*/
static void floorFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	double rVal = 0.0;
	i64 iVal = 0;

	assert(argc == 1);
	switch (sqlite3_value_type(argv[0]) ) {
	case SQLITE_INTEGER: {
		iVal = sqlite3_value_int64(argv[0]);
		sqlite3_result_int64(context, iVal);
		break;
	}
	case SQLITE_NULL: {
		sqlite3_result_null(context);
		break;
	}
	default: {
		rVal = sqlite3_value_double(argv[0]);
		sqlite3_result_int64(context, (i64)floor(rVal));
		break;
	}
	}
}


/*
** Given a string (s) in the first argument and an integer (n) in the second returns the
** string that constains s contatenated n times
*/
static void replicateFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	unsigned char * z;      /* input string */
	unsigned char * zo;     /* result string */
	i64 iCount;             /* times to repeat */
	i64 nLen;               /* length of the input string (no multibyte considerations) */
	i64 nTLen;              /* length of the result string (no multibyte considerations) */
	i64 i = 0;

	if (argc != 2 || SQLITE_NULL == sqlite3_value_type(argv[0]) )
		return;

	iCount = sqlite3_value_int64(argv[1]);

	if (iCount < 0) {
		sqlite3_result_error(context, "domain error", -1);
	}else{

		nLen = sqlite3_value_bytes(argv[0]);
		nTLen = nLen * iCount;
		z = (unsigned char *)sqlite3_malloc(nTLen + 1);
		zo = (unsigned char *)sqlite3_malloc(nLen + 1);
		if (!z || !zo) {
			sqlite3_result_error_nomem(context);
			if (z) sqlite3_free(z);
			if (zo) sqlite3_free(zo);
			return;
		}
		strcpy((char *)zo, (char *)sqlite3_value_text(argv[0]));

		for (i = 0; i < iCount; ++i) {
			strcpy((char *)(z + i * nLen), (char *)zo);
		}

		sqlite3_result_text(context, (char *)z, -1, SQLITE_TRANSIENT);
		sqlite3_free(z);
		sqlite3_free(zo);
	}
}


/*
** Some systems (win32 among others) don't have an isblank function, this will emulate it.
** This function is not UFT-8 safe since it only analyses a byte character.
*/
#ifndef HAVE_ISBLANK
int isblank(char c)
{
	return(' ' == c || '\t' == c);
}


#endif

static void properFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const unsigned char * z;    /* input string */
	unsigned char * zo;         /* output string */
	unsigned char * zt;         /* iterator */
	char r;
	int c = 1;

	assert(argc == 1);
	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}

	z = sqlite3_value_text(argv[0]);
	zo = (unsigned char *)sqlite3StrDup((char *)z);
	if (!zo) {
		sqlite3_result_error_nomem(context);
		return;
	}
	zt = zo;

	while ( (r = *(z++)) != 0) {
		if (isblank(r) ) {
			c = 1;
		}else{
			if (c == 1) {
				r = toupper(r);
			}else{
				r = tolower(r);
			}
			c = 0;
		}
		*(zt++) = r;
	}
	*zt = '\0';

	sqlite3_result_text(context, (char *)zo, -1, SQLITE_TRANSIENT);
	sqlite3_free(zo);
}


/*
** given an input string (s) and an integer (n) adds spaces at the begining of  s
** until it has a length of n characters.
** When s has a length >=n it's a NOP
** padl(NULL) = NULL
*/
static void padlFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	i64 ilen;           /* length to pad to */
	i64 zl;             /* length of the input string (UTF-8 chars) */
	int i = 0;
	const char * zi;    /* input string */
	char * zo;          /* output string */
	char * zt;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		zi = (char *)sqlite3_value_text(argv[0]);
		ilen = sqlite3_value_int64(argv[1]);
		/* check domain */
		if (ilen < 0) {
			sqlite3_result_error(context, "domain error", -1);
			return;
		}
		zl = sqlite3Utf8CharLen(zi, -1);
		if (zl >= ilen) {
			/* string is longer than the requested pad length, return the same string (dup it) */
			zo = sqlite3StrDup(zi);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		}else{
			zo = (char *)sqlite3_malloc(strlen(zi) + ilen - zl + 1);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			zt = zo;
			for (i = 1; i + zl <= ilen; ++i) {
				*(zt++) = ' ';
			}
			/* no need to take UTF-8 into consideration here */
			strcpy(zt, zi);
		}
		sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		sqlite3_free(zo);
	}
}


/*
** given an input string (s) and an integer (n) appends spaces at the end of  s
** until it has a length of n characters.
** When s has a length >=n it's a NOP
** padl(NULL) = NULL
*/
static void padrFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	i64 ilen;           /* length to pad to */
	i64 zl;             /* length of the input string (UTF-8 chars) */
	i64 zll;            /* length of the input string (bytes) */
	int i = 0;
	const char * zi;    /* input string */
	char * zo;          /* output string */
	char * zt;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		zi = (char *)sqlite3_value_text(argv[0]);
		ilen = sqlite3_value_int64(argv[1]);
		/* check domain */
		if (ilen < 0) {
			sqlite3_result_error(context, "domain error", -1);
			return;
		}
		zl = sqlite3Utf8CharLen(zi, -1);
		if (zl >= ilen) {
			/* string is longer than the requested pad length, return the same string (dup it) */
			zo = sqlite3StrDup(zi);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		}else{
			zll = strlen(zi);
			zo = (char *)sqlite3_malloc(zll + ilen - zl + 1);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			zt = strcpy(zo, zi) + zll;
			for (i = 1; i + zl <= ilen; ++i) {
				*(zt++) = ' ';
			}
			*zt = '\0';
		}
		sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		sqlite3_free(zo);
	}
}


/*
** given an input string (s) and an integer (n) appends spaces at the end of  s
** and adds spaces at the begining of s until it has a length of n characters.
** Tries to add has many characters at the left as at the right.
** When s has a length >=n it's a NOP
** padl(NULL) = NULL
*/
static void padcFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	i64 ilen;           /* length to pad to */
	i64 zl;             /* length of the input string (UTF-8 chars) */
	i64 zll;            /* length of the input string (bytes) */
	int i = 0;
	const char * zi;    /* input string */
	char * zo;          /* output string */
	char * zt;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		zi = (char *)sqlite3_value_text(argv[0]);
		ilen = sqlite3_value_int64(argv[1]);
		/* check domain */
		if (ilen < 0) {
			sqlite3_result_error(context, "domain error", -1);
			return;
		}
		zl = sqlite3Utf8CharLen(zi, -1);
		if (zl >= ilen) {
			/* string is longer than the requested pad length, return the same string (dup it) */
			zo = sqlite3StrDup(zi);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		}else{
			zll = strlen(zi);
			zo = (char *)sqlite3_malloc(zll + ilen - zl + 1);
			if (!zo) {
				sqlite3_result_error_nomem(context);
				return;
			}
			zt = zo;
			for (i = 1; 2 * i + zl <= ilen; ++i) {
				*(zt++) = ' ';
			}
			strcpy(zt, zi);
			zt += zll;
			for (; i + zl <= ilen; ++i) {
				*(zt++) = ' ';
			}
			*zt = '\0';
		}
		sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		sqlite3_free(zo);
	}
}


/*
** given 2 string (s1,s2) returns the string s1 with the characters NOT in s2 removed
** assumes strings are UTF-8 encoded
*/
static void strfilterFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * zi1;       /* first parameter string (searched string) */
	const char * zi2;       /* second parameter string (vcontains valid characters) */
	const char * z1;
	const char * z21;
	const char * z22;
	char * zo;            /* output string */
	char * zot;
	int c1 = 0;
	int c2 = 0;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
		sqlite3_result_null(context);
	}else{
		zi1 = (char *)sqlite3_value_text(argv[0]);
		zi2 = (char *)sqlite3_value_text(argv[1]);
		/*
		** maybe I could allocate less, but that would imply 2 passes, rather waste
		** (possibly) some memory
		*/
		zo = (char *)sqlite3_malloc(strlen(zi1) + 1);
		if (!zo) {
			sqlite3_result_error_nomem(context);
			return;
		}
		zot = zo;
		z1 = zi1;
		while ( (c1 = sqliteCharVal((unsigned char *)z1)) != 0) {
			z21 = zi2;
			while ( (c2 = sqliteCharVal((unsigned char *)z21)) != 0 && c2 != c1) {
				sqliteNextChar(z21);
			}
			if (c2 != 0) {
				z22 = z21;
				sqliteNextChar(z22);
				strncpy(zot, z21, z22 - z21);
				zot += z22 - z21;
			}
			sqliteNextChar(z1);
		}
		*zot = '\0';

		sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
		sqlite3_free(zo);
	}
}


/*
** Given a string z1, retutns the (0 based) index of it's first occurence
** in z2 after the first s characters.
** Returns -1 when there isn't a match.
** updates p to point to the character where the match occured.
** This is an auxiliary function.
*/
static int _substr(const char * z1, const char * z2, int s, const char ** p)
{
	int c = 0;
	int rVal = -1;
	const char * zt1;
	const char * zt2;
	int c1, c2;

	if ('\0' == *z1) {
		return -1;
	}

	while ( (sqliteCharVal((unsigned char *)z2) != 0) && (c++) < s) {
		sqliteNextChar(z2);
	}

	c = 0;
	while ( (sqliteCharVal((unsigned char *)z2)) != 0) {
		zt1 = z1;
		zt2 = z2;

		do {
			c1 = sqliteCharVal((unsigned char *)zt1);
			c2 = sqliteCharVal((unsigned char *)zt2);
			sqliteNextChar(zt1);
			sqliteNextChar(zt2);
		} while (c1 == c2 && c1 != 0 && c2 != 0);

		if (c1 == 0) {
			rVal = c;
			break;
		}

		sqliteNextChar(z2);
		++c;
	}
	if (p) {
		*p = z2;
	}
	return rVal >= 0 ? rVal + s : rVal;
}


/*
** given 2 input strings (s1,s2) and an integer (n) searches from the nth character
** for the string s1. Returns the position where the match occured.
** Characters are counted from 1.
** 0 is returned when no match occurs.
*/

static void charindexFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const u8 * z1;          /* s1 string */
	u8 * z2;                /* s2 string */
	int s = 0;
	int rVal = 0;

	assert(argc == 3 || argc == 2);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) || SQLITE_NULL == sqlite3_value_type(argv[1])) {
		sqlite3_result_null(context);
		return;
	}

	z1 = sqlite3_value_text(argv[0]);
	if (z1 == 0) return;
	z2 = (u8 *)sqlite3_value_text(argv[1]);
	if (argc == 3) {
		s = sqlite3_value_int(argv[2]) - 1;
		if (s < 0) {
			s = 0;
		}
	}else{
		s = 0;
	}

	rVal = _substr((char *)z1, (char *)z2, s, NULL);
	sqlite3_result_int(context, rVal + 1);
}


/*
** given a string (s) and an integer (n) returns the n leftmost (UTF-8) characters
** if the string has a length<=n or is NULL this function is NOP
*/
static void leftFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	int c = 0;
	int cc = 0;
	int l = 0;
	const unsigned char * z;    /* input string */
	const unsigned char * zt;
	unsigned char * rz;         /* output string */

	assert(argc == 2);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) || SQLITE_NULL == sqlite3_value_type(argv[1])) {
		sqlite3_result_null(context);
		return;
	}

	z = sqlite3_value_text(argv[0]);
	l = sqlite3_value_int(argv[1]);
	zt = z;

	while (sqliteCharVal(zt) && c++ < l)
		sqliteNextChar(zt);

	cc = zt - z;

	rz = (unsigned char *)sqlite3_malloc(zt - z + 1);
	if (!rz) {
		sqlite3_result_error_nomem(context);
		return;
	}
	strncpy((char *)rz, (char *)z, zt - z);
	*(rz + cc) = '\0';
	sqlite3_result_text(context, (char *)rz, -1, SQLITE_TRANSIENT);
	sqlite3_free(rz);
}


/*
** given a string (s) and an integer (n) returns the n rightmost (UTF-8) characters
** if the string has a length<=n or is NULL this function is NOP
*/
static void rightFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	int l = 0;
	int c = 0;
	int cc = 0;
	const char * z;
	const char * zt;
	const char * ze;
	char * rz;

	assert(argc == 2);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) || SQLITE_NULL == sqlite3_value_type(argv[1])) {
		sqlite3_result_null(context);
		return;
	}

	z = (char *)sqlite3_value_text(argv[0]);
	l = sqlite3_value_int(argv[1]);
	zt = z;

	while (sqliteCharVal((unsigned char *)zt) != 0) {
		sqliteNextChar(zt);
		++c;
	}

	ze = zt;
	zt = z;

	cc = c - l;
	if (cc < 0)
		cc = 0;

	while (cc-- > 0) {
		sqliteNextChar(zt);
	}

	rz = (char *)sqlite3_malloc(ze - zt + 1);
	if (!rz) {
		sqlite3_result_error_nomem(context);
		return;
	}
	strcpy((char *)rz, (char *)(zt));
	sqlite3_result_text(context, (char *)rz, -1, SQLITE_TRANSIENT);
	sqlite3_free(rz);
}


#ifndef HAVE_TRIM
/*
** removes the whitespaces at the begining of a string.
*/
const char * ltrim(const char * s)
{
	while (*s == ' ')
		++s;
	return s;
}


/*
** removes the whitespaces at the end of a string.
** !mutates the input string!
*/
void rtrim(char * s)
{
	char * ss = s + strlen(s) - 1;

	while (ss >= s && *ss == ' ')
		--ss;
	*(ss + 1) = '\0';
}


/*
**  Removes the whitespace at the begining of a string
*/
static void ltrimFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * z;

	assert(argc == 1);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}
	z = sqlite3_value_text(argv[0]);
	sqlite3_result_text(context, ltrim(z), -1, SQLITE_TRANSIENT);
}


/*
**  Removes the whitespace at the end of a string
*/
static void rtrimFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * z;
	char * rz;

	/* try not to change data in argv */

	assert(argc == 1);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}
	z = sqlite3_value_text(argv[0]);
	rz = sqlite3StrDup(z);
	rtrim(rz);
	sqlite3_result_text(context, rz, -1, SQLITE_TRANSIENT);
	sqlite3_free(rz);
}


/*
**  Removes the whitespace at the begining and end of a string
*/
static void trimFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * z;
	char * rz;

	/* try not to change data in argv */

	assert(argc == 1);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}
	z = sqlite3_value_text(argv[0]);
	rz = sqlite3StrDup(z);
	rtrim(rz);
	sqlite3_result_text(context, ltrim(rz), -1, SQLITE_TRANSIENT);
	sqlite3_free(rz);
}


#endif

/*
** given a pointer to a string s1, the length of that string (l1), a new string (s2)
** and it's length (l2) appends s2 to s1.
** All lengths in bytes.
** This is just an auxiliary function
*/
// static void _append(char **s1, int l1, const char *s2, int l2){
//   *s1 = realloc(*s1, (l1+l2+1)*sizeof(char));
//   strncpy((*s1)+l1, s2, l2);
//   *(*(s1)+l1+l2) = '\0';
// }

#ifndef HAVE_TRIM

/*
** given strings s, s1 and s2 replaces occurrences of s1 in s by s2
*/
static void replaceFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * z1;    /* string s (first parameter) */
	const char * z2;    /* string s1 (second parameter) string to look for */
	const char * z3;    /* string s2 (third parameter) string to replace occurrences of s1 with */
	int lz1;
	int lz2;
	int lz3;
	int lzo = 0;
	char * zo = 0;
	int ret = 0;
	const char * zt1;
	const char * zt2;

	assert(3 == argc);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}

	z1 = sqlite3_value_text(argv[0]);
	z2 = sqlite3_value_text(argv[1]);
	z3 = sqlite3_value_text(argv[2]);
	/* handle possible null values */
	if (0 == z2) {
		z2 = "";
	}
	if (0 == z3) {
		z3 = "";
	}

	lz1 = strlen(z1);
	lz2 = strlen(z2);
	lz3 = strlen(z3);

#  if 0
	/* special case when z2 is empty (or null) nothing will be changed */
	if (0 == lz2) {
		sqlite3_result_text(context, z1, -1, SQLITE_TRANSIENT);
		return;
	}
#  endif

	zt1 = z1;
	zt2 = z1;

	while (1) {
		ret = _substr(z2, zt1, 0, &zt2);

		if (ret < 0)
			break;

		_append(&zo, lzo, zt1, zt2 - zt1);
		lzo += zt2 - zt1;
		_append(&zo, lzo, z3, lz3);
		lzo += lz3;

		zt1 = zt2 + lz2;
	}
	_append(&zo, lzo, zt1, lz1 - (zt1 - z1));
	sqlite3_result_text(context, zo, -1, SQLITE_TRANSIENT);
	sqlite3_free(zo);
}


#endif

/*
** given a string returns the same string but with the characters in reverse order
*/
static void reverseFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	const char * z;
	const char * zt;
	char * rz;
	char * rzt;
	int l = 0;
	int i = 0;

	assert(1 == argc);

	if (SQLITE_NULL == sqlite3_value_type(argv[0]) ) {
		sqlite3_result_null(context);
		return;
	}
	z = (char *)sqlite3_value_text(argv[0]);
	l = strlen(z);
	rz = (char *)sqlite3_malloc(l + 1);
	if (!rz) {
		sqlite3_result_error_nomem(context);
		return;
	}
	rzt = rz + l;
	*(rzt--) = '\0';

	zt = z;
	while (sqliteCharVal((unsigned char *)zt) != 0) {
		z = zt;
		sqliteNextChar(zt);
		for (i = 1; zt - i >= z; ++i) {
			*(rzt--) = *(zt - i);
		}
	}

	sqlite3_result_text(context, rz, -1, SQLITE_TRANSIENT);
	sqlite3_free(rz);
}


/*
** An instance of the following structure holds the context of a
** stdev() or variance() aggregate computation.
** implementaion of http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Algorithm_II
** less prone to rounding errors
*/
typedef struct StdevCtx StdevCtx;
struct StdevCtx
{
	double rM;
	double rS;
	i64 cnt;        /* number of elements */
};

/*
** An instance of the following structure holds the context of a
** mode() or median() aggregate computation.
** Depends on structures defined in map.c (see map & map)
** These aggregate functions only work for integers and floats although
** they could be made to work for strings. This is usually considered meaningless.
** Only usuall order (for median), no use of collation functions (would this even make sense?)
*/
typedef struct ModeCtx ModeCtx;
struct ModeCtx
{
	i64 riM;            /* integer value found so far */
	double rdM;         /* double value found so far */
	i64 cnt;            /* number of elements so far */
	double pcnt;        /* number of elements smaller than a percentile */
	i64 mcnt;           /* maximum number of occurrences (for mode) */
	i64 mn;             /* number of occurrences (for mode and percentiles) */
	i64 is_double;      /* whether the computation is being done for doubles (>0) or integers (=0) */
	map * m;            /* map structure used for the computation */
	int done;           /* whether the answer has been found */
};

/*
** called for each value received during a calculation of stdev or variance
*/
static void varianceStep(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	StdevCtx * p;

	double delta;
	double x;

	assert(argc == 1);
	p = (StdevCtx *)sqlite3_aggregate_context(context, sizeof(*p));
	/* only consider non-null values */
	if (SQLITE_NULL != sqlite3_value_numeric_type(argv[0]) ) {
		p->cnt++;
		x = sqlite3_value_double(argv[0]);
		delta = (x - p->rM);
		p->rM += delta / p->cnt;
		p->rS += delta * (x - p->rM);
	}
}


/*
** called for each value received during a calculation of mode of median
*/
static void modeStep(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	ModeCtx * p;
	i64 xi = 0;
	double xd = 0.0;
	i64 * iptr;
	double * dptr;
	int type;

	assert(argc == 1);
	type = sqlite3_value_numeric_type(argv[0]);

	if (type == SQLITE_NULL)
		return;

	p = (ModeCtx *)sqlite3_aggregate_context(context, sizeof(*p));

	if (0 == (p->m) ) {
		p->m = (map *)calloc(1, sizeof(map));
		if (type == SQLITE_INTEGER) {
			/* map will be used for integers */
			*(p->m) = map_make(int_cmp);
			p->is_double = 0;
		}else{
			p->is_double = 1;
			/* map will be used for doubles */
			*(p->m) = map_make(double_cmp);
		}
	}

	++(p->cnt);

	if (0 == p->is_double) {
		xi = sqlite3_value_int64(argv[0]);
		iptr = (i64 *)calloc(1, sizeof(i64));
		*iptr = xi;
		map_insert(p->m, iptr);
	}else{
		xd = sqlite3_value_double(argv[0]);
		dptr = (double *)calloc(1, sizeof(double));
		*dptr = xd;
		map_insert(p->m, dptr);
	}
}


/*
**  Auxiliary function that iterates all elements in a map and finds the mode
**  (most frequent value)
*/
static void modeIterate(void * e, i64 c, void * pp)
{
	i64 ei;
	double ed;
	ModeCtx * p = (ModeCtx *)pp;

	if (0 == p->is_double) {
		ei = *(int *)(e);

		if (p->mcnt == c) {
			++p->mn;
		}else if (p->mcnt < c) {
			p->riM = ei;
			p->mcnt = c;
			p->mn = 1;
		}
	}else{
		ed = *(double *)(e);

		if (p->mcnt == c) {
			++p->mn;
		}else if (p->mcnt < c) {
			p->rdM = ed;
			p->mcnt = c;
			p->mn = 1;
		}
	}
}


/*
**  Auxiliary function that iterates all elements in a map and finds the median
**  (the value such that the number of elements smaller is equal the the number of
**  elements larger)
*/
static void medianIterate(void * e, i64 c, void * pp)
{
	i64 ei;
	double ed;
	double iL;
	double iR;
	int il;
	int ir;
	ModeCtx * p = (ModeCtx *)pp;

	if (p->done > 0)
		return;

	iL = p->pcnt;
	iR = p->cnt - p->pcnt;
	il = p->mcnt + c;
	ir = p->cnt - p->mcnt;

	if (il >= iL) {
		if (ir >= iR) {
			++p->mn;
			if (0 == p->is_double) {
				ei = *(int *)(e);
				p->riM += ei;
			}else{
				ed = *(double *)(e);
				p->rdM += ed;
			}
		}else{
			p->done = 1;
		}
	}
	p->mcnt += c;
}


/*
** Returns the mode value
*/
static void modeFinalize(sqlite3_context * context)
{
	ModeCtx * p;

	p = (ModeCtx *)sqlite3_aggregate_context(context, 0);
	if (p && p->m) {
		map_iterate(p->m, modeIterate, p);
		map_destroy(p->m);
		free(p->m);

		if (1 == p->mn) {
			if (0 == p->is_double)
				sqlite3_result_int64(context, p->riM);
			else
				sqlite3_result_double(context, p->rdM);
		}
	}
}


/*
** auxiliary function for percentiles
*/
static void _medianFinalize(sqlite3_context * context)
{
	ModeCtx * p;

	p = (ModeCtx *)sqlite3_aggregate_context(context, 0);
	if (p && p->m) {
		p->done = 0;
		map_iterate(p->m, medianIterate, p);
		map_destroy(p->m);
		free(p->m);

		if (0 == p->is_double)
			if (1 == p->mn)
				sqlite3_result_int64(context, p->riM);
			else
				sqlite3_result_double(context, p->riM * 1.0 / p->mn);
		else
			sqlite3_result_double(context, p->rdM / p->mn);
	}
}


/*
** Returns the median value
*/
static void medianFinalize(sqlite3_context * context)
{
	ModeCtx * p;

	p = (ModeCtx *)sqlite3_aggregate_context(context, 0);
	if (p != 0) {
		p->pcnt = (p->cnt) / 2.0;
		_medianFinalize(context);
	}
}


/*
** Returns the lower_quartile value
*/
static void lower_quartileFinalize(sqlite3_context * context)
{
	ModeCtx * p;

	p = (ModeCtx *)sqlite3_aggregate_context(context, 0);
	if (p != 0) {
		p->pcnt = (p->cnt) / 4.0;
		_medianFinalize(context);
	}
}


/*
** Returns the upper_quartile value
*/
static void upper_quartileFinalize(sqlite3_context * context)
{
	ModeCtx * p;

	p = (ModeCtx *)sqlite3_aggregate_context(context, 0);
	if (p != 0) {
		p->pcnt = (p->cnt) * 3 / 4.0;
		_medianFinalize(context);
	}
}


/*
** Returns the stdev value
*/
static void stdevFinalize(sqlite3_context * context)
{
	StdevCtx * p;

	p = (StdevCtx *)sqlite3_aggregate_context(context, 0);
	if (p && p->cnt > 1) {
		sqlite3_result_double(context, sqrt(p->rS / (p->cnt - 1)));
	}else{
		sqlite3_result_double(context, 0.0);
	}
}


/*
** Returns the variance value
*/
static void varianceFinalize(sqlite3_context * context)
{
	StdevCtx * p;

	p = (StdevCtx *)sqlite3_aggregate_context(context, 0);
	if (p && p->cnt > 1) {
		sqlite3_result_double(context, p->rS / (p->cnt - 1));
	}else{
		sqlite3_result_double(context, 0.0);
	}
}


#ifdef SQLITE_SOUNDEX

/* relicoder factored code */
/*
** Calculates the soundex value of a string
*/

static void soundex(const u8 * zIn, char * zResult)
{
	int i, j;
	static const unsigned char iCode[] = {
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 0, 1, 2, 3, 0, 1, 2, 0, 0, 2, 2, 4, 5, 5, 0,
		1, 2, 6, 2, 3, 0, 1, 0, 2, 0, 2, 0, 0, 0, 0, 0,
		0, 0, 1, 2, 3, 0, 1, 2, 0, 0, 2, 2, 4, 5, 5, 0,
		1, 2, 6, 2, 3, 0, 1, 0, 2, 0, 2, 0, 0, 0, 0, 0,
	};

	for (i = 0; zIn[i] && !isalpha(zIn[i]); i++) {
	}
	if (zIn[i]) {
		zResult[0] = toupper(zIn[i]);
		for (j = 1; j < 4 && zIn[i]; i++) {
			int code = iCode[zIn[i] & 0x7f];
			if (code > 0) {
				zResult[j++] = code + '0';
			}
		}
		while (j < 4) {
			zResult[j++] = '0';
		}
		zResult[j] = 0;
	}else{
		strcpy(zResult, "?000");
	}
}


/*
** computes the number of different characters between the soundex value fo 2 strings
*/
static void differenceFunc(sqlite3_context * context, int argc, sqlite3_value ** argv)
{
	char zResult1[8];
	char zResult2[8];
	char * zR1 = zResult1;
	char * zR2 = zResult2;
	int rVal = 0;
	int i = 0;
	const u8 * zIn1;
	const u8 * zIn2;

	assert(argc == 2);

	if (sqlite3_value_type(argv[0]) == SQLITE_NULL || sqlite3_value_type(argv[1]) == SQLITE_NULL) {
		sqlite3_result_null(context);
		return;
	}

	zIn1 = (u8 *)sqlite3_value_text(argv[0]);
	zIn2 = (u8 *)sqlite3_value_text(argv[1]);

	soundex(zIn1, zR1);
	soundex(zIn2, zR2);

	for (i = 0; i < 4; ++i) {
		if (sqliteCharVal((unsigned char *)zR1) == sqliteCharVal((unsigned char *)zR2) )
			++rVal;
		sqliteNextChar(zR1);
		sqliteNextChar(zR2);
	}
	sqlite3_result_int(context, rVal);
}


#endif

/*
** This function registered all of the above C functions as SQL
** functions.  This should be the only routine in this file with
** external linkage.
*/
int RegisterExtensionFunctions(sqlite3 * db)
{
	static const struct FuncDef
	{
		char * zName;
		signed char nArg;
		u8 argType;         /* 0: none.  1: db  2: (-1) */
		u8 eTextRep;        /* 1: UTF-16.  0: UTF-8 */
		u8 needCollSeq;
		void (* xFunc)(sqlite3_context *, int, sqlite3_value **);
	} aFuncs[] = {
		/* math.h */
		{ "acos",		1,		   0,		  SQLITE_UTF8,		   0,		  acosFunc																					},
		{ "asin",		1,		   0,		  SQLITE_UTF8,		   0,		  asinFunc																					},
		{ "atan",		1,		   0,		  SQLITE_UTF8,		   0,		  atanFunc																					},
		{ "atn2",		2,		   0,		  SQLITE_UTF8,		   0,		  atn2Func																					},
		/* XXX alias */
		{ "atan2",		2,		   0,		  SQLITE_UTF8,		   0,		  atn2Func																					},
		{ "acosh",		1,		   0,		  SQLITE_UTF8,		   0,		  acoshFunc																					},
		{ "asinh",		1,		   0,		  SQLITE_UTF8,		   0,		  asinhFunc																					},
		{ "atanh",		1,		   0,		  SQLITE_UTF8,		   0,		  atanhFunc																					},

		{ "difference", 2,		   0,		  SQLITE_UTF8,		   0,		  differenceFunc																			},
		{ "degrees",	1,		   0,		  SQLITE_UTF8,		   0,		  rad2degFunc																				},
		{ "radians",	1,		   0,		  SQLITE_UTF8,		   0,		  deg2radFunc																				},

		{ "cos",		1,		   0,		  SQLITE_UTF8,		   0,		  cosFunc																					},
		{ "sin",		1,		   0,		  SQLITE_UTF8,		   0,		  sinFunc																					},
		{ "tan",		1,		   0,		  SQLITE_UTF8,		   0,		  tanFunc																					},
		{ "cot",		1,		   0,		  SQLITE_UTF8,		   0,		  cotFunc																					},
		{ "cosh",		1,		   0,		  SQLITE_UTF8,		   0,		  coshFunc																					},
		{ "sinh",		1,		   0,		  SQLITE_UTF8,		   0,		  sinhFunc																					},
		{ "tanh",		1,		   0,		  SQLITE_UTF8,		   0,		  tanhFunc																					},
		{ "coth",		1,		   0,		  SQLITE_UTF8,		   0,		  cothFunc																					},

		{ "exp",		1,		   0,		  SQLITE_UTF8,		   0,		  expFunc																					},
		{ "log",		1,		   0,		  SQLITE_UTF8,		   0,		  logFunc																					},
		{ "log10",		1,		   0,		  SQLITE_UTF8,		   0,		  log10Func																					},
		{ "power",		2,		   0,		  SQLITE_UTF8,		   0,		  powerFunc																					},
		{ "sign",		1,		   0,		  SQLITE_UTF8,		   0,		  signFunc																					},
		{ "sqrt",		1,		   0,		  SQLITE_UTF8,		   0,		  sqrtFunc																					},
		{ "square",		1,		   0,		  SQLITE_UTF8,		   0,		  squareFunc																				},

		{ "ceil",		1,		   0,		  SQLITE_UTF8,		   0,		  ceilFunc																					},
		{ "floor",		1,		   0,		  SQLITE_UTF8,		   0,		  floorFunc																					},

		{ "pi",			0,		   0,		  SQLITE_UTF8,		   1,		  piFunc																					},


		/* string */
		{ "replicate",	2,		   0,		  SQLITE_UTF8,		   0,		  replicateFunc																				},
		{ "charindex",	2,		   0,		  SQLITE_UTF8,		   0,		  charindexFunc																				},
		{ "charindex",	3,		   0,		  SQLITE_UTF8,		   0,		  charindexFunc																				},
		{ "leftstr",	2,		   0,		  SQLITE_UTF8,		   0,		  leftFunc																					},
		{ "rightstr",	2,		   0,		  SQLITE_UTF8,		   0,		  rightFunc																					},
#ifndef HAVE_TRIM
		{ "ltrim",		1,		   0,		  SQLITE_UTF8,		   0,		  ltrimFunc																					},
		{ "rtrim",		1,		   0,		  SQLITE_UTF8,		   0,		  rtrimFunc																					},
		{ "trim",		1,		   0,		  SQLITE_UTF8,		   0,		  trimFunc																					},
		{ "replace",	3,		   0,		  SQLITE_UTF8,		   0,		  replaceFunc																				},
#endif
		{ "reverse",	1,		   0,		  SQLITE_UTF8,		   0,		  reverseFunc																				},
		{ "proper",		1,		   0,		  SQLITE_UTF8,		   0,		  properFunc																				},
		{ "padl",		2,		   0,		  SQLITE_UTF8,		   0,		  padlFunc																					},
		{ "padr",		2,		   0,		  SQLITE_UTF8,		   0,		  padrFunc																					},
		{ "padc",		2,		   0,		  SQLITE_UTF8,		   0,		  padcFunc																					},
		{ "strfilter",	2,		   0,		  SQLITE_UTF8,		   0,		  strfilterFunc																				},

	};
	/* Aggregate functions */
	static const struct FuncDefAgg
	{
		char * zName;
		signed char nArg;
		u8 argType;
		u8 needCollSeq;
		void (* xStep)(sqlite3_context *, int, sqlite3_value **);
		void (* xFinalize)(sqlite3_context *);
	} aAggs[] = {
		{ "stdev",			1,			 0,			  0,			varianceStep,			  stdevFinalize															},
		{ "variance",		1,			 0,			  0,			varianceStep,			  varianceFinalize														},
		{ "mode",			1,			 0,			  0,			modeStep,				  modeFinalize															},
		{ "median",			1,			 0,			  0,			modeStep,				  medianFinalize														},
		{ "lower_quartile", 1,			 0,			  0,			modeStep,				  lower_quartileFinalize												},
		{ "upper_quartile", 1,			 0,			  0,			modeStep,				  upper_quartileFinalize												},
	};
	int i;

	for (i = 0; i < sizeof(aFuncs) / sizeof(aFuncs[0]); i++) {
		void * pArg = 0;
		switch (aFuncs[i].argType) {
		case 1: pArg = db; break;
		case 2: pArg = (void *)(-1); break;
		}
		//sqlite3CreateFunc
		/* LMH no error checking */
		sqlite3_create_function(db, aFuncs[i].zName, aFuncs[i].nArg,
			aFuncs[i].eTextRep, pArg, aFuncs[i].xFunc, 0, 0);
#if 0
		if (aFuncs[i].needCollSeq) {
			struct FuncDef * pFunc = sqlite3FindFunction(db, aFuncs[i].zName,
				strlen(aFuncs[i].zName), aFuncs[i].nArg, aFuncs[i].eTextRep, 0);
			if (pFunc && aFuncs[i].needCollSeq) {
				pFunc->needCollSeq = 1;
			}
		}
#endif
	}

	for (i = 0; i < sizeof(aAggs) / sizeof(aAggs[0]); i++) {
		void * pArg = 0;
		switch (aAggs[i].argType) {
		case 1: pArg = db; break;
		case 2: pArg = (void *)(-1); break;
		}
		//sqlite3CreateFunc
		/* LMH no error checking */
		sqlite3_create_function(db, aAggs[i].zName, aAggs[i].nArg, SQLITE_UTF8,
			pArg, 0, aAggs[i].xStep, aAggs[i].xFinalize);
#if 0
		if (aAggs[i].needCollSeq) {
			struct FuncDefAgg * pFunc = sqlite3FindFunction(db, aAggs[i].zName,
				strlen(aAggs[i].zName), aAggs[i].nArg, SQLITE_UTF8, 0);
			if (pFunc && aAggs[i].needCollSeq) {
				pFunc->needCollSeq = 1;
			}
		}
#endif
	}
	return 0;
}


#ifdef COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE
int sqlite3_extension_init(
                           sqlite3 * db, char ** pzErrMsg, const sqlite3_api_routines * pApi)
{
	SQLITE_EXTENSION_INIT2(pApi);
	RegisterExtensionFunctions(db);
	sqlite3_my_extension_init(db, pzErrMsg, pApi);
	return 0;
}


#endif /* COMPILE_SQLITE_EXTENSIONS_AS_LOADABLE_MODULE */

map map_make(cmp_func cmp)
{
	map r;

	r.cmp = cmp;
	r.base = 0;

	return r;
}


void * xcalloc(size_t nmemb, size_t size, char * s)
{
	void * ret = calloc(nmemb, size);

	return ret;
}


void xfree(void * p)
{
	free(p);
}


void node_insert(node ** n, cmp_func cmp, void * e)
{
	int c;
	node * nn;

	if (*n == 0) {
		nn = (node *)xcalloc(1, sizeof(node), "for node");
		nn->data = e;
		nn->count = 1;
		*n = nn;
	}else{
		c = cmp((*n)->data, e);
		if (0 == c) {
			++((*n)->count);
			xfree(e);
		}else if (c > 0) {
			/* put it right here */
			node_insert(&((*n)->l), cmp, e);
		}else{
			node_insert(&((*n)->r), cmp, e);
		}
	}
}


void map_insert(map * m, void * e)
{
	node_insert(&(m->base), m->cmp, e);
}


void node_iterate(node * n, map_iterator iter, void * p)
{
	if (n) {
		if (n->l)
			node_iterate(n->l, iter, p);
		iter(n->data, n->count, p);
		if (n->r)
			node_iterate(n->r, iter, p);
	}
}


void map_iterate(map * m, map_iterator iter, void * p)
{
	node_iterate(m->base, iter, p);
}


void node_destroy(node * n)
{
	if (0 != n) {
		xfree(n->data);
		if (n->l)
			node_destroy(n->l);
		if (n->r)
			node_destroy(n->r);

		xfree(n);
	}
}


void map_destroy(map * m)
{
	node_destroy(m->base);
}


int int_cmp(const void * a, const void * b)
{
	int64_t aa = *(int64_t *)(a);
	int64_t bb = *(int64_t *)(b);

	/* printf("cmp %d <=> %d\n",aa,bb); */
	if (aa == bb)
		return 0;
	else if (aa < bb)
		return -1;
	else
		return 1;
}


int double_cmp(const void * a, const void * b)
{
	double aa = *(double *)(a);
	double bb = *(double *)(b);

	/* printf("cmp %d <=> %d\n",aa,bb); */
	if (aa == bb)
		return 0;
	else if (aa < bb)
		return -1;
	else
		return 1;
}


void print_elem(void * e, int64_t c, void * p)
{
	int ee = *(int *)(e);

	printf("%d => %lld\n", ee, c);
}


