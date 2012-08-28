/*
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

#include <sqlite3ext.h>
SQLITE_EXTENSION_INIT1

/*
** Demo: The half() SQL function returns half of its input value.
*/
static void halfFunc(
                     sqlite3_context * context,
                     int argc,
                     sqlite3_value ** argv
                     )
{
	sqlite3_result_double(context, 0.5 * sqlite3_value_double(argv[0]));
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
	x12;
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


// add an init function to make python extension happy
void init_vt_sqlite3_ext()
{
}

/* SQLite invokes this routine once when it loads the extension.
** Create new functions, collating sequences, and virtual table
** modules here.  This is usually the only exported symbol in
** the shared library.
*/
int sqlite3_extension_init(
                           sqlite3 * db,
                           char ** pzErrMsg,
                           const sqlite3_api_routines * pApi
                           )
{
	SQLITE_EXTENSION_INIT2(pApi)
	sqlite3_create_function(db, "half", 1, SQLITE_ANY, 0, halfFunc, 0, 0);
	//http://www.sqlite.org/c3ref/create_function.html
	//The first parameter is the database connection to which the SQL function is to be added
	//The second parameter is the name of the SQL function to be created or redefined.
	//The third parameter (nArg) is the number of arguments that the SQL function or aggregate takes.
	//The fourth parameter, eTextRep, specifies what text encoding this SQL function prefers for its parameters.
	//The fifth parameter is an arbitrary pointer.
	//The sixth, seventh and eighth parameters, xFunc, xStep and xFinal, are pointers to C-language functions that implement the SQL function or aggregate.
	sqlite3_create_function(db, "HWE_exact", -1, SQLITE_ANY, 0, hwe_exact, 0, 0);
	sqlite3_create_function(db, "Fisher_exact", 4, SQLITE_ANY, 0, fisher_exact, 0, 0);
	return 0;
}


