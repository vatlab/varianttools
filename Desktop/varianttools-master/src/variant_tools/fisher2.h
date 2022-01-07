/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 1998-2001   The R Development Core Team.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*  fisher2.c: this is fexact.c, part of the R
    (http://cran.r-project.org) distribution
    (package ctest).
    I have just made a few changes to make it easier
    to use it from the rest of my code (adjusted10.cpp and
    adj-fast.cpp).
     Ramón Díaz-Uriarte (rdiaz@cnio.es). May 2003.
    http://ligarto.org/rdiaz/
 */

#ifndef GUARD_fisher2
#define GUARD_fisher2

/*Constant.h*/
#ifndef M_PI
#  define M_PI 3.141592653589793238462643383279502884197169399375
#endif

#define PI             M_PI
#define SINGLE_EPS     FLT_EPSILON
#define SINGLE_BASE    FLT_RADIX
#define SINGLE_XMIN    FLT_MIN
#define SINGLE_XMAX    FLT_MAX
#define DOUBLE_DIGITS  DBL_MANT_DIG
#define DOUBLE_EPS     DBL_EPSILON
#define DOUBLE_XMAX    DBL_MAX
#define DOUBLE_XMIN    DBL_MIN

/*Memory.h
#ifdef  __cplusplus
extern "C" {
#endif

char * vmaxget(void);

void    vmaxset(char *);

void    R_gc(void);

char * R_alloc(long, int);

char * S_alloc(long, int);

char * S_realloc(char *, long, long, int);

#ifdef  __cplusplus
}
#endif
*/
/*boolean.h*/
#undef FALSE
#undef TRUE

#ifdef  __cplusplus
extern "C" {
#endif
typedef enum { FALSE = 0, TRUE /*, MAYBE */
} Rboolean;

#ifdef  __cplusplus
}
#endif

/*fisher2.h*/
#ifdef  __cplusplus
extern "C" {
#endif
void fexact(int * nrow, int * ncol, double * table, int * ldtabl,
	double * expect, double * percnt, double * emin, double * prt,
	double * pre, /* new in C : */ int * workspace);
#ifdef  __cplusplus
}
#endif
#endif
