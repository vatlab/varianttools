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
** The half() SQL function returns half of its input value.
*/
static void halfFunc(
  sqlite3_context *context,
  int argc,
  sqlite3_value **argv
){
  sqlite3_result_double(context, 0.5*sqlite3_value_double(argv[0]));
}

/* SQLite invokes this routine once when it loads the extension.
** Create new functions, collating sequences, and virtual table
** modules here.  This is usually the only exported symbol in
** the shared library.
*/
int sqlite3_extension_init(
  sqlite3 *db,
  char **pzErrMsg,
  const sqlite3_api_routines *pApi
)
{
  SQLITE_EXTENSION_INIT2(pApi)
  sqlite3_create_function(db, "half", 1, SQLITE_ANY, 0, halfFunc, 0, 0);
  return 0;
}
