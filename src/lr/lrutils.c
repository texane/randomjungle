/*
  Logistic Regression using Truncated Iteratively Re-weighted Least Squares
  (includes several programs)
  Copyright (C) 2005  Paul Komarek

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

  Author: Paul Komarek, komarek@cmu.edu
  Alternate contact: Andrew Moore, awm@cs.cmu.edu

*/


/*
   File:        lrutils.c
   Author:      Paul Komarek
   Created:     Wed May 18 14:44:23 EDT 2005
   Description: Utilities

   Copyright 2005, The Auton Lab, CMU
*/


#include <stdio.h>
#include <stdlib.h>

#include "am_string.h"
#include "amdyv.h"

#include "lrutils.h"




/****************************************************************************/
/* ivec utilities                                                           */
/****************************************************************************/

/* Returns the number of times val occurs in iv. */
int count_in_ivec( ivec *iv, int val)
{
  int i, count;
  count = 0;
  for (i=0; i<ivec_size( iv); ++i) if (ivec_ref( iv, i) == val) count += 1;
  return count;
}



/****************************************************************************/
/* dyv utilities                                                            */
/****************************************************************************/

void dyv_write( PFILE *f, const dyv *dv)
{
  int size, i;
  size = dyv_size( dv);
  pfprintf( f, "%d\n", size);
  for (i=0; i<size; ++i) pfprintf( f, "%.16f\n", dyv_ref(dv,i));
  return;
}

void dyv_write_fname( char *filename, const dyv *dv)
{
  PFILE *f;
  f = sure_pfopen( filename, "wb");
  dyv_write( f, dv);
  pfclose( f);
  return;
}

dyv *mk_dyv_read( PFILE *f)
{
  int i, size, lineno;
  double val;
  char line[101];
  dyv *dv;

  lineno = 1;
  line[100] = '\0';

  /* Read size and make dyv. */
    if (pfeof(f)) {
      my_errorf( "mk_dyv_read: unexpected end-of-file while reading size,\n"
                 "after line %d of file", lineno);
    }
  if (pfgets( line, 100, f) == NULL) {
    my_errorf( "mk_dyv_read: failed to read line %d from the passed stream.",
               lineno);
  }
  else lineno++;
  size = atoi( line);
  dv = mk_dyv( size);


  /* Read values. */
  for (i=0; i<size; ++i) {
    if (pfeof(f)) {
      my_errorf( "mk_dyv_read: unexpected end-of-file while reading %d vals,\n"
                 "after line %d of file (after the %dth value)",
                 size, lineno, lineno-1);
    }
    if (pfgets( line, 100, f) == NULL) {
      my_errorf( "mk_dyv_read: failed to read line %d from the passed stream.",
                 lineno);
    }
    else lineno++;

    val = atof( line);
    dyv_set( dv, i, val);
  }

  return dv;
}

dyv *mk_dyv_read_fname( char *filename)
{
  PFILE *f;
  dyv *dv;
  f = sure_pfopen( filename, "rb");
  dv = mk_dyv_read( f);
  pfclose( f);
  return dv;
}
