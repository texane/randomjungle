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
   File:        lrutils.h
   Author:      Paul Komarek
   Created:     Wed May 18 14:46:13 EDT 2005
   Description: Utilities

   Copyright 2005, The Auton Lab, CMU
*/



#ifndef LRUTILS
#define LRUTILS



#include <stdio.h>

#include "file.h"
#include "amdyv.h"


/****************************************************************************/
/* ivec utilities                                                           */
/****************************************************************************/

int count_in_ivec( ivec *iv, int val);



/****************************************************************************/
/* dyv utilities                                                            */
/****************************************************************************/

void dyv_write( PFILE *f, const dyv *dv);
void dyv_write_fname( char *filename, const dyv *dv);
dyv *mk_dyv_read( PFILE *f);
dyv *mk_dyv_read_fname( char *filename);



#endif
