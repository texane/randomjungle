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
   File:        precs.h
   Author:      Andrew W. Moore
   Created:     Nov 2000
   Description: Header for An array of precs

   Copyright 2000, Andrew W. Moore
*/


#ifndef precs_H
#define precs_H

#include "genarray.h"
#include "prec.h"

typedef struct precs
{
  generic_array *genarr;
} precs;

precs *mk_empty_precs(void);
void add_to_precs(precs *array,prec *element);

/* Does NOT copy in the element, merely its pointer, so after calling, the
   user must NOT access or directly free element ever again */
void add_pointer_to_precs(precs *array,prec *element);

int precs_size(precs *array);
prec *precs_ref(precs *array,int index);
void fprintf_precs(FILE *s,char *m1,precs *array,char *m2);
void pprecs(precs *array);
void free_precs(precs *array);
precs *mk_copy_precs(precs *array);

/* Increases the length of sps by 1, puts a COPY of element in as the
   i'th element, and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one higher index.
*/
void insert_in_precs(precs *sps,int i,prec *rl);

/* Decreases the length of sps by 1, removes the i'th element
   and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one lower index.
*/
void precs_remove(precs *sps,int i);

/* The entry contents of sps[i] are AMFREED. The pointer to sp
   goes in its place with NO COPYING. After this call sp should never again
   be referenced directly or freed directly. */
void precs_set_pointer(precs *sps,int i,prec *sp);


#endif /* #ifndef precs_H */
