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
   File:        precs.c
   Author:      Andrew W. Moore
   Created:     Nov 2000
   Description: An array of precs

   Copyright 2000, Andrew W. Moore
*/

#include "am_string.h"
#include "precs.h"

void void_free_prec(void *data)
{
  free_prec((prec *)data);
}

void *void_mk_copy_prec(void *data)
{
  return (void *) mk_copy_prec((prec *)data);
}

void void_fprintf_prec(FILE *s,char *m1,void *data,char *m2)
{
  my_error("Don't call me");
}

precs *mk_empty_precs(void)
{
  precs *array = AM_MALLOC(precs);
  array -> genarr = mk_empty_generic_array(void_free_prec,
					   void_mk_copy_prec,
					   void_fprintf_prec);
  return array;
}

void add_to_precs(precs *array,prec *element)
{
  add_to_generic_array(array->genarr,(void *)element);
}

void add_pointer_to_precs(precs *array,prec *element)
{
  add_pointer_to_generic_array(array->genarr,(void *)element);
}

int precs_size(precs *array)
{
  return(generic_array_size(array->genarr));
}

prec *precs_ref(precs *array,int index)
{
  return (prec *) generic_array_ref(array->genarr,index);
}

void fprintf_precs(FILE *s,char *m1,precs *array,char *m2)
{
  if ( precs_size(array) == 0 )
    fprintf(s,"%s = <empty precs>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < precs_size(array) ; i++ )
    {
      prec *q = precs_ref(array,i);
      char *buff = mk_printf("%s[%d]",m1,i);
      fprintf_prec(s,buff,q,m2);
      free_string(buff);
    }
  }
}

void free_precs(precs *array)
{
  free_generic_array(array->genarr);
  AM_FREE(array,precs);
}

precs *mk_copy_precs(precs *array)
{
  precs *new_precs = AM_MALLOC(precs);
  new_precs -> genarr = mk_copy_generic_array(array->genarr);
  return new_precs;
}

/* Increases the length of qs by 1, puts a COPY of element in as the
   i'th element, and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one higher index.
*/
void insert_in_precs(precs *qs,int i,prec *rl)
{
  insert_in_generic_array(qs->genarr,i,(void *)rl);
}

/* Decreases the length of qs by 1, removes the i'th element
   and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one lower index.
*/
void precs_remove(precs *qs,int i)
{
  generic_array_remove(qs->genarr,i);
}

/* The entry contents of qs[i] are AMFREED. The pointer to q
   goes in its place with NO COPYING. After this call q should never again
   be referenced directly or freed directly. */
void precs_set_pointer(precs *qs,int i,prec *q)
{
  generic_array_set_pointer(qs->genarr,i,(void *)q);
}

