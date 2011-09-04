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
   File:        amdyv_array.c
   Author:      Andrew W. Moore
   Created:     Tue Jun 15 12:31:18 EDT 2004
   Description: Extensions and advanced amdm stuff (no direct dym data access)

   Copyright 1996, Schenley Park Research

   This file contains advanced utility functions involving dyv_arrays, dyvs and
   ivecs. It never accesses the data structures directly, so if the
   underlying representation of dyms and dyvs changes these won't need to.

   The prototypes of these functions are declared in amdyv_array.h
*/

#include "amdyv_array.h"


/***** Now we'll play with dyv_arrays which are adjustable length
       arrays of dyvs
*****/

#define INITIAL_DYV_ARRAY_SIZE 10

dyv_array *mk_empty_dyv_array(void)
{
  dyv_array *da = AM_MALLOC(dyv_array);
  da -> size = 0;
  da -> array_size = INITIAL_DYV_ARRAY_SIZE;
  da -> array = AM_MALLOC_ARRAY(dyv_ptr,da->array_size);
  return(da);
}

void add_to_dyv_array(dyv_array *da, const dyv *dv)
/*
     Assume dyv_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of dv.
*/
{
  if ( da -> size == da -> array_size )
  {
    int new_size = 2 + 2 * da->array_size;
    dyv **new_array = AM_MALLOC_ARRAY(dyv_ptr,new_size);
    int i;
    for ( i = 0 ; i < da -> array_size ; i++ )
      new_array[i] = da->array[i];
    AM_FREE_ARRAY(da->array,dyv_ptr,da->array_size);
    da -> array = new_array;
    da -> array_size = new_size;
  }
  da->array[da->size] = (dv==NULL) ? NULL : mk_copy_dyv(dv);
  da->size += 1;
}

dyv_array *mk_const_dyv_array(dyv *base_vec, int size){
  int i;
  dyv_array *result = mk_empty_dyv_array();

  for (i=0; i<size; i++) {
    add_to_dyv_array(result, base_vec);
  }
  return result;
}

dyv_array *mk_rectangular_dyv_array( int numdyvs, int dyvlen)
{
  int i;
  dyv_array *dva;
  dyv*       temp;

  dva  = mk_dyv_array( numdyvs);
  temp = mk_dyv( dyvlen);

  for (i=0; i<numdyvs; ++i) dyv_array_set( dva, i, temp);

  free_dyv(temp);

  return dva;
}

dyv_array *mk_zero_dyv_array(int size){
  dyv *temp_dyv = mk_dyv(0);
  dyv_array *result = mk_const_dyv_array(temp_dyv, size);
  free_dyv(temp_dyv);
  return result;
}

int dyv_array_size(const dyv_array *da)
{
  return(da->size);
}

void dyv_array_set(dyv_array *dva, int idx, const dyv *dv)
{
  if ((idx < 0) || (dva == NULL) || (idx >= dva->size))
        my_error("dyv_array_set: called with incompatible arguments");
  if (dva->array[idx] != NULL)
        free_dyv(dva->array[idx]);
  dva->array[idx] = (dv == NULL) ? NULL : mk_copy_dyv(dv);
}

/* Use this function at your peril! */
void dyv_array_set_no_copy( dyv_array *dva, int idx, dyv *dv)
{
  if ((idx < 0) || (dva == NULL) || (idx >= dva->size)) {
    my_errorf( "dyv_array_set_no_copy: index %d is out of bounds [0,%d]\n",
	       idx, dva->size);
  }
  if (dva->array[idx] != NULL) free_dyv(dva->array[idx]);
  dva->array[idx] = dv;
  return;
}

void fprintf_dyv_array(FILE *s, const char *m1, const dyv_array *da,
                       const char *m2)
{
  if ( dyv_array_size(da) == 0 )
    fprintf(s,"%s = <dyv_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%s[%2d]",m1,i);
      fprintf_dyv(s,buff,dyv_array_ref(da,i),m2);
    }
  }
}

void free_dyv_array(dyv_array *da)
{
  int i;
  for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    if ( da->array[i] != NULL )
      free_dyv(da->array[i]);
  AM_FREE_ARRAY(da->array,dyv_ptr,da->array_size);
  AM_FREE(da,dyv_array);
}

dyv_array *mk_copy_dyv_array(const dyv_array *da)
{
  dyv_array *new_ar = mk_empty_dyv_array();
  int i;

  for ( i = 0 ; i < dyv_array_size(da) ; i++ )
    add_to_dyv_array(new_ar,dyv_array_ref(da,i));

  return(new_ar);
}

void dyv_array_remove(dyv_array *dva,int idx)
{
  int i;
  dyv *dv = dyv_array_ref(dva,idx);
  if ( dv != NULL ) free_dyv(dv);
  for ( i = idx ; i < dva->size-1 ; i++ )
    dva->array[i] = dva->array[i+1];
  dva->array[dva->size-1] = NULL;
  dva->size -= 1;
}

dyv_array *mk_array_of_zero_length_dyvs(int size)
{
  dyv_array *dva = mk_empty_dyv_array();
  dyv *dv = mk_dyv(0);
  int i;

  for (i = 0; i < size; i++)
        add_to_dyv_array(dva, dv);
  free_dyv(dv);
  return(dva);
}

dyv_array *mk_array_of_null_dyvs(int size)
{
  dyv_array *dva = mk_empty_dyv_array();
  int i;

  for (i = 0; i < size; i++)
    add_to_dyv_array(dva,NULL);

  return(dva);
}

dyv_array *mk_dyv_array( int size)
{
  int i;
  dyv_array *da;
 
  da = AM_MALLOC( dyv_array);
  da->size = size;
  da->array_size = size;
  da->array = AM_MALLOC_ARRAY( dyv_ptr, size);

  for (i=0; i<size; ++i) da->array[i] = NULL;

  return(da);
}

dyv_array *mk_dyv_array_subset( dyv_array *da, ivec *indices)
{
  int dasize, i, idx;
  dyv *dv;
  dyv_array *subda;

  dasize = ivec_size( indices);
  subda = mk_dyv_array( dasize);
  for (i=0; i<dasize; ++i) {
    idx = ivec_ref( indices, i);
    dv = dyv_array_ref( da, idx); /* Get from da[index]. */
    dyv_array_set( subda, i, dv);   /* Store in subda[i].  */
  }

  return subda;
}
