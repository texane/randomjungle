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
   File:        genarray.c
   Author:      Andrew W. Moore
   Created:     Sat Jun 26 03:22:44 EDT 1999
   Description: Generic Array of Objects

   Copyright 1999, Schenley Park Research
*/

#include "genarray.h"
#include "amma.h"

/***** Now we'll play with generic_arrays which are adjustable length
       arrays of generics
*****/

#define INITIAL_GENERIC_ARRAY_SIZE 10

generic_array *mk_empty_generic_array(
  void (*free_data)(void *data),
  void * (*mk_copy_data)(void *data),
  void (*fprintf_data)(FILE *s,char *m1,void *data,char *m2) /* May be NULL */
)
{
  generic_array *ga = AM_MALLOC(generic_array);
  ga -> size = 0;
  ga -> array_size = INITIAL_GENERIC_ARRAY_SIZE;
  ga -> array = AM_MALLOC_ARRAY(void *,ga->array_size);

  ga -> free_data = free_data;
  ga -> mk_copy_data = mk_copy_data;
  ga -> fprintf_data = fprintf_data;
  
  return(ga);
}

/*
   Assume generic_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is the same pointer to the same thing. After
   calling this, thing should never be directly referenced or 
   directly freed again.
*/
void add_pointer_to_generic_array(generic_array *ga,void *thing)
{
  if ( ga -> size == ga -> array_size )
  {
    int new_size = 2 + 2 * ga->array_size;
    void **new_array = AM_MALLOC_ARRAY(void *,new_size);
    int i;
    for ( i = 0 ; i < ga -> array_size ; i++ )
      new_array[i] = ga->array[i];
    AM_FREE_ARRAY(ga->array,void *,ga->array_size);
    ga -> array = new_array;
    ga -> array_size = new_size;
  }
  ga->array[ga->size] = thing;
  ga->size += 1;
}

/*
   Assume generic_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of thing.
*/
void add_to_generic_array(generic_array *ga,void *thing)
{
  add_pointer_to_generic_array(ga,
			       (thing==NULL) ? NULL : ga->mk_copy_data(thing));
}

int generic_array_size(generic_array *ga)
{
  return(ga->size);
}

void *generic_array_ref(generic_array *ga,int idx)
/*
   Returns a pointer (not a copy) to the index'th element stored in
   the generic_array. Error if index < 0 or index >= size
*/
{
  void *result;
  if ( idx < 0 || idx >= generic_array_size(ga) )
  {
    result = NULL;
    my_error("generic_array_ref");
  }
  else
    result = ga->array[idx];
  return(result);
}

void generic_array_set(generic_array *ga,int idx,void *element)
{
  if ( idx < 0 || idx >= generic_array_size(ga) )
    my_error("generic_array_set");
  else 
  {
    if ( ga->array[idx] != NULL ) ga->free_data(ga->array[idx]);
    ga->array[idx] = (element == NULL) ? NULL : ga->mk_copy_data(element);
  }
}
  
void fprintf_generic_array(FILE *s,char *m1,generic_array *ga,char *m2)
{
  if ( generic_array_size(ga) == 0 )
    fprintf(s,"%s = <generic_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for ( i = 0 ; i < generic_array_size(ga) ; i++ )
    {
      char buff[5000];
      void *element = generic_array_ref(ga,i);
      sprintf(buff,"%s[%2d]",m1,i);
      if ( element == NULL )
	fprintf(s,"%s = NULL%s",buff,m2);
      else
	ga->fprintf_data(s,buff,element,m2);
    }
  }
}

void free_generic_array(generic_array *ga)
{
  int i;
  for ( i = 0 ; i < generic_array_size(ga) ; i++ )
    if ( ga->array[i] != NULL ) ga->free_data(ga->array[i]);
  AM_FREE_ARRAY(ga->array,void *,ga->array_size);
  AM_FREE(ga,generic_array);
}

void generic_array_free_without_freeing_contents(generic_array *ga)
{
  AM_FREE_ARRAY(ga->array,void *,ga->array_size);
  AM_FREE(ga,generic_array);
}

generic_array *mk_copy_generic_array(generic_array *ga)
{
  generic_array *new_ar = 
    mk_empty_generic_array(ga->free_data,ga->mk_copy_data,ga->fprintf_data);

  int i;

  for ( i = 0 ; i < generic_array_size(ga) ; i++ )
    add_to_generic_array(new_ar,generic_array_ref(ga,i));

  return(new_ar);
}

/* Increases the length of ga by 1, puts a COPY of element in as the
   i'th element, and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one higher index.
*/
void insert_in_generic_array(generic_array *ga,int i,void *element)
{
  void *new_ptr;
  int j;
  
  if ( i < 0 || i > ga->size ) my_error("insert_in_generic_array: bad index");

  add_to_generic_array(ga,element);
  new_ptr = ga->array[ga->size-1];
  for ( j = ga->size-2 ; j >= i ; j-- )
    ga->array[j+1] = ga->array[j];
  ga->array[i] = new_ptr;
}

/* Decreases the length of ga by 1, removes the i'th element
   and any elements to the right of the i'th element (i.e.
   any elements with index >= i) are moved to having one lower index.
*/
void generic_array_remove(generic_array *ga,int i)
{
  int j;
  if ( i < 0 || i >= ga->size ) my_error("generic_array_remove: bad index");
  
  if ( ga->array[i] != NULL ) ga->free_data(ga->array[i]);

  for ( j = i ; j < ga->size-1 ; j++ )
    ga->array[j] = ga->array[j+1];
  ga->array[ga->size-1] = NULL;

  ga -> size -= 1;
}

/* The entry contents of ga[i] are AMFREED. The pointer to sp
   goes in its place with NO COPYING. After this call sp should never again
   be referenced directly or freed directly. */
void generic_array_set_pointer(generic_array *ga,int i,void *item)
{
  if ( ga->array[i] != NULL ) ga->free_data(ga->array[i]);
  ga->array[i] = item;
}

/* The pointer to the contents of ga[i] are overwritten with NULL.
   The previous contents are NOT freed, so someone else must have taken
   care of (or be about to take care of) that */
void generic_array_null_without_freeing(generic_array *genarr,int i)
{
  genarr->array[i] = NULL;
}

/* Only use this free function if
     (a) All items in genarr were added with set_pointer or
         add_pointer_to

         AND

     (b) Someone else is looking after the pointers to the contents
         and will eventually free them
*/
void free_generic_array_not_contents(generic_array *genarr)
{
  int i;
  for ( i = 0 ; i < generic_array_size(genarr) ; i++ )
    generic_array_null_without_freeing(genarr,i);
  free_generic_array(genarr);
}
