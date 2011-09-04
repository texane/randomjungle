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


/* *
   File:        amiv.c
   Author:      Andrew W. Moore
   Created:     Sat Apr  8 18:48:26 EDT 1995
   Updated:     29 March 97
   Description: Integer and Boolean Dynamic vectors

   Copyright 1997, Schenley Park Research
*/

#include "amiv.h"
#include "am_string.h"
#include "am_string_array.h"
#include "amma.h"

#define IVEC_CODE 20541

int Ivecs_mallocked = 0;
int Ivecs_freed = 0;

/* Headers for private functions. Are these all meant to be private? */
int ivec_arg_extreme(const ivec *iv,bool lowest);
int ivec_extreme(const ivec *iv,bool lowest);

ivec *mk_ivec(int size)
{
  ivec *result = AM_MALLOC(ivec);
  if ( size < 0 ) my_error("mk_ivec : size < 0 illegal");
  result -> ivec_code = IVEC_CODE;
  result -> array_size = size;
  result -> size = size;
  result -> iarr = am_malloc_ints(size);
  Ivecs_mallocked += 1;
  return(result);
}

/* Acts as though you created a ivec of size 0,        */
/* but actually allocates an array of size "capacity". */
/* This is very useful if you want to use add_to_ivec  */
/* and happen to have a reasonable upper bound to the  */
/* number of elements you want to add.                 */
ivec *mk_empty_ivec(int capacity) {
  ivec *result = AM_MALLOC(ivec);
  if ( capacity < 0 ) my_error("mk_empty_ivec : capacity < 0 illegal");
  result -> ivec_code = IVEC_CODE;
  result -> array_size = capacity;
  result -> size = 0;
  result -> iarr = am_malloc_ints(capacity);
  Ivecs_mallocked += 1;
  return(result);
}


ivec *mk_ivec_x( int size, ...)
{
  /* Warning: no type checking can be done by the compiler.  You *must*
     send the values as integers for this to work correctly. */
  int i, ival;
  va_list argptr;
  ivec *iv;
  
  iv = mk_ivec(size);

  va_start( argptr, size);
  for (i=0; i<size; ++i) {
    ival = va_arg( argptr, int);
    ivec_set( iv, i, ival);
  }
  va_end(argptr);

  return iv;
}

void free_ivec(ivec *iv)
{
  iv -> ivec_code = 7777;

  am_free_ints(iv->iarr,iv->array_size);
  AM_FREE(iv,ivec);

  Ivecs_freed += 1;
}

void fprintf_ivec(FILE *s,char *m1, const ivec *iv,char *m2)
{
  if ( iv == NULL )
    fprintf(s,"%s = (ivec *)NULL%s",m1,m2);
  else
  {
    char buff[256];

    /* Modified so that we don't cause a seg fault when m1 is bigger than
       the buffer - WKW*/
    if( strlen(m1) <= 255 )
      {
	strncpy(buff,m1,strlen(m1)+1);
      }
    else
      {
	strncpy(buff,m1,sizeof(buff));
	buff[255] = '\0';
      }
    fprintf_ints(s,buff,iv->iarr,iv->size,m2);
  }
}

void copy_iarr_to_ivec(int *iarr,int size,ivec *r_iv)
{
  copy_ints(iarr,r_iv->iarr,size);
}

void copy_ivec_to_iarr(ivec *iv, int *iarr)
{
  copy_ints(iv->iarr,iarr,iv->size);
}

int *mk_iarr_from_ivec(ivec *iv)
{
  int *result;
  result = am_malloc_ints(iv->size);
  copy_ivec_to_iarr(iv,result);
  return(result);
}

/* Makes an ivec of size end - start:
   { start , start+1 , .... end-2 , end-1 } */
ivec *mk_sequence_ivec(int start_value,int end_value)
{
  int size = end_value - start_value;
  ivec *iv = mk_ivec(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    ivec_set(iv,i,start_value+i);
  return iv;
}

/*
   Allocates and returns an ivec of size size in which ivec[i] = i
   ivec[0] = 0
   ivec[1] = 1
    .
    .
   ivec[size-1] = size-1
*/
ivec *mk_identity_ivec(int size)
{
  return mk_sequence_ivec(0,size);
}

void shuffle_ivec(ivec *iv)
/*
   A random permutation of iv is returned
*/
{
  int size = ivec_size(iv);
  int i;
  for ( i = 0 ; i < size-1 ; i++ )
  {
    int j = int_random(size - i);
    if ( j > 0 )
    {
      int swap_me_1 = ivec_ref(iv,i);
      int swap_me_2 = ivec_ref(iv,i+j);
      ivec_set(iv,i,swap_me_2);
      ivec_set(iv,i+j,swap_me_1);
    }
  }
}

void constant_ivec(ivec *iv,int value)
{
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    ivec_set(iv,i,value);
}

ivec *mk_constant_ivec(int size,int value)
{
  ivec *iv = mk_ivec(size);
  constant_ivec(iv,value);
  return(iv);
}

void zero_ivec(ivec *iv)
{
  constant_ivec(iv,0);
}

ivec *mk_zero_ivec(int size)
{
  return(mk_constant_ivec(size,0));
}

void ivec_scalar_plus(const ivec *iv, int addme, ivec *riv)
{
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    ivec_set(riv,i,ivec_ref(iv,i) + addme);
}

ivec *mk_ivec_scalar_plus(ivec *iv,int delta)
{
  ivec *result = mk_ivec(ivec_size(iv));
  ivec_scalar_plus(iv,delta,result);
  return(result);
}

void ivec_scalar_mult(const ivec *iv,int scale,ivec *riv)
{
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    ivec_set(riv,i,ivec_ref(iv,i) * scale);
}

ivec *mk_ivec_scalar_mult(ivec *iv,int scale)
{
  ivec *result = mk_ivec(ivec_size(iv));
  ivec_scalar_mult(iv,scale,result);
  return(result);
}

void copy_ivec(const ivec *iv,ivec *r_iv)
{
  int i, size;
  size = iv->size;
  for (i=0; i<size; ++i) ivec_set(r_iv, i, ivec_ref(iv, i));
  return;
  /* ivec_scalar_mult(iv,1,r_iv); */
}

ivec *mk_copy_ivec(const ivec *iv)
{
  ivec *result = mk_ivec(ivec_size(iv));
  copy_ivec(iv,result);
  return(result);
}

void copy_ivec_subset(ivec *iv, int start, int size, ivec *r_iv)
{
  int i;
  for (i = start; i < start + size; i++)
    ivec_set(r_iv, i - start, ivec_ref(iv, i));
}

ivec *mk_copy_ivec_subset(ivec *iv, int start, int size)
{
  ivec *result = mk_ivec(size);
  copy_ivec_subset(iv, start, size, result);
  return (result);
}

int num_of_given_value(ivec *iv,int value)
{
  int i;
  int result = 0;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    if ( ivec_ref(iv,i) == value ) result += 1;

  return(result);
}

int num_zero_entries(ivec *iv)
{
  return(num_of_given_value(iv,0));
}

int num_nonzero_entries(ivec *iv)
{
  return(ivec_size(iv) - num_zero_entries(iv));
}

int ivec_arg_extreme(const ivec *iv,bool lowest)
{
  int extreme_index;
  int extremest_value;
  int i;

  if ( ivec_size(iv) <= 0 )
    my_error("Can't find min or max of empty ivec");
  
  extreme_index = 0;
  extremest_value = ivec_ref(iv,0);

  for ( i = 1 ; i < ivec_size(iv) ; i++ )
  {
    int v = ivec_ref(iv,i);
    if ( (lowest && v < extremest_value) || (!lowest && v > extremest_value) ) 
    {
      extremest_value = v;
      extreme_index = i;
    }
  }

  return(extreme_index);
}

int ivec_extreme(const ivec *iv,bool lowest)
{
  return ivec_ref(iv,ivec_arg_extreme(iv,lowest));
}

int ivec_min(const ivec *iv)
{
  return ivec_ref(iv,ivec_arg_extreme(iv,TRUE));
}

int ivec_max(const ivec *iv)
{
  return ivec_ref(iv,ivec_arg_extreme(iv,FALSE));
}

int ivec_argmin(const ivec *iv)
{
  return ivec_arg_extreme(iv,TRUE);
}

int ivec_argmax(const ivec *iv)
{
  return ivec_arg_extreme(iv,FALSE);
}

/* Sensible if args are NULL. False if different size */
bool ivec_equal(const ivec *x1,const ivec *x2)
{
  bool result = TRUE;

  if ( EQ_PTR(x1,x2) )
    result = TRUE;
  else if ( x1 == NULL || x2 == NULL )
    result = FALSE;
  else if ( ivec_size(x1) != ivec_size(x2) )
    result = FALSE;
  else
  {
    int i;
    for ( i = 0 ; result && i < ivec_size(x1) ; i++ ) 
      result = result && ivec_ref(x1,i) == ivec_ref(x2,i);
  }
  return(result);
}

/**** Removal functions on ivecs ****/

/* Reduces the size of iv by one.
   ivec_ref(iv,index) disappears.
   Everything to the right of ivec_ref(iv,index) is copied one to the left.

   Formally: Let ivold be the ivec value beore calling this function
             Let ivnew be the ivec value after calling this function.

PRE: ivec_size(ivold) > 0
     0 <= index < ivec_size(ivold)

POST: ivec_size(ivnew) = ivec_size(ivold)-1
      for j = 0 , 1, 2 ... index-1  : 
         ivec_ref(ivnew,j) == ivec_ref(ivold,j)

      for j = index , index+1 , ... ivec_size(ivnew)-1:
         ivec_ref(ivnew,j) == ivec_ref(ivold,j+1)
*/
void ivec_remove(ivec *iv,int idx)
{
  int i;
  int isize = ivec_size(iv);

#ifndef AMFAST
  if ( isize <= 0 ) my_error("ivec_remove: empty ivec");
  if ( idx < 0 || idx >= isize )
    my_error("ivec_remove: bad index");
#endif /* #ifndef AMFAST */

  for ( i = idx ; i < isize - 1 ; i++ )
    ivec_set(iv,i,ivec_ref(iv,i+1));
  iv -> size -= 1;
}

int ivec_sum(const ivec *iv)
{
  int result = 0;
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    result += ivec_ref(iv,i);
  return(result);
}

bool equal_ivecs(const ivec *iv1, const ivec *iv2)
{
  bool result = TRUE;
  int i;

  if ( ivec_size(iv1) != ivec_size(iv2) )
    my_error("qual_ivecs wrong sizes");

  for ( i = 0 ; i < ivec_size(iv1) && result ; i++ )
    result = ivec_ref(iv1,i) == ivec_ref(iv2,i);
  return(result);
}

void add_to_ivec(ivec *iv,int new_val)
{
  int size;


  if ( iv->array_size < 0 )
    my_error("An ivec should never have a negative array_size if\n"
	     "it was constructed using legal API calls");

  if ( iv->size > iv->array_size )
    my_error("An ivec should never have size greater than array_size if\n"
    "it was constructed using legal API calls");

  size = iv->size;

  if ( size == iv->array_size )
  {
    int new_array_size = 2 * size + 10;
    int *iarr_new = AM_MALLOC_ARRAY(int,new_array_size);
    int i;

    for ( i = 0 ; i < size ; i++ )
      iarr_new[i] = iv->iarr[i];
    AM_FREE_ARRAY(iv->iarr,int,size);
    iv -> iarr = iarr_new;
    iv -> array_size = new_array_size;
  }
  iv->iarr[size] = new_val;
  iv -> size += 1;
}


void add_to_ivec_unique(ivec *iv,int val) {
  printf("AWM comments: You should probably be using sivecs\n");
  if (!is_in_ivec(iv, val)) add_to_ivec(iv, val);
}

/* Increases iv in length by 1 and shifts all elements
   with original index greater or equal to index one to the
   right and inserts val at index. */
void ivec_insert(ivec *iv,int idx,int val)
{
  int i;
  add_to_ivec(iv,0);
  for ( i = ivec_size(iv)-1 ; i > idx ; i-- )
    ivec_set(iv,i,ivec_ref(iv,i-1));
  ivec_set(iv,idx,val);
}

/* Creates a new ivec with the given element inserted.
   This is more efficient than copying and inserting. */
ivec *mk_ivec_insert( ivec *iv, int idx, int val)
{
  int i, newsize;
  ivec *result;

  newsize = ivec_size( iv) + 1;
  result = mk_ivec( newsize);

#ifndef AMFAST
  if (idx >= newsize || idx < 0) {
    my_errorf( "mk_ivec_insert: idx=%d is out-of-bounds [0,%d]",
               idx, newsize-1);
  }
#endif

  for (i=0; i<idx; ++i) ivec_set( result, i, ivec_ref( iv, i));
  ivec_set( result, idx, val);
  for (i=idx+1; i<newsize; ++i) ivec_set( result, i, ivec_ref( iv, i-1));
  
  return result;
}



/* Find least index i such that value = ivec_ref(iv,i).
  If not found, returns -1
*/
int find_index_in_ivec(const ivec *iv, int value)
{
  int result = -1;
  int i;

  for ( i = 0 ; i < ivec_size(iv) && result < 0 ; i++ )
    if (value == ivec_ref(iv,i)) 
      result = i;
  return(result);
}

/* make and return an ivec containing the indices of iv that have value */
ivec *mk_indices_in_ivec(const ivec *iv, int value)
{
  int i;
  ivec *res = mk_ivec(0);
  for (i=0;i<ivec_size(iv);i++) if (ivec_ref(iv,i)==value) add_to_ivec(res,i);
  return res;
}

int find_index_in_ivec_hint(const ivec *iv, int value, int hint)
{
  int i;
  int sign = -1;
  int sz = ivec_size(iv);

  for (i=1; i<=sz; i++) {
	int disp = i/2;
	int idx = (hint + sign*disp);
    /* make sure idx is in bounds. Adding sz will take care of negative idx.*/
    idx = (idx + sz) % sz;
    if (value == ivec_ref(iv, idx)) 
      return idx;
    sign *= -1;
  }
  return -1;
}

/* Finds leftmost index such that siv[index] > val. If none
   such, returns ivec_size(siv) */
int index_in_sorted_ivec(ivec *siv,int val)
{
  int result = -1;
  int i;
  for ( i = 0 ; i < ivec_size(siv) && result < 0 ; i++ )
    if ( ivec_ref(siv,i) > val )
      result = i;
  if ( result < 0 ) result = ivec_size(siv);
  return(result);
}

int find_in_sorted_ivec(ivec *siv, int val)
{
  int begin = 0;
  int end = ivec_size(siv);
  int cur = (begin+end)/2;
  while (end > begin)
  {
    if (ivec_ref(siv,cur) == val) return cur;
    else if (ivec_ref(siv,cur) < val) begin = cur+1;
    else                              end = cur;
    cur = (begin+end)/2;
  }
  return -1;
}

void add_to_sorted_ivec(ivec *siv, int val)
{
  int i;
  add_to_ivec(siv, 0);
  for (i = ivec_size(siv) - 2; i >= 0 && val < ivec_ref(siv, i); i--)
    ivec_set(siv, i + 1, ivec_ref(siv, i));
  ivec_set(siv, i + 1, val);
}

bool is_in_ivec( const ivec *iv,int value)
{
  return(find_index_in_ivec(iv,value) >= 0);
}

/* Returns list of indices at which value appears. */
ivec *mk_is_in_ivec( ivec *iv, int value)
{
  int size, num, i;
  ivec *where, *r_iv;

  /* Find locations of value. */
  size = ivec_size( iv);
  where = mk_ivec( size);
  num = 0;
  for (i=0; i<size; ++i) {
    if (value == ivec_ref( iv, i)) {
      ivec_set( where, num, i);
      num += 1;
    }
  }

  /* Copy useful part of where into r_iv. */
  r_iv = mk_ivec( num);
  for (i=0; i<num; ++i) ivec_set( r_iv, i, ivec_ref( where, i));
  free_ivec( where);

  return r_iv;
}

/* return the number of elements appearing in both ivecs */
int count_ivec_overlap(ivec *a, ivec *b)
{
  int i;
  int result = 0;

  for (i=0;i<ivec_size(a);i++) if (is_in_ivec(b,ivec_ref(a,i))) result++;

  return result;
}

/* It is fine for iv and r_iv to occupy the same memory
   (see amar.c:sort_ints). */
void ivec_sort(const ivec *iv,ivec *r_iv)
{
  copy_ivec(iv,r_iv);
  sort_ints(r_iv->iarr,ivec_size(r_iv),r_iv->iarr);
}

ivec *mk_ivec_sort(const ivec *iv)
{
  ivec *result = mk_ivec(ivec_size(iv));
  ivec_sort(iv,result);
  return result;
}

/*
 Creates a ivec of indices such that indices[i] is the origional
 location (in the unsorted iv) of the ith smallest value.
 Used when you want the location of the sorted values instead of 
 the sorted vector.
*/
ivec *mk_sorted_ivec_indices(ivec *iv)
{
  ivec* indices;
  int size = ivec_size(iv);
  int i;
  int* iarr_ind;
  int* iarr;

  iarr = mk_iarr_from_ivec(iv);
  iarr_ind = am_malloc_ints(size);

  indices_sort_integers(iarr,size,iarr_ind);
  indices = mk_ivec(size);

  for(i=0;i<size;i++)
  {
    ivec_set(indices,i,iarr_ind[i]);
  }

  am_free_ints(iarr_ind,size);
  am_free_ints(iarr,size);

  return indices;
}

string_array *mk_string_array_from_ivec(ivec *iv)
{
  string_array *sa;
  
  if ( ivec_size(iv) == 0 )
  {
    sa = mk_string_array(1);
    string_array_set(sa,0,"empty");
  }
  else
  {
    int i;
    sa = mk_string_array(ivec_size(iv));
    for ( i = 0 ; i < ivec_size(iv) ; i++ )
    {
      char buff[100];
      sprintf(buff,"%d",ivec_ref(iv,i));
      string_array_set(sa,i,buff);
    }
  }
  return(sa);
}


/***** ivec_array time! ****/

#define INITIAL_ivec_ARRAY_SIZE 10

ivec_array *mk_ivec_array(int size)
{
  int i;
  ivec_array *iva;
  iva = AM_MALLOC( ivec_array);
  iva->size = size;
  iva->array_size = int_max( size, INITIAL_ivec_ARRAY_SIZE);
  iva->array = AM_MALLOC_ARRAY( ivec_ptr, iva->array_size);
  for (i = 0; i < iva->array_size; i++) iva->array[i] = NULL;
  return iva;
}

ivec_array *mk_empty_ivec_array(void)
{
  ivec_array *ivecarr = AM_MALLOC(ivec_array);
  ivecarr -> size = 0;
  ivecarr -> array_size = INITIAL_ivec_ARRAY_SIZE;
  ivecarr -> array = AM_MALLOC_ARRAY(ivec_ptr,ivecarr->array_size);
  return(ivecarr);
}

ivec_array *mk_const_ivec_array(ivec *base_vec, int size){
  int i;
  ivec_array *result = mk_empty_ivec_array();
  for (i=0; i<size; i++) {
    add_to_ivec_array(result, base_vec);
  }
  return result;
}

ivec_array *mk_zero_ivec_array(int size){
  ivec *temp_vec = mk_ivec(0);
  ivec_array *result = mk_const_ivec_array(temp_vec, size);
  free_ivec(temp_vec);
  return result;
}

/* does the same thing as mk_zero_ivec_array */
ivec_array *mk_array_of_zero_length_ivecs(int size)
{
  ivec_array *iva = mk_empty_ivec_array();
  ivec *iv = mk_ivec(0);
  int i;

  for (i = 0; i < size; i++)
        add_to_ivec_array(iva, iv);
  free_ivec(iv);
  return(iva);
}

void add_to_ivec_array(ivec_array *ivecarr,const ivec *this_ivec)
/*
   Assume ivec_array is previously of size n. After this it is of size
   n+1, and the n+1'th element is a COPY of this_ivec.
*/
{
  if ( ivecarr -> size == ivecarr -> array_size )
  {
    int new_size = 2 + 2 * ivecarr->array_size;
    ivec **new_array = AM_MALLOC_ARRAY(ivec_ptr,new_size);
    int i;
    for ( i = 0 ; i < ivecarr -> array_size ; i++ )
      new_array[i] = ivecarr->array[i];
    AM_FREE_ARRAY(ivecarr->array,ivec_ptr,ivecarr->array_size);
    ivecarr -> array = new_array;
    ivecarr -> array_size = new_size;
  }
  ivecarr->array[ivecarr->size] = (this_ivec==NULL) ? NULL : mk_copy_ivec(this_ivec);
  ivecarr->size += 1;
}

int ivec_array_size(const ivec_array *ivecarr)
{
  return(ivecarr->size);
}

/* Returns the sum of all ivec_size(...) values of all ivecs
   in iva */
int sum_of_ivec_array_sizes(ivec_array *iva)
{
  int result = 0;
  int i;
  for ( i = 0 ; i < ivec_array_size(iva) ; i++ )
    result += ivec_size(ivec_array_ref(iva,i));
  return result;
}

void ivec_array_set(ivec_array *iva, int idx, const ivec *iv)
{
  if ((idx < 0) || (iva == NULL) || (idx >= iva->size))
        my_error("ivec_array_set: called with incompatible arguments");
  if (iva->array[idx] != NULL)
        free_ivec(iva->array[idx]);
  iva->array[idx] = (iv == NULL) ? NULL : mk_copy_ivec(iv);
}

void fprintf_ivec_array(FILE *s,char *m1,ivec_array *ivecarr,char *m2)
{
  if ( ivec_array_size(ivecarr) == 0 )
    fprintf(s,"%s = <ivec_array with zero entries>%s",m1,m2);
  else
  {
    int i;
    for (i=0; i<ivec_array_size(ivecarr); i++){
      char buff[100];
      ivec *iv = ivec_array_ref(ivecarr,i);
      sprintf(buff,"%s[%2d]",m1,i);
      if ( iv == NULL )
        fprintf(s,"%s = NULL%s",buff,m2);
      else
        fprintf_ivec(s,buff,iv,m2);
    }
  }
}

void free_ivec_array(ivec_array *ivecarr)
{
  int i;
  for ( i = 0 ; i < ivec_array_size(ivecarr) ; i++ )
    if ( ivecarr->array[i] != NULL )
      free_ivec(ivecarr->array[i]);
  AM_FREE_ARRAY(ivecarr->array,ivec_ptr,ivecarr->array_size);
  AM_FREE(ivecarr,ivec_array);
}

ivec_array *mk_copy_ivec_array(const ivec_array *ivecarr)
{
  ivec_array *new_ar = mk_empty_ivec_array();
  int i;

  for ( i = 0 ; i < ivec_array_size(ivecarr) ; i++ )
    add_to_ivec_array(new_ar,ivec_array_ref(ivecarr,i));

  return(new_ar);
}

bool ivec_array_equal(ivec_array *iva1,ivec_array *iva2)
{
  int size = ivec_array_size(iva1);
  bool result = size == ivec_array_size(iva2);
  int i;
  for ( i = 0 ; result && i < size ; i++ )
    result = ivec_equal(ivec_array_ref(iva1,i),ivec_array_ref(iva2,i));
  return result;
}

void ivec_array_remove(ivec_array *iva,int idx)
{
  int i;
  ivec *iv = ivec_array_ref(iva,idx);
  if ( iv != NULL ) free_ivec(iv);
  for ( i = idx ; i < iva->size-1 ; i++ )
    iva->array[i] = iva->array[i+1];
  iva->array[iva->size-1] = NULL;
  iva->size -= 1;
}

/* Returns an ivec_array of length equal to ivec_size(rows)
   in which result[i] = iva[rows[i]] */
ivec_array *mk_ivec_array_subset(ivec_array *iva,ivec *rows)
{
  int size, row, i;
  ivec *iv, *myrows;
  ivec_array *result;

  if ( rows == NULL) myrows = mk_identity_ivec( ivec_array_size( iva));
  else myrows = rows;
  size = ivec_size( myrows);
  result = mk_ivec_array( size);

  for (i=0; i < size; i++) {
    row = ivec_ref( myrows, i);
    iv = ivec_array_ref( iva, row);
    ivec_array_set( result, i, iv);
  }

  if (rows == NULL) free_ivec( myrows);
  return result;
}


/* Returns the max value in any of the ivecs in iva.
   PRE: Contains at least one non-zero-length ivec */
int ivec_array_max_value(ivec_array *iva)
{
  int result = -1;
  bool started = FALSE;
  int i;

  for ( i = 0 ; i < ivec_array_size(iva) ; i++ )
  {
    ivec *iv = ivec_array_ref(iva,i);
    if ( ivec_size(iv) > 0 )
    {  
      int m = ivec_max(iv);
      if ( started )
	result = int_max(result,m);
      else
      {
	result = m;
	started = TRUE;
      }
    }
  }

  if ( !started )
    my_error("ivec_array_max_value: no entries or all zero-length entries");

  return result;
}

/* Returns ivec of size ivec_size(rows) in which
    result[i] = x[rows[i]] */
ivec *mk_ivec_subset(ivec *x,ivec *rows)
{
  int size = ivec_size(rows);
  ivec *y = mk_ivec(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    ivec_set(y,i,ivec_ref(x,ivec_ref(rows,i)));
  return y;
}



/* Returns the minimum value in sivec. 
   Time cost: Constant */
int sivec_max(const ivec *siv)
{
  return ivec_ref(siv,ivec_size(siv)-1);
}



void sivres_update(sivres *sr,int value)
{
  if ( sr->intersection == NULL )
    sr -> intersection_size += 1;
  else
  {
    my_assert(ivec_size(sr->intersection) == 0 || 
	     ivec_ref(sr->intersection,ivec_size(sr->intersection)-1) < value);
    add_to_ivec(sr->intersection,value);
  }
}

void sivseg_cut(sivseg *old,sivseg *left,sivseg *right,
		int left_hi,int right_lo)
{
  left -> iv = old -> iv;
  left -> lo = sivseg_lo(old);
  left -> hi = left_hi;

  right -> iv = old -> iv;
  right -> lo = right_lo;
  right -> hi = sivseg_hi(old);
}

int sivseg_value_to_index(sivseg *ss,int value,bool *r_exists)
{
  int lo = sivseg_lo(ss);
  int hi = sivseg_hi(ss);
  int loval = sivseg_index_to_value(ss,lo);
  int hival = sivseg_index_to_value(ss,hi);
  int result = -7777;

  if ( value < loval )
  {
    result = lo-1;
    *r_exists = FALSE;
  }
  else if ( value > hival )
  {
    result = hi;
    *r_exists = FALSE;
  }
  else if ( value == loval )
  {
    result = lo;
    *r_exists = TRUE;
  }
  else if ( value == hival )
  {
    result = hi;
    *r_exists = TRUE;
  }
  else
  {
    bool found = FALSE;

    while ( hi > lo+1 && !found )
    {
      int mid = (lo + hi) / 2;
      int midval = sivseg_index_to_value(ss,mid);
   
      my_assert(mid > lo);
      my_assert(mid < hi);
      my_assert(loval < value);
      my_assert(hival > value);

      if ( midval == value )
      {
	found = TRUE;
	result = mid;
      }
      else if ( midval < value )
      {
	lo = mid;
	loval = midval;
      }
      else
      {
	hi = mid;
	hival = midval;
      }
    }

    if ( found )
      *r_exists = TRUE;
    else
    {
      *r_exists = FALSE;
      result = lo;
    }
  }

  return result;
}



void sivseg_search(sivseg *a,sivseg *b,sivres *sr)
{
  if ( sivseg_size(a) == 0 || sivseg_size(b) == 0 )
  {
    /* do nothing */
  }
  else if ( sivseg_lo_value(a) > sivseg_hi_value(b) )
  {
    /* do nothing */
  }
  else if ( sivseg_lo_value(b) > sivseg_hi_value(a) )
  {
    /* do nothing */
  }
  else
  {
    if ( sivseg_size(a) > sivseg_size(b) )
    {
      sivseg *c = a;
      a = b;
      b = c;
    }

    if ( sivseg_size(a) <= 2 )
    {
      int j;
      for ( j = sivseg_lo(a) ; j <= sivseg_hi(a) ; j++ )
      {
	bool exists;
	int value = sivseg_index_to_value(a,j);
	(void) sivseg_value_to_index(b,value,&exists);
	if ( exists )
	  sivres_update(sr,value);
      }
    }
    else
    {
      int mid = ( sivseg_lo(a) + sivseg_hi(a) ) / 2;
      int midval = sivseg_index_to_value(a,mid);
      bool exists;
      int bindex = sivseg_value_to_index(b,midval,&exists);
      sivseg aleft[1];
      sivseg aright[1];
      sivseg bleft[1];
      sivseg bright[1];

      sivseg_cut(a,aleft,aright,mid-1,mid+1);
      sivseg_cut(b,bleft,bright,(exists)?bindex-1:bindex,bindex+1);

      if ( sivseg_hi(bleft) >= sivseg_lo(bleft) ) {
	sivseg_search(aleft,bleft,sr);
      }

      if ( exists ) sivres_update(sr,midval);

      if ( sivseg_hi(bright) >= sivseg_lo(bright) )
	sivseg_search(aright,bright,sr);
    }
  }

  return;
}

/* define subset(siv,lo,hi) == { siv[lo] , siv[lo+1] , ... siv[hi] }
   thus size subset = hi - lo + 1.

   This function returns TRUE if and only if value is in subset(siv,lo,hi)
*/
bool is_in_sivec_between(const ivec *siv,int lo,int hi,int value)
{
  bool result;
  int vlo = ivec_ref(siv,lo);
  int vhi = ivec_ref(siv,hi);

  if ( vlo == value || vhi == value )
    result = TRUE;
  else if ( lo == hi )
    result = FALSE;
  else if ( value < vlo )
    result = FALSE;
  else if ( value > vhi )
    result = FALSE;
  else
  {
    result = FALSE;
    while ( !result && hi > lo+1 )
    {
      int mid = (lo + hi) / 2;
      int vmid = ivec_ref(siv,mid);
      if ( vmid == value )
	result = TRUE;
      else if ( vmid < value )
      {
	lo = mid;
	vlo = vmid;
      }
      else
      {
	hi = mid;
	vhi = vmid;
      }
    }
  }

  return result;
}
    
ivec *mk_new_sivec_intersection(const ivec *a,const ivec *b)
{
  sivres sr[1];
  sivseg ass[1];
  sivseg bss[1];

  sr -> intersection = mk_ivec(0);
  sr -> intersection_size = -77;

  ass -> iv = (ivec *) a;
  ass -> lo = 0;
  ass -> hi = ivec_size(a)-1;

  bss -> iv = (ivec *) b;
  bss -> lo = 0;
  bss -> hi = ivec_size(b)-1;

  sivseg_search(ass,bss,sr);

  my_assert(sr->intersection_size == -77);

  return sr -> intersection;
}

ivec *mk_sivec_intersection_basic(const ivec *siva,const ivec *sivb)
{
  ivec *bigresult, *result;
  int ai = 0;
  int bi = 0;
  int asize, bsize, nextidx;

  nextidx = 0;
  asize = ivec_size( siva);
  bsize = ivec_size( sivb);
  bigresult = mk_ivec( int_min( asize, bsize));

  while ( ai < ivec_size(siva) && bi < ivec_size(sivb) ) {
    while ( ai < ivec_size(siva) && ivec_ref(siva,ai) < ivec_ref(sivb,bi) ) {
      ai += 1;
    }
    if ( ai < ivec_size(siva) ) {
      while ( bi < ivec_size(sivb) && ivec_ref(sivb,bi) < ivec_ref(siva,ai) ) {
        bi += 1;
      }
      if ( bi < ivec_size(sivb) && ivec_ref(siva,ai) == ivec_ref(sivb,bi) ) {
	ivec_set( bigresult, nextidx, ivec_ref( siva, ai));
	nextidx += 1;
	ai += 1;
	bi += 1;
      }
    }
  }

  result = mk_copy_ivec_subset( bigresult, 0, nextidx);
  free_ivec( bigresult);

  return result;
}

ivec *mk_sivec_intersection(const ivec *siva,const ivec *sivb)
{
  int asize = ivec_size(siva);
  int bsize = ivec_size(sivb);
  int max_size = int_max(asize,bsize);
  int min_size = int_min(asize,bsize);
  ivec *result;

  if ( min_size * 100 < max_size ) {
    result = mk_new_sivec_intersection(siva,sivb);
  }
  else result = mk_sivec_intersection_basic(siva,sivb);
  return result;
}


/* define subset(siv,lo,hi) == { siv[lo] , siv[lo+1] , ... siv[hi] }
   thus size subset = hi - lo + 1.

   This function returns TRUE if and only if value is in subset(siv,lo,hi)
*/
int sio(const ivec *a,int la,int ha,const ivec *b,int lb,int hb)
{
  int result;

  if ( la == ha )
  {
    int va = ivec_ref(a,la);
    if ( lb == hb )
      result = ( va == ivec_ref(b,lb) ) ? 1 : 0;
    else
      result = ( is_in_sivec_between(b,lb,hb,va) ) ? 1 : 0;
  }
  else if ( lb == hb )
  {
    int vb = ivec_ref(b,lb);
    result = ( is_in_sivec_between(a,la,ha,vb) ) ? 1 : 0;
  }
  else
  {
    int vla = ivec_ref(a,la);
    int vlb = ivec_ref(b,lb);
    int vha = ivec_ref(a,ha);
    int vhb = ivec_ref(b,hb);
    
    if ( vha < vlb )
      result = 0;
    else if ( vhb < vla )
      result = 0;
    else
    {
      int ma = (la + ha) / 2;
      int mb = (lb + hb) / 2;
      result = sio(a,la,ma,b,lb,mb) + 
	       sio(a,la,ma,b,mb+1,hb) + 
               sio(a,ma+1,ha,b,lb,mb) + 
               sio(a,ma+1,ha,b,mb+1,hb);
    }
  }

  return result;
}

int size_of_sivec_intersection(const ivec *siva,const ivec *sivb)
{
  int result;

  if ( ivec_size(siva) == 0 ) result = 0;
  else if ( ivec_size(sivb) == 0 ) result = 0;
  else {
    result = sio(siva,0,ivec_size(siva)-1,sivb,0,ivec_size(sivb)-1);
  }

  return result;
}



/* Inverts ivec f, writing magic to values of inverse that don't occur in f.
   Passed ivec must be non-negative -- this is an unnecessary restriction,
   but 1) it seems like the usual case for calls to this function, and
   2) it avoid making the caller pass a "known good" magic value to use
   if the f is into and not onto [0,max(f)] */
ivec *mk_invert_nonneg_ivec( ivec *f)
{
  int size, max, i, val;
  ivec *finv;
  size  = ivec_size( f);
  max   = ivec_max( f);
  finv  = mk_constant_ivec( max+1, -1);
  /* Invert f. */
  for (i=0; i<size; ++i) {
    val = ivec_ref( f, i);
    ivec_set( finv, val, i);
  }
  return finv;
}

/* Reproducibly partition dataset rows for cross-validation. */
void make_kfold_rows( ivec *train_and_test_rows, int num_rows, int num_folds,
                      int fold_num, ivec **r_train_rows, ivec **r_test_rows)
{
  int save_seed = int_random(300000);
  ivec *srows = (train_and_test_rows==NULL) ? mk_identity_ivec(num_rows) :
                                              mk_copy_ivec(train_and_test_rows);
  int srows_size = ivec_size(srows);
  int start_i = (int)floor(fold_num * num_rows / (double) num_folds);
  int end_i = (int)floor((fold_num+1) * num_rows / (double) num_folds);
  int i;

  *r_train_rows = mk_ivec(srows_size - (end_i - start_i));
  *r_test_rows = mk_ivec(end_i - start_i);
  am_srand(12345);
  shuffle_ivec(srows);

  for ( i = 0 ; i < srows_size ; i++ )
  {
    ivec *update_me = (i >= start_i && i < end_i) ? *r_test_rows : *r_train_rows;
    int update_index = (i < start_i) ? i :
                       (i < end_i) ? i - start_i : i - (end_i - start_i);
    ivec_set(update_me,update_index,ivec_ref(srows,i));
  }
  free_ivec(srows);
  am_srand(save_seed);
  ivec_sort(*r_train_rows,*r_train_rows);
  ivec_sort(*r_test_rows,*r_test_rows);

  return;
}



ivec *mk_find_all_in_ivec( const ivec *iv, int val)
{
  int size, i, found;
  ivec *tmpiv, *result;

  size = ivec_size( iv);
  tmpiv = mk_ivec( size);

  found = 0;
  for (i=0; i<size; ++i) if (ivec_ref(iv,i)==val) ivec_set( tmpiv, found++, i);

  /* Copy tmpiv into properly-sized ivec. */
  result = mk_ivec( found);
  for (i=0; i<found; ++i) ivec_set( result, i, ivec_ref(tmpiv,i));
  free_ivec( tmpiv);

  return result;
}
