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



/* !!!!
   File:        amar.c
   Author:      Andrew W. Moore
   Created:     Sat Sep 12 14:53:13 EDT 1992
   Description: Obvious operations on 1-d arrays

   Copyright (C) 1992, Andrew W. Moore
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amar.h"      /* Obvious operations on 1-d arrays */

/* Private prototypes */
void os_sort_ints(int *farr, int size, int *r_farr);

#ifdef WIN32
int __cdecl double_comp(const void *r1_ptr, const void *r2_ptr);
int __cdecl compare_isf_elt(const void *i_ptr, const void *j_ptr);
#else
int double_comp(const void *r1_ptr, const void *r2_ptr);
int compare_isf_elt(const void *i_ptr, const void *j_ptr);
#endif

void realnums_permute(double *farr,int *indices,double *r_farr,int size);
void ints_permute(int *iarr,int *indices,int *r_iarr,int size);
void my_qs_2(double *farr, int start, int *indices);
void my_qs(double *farr, int size, int start, int end, int *indices);
bool indices_are_sorted(double *farr,int *indices,int size); /* This maybe should be public? */
void spr_indices_sort_realnums(double *farr, int size, int *indices);
void os_sort_realnums(double *farr, int size, double *r_farr);
int ints_argextreme(int *ints, int size, bool lowest);
int doubles_argextreme(double *doubles, int size, bool lowest);
void spr_sort_realnums(double *farr, int size, double *r_farr);
void os_indices_sort_realnums(double *farr, int size, int *indices);
int num_different_values_in_sorted_farr(double *farr,int size); /* Useful. Should be public. */
int frequency_of_most_common_value_in_sorted_farr(double *farr,int size);
bool avoid_os_sort_given_subset(double *subset,int subsize);
bool avoid_os_sort(double *farr,int size);
bool avoid_os_int_sort(int *iarr,int size);
void os_indices_sort_integers(int *iarr, int size, int *indices);
double *mk_farr_from_iarr(int *iarr,int size); /* Should be public */
void spr_indices_sort_integers(int *iarr,int size,int *indices);
void spr_sort_ints(int *iarr, int size, int *r_iarr);



void copy_ints(int *src, int *dest, int size)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    dest[i] = src[i];
}

void fprintf_ints(FILE *s, char *message1, int *ints, int size, char *message2)
{
  int i;
  fprintf(s,"%s = ",message1);

  if ( ints == NULL && size > 0 )
    fprintf(s,"(int *) NULL");
  else
  {
    fprintf(s,"{ ");
    for ( i = 0 ; i < size ; i++ )
      fprintf(s,"%d %s",ints[i],(i==size-1) ? "" : ", ");
    fprintf(s,"}");
  }

  fprintf(s," %s",message2);
}

void printf_ints(char *message1, int *ints, int size, char *message2)
{
  fprintf_ints(stdout,message1,ints,size,message2);
}

int ints_extreme(int *ints, int size, bool lowest)
{
  int result;
  if ( size <= 0 )
  {
    result = 0;
    my_error("ints_extreme: zero (or -ve) size");
  }
  else
  {
    int i;
    result = ints[0];
    for ( i = 1 ; i < size ; i++ )
      result = (lowest) ? int_min(result,ints[i]) : int_max(result,ints[i]);
  }
  return(result);
}

int ints_min(int *ints, int size)
{
  return(ints_extreme(ints,size,TRUE));
}

int ints_max(int *ints, int size)
{
  return(ints_extreme(ints,size,FALSE));
}

/*
* Frank Dellaert note (6/23/97)
* qsort expects cdecl, but the mblc DLL is compiled using
* __stdcall calling conventions
*/
#ifdef WIN32
int __cdecl int_comp(const void *r1_ptr, const void *r2_ptr)
#else
int int_comp(const void *r1_ptr, const void *r2_ptr)
#endif
{
  int *r1 = (int *) r1_ptr;
  int *r2 = (int *) r2_ptr;
  if ( *r1 < *r2 )
    return(-1);
  else if ( *r1 > *r2 )
    return(1);
  else
    return(0);
}

void os_sort_ints(int *farr, int size, int *r_farr)
/*
   It is fine for farr to be the same memory as r_farr
*/
{
  copy_ints(farr,r_farr,size);
  qsort((char *)r_farr,
        size,
        sizeof(int),
	int_comp
       );
}

void copy_realnums(double *src, double *dest, int size)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    dest[i] = src[i];
}

void fprintf_realnums(FILE *s, char *message1, double *doubles, int size, char *message2)
{
  int i;
  fprintf(s,"%s ",message1);

  if ( doubles == NULL && size > 0 )
    fprintf(s,"(double *) NULL");
  else
  {
    fprintf(s,"{ ");
    for ( i = 0 ; i < size ; i++ )
      fprintf(s,"%g %s",doubles[i],(i==size-1) ? "" : ", ");
    fprintf(s,"}");
  }

  fprintf(s," %s",message2);
}

void printf_realnums(char *message1, double *doubles, int size, char *message2)
{
  fprintf_realnums(stdout,message1,doubles,size,message2);
}

void set_realnums_constant(double *doubles, int size, double val)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    doubles[i] = val;
}

double doubles_extreme(double *doubles, int size, bool lowest)
{
  double result;
  if ( size <= 0 )
  {
    result = 0.0;
    my_error("doubles_extreme: zero (or -ve) size");
  }
  else
  {
    int i;
    result = doubles[0];
    for ( i = 1 ; i < size ; i++ )
      result = (lowest) ? 
               real_min(result,doubles[i]) : real_max(result,doubles[i]);
  }
  return(result);
}

double doubles_min(double *doubles, int size)
{
  return(doubles_extreme(doubles,size,TRUE));
}

double doubles_max(double *doubles, int size)
{
  return(doubles_extreme(doubles,size,FALSE));
}

double doubles_inner_product(double *f1, double *f2, int size)
{
  int i;
  double result = 0.0;
  for ( i = 0 ; i < size ; i++ )
    result += f1[i] * f2[i];
  return(result);
}

double doubles_magnitude_sqd(double *f, int size)
{
  return(doubles_inner_product(f,f,size));
}

void doubles_scalar_multiply(double *f, int size, double alpha, double *res_f)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
   res_f[i] = alpha * f[i];
}

void doubles_add(double *f1, double *f2, int size, double *res_f)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
   res_f[i] = f1[i] + f2[i];
}

void doubles_subtract(double *f1, double *f2, int size, double *res_f)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
   res_f[i] = f1[i] - f2[i];
}

double doubles_scaled_dsqd(double *f1, double *f2, double *scales, int size)
/* Returns scales . ( f1 - f2 ) (. = dot-product) */
{
  double result = 0.0;
  int i;
  for ( i = 0 ; i < size ; i++ )
  {
    double d = scales[i] * (f1[i] - f2[i]);
    result += d * d;
  }
  return(result);
}

double doubles_antiscaled_dsqd(double *f1, double *f2, double *scales, int size)
/* Returns SUM ( f1[i] - f2[i] ) / scales[i] */
{
  double result = 0.0;
  int i;
  for ( i = 0 ; i < size ; i++ )
  {
    double d = (f1[i] - f2[i]) / scales[i];
    result += d * d;
  }
  return(result);
}

void doubles_reciprocal(double *farr, int size, double *res_farr)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
  {
    if ( fabs(farr[i]) < 1e-20 )
      printf("doubles_reciprocal about to hit a near-zero\n");

    res_farr[i] = 1.0 / farr[i];
  }
}

void doubles_scalar_add(double *farr1, double *farr2, int size, double alpha, double *res_farr)
/*
   Any of the input arrays may share memory.
   res_farr := farr1 + alpha * farr2
*/
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    res_farr[i] = farr1[i] + alpha * farr2[i];
}

/* 2d_arrays. The entry at the ith row and jth colum is index as
   x[i][j]. The 2d arrays are implemented as an array of pointers
   to arrays of doubles. indexing is [0..rows-1] [0..columns-1] inclusive
*/

double **am_malloc_2d_realnums(int rows, int cols)
{
  double **result = AM_MALLOC_ARRAY(double_ptr,rows);
  int i;

  for ( i = 0 ; i < rows ; i++ )
    result[i] = am_malloc_realnums(cols);

  return(result);
}

void doubles_2d_scalar_add(double **tdarr1, double **tdarr2, int rows,
 int cols, double alpha, double **res_tdarr)
/*
   Any of the input arrays may share memory.
   res_tdarr := tdarr1 + alpha * tdarr2
*/
{
  int i;
  for ( i = 0 ; i < rows ; i++ )
    doubles_scalar_add(tdarr1[i],tdarr2[i],cols,alpha,res_tdarr[i]);
}

void fprintf_2d_realnums(FILE *s, char *m1, double **tdarr,
 int rows, int cols, char *m2)
{
  char buff[100];
  int i;

  for ( i = 0 ; i < rows ; i++ )
  {
    sprintf(buff,"%s[%d] = ",m1,i);
    fprintf_realnums(s,buff,tdarr[i],cols,m2);
  }
}

void free_2d_realnums(double **tdarr, int rows, int cols)
{
  int i;
  for ( i = 0 ; i < rows ; i++ )
    AM_FREE_ARRAY(tdarr[i],double,cols);
  AM_FREE_ARRAY(tdarr,double_ptr,rows);
}

void am_free_2d_realnums(double **tdarr, int rows, int cols)
/* Exactly the same as "free_2d_realnums". Added for naming consistency. */
{
  free_2d_realnums(tdarr,rows,cols);
}

void set_2d_realnums_constant(double **tdarr, int rows, int cols, double val)
{
  int i;
  for ( i = 0 ; i < rows ; i++ )
    set_realnums_constant(tdarr[i],cols,val);
}

void copy_2d_realnums(double **tdarr1, double **tdarr2, int rows, int cols)
{
  int i;
  for ( i = 0 ; i < rows ; i++ )
    copy_realnums(tdarr1[i],tdarr2[i],cols);
}

double doubles_sum(double *farr, int size)
{
  double result = 0.0;
  int i;
  for ( i = 0 ; i < size ; i++ )
    result += farr[i];
  return(result);
}

double doubles_mean(double *farr, int size)
{
  return(doubles_sum(farr,size) / int_max(size,1));
}

double doubles_sdev(double *farr, int size)
{
  double mean = doubles_mean(farr,size);
  int i;
  double sum_sqs = 0.0;
  for ( i = 0 ; i < size ; i++ )
    sum_sqs += real_square(farr[i] - mean);

  return(sqrt(sum_sqs / int_max(size,1)));
}

void random_shuffle_ints(int *iarr, int size)
{
  int i;
  for ( i = 0 ; i < size-1 ; i++ )
  {
    int j = i + int_random(size-i);
    int temp = iarr[i];
    iarr[i] = iarr[j];
    iarr[j] = temp;
  }
}


/* see __cdecl note at int_comp */
#ifdef WIN32
int __cdecl double_comp(const void *r1_ptr, const void *r2_ptr)
#else
int double_comp(const void *r1_ptr, const void *r2_ptr)
#endif
{
  double *r1 = (double *) r1_ptr;
  double *r2 = (double *) r2_ptr;
  if ( *r1 < *r2 )
    return(-1);
  else if ( *r1 > *r2 )
    return(1);
  else
    return(0);
}

/* PRE:  farr and r_farr must be different bits of memory.
         indices must be a permutation of 0,1,...size-1 
   POST: r_farr[i] = farr[indices[i]] */
void realnums_permute(double *farr,int *indices,double *r_farr,int size)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    r_farr[i] = farr[indices[i]];
}

/* PRE:  iarr and r_iarr must be different bits of memory.
         indices must be a permutation of 0,1,...size-1 
   POST: r_iarr[i] = iarr[indices[i]] */
void ints_permute(int *iarr,int *indices,int *r_iarr,int size)
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    r_iarr[i] = iarr[indices[i]];
}

/*********** SPR's own QuickSort implementation.

    AWM: It would usually be insane to use anything except the OS's
    built in sort but it so happens that microsoft's qsort goes
    slowly when there are lots of duplicates. The following implementation
    behaves very well with many duplicates. 

    Later on you'll see that our bottom level sorting chooses whether or
    not to use the OS sort according to whether there are many duplicates
    or not
*/

#define my_qs_swap(indices,i,j) \
{ \
  int temp = indices[i]; \
  indices[i] = indices[j]; \
  indices[j] = temp; \
}

void my_qs_2(double *farr, int start, int *indices)
{
  double x0 = farr[indices[start+0]];
  double x1 = farr[indices[start+1]];
  if ( x1 < x0 )
    my_qs_swap(indices,start,start+1);
}
  
void my_qs(double *farr, int size, int start, int end, int *indices)
{
  if ( end <= start+1 )
  {
    /* skip */
  }
  else if ( end == start + 2 )
  {
    my_qs_2(farr,start,indices);
  }
  else if ( end == start + 3 )
  {
    double x0 = farr[indices[start+0]];
    double x1 = farr[indices[start+1]];
    double x2 = farr[indices[start+2]];
    if ( x0 > x1 || x0 > x2 )
    {
      if ( x1 < x2 )
      {
        my_qs_swap(indices,start,start+1);
        my_qs_2(farr,start+1,indices);
      }
      else
      {
        my_qs_swap(indices,start,start+2);
        my_qs_2(farr,start+1,indices);
      }
    }
    else
      my_qs_2(farr,start+1,indices);
  }
  else
  {
    int pivot_index = start + int_random(end-start);
    double pivot_value = farr[indices[pivot_index]];
    int i = start;
    int j = start;
    int k = end;

    while ( j < k )
    {
      /* Assert: indices[start] ... indices[i-1] point to values < pivot_value
                 indices[i] .. indices[j-1] point to values == pivot_value
                 indices[k] .. indices[end-1] point to values > pivot_value
           indices[start] .. indices[end-1] are a permutation of the indices on the
           way in */

      while ( j < k && farr[indices[j]] <= pivot_value )
      {
        while ( j < k && farr[indices[j]] < pivot_value )
        {
          if ( i < j )
            my_qs_swap(indices,i,j);
          i += 1;
          j += 1;
        }

        while ( j < k && farr[indices[j]] == pivot_value )
        {
          if ( i < j )
            my_qs_swap(indices,i,j);
          j += 1;
        }

      }

      while ( j < k && farr[indices[k-1]] > pivot_value )
        k -= 1;

      if ( j < k && farr[indices[k-1]] == pivot_value )
      {
        my_qs_swap(indices,j,k-1);
        j += 1;
        if ( j < k )
          k -= 1;
      }
      else if ( j < k )
      {
        if ( i == j )
        {
          my_qs_swap(indices,i,k-1);
          i += 1;
          j += 1;
          if ( j < k )
            k -= 1;
        }
        else
        {
          my_qs_swap(indices,j,k-1);
          my_qs_swap(indices,i,j);
          i += 1;
          j += 1;
          if ( j < k )
            k -= 1;
        }
      }

    }

    my_qs(farr,size,start,i,indices);
    my_qs(farr,size,k,end,indices);
  }
}

bool indices_are_sorted(double *farr,int *indices,int size)
{ /* XXX: This would be more efficient if the && were removed, and the "sorted =" were
     replaced with sorted *=     */
  bool sorted = TRUE;
  int i;
  for ( i = 0 ; sorted && i < size-1 ; i++ )
    sorted = farr[indices[i]] <= farr[indices[i+1]];
  return sorted;
}

void spr_indices_sort_realnums(double *farr, int size, int *indices)
/*
   PRE: farr and "indices" both have "size" elements.
        indices's contents irrelevant.
   POST: indices contains sorted indiceses into farr, so that
         forall j, indices[j] is the j'th smallest element of farr.

         thus forall i,j, (i < j) => farr[indices[i]] <= farr[indices[j]]
          and indices contains a permutation of [0 ...size-1]
*/
{
  int i;
  for ( i = 0 ; i < size ; i++ )
    indices[i] = i;

  /* All the following am_srand stuff is because the above 
     implementation of quicksort randomly selects its pivots and
     (a) I want to do the same pivot selection each time I sort
     the same array and (b) I want the random number generator
     state to be the same after a sort no matter what kind of sort
     I used. */
      
  push_current_am_srand_state();
  am_srand(575297);
  my_qs(farr,size,0,size,indices);
  pop_current_am_srand_state();

#ifndef AMFAST
  if ( !indices_are_sorted(farr,indices,size) )
    my_error("indices_sort_realnums failed");
#endif
}

void spr_sort_realnums(double *farr, int size, double *r_farr)
/*
   It is fine for farr to be the same memory as r_farr
*/
{
  int *indices = AM_MALLOC_ARRAY(int,size);
 
  spr_indices_sort_realnums(farr,size,indices);

  if ( farr == r_farr )
  {
    double *temp = AM_MALLOC_ARRAY(double,size);
    realnums_permute(farr,indices,temp,size);
    copy_realnums(temp,r_farr,size);
    AM_FREE_ARRAY(temp,double,size);
  }
  else
    realnums_permute(farr,indices,r_farr,size);

  AM_FREE_ARRAY(indices,int,size);
}

void os_sort_realnums(double *farr, int size, double *r_farr)
/*
   It is fine for farr to be the same memory as r_farr
*/
{
  copy_realnums(farr,r_farr,size);
  qsort((char *)r_farr,
        size,
        sizeof(double),
	double_comp
       );
}

double doubles_median(double *farr, int size)
{
  double *f2 = am_malloc_realnums(size);
  int med_index = size/2;
  double result;

  sort_realnums(farr,size,f2);

  if ( (size%2) == 0 )
    result = f2[med_index];
  else
    result = (f2[med_index] + f2[med_index+1])/2.0;

  AM_FREE_ARRAY(f2,double,size);

  return(result);
}

double *am_malloc_realnums(
    int size
  )
{
  double *result = AM_MALLOC_ARRAY(double,size);
#ifndef AMFAST
  set_realnums_constant(result,size,-7.7777e27);
#endif
  return(result);
}

void am_free_realnums(
  double *farr,
  int size
  )
{
  AM_FREE_ARRAY(farr,double,size);
}

int *am_malloc_ints(
    int size
  )
{
  int *result = AM_MALLOC_ARRAY(int,size);
  return(result);
}

void am_free_ints(
  int *iarr,
  int size
  )
{
  AM_FREE_ARRAY(iarr,int,size);
}

long *am_malloc_longs(
    int size
  )
{
  long *result = AM_MALLOC_ARRAY(long,size);
  return(result);
}

void am_free_longs(
  long *larr,
  int size
  )
{
  AM_FREE_ARRAY(larr,long,size);
}

bool *am_malloc_bools(
    int size
  )
{
  bool *result = AM_MALLOC_ARRAY(bool,size);
  return(result);
}

void am_free_bools(
  bool *barr,
  int size
  )
{
  AM_FREE_ARRAY(barr,bool,size);
}

int ints_argextreme(int *ints, int size, bool lowest)
{
  int i, r=0;
  int result;

  if (size <= 0)
    my_error("int_argextreme: zero (or -ve) size");
 
  result=ints[0];
  for (i=1; i<size; i++) {
    if ( lowest ? (ints[i] < result) : (ints[i] > result)) {
      r=i;
      result=ints[i];
    }
  }
  return r;
}

int ints_argmin(int *ints, int size)
{
  return ints_argextreme(ints,size,TRUE);
}

int ints_argmax(int *ints, int size)
{
  return ints_argextreme(ints,size,FALSE);
}

int doubles_argextreme(double *doubles, int size, bool lowest)
{
  int i, r=0;
  double result;
  if (size <= 0)
    my_error("double_argextreme: zero (or -ve) size");

  result = doubles[0];
  for (i=1; i<size; i++) {
    if ( lowest ? (doubles[i] < result) : (doubles[i] > result)) {
      r=i;
      result=doubles[i];
    }
  }
  return r;
}

int doubles_argmin(double *doubles, int size)
{
  return doubles_argextreme(doubles, size, TRUE);
}

int doubles_argmax(double *doubles, int size)
{
  return doubles_argextreme(doubles, size, FALSE);
}

/* indices_sort_realnums Courtesy of Justin Boyan.... */

typedef struct isf_elt_struct { int ix; double v; } isf_elt;

/* see __cdecl note at int_comp */
#ifdef WIN32
int __cdecl compare_isf_elt(const void *i_ptr, const void *j_ptr)
#else
int compare_isf_elt(const void *i_ptr, const void *j_ptr)
#endif
{
  isf_elt *i = (isf_elt *) i_ptr;
  isf_elt *j = (isf_elt *) j_ptr;
  return (i->v < j->v) ? -1 : (i->v > j->v) ? 1 : 0;
}

void os_indices_sort_realnums(double *farr, int size, int *indices)
/*
   PRE: farr and "indices" both have "size" elements.
        indices's contents irrelevant.
   POST: indices contains sorted indiceses into farr, so that
         forall j, indices[j] is the j'th smallest element of farr.

         thus forall i,j, (i < j) => farr[indices[i]] <= farr[indices[j]]
          and indices contains a permutation of [0 ...size-1]
*/
{
  int i;
  isf_elt *isf = AM_MALLOC_ARRAY(isf_elt, size);

  for (i=0; i<size; i++) {
    isf[i].ix = i;
    isf[i].v = farr[i];
  }
  qsort((char *)isf, size, sizeof(isf_elt), compare_isf_elt);
  for (i=0; i<size; i++) {
    indices[i] = isf[i].ix;
  }
  AM_FREE_ARRAY(isf, isf_elt, size);
}

int num_different_values_in_sorted_farr(double *farr,int size)
{
  int result = 1;
  int i;
  for ( i = 1 ; i < size ; i++ )
  {
    if ( farr[i] > farr[i-1] )
      result += 1;
  }
  return result;
}

int frequency_of_most_common_value_in_sorted_farr(double *farr,int size)
{
  int result = 1;
  int i;
  int count = 1;
  for ( i = 1 ; i < size ; i++ )
  {
    if ( farr[i] == farr[i-1] )
    {
      count += 1;
      result = int_max(result,count);
    }
    else
      count = 1;
  }
  return result;
}

bool avoid_os_sort_given_subset(double *subset,int subsize)
{
  bool result;
  os_sort_realnums(subset,subsize,subset);
  result = num_different_values_in_sorted_farr(subset,subsize) < 40 ||
           frequency_of_most_common_value_in_sorted_farr(subset,subsize) > subsize/2;
  if ( result )
    printf("Avoiding Microsoft's sort because of poor performance with duplicates\n");

  return result;
}

bool avoid_os_sort(double *farr,int size)
{
  bool result;
#ifndef PC_PLATFORM
  result = FALSE;
#else
  if ( size < 5000 )
    result = FALSE;
  else
  {
    int subsize = 200;
    double *subset = AM_MALLOC_ARRAY(double,subsize);
    int i;
    for ( i = 0 ; i < subsize ; i++ )
      subset[i] = farr[int_random(size)];
    result = avoid_os_sort_given_subset(subset,subsize);
    AM_FREE_ARRAY(subset,double,subsize);
  }
#endif

  return result;
}

bool avoid_os_int_sort(int *iarr,int size)
{
  bool result;
#ifndef PC_PLATFORM
  result = FALSE;
#else
  if ( size < 5000 )
    result = FALSE;
  else
  {
    int subsize = 200;
    double *subset = AM_MALLOC_ARRAY(double,subsize);
    int i;
    for ( i = 0 ; i < subsize ; i++ )
      subset[i] = (double) iarr[int_random(size)];
    result = avoid_os_sort_given_subset(subset,subsize);
    AM_FREE_ARRAY(subset,double,subsize);
  }
#endif
  return result;
}

/* by Artur ----------------------------------------------------------------- */
void os_indices_sort_integers(int *iarr, int size, int *indices)
/*
   PRE: iarr and "indices" both have "size" elements.
        indices's contents irrelevant.
   POST: indices contains sorted indiceses into iarr, so that
         forall j, indices[j] is the j'th smallest element of iarr.

         thus forall i,j, (i < j) => iarr[indices[i]] <= iarr[indices[j]]
          and indices contains a permutation of [0 ...size-1]
*/
{
  int i;
  isf_elt *isf = AM_MALLOC_ARRAY(isf_elt, size);
  for (i=0; i<size; i++) {
    isf[i].ix = i;
    isf[i].v = iarr[i];
  }
  qsort((char *)isf, size, sizeof(isf_elt), compare_isf_elt);
  for (i=0; i<size; i++) {
    indices[i] = isf[i].ix;
  }
  AM_FREE_ARRAY(isf, isf_elt, size);
}

double *mk_farr_from_iarr(int *iarr,int size)
{
  double *farr = AM_MALLOC_ARRAY(double,size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    farr[i] = (double) iarr[i];
  return farr;
}

void spr_indices_sort_integers(int *iarr,int size,int *indices)
{
  double *farr = mk_farr_from_iarr(iarr,size);
  spr_indices_sort_realnums(farr,size,indices);
  AM_FREE_ARRAY(farr,double,size);
}

void spr_sort_ints(int *iarr, int size, int *r_iarr)
/*
   It is fine for iarr to be the same memory as r_iarr
*/
{
  int *indices = AM_MALLOC_ARRAY(int,size);
  double *farr = mk_farr_from_iarr(iarr,size);
  spr_indices_sort_realnums(farr,size,indices);

  if ( iarr == r_iarr )
  {
    int *temp = AM_MALLOC_ARRAY(int,size);
    ints_permute(iarr,indices,temp,size);
    copy_ints(temp,r_iarr,size);
    AM_FREE_ARRAY(temp,int,size);
  }
  else
    ints_permute(iarr,indices,r_iarr,size);

  AM_FREE_ARRAY(indices,int,size);
  AM_FREE_ARRAY(farr,double,size);
}

/*
   PRE: farr and "indices" both have "size" elements.
        indices's contents irrelevant.
   POST: indices contains sorted indiceses into farr, so that
         forall j, indices[j] is the j'th smallest element of farr.

         thus forall i,j, (i < j) => farr[indices[i]] <= farr[indices[j]]
          and indices contains a permutation of [0 ...size-1]
*/
void indices_sort_realnums(double *farr, int size, int *indices)
{
  if ( avoid_os_sort(farr,size) )
    spr_indices_sort_realnums(farr,size,indices);
  else
    os_indices_sort_realnums(farr,size,indices);
}

void sort_realnums(double *farr, int size, double *r_farr)
{
  if ( avoid_os_sort(farr,size) )
    spr_sort_realnums(farr,size,r_farr);
  else
    os_sort_realnums(farr,size,r_farr);
}

/*
   PRE: farr and "indices" both have "size" elements.
        indices's contents irrelevant.
   POST: indices contains sorted indiceses into farr, so that
         forall j, indices[j] is the j'th smallest element of farr.

         thus forall i,j, (i < j) => farr[indices[i]] <= farr[indices[j]]
          and indices contains a permutation of [0 ...size-1]
*/
void indices_sort_integers(int *iarr, int size, int *indices)
{
  if ( avoid_os_int_sort(iarr,size) )
    spr_indices_sort_integers(iarr,size,indices);
  else
    os_indices_sort_integers(iarr,size,indices);
}

void sort_ints(int *iarr, int size, int *r_iarr)
{
  if ( avoid_os_int_sort(iarr,size) )
    spr_sort_ints(iarr,size,r_iarr);
  else
    os_sort_ints(iarr,size,r_iarr);
}







