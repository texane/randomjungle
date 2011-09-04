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
   File:        amdyv.c
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:13 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Header for Dynamically allocated and deallocated vectors

   Copyright 1996, Schenley Park Research
*/

#include <stdio.h>
#include <math.h>

#include "amdyv.h"     /* Dynamically allocated and deallocated vectors */
#include "amma.h"      /* Fast, non-fragmenting, memory management */
#include "amar.h"      /* Obvious operations on 1-d arrays */
#include "am_string.h"
#include "am_string_array.h"



int Dyvs_mallocked = 0;
int Dyvs_freed = 0;

dyv *mk_dyv(int size)
{
  dyv *result = AM_MALLOC(dyv);
  result -> dyv_code = DYV_CODE;
  result -> array_size = size;
  result -> size = size;
  result -> farr = am_malloc_realnums(size);
  Dyvs_mallocked += 1;
  return(result);
}


/* Acts as though you created a dyv of size 0,         */
/* but actually allocates an array of size "capacity". */
/* This is very useful if you want to use add_to_dyv   */
/* and happen to have a reasonable upper bound to the  */
/* number of elements you want to add.                 */
dyv *mk_empty_dyv(int capacity) {
  dyv *result = AM_MALLOC(dyv);
  if ( capacity < 0 ) my_error("mk_empty_ivec : capacity < 0 illegal");
  result -> dyv_code = DYV_CODE;
  result -> array_size = capacity;
  result -> size = 0;
  result -> farr = am_malloc_realnums(capacity);
  Dyvs_mallocked += 1;
  return(result);
}


dyv *mk_dyv_x( int size, ...)
{
  /* Warning: no type checking can be done by the compiler.  You *must*
     send the values as doubles for this to work correctly. */
  int i;
  double val;
  va_list argptr;
  dyv *dv;
  
  dv = mk_dyv( size);

  va_start( argptr, size);
  for (i=0; i<size; ++i) {
    val = va_arg( argptr, double);
    dyv_set( dv, i, val);
  }
  va_end(argptr);

  return dv;
}

double dyv_sum(const dyv *dv)
{
  double result = 0.0;
  int i;
  for ( i = 0 ; i < dyv_size(dv) ; i++ )
    result += dyv_ref(dv,i);
  return(result);
}

double dyv_product(const dyv *dv)
{
  double result = 1.0;
  int i;
  for ( i = 0 ; i < dyv_size(dv) ; i++ )
    result *= dyv_ref(dv,i);
  return(result);
}

/* Returns 1 if all elements are 0.0 or 1.0.  Note that 0.0 and 1.0 are
   perfectly representable in IEEE 754. */
int dyv_is_binary( const dyv *dv)
{
  int size, i;
  double val;
  size = dyv_size( dv);
  for (i=0; i<size; ++i) {
    val = dyv_ref( dv, i);
    if (val != 0.0 && val != 1.0) return 0;
  }
  return 1;
}

void free_dyv(dyv *d)
{
  d -> dyv_code = 7777;

  am_free_realnums(d->farr,d->array_size);
  AM_FREE(d,dyv);

  Dyvs_freed += 1;
}

void add_to_dyv(dyv *d,double new_val)
{
  if ( d->array_size < 0 || d->size > d->array_size )
    my_error("dyv size or array size has got muddles. Talk to AWM");

  if ( d->size == d->array_size )
  {
    int new_array_size = 2 * d->size + 2;
    double *farr_new = AM_MALLOC_ARRAY(double,new_array_size);
    int i;
    for ( i = 0 ; i < d->size ; i++ )
      farr_new[i] = d->farr[i];
    AM_FREE_ARRAY(d->farr,double,d->size);
    d -> farr = farr_new;
    d -> array_size = new_array_size;
  }
  d->farr[d->size] = new_val;
  d -> size += 1;
}

double *mk_farr_from_dyv(const dyv *d)
{
  double *result;
  result = am_malloc_realnums(d->size);
  copy_dyv_to_farr(d,result);
  return(result);
}

void copy_dyv_to_farr(const dyv *d, double *farr)
{
  copy_realnums(d->farr,farr,d->size);
}
  
void copy_farr_to_dyv(double *farr,int size,dyv *r_d)
{
  copy_realnums(farr,r_d->farr,size);
}

dyv *mk_dyv_from_farr(double *farr,int size)
{
  dyv *result = mk_dyv(size);
  copy_farr_to_dyv(farr,size,result);
  return(result);
}

void constant_dyv(dyv *r_d,double v)
{
  set_realnums_constant(r_d->farr,r_d->size,v);
}

void zero_dyv(dyv *r_d)
{
  constant_dyv(r_d,0.0);
}

dyv *mk_constant_dyv(int size,double v)
{
  dyv *result = mk_dyv(size);
  constant_dyv(result,v);
  return(result);
}

dyv *mk_zero_dyv(int size)
{
  dyv *result = mk_dyv(size);
  zero_dyv(result);
  return(result);
}

/********* Standard operations of dyvs *********/

void fprint_dyv_csv(FILE *s,dyv *x)
{
  int i;
  int size = dyv_size(x);

  for ( i = 0 ; i < size ; i++ )
    fprintf(s,"%g%s",dyv_ref(x,i),(i==size-1)?"\n":",");
}

void fprintf_dyv(FILE *s, const char *m1, const dyv *d, const char *m2)
{
  if ( d == NULL )
    fprintf(s,"%s = (dyv *)NULL%s",m1,m2);
  else if ( d->dyv_code != DYV_CODE )
  {
    fprintf(stderr,"fprintf_dyv(s,\"%s\",d,\"\\n\"\n",m1);
    my_error("fprintf_dyv called with a non-allocated dyv (DYnamic Vector)");
  }
  else if ( d->size <= 0 )
    fprintf(s,"%s = <Dyv of size %d>%s",m1,d->size,m2);
  else
  {
    int i;
    buftab bt;
    int cols = 1;

    init_buftab(&bt,d->size,cols + 4);

    for ( i = 0 ; i < d->size ; i++ )
    {
      char buff[100];
      set_buftab(&bt,i,0,(i == (d->size-1)/2) ? m1 : "");
      set_buftab(&bt,i,1,(i == (d->size-1)/2) ? "=" : "");
      set_buftab(&bt,i,2,"(");

      sprintf(buff," %g ",d->farr[i]);
      set_buftab(&bt,i,3,buff);

      set_buftab(&bt,i,3+cols,")");
    }

    fprint_buftab(s,&bt);
    free_buftab_contents(&bt);
  }
  fprintf(s,"\n");
}

/* Reduces the size of d by one.
   dyv_ref(d,index) disappears.
   Everything to the right of dyv_ref(d,index) is copied one to the left.

   Formally: Let dold be the dyv value beore calling this function
             Let dnew be the dyv value after calling this function.

PRE: dyv_size(dold) > 0
     0 <= index < dyv_size(dold)

POST: dyv_size(dnew) = dyv_size(dold)-1
      for j = 0 , 1, 2 ... index-1  : 
         dyv_ref(dnew,j) == dyv_ref(dold,j)

      for j = index , index+1 , ... dyv_size(dnew)-1:
         dyv_ref(dnew,j) == dyv_ref(dold,j+1)
*/
void dyv_remove(dyv *d,int index)
{
  int i;
  int dsize = dyv_size(d);

#ifndef AMFAST
  if ( dsize <= 0 ) my_error("dyv_remove: empty dyv");
  if ( index < 0 || index >= dsize ) my_error("dyv_remove: bad index");
#endif /* #ifndef AMFAST */

  for ( i = index ; i < dsize - 1 ; i++ )
    dyv_set(d,i,dyv_ref(d,i+1));
  d -> size -= 1;
}

void dyv_mult(const dyv *d1, const dyv *d2, dyv *rd)
{
  int i;
  int n = d1->size;
  for (i=0;i<n;i++) rd->farr[i] = d1->farr[i] * d2->farr[i];
}

dyv *mk_dyv_mult(const dyv *d1, const dyv *d2)
{
  dyv *result;
  result = mk_dyv(d1->size);
  dyv_mult(d1,d2,result);
  return result;
}

void dyv_scalar_mult(const dyv *d, double alpha, dyv *r_d)
{
  int i, n = d->size;
  for ( i = 0 ; i < n ; i++ )
    r_d -> farr[i] = d->farr[i] * alpha;
}

dyv *mk_dyv_scalar_mult(const dyv *d, double alpha)
{
  dyv *result;
  result = mk_dyv(d->size);
  dyv_scalar_mult(d,alpha,result);
  return(result);
}

void dyv_scalar_add(dyv *d, double alpha, dyv *r_d)
{
  int i;
  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d->farr[i] + alpha;
}

dyv *mk_dyv_scalar_add(dyv *d,double alpha)
{
  dyv *result;
  result = mk_dyv(d->size);
  dyv_scalar_add(d,alpha,result);
  return(result);
}

void dyv_madd( double factor, dyv *dv1, dyv *dv2, dyv *r_dv)
{
  /* This functions hopes to take advantage of the popular madd instruction. */
  int i;
  double v1, v2, result;

  for (i=0; i<dyv_size(dv1); ++i) {
    v1 = dyv_ref( dv1, i);
    v2 = dyv_ref( dv2, i);
    result = factor*v1 + v2;
    dyv_set( r_dv, i, result);
  }

  return;
}

dyv *mk_dyv_madd( double factor, dyv *dv1, dyv *dv2)
{
  dyv *result;
  result = mk_dyv( dyv_size(dv1));
  dyv_madd( factor, dv1, dv2, result);
  return result;
}


void copy_dyv(const dyv *d, dyv *r_d)
{
  dyv_scalar_mult(d,1.0,r_d);
}
    
dyv *mk_copy_dyv(const dyv *d)
{
  return(mk_dyv_scalar_mult(d,1.0));
}

/* Returns dyv of size ivec_size(rows) in which
    result[i] = x[rows[i]] */
dyv *mk_dyv_subset( const dyv *x, const ivec *rows)
{
  int size, i;
  dyv *y;

  size = ivec_size( rows);
  y = mk_dyv(size);
  
  for (i=0; i<size; i++) dyv_set( y, i, dyv_ref(x, ivec_ref(rows,i)));

  return y;
}

void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d)
{
  int i;
  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_dyv(stderr,"d_1",d_1,"\n");
    fprintf_dyv(stderr,"d_2",d_2,"\n");
    my_error("dyv_plus: dyvs (DYnamic Vectors) different shape");
  }

  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d_1->farr[i] + d_2 -> farr[i];
}

dyv *mk_dyv_plus(const dyv *a,const dyv *b)
{
  dyv *result = mk_dyv(a->size);
  dyv_plus(a,b,result);
  return(result);
}

void dyv_subtract(const dyv *d_1,const dyv *d_2,dyv *r_d)
{
  int i;
  if ( d_1 -> size != d_2 -> size )
  {
    fprintf_dyv(stderr,"d_1",d_1,"\n");
    fprintf_dyv(stderr,"d_2",d_2,"\n");
    my_error("dyv_subtract: dyvs (DYnamic Vectors) different shape");
  }

  for ( i = 0 ; i < r_d -> size ; i++ )
    r_d -> farr[i] = d_1->farr[i] - d_2 -> farr[i];
}

dyv *mk_dyv_subtract(const dyv *a,const dyv *b)
{
  dyv *result; 
  result = mk_dyv(a->size);
  dyv_subtract(a,b,result);
  return(result);
}

void dyv_abs( const dyv *dv, dyv *r_dv)
{
  int i;
  for (i=0; i<dyv_size( dv); ++i) dyv_set( r_dv, i, fabs(dyv_ref( dv, i)));
  return;
}

dyv *mk_dyv_abs( const dyv *dv)
{
  dyv *absdv;
  absdv = mk_dyv( dyv_size( dv));
  dyv_abs( dv, absdv);
  return absdv;
}

/***** More complex operations ******/

double dyv_scalar_product(const dyv *a,const dyv *b)
{
  int i,n= a->size;
  double result = 0.0;
  if ( b -> size != n)
  {
    fprintf(stderr,
            "dyv_scalar_product: sizes differ (e.g. model vec different size\n"
            "    than data point)\n");
    fprintf_dyv(stderr,"a",a,"\n");
    fprintf_dyv(stderr,"b",b,"\n");
    my_error("dyv_scalar_product: sizes differ (e.g. model vec different\n"
             "    size than data point)\n");
  }

  for ( i = 0 ; i < n ; i++ )
    result += a->farr[i] * b->farr[i];
  return(result);
}

double dyv_magnitude(const dyv *a)
{
  double result = sqrt(dyv_scalar_product(a,a));
  return(result);
}

double dyv_mean(const dyv *dv)
{
  return(doubles_mean(dv->farr,dv->size));
}

double dyv_sdev(const dyv *dv)
{
  double sum_sq = 0.0;
  int i;
  double mean = dyv_mean(dv);
  double result;

  for ( i = 0 ; i < dv->size ; i++ )
    sum_sq += real_square(dv->farr[i] - mean);
    
  result = sqrt(sum_sq / int_max(dv->size-1,1));

  return(result);
}

double dyv_min(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_min: empty dyv");
  return(doubles_min(dv->farr,dv->size));
}

double dyv_max(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_max: empty dyv");
  return(doubles_max(dv->farr,dv->size));
}

int dyv_argmin(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_argmin: empty dyv");
  return(doubles_argmin(dv->farr,dv->size));
}

int dyv_argmax(const dyv *dv)
{
  if ( dyv_size(dv) < 1 )
    my_error("dyv_argmax: empty dyv");
  return(doubles_argmax(dv->farr,dv->size));
}

/* Returns first index i such that dv[i] >= key.  If the last element of
   dv is strictly less than the key, then dyv_size(dv) is returned. */
int index_in_sorted_dyv( dyv *dv, double key)
{
  int idxlo, idxhi, idx;
  double val;

  /* idxhi is always >= desired idx.  idxlo is always <= desired idx. */

  idxlo = 0;
  idxhi = dyv_size( dv);
  if (idxhi == 0) return 0;

  while (1) {
    /* Get next value. */
    idx = (idxlo + idxhi) / 2;
    val = dyv_ref( dv, idx);

    /* Adjust bounds. */
    if (val >= key) idxhi = idx ;  /* idx is too high. */
    else idxlo = idx + 1;          /* idx is too low.  */

    /* If the boundaries meet, we're done. */
    if (idxhi == idxlo) break;
  }

  return idxlo;
}


/* We're not quite sure what this code computed.  On the bright
   side, since it isn't commented we can easily comment the
   entire function out. */
/* 
int index_in_sorted_dyv(dyv *d,double t){
  int i1 = dyv_size(d), i0 = 0;
  int i = (i1-i0)>>1;
  double n2 = dyv_ref(d,i), n1 = (i)? dyv_ref(d,i-1):-1;
  while((i1-i0)>1 && !(n2>=t && (!i || n1<t))){
    if(n2>t){
      i1 = i;
      i = i0+((i-i0)>>1);
    } else {
      i0 = i+1;
      i = i+((i1-i)>>1);
    }
  }
  return i;
}
*/

void dyv_sort(const dyv *dv,dyv *r_dv)
{
  int size;
  double *farr;

  size = dyv_size(dv);
  farr = mk_farr_from_dyv(dv);
  sort_realnums(farr,size,farr);
  copy_farr_to_dyv(farr,size,r_dv);
  am_free_realnums(farr,size);
}


dyv *mk_dyv_sort(const dyv *dv)
{
  dyv *result;
  result = mk_dyv(dyv_size(dv));
  dyv_sort(dv,result);
  return(result);
}

/*
 Creates a dyv of indices such that indices[i] is the origional
 location (in the unsorted dv) of the ith smallest value.
 Used when you want the location of the sorted values instead of
 the sorted vector itself.
*/
dyv *mk_sorted_dyv_indices(dyv *dv)
{
  dyv* indices;
  int size = dyv_size(dv);
  int i;
  int* iarr; 
  double* farr;

  farr = mk_farr_from_dyv(dv);
  iarr = am_malloc_ints(size); 

  indices_sort_realnums(farr,size,iarr);
  indices = mk_dyv(size);

  for(i=0;i<size;i++)
  {
    dyv_set(indices,i,iarr[i]);
  }

  am_free_realnums(farr,size);
  am_free_ints(iarr,size);

  return indices;
}

ivec *mk_ivec_sorted_dyv_indices(dyv *dv)
{
  ivec* indices;
  int size = dyv_size(dv);
  double* farr;

  farr = mk_farr_from_dyv(dv);

  indices = mk_ivec(size);
  indices_sort_realnums(farr,size,indices->iarr);

  am_free_realnums(farr,size);

  return indices;
}


/* Increases dv in length by 1 and shifts all elements
   with original index greater or equal to index one to the
   right and inserts val at index. */
void dyv_insert(dyv *dv,int index,double val)
{
  int i;
  add_to_dyv(dv,0);
  for ( i = dyv_size(dv)-1 ; i > index ; i-- )
    dyv_set(dv,i,dyv_ref(dv,i-1));
  dyv_set(dv,index,val);
}

void fprintf_oneline_dyv(FILE *s,const char *m1, const dyv *d, const char *m2)
{
  int i;
  fprintf(s,"%s ",m1);
  for ( i = 0 ; i < dyv_size(d) ; i++ )
    fprintf(s,"%8g%s",dyv_ref(d,i),(i==dyv_size(d)-1) ? "" : " ");
  fprintf(s,"%s",m2);
}






void indices_of_sorted_dyv(const dyv *dv,ivec *iv)
/*
   NOTE: ivec structure (integer vectors) defined in sortind.ch
   PRE: dv and iv must be same size. iv's contents irrelevant.
   POST: iv contains sorted indexes into dv, so that
         forall iv[j] is the j'th smallest element of dv.

         thus forall i,j, (i < j) => dyv_ref(dv,iv[i]) <= dyv_ref(dv,iv[j])
          and iv contains a permutation of [0 ... iv->size-1]
*/
{
  int *iarr = am_malloc_ints(dyv_size(dv));
  indices_sort_realnums(dv->farr,dv->size,iarr);
  copy_iarr_to_ivec(iarr,dv->size,iv);
  am_free_ints(iarr,dv->size);
}

ivec *mk_indices_of_sorted_dyv(const dyv *dv)
{
  ivec *iv = mk_ivec(dyv_size(dv));
  indices_of_sorted_dyv(dv,iv);
  return(iv);
}







double dyv_partial_sum( const dyv *dv, const ivec *indices)
{
  double sum;
  sum = 0.0;
  
# ifdef DEBUG
  int i, imax;
  imax = ivec_size(indices);
  for (i=0; i<imax; ++i) {
    if (ivec_ref(indices, i) >= dyv_size(dv)) {
      my_errorf( "dyv_partial_sum: ivec index %d at %d greater than dyv "
                 "size %d\n"
                 "    (e.g. data point has more columns than model vector)\n",
                 ivec_ref(indices, i), i, dyv_size(dv));
    }
    sum += dyv_ref( dv, ivec_ref(indices, i));
  }
# else /* NON-DEBUG */
  int *iarr, *imax;
  double *farr;
  iarr = indices->iarr;
  imax = iarr + indices->size;
  farr = dv->farr;
  while (iarr < imax) sum += farr[*iarr++];
# endif
  
  return sum;
}

dyv *mk_dyv_from_ivec(const ivec *iv)
{
  dyv *d = mk_dyv(ivec_size(iv));
  int i;
  for ( i = 0 ; i < ivec_size(iv) ; i++ )
    dyv_set(d,i, ivec_ref(iv,i));

  return(d);
}

ivec *mk_ivec_from_dyv( const dyv *dv)
{
  int size, i;
  ivec *iv;
  size = dyv_size( dv);
  iv = mk_ivec( size);
  for (i=0; i<size; ++i) ivec_set( iv, i, (int) dyv_ref(dv,i));
  return iv;
}

/* Writes dv[a] to dv[b-1] into r_dyv.  If a or b is negative, then
   it is replaced by by len+a or len+b. Returns the number of elements
   written.  This is more-or-less equivalent to Python slice syntax.
*/
int dyv_slice( const dyv *dv, int a, int b, dyv *r_dv)
{
  int dvlen, lb, ub, r_dvlen, srcidx, dstidx;
  double dval;

  /* Set lower- and (strict) upper-bounds. */
  dvlen = dyv_size( dv);
  if (a < 0) lb = int_max(0, dvlen+a);
  else lb = int_min(dvlen, a);
  if (b < 0) ub = int_max(0, dvlen+b);
  else ub = int_min(dvlen, b);
  r_dvlen = ub - lb;

  /* Copy elements. */
  for (srcidx=lb,dstidx=0;  srcidx<ub;  ++srcidx,++dstidx) {
    dval = dyv_ref( dv, srcidx);
    dyv_set( r_dv, dstidx, dval);
  }

  return r_dvlen;
}

/* If b >= a >= 0, creates "slice" of dv on [a,b-1].  If a or b < 0, then
   it is replaced by len+a or len+b.  Indices which exceed bounds are
   replaced by appropriate indices.  Therefore this function always succeeds,
   sometimes returning an empty dyv.
*/
dyv *mk_dyv_slice( const dyv *dv, int a, int b)
{
  int dvlen, lb, ub, r_dvlen;
  dyv *r_dv;

  /* Set lower- and (strict) upper-bounds. */
  dvlen = dyv_size( dv);
  if (a < 0) lb = int_max(0, dvlen+a);
  else lb = int_min(dvlen, a);
  if (b < 0) ub = int_max(0, dvlen+b);
  else ub = int_min(dvlen, b);
  r_dvlen = ub - lb;

  r_dv = mk_dyv( r_dvlen);
  dyv_slice( dv, lb, ub, r_dv);
  return r_dv;
}
