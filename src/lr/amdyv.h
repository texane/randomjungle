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
   File:        amdyv.h
   Author:      Andrew W. Moore
   Created:     Thu Sep 15 21:01:12 EDT 1994
   Updated:     amdm was split into amdyv, amdym and svd by Frank Dellaert, Aug 14 1997
   Description: Header for Dynamically allocated and deallocated vectors

   Copyright 1996, Schenley Park Research
*/

#ifndef AMDYV_H
#define AMDYV_H

#include "standard.h"
#include "ambs.h"
#include "amiv.h"

/*
%%%%%% DYVs (DYnamic Vectors)

dyvs are just one dimensional matrices. You can conveniently allocate
and free them too using the same conventions as for dynamic matrices.

dyv *mk_dyv(size) CREATES and returns one of the fellows.

int dyv_size(dyv *dv) returns the number of elements in the dyv

double dyv_ref(dv,i) returns the value of the i'th element.
  Indexing begins at 0: 0 <= i < dyv_size(dv).

All the numerical operations on dym's have counterpart numerical
operations on dyv's:

void constant_dyv(dyv *r_d,double v);
void zero_dyv(dyv *r_d);
dyv *mk_constant_dyv(int size,double v);
dyv *mk_zero_dyv(int size);
void dyv_scalar_mult(dyv *d, double alpha, dyv *r_d);
dyv *mk_dyv_scalar_mult(dyv *d,double alpha);
void dyv_scalar_add(dyv *d, double alpha, dyv *r_d);
dyv *mk_dyv_scalar_add(dyv *d,double alpha);
void copy_dyv(const dyv *d, dyv *r_d);
dyv *mk_copy_dyv(dyv *d);
void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d);
dyv *mk_dyv_plus(dyv *a,dyv *b);
void dyv_subtract(dyv *d_1,dyv *d_2,dyv *r_d);
dyv *mk_dyv_subtract(dyv *a,dyv *b);
void dyv_sort(dyv *dv,dyv *r_dv);  It is fine, as always, if dy and r_dv are
                                   the same.
dyv *mk_dyv_sort(dyv *dv);

%%%%%%%%% Making small dyvs 
dyv *mk_dyv_1(double x) makes a 1-element dyv containing x as its only element
dyv *mk_dyv_2(double x,double y) 
         makes a 2-element dyv containing x as its 0th-indexed element
                                          y as its 1-index element
dyv *mk_dyv_3(double x,double y , double z) .... obvious.

%%%%%%%%% Complex vector operations

double dyv_scalar_product(dyv *a,dyv *b);
Returns a . b

double dyv_dsqd(dyv *a,dyv *b)
Returns (a - b).(a - b)

*/


typedef struct dyv_struct
{
  int dyv_code;
  int array_size;
  int size;
  double *farr;
} dyv, *dyv_ptr;

dyv *mk_dyv(int size);

/* Acts as though you created a dyv of size 0,         */
/* but actually allocates an array of size "capacity". */
/* This is very useful if you want to use add_to_dyv   */
/* and happen to have a reasonable upper bound to the  */
/* number of elements you want to add.                 */
dyv *mk_empty_dyv(int capacity);

/* Note that this clears all the old data your dyv had */
void dyv_destructive_resize(dyv* dv, int size);

/* Warning: no type checking can be done by the compiler.  You *must*
   send the values as doubles for this to work correctly. */
dyv *mk_dyv_x(int size, ...);

void free_dyv(dyv *d);
int dyv_num_bytes(dyv *d);

void free_dyv_NaN_ok(dyv *d);

void dyv_malloc_report(void);

int safe_dyv_size(const dyv *d);
double safe_dyv_ref(const dyv *d, int i);
void safe_dyv_set(dyv *d,int i,double value);
void safe_dyv_increment(dyv *d,int i,double value);

#define dyv_size(d) ((d)->size)
#define dyv_ref(d,i) ((d)->farr[i])
#define dyv_set(d,i,v) ((d)->farr[i] = (v))
#define dyv_increment(d,i,v) ((d)->farr[i] += (v))


/*Added by Dan: Something I've wanted for a LONG time!*/
#define dyv_array_ref_ref(d,i,j) dyv_ref(dyv_array_ref(d,i),j)
#define dyv_array_ref_set(d,i,j,x) dyv_set(dyv_array_ref(d,i),j,x)
#define dyv_array_ref_size(d,i) dyv_size(dyv_array_ref(d,i))

#define add_to_dyv_array_ref(da,i,x) add_to_dyv(dyv_array_ref(da,i),x)

void dyv_increase_length(dyv *d,int extra_size);
void dyv_decrease_length(dyv *d,int new_length);
void add_to_dyv(dyv *d,double new_val);

void dyv_remove(dyv *d,int index); /* Reduces size by 1, removes index'th 
                                      element. All elements to right of
                                      delete point copied one to left.
                                      See comments in amdm.c more details */
void dyv_remove_last_element(dyv *d); /* Reduce size by 1, remove 
                                        last rightmost elt */

/* Increases dv in length by 1 and shifts all elements
   with original index greater or equal to index one to the
   right and inserts val at index. */
void dyv_insert(dyv *dv,int index,double val);

void copy_dyv_to_farr(const dyv *d, double *farr);

double *mk_farr_from_dyv(const dyv *d);

void copy_farr_to_dyv(double *farr,int size,dyv *r_d);

dyv *mk_dyv_from_farr(double *farr,int size);

void copy_dyv_to_tdarr_row(dyv *dv,double **tdarr,int row);

void copy_dyv_to_tdarr_col(dyv *dv,double **tdarr,int col);

void copy_tdarr_row_to_dyv(double **tdarr,dyv *dv,int row);

dyv *mk_dyv_from_tdarr_row(double **tdarr,int row,int tdarr_cols);

void copy_tdarr_col_to_dyv(double **tdarr,dyv *dv,int col);

dyv *mk_dyv_from_tdarr_col(double **tdarr,int col,int tdarr_rows);


void constant_dyv(dyv *r_d,double v);

void zero_dyv(dyv *r_d);

bool zero_dyvp(dyv *d);

dyv *mk_constant_dyv(int size,double v);

dyv *mk_zero_dyv(int size);

void dyv_mult(const dyv *d1, const dyv *d2, dyv *rd);
dyv *mk_dyv_mult(const dyv *d1, const dyv *d2);

void dyv_scalar_mult(const dyv *d, double alpha, dyv *r_d);

dyv *mk_dyv_scalar_mult(const dyv *d,double alpha);

void dyv_scalar_add(dyv *d, double alpha, dyv *r_d);

dyv *mk_dyv_scalar_add(dyv *d,double alpha);

void dyv_madd( double factor, dyv *dv1, dyv *dv2, dyv *r_dv);

dyv *mk_dyv_madd( double factor, dyv *dv1, dyv *dv2); 

void copy_dyv(const dyv *d, dyv *r_d);

dyv *mk_copy_dyv(const dyv *d);

dyv *mk_dyv_subset( const dyv *x, const ivec *rows);

void dyv_plus(const dyv *d_1, const dyv *d_2, dyv *r_d);

dyv *mk_dyv_plus(const dyv *a,const dyv *b);

void dyv_subtract(const dyv *d_1,const dyv *d_2,dyv *r_d);

dyv *mk_dyv_subtract(const dyv *a,const dyv *b);

void dyv_abs( const dyv *dv, dyv *r_dv);

dyv *mk_dyv_abs( const dyv *dv);

double dyv_scalar_product(const dyv *a,const dyv *b);

double dyv_dsqd(const dyv *a,const dyv *b);

double pnorm( dyv *v, double p);

double dyv_magnitude(const dyv *a);

int paul_index_in_sorted_dyv( dyv *dv, double key);
int index_in_sorted_dyv(dyv *a,double x);

void dyv_sort(const dyv *dv,dyv *r_dv);

dyv *mk_dyv_sort(const dyv *dv);

/*
 Creates a dyv of indices such that indices[i] is the origional
 location (in the unsorted dv) of the ith smallest value.
 Used when you want the location of the sorted values instead of
 the sorted vector itself.
*/
dyv *mk_sorted_dyv_indices(dyv *dv);

ivec *mk_ivec_sorted_dyv_indices(dyv *dv);


double dyv_sum(const dyv *dv);
double dyv_product(const dyv *dv);

double dyv_mean(const dyv *dv);
  /* Mean of all elements in dv */

double dyv_median(const dyv *dv); /* added by Artur Dubrawski on Aug 02 1996, efficiented by AWM */
  /* Median of all elements in dv */

double dyv_sdev(const dyv *dv);
  /* Sdev of all elements in dv */

double dyv_min(const dyv *dv);
  /* Min of all elements in dv. ERROR if dv sized zero */

double dyv_max(const dyv *dv);
  /* Max of all elements in dv. ERROR if dv sized zero */

int dyv_argmin(const dyv *dv);

int dyv_argmax(const dyv *dv);

void fprintf_dyv(FILE *s, const char *m1, const dyv *d, const char *m2);


/* Returns TRUE if any of x's elements are NaN */
/* Declared but not defined - sir 8/6/2000
bool dyv_isnan(dyv *x);*/

/* Returns TRUE if any elements are NaN or Inf */
bool dyv_is_ill_defined(dyv *x);

/* Returns 1 if all elements are 0.0 or 1.0.  Note that 0.0 and 1.0 are
   perfectly representable in IEEE 754. */
int dyv_is_binary( const dyv *dv);


/* returns TRUE if dyvs are equal */
bool dyv_equal(const dyv *dx, const dyv *dy);

void fprintf_oneline_dyv(FILE *s,const char *m1, const dyv *d, const char *m2);

/* After calling, ivec_ref(iv,k) points to the k'th smallest
   element in dv (ties handles arbitrarily) */
void indices_of_sorted_dyv(const dyv *dv,ivec *iv);

dyv *mk_sorted_dyv(dyv *x);

/* After calling, ivec_ref(result,k) points to the k'th smallest
   element in dv (ties handles arbitrarily) */
ivec *mk_indices_of_sorted_dyv(const dyv *dv);

ivec *mk_indices_of_sorted_ivec(ivec *v); /* Artur */





double dyv_partial_sum( const dyv *dv, const ivec *indices);
dyv *mk_dyv_from_ivec(const ivec *iv);
ivec *mk_ivec_from_dyv( const dyv *dv);

/* If b >= a >= 0, creates "slice" of dv on [a,b-1].  If a or b < 0, then
   it is replaced by len+a or len+b.  Indices which exceed bounds are
   replaced by appropriate indices.  Therefore this function always succeeds,
   sometimes returning an empty dyv.
*/
dyv *mk_dyv_slice( const dyv *dv, int a, int b);



#define DYV_CODE 4509

#endif /* #ifndef AMDYV_H */
