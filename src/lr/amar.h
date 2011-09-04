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
   File:        amar.h
   Author:      Andrew W. Moore
   Created:     Sat Sep 12 14:53:01 EDT 1992
   Description: Header for Obvious operations on 1-d arrays

   Copyright (C) 1992, Andrew W. Moore
*/

#ifndef AMAR_H
#define AMAR_H

#include "standard.h"
#include "ambs.h"
#include "amma.h"

void copy_ints(int *src, int *dest, int size);

void fprintf_ints(FILE *s, char *msg1, int *ints, int size, char *msg2);

void printf_ints(char *msg1, int *ints, int size, char *msg2);

void set_ints_constant(int *ints, int size, int val);

int ints_extreme(int *ints, int size, bool lowest);

int ints_min(int *ints, int size);

int ints_max(int *ints, int size);

int int_comp(const void *r1_ptr, const void *r2_ptr);

void sort_ints(int *farr, int size, int *r_farr);

void copy_longs(long *src, long *dest, int size);

void fprintf_longs(FILE *s, char *msg1, long *longs, int size, char *msg2);

void printf_longs(char *msg1, long *longs, int size, char *msg2);

void set_longs_constant(long *longs, int size, long val);

long longs_extreme(long *longs, int size, bool lowest);

long longs_min(long *longs, int size);

long longs_max(long *longs, int size);

int  long_comp(const void *r1_ptr, const void *r2_ptr);

void sort_longs(long *farr, int size, long *r_farr);

void copy_bools(bool *src, bool *dest, int size);

void fprintf_bools(FILE *s, char *msg1, bool *bools, int size, char *msg2);

void printf_bools(char *msg1, bool *bools, int size, char *msg2);

void set_bools_constant(bool *bools, int size, bool val);

bool bools_and(bool *bools, int size);

bool bools_or(bool *bools, int size);

void copy_realnums(double *src, double *dest, int size);

void fprintf_realnums(FILE *s, char *msg1, double *doubles, int size, char *msg2);

void printf_realnums(char *msg1, double *doubles, int size, char *msg2);

void set_realnums_constant(double *doubles, int size, double val);

double doubles_extreme(double *doubles, int size, bool lowest);

double doubles_min(double *doubles, int size);

double doubles_max(double *doubles, int size);

double doubles_inner_product(double *f1, double *f2, int size);

double doubles_magnitude_sqd(double *f, int size);

void doubles_scalar_multiply(double *f, int size, double alpha, double *res_f);

void doubles_add(double *f1, double *f2, int size, double *res_f);

void doubles_subtract(double *f1, double *f2, int size, double *res_f);

double doubles_scaled_dsqd(double *f1, double *f2, double *scales, int size);

double doubles_antiscaled_dsqd(double *f1, double *f2, double *scales, int size);

void doubles_reciprocal(double *farr, int size, double *res_farr);

void doubles_scalar_add(double *farr1, double *farr2, int size, double alpha, double *res_farr);

double **am_malloc_2d_realnums(int rows, int cols);
/* Declared but not defined - sir 8/6/2000
int **am_malloc_2d_ints(int rows, int cols); */

void doubles_2d_scalar_add(double **tdarr1, double **tdarr2, int rows,
 int cols, double alpha, double **res_tdarr);

void fprintf_2d_realnums(FILE *s, char *m1, double **tdarr,
 int rows, int cols, char *m2);

void free_2d_realnums(double **tdarr, int rows, int cols);

void am_free_2d_realnums(double **tdarr, int rows, int cols);
/* Declared but not defined - sir 8/6/2000
void am_free_2d_ints(int **itdarr, int rows, int cols);*/

void set_2d_realnums_constant(double **tdarr, int rows, int cols, double val);

void copy_2d_realnums(double **tdarr1, double **tdarr2, int rows, int cols);

double doubles_sum(double *farr, int size);

double doubles_mean(double *farr, int size);

double doubles_sdev(double *farr, int size);

void random_shuffle_ints(int *iarr, int size);

void sort_realnums(double *farr, int size, double *r_farr);
/*
   It is fine for farr to be the same memory as r_farr
*/

double doubles_median(double *farr, int size);

double *am_malloc_realnums(
    int size
  );

void am_free_realnums(
  double *farr,
  int size
  );

int *am_malloc_ints(
    int size
  );

void am_free_ints(
  int *iarr,
  int size
  );

long *am_malloc_longs(
    int size
  );

void am_free_longs(
  long *larr,
  int size
  );

bool *am_malloc_bools(
    int size
  );

void am_free_bools(
  bool *barr,
  int size
  );

int doubles_argmax(double *doubles, int size);			    
int doubles_argmin(double *doubles, int size);			    
int ints_argmax(int *ints, int size);			    
int ints_argmin(int *ints, int size);

void indices_sort_realnums(double *farr, int size, int *indices);

void indices_sort_integers(int *iarr, int size, int *indices); /*  Artur */

#endif /* #ifdef AMAR_H */
