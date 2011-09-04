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


#ifndef _H_AMDYM_H
#define _H_AMDYM_H



#include "file.h"
#include "amiv.h"
#include "amdyv.h"
#include "amdyv_array.h"



typedef struct dym_struct
{
  int dym_code;
  int rows;
  int cols;
  int rows_allocated;
  double **tdarr;
} dym, *dym_ptr;



dym* mk_dym(int rows, int col);
dym *mk_copy_dym(const dym *d);
void free_dym(dym *d);
void fprintf_dym(FILE *s, const char *m1, const dym *d, const char *m2);


dym *mk_read_dym(const char* ascii_file_name);
void mk_read_dym_for_csv(const char* filename, dym **factors, dyv **outputs);
int save_dym(PFILE* s, dym* d);

int dym_rows(const dym *d);
int dym_cols(const dym *d);

void dym_scalar_mult(const dym *d, double alpha, dym *r_d);
dym *mk_dym_scalar_mult(const dym *d,double alpha);
dym *mk_copy_dym(const dym *d);
void copy_dym(dym *d, dym *r_d);

void dym_times_dyv(const dym *a, const dyv *b, dyv *result);
dyv *mk_dym_times_dyv(const dym *a, const dyv *b);

void dym_transpose( const dym *d, dym *r_d);
dym *mk_dym_transpose( const dym *a);

void dym_transpose_times_dyv(dym *p,dyv *v, dyv *r_dyv);
dyv *mk_dym_transpose_times_dyv(dym *p,dyv *v);

dym *mk_dym_subset( const dym *x, const ivec *rows, const ivec *cols);
dyv *mk_dyv_from_dym_column( const dym *dm, int col);
dyv *mk_dyv_from_dym_row( const dym *dm, int row);
ivec *mk_ivec_from_binary_dym_column( const dym *dm, int col);
dyv_array *mk_dyv_array_from_dym( const dym *dm);
dyv_array *mk_dyv_array_from_part_of_dym( const dym *dm, const ivec *rows,
                                              const ivec *atts);


#define dym_ref(d,i,j) ((d)->tdarr[i][j])
#define dym_set(d,i,j,v) (d)->tdarr[i][j] = (v)
#define dym_increment(d,i,j,v) (d)->tdarr[i][j] += (v)




#endif /* H_AMDYM_H */
