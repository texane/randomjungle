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
   File:        amdyv_array.h
   Author:      Andrew W. Moore
   Created:     June 15, 2004
   Description: Extensions and advanced amdm stuff

   Copyright 1996, Schenley Park Research

   This file contains utility functions involving dyv_arrays, dyvs and ivecs.

   The prototypes of these functions used to be declared in amdmex.h.
   Now it's amdyv_array.h
*/

#ifndef AMDYV_ARRAY_H
#define AMDYV_ARRAY_H

#include "standard.h"
#include "ambs.h"
#include "amiv.h"
#include "amdyv.h"
#include "am_string_array.h"

typedef struct dyv_array
{
  int size;
  int array_size;
  dyv **array;
} dyv_array;



dyv_array *mk_empty_dyv_array(void);
void add_to_dyv_array(dyv_array *da, const dyv *dv);

dyv_array *mk_const_dyv_array(dyv *base_vec, int size);

dyv_array *mk_rectangular_dyv_array( int numdyvs, int dyvlen);

/* Creates an dyv array with size entries each composed of an dyv of size 0  [e.g. mk_dyv(0)] */
dyv_array *mk_zero_dyv_array(int size);

int dyv_array_size(const dyv_array *da);
dyv *safe_dyv_array_ref(const dyv_array *da,int idx);

dyv_array *mk_array_of_zero_length_dyvs(int size);
dyv_array *mk_array_of_null_dyvs(int size);

dyv_array *mk_dyv_array( int size);
dyv_array *mk_dyv_array_subset( dyv_array *da, ivec *indices);

#define dyv_array_ref(dva,i) ((dva)->array[i])

void dyv_array_set(dyv_array *dva, int idx, const dyv *dv);
void dyv_array_set_no_copy( dyv_array *dva, int idx, dyv *dv);
void fprintf_dyv_array(FILE *s, const char *m1, const dyv_array *da, const char *m2);

void free_dyv_array(dyv_array *da);
dyv_array *mk_copy_dyv_array(const dyv_array *da);
void dyv_array_remove(dyv_array *dva,int idx);



#endif
