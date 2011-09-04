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


/* This notice hereby gives fair warning that I'm going to remove this
   file at some point -- Pat */

#error If you need this, please contact pgunn@cs.cmu.edu
#ifdef NEVER

/*

amiv analogue for void*'s.  
Author: Scott Davies, pattern-matched off of Andrew's amiv code.

Insert obligatory "junk like this is why we should be using C++ 
with STL instead" whinge here.

*/

#ifndef AMVPV_H
#define AMVPV_H

#include "amma.h"

typedef struct vpvec_struct
{
  int vpvec_code;
  int array_size;
  int size;
  void** vparr;
} vpvec;

void* safe_vpvec_ref(const vpvec* v, int i);
void safe_vpvec_set(vpvec* v, int i, void* value);
int safe_vpvec_size(const vpvec* v);

#ifdef AMFAST

#define vpvec_ref(v, i) ((v)->vparr[i])
#define vpvec_set(v,i,val) ((v)->vparr[i] = (val))
#define vpvec_size(v) ((v)->size)

#else

#define vpvec_ref(v, i) (safe_vpvec_ref(v, i))
#define vpvec_set(v, i, val) (safe_vpvec_set(v, i, val))
#define vpvec_size(v) (safe_vpvec_size(v))

#endif

vpvec* mk_vpvec(int size);
void free_vpvec(vpvec* v);
void add_to_vpvec(vpvec* v, void* val);

void copy_vpvec(vpvec* v, vpvec* r_v);
/* Note: *not* a deep copy! */

vpvec* mk_copy_vpvec(vpvec* v);

#endif

#endif
