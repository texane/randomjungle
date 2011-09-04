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


int dummy = -1;

#ifdef NEVER

#include "amvpv.h"

void check_vpvec_code(const vpvec* v, char* name);
void check_vpvec_access(const vpvec* v, int i, char* name);
#define VPVEC_CODE 67234

int Vpvecs_mallocked = 0;
int Vpvecs_freed = 0;

#ifdef AMFAST

#define NOTHING_TO_DO

#define check_vpvec_code(v, name) NOTHING_TO_DO

#else

void check_vpvec_code(const vpvec* v, char* name)
{
  if (v == NULL)
    {
      fprintf(stderr, "NULL vpvec passed in operation %s\n", name);
      my_error("vpvec data structure");
    }
  if (v->vpvec_code != VPVEC_CODE)
    {
      fprintf(stderr, "Attempt to access a non-allocated void* vector\n");
      fprintf(stderr, "This is in the operation %s\n", name);
      my_error("vpvec data structure error");
    }
  if (v->array_size < v->size)
    {
      my_error("check_vpvec_code: array_size and size muddled");
    }
}

#endif

#ifdef AMFAST

#define check_vpvec_access(v, i, name) NOTHING_TO_DO

#else

void check_vpvec_access(const vpvec* v, int i, char* name)
{
  check_vpvec_code(v, name);
  if (i < 0 || i >= v->size)
    {
      fprintf(stderr,"In operation \"%s\"\n",name);
      fprintf(stderr,"the vpvec (void* vector) has size = %d\n",v->size);
      fprintf(stderr,"You tried to use index i=%d\n",i);
      my_error("check_vpvec_access");
    }
}

#endif

vpvec* mk_vpvec(int size)
{
  vpvec* result = AM_MALLOC(vpvec);
  if (size < 0) my_error("mk_vpvec: size < 0 illegal");
  result->vpvec_code = VPVEC_CODE;
  result->array_size = size;
  result->size = size;
  result->vparr = AM_MALLOC_ARRAY(void*, size);
  Vpvecs_mallocked += 1;
  return result;
}

void free_vpvec(vpvec* v)
{
  check_vpvec_code(v, "free_vpvec");
  v->vpvec_code = 7777;
  AM_FREE_ARRAY(v->vparr, void*, v->array_size);
  AM_FREE(v, vpvec);
  Vpvecs_freed += 1;
}

void* safe_vpvec_ref(const vpvec* v, int i)
{
  check_vpvec_access(v, i, "vpvec_ref");
  return (v->vparr[i]);
}

void safe_vpvec_set(vpvec* v, int i, void* value)
{
  check_vpvec_access(v, i, "vpvec_set");
  v->vparr[i] = value;
}

int safe_vpvec_size(const vpvec* v)
{
  check_vpvec_code(v, "vpvec_size");
  return (v->size);
}

void add_to_vpvec(vpvec* v, void* new_val)
{
  int size;
  check_vpvec_code(v, "add_to_vpvec");
#ifndef AMFAST
  if (v->array_size < 0 || v->size > v->array_size)
    {
      my_error("vpvec size or array size is screwed.");
    }
#endif
  size = v->size;
  if (size == v->array_size)
    {
      int new_array_size = 2 * v->size + 2;
      void** vparr_new = AM_MALLOC_ARRAY(void*, new_array_size);
      int i;
      for (i = 0; i < size; i++)
	vparr_new[i] = v->vparr[i];
      AM_FREE_ARRAY(v->vparr, void*, size);
      v->vparr = vparr_new;
      v->array_size = new_array_size;
    }
  v->vparr[size] = new_val;
  v->size += 1;
}

void copy_vpvec(vpvec* v, vpvec* r_v)
{
  int i;
  check_vpvec_code(v, "copy_vpvec");
  check_vpvec_code(r_v, "copy_vpvec");
#ifndef AMFAST
  if (v->size != r_v->size)
    my_error("copy_vpvec: sizes don't match!");
#endif
  for (i = 0; i < v->size; i++)
    {
      r_v->vparr[i] = v->vparr[i];
    }
}

vpvec* mk_copy_vpvec(vpvec* v)
{
  int sz = v->size;
  vpvec* result = mk_vpvec(sz);
  copy_vpvec(v, result);
  return result;
}

#endif
