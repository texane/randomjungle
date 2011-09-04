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
   File:        prec.c
   Author:      Andrew W. Moore
   Created:     Mon Feb 11 10:43:10 EST 2002
   Description: Record in a sparse database

   Copyright 2002, The Auton Lab, CMU
*/

#include "amiv.h"
#include "am_string.h"
#include "prec.h"

prec *mk_prec(double activation,ivec *factors)
{
  prec *p = AM_MALLOC(prec);
  p -> activation = activation;
  p -> factors = mk_copy_ivec(factors);
  return p;
}

void free_prec(prec *p)
{
  free_ivec(p->factors);
  AM_FREE(p,prec);
}

prec *mk_copy_prec(prec *old)
{
  prec *p = AM_MALLOC(prec);
  p -> activation = old->activation;
  p -> factors = mk_copy_ivec(old->factors);
  return p;
}

double prec_activation(prec *p)
{
  return p->activation;
}

ivec *prec_factors(prec *p)
{
  return p->factors;
}

void fprintf_prec(FILE *s,char *m1,prec *p,char *m2)
{
  char *buff;
  fprintf(s,"%s -> activation = %g%s\n",m1,prec_activation(p),m2);
  buff = mk_printf("%s -> factors",m1);
  fprintf_ivec(s,buff,prec_factors(p),m2);
  free_string(buff);
}

