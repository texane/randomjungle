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
   File:        prec.h
   Author:      Andrew W. Moore
   Created:     Mon Feb 11 10:43:10 EST 2002
   Description: Header for Record in a sparse database

   Copyright 2002, The Auton Lab, CMU
*/


#ifndef PREC_H
#define PREC_H
#include "amiv.h"

typedef struct prec
{
  double activation;
  ivec *factors;
} prec;

prec *mk_prec(double activation,ivec *factors);

void free_prec(prec *p);

prec *mk_copy_prec(prec *old);

void fprintf_prec(FILE *s,char *m1,prec *p,char *m2);

double prec_activation(prec *p);

ivec *prec_factors(prec *p);


#endif /* #ifndef PREC_H */
