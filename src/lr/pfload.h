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
   File:        pfload.h
   Author:      Andrew W. Moore
   Created:     Mon Feb 11 10:50:18 EST 2002
   Description: Header for Load in a sparse  database

   Copyright 2002, The Auton Lab, CMU
*/


#ifndef PFLOAD_H
#define PFLOAD_H

#include "amiv.h"
#include "am_string_array.h"
#include "precs.h"

precs *mk_precs_from_filename(char *fname);

void ps_save_factors( char *fname, precs *ps);
void ps_save_activations( char *fname, precs *ps);
void ps_save( char *fname, precs *ps);

int precs_num_rows(precs *ps);

double precs_row_to_activation(precs *ps,int row);

ivec *precs_row_to_factors(precs *ps,int row);

int precs_num_factors(precs *ps);

ivec *mk_factors_from_string_array(string_array *sa,int line);

#endif /* #ifndef PFLOAD_H */
