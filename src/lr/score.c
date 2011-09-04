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
   File:        score.c
   Author:      Paul Komarek
   Created:     Wed May 18 14:26:10 EDT 2005
   Description: Scoring functions for classifiers.

   Copyright 2005, The Auton Lab, CMU
*/


#include "amiv.h"
#include "amdyv.h"
#include "am_string.h"

#include "lrutils.h"
#include "score.h"



/************************************************************************/
/* ROC CURVES                                                           */
/************************************************************************/

double mk_roc_curve( ivec *outputs, dyv *predicts, ivec **x, ivec **y)
{
  /* Returns auc and stores roc curve in x and y.  Set x or y to NULL
     to prevent writing to either.  If x or y isn't NULL, then x and y
     will point to ivecs when this function returns. */

  bool write_curve;
  unsigned long int numrows, num_wrong, num_right, i, rank, row;
  unsigned long int tot_right, tot_wrong;
  double auc, unitarea;
  ivec *indexes, *order, *ocopy;
  dyv *pcopy;
  char *fstr;

  numrows = dyv_size( predicts);

  /* Create randomized copy of predicts, so that ties are broken randomly.
     Note that both predicts and outputs need to be randomized in
     parallel.  Note that we are re-seeding the random number generator
     in the same way every time, for consitency between identical learners
     and experiments. */
  order = mk_identity_ivec( numrows);

  push_current_am_srand_state();
  am_srand(8675309);
  shuffle_ivec( order);
  pop_current_am_srand_state();

  pcopy = mk_dyv_subset( predicts, order);
  ocopy = mk_ivec_subset( outputs, order);
  free_ivec( order);

  /* Allocate return vectors. */
  write_curve = ((x != NULL) && (y != NULL));
  if (write_curve) {
    *x = mk_ivec( numrows);
    *y = mk_ivec( numrows);
  }

  /* Compute unit area for a grid cell as if we gridded the unit square. */
  tot_wrong = count_in_ivec( outputs, 0);
  tot_right = numrows - tot_wrong;
  unitarea = 1.0 / (double) tot_wrong / (double) tot_right;

  /* Compute points on ROC curve, and AUC. */
  indexes = mk_indices_of_sorted_dyv( pcopy);
  free_dyv( pcopy);
  num_wrong = 0;
  num_right = 0;
  auc = 0.0;
  for (i = 0; i < numrows; i++) {
    /* Get predictions in decreasing order. */
    rank   = numrows - i - 1;
    row    = ivec_ref( indexes, rank);

    if (ivec_ref(ocopy, row)==1) {
      num_right++;
    }
    else {
      num_wrong++;
      /* Increment area by curve height = num_right. */
      auc += num_right * unitarea;
    }

    if (write_curve) {
      ivec_set( *x, i, num_wrong);
      ivec_set( *y, i, num_right);
    }
  }

  free_ivec( indexes);
  free_ivec( ocopy);

  /* Potential error message. */
  fstr = mk_copy_string( "WARNING: mk_roc_curve: %s = 0, which means I was "
                         "passed\nno %s instances.  This may be no big deal "
                         "if this is a fold.  In that\ncase, it might be "
                         "possibe to avoid the problem by decreasing the "
                         "number\n of folds. Setting AUC to 1.0.\n\n");

  /* There exists a perspective from which num_right==0 => auc=1.0.
     The same goes for num_wrong==0. */
  if (num_right == 0) {
    auc = 1.0;
    /* DO NOT use my_warning(), because it will wait for a key press. */
    fprintf( stderr, fstr, "num_right", "positive");
  }
  if (num_wrong == 0) {
    auc = 1.0;
    /* DO NOT use my_warning(), because it will wait for a key press. */
    fprintf( stderr, fstr, "num_wrong", "negative");
  }
  free_string(fstr);

  return auc;
}
