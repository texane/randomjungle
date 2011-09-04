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
   File:        lr.c
   Author:      Paul Komarek
   Created:     Thu Mar  6 12:14:39 EST 2003
   Description: Logistic regression implmentation.

   Copyright 2003, The Auton Lab, CMU
*/

#include <stdio.h>
#include <float.h>
#include <math.h>


#include "amma.h"
#include "amdyv.h"
#include "amdym.h"
#include "amdyv_array.h"
#include "spardat.h"

#include "lrutils.h"
#include "lin_conjgrad.h"
#include "lr.h"



/***********************************************************************/
/* Option Parsing                                                      */
/***********************************************************************/

lr_options *mk_lr_options(void)
{
  lr_options *opts;
  opts = AM_MALLOC( lr_options);

  /* Options which are always available. */
  opts->rrlambda   = 10.0;

  /*   Termination criteria for lr iterations. */
  opts->lreps      = 0.05;
  opts->lrmax      = 30;

  /* cg options are available in conjuagate gradient runs. */
  opts->cgbinit    = 1;
  opts->cgdeveps   = 0.005;    /* suggestion: 0.005 */
  opts->cgeps      = 0.000;  /* multiplied by initial CG rsqr. */
  opts->cgmax      = 200;

  opts->cgwindow   = 3;      /* Number of bad iterations allowed. */
  opts->cgdecay    = 1000.0; /* Factor worse than best-seen that is allowed. */
  
  return opts;
}

lr_options *mk_copy_lr_options( lr_options *opts)
{
  lr_options *newopts;
  newopts = AM_MALLOC( lr_options);

  newopts->rrlambda        = opts->rrlambda;
  newopts->lreps           = opts->lreps;
  newopts->lrmax           = opts->lrmax;
  newopts->cgbinit         = opts->cgbinit;
  newopts->cgdeveps        = opts->cgdeveps;
  newopts->cgeps           = opts->cgeps;
  newopts->cgmax           = opts->cgmax;
  newopts->cgwindow        = opts->cgwindow;
  newopts->cgdecay         = opts->cgdecay;

  return newopts;
}

void free_lr_options( lr_options *opts)
{
  if (opts != NULL) AM_FREE( opts, lr_options);
  return;
}

void parse_lr_options( lr_options *opts, int argc, char **argv)
{
  /* Create an options struct with mk_options(), then call this function
     to get the options from the command line. */

  opts->rrlambda = double_from_args( "rrlambda", argc, argv, opts->rrlambda);

  opts->lreps = double_from_args( "lreps", argc, argv, opts->lreps);
  opts->lrmax = int_from_args( "lrmax", argc, argv, opts->lrmax);

  /* check_lr_options() ensures not both cgdeveps and cgeps are specified
     on the command line. */
  opts->cgdeveps = double_from_args( "cgdeveps", argc, argv, opts->cgdeveps);
  if (opts->cgdeveps != 0.0) {
    opts->cgeps = 0.0;
    opts->cgbinit = TRUE;
  }
  opts->cgeps = double_from_args( "cgeps", argc, argv, opts->cgeps);
  if (opts->cgeps != 0.0) opts->cgdeveps = 0.0;

  opts->cgmax = int_from_args( "cgmax", argc, argv, opts->cgmax);

  opts->cgwindow = int_from_args( "cgwindow", argc, argv, opts->cgwindow);
  opts->cgdecay  = double_from_args( "cgdecay", argc, argv, opts->cgdecay);

  return;
}

void check_lr_options( lr_options *opts, int argc, char **argv)
{
  /* Checks that options are sane.  Aborts if options are not sane. */
  int use_cgdeveps, use_cgeps;
  char *co = "check_lr_options";

  if (opts == NULL) my_errorf( "%s: opts==NULL", co);

  /* Check that the Ridge Regression parameter is non-negative. */
  if (opts->rrlambda < 0.0) {
    my_errorf( "%s: rrlambda(=%g) must be non-negative.", co, opts->rrlambda);
  }

  /* Check lreps and lrmax. */
  if (opts->lreps < 1e-10) my_errorf( "%s: lreps(=%g) < 1e-10 is unreasonable",
                                     co, opts->lreps);
  if (opts->lrmax < 0) my_errorf( "%s: lrmax(=%d) < 0 is unreasonable",
                                     co, opts->lrmax);

  use_cgdeveps = (index_of_arg( "cgdeveps", argc, argv) > 0);
  use_cgeps = (index_of_arg( "cgeps", argc, argv) > 0);

  /* Check cgdeveps, cgeps and cgmax. */
  if ( use_cgeps && use_cgdeveps) {
    my_errorf( "%s: Cannot specify both cgdeveps and cgeps", co);
  }
  if (use_cgdeveps && opts->cgdeveps < 1e-10) {
    my_errorf( "%s: cgdeveps(=%g) < 1e-10 is unreasonable",
               co, opts->cgdeveps);
  }
  if (use_cgeps && opts->cgeps < 1e-10) {
    my_errorf( "%s: cgeps(=%g) < 1e-10 is unreasonable", co, opts->cgeps);
  }
  if (opts->cgmax < 0) my_errorf( "%s: cgmax(=%d) < 0 is unreasonable",
                                     co, opts->cgmax);
  /* Check cgwindow and cgdecay. */
  if (opts->cgwindow < 0) my_errorf( "%s: cgwindow(=%d) < 0 is unreasonable",
                                     co, opts->cgwindow);
  if (opts->cgdecay < 1.0) my_errorf( "%s: cgdecay(=%g) < 1.0 is unreasonable",
                                     co, opts->cgdecay);

  return;
}



/***********************************************************************/
/* LR_STATE STRUCT                                                     */
/***********************************************************************/

lr_state *mk_lr_state( lr_train *lrt, lr_options *opts)
{
  lr_state *lrs;
  lrs = AM_MALLOC( lr_state);

  lrs->b0         = 0.0;

  lrs->b          = mk_zero_dyv(lrt->numatts-1);
  lrs->n          = mk_zero_dyv(lrt->numrows);
  lrs->u          = mk_zero_dyv(lrt->numrows);
  lrs->w          = mk_zero_dyv(lrt->numrows);
  lrs->z          = mk_zero_dyv(lrt->numrows);

	lrs->converged  = 0;

  return lrs;
}

lr_state *mk_copy_lr_state( lr_state *lrs)
{
  lr_state *lrscopy;

  lrscopy = AM_MALLOC( lr_state);

  lrscopy->b0          = lrs->b0;
  lrscopy->b           = mk_copy_dyv( lrs->b);
  lrscopy->n           = mk_copy_dyv( lrs->n);
  lrscopy->u           = mk_copy_dyv( lrs->u);
  lrscopy->w           = mk_copy_dyv( lrs->w);
  lrscopy->z           = mk_copy_dyv( lrs->z);
  lrscopy->converged   = lrs->converged;


  return lrscopy;
}

void fprintf_lr_state( FILE *f, char *pre, lr_state *lrs)
{
  fprintf( f, "%sb0: %g\n", pre, lrs->b0);
  fprintf( f, "%s", pre);
  fprintf_oneline_dyv( f, "b:", lrs->b, "\n");
  fprintf( f, "%s", pre);
  fprintf_oneline_dyv( f, "n:", lrs->n, "\n");
  fprintf( f, "%s", pre);
  fprintf_oneline_dyv( f, "u:", lrs->u, "\n");
  fprintf( f, "%s", pre);
  fprintf_oneline_dyv( f, "w:", lrs->w, "\n");
  fprintf( f, "%s", pre);
  fprintf_oneline_dyv( f, "z:", lrs->z, "\n");
  fprintf( f, "%s", pre);
  fprintf( f, "%sconverged: %i\n", pre, lrs->converged);

  return;
}

/* Copies initb to lrs->b. */
void lr_state_overwrite_b( lr_state *lrs, dyv *initb)
{
  copy_dyv( initb, lrs->b);
  return;
}

void free_lr_state( lr_state *lrs)
{
  if (lrs != NULL) {
    if (lrs->b != NULL) free_dyv( lrs->b);
    if (lrs->n != NULL) free_dyv( lrs->n);
    if (lrs->u != NULL) free_dyv( lrs->u);
    if (lrs->w != NULL) free_dyv( lrs->w);
    if (lrs->z != NULL) free_dyv( lrs->z);
    AM_FREE( lrs, lr_state);
  }
  return;
}

/***********************************************************************/
/* LR_STATEARR STRUCT                                                  */
/***********************************************************************/

lr_statearr *mk_array_of_null_lr_states( int size)
{
  int i;
  lr_statearr *lrsarr;
  lrsarr = AM_MALLOC( lr_statearr);
  lrsarr->size = size;
  lrsarr->arr = AM_MALLOC_ARRAY( lr_state *, size);
  for (i=0; i<size; ++i) lrsarr->arr[i] = NULL;
  return lrsarr;
}

lr_state *lr_statearr_ref( lr_statearr *lrsarr, int index)
{
#ifndef AMFAST
  if (index < 0 || index > lrsarr->size) {
    my_errorf( "lr_statearr_ref: illegal index %d not within [%d,%d]",
               index, 0, lrsarr->size-1);
  }
#endif
  return lrsarr->arr[index];
}

/* Copies lr_state. */
void lr_statearr_set( lr_statearr *lrsarr, int index, lr_state *lrs)
{
#ifndef AMFAST
  if (index < 0 || index > lrsarr->size) {
    my_errorf( "lr_statearr_set: illegal index %d not within [%d,%d]",
               index, 0, lrsarr->size-1);
  }
#endif
  lrsarr->arr[index] = mk_copy_lr_state( lrs);
  return;
}

void free_lr_statearr( lr_statearr *lrsarr)
{
  int i;
  if (lrsarr != NULL) {
    if (lrsarr->arr != NULL) {
      for (i=0; i < lrsarr->size; ++i) {
        if (lrsarr->arr[i] != NULL) free_lr_state( lrsarr->arr[i]);
      }
      AM_FREE_ARRAY( lrsarr->arr, lr_state *, lrsarr->size);
    }
    AM_FREE( lrsarr, lr_statearr);
  }
}

/***********************************************************************/
/* LR_TRAIN STRUCT                                                     */
/***********************************************************************/

lr_train *mk_lr_train_from_dym( dym *factors, dyv *outputs, lr_options *opts)
{
  /* Set rows to NULL if you want all rows from ds to be used. */

  int numrows, numatts;
  lr_train *lrt;

  numrows = dym_rows( factors);
  numatts = dym_cols( factors)+1; /* Number of factors including constant. */

  /* Create lr lrt structure. */
  lrt = AM_MALLOC(lr_train);

  /* Copy in opts. */
  lrt->opts = mk_copy_lr_options( opts);

  /* Assign factors and outputs into lr structure. */
  lrt->X = NULL;
  lrt->M = factors;

  /* Outputs. */
  lrt->y = mk_copy_dyv( outputs);
  if (!dyv_is_binary( outputs)) {
    my_error( "mk_lr_train: Error: outputs are not binary.\n");
  }

  /* Set log likelihood of saturated model. */
  lrt->likesat = 0.0;

  /* Initialize remainder of lr struct */
  lrt->numatts = numatts;
  lrt->numrows = numrows;

  /* Create lr_state member. */
  lrt->lrs = mk_lr_state( lrt, opts);

  /* Now that the structure is complete, update n and u to prepare for
     iterations. */
  lr_train_update_n(lrt);
  lr_train_update_u(lrt);

  return lrt;
}

lr_train *mk_lr_train_from_spardat( spardat *X, lr_options *opts)
{
  /* Copies spardat X. */
  int numrows, numatts;
  lr_train *lrt;

  numrows = spardat_num_rows(X);
  numatts = spardat_num_atts(X)+1; /* Add 1 for constant att. */
  lrt = AM_MALLOC(lr_train);

  /* Copy in opts. */
  lrt->opts = mk_copy_lr_options( opts);

  /* Do not make a copy of the caller's spardat; too expensive. */
  lrt->X = X;
  lrt->M = NULL;

  /* Initialize reaminder of lr struct */
  lrt->numatts    = numatts;
  lrt->numrows    = numrows;

  /* Futz with 0-1 probabilities. */
  lrt->y = mk_dyv_from_ivec( X->row_to_outval);

  /* Set log likelihood of saturated model. */
  lrt->likesat = 0.0;

  /* Make lr_state member. */
  lrt->lrs = mk_lr_state( lrt, opts);

  /* Now that the structure is complete, update n and u to prepare for
     iterations. */
  lr_train_update_n(lrt);
  lr_train_update_u(lrt);

  return lrt;
}

lr_train *mk_copy_lr_train( const lr_train *source)
{
  lr_train *dest;

  dest = AM_MALLOC(lr_train);
  dest->opts = mk_copy_lr_options( source->opts);

  /* Don't copy spardat or factors dym.  Just keep pointer. */
  dest->X = source->X;
  dest->M = source->M;

  dest->numatts    = source->numatts;
  dest->numrows    = source->numrows;
  dest->y          = mk_copy_dyv(source->y);
  dest->likesat    = source->likesat;
  dest->lrs        = mk_copy_lr_state( source->lrs);
  return dest;
}

void free_lr_train( lr_train *lrt)
{
  if (lrt != NULL) {
    if (lrt->lrs != NULL)      free_lr_state( lrt->lrs);
    if (lrt->y != NULL)        free_dyv( lrt->y);
    if (lrt->opts != NULL)     free_lr_options( lrt->opts);
    AM_FREE(lrt, lr_train);
  }
  return;
}

void fprintf_lr_train( FILE *f, char *pre, lr_train *lrt)
{
  int numatts, numrows;

  numatts = lrt->numatts;
  numrows = lrt->numrows;

  fprintf( f, "%snumatts: %d\n", pre, lrt->numatts);
  fprintf( f, "%snumrows: %d\n", pre, lrt->numrows);

  /* Print lr_state member. */
  if (numatts < 15 && numrows < 500) fprintf_lr_state( f, pre, lrt->lrs);

  return;
}

/* Copies initb to lr->b. */
void lr_train_overwrite_b( lr_train *lrt, dyv *initb)
{
  lr_state_overwrite_b( lrt->lrs, initb);
  return;
}



/***********************************************************************/
/* LR ITERATIONS                                                       */
/***********************************************************************/

void lr_train_update_n( lr_train *lrt)
{
  /* X,b -> n, or M,b -> n */
  /* assumes that X is binary, but M can be any matrix. */
  if ( lrt->X != NULL) lr_compute_n_from_spardat( lrt->X, lrt_b0_ref(lrt),
                                                  lrt_b_ref(lrt),
                                                  lrt_n_ref(lrt));
  else lr_compute_n_from_dym( lrt->M,  lrt_b0_ref(lrt), lrt_b_ref(lrt),
                              lrt_n_ref(lrt));
  return;
}

void lr_train_update_u( lr_train *lrt)
{
  /* n -> u */
  lr_compute_u_from_n( lrt_n_ref(lrt), lrt_u_ref(lrt));
  return;
}

void lr_train_update_w( lr_train *lrt)
{
  /* u -> w */
  int i;
  double ui, val;
  for (i=0; i < lrt->numrows; ++i) {
    ui  = dyv_ref( lrt_u_ref(lrt), i);
    val = ui * (1-ui);
    dyv_set( lrt_w_ref(lrt), i, val);
  }
  return;
}

void lr_train_update_z( lr_train *lrt)
{
  /* y,n,u,w -> z */
  int i;
  double yi, ni, ui, wi, val;

  for (i=0; i < lrt->numrows; ++i) {
    yi  = dyv_ref( lrt->y, i);
    ni  = dyv_ref( lrt_n_ref(lrt), i);
    ui  = dyv_ref( lrt_u_ref(lrt), i);
    wi  = dyv_ref( lrt_w_ref(lrt), i);
    val = ni + (yi-ui) / wi;
    
#ifndef AMFAST
    if (!am_isnum( val)) {
      my_errorf( "lr_train_update_z: NaN or Inf problem: val is %f.\n"
		 "Inputs: i=%d, yi=%f, ni=%f, ui=%f, wi=%f\n",
		 val, i, yi, ni, ui, wi);
    }
#endif

    dyv_set( lrt_z_ref(lrt), i, val);
  }
  return;
}

int lr_train_update_b( lr_train *lrt)
{
  /* X,w,z -> b */
  /*
                   [1t]                [1t]
    Compute b = (( [--] W [1|X])^-1) * [--] W z, where W = diag(w).
                   [Xt]                [Xt]
  */
  int numatts, i, iters;
  double cgeps, cgdeveps, val;
  dyv *B, *initb;

  numatts = lrt->numatts;

  /* We are now using initial CG residuaal for scaling cgeps.
     This is best done inside mk_lr_cgresult(). */
  /* cgeps = lrt->numatts * lrt->opts->cgeps; */
  cgeps = lrt->opts->cgeps;
  cgdeveps = lrt->opts->cgdeveps;

  /* Create initb. */
  initb = NULL;
  if (lrt->opts->cgbinit) {
    initb = mk_dyv( numatts);
    dyv_set( initb, 0, lrt_b0_ref(lrt));
    for (i=1; i<numatts; ++i) {
      val = dyv_ref( lrt_b_ref(lrt), i-1);
      dyv_set( initb, i, val);
    }
  }

  B = mk_lr_update_b_conjugate_gradient_helper( lrt, cgeps, cgdeveps,
                                                lrt->opts->cgmax, &iters,
                                                initb);

  if (initb != NULL) free_dyv( initb);

  /* Break newb into ( b0, b ). */
  lrt_b0_set(lrt, dyv_ref( B, 0));
  for (i=1; i<numatts; ++i) {
    val = dyv_ref( B, i);
    dyv_set( lrt_b_ref(lrt), i-1, val);
  }

  free_dyv( B);

  /* Hitting cgmax is considered a failure. */
  if ( iters > lrt->opts->cgmax) return -2;

  return 1;
}

int lr_train_iterate( lr_train *lrt)
{
  int rval;

  lr_train_update_w(lrt);
  lr_train_update_z(lrt);
  rval = lr_train_update_b(lrt);

  if (rval != 0) {
    /* So long as we don't have a fatal error, update n and u. */
    lr_train_update_n(lrt);
    lr_train_update_u(lrt);
  }

  return rval;
}

/***********************************************************************/
/* METRICS                                                             */
/***********************************************************************/

double lr_log_likelihood_basic( dyv *y, dyv *u)
{
  /* Compute log likelihood L(b) = Sum( yi*ln(ui) + (1-yi)*ln(1-ui) ). */
  /* Note that this falls apart if u == 0.0 or u == 1.0, which it should
     never be for the logit. */
  int numrows, row;
  double sum, val, ui, yi;

  numrows = dyv_size( y);

  /* Compute log likelihood. */
  sum = 0.0;
  for (row=0; row<numrows; ++row) {
    ui  = dyv_ref( u, row);
    yi  = dyv_ref( y, row);
    val = yi*log(ui) + (1.0 - yi)*log(1.0 - ui);
    sum += val;
  }

  /* Done. */
  return sum;
}

double lr_log_likelihood_from_deviance( double deviance, double likesat)
{
  return likesat - (0.5 * deviance);
}

double lr_deviance_from_log_likelihood( double likelihood, double likesat)
{
  return -2.0 * ( likelihood - likesat);
}

double lr_deviance_basic( dyv *y, dyv *u)
{
  int i;
  double yi, ui, sum, deviance;

  /* deviance = -2*sum(yi*ln(yi/ui) + (1-yi)ln((1-yi)/(1-ui))),
     but with binary yi we have to compute the terms conditionally. */
  sum = 0.0;
  for (i=0; i < dyv_size(y); ++i) {
    yi = dyv_ref( y, i);
    ui = dyv_ref( u, i);

    if (yi != 0) sum += log(ui);
    else sum += log(1-ui);

    /* Stop summing after overflow. */
    if (sum < -FLT_MAX) return FLT_MAX;
  }

  deviance = -2*sum;
  return deviance;
}

double lr_deviance_from_spardat_b( const spardat *X, dyv *y, double b0,
                                   dyv *b)
{
  int numrows;
  double dev;
  dyv *n, *u;

  numrows = spardat_num_rows( X);
  n = mk_dyv( numrows);
  u = mk_dyv( numrows);

  lr_compute_n_from_spardat( X, b0, b, n);
  lr_compute_u_from_n( n, u);
  dev = lr_deviance_basic( y, u);

  free_dyv( n);
  free_dyv( u);
  return dev;
}

double lr_deviance_from_dym_b( const dym *M, dyv *y, double b0, dyv *b)
{
  int numrows;
  double dev;
  dyv *n, *u;

  numrows = dym_rows( M);
  n = mk_dyv( numrows);
  u = mk_dyv( numrows);

  lr_compute_n_from_dym( M, b0, b, n);
  lr_compute_u_from_n( n, u);
  dev = lr_deviance_basic( y, u);

  free_dyv( n);
  free_dyv( u);
  return dev;
}

double lr_train_deviance( lr_train *lrt)
{
  double likelihood, dev;
  likelihood = lr_log_likelihood_basic( lrt->y, lrt_u_ref(lrt));
  dev = lr_deviance_from_log_likelihood( likelihood, lrt->likesat);
  return dev;
}

/***********************************************************************/
/* LR_PREDICT                                                          */
/***********************************************************************/

lr_predict *mk_lr_predict( double b0, dyv *b)
{
  lr_predict *lrp;
  lrp = AM_MALLOC( lr_predict);
  lrp->b0       = b0;
  lrp->b        = mk_copy_dyv( b);
  return lrp;
}

lr_predict *mk_copy_lr_predict( lr_predict *lrp)
{
  lr_predict *lrpcopy;
  lrpcopy = mk_lr_predict( lrp->b0, lrp->b);
  return lrpcopy;
}

void free_lr_predict( lr_predict *lrp)
{
  if (lrp != NULL) {
    if (lrp->b   != NULL) free_dyv( lrp->b);
    AM_FREE( lrp, lr_predict);
  }
  return;
}

double lr_predict_predict( ivec *posatts, dyv *attvals, lr_predict *lrp)
{
  return lr_prediction( lrp->b0, lrp->b, posatts, attvals);
}

/***********************************************************************/
/* UTILITY                                                             */
/***********************************************************************/

double lr_prediction( double b0, dyv *b, ivec *posatts, dyv *attvals)
{
  /* posatts should be a sivec. */
  /* e^n / (1+e^n), n = b0 + lrb->b . posatts. */
  double n, en, result;

  /* Compute n = <lrb->b, posatts>. */
  n = b0;
  if (posatts != NULL) n += dyv_partial_sum( b, posatts);
  else n += dyv_scalar_product( b, attvals);

  /* Compute model value, u = e^n / (1.0 + e^n) */
  en = exp(n);
  result = en / (1.0 + en);

  return result;
}

void lr_compute_n_from_spardat( const spardat *X, double b0, dyv *b, dyv *n)
{
  int numrows, row;
  ivec *posatts;
  double sum;

  numrows = spardat_num_rows( X);
  for (row=0; row < numrows; ++row) {
    posatts = spardat_row_to_posatts( X, row);
    /* Remember that at one time we made a copy of posatts because
       of a very mysterious and hardwared/compiler-looking bug. */
    sum = dyv_partial_sum( b, posatts);
    sum += b0;
    dyv_set( n, row, sum);
  }
  return;
}

void lr_compute_n_from_dym( const dym *M, double b0, dyv *b, dyv *n)
{
  int numrows, numgood, row, j;
  double sum;

  numrows = dym_rows( M);
  numgood = dyv_size(b);
  for (row=0; row < numrows; ++row) {
    sum = 0.0;
    for (j=0; j<numgood; ++j) sum += dym_ref( M, row, j) * dyv_ref( b, j);
    sum += b0;
    dyv_set( n, row, sum);
  }
  return;
}

void lr_compute_u_from_n( dyv *n, dyv *u)
{
  int numrows, i;
  double en, val, ni;
  numrows = dyv_size( n);

  for (i=0; i < numrows; ++i) {
    ni  = dyv_ref( n, i);
    en = exp(ni);
    val = en / (1.0 + en);
    dyv_set( u, i, val);
  }

  return;
}

dyv *mk_lr_XtWXv_dyv( const lr_train *lrt, const dyv *v)
{
  /* Compute [1t]
             [--] W [1|X] v
             [Xt]            */

  double v0, cterm;
  dyv *subv, *Xv, *XtWXv;

  /* Split v into v0=v[0] and subv=v[1:] */
  v0 = dyv_ref( v, 0);
  subv = mk_dyv_slice( v, 1, dyv_size( v));

  /* Compute [1|X] v. */
  if (lrt->X != NULL) Xv = mk_spardat_times_dyv( lrt->X, subv);
  else Xv = mk_dym_times_dyv( lrt->M, subv);
  dyv_scalar_add( Xv, v0, Xv);
  free_dyv( subv);

  /* Compute W [1|X] v. */
  dyv_mult( Xv, lrt_w_ref(lrt), Xv);    /* Xv now stores WXv. */

  /* Compute Xt W [1|X] v and  1t W [1|X] v separately. Both get stored in
     XtWXv. */
  if (lrt->X != NULL) XtWXv = mk_spardat_transpose_times_dyv( lrt->X, Xv);
  else XtWXv = mk_dym_transpose_times_dyv( lrt->M, Xv);
  cterm = dyv_sum( Xv);
  dyv_insert( XtWXv, 0, cterm);

  free_dyv( Xv);

  return XtWXv;
}

dyv *mk_lr_XtWz_dyv( const lr_train *lrt)
{
  /* XtWz = Xt v = [ r_i ], v = w * z, elementwise,
     r_i = dyv_partial_sum( v, posrows_i) */
  double cterm;
  dyv *Wz, *XtWz;

  /*
            [1t]    
    Compute [--] W z
            [Xt]
  */

  /* Compute Wz. */
  Wz = mk_dyv_mult( lrt_w_ref(lrt), lrt_z_ref(lrt));

  /* Compute XtWz. */
  if (lrt->X != NULL) XtWz = mk_spardat_transpose_times_dyv( lrt->X, Wz);
  else XtWz = mk_dym_transpose_times_dyv( lrt->M, Wz);

  /* Insert 1t Wz at beginning of XtWz. */
  cterm = dyv_sum( Wz);
  dyv_insert( XtWz, 0, cterm);
  free_dyv( Wz);

  return XtWz;
}

/***********************************************************************/
/* CONJGRAD HELPERS                                                    */
/***********************************************************************/

void *lr_cg_mk_copy_userdata( const void *userdata)
{
  return (void *) mk_copy_lr_train( (lr_train *) userdata);
}

void lr_cg_free_userdata( void *userdata)
{
  free_lr_train( (lr_train *) userdata);
  return;
}

void lr_cg_multA( const dyv *v, dyv *result, void *userdata)
{
  double lambda;
  lr_train *lrt;
  dyv *Av, *lv;

  lrt = (lr_train *) userdata;

  /* Do sparse matrix-vector multiply. */
  Av = mk_lr_XtWXv_dyv( lrt, v);

  lambda = lrt->opts->rrlambda;
  if (lambda > 0.0) {
    /* Add Ridge Regression term. */
    lv = mk_dyv_scalar_mult( v, lambda);
    dyv_set( lv, 0, 0.0);  /* Don't penalize constant term. */
    dyv_plus( Av, lv, result);
    free_dyv( lv);
  }
  else {
    /* Don't do Ridge Regression. */
    copy_dyv( Av, result);
  }

  free_dyv( Av);
  return;
}

double lr_deviance_from_cg( lr_train *lrt, conjgrad *cg)
{
  int numrows;
  double cgb0, likelihood, dev;
  dyv *cgb, *cgn, *cgu;

  /* Get beta. */
  cgb = mk_copy_dyv( conjgrad_x_ref( cg)); /* good params */
  cgb0 = dyv_ref( cgb, 0);
  dyv_remove( cgb, 0);

  numrows = lrt->numrows;
  cgn = mk_dyv( numrows);
  cgu = mk_dyv( numrows);

  /* Compute u and n. */
  if (lrt->X != NULL) lr_compute_n_from_spardat( lrt->X, cgb0, cgb, cgn);
  else lr_compute_n_from_dym( lrt->M, cgb0, cgb, cgn);
  free_dyv( cgb);
  lr_compute_u_from_n( cgn, cgu);
  free_dyv( cgn);

  /* Compute likelihood and deviance. */
  likelihood = lr_log_likelihood_basic( lrt->y, cgu);
  free_dyv( cgu);

  dev = lr_deviance_from_log_likelihood( likelihood, lrt->likesat);

  return dev;
}

dyv *mk_lr_cgresult_cgeps( lr_train *lrt, double unscaled_cgeps,
                           int maxiters, conjgrad *cg)
{
  int iters, bestiter, window;
  double rsqr, bestrsqr, decay, cgeps, decthresh;
  dyv *x, *result, *rsqrhist;
  dyv_array *paramhist;

  /* Initialize paramters. */
  rsqrhist    = mk_constant_dyv( maxiters, FLT_MAX);
  paramhist  = mk_array_of_null_dyvs( maxiters);
  bestrsqr   = FLT_MAX;
  window     = lrt->opts->cgwindow;
  decay      = lrt->opts->cgdecay;
  decthresh  = FLT_MAX;

  /* Scale cgeps. */
  rsqr = sqrt(dyv_scalar_product( cg->cgs->r, cg->cgs->r));
  if (Verbosity >= 2) printf( "    CGINITIAL RSQR: %g\n", rsqr);
  cgeps = unscaled_cgeps * rsqr;

  /* Store initial position in history. */
  iters = 0;
  dyv_set( rsqrhist, iters, rsqr);
  dyv_array_set( paramhist, iters, conjgrad_x_ref( cg));
  iters += 1;

  /* Abort iterations if rsqr gets too small for calcs to proceed. */
  while (rsqr >= cgeps) {
    if (Verbosity > 3) {
      fprintf_oneline_dyv( stdout, "    CG POS:", cg->cgs->x, "\n");
    }

    /* Non-epsilon termination conditions. */
    if (iters >= maxiters) break;
    if (window <= 0) break;
    if (rsqr > decthresh) break;

    /* Iterate. */
    cgiter( cg);

    /* CG resisdual Euclidean norm. */
    rsqr = dyv_magnitude( conjgrad_r_ref( cg));

    /* Store history. */
    dyv_set( rsqrhist, iters, rsqr);
    dyv_array_set( paramhist, iters, conjgrad_x_ref( cg));
    if (Verbosity >= 2) printf( "    CGEPS RSQR: %g\n", rsqr);


    /* Update records. */
    if (rsqr <= bestrsqr) {
      bestrsqr  = rsqr;
      window    = lrt->opts->cgwindow;
      decthresh = decay * bestrsqr;
    }
    else window -= 1;

    /* Count number of iters. */
    iters += 1;
  }

  /* Select parameters. */
  /* CG residual: use last iteration's parameter vector. */
  /* x = conjgrad_x_ref( cg); */
  /* Get best params from paramhist. */
  bestiter = dyv_argmin( rsqrhist);
  x = dyv_array_ref( paramhist, bestiter);
  if (x == NULL) {
    my_errorf( "mk_lr_cgresult_cgeps: NULL param vec %d", bestiter);
  }

  if (Verbosity >= 2) {
    rsqr = sqrt(dyv_scalar_product( cg->cgs->r, cg->cgs->r));
    printf( "    CGFINAL RSQR: %g\n", rsqr);
  }

  result = mk_copy_dyv( x);
  free_dyv_array( paramhist);
  free_dyv( rsqrhist);
  return result;

}

dyv *mk_lr_cgresult_cgdeveps( lr_train *lrt, double cgdeveps,
                              int maxiters, conjgrad *cg)
{
  int iters, bestiter, window;
  double dev, olddev, bestdev, rsqr, decay, decthresh;
  dyv *devhist, *x, *result;
  dyv_array *paramhist;

  /* Run conjugate gradient. */
  devhist    = mk_constant_dyv( maxiters, FLT_MAX);
  paramhist  = mk_array_of_null_dyvs( maxiters);
  dev        = -FLT_MAX;
  bestdev    = FLT_MAX;
  window     = lrt->opts->cgwindow;
  decay      = lrt->opts->cgdecay;
  decthresh  = FLT_MAX;

  /* Scale cgeps. */
  rsqr = sqrt(dyv_scalar_product( cg->cgs->r, cg->cgs->r));

  /* Store initial position in history. */
  iters = 0;
  dev = lr_deviance_from_cg( lrt, cg);
  dyv_set( devhist, iters, dev);
  dyv_array_set( paramhist, iters, conjgrad_x_ref( cg));
  iters += 1;

  /* Abort the iters if rsqr gets too small for calcs to proceed. */
  while (rsqr > 1e-300) {
    if (Verbosity > 3) {
      fprintf_oneline_dyv( stdout, "    CG POS:", cg->cgs->x, "\n");
    }

    /* Non-deviance termination criteria. */
    if (iters > maxiters) break; /* Strict, since we start with iters=1. */
    if (window <= 0) break;
    if (dev > decthresh) break;

    /* Iterate. */
    olddev  = dev;
    cgiter( cg);

    /* Relative difference of deviance. */
    dev = lr_deviance_from_cg( lrt, cg);
    if (dev <= bestdev) {
      bestdev   = dev;
      window    = lrt->opts->cgwindow;
      decthresh = decay * bestdev;
    }
    else window -= 1;

    /* Store history. */
    dyv_set( devhist, iters, dev);
    dyv_array_set( paramhist, iters,
                   conjgrad_x_ref( cg) /* good params */);
    if (Verbosity >= 2) printf( "CG DEVIANCE: %g\n", dev);

    /* Terminate on rel diff of deviance. */
    if (fabs(olddev-dev) < dev*cgdeveps) break;

    /* Count number of iters. */
    iters += 1;

    /* We must calculate rsqr for the while-loop condition. */
    rsqr = dyv_magnitude( conjgrad_r_ref( cg));
  }

  /* Select parameters. */
  /* Get best params from paramhist. */
  bestiter = dyv_argmin( devhist);
  x = dyv_array_ref( paramhist, bestiter);
  if (x == NULL) {
    my_errorf( "mk_lr_cgresult_cgdeveps: NULL param vec %d", bestiter);
  }

  result = mk_copy_dyv( x);
  free_dyv_array( paramhist);
  free_dyv( devhist);
  return result;
}

/* Run conjugate gradient. */
dyv *mk_lr_cgresult( lr_train *lrt, double unscaled_cgeps, double cgdeveps,
                     int maxiters, conjgrad *cg)
{
  dyv *result;

  if (cgdeveps > 0.0) {
    result = mk_lr_cgresult_cgdeveps( lrt, cgdeveps, maxiters, cg);
  }
  else {
    result = mk_lr_cgresult_cgeps( lrt, unscaled_cgeps, maxiters, cg);
  }

  if (Verbosity > 3) {
    fprintf_oneline_dyv( stdout, "    CG POS:", result, "\n");
  }

  return result;
}

void diag_precond( const dyv *v, dyv *result, void *userdata)
{
  /* Get diagonal     ( [1t]         )
                  diag( [--] W [1|X] )  = [ m_ii = Sum(x_ki^2 * w_k over k) ]
                      ( [Xt]         )
     In the sparse case, X is binary and x_ki^2 == x_ki, and the
     diagonal is [ m_ii = Sum(w_k over posrows_i) ].
     Preconditioning matrix is the diagonal matrix.  Multiply inverse
     of this matrix time v, which is an element-wise product.
  */
  int colidx;
  double divisor, val;
  ivec *posrows;
  dyv *w;
  lr_train *lrt;

  lrt = (lr_train *) userdata;

  if (lrt->X == NULL) {
    my_error( "diag_precond: dense problems not yet supported.");
  }

  w = lrt_w_ref( lrt);
  val = dyv_ref( v, 0);
  dyv_set( result, 0, val / dyv_sum( w));


  for (colidx=1; colidx < lrt->numatts; ++colidx) {
    posrows = spardat_attnum_to_posrows( lrt->X, colidx-1);
    divisor = dyv_partial_sum( w, posrows);
    val = dyv_ref( v, colidx);
    dyv_set( result, colidx, val / divisor);
  }

  return;
}

dyv *mk_lr_update_b_conjugate_gradient_helper( lr_train *lrt, double cgeps,
                                               double cgdeveps,
                                               int maxiters, int *iters,
                                               dyv *initx)
{
  int numatts;
  dyv *B, *x;
  cgopts *cgo;
  conjgrad *cg;

  numatts = lrt->numatts;

  B = mk_lr_XtWz_dyv( lrt);
  cgo = mk_cgopts_qspd( numatts, maxiters, -1.0 /* eps for runcg */,
                        B, initx,
                        (void *) lrt,
                        lr_cg_mk_copy_userdata,
                        lr_cg_free_userdata,
                        lr_cg_multA);
  free_dyv( B);

  /* Set up preconditioning. */
  /* set_cgopts_multMinv( cgo, diag_precond); */

  /* Remainder of setup work. */
  cg  = mk_conjgrad_from_cgopts( cgo);
  free_cgopts( cgo);

  /* Run conjugate gradient. */
  /* Course-grained method: runcg( cg); x = conjgrad_x_ref( cg); */
  /* Fine-grained method: */
  x = mk_lr_cgresult( lrt, cgeps, cgdeveps, maxiters, cg);

  /* Print iteration information. */
  *iters = cg->cgs->iterations;
  if (Verbosity >= 3) printf( "CG iterations=%d\n", cg->cgs->iterations);

  /* Done. */
  free_conjgrad( cg);
  return x;
}

/***********************************************************************/
/* LOGISTIC START                                                      */
/***********************************************************************/

int lr_deviance_test( lr_train *lrt, double epsilon,
                       double olddev, double *dev)
{
  /* Compute deviance */
  *dev = lr_train_deviance( lrt);

  /* Relative difference is small: |a-b|/a < epsilon  */
  if (fabs(*dev-olddev) < epsilon * *dev) return 1;

  /* If we reach this point, then we have not yet converged. */
  return 0;
}

/* Exactly one of X and ds should be NULL. */
lr_train *mk_lr_train( spardat *X, dym *factors, dyv *outputs,
                       dyv *initb, lr_options *opts)
{
  /* initb is copied into lr->b. */
  int converge, rc;
  int numiters, bestiter;
  double dev, olddev;
  dyv *devhist;
  lr_train *lrt;
  lr_state *bestlrs;
  lr_statearr *lrsarr;

  /* Create lr_train struct. */
  if (X != NULL) lrt = mk_lr_train_from_spardat( X, opts);
  else lrt = mk_lr_train_from_dym( factors, outputs, opts);

  /* Set initial value of model parameters, if desired. */
  if (initb != NULL) lr_train_overwrite_b( lrt, initb);

  /* Initialize our loop state */
  dev = -1000.0;
  lrsarr = mk_array_of_null_lr_states( opts->lrmax);
  devhist = mk_constant_dyv( opts->lrmax, FLT_MAX);

  /* START OF IRLS ITERATIONS */
  /* Iterate until the change in deviance is relatively small. */
  for (numiters=0; numiters < opts->lrmax; ++numiters) {

    /* Update olddev and iterate. */
    olddev = dev;
    rc = lr_train_iterate(lrt);

    /* Test for convergence. */
    lr_statearr_set( lrsarr, numiters, lrt->lrs);
    converge = lr_deviance_test( lrt, opts->lreps, olddev, &dev);
    dyv_set( devhist, numiters, dev);

    /* Print stuff. */
    if (Verbosity >= 1) printf( ".");
    if (Verbosity >= 3) {
      printf( "LR ITER %d: likesat: %g, likelihood: %g, deviance: %g\n",
	      numiters, lrt->likesat,
              lr_log_likelihood_from_deviance( dev, lrt->likesat), dev);
    }
    if (Verbosity >= 5) {
      /* Print all or most extreme attributes. */
        printf( "  Params, b0: %g\n", lrt->lrs->b0);
        fprintf_oneline_dyv( stdout, "  Params, b:", lrt->lrs->b, "\n");
    }

    if (converge) break;
    else if (rc == -2) break; /* Exceeded cgmax. */
    else if (am_isnan(dev)) break;
  }
  /* END OF ITERATIONS */

  /* Check state history for best holdout performance. */
  bestiter = dyv_argmin( devhist);
  bestlrs  = lr_statearr_ref( lrsarr, bestiter);
  free_lr_state( lrt->lrs);
  lrt->lrs = mk_copy_lr_state( bestlrs);
	if (converge) lrt->lrs->converged = converge;
  if (Verbosity == 1) printf( "\n");
  if (Verbosity >= 2) {
    printf( "CHOOSING ITERATION %d WITH DEVIANCE %g\n",
            bestiter, dyv_ref( devhist, bestiter));
  }
  if (Verbosity >= 2) {
    fprintf_oneline_dyv( stdout, "  devhist:", devhist, "\n");
  }

  /* Free state history. */
  free_lr_statearr( lrsarr);
  free_dyv( devhist);

  /* Done. */
  return lrt;
}

/***********************************************************************/
/* INPUT/OUTPUT                                                        */
/***********************************************************************/

void out_lr_predict( PFILE *f, lr_predict *lrp)
{
  int nump, i;
  double val;
  dyv *dv;

  nump = dyv_size( lrp->b) + 1;

  /* Copy b0, b into a single dyv. */
  dv = mk_dyv( nump);
  dyv_set( dv, 0, lrp->b0);
  for (i=1; i<nump; ++i) {
    val = dyv_ref( lrp->b, i-1);
    dyv_set( dv, i, val);
  }

  dyv_write( f, dv);

  free_dyv( dv);
  return;
}

lr_predict *mk_in_lr_predict( PFILE *f)
{
  int i, size;
  double val;
  dyv *dv, *b;
  lr_predict *lrp;

  lrp = AM_MALLOC( lr_predict);

  dv = mk_dyv_read( f);
  size = dyv_size( dv);

  lrp->b0 = dyv_ref( dv, 0);

  b = mk_dyv( size-1);
  for (i=1; i<size; ++i) {
    val = dyv_ref( dv, i);
    dyv_set( b, i-1, val);
  }
  lrp->b = b;

  free_dyv( dv);

  return lrp;
}
