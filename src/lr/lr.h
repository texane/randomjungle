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
   File:        lr.h
   Author:      Paul Komarek
   Created:     Wed Mar  5 18:01:21 EST 2003
   Description: Logistic regression implmentation.

   Copyright 2003, The Auton Lab, CMU
*/

#ifndef LR_H
#define LR_H

#ifdef __cplusplus 
extern "C" {
#endif

#include "amdyv.h"
#include "amdym.h"
#include "spardat.h"

#include "lin_conjgrad.h"




/***********************************************************************/
/* LR_OPTIONS STRUCT                                                   */
/***********************************************************************/

typedef struct lr_options_struct {
  double rrlambda;
  double lreps;
  int    lrmax;

  double cgdeveps;
  double cgeps;
  int    cgmax;
  int    cgbinit;     /* Set automatically. */

  int    cgwindow;    /* Number of bad iterations allowed. */
  double cgdecay;     /* Factor worse than best-seen so far that is allowed. */
} lr_options;

lr_options *mk_lr_options(void);
lr_options *mk_copy_lr_options( lr_options *opts);
void free_lr_options( lr_options *opts);
void parse_lr_options( lr_options *opts, int argc, char **argv);
void check_lr_options( lr_options *opts, int argc, char **argv);
char *mk_string_from_lr_options( lr_options *opts, int argc, char **argv);
void fprintf_lr_options( FILE *f, char *pre, lr_options *opts,
			 int argc, char **argv, char *post);
void write_lr_options( lr_options *opts, int argc, char **argv);

/***********************************************************************/
/* LR_STATE                                                            */
/***********************************************************************/

typedef struct lr_state_struct {
  double  b0;    /* Constant factor. */
  dyv     *b;    /* (beta) Current estimates of regression coeffs.  Note
		    that b[1] is beta[1] in most descriptions of logistic
		    model.  The indices of b correspond to those of
		    X or M, whichever is not NULL. */

  /* The stuff below is used during iterations. */
  dyv     *n;    /* (eta)  Predicted values in linearized model. */
  dyv     *u;    /* (mu)   Predicted values. */
  dyv     *w;    /* Weights. */
  dyv     *z;    /* Adjusted outputs in linearized model. */
	int converged; /* does lr iterations converged */
} lr_state;


/* Prototype for mk_lr_state() occurs after lr_train is typedef'd. */
/* lr_state *mk_lr_state( lr_train *data, lr_options *opts); */
lr_state *mk_copy_lr_state( lr_state *lrs);
void fprintf_lr_state_brief( FILE *f, char *pre, lr_state *lrs);
void fprintf_lr_state( FILE *f, char *pre, lr_state *lrs);
void lr_state_overwrite_b( lr_state *lrs, dyv *initb);
void free_lr_state( lr_state *lrs);

/***********************************************************************/
/* LR_STATEARR STRUCT                                                   */
/***********************************************************************/

typedef struct lr_statearr_struct {
  int size;
  lr_state **arr;
} lr_statearr;

lr_statearr *mk_array_of_null_lr_states( int size);
lr_state *lr_statearr_ref( lr_statearr *lrsarr, int index);
void lr_statearr_set( lr_statearr *lrsarr, int index, lr_state *lrs);
void free_lr_statearr( lr_statearr *lrsarr);

/***********************************************************************/
/* LR_TRAIN STRUCT                                                     */
/***********************************************************************/

typedef struct lr_train_struct {
  /* Parameters. */
  lr_options *opts;

  /* One of the two following fields should be NULL. */
  const spardat *X;    /* sparse form of design matrix and outputs */
  dym   *M;            /* dense form of design matrix */
  int numatts;   /* Number of factors including constant factor. */
  int numrows;
  dyv     *y;     /* outputs as doubles */
  double likesat; /* likelihood of saturated model, i.e. with u=y. */

  lr_state *lrs;  /* LR state structure, used during iterations. */
} lr_train;

#define lrt_b0_ref(lrt)              ((lrt)->lrs->b0)
#define lrt_b_ref(lrt)               ((lrt)->lrs->b)
#define lrt_n_ref(lrt)               ((lrt)->lrs->n)
#define lrt_u_ref(lrt)               ((lrt)->lrs->u)
#define lrt_w_ref(lrt)               ((lrt)->lrs->w)
#define lrt_z_ref(lrt)               ((lrt)->lrs->z)

#define lrt_b0_set(lrt,val)          ((lrt)->lrs->b0 = val)
#define lrt_b_set(lrt,dv)            ((lrt)->lrs->b = dv)
#define lrt_n_set(lrt,dv)            ((lrt)->lrs->n = dv)
#define lrt_u_set(lrt,dv)            ((lrt)->lrs->u = dv)
#define lrt_w_set(lrt,dv)            ((lrt)->lrs->w = dv)
#define lrt_z_set(lrt,dv)            ((lrt)->lrs->z = dv)
#define lrt_converged_set(lrt,conv)            ((lrt)->lrs->converged = conv)

/* Defined here because it depends on lr_train. */
lr_state *mk_lr_state( lr_train *data, lr_options *opts);


lr_train *mk_lr_train_from_dym( dym *factors, dyv *outputs, lr_options *opts);
lr_train *mk_lr_train_from_spardat( spardat *X, lr_options *opts);
lr_train *mk_copy_lr_train( const lr_train *source);
void free_lr_train( lr_train *lrt);
void fprintf_lr_train( FILE *f, char *pre, lr_train *lrt);
void lr_train_overwrite_b( lr_train *lrt, dyv *initb);
int lr_train_iterate( lr_train *lrt);

double lr_train_deviance( lr_train *lrt);

void lr_train_split_b( dyv *b, lr_train *lrt);
void lr_train_join_b( lr_train *lrt, double b0, dyv *b);

void lr_train_update_w( lr_train *lrt);
void lr_train_update_z( lr_train *lrt);
int  lr_train_update_b( lr_train *lrt);
void lr_train_update_n( lr_train *lrt);
void lr_train_update_u( lr_train *lrt);
int lr_train_iterate( lr_train *lrt);

/***********************************************************************/
/* LR_PREDICT                                                          */
/***********************************************************************/

typedef struct lr_predict_struct {
  double b0;
  dyv *b;
} lr_predict;

lr_predict *mk_lr_predict( double b0, dyv *b);

lr_predict *mk_copy_lr_predict( lr_predict *lrp);
void free_lr_predict( lr_predict *lrp);
double lr_predict_predict( ivec *posatts, dyv *attvals, lr_predict *lrp);

/***********************************************************************/
/* UTILITY                                                             */
/***********************************************************************/

/* Predictions. */
double lr_prediction( double b0, dyv *b, ivec *posatts, dyv *attvals);

/* Computing n (eta). */
void lr_compute_n_from_spardat( const spardat *X, double b0, dyv *b, dyv *n);
void lr_compute_n_from_dym( const dym *M, double b0, dyv *b, dyv *n);
void lr_n_from_spardat( const spardat *X, double b0, dyv *b, dyv *n);
void lr_n_from_dym( const dym *M, double b0, dyv *b, dyv *n);

/* Computing u (mu). */
void lr_compute_u_from_n( dyv *n, dyv *u);

/* Computing b (beta). */
dyv *mk_lr_update_b_conjugate_gradient_helper( lr_train *lrt, double cgeps,
                                               double cgdeveps,
                                               int maxiters, int *iters,
                                               dyv *initx);

/* Likelihood and Deviance. */
double lr_log_likelihood_basic( dyv *y, dyv *u);
double lr_log_likelihood_from_deviance( double deviance, double likesat);
double lr_deviance_from_log_likelihood( double likelihood, double likesat);

double lr_deviance_basic( dyv *y, dyv *u);
double lr_deviance_from_spardat_b( const spardat *X, dyv *y, double b0,
                                   dyv *b);
double lr_deviance_from_dym_b( const dym *M, dyv *y, double b0, dyv *b);
double lr_deviance_from_cg( lr_train *lrt, conjgrad *cg);

/* Exactly one of posatts and attvals should be NULL. */
dyv *mk_lr_XtWXv_dyv( const lr_train *lrt, const dyv *v);
dyv *mk_lr_XtWz_dyv( const lr_train *lrt);

/***********************************************************************/
/* CONJGRAD HELPERS                                                    */
/***********************************************************************/

/* Copy [b0, b1, ..., bn] into lrt->b0 and lrt->b. */
void lr_copy_full_b_to_lr_train( dyv *sourceb, lr_train *lrt);
void lr_copy_cgs_x_to_lr_train( cgstate *cgs, lr_train *lrt);
void *lr_cg_mk_copy_userdata( const void *userdata);
void lr_cg_free_userdata( void *userdata);
void lr_cg_multA( const dyv *v, dyv *result, void *userdata);

/***********************************************************************/
/* LR LEARN AND PREDICT                                                */
/***********************************************************************/

/* Exactly one of X and ds should be NULL. */
lr_train *mk_lr_train( spardat *X, dym *factors, dyv *outputs,
                       dyv *initb, lr_options *opts);

/***********************************************************************/
/* INOUT                                                               */
/***********************************************************************/

void out_lr_predict( PFILE *f, lr_predict *lrp);
lr_predict *mk_in_lr_predict( PFILE *f);

#ifdef __cplusplus
}
#endif

#endif
