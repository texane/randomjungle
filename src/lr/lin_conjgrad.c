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
   File:        conjgrad.c
   Author:      Paul Komarek
   Created:     Wed Oct  2 13:56:31 EDT 2002
   Description: Linear conjugate gradient implementation.

   Copyright 2002, The Auton Lab, CMU
*/

#include <math.h>
#include "ambs.h"
#include "amma.h"
#include "amdyv.h"

#include "lin_conjgrad.h"



/* Update negative gradient, i.e.  the CG linear system residual.
   Note that cgs is a cgstate object, and the result is stored in
   cgs->r. */
/* See Nash and Sofer, "Linear and Nonlinear Programming", Section
   12.2, p. 385; and exercise 7 on p.390.  This produces equivalent
   results to cg_qspd_neggradf. */
#define UPDATE_NEGGRADF(cgs) \
  (dyv_madd( -(cgs)->alpha, (cgs)->Ap, (cgs)->r, (cgs)->r))



/*********************************************************************/
/* CONJGRAD                                                          */
/*********************************************************************/

conjgrad *mk_conjgrad( const cgopts *cgo, const cgstate *cgs)
{
  conjgrad *cg;
  cg = AM_MALLOC( conjgrad);
  cg->cgo = mk_copy_cgopts( cgo);
  cg->cgs = mk_copy_cgstate( cgs);
  return cg;
}

conjgrad *mk_conjgrad_from_cgopts( const cgopts *cgo)
{
  conjgrad *cg;
  cgstate *cgs;
  cgs = mk_cgstate( cgo);
  cg  = mk_conjgrad( cgo, cgs);
  free_cgstate( cgs);
  return cg;
}

conjgrad *mk_copy_conjgrad( const conjgrad *cg)
{
  conjgrad *newcg;
  newcg = mk_conjgrad( cg->cgo, cg->cgs);
  return newcg;
}

void free_conjgrad( conjgrad *cg)
{
  if (cg != NULL) {
    if (cg->cgo != NULL) free_cgopts( cg->cgo);
    if (cg->cgs != NULL) free_cgstate( cg->cgs);
    AM_FREE( cg, conjgrad);
  }
  return;
}

dyv *conjgrad_x_ref( conjgrad *cg)
{
  return cg->cgs->x;
}

dyv *conjgrad_r_ref( conjgrad *cg)
{
  return cg->cgs->r;
}

/*********************************************************************/
/* CGOPTS                                                            */
/*********************************************************************/

cgopts *mk_empty_cgopts( void)
{
  cgopts *cgo;

  cgo = AM_MALLOC( cgopts);

  cgo->userdata           = NULL;
  cgo->mk_copy_userdata   = NULL;
  cgo->free_userdata      = NULL;

  cgo->size               = -1;
  cgo->maxiter            = -1;
  cgo->smallenough        = -1.0;

  cgo->initx              = NULL;
  cgo->initr              = NULL;

  cgo->multMinv           = NULL;

  cgo->multA              = NULL;
  cgo->evalf              = NULL;

  return cgo;
}

cgopts *mk_copy_cgopts( const cgopts *cgo)
{
  cgopts *newcgo;

  newcgo = mk_empty_cgopts();

  if (cgo->userdata != NULL) {
    newcgo->userdata         = cgo->mk_copy_userdata( cgo->userdata);
    newcgo->mk_copy_userdata = cgo->mk_copy_userdata;
    newcgo->free_userdata    = cgo->free_userdata;
  }
  else {
    newcgo->userdata         = NULL;
    newcgo->mk_copy_userdata = NULL;
    newcgo->free_userdata    = NULL;
  }

  newcgo->size               = cgo->size;
  newcgo->maxiter            = cgo->maxiter;
  newcgo->smallenough        = cgo->smallenough;

  newcgo->initx = (cgo->initx==NULL) ? NULL : mk_copy_dyv( cgo->initx);
  newcgo->initr = (cgo->initr==NULL) ? NULL : mk_copy_dyv( cgo->initr);

  newcgo->multMinv           = cgo->multMinv;

  newcgo->multA              = cgo->multA;
  newcgo->evalf              = cgo->evalf;

  return newcgo;
}

void free_cgopts( cgopts *cgo)
{
  if (cgo != NULL) {
    if (cgo->userdata != NULL && cgo->free_userdata != NULL) {
      cgo->free_userdata( cgo->userdata);
    }

    if (cgo->initx != NULL) free_dyv( cgo->initx);
    if (cgo->initr != NULL) free_dyv( cgo->initr);

    AM_FREE( cgo, cgopts);
  }
  return;
}

/*********************************************************************/
/* CGSTATE                                                           */
/*********************************************************************/

cgstate *mk_empty_cgstate( void)
{
  cgstate *cgs;
  cgs = AM_MALLOC( cgstate);
  cgs->alpha = 1.0;
  cgs->p  = NULL;
  cgs->Ap = NULL;
  cgs->x  = NULL;
  cgs->r  = NULL;
  cgs->z  = NULL;
  cgs->rz = 0.0;
  cgs->oldrz = 0.0;
  cgs->iterations = 0;
  return cgs;
}

cgstate *mk_cgstate( const cgopts *cgo)
{
  cgstate *cgs;
  cgs = mk_empty_cgstate();

  /* Initialize alpha. */
  cgs->alpha = 1.0;

  /* Initialize p. */
  cgs->p = mk_zero_dyv( cgo->size);

  /* Initialize Ap. */
  cgs->Ap = mk_dyv( cgo->size);
  cgo->multA( cgs->p, cgs->Ap, cgo->userdata);

  /* Initialize x. */
  if (cgo->initx != NULL) cgs->x = mk_copy_dyv( cgo->initx);
  else cgs->x = mk_zero_dyv( cgo->size);

  /* Initialize r. */
  if (cgo->initr != NULL) cgs->r = mk_copy_dyv( cgo->initr);
  else cgs->r = mk_zero_dyv( cgo->size);

  /* Compute r using p, Ap, x, and r as needed. */
  UPDATE_NEGGRADF( cgs);

  /* Compute z from r. */
  cgs->z = mk_zero_dyv( cgo->size);
  if (cgo->multMinv != NULL) cgo->multMinv( cgs->r, cgs->z, cgo->userdata);
  else copy_dyv( cgs->r, cgs->z);

  /* Compute size of gradient at x. */
  cgs->rz = dyv_scalar_product( cgs->r, cgs->z);
  cgs->oldrz = 0.0;

  /* Initialize iterations. */
  cgs->iterations = 0;

  return cgs;
}

cgstate *mk_copy_cgstate( const cgstate *cgs)
{
  cgstate *newcgs;
  newcgs = AM_MALLOC( cgstate);
  newcgs->alpha = cgs->alpha;
  newcgs->x  = (cgs->x  == NULL) ? NULL : mk_copy_dyv( cgs->x);
  newcgs->p  = (cgs->p  == NULL) ? NULL : mk_copy_dyv( cgs->p);
  newcgs->Ap = (cgs->Ap == NULL) ? NULL : mk_copy_dyv( cgs->Ap);
  newcgs->r  = (cgs->r  == NULL) ? NULL : mk_copy_dyv( cgs->r);
  newcgs->z  = (cgs->z  == NULL) ? NULL : mk_copy_dyv( cgs->z);
  newcgs->rz = cgs->rz;
  newcgs->oldrz = cgs->oldrz;
  newcgs->iterations = cgs->iterations;
  return newcgs;
}

void free_cgstate( cgstate *cgs)
{
  if (cgs != NULL) {
    if (cgs->x  != NULL) free_dyv( cgs->x);
    if (cgs->p  != NULL) free_dyv( cgs->p);
    if (cgs->Ap != NULL) free_dyv( cgs->Ap);
    if (cgs->r  != NULL) free_dyv( cgs->r);
    if (cgs->z  != NULL) free_dyv( cgs->z);
    AM_FREE( cgs, cgstate);
  }
  return;
}

dyv *cgstate_x_ref( cgstate *cgs)
{
  return cgs->x;
}

void fprintf_cgstate( FILE *f, char *pre, const cgstate *cgs, char *post)
{
  fprintf( f, "%scgs: %p%s", pre, (void *) cgs, post);
  if (cgs != NULL) {
    fprintf( f, "%scgs->alpha: %f%s", pre, cgs->alpha, post);
    fprintf( f, "%scgs->p:     %p%s", pre, (void *) cgs->p, post);
    fprintf( f, "%scgs->Ap:    %p%s", pre, (void *) cgs->Ap, post);
    fprintf( f, "%scgs->x:     %p%s", pre, (void *) cgs->x, post);
    fprintf( f, "%scgs->r:     %p%s", pre, (void *) cgs->r, post);
    fprintf( f, "%scgs->z:     %p%s", pre, (void *) cgs->z, post);
    fprintf( f, "%scgs->rz:    %f%s", pre, cgs->rz, post);
    fprintf( f, "%scgs->oldrz: %f%s", pre, cgs->oldrz, post);
    fprintf( f, "%scgs->iterations: %d%s", pre, cgs->iterations, post);
  }
  return;
}

void fprintf_cgstate_verbose( FILE *f, char *pre, const cgstate *cgs,
			      char *post)
{
  fprintf( f, "%scgs: %p%s", pre, (void *) cgs, post);
  if (cgs != NULL) {
    fprintf( f, "%scgs->alpha: %f%s", pre, cgs->alpha, post);
    fprintf( f, "%scgs->p: ", pre);
    fprintf_oneline_dyv( f, "", cgs->p, post);

    fprintf( f, "%scgs->Ap:", pre);
    if (cgs->Ap == NULL) fprintf( f, "NULL\n");
    else fprintf_oneline_dyv( f, "", cgs->Ap, post);

    fprintf( f, "%scgs->x: ", pre);
    fprintf_oneline_dyv( f, "", cgs->x, post);
    fprintf( f, "%scgs->r: ", pre);
    fprintf_oneline_dyv( f, "", cgs->r, post);
    fprintf( f, "%scgs->z: ", pre);
    fprintf_oneline_dyv( f, "", cgs->z, post);
    fprintf( f, "%scgs->rz:    %f%s", pre, cgs->rz, post);
    fprintf( f, "%scgs->oldrz: %f%s", pre, cgs->oldrz, post);
    fprintf( f, "%scgs->iterations: %d%s", pre, cgs->iterations, post);
  }
  return;
}

/***********************************************************************/
/* CONJUGATE GRADIENT ITERATIONS                                       */
/***********************************************************************/

void runcg( conjgrad *cg)
{
  double rsqr;
  cgstate *cgs;

  /* Do conjugate gradient iterations. */
  cgs = cg->cgs;
  rsqr = -1000.0;
  while (1) {
    if (cgs->iterations > cg->cgo->maxiter) break;
    cgiter( cg);
    rsqr = dyv_scalar_product( cgs->r, cgs->r);
    if (sqrt(rsqr) < cg->cgo->smallenough) break;
  }

  return;
}

void cgiter( conjgrad *cg)
{
  double beta;
  cgstate *cgs;
  cgopts *cgo;

  cgs = cg->cgs;
  cgo = cg->cgo;

  /* This implementation is roughly the same as Nash and Sofer, "Linear
     and Nonlinear Programming", section 12.2, p. 385. */

  /* Check that residual is nonzero.  The user may be calling this
     function with their own termination criteria, but calculations
     fail if the residual is zero.  A near-zero residual is likely to
     cause problems, too, but we shouldn't second-guess the caller. */
  if (dyv_magnitude(cgs->r) == 0) {
    my_error( "cgiter: Error: called with zero residual: ||cgs->r|| = 0.");
  }

  /* Compute beta for search direction update. */
  if (cgs->iterations == 0) beta = 0.0;
  else beta = cgs->rz / cgs->oldrz;

  /* Compute next search direction. */
  dyv_madd( beta, cgs->p, cgs->z, cgs->p);           /* beta*p + z -> p */

  /* Compute Ap for functions with constant Hessian (constant, linear,
     or quadratic).  This is used by the incremental neggradf update. */
  cgo->multA( cgs->p, cgs->Ap, cgo->userdata);

  /* Compute line search coefficient exactly. */
  cgs->alpha = cgs->rz / dyv_scalar_product(cgs->p, cgs->Ap);

  /* Compute new test point x.  In terms of indexing our vectors, this
     is where the 'next' iteration starts. */
  dyv_madd( cgs->alpha, cgs->p, cgs->x, cgs->x);      /* alpha*p + x -> x */

  /* Compute new gradient at new x (quadratic residual). */
  UPDATE_NEGGRADF(cgs);

  /* Compute preconditioning vector from preconditioning matrix. */
  if (cgo->multMinv != NULL) cgo->multMinv( cgs->r, cgs->z, cgo->userdata);
  else copy_dyv( cgs->r, cgs->z);

  /* Compute size of gradient at x. */
  cgs->oldrz = cgs->rz;
  cgs->rz = dyv_scalar_product( cgs->r, cgs->z);

  cgs->iterations += 1;

  return;
}

/***************************************************************************/
/* CGOPTS/CGSTATE INITIALIZATION                                           */
/***************************************************************************/

void set_cgopts_initx( cgopts *cgo, const dyv *initx)
{
  if (cgo->initx != NULL) free_dyv( cgo->initx);
  if (initx == NULL) cgo->initx = NULL;
  else cgo->initx = mk_copy_dyv( initx);
  return;
}

void set_cgopts_initr( cgopts *cgo, const dyv *initr)
{
  if (cgo->initr != NULL) free_dyv( cgo->initr);
  if (initr == NULL) cgo->initr = NULL;
  else cgo->initr = mk_copy_dyv( initr);
  return;
}

/***************************************************************************/
/* PRECONDITIONING                                                         */
/***************************************************************************/

void set_cgopts_multMinv( cgopts *cgo, void (*multMinv)(const dyv *v, dyv *r,
							void *userdata))
{
  cgo->multMinv = multMinv;
  return;
}

/***************************************************************************/
/* FUNCTION AND GRADIENT EVALUATION                                        */
/***************************************************************************/

void set_cgopts_evalf( cgopts *cgo, double (*evalf)(dyv *x, void *userdata))
{
  cgo->evalf = evalf;
  return;
}				


/***************************************************************************/
/* MISCELLANEOUS OPTIONS                                                   */
/***************************************************************************/

void set_cgopts_multA( cgopts *cgo, void (*multA)(const dyv *v, dyv *result,
						  void *userdata))
{
  cgo->multA = multA;
  return;
}

/***************************************************************************/
/* LINEAR SYSTEMS                                                          */
/***************************************************************************/

/* Copies userdata. */
/* Create cgopts struct with external multA function, for instance
   for sparse matrices. */
cgopts *mk_cgopts_qspd( int size, int maxiter, double smallenough,
			const dyv *b, /* only used to initialize r */
			const dyv *initx,
			void *userdata,
			void *(*mk_copy_userdata)(const void *userdata),
			void (*free_userdata)(void *userdata),
			void (*multA)(const dyv *v, dyv *r, void *userdata))
{
  dyv *Ax, *initr;
  cgopts *cgo;
  cgo = mk_empty_cgopts();

  /* userdata stuff. */
  if (userdata != NULL) {
    cgo->userdata         = (void *) mk_copy_userdata( userdata);
    cgo->mk_copy_userdata = mk_copy_userdata;
    cgo->free_userdata    = free_userdata;
  }

  /* Algorithm parameters. */
  cgo->size = size;
  cgo->maxiter = maxiter;
  cgo->smallenough = smallenough;

  /* Set multA. */
  set_cgopts_multA( cgo, multA);


  set_cgopts_initx( cgo, initx);
  if (initx != NULL) {   /* initr = b-Ax. */
    Ax = mk_dyv( cgo->size);
    cgo->multA( initx, Ax, cgo->userdata);
    initr = mk_dyv_subtract( b, Ax);
    free_dyv( Ax);
  }
  else initr = mk_copy_dyv( b);
  set_cgopts_initr( cgo, initr);
  free_dyv( initr);

  return cgo;
}
