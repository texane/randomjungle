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
   File:        conjgrad.h
   Author:      Paul Komarek
   Created:     Wed Oct 30 18:27:21 EST 2002
   Description: Linear conjugate gradient ipmlementation.

   Copyright 2002, The Auton Lab, CMU
*/

#ifndef __CONJGRAD_H__
#define __CONJGRAD_H__



#include "ambs.h"
#include "amdyv.h"



/*********************************************************************/
/* STRUCTURES                                                        */
/*********************************************************************/

typedef struct conjgrad_struct {
  struct cgopts_struct  *cgo;
  struct cgstate_struct *cgs;
} conjgrad;

typedef struct cgstate_struct {
  double alpha;  /* line-search coefficient */
  dyv    *x;     /* current guess at minimizer */
  dyv    *p;     /* current search direction */
  dyv    *Ap;    /* change in the gradient per unit step in direction p */
  dyv    *r;     /* -grad f (== b-Ax for quadratic f) */
  dyv    *z;     /* M^(-1).r, where M is preconditioning matrix. */
  double rz;     /* r.z */
  double oldrz;  /* r[i-1].z[i-1]. */
  
  int iterations;
} cgstate;

typedef struct cgopts_struct {
  void *userdata;
  void *(*mk_copy_userdata)( const void *userdata);
  void (*free_userdata)( void *userdata);
  
  /* Algorithm parameters. */
  int size;            /* Number of variables. */
  int maxiter;
  double smallenough;  /* Residual bound for cg iterations. */
  
  /* Initial state. */
  dyv *initx;
  dyv *initr; /* initr might not do what you think it does; see mk_cgstate()
		 for clarification. */
  
  /* Additional functions which may or may not be required depending upon
     how you initialized cgopts, your choice of linesearch, whether or not
     you call fprintf_cg_linesearch(), etc.  Some of these can be chosen from
     helpers provided further below. */

 /* Function to multiply constant Hessian matrix by a vector. */
  void (*multA)(const dyv *v, dyv *result, void *userdata);

  void (*multMinv)(const dyv *v, dyv *result, void *userdata);

  double (*evalf)(dyv *x, void *userdata);                       /* f(x)     */
  
} cgopts;


/*********************************************************************/
/* CONJGRAD                                                          */
/*********************************************************************/

conjgrad *mk_conjgrad( const cgopts *cgo, const cgstate *cgs);
conjgrad *mk_conjgrad_from_cgopts( const cgopts *cgo);
conjgrad *mk_copy_conjgrad( const conjgrad *cg);
void free_conjgrad( conjgrad *cg);
dyv *conjgrad_x_ref( conjgrad *cg);
dyv *conjgrad_r_ref( conjgrad *cg);

/*********************************************************************/
/* CGOPTS                                                            */
/*********************************************************************/

cgopts *mk_empty_cgopts( void);
cgopts *mk_copy_cgopts( const cgopts *cgo);
void free_cgopts( cgopts *cgo);

/*********************************************************************/
/* CGSTATE                                                           */
/*********************************************************************/

cgstate *mk_empty_cgstate( void);
cgstate *mk_cgstate( const cgopts *cgo);
cgstate *mk_copy_cgstate( const cgstate *cgs);
void free_cgstate( cgstate *cgs);
dyv *cgstate_x_ref( cgstate *cgs);
void fprintf_cgstate( FILE *f, char *pre, const cgstate *cgs, char *post);
void fprintf_cgstate_verbose( FILE *f, char *pre, const cgstate *cgs,
			      char *post);

/***********************************************************************/
/* CONJUGATE GRADIENT ITERATIONS                                       */
/***********************************************************************/

void runcg( conjgrad *cg);
void cgiter( conjgrad *cg);

/***************************************************************************/
/* CGOPTS/CGSTATE INITIALIZATION                                           */
/***************************************************************************/

void set_cgopts_initx( cgopts *cgo, const dyv *initx);
void set_cgopts_initr( cgopts *cgo, const dyv *initr);

/***************************************************************************/
/* PRECONDITIONING                                                         */
/***************************************************************************/

void set_cgopts_multMinv( cgopts *cgo, void (*multMinv)(const dyv *v, dyv *r,
							void *userdata));

/***************************************************************************/
/* FUNCTION AND GRADIENT EVALUATION                                        */
/***************************************************************************/

void set_cgopts_evalf( cgopts *cgo, double (*evalf)(dyv *x, void *userdata));

/***************************************************************************/
/* MISCELLANEOUS OPTIONS                                                   */
/***************************************************************************/

void set_cgopts_multA( cgopts *cgo, void (*multA)(const dyv *v, dyv *result,
						  void *userdata));

/***************************************************************************/
/* LINEAR SYSTEMS                                                          */
/***************************************************************************/

/* Create cgopts struct with external multA function, for instance
   for sparse matrices.  Copies userdata. */
cgopts *mk_cgopts_qspd( int size, int maxiter, double smallenough,
			const dyv *b,
			const dyv *initx,
			void *userdata,
			void *(*mk_copy_userdata)(const void *userdata),
			void (*free_userdata)(void *userdata),
			void (*multA)(const dyv *v, dyv *r, void *userdata));


#endif
