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
   File:        amma.h
   Author:      Rahul Kalaskar
   Created:     Fri May 20 2005
   Description: Stripped down verion of Auton memory allocator

*/

#ifndef _AMMA_LITE_H_
#define _AMMA_LITE_H_

#include<stdlib.h> /* malloc, free */

#define AM_MALLOC(t) ((t *) malloc(sizeof(t)))
#define AM_MALLOC_ARRAY(t,len) ((t *) malloc((len) * sizeof(t)))

#define AM_FREE(thing,type) free(thing)
#define AM_FREE_ARRAY(thing,type,len) free(thing)

#define am_malloc(len) ((char*) malloc(len))
#define am_free(thing,len) free(thing)






#endif 
