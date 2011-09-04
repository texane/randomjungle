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


/* **
   File:            file.c
   Author:          Paul Komarek
   Date:            Thu Jul  7 14:17:51 EDT 2005
   Description:     File manipulation.
*/



#ifndef FILE_H
#define FILE_H


#include <stdio.h>

#if USE_ZLIB == 1
#include <zlib.h>
#endif

#include "amma.h"




/*  As of 2005-07-08, zlib's gzopen() will always compress files
   opened with mode "w".  The structure below is our replacement for
   the basic file pointer, to allow zlib integration yet still provide
   creation of uncompressed files.
*/
typedef struct pfile_struct {
  char *fname;
  char *mode;
  int gzip;
  void *fp;
} pfile;

typedef pfile PFILE;


/* This macro makes it easy to hide the zlib functions when zlib is
   not available. */
#if USE_ZLIB == 1
#  define SELECT_PFUNCTION(f, gzcmd, cmd) ((f->gzip) ? (gzcmd) : (cmd))
#else
#  define SELECT_PFUNCTION(f, gzcmd, cmd) (cmd)
#endif


/* This macro retrieves the file pointer and casts it appropriately. */
#if USE_ZLIB == 1
#  define PZFP(f) ((gzFile) (f->fp))
#else
/* gzFile does not exist since zlib.h not included. */
#  define PZFP(f) ((void *) (f->fp))
#endif

# define PFP(f) ((FILE *) (f->fp))



/**************************************************************************/
/* OVERRIDDEN FILE FUNCTIONS                                              */
/**************************************************************************/

PFILE *pfopen( const char *fname, const char *mode);
int pfclose( PFILE *pf);
int pfgetc( PFILE *f);
char *pfgets( char *buf, int size, PFILE *f);
int pfflush( PFILE *f);
int pfeof( PFILE *f);
void prewind( PFILE *f);
int pfprintf( PFILE *f, const char *format, ...);



/**************************************************************************/
/* VALUE-ADDED FILE FUNCTIONS                                             */
/**************************************************************************/

PFILE *safe_pfopen( const char *fname, const char *mode);
PFILE *sure_pfopen( const char *filename, const char *mode);
void sure_pfclose( PFILE *f, const char *filename);



#endif /* FILE_H */
