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



#include <stdio.h>

#if USE_ZLIB == 1
#include <zlib.h>
#endif

#include "ambs.h"
#include "file.h"




/**************************************************************************/
/* BASIC pfile FUNCTIONS, consider them private.                          */
/**************************************************************************/

pfile *mk_pfile( const char *fname, const char *mode, int retry)
{
  int namelen;
  pfile *pf;
  pf = AM_MALLOC( pfile);

  /* Check name and mode. */
  namelen = strlen( fname);
  if (namelen <= 0) {
    my_errorf( "mk_pfile: fname '%s' appears unreasonable, length is %d.",
               namelen);
  }
  if (strlen(mode) <= 0) {
    my_errorf( "mk_pfile: fname '%s' appears unreasonable, length is %d.",
               strlen(mode));
  }

  /* Fill in structure. */
  pf->fname = AM_MALLOC_ARRAY( char, strlen(fname)+1);
  strcpy( pf->fname, fname);

  pf->mode = AM_MALLOC_ARRAY( char, strlen(fname)+1);
  strcpy( pf->mode, mode);


  /* Set gzip flag.  We cannot call am_string's string_has_suffix(),
     so we hack a .gz test to deterimines compression-ness. */
  pf->gzip = 0;
  if (fname[namelen-1] == 'z' \
      && fname[namelen-2] == 'g' \
      && fname[namelen-3] == '.') {
    /* Ends in .gz. */
#if USE_ZLIB == 1
    pf->gzip = 1;
#else
    fprintf( stderr,
             "\n\n"
             "************************************************************\n"
             "WARNING:  zlib support is not enabled in this build.\n"
             "We will still try to %s the file '%s',\n"
             "but the results may be unexpected.\n"
             "************************************************************\n"
             "\n\n",
             strchr(pf->mode, 'r') == NULL ? "write to" : "read from",
             pf->fname);
    pf->gzip = 0;
#endif

  } /* end if(fname ends in .gz) {...} */


  /* File pointer, filled in by pfopen(). */
  pf->fp = NULL;

  return pf;
}

void free_pfile( pfile *pf)
{
  if (pf != NULL) {
    if (pf->fname != NULL) AM_FREE_ARRAY( pf->fname, char, strlen(fname)+1);
    if (pf->mode != NULL) AM_FREE_ARRAY( pf->mode, char, strlen(fname)+1);
    AM_FREE( pf, pfile);
  }
  return;
}




/**************************************************************************/
/* OVERRIDDEN FILE FUNCTIONS                                              */
/**************************************************************************/

PFILE *pfopen( const char *fname, const char *mode)
{
  pfile *pf;
  pf = mk_pfile( fname, mode, 0);
  pf->fp = (void *) SELECT_PFUNCTION( pf,
                                      gzopen( fname, mode),
                                      fopen( fname, mode));
  return pf;
}

int pfclose( PFILE *pf)
{
  int rc;

  rc = -1;

  if (pf == NULL) my_error( "pfclose: pf == NULL.");
  else if (pf->fp == NULL) my_error( "pfclose: pf->fp == NULL.");
  else {
    rc = SELECT_PFUNCTION( pf, gzclose(PZFP(pf)), fclose(PFP(pf)));
    if (rc == 0) free_pfile( pf);
  }

  return rc;
}

int pfgetc( PFILE *f)
{
  return SELECT_PFUNCTION( f, gzgetc(PZFP(f)), fgetc(PFP(f)));
}

char *pfgets( char *buf, int size, PFILE *f)
{
  return SELECT_PFUNCTION( f,
                           gzgets(PZFP(f), buf, size),
                           fgets(buf, size, PFP(f)));
}

int pfflush( PFILE *f)
{
  return SELECT_PFUNCTION( f, gzflush(PZFP(f), Z_SYNC_FLUSH), fflush(PFP(f)));
}

int pfeof( PFILE *f)
{
  return SELECT_PFUNCTION( f, gzeof(PZFP(f)), feof(PFP(f)));
}

void prewind( PFILE *f)
{
  SELECT_PFUNCTION( f, gzrewind(PZFP(f)), rewind(PFP(f)));
  return;
}

int pfprintf( PFILE *f, const char *format, ...)
{
  int count, bufsize, rc;
  char *buf;
  va_list ap;

  /* zlib does not provide a gzfprintf() function according to the
     online docs as of 2005-07-08. */

  /* WARNING: beware of vsnprintf() for glibc earlier than 2.1.  Their
     return codes did not necessarily mean the same thing as the
     C99 standard specifies.  Our code below will die unncessarily
     but harmlessly if this occurs, if the string requires more space
     than our initial buffer provides.  We are unlikely to find
     a pre-2.1 glibc platform anymore. 
  */

  /* WARNING: Microsoft has chosen not to support C99 as of
     2005-07-19.  They do have vsnprintf(), but only by the name
     _vsnprintf().  Unfortunately, they do not use the C99 return
     value in their _vsnprintf().  Thus our code may die unnecessarily
     but harmlessly under Windows under the same circumstances as
     it would under pre-2.1 glibc.
  */

  /* Print into buffer.  If vsnprintf() needs more space, we
     re-allocate and try again (see WARNINGs, above).  Because of
     Windows (see WARNING, above), we should make this buffer somewhat
     larger than we would like. */
  bufsize = 16 * 4096;
  buf = AM_MALLOC_ARRAY( char, bufsize);

  va_start(ap, format);
#ifndef WIN32
  count = vsnprintf( buf, bufsize, format, ap);
#else
  count = _vsnprintf( buf, bufsize, format, ap);
#endif
  va_end(ap);

  /* Check for overflow or error. */
  if (count < 0 ) {
    my_errorf( "pfprintf: encountered error %d from vsnprintf.", count);
  }
  else if (count >= bufsize) {
    /* Re-allocate. */
    AM_FREE_ARRAY( buf, char, bufsize);
    bufsize = count+1;
    buf = AM_MALLOC_ARRAY( char, bufsize);

    va_start(ap, format);
#ifndef WIN32
    count = vsnprintf( buf, bufsize, format, ap);
#else
    count = _vsnprintf( buf, bufsize, format, ap);
#endif
    va_end(ap);
    
    if (count < 0 ) {
      my_errorf( "pfprintf: encountered error %d from vsnprintf on the second "
                 "pass.", count);
    }
    else if (count >= bufsize) my_error( "pfprintf: internal error.");
  }

  /* We can finally print the string. */
  rc = SELECT_PFUNCTION( f,
                         gzprintf(PZFP(f), "%s", buf),
                         fprintf(PFP(f), "%s", buf));

  AM_FREE_ARRAY( buf, char, bufsize);

  return rc;
}



/**************************************************************************/
/* VALUE-ADDED FILE FUNCTIONS                                             */
/**************************************************************************/

PFILE *safe_pfopen( const char *fname, const char *access)
{
  PFILE *s;

  s = pfopen( fname, access);
  if ( s->fp == NULL ) {
    perror( "safe_fopen");
    printf( "Could not fopen('%s','%s')\n", fname, access);
    my_error( "safe_fopen() failed");
  }

  return s;
}

PFILE *sure_pfopen( const char *filename, const char *mode)
{
  PFILE *f;
  while (1) {
    f = pfopen( filename, mode);
    if (f->fp != NULL) break;
    fprintf( stderr, "sure_fopen: cannot open '%s' with mode '%s': ",
             filename, mode);
    perror( "");
    fprintf( stderr, "  Press return to try again, ctrl-c to abort.\n");
    fgetc( stdin);
  }
  return f;
}

void sure_pfclose( PFILE *f, const char *filename)
{
  while ( pfclose(f)) {
    fprintf( stderr, "sure_fclose: failed to close file '%s'", filename);
    perror( "");
    fprintf( stderr, "  Until the file is closed successfully, data may "
	     "be lost.\n");
    fprintf( stderr, "  Press return to try again, ctrl-c to abort.\n");
    fgetc( stdin);
  }
  return;
}
