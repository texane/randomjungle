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
   File:        standard.h
   Author:      Andrew W. Moore
   Created:     Wed Jun  7 15:29:06 EDT 1995
   Description: Includer of standard libraries

   Copyright (C) 1995, Andrew W. Moore
*/

#ifndef STANDARD_H
#define STANDARD_H

#ifndef UNIX_TTY_PLATFORM
#ifdef WIN32
/* Precisely one of the following constants must be #defined and the
   others must not be #defined. */
/* #define UNIX_TTY_PLATFORM */
 /*#define UNIX_XW_PLATFORM*/
/* #define PC_TTY_PLATFORM */
#define PC_TTY_PLATFORM
#else
#define UNIX_XW_PLATFORM
#endif
#endif

/* 
   UNIX_TTY_PLATFORM -- this code will run under Unix and not produce
                        any graphics.  Any graphics calls will be
                        ignored, or possibly produce a printed warning
   UNIX_XW_PLATFORM  -- this code will run under Unix and is free to
                        produce graphics using damut/amgr.h graphics
                        routines
   PC_TTY_PLATFORM   -- this code will run on PCs (compiled by Visual C++)
                        and must not produce any graphics.  Any
			graphics calls will be ignored, or possibly
			produce a printed warning, or possibly pop up
                        graphics windows, though all control will remain
                        with the standard input.
   PC_MVIS_PLATFORM  -- this code will run on PCs, compiled by Visual
                        C++ as part of a single-document project.  It
                        requires the use of Mary's EYE files.  There
                        is no stdio/stderr output or input.  User
                        communication uses expos, apicts, and aform
			interface functions.  The programmer is free
			to produce graphics using damut/amgr.h graphics
                        routines.
*/

/* The following #defined constants are created as a function of which
   one of the above list is defined.  */
#ifdef UNIX_TTY_PLATFORM
#define UNIX_PLATFORM
#define TTY_PLATFORM
#endif /* UNIX_TTY_PLATFORM */

#ifdef PC_TTY_PLATFORM
#define TTY_PLATFORM
#define PC_PLATFORM
#endif /* PC_TTY_PLATFORM */

#ifdef UNIX_XW_PLATFORM
#define UNIX_PLATFORM
#endif /* UNIX_XW_PLATFORM */

#ifdef PC_MVIS_PLATFORM
#define PC_PLATFORM
#endif /* PC_MVIS_PLATFORM */

/* Use for exporting symbols from dlls.  Needed on MSVC.  Also adding code
   for gcc so only export needed symbols where possible. */
#ifdef PC_PLATFORM
   /* define flag to use for exporting symbols from a dll on windows*/
#define AMEXPORT __declspec(dllexport)
#elif __GNUC__ > 3 ||				\
  (__GNUC__ == 3 && (__GNUC_MINOR__ > 3 ||      \
		     (__GNUC_MINOR__ == 3)))
# define AMEXPORT __attribute__ ((visibility("default")))
#else
# define AMEXPORT 
#endif

#ifdef _MSC_VER
   /* !!!!  Comment out the pragma statement on non-Visual C++ platforms.
   I could not put it within a #ifdef PC_MVIS_PLATFORM because xdamut
   and xambl no longer set this in their Build Settings. 

   The pragma statement should disable the level 4 warning: "unreferenced 
   inline function has been removed" */
#pragma warning( disable : 4514)
   /* Fixing some microsoft brain-damage here. */
#define snprintf _snprintf
#endif

/* Directory separation string  e.g. "/" for Unix, "\" for DOS */
#ifdef UNIX_PLATFORM
#define DIRSEP "/"
#else
#define DIRSEP "\\"
#endif /* UNIX_PLATFORM */


#ifdef UNIX_PLATFORM
#ifndef _SYS_TIMES_H
#include <sys/times.h>
#ifndef _SYS_TIMES_H
#define _SYS_TIMES_H 1
#endif /* sys_times_h sanity check */
#endif /* _SYS_TIMES_H */
#endif /* UNIX_PLATFORM */

#ifdef UNIX_PLATFORM
#include <sys/time.h>
#include <signal.h>
#include <unistd.h>
#endif

#ifdef PC_PLATFORM
#include <winsock.h>
#endif

#ifndef STDIO_H
#include <stdio.h>
#ifndef STDIO_H
#define STDIO_H
#endif /* STDIO_H inner */
#endif /* STDIO_H outer */

#ifndef STRING_H
#include <string.h>
#ifndef STRING_H
#define STRING_H
#endif /* STRING_H inner */
#endif /* STRING_H outer */

#ifndef STDLIB_H
#include <stdlib.h>
#ifndef STDLIB_H
#define STDLIB_H
#endif /* STDLIB_H inner */
#endif /* STDLIB_H outer */

#ifndef TIME_H 
#include <time.h>
#ifndef TIME_H 
#define TIME_H
#endif /* TIME_H inner */
#endif /* TIME_H outer */

#ifndef TYPES_H 
#include <sys/types.h>
#ifndef TYPES_H 
#define TYPES_H 
#endif /* TYPES_H inner */
#endif /* TYPES_H outer */

#ifndef SYSTIMEB_H 
#include <sys/timeb.h>
#ifndef SYSTIMEB_H
#define SYSTIMEB_H
#endif /* inner */
#endif /* outer */

#ifndef LIMITS_H
#include <limits.h>
#ifndef LIMITS_H
#define LIMITS_H
#endif /* inner */
#endif /* outer */

#ifdef OLD_SPR
#ifndef ADGUI_H
#include "adgui.h"
#ifndef ADGUI_H
#define ADGUI_H
#endif /* inner */
#endif /* outer */
#endif

#ifndef MATH_H
#ifdef PC_PLATFORM
#define _USE_MATH_DEFINES 1
#endif /* PC_PLATFORM */
#include <math.h>
#ifndef MATH_H
#define MATH_H
#endif /* MATH_H inner */
#endif /* MATH_H outer */

#ifdef PC_PLATFORM
#ifdef NDEBUG
#ifndef AMFAST
#define AMFAST
#endif
#endif
#endif

#ifndef STDARG_H
#include <stdarg.h>
#ifndef STDARG_H
#define STDARG_H
#endif
#endif

#ifdef PC_PLATFORM
#define strcasecmp _stricmp
#endif

#endif /* standrd_h */
