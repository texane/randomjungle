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


/** @file ambs.h
 @verbatim
   File:        ambs.h
   Author:      Andrew W. Moore
   Created:     Sat Sep 12 15:22:40 EDT 1992
   Updated:     June 25 97
   Description: Basic

   Copyright 1996, Schenley Park Research
@endverbatim

ambs.c and ambs.h contain our set of basic utility functions--- the kind
of set of utilities that most C programming projects have. If you're working
on the Auton project feel free but not obliged to use them. Try to use
them in preference to your own version of the same thing if possible. If
something's missing, please feel free to add it, and then email awm+lab@cs
to tell us.

This file defines bool to be the boolean datatype, (with a typedef to int)
and #defines TRUE and FALSE to 1 and 0 respectively. It does a couple of
other convenient #defines too.
*/

#ifndef AMBS_H
#define AMBS_H

#include "standard.h"

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


#ifndef BOOL_DEFINED
typedef unsigned char auton_u_char;
#ifndef __cplusplus 
#define bool auton_u_char
#endif
#define BOOL_DEFINED
#endif

#ifndef DOUBLE_PTR_DEFINED
typedef double *double_ptr;
#define DOUBLE_PTR_DEFINED
#endif

#ifndef CHAR_PTR_DEFINED
typedef char *char_ptr;
#define CHAR_PTR_DEFINED
#endif

#ifndef CHAR_PTR_PTR_DEFINED
typedef char_ptr *char_ptr_ptr;
#define CHAR_PTR_PTR_DEFINED
#endif

#ifndef PI
#define PI 3.1415926535897932384626
#endif
#define ROOT_TWO 1.4142135623730952
#define ROOT_OF_TWO_PI 2.506628274631000242

#ifndef LOG_2
#define LOG_2 0.693147180559945286
#endif
#ifndef INVLOG_2
#define INVLOG_2 1.442695040888963407
#endif

#define LOG_TWO_PI 1.837877066409345339

#define LOG_ROOT_TWO_PI 0.918938533204672670

#define EQ_PTR(p1,p2) (((long) p1) == ((long) p2))

extern double Verbosity;

extern bool htmlp;
/**
Use this if you wish to print out diagnostics in development mode, but not
in normal use. Conditionally print out if Verbosity is greater than a certain
value. Informally, the user should expect that with Verbosity 0.0, only
useful stuff appears, Verbosity = 1.0 adds a bit more, Verbosity = 10.0
prints a fair amount and occasionally even stops and waits for a key,
Verbosity = 100.0 prints all kinds of rubbish.
*/

extern bool SuppressMessages;

/** specify the user interactin mode.  This determines the 
   behavior of functions that ask for user input: ex. wait_for_key().
   The default mode is kShell.
*/
typedef enum ui_mode {
    kUI_Shell = 0,  /**< default ui mode. user interaction happens at the shell.*/
    kUI_NoInput,    /**< skip user interaction requests (ex. for batch runs)*/
    kUI_GUIwx,      /**< special messages will be output that the applic_wxgui 
		      knows how to interpret for user interaction*/
    kNumUIModes
} ui_mode;

void Global_set_ui_mode(ui_mode mode);
int Global_get_ui_mode(void);

/** specify the name of the logfile to write log messages to.  Can
    also be set via the AM_LOGFILE env var.
*/
bool set_logfile_name(const char *str);

/** Write a message to the logfile (specified by AM_LOGFILE env var
    or by calling set_logfile_name().  All of the below info, warn,
    error functions will also write to the logfile.
*/
void my_logmsg(const char *format, ...);

/** Output a message, and then wait for a key
*/
void my_info(char *string);
void my_infof(char *format, ...);

/** Signal a warning with the message, and then wait for a key
*/
void my_warning(char *string);
void my_warning_no_wait(char *string);

/** Signal an *** Auton Error: with the message, and then stop. Actually
    waits for a key before exitting so that in the debugger you get hit
    cntrl-c and examine the program stack
*/
void my_error(const char *string);

/** You should should use my_assert( boolean statement ) to assert
   facts that you want checked for code correctness. */
/** Make sure the assert statement has no side effects that you need
   because it won't be executed in fast mode */

/** If you want the check to also happen in FAST mode, call this instead... */
void my_assert_always(bool b);

/** Make sure the assert statement has no side effects that you need
   because it won't be executed in fast mode */

/** Uncomment the next line if you want the program to check ALL assertions
   even in fast mode */

/** #define ALWAYS_ASSERT */

/** In your code you can assume that if DOING_ASSERTS is defined then
   my_assert will execute and that the person who compiled the code
   wants all assertions checked.

   If DOING_ASSERTS is undefined then the person who compiled does
   not want assertions checked. */

#ifdef ALWAYS_ASSERT
#define my_assert(b) my_assert_always(b)
#define DOING_ASSERTS
#else
#ifdef AMFAST
#define my_assert(b) /**< do nothing */
#else
#define my_assert(b) my_assert_always(b)
#define DOING_ASSERTS
#endif
#endif

/** like my_error, but accepts printf-like params */
void my_errorf(const char *format, ...);

/** Waits for user to hit return before continuing. Also gives user
   the option of doing some other things (including turning off
   future wait_for_keys(). In a windows environment, MIGHT do a
   dialog box instead of a console print */
void wait_for_key(void);
/** This version can't be turned off */
void really_wait_for_key(void);

/** Show a message to the user and wait for a text response */
/** (Use this in place of mk_string_from_user (am_string.c) if
   your code will be used by applic. 
*/
char *mk_get_user_input(const char *message, const char *def_reply);

/** Prints the message where the user will see it. In a windows
   environment might write this into a status window instead of
   the regular console. */
void status_report(char *message);

char *safe_malloc(unsigned size);
/**
  Very very basic malloc. Please don't use this. Use am_malloc defined
  in amma.h
*/

bool eq_string(const char *s1, const char *s2);
bool eq_string_with_length(const char *s1, const char *s2, int n);
/**
   Returns TRUE if and only if the strings contain the same characters.
   The first reads up to a NULL character.  The 2nd reads the 1st n chars
*/

double pd_random_double(void);

/** Uses the time to set a random seed for the random() functions above */
void am_randomize(void);

/** Following seeds the random number generator */
void am_srand(int seed);

void push_current_am_srand_state(void);
void pop_current_am_srand_state(void);

/**  Returns a random integer uniformly >= 0 and < n. */
int int_random(int n);

double range_random(double lo, double hi);
/**
   Returns a random double x uniformly s.t. lo <= x <= hi
*/

/** Generate a gaussian random number mean = 0, variance
   (= sdev) = 1.0 */
double gen_gauss(void);


/** @name Basic Type Utilities
   Following functions have obvious properties. In some cases standard macro
   implementations exists. You can use these if you fear strange behavior
   of like to have things typechecked properly, or you fear they might
   not be implemented in all c standard libraries
*/
/*@{*/
double real_min(double x, double y);
double real_max(double x, double y);
double real_abs(double x);
double real_square(double x);
double real_cube(double x);

int int_min(int x, int y);
int int_max(int x, int y);
int int_abs(int x);
int int_square(int x);
int int_cube(int x);

long long_min(long x, long y);
long long_max(long x, long y);
long long_abs(long x);
long long_square(long x);
long long_cube(long x);

bool am_isinf( double val);
bool am_isnan( double val);
bool am_isnum( double val);
double am_copysign(double x, double y);

char *input_string(char *mess, char *string, int string_size);
double input_realnum(char *mess);
int input_int(char *mess);
bool input_bool(char *mess);
/*@}*/


/** @name Command Line Parse Utilities
   The following functions are very useful from extracting user options
   from the command-line argc and argv.

   The first simply tells you whether a string appeared on the command
   line existed, and if so what its index in argv[] is. 

    index_of_arg(string,argc,argv) returns -1 if arg doesn't exist, returns
    index otherwise.

    string_from_args(key,argc,argv,default_string)
      searches for key on the command line. If it finds it, and
      theres another string to the right of it on the command line,
      returns that string (a pointer to the string...allocates no memory).
      If it doesn't find it, returns the default.

     {int,double,bool}_from_args are similar.

    EXAMPLE:
<pre>
     {
       int size = int_from_args("size",argc,argv,20);
       double weight = double_from_args("weight",argc,argv,12.4);
     }
</pre>
     if the program was ran with

       prog weight 404

     then after that code segment, size would be 20 and weight would be 404.0

  Details:
     These functions permit a leading dash in front of an argument, e.g.

       prog -weight 404

     would have had the same effect.

  
     If the string arghelp appears on the command line, then all these functions
     tell the user what key they expect, what type they are looking for, and
     what the default value is.
*/ 
/*@{*/    
bool arghelp_wanted(int argc, char *argv[]);
int index_of_arg(const char *opt_string, int argc, char *argv[]);

char *mk_string_from_args(const char *key,int argc,char *argv[],
			  char *default_value);

char *string_from_args(const char *key,int argc,char *argv[],
		       char *default_value);

char *string_from_args_insist(char *key,int argc,char *argv[]);

double double_from_args(const char *key, int argc, char *argv[], double default_value);
int int_from_args(const char *key, int argc, char *argv[], int default_value);
bool bool_from_args(const char *key, int argc, char *argv[], bool default_value);

/*@}*/    

/**
   rounds a double and turns it into an int.
   irint() the library function doesn't exists everywhere, so use this instead
*/
int my_irint(double x);

/** @name Printing of Basic Types
   The following functions print m1 then the argumnet and m2 to the requested
   stream. What use are they? They can be useful in semi-automated contruction
   of fprintf functions for big user-defined structures. 
*/
/*@{*/
void fprintf_int(FILE *s,char *m1,int x,char *m2);
void fprintf_realnum(FILE *s,char *m1,double x,char *m2);
void fprintf_float(FILE *s,char *m1,double x,char *m2);
void fprintf_double(FILE *s,char *m1,double x,char *m2);
void fprintf_bool(FILE *s,char *m1,bool x,char *m2);
void fprintf_string(FILE *s,char *m1,char *x,char *m2);
/*@}*/

/** @name Trivial string functions
*/
/*@{*/
bool char_is_in(char *s,char c);
int num_of_char_in_string(const char *s, char c);
int index_of_char(char *s,char c);
/*@}*/

FILE *am_fopen(char *filename,char *mode);



bool is_all_digits(const char *string);
bool is_a_number(const char *string);
bool bool_from_string(char *s,bool *r_ok);
int int_from_string(char *s,bool *r_ok);
double double_from_string(char *s,bool *r_ok);

extern void sensible_limits(
    double xlo,
    double xhi,
    double *res_lo,
    double *res_hi,
    double *res_delta
  );

int next_highest_power_of_two(int n);

double roundest_number_between(double lo,double hi);

/**convenience function for generating a valid tmp filename*/
char *mk_unique_tmpfile_path(const char *filename_prefix);

/** This function mirrors the glibc (GNU/Linux) counterpart (well, mirrors
   the API but not the security, etc).  See man 3 mktemp for details.
   NOTE: you *must* call this with a writeable character array. 
   Example pattern: char pattern[] = "/tmp/foo.XXXXXX" */
char *am_mktemp( char pattern[]);

/** get a valid path to use when making temp files. (windows or unix)*/
char *mk_valid_temp_dir_string(void);

void my_breakpoint(void);

/** Returns TRUE if and only if x is NaN or Inf (i.e. returns FALSE
   if and only if x is a completely legal number) */         
bool is_ill_defined(double x);

bool doubles_very_close(double a,double b);

/* The following segment wraps up the ability to new printf-style
   functions in a easy-to-use way.

   Provide you have #included the following two macros you can, for
   example do the following kind of definition... 

void rprint(report *rep,char *s)
{
  rprint_with_type(rep,s,WORDS_REPORT_TYPE);
}

void rprintf(report *rep,const char *format, ...)
{
  VA_MAGIC_INTO_BIG_ARRAY
  rprint(rep,big_array);
}

  And then if a programmer calls

  rprintf(rep,"%s = %g","pi",3.14159);

   it's as though the programmer had really called

  rprint(rep,"pi = 3.14159");

  This is fully portable C. 

  Note the dangerous buffer overflow possibility. Does anyone know how
  to avoid this problem? It would need an underlying function that could
  tell us how long the big_array is going to have to be so we could 
  allocate that amount of memory
*/

#define BIG_ARRAY_SIZE 50000

#define VA_MAGIC_INTO_BIG_ARRAY \
  static char big_array[BIG_ARRAY_SIZE]; \
  va_list args; \
  va_start(args, format); \
  if ( 0 > vsprintf(big_array,format,args) ) \
    my_error("vsprintf returned -ve number"); \
  va_end(args); 


/** Returns the time in seconds since some date in 1970.
   These days its actually implemented in h/utils/am_time.h, but 
   needs to be declared here to prevent vast swathes of warnings */
int global_time(void);

bool is_power_of_two(int x);

/** this is needed by both fprintf_dyv and fprintf_dym, which
   are now in separate files amdyv.c and amdym.c
   */
typedef struct buftab_struct
{
  int rows,cols;
  char ***strings;
} buftab;

char *bufstr(buftab *bt,int i,int j);
void fprint_buftab( FILE *s, buftab *bt);
void init_buftab( buftab *bt, int rows, int cols);
void free_buftab_contents(buftab *bt);
void set_buftab( buftab *bt, int i, int j, const char *str);

/**-----------------------------------------------------------------*/


#endif /** AMBS_H */
