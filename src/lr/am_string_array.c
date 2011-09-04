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


/* 	File: am_string_array.c 
 * 	Author(s):
 * 	Date:
 * 	Purpose:
 */

#include "standard.h"
#include "ambs.h"
#include "amiv.h"
#include "am_string_array.h"
#include "amma.h"
#include "am_string.h"
#include "amdyv_array.h"


int String_Arrays_mallocked = 0;
int String_Arrays_freed = 0;


#define STRING_ARRAY_CODE 20542

string_array *mk_string_array(int size)
{
  string_array *result = AM_MALLOC(string_array);
  int i;

  result -> string_array_code = STRING_ARRAY_CODE;
  result -> size = size;
  result -> sarr_size = size;
  result -> sarr = AM_MALLOC_ARRAY(char_ptr,size);

  for ( i = 0 ; i < size ; i++ )
    result->sarr[i] = mk_copy_string("<UnDeFiNeD>");
  String_Arrays_mallocked += 1;
  return(result);
}

void free_string_array(string_array *sar)
{
  int i;
  sar -> string_array_code = 7777;

  for ( i = 0 ; i < sar->size ; i++ )
    if ( sar->sarr[i] != NULL )
      am_free_string(sar->sarr[i]);

  AM_FREE_ARRAY(sar->sarr,char_ptr,sar->sarr_size);
  AM_FREE(sar,string_array);

  String_Arrays_freed += 1;
}

string_array *mk_copy_string_array(const string_array *sa)
{
  string_array *nsa = mk_string_array(string_array_size(sa));
  int i;
  for ( i = 0 ; i < string_array_size(sa) ; i++ )
    string_array_set(nsa,i,string_array_ref(sa,i));
  return(nsa);
}


string_array *mk_string_array_x( int size, ...)
{
  /* XXX: no type checking */
  int i;
  char *s;
  va_list argptr;
  string_array *sa;

  sa = mk_string_array( size);

  va_start( argptr, size);
  for (i=0; i<size; ++i) {
    s = va_arg( argptr, char *);
    string_array_set( sa, i, s); /* Copies string into string array. */
  }
  va_end(argptr);

  return sa;
}


char **mk_array_from_string_array( string_array *sa)
{
  int size, i;
  char **array, *s;

  size = string_array_size( sa);
  array = AM_MALLOC_ARRAY( char *, size);
  for (i=0; i<size; ++i) {
    s = string_array_ref( sa, i);
    array[i] = mk_copy_string( s);
  }

  return array;
}

void fprintf_string_array(FILE *s,const char *m1,const string_array *sar,
                          const char *m2)
{
  int i;

  if (sar == NULL)
    fprintf(s,"%s = <NULL string array>%s",m1,m2);
  else
  {
    if ( string_array_size(sar) == 0 )
      fprintf(s,"%s = <empty string array>%s",m1,m2);
    else
    {
      for ( i = 0 ; i < string_array_size(sar) ; i++ )
	fprintf(s,"%s[%3d] = %s%s",
		m1,i,
		(string_array_ref(sar,i)==NULL) ? "NULL" : string_array_ref(sar,i),
		m2
		);
    }
  }
}

char *string_array_ref(const string_array *sar, int i)
{
  return(sar->sarr[i]);
}

void string_array_set(string_array *sar,int i,char *value)
/* value is COPIED in */
{
  if ( sar->sarr[i] != NULL )
    am_free_string(sar->sarr[i]);
  sar->sarr[i] = mk_copy_string(value);
}

int string_array_size(const string_array *sar)
{
  return(sar->size);
}

/* only for people who really want to save time allocating memory.
   after calling this, you should forget about the memory in string without
   freeing it.
 */
void add_to_string_array_no_copy(string_array *sa, char *string)
{
  if ( sa -> size == sa -> sarr_size )
  {
    int new_sarr_size = 2 + 2 * sa->sarr_size;
    char **new_sarr = AM_MALLOC_ARRAY(char_ptr,new_sarr_size);
    int i;

    for ( i = 0 ; i < sa->size ; i++ ) new_sarr[i] = sa->sarr[i];

    AM_FREE_ARRAY(sa->sarr,char_ptr,sa->sarr_size);
    sa -> sarr_size = new_sarr_size;
    sa -> sarr = new_sarr;
  }

  sa -> size += 1;
  sa -> sarr[sa->size-1] = string;
}

void add_to_string_array(string_array *sa, const char *string)
{
  if ( sa -> size == sa -> sarr_size )
    {
    int new_sarr_size = 2 + 2 * sa->sarr_size;
    char **new_sarr = AM_MALLOC_ARRAY(char_ptr,new_sarr_size);
    int i;

    for ( i = 0 ; i < sa->size ; i++ )
      new_sarr[i] = sa->sarr[i];

    AM_FREE_ARRAY(sa->sarr,char_ptr,sa->sarr_size);
    sa -> sarr_size = new_sarr_size;
    sa -> sarr = new_sarr;
    }

  sa -> size += 1;
  sa -> sarr[sa->size-1] = (string==NULL) ? NULL : mk_copy_string(string);
}

void string_array_add(string_array *sa,char *string)
{
  add_to_string_array(sa,string);
}

void string_array_remove(string_array* sa, int idx)
{
  int i;
  if ( string_array_size(sa) <= 0 ) my_error("string_array_remove: empty string_array");
  if ( idx < 0 || idx >= string_array_size(sa) )
    my_error("string_array_remove: bad index");

  /* WKW - New version shuffles the pointers over by one instead
     of using string_array_set, which copies the (i+1)th element into the
     ith spot and then frees the (i+1)th element.  This is faster and
     it avoids some memory bugs in the old version */
  if( sa->sarr[idx] != NULL )
    {
    am_free_string(sa->sarr[idx]);
    sa->sarr[idx]=NULL;
    }

  for( i = idx ; i < (sa->size - 1); i++ )
    sa->sarr[i] = sa->sarr[i+1];

  sa->sarr[sa->size-1] = NULL; /* Might help catch some weird errors */
  sa -> size -= 1;
}


string_array *mk_broken_string_using_seppers(const char *string, const char *seppers)
{
  string_array *sa = mk_string_array(0);
  int i = 0;
  while ( string != NULL && string[i] != '\0' )
    {
    int j = 0;
    while ( is_sepper(string[i],seppers) )
      i += 1;

    if ( string[i] != '\0' )
      {
      int part_size,backslashcount = 0;
      char *part,stopchar = ' ';
      int k;

      while ( string[i+j] != '\0' )
        {
	if(stopchar == ' ')
	  {
	  if(is_sepper(string[i+j],seppers))
            break;
	  else if((string[i+j]=='\"') && !(backslashcount%2))
            stopchar = '\"';
	  }
	else if(stopchar == '\"')
	  {
	  /* bug fix? this used to say stopchar = '\n' which made it put the
             whole rest of the line into one string once it had seen a double
             quote.  Now it only includes up to the next double quote. 8/24/99 JGS */
	  if((string[i+j] == '\"') && !(backslashcount %2))
            stopchar = ' ';
	  }

	if (string[i+j] == '\\') backslashcount++;
	else                     backslashcount=0;
	
        j++;
        }

      part_size = j+1;
      part = AM_MALLOC_ARRAY(char,part_size);

      for ( k = 0 ; k < j ; k++ )
        part[k] = string[i+k];
      if ( k != part_size-1 ) my_error("oaibxwibxpwibx");
      part[k] = '\0';
      string_array_add(sa,part);
      AM_FREE_ARRAY(part,char,part_size);
      }
    i = i+j;
    }
  return(sa);
}

string_array *mk_broken_string_using_seppers_only(char *string,char *seppers)
{
  string_array *sa = mk_string_array(0);
  char *p = strpbrk(string,seppers);
  char c;
  while(p){
    c = *p;
    *p = '\0';
    add_to_string_array(sa,string);
    *p = c;
    string = p+1;
    p = strpbrk(string,seppers);
  }
  add_to_string_array(sa,string);
  return sa;
}

string_array *mk_broken_string(const char *string)
{
  string_array *result = mk_broken_string_using_seppers(string,NULL);
  return(result);
}

string_array *mk_broken_quoteless_string(char *string)
{
  char *quoteless = mk_quoteless_string(string);
  string_array *result = mk_broken_string_using_seppers(quoteless,NULL);
  free_string(quoteless);
  return(result);
}

string_array *mk_string_array_subset(string_array *sa,ivec *indices)
{
  int size = ivec_size(indices);
  string_array *result = mk_string_array(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    string_array_set(result,i,string_array_ref(sa,ivec_ref(indices,i)));
  return result;
}


int find_index_in_string_array(const string_array *sa, const char *string)
{
  int result = -1;
  int i;
  for ( i = 0 ; i < string_array_size(sa) && result < 0 ; i++ )
    if ( eq_string(string,string_array_ref(sa,i)) )
      result = i;
  return(result);
}

bool string_array_member(string_array *sa,char *string)
{
  return(find_index_in_string_array(sa,string)>=0);
}



/************* NEW LINE PARSING CODE *************/

/* If line_format is WHITESPACE then the line is read SPACE_STYLE
   if line_format is COMMA      then the line is read COMMA_STYLE
   if line_format is BAR        then the line is read BAR_STYLE
   if lineformat is  ANY        then

        count the number of unquoted commas on the line (n_comma)
        and count the number of unquoted bars on the line (n_bar)

        if ( n_comma == 0 && n_bar == 0 ) use SPACE_FORMAT
        if ( n_comma >= n_bar )           use COMMA_FORMAT
        if ( n_bar > n_comma )            use BAR_FORMAT

   The line parser runs through a finite state machine. On
   each character it looks at the character type:

     S Space       - The character is the ' ' char
     C Comma       - The character is the ',' char
     A SingleQuote - The character is the '\'' char
     Q DoubleQuote - The character is the '\"' char
     T Token       - The character is something other than the above
     
   The line parser is building up an array of tokens. It begins with
   an empty array of tokens. It has a current token being built. It begins
   with the current token empty. After each character is read, it performs
   one of the following actions:

     ADD   Add the curent token to the array. Set the current token to empty
     PUSH  Put the current character at the end of the current token
     NIL   Do nothing
     DASH  Put a dash character at the end of the current token
     DP    Put a dash, then the current character at end of token
     UNKN  Add the UNKNOWN_STRING to the array. Clear current token


  COMMA_STYLE parsing:

       All whitespace to immediate left and right of commas is removed.
       All other contiguous blocks of whitespace are replaced with - symbols
         (outside quotes, N contiguous spaces are replaced with one -.
          inside quotes, N contiguous spaces are replaced with N -'s)
       The resulting tokens between commas are used.
       Empty string between commas is turned into UNKNOWN STRING

  BAR_STYLE parsing:
       Just like commas, except use bar (|) symbol instead of , symbol

  SPACE_STYLE parsing:

       All whitespace inside quotes are turned to dashes
       All other CONTIGUOUS blocks of whitespace are collapsed to one space
       Then the resulting tokens between whitespaces are used.
*/

/* Returns TRUE iff all characters in line_string are c (note special
   case: returns TRUE if string hads length zero) */
bool line_contains_only_character_c(char *line_string,char c)
{
  bool result = TRUE;
  int i;
  for ( i = 0 ; result && line_string[i] != '\0' ; i++ )
    result = line_string[i] == c;
  return result;
}

/* Returns TRUE iff all characters in line_string are '-' (note special
   case: returns TRUE if string hads length zero) */
bool line_contains_only_dashes(char *line_string)
{
  return line_contains_only_character_c(line_string,'-');
}

/* A line is interesting if its not all white space and
the leftmost non whitespace character isnt # */
bool line_string_is_interesting(char *line_string)
{
  int i;
  char first_non_whitespace = ' ';
  char second_non_whitespace = ' ';
  bool result;

  for ( i = 0 ; first_non_whitespace == ' ' && line_string[i] != '\0' ; i++ )
  {
    if ( line_string[i] != ' ' && line_string[i] != '\t' && 
         line_string[i] != '\r' )
      first_non_whitespace = line_string[i];
    if (first_non_whitespace != '\0') second_non_whitespace = line_string[i+1];
  }
  result = ( first_non_whitespace != ' ' );
  /* we allow the special sequence '##' to be a "machine readable comment */
  if ((first_non_whitespace == '#') && (second_non_whitespace != '#'))
    result = FALSE;

  if ( result && line_contains_only_dashes(line_string) )
    result = FALSE;

  return(result);
}

/* Searches the file for the next line that isn't all whitespace and
   that doesn't have # as its first non-whitespace character. 

   If no-such line before file-end, returns NULL */
char *mk_next_interesting_line_string(PFILE *s,int *line_number)
{
  char *line_string = NULL;
  bool finished = FALSE;

  while ( !finished )
  {
    line_string = mk_string_from_line(s);
    *line_number += 1;
    if ( line_string == NULL )
      finished = TRUE;
    else
      finished = line_string_is_interesting(line_string);

    if ( !finished && line_string != NULL )
    {
      free_string(line_string);
      line_string = NULL;
    }
  }

  return(line_string);
}

/* As above excepts breaks resulting line into a string array of tokens... */
string_array *mk_next_interesting_line(PFILE *s,int *line_number)
{
  char *str = mk_next_interesting_line_string(s,line_number);
  string_array *sa = (str == NULL) ? NULL : mk_broken_string(str);
  if ( str != NULL ) free_string(str);
  return(sa);
}

bool contains_a_number(string_array *sa)
{
  bool result = FALSE;
  int i;
  for ( i = 0 ; !result && i < string_array_size(sa) ; i++ )
    result = is_a_number(string_array_ref(sa,i));
  return(result);
}

#define UQ_START    0
#define UQ_MIDDLE   1
#define UQ_INSIDE_Q 2
#define UQ_INSIDE_A 3
#define UQ_STOP     4

void line_count_unquoted(char *string,
			 int *r_num_unquoted_commas,
			 int *r_num_unquoted_bars)
{
  int state = UQ_START;
  int i = 0;

  *r_num_unquoted_commas = 0;
  *r_num_unquoted_bars = 0;

  while ( state != UQ_STOP )
  {
    char c = string[i];
    if ( c == '\0' )
      state = UQ_STOP;
    else
    {
      switch ( state )
      {
        case UQ_START:
          if ( c == ' ' ) state = UQ_START;
          else if ( c == '\"' ) state = UQ_INSIDE_Q;
          else if ( c == '\'' ) state = UQ_INSIDE_A;
          else if ( c == ',' )
          {
	    *r_num_unquoted_commas += 1;
	    state = UQ_START;
	  }
          else if ( c == '|' )
          {
	    *r_num_unquoted_bars += 1;
	    state = UQ_START;
	  }
          else state = UQ_MIDDLE;
        break;
        case UQ_MIDDLE:
          if ( c == ' ' ) state = UQ_START;
          else if ( c == ',' )
          {
	    *r_num_unquoted_commas += 1;
	    state = UQ_START;
	  }
          else if ( c == '|' )
          {
	    *r_num_unquoted_bars += 1;
	    state = UQ_START;
	  }
          else state = UQ_MIDDLE;
        break;
        case UQ_INSIDE_A:
          if ( c == '\'' ) state = UQ_START;
          else state = UQ_INSIDE_A;
        break;
        case UQ_INSIDE_Q:
          if ( c == '\"' ) state = UQ_START;
          else state = UQ_INSIDE_Q;
        break;
        default: my_error("wiudbiuwb"); break;
      }
    }
    i += 1;
  }
}

string_array *mk_parse_data_line(char *string,int line_format)
{
  char separator;

  switch ( line_format )
  {
    case WHITESPACE_FORMAT: separator = ' '; break;
    case COMMA_FORMAT: separator = ','; break;
    case BAR_FORMAT: separator = '|'; break;
    case AUTO_FORMAT:
    {
      int num_unquoted_commas;
      int num_unquoted_bars;
      line_count_unquoted(string,&num_unquoted_commas,&num_unquoted_bars);
      if ( num_unquoted_commas == 0 && num_unquoted_bars == 0 )
	separator = ' ';
      else if ( num_unquoted_commas >= num_unquoted_bars )
	separator = ',';
      else
	separator = '|';

      break;
    }
    default:
      separator = '?';
      my_error("Unknown line_format");
      break;
  }

  return mk_string_array_from_parse(string,separator);
}

string_array *mk_next_tokens(PFILE *s,int *line_number,int line_format)
{
  char *line_string = mk_next_interesting_line_string(s,line_number);
  string_array *tokens;
  
  if ( line_string == NULL )
    tokens = NULL;
  else
    tokens = mk_parse_data_line(line_string,line_format);

  if ( line_string != NULL )
    free_string(line_string);

  return(tokens);
}

string_array *mk_default_attribute_names(int num_atts)
{
  string_array *sa = mk_string_array(num_atts);
  int i;
  for ( i = 0 ; i < num_atts ; i++ )
  {
    char buff[100];
    sprintf(buff,"x%d",i+1);
    string_array_set(sa,i,buff);
  }
  return(sa);
} 

/* private struct */
typedef struct parse_state
{
  int id;
  int s_action;  int s_next;
  int c_action;  int c_next;
  int a_action;  int a_next;
  int q_action;  int q_next;
  int t_action;  int t_next;
  int end_action;
} parse_state;

#define CGO  0
#define C1   1
#define C2   2
#define CQS  3
#define CQ   4
#define CAS  5
#define CA   6
#define SGO  7
#define S1   8
#define SQST 9
#define SQ   10
#define SAST 11
#define SA   12

#define ADD  0
#define PUSH 1
#define NIL  2
#define DASH 3
#define DP   4
#define UNKN 5
parse_state Parse_array[] =
/* State   S(act,next)  C(act,next)  A(act,next)  Q(act,next)  T(act,next) END(act)*/
{{ CGO ,   NIL ,CGO  ,  UNKN,CGO  ,  NIL ,CAS  ,  NIL ,CQS  ,  PUSH,C1   , UNKN},
 { C1  ,   NIL ,C2   ,  ADD ,CGO  ,  PUSH,C1   ,  PUSH,C1   ,  PUSH,C1   , ADD },
 { C2  ,   NIL ,C2   ,  ADD ,CGO  ,  DP  ,C1   ,  DP  ,C1   ,  DP  ,C1   , ADD },
 { CQS ,   DASH,CQ   ,  PUSH,CQ   ,  PUSH,CQ   ,  NIL ,CGO  ,  PUSH,CQ   , UNKN},
 { CQ  ,   DASH,CQ   ,  PUSH,CQ   ,  PUSH,CQ   ,  NIL ,C1   ,  PUSH,CQ   , ADD },
 { CAS ,   DASH,CA   ,  PUSH,CA   ,  NIL ,CGO  ,  PUSH,CA   ,  PUSH,CA   , UNKN},
 { CA  ,   DASH,CA   ,  PUSH,CA   ,  NIL ,C1   ,  PUSH,CA   ,  PUSH,CA   , ADD },
 { SGO ,   NIL ,SGO  ,  PUSH,S1   ,  NIL ,SAST ,  NIL ,SQST ,  PUSH,S1   , NIL },
 { S1  ,   ADD ,SGO  ,  PUSH,S1   ,  PUSH,S1   ,  PUSH,S1   ,  PUSH,S1   , ADD },
 { SQST,   DASH,SQ   ,  PUSH,SQ   ,  PUSH,SQ   ,  UNKN,SGO  ,  PUSH,SQ   , UNKN},
 { SQ  ,   DASH,SQ   ,  PUSH,SQ   ,  PUSH,SQ   ,  ADD ,SGO  ,  PUSH,SQ   , ADD },
 { SAST,   DASH,SA   ,  PUSH,SA   ,  UNKN,SA   ,  PUSH,SA   ,  PUSH,SA   , UNKN},
 { SA  ,   DASH,SA   ,  PUSH,SA   ,  ADD ,SA   ,  PUSH,SA   ,  PUSH,SA   , ADD }};

string_array *mk_string_array_from_parse(char *s,char separator)
{
  bool space_separated = separator == ' ';
  int state = (space_separated) ? SGO : CGO;
  bool finished = FALSE;
  int s_ptr = 0;
  string_array *tokarray = mk_string_array(0);
  int currtok_size = strlen(s) + 1;
  char *currtok = AM_MALLOC_ARRAY(char,currtok_size);
  int currtok_ptr = 0;

  while ( !finished )
  {
    parse_state *ps = &(Parse_array[state]);
    char c = s[s_ptr];
    int action;
    int next;

    if ( state != ps->id ) my_error("Parse_array misspecified");

    if ( c == '\0' )
    {
      finished = TRUE;
      next = -1;
      action = ps->end_action;
    }
    else if ( c == ' '  ) { action = ps->s_action ; next = ps->s_next; }
    else if ( c == '\t'  ) { action = ps->s_action ; next = ps->s_next; }
    else if ( c == separator ) { action = ps->c_action ; next = ps->c_next; }
    else if ( c == '\'' ) { action = ps->a_action ; next = ps->a_next; }
    else if ( c == '\"' ) { action = ps->q_action ; next = ps->q_next; }
    else                  { action = ps->t_action ; next = ps->t_next; }

    switch ( action )
    {
      case ADD :
        currtok[currtok_ptr] = '\0';
        add_to_string_array(tokarray,currtok);
        currtok_ptr = 0;
      break;
      case PUSH:
        currtok[currtok_ptr] = c;
        currtok_ptr += 1;
      break;
      case NIL :
        /* skip */
      break;
      case DASH:
        currtok[currtok_ptr] = '_';
        currtok_ptr += 1;
      break;
      case DP  :
        currtok[currtok_ptr] = '_';
        currtok_ptr += 1;
        currtok[currtok_ptr] = c;
        currtok_ptr += 1;
      break;
      case UNKN:
        add_to_string_array(tokarray,"?");
        currtok_ptr = 0;
      break;
      default: my_error("ljdnlkjs"); break;
    }

    state = next;
    s_ptr += 1;
  }

  AM_FREE_ARRAY(currtok,char,currtok_size);
  return(tokarray);
}

string_array *mk_string_array_from_argc_argv(int argc,char *argv[])
{
  string_array *sa = mk_string_array(argc);
  int i;
  for ( i = 0 ; i < argc ; i++ )
    string_array_set(sa,i,argv[i]);
  return(sa);
}

/* Makes and returns a string array of given size in which
   every entry contains a copy of the string stored in value */
string_array *mk_constant_string_array(int size,char *value)
{
  string_array *sa = mk_string_array(size);
  int i;
  for ( i = 0 ; i < size ; i++ )
    string_array_set(sa,i,value);
  return sa;
}

/* Parse a boolean from a string array. */
bool bool_from_argarray( char *key, string_array *argarray, bool defval)
{
  bool result;
  int fakeargc, i;
  char **fakeargv;

  fakeargc = string_array_size( argarray);
  fakeargv = mk_array_from_string_array( argarray);
  result = bool_from_args( key, fakeargc, fakeargv, defval);

  for( i=0; i<fakeargc; ++i) free_string( fakeargv[i]);
  AM_FREE_ARRAY( fakeargv, char *, fakeargc);

  return result;
}

/* Parse a int from a string array. */
int int_from_argarray( char *key, string_array *argarray, int defval)
{
  int result, fakeargc, i;
  char **fakeargv;

  fakeargc = string_array_size( argarray);
  fakeargv = mk_array_from_string_array( argarray);
  result = int_from_args( key, fakeargc, fakeargv, defval);

  for( i=0; i<fakeargc; ++i) free_string( fakeargv[i]);
  AM_FREE_ARRAY( fakeargv, char *, fakeargc);

  return result;
}

/* Parse a double from a string array. */
double double_from_argarray( char *key, string_array *argarray,
                             double defval)
{
  int fakeargc, i;
  double result;
  char **fakeargv;

  fakeargc = string_array_size( argarray);
  fakeargv = mk_array_from_string_array( argarray);
  result = double_from_args( key, fakeargc, fakeargv, defval);

  for( i=0; i<fakeargc; ++i) free_string( fakeargv[i]);
  AM_FREE_ARRAY( fakeargv, char *, fakeargc);

  return result;
}

/* Parse a string from a string array. */
char *mk_string_from_argarray( char *key, string_array *argarray,
			       char *defval)
{
  int size, idx;
  char *result;

  size = string_array_size( argarray);
  idx = find_index_in_string_array( argarray, key);

  if (idx == size-1) {
    my_errorf( "mk_string_from_argarray: key '%s' found at end of arguments\n",
	       "but we require a value after the key.",
	       key);
  }

  if (idx < 0) result = mk_copy_string( defval);
  else result = mk_copy_string( string_array_ref( argarray, idx+1));

  return result;
}

/* Split argument array into our options and extra options.  This
   allows one project to parse its options, then pass the unparsed
   options to another project's parse function. */
void mk_split_option_array( string_array *keys_with_args,
                            string_array *keys_without_args,
                            string_array *argarray,
                            string_array **my_opts,
                            string_array **extra_opts)
{
  int i, numargs;
  char *s;
  string_array *sa_mine, *sa_extra;
  
  numargs = string_array_size( argarray);

  /* Arrays to hold our args and extra args. */
  sa_mine = mk_string_array( 0);
  sa_extra = mk_string_array( 0);

  /* Sort args between ours (i.e. one of the keys listed above), and not
     ours.  Note that an extra options which takes an arg will be copied
     in-order to the sa_extra array, unless the arg val is the same as
     one of our keys. */
  for (i=0; i<numargs; ++i) {
    /* Token. */
    s = string_array_ref( argarray, i);

    if (string_array_member( keys_without_args, s)) {
      /* Token is lonely. */
      add_to_string_array( sa_mine, s);
    }

    else if (string_array_member( keys_with_args, s)) {
      /* Token should be followed by argument. */
      add_to_string_array( sa_mine, s);
      if (i+1 == numargs) {
        /* Print error if argument is missing. */
        my_errorf( "mk_active_train_option_arrays:\n"
                   "option '%s' is missing its required argument\n",
                   s);
      }
      s = string_array_ref( argarray, i+1);
      add_to_string_array( sa_mine, s);

      /* Increment i an extra time. */
      i = i + 1;
    }
    else {
      /* Extra arg. */
      add_to_string_array( sa_extra, s);
    }
  }

  /* Assign arg arrays to my_opts and extra_opts. */
  *my_opts = sa_mine;
  *extra_opts = sa_extra;

  return;
}



/*----------------------------------------------------------------*/

string_array *mk_split_string(const char *str, const char *delimiter)
{
    string_array *ret = mk_string_array(0);
    const char *begin = str;
    const char *c = str;

    if (!str || !delimiter || strlen(delimiter) == 0) return ret;
        
    while ((c = strstr(begin, delimiter)))
    {
	if (c - begin == 0) /*two delimiters in a row.  Add an empty string*/
	{
	    string_array_add(ret, "");
	}
	else /*add the sub string*/
	{
	    char *sub_str = AM_MALLOC_ARRAY(char, c-begin+1);
	    strncpy(sub_str, begin, c-begin);
	    sub_str[c-begin] = '\0';
	    
	    string_array_add(ret, sub_str);
	    AM_FREE_ARRAY(sub_str, char, c-begin+1);
	}
	c+= strlen(delimiter);
	begin = c;
    }
    if (!begin) /*str ended in a delimeter.  Add an empty string*/
    {
	string_array_add(ret, "");
    }
    else /*add the last string*/
    {
	char *tmp = mk_copy_string(begin);
	string_array_add(ret, tmp);
	free_string(tmp);
    }

    return ret;
}


char *mk_join_string_array(const string_array *sa, const char *delimiter)
{
    char *ret = mk_copy_string("");
    char *tmp;
    int i;
    int size = string_array_size(sa);
    for (i = 0; i < size; i++)
    {
	const char *delim = (i == 0) ? "" : delimiter;

	tmp = mk_printf("%s%s%s", ret, delim, string_array_ref(sa, i));
	free_string(ret);
	ret = tmp;
    }
    return ret;
}

string_array *mk_append_string_arrays(string_array *a,string_array *b)
{
  string_array *sa = mk_copy_string_array(a);
  int i;
  for ( i = 0 ; i < string_array_size(b) ; i++ )
    add_to_string_array(sa,string_array_ref(b,i));
  return sa;
}

