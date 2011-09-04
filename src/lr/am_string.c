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



/* am_string.c
   File:         amstr.c
   Author:       Andrew Moore
   Created:      June 25 1997 from amma.c by Frank Dellaert
   Description:  String related functions.

   Copyright (C) Andrew W. Moore, 1992
*/

#include <ctype.h>

#include "am_string.h"
#include "am_string_array.h"
#include "amma.h"
#include "file.h"

/* Used to remember messages about warnings. */
bool Control_m_message_delivered = FALSE;
bool Lone_control_m_message_delivered = FALSE;

/* Proto for private functions */
char *mk_amstr_substring(char *s,int start,int end);
bool old_string_pattern_matches(char *pattern, char *string);


char *mk_copy_string(const char *s)
{
  char *newv;

  if ( s == NULL )
  {
    newv = NULL;
  }
  else
  {
    newv = AM_MALLOC_ARRAY(char,1+strlen(s));
    sprintf(newv,"%s",s);
  }

  return(newv);
}

/*copy first n chars*/
char *mk_copy_string_n(const char *s, int n)
{
    char *newv;

    if (s == NULL)
    {
	newv = NULL;
    }
    else
    {
	int len = strlen(s);
	int num = n;
	if (num > len) num = len;

	newv = AM_MALLOC_ARRAY(char, 1+num);
	strncpy(newv, s, num);
	newv[num] = '\0';
    }

    return(newv);
}


#define MK_PRINTF_BUFFER_SIZE 100000

char *mk_printf(char *fmt, ...)
{
  va_list ap;
  char buff[MK_PRINTF_BUFFER_SIZE];

  va_start(ap, fmt);
  vsprintf(buff, fmt, ap);
  va_end(ap);
  if (strlen(buff) >= MK_PRINTF_BUFFER_SIZE)
  {
    printf("string length is %ld\n", (unsigned long) strlen(buff));
    my_error("mk_printf: string length exceeded static buffer.  Impossible to find error will occur if this is allowed to happen\n");
  }
  return (mk_copy_string(buff));
}

char *make_copy_string(const char *s)
{
  return(mk_copy_string(s));
}

int find_char_in_string(char *s,char c){
  char *p = strchr(s,c);
  return (p)? (p-s):-1;
}

void free_string(char *s)
{
  if ( s == NULL )
    my_error("am_free_string NULL s");
  AM_FREE_ARRAY(s,char,(1 + strlen(s)));
}

void am_free_string(char *s)
{
  free_string(s);
}

char *mk_concat_strings(char *s1, char *s2)
{
  char *s = AM_MALLOC_ARRAY(char,1+strlen(s1)+strlen(s2));
  sprintf(s,"%s%s",s1,s2);
  return s;
}

/* string_has_suffix("plop.ps",".ps") == TRUE
   string_has_suffix("plop.ps","ps") == TRUE
   string_has_suffix("plops",".ps") == FALSE
   string_has_suffix("plops","ps") == TRUE
   string_has_suffix("plop.ps.hello",".ps") == FALSE */
bool string_has_suffix(char *string,char *suffix)
{
  int suffix_length = (int) strlen(suffix);
  int string_length = (int) strlen(string);
  int matches = string_length >= suffix_length;
  int i;
  for ( i = 0 ; matches && i < suffix_length ; i++ )
    matches = string[string_length - suffix_length + i] == suffix[i];
  return matches;
}

/* Returns -1 if c does not appear in s
   Else returns smallest i such that s[i] == c */
int find_index_of_char_in_string(char *s,char c)
{
  int result = -1;
  int i;
  for ( i = 0 ; result < 0 && s[i] != '\0' ; i++ )
  {
    if ( c == s[i] )
      result = i;
  }
  return result;
}

int find_last_index_of_char_in_string(char *s,char c)
{
  int i = strlen(s) - 1;
  while(i > -1)
    {
    if ( c == s[i] )
      return i;
    i--;
    }
  return i;
}

/* Returns -1 if no member of s appears in stops.
   Else returns smallest i such that s[i] appears in stops */
int find_index_in_string(char *s,char *stops)
{
  int result = -1;
  int i;
  for ( i = 0 ; result < 0 && s[i] != '\0' ; i++ )
  {
    if ( 0 <= find_index_of_char_in_string(stops,s[i]) )
      result = i;
  }
  return result;
}

/* Makes a new null-terminated string using character
   s[start] s[start+1] ... s[end-1]
   strlen(result) = end - start

   PRE: 0 <= start < end <= strlen(s) 
*/
char *mk_substring(char *s,int start,int end)
{
  char *result;
  int i;
  if ( start < 0 || start >= end || end > (int)strlen(s) )
    my_error("mk_substring bad start or end");

  result = AM_MALLOC_ARRAY(char,end-start+1); /* +1 for '\0' */

  for ( i = 0 ; i < end-start ; i++ )
    result[i] = s[start+i];
  result[end-start] = '\0';
  return result;
}

/* Just like mk_substring except returns an empty string ("")
   if start >= end, and makes start at least 0 and end at most
   strlen(s) before going into action. Always succeeds. */
char *mk_friendly_substring(char *s,int start,int end)
{
  char *result;

  start = int_max(0,start);
  end = int_min((int)strlen(s),end);

  if ( end <= start )
    result = mk_copy_string("");
  else
    result = mk_substring(s,start,end);

  return result;
}

char upcase_char(char c)
{ /* Why not just use toupper() ? */
  return (c >= 'a' && c <= 'z') ? (c - 'a' + 'A') : c;
}

char *mk_quoteless_string(char *s)
{
  return mk_string_without_character(s,'\"');
}

char *mk_string_without_character(char *string,char c)
{
  int string_num_chars = (int)strlen(string);
  int string_num_bytes = string_num_chars + 1;
  char *temp = AM_MALLOC_ARRAY(char,string_num_bytes);
  int si = 0;
  int ti = 0;
  char *result;

  for ( si = 0 ; si < string_num_chars ; si++ )
  {
    if ( string[si] != c )
    {
      temp[ti] = string[si];
      ti += 1;
    }
  }

  temp[ti] = '\0';
  result = mk_copy_string(temp);
  AM_FREE_ARRAY(temp,char,string_num_bytes);
  return result;
}

bool is_sepper(char c, const char *seppers)
{
  bool result = c == ' ' || c == '\n' || c == '\r' || c == '\t';
  if ( !result && seppers != NULL )
  {
    int i;
    for ( i = 0 ; !result && seppers[i] != '\0' ; i++ )
      if ( seppers[i] == c ) result = TRUE;
  }
  return(result);
}

void add_char_to_dstring_s_p(char **ds_ref, char c,
			     int *pos_ref, int *power_ref)
{
  char *new_ds;
  int new_pow = *power_ref;

  if (*pos_ref >= *power_ref)
  {
    while (*pos_ref >= new_pow)
      if (*power_ref == 0)
	new_pow = 32;
      else
	new_pow += new_pow;
    new_ds = AM_MALLOC_ARRAY(char, new_pow);
    if (*ds_ref != NULL)
    {
      memcpy(new_ds, *ds_ref, *power_ref);
      AM_FREE_ARRAY(*ds_ref, char, *power_ref);
    }
    *ds_ref = new_ds;
    *power_ref = new_pow;
  }
  (*ds_ref)[(*pos_ref)++] = c;
}

char *mk_string_from_line( PFILE *f)
{
  unsigned char c;
  int bufsize, pos;
  char *buf, *s;
  
  bufsize = 16 * 4096;
  buf = AM_MALLOC_ARRAY( char, bufsize);
  buf[4095] = '\0';

  /* Newlines:
     UNIX/C specification: \n
     DOS/Windows:          \r\n
     Mac OS9 and earlier:  \r
     Mac OSX:              \n
  */

  /* Copy file chars into buf until hit non-printable (np).  If np ==
     '\r' (CR), then we read one more char (omc).  If omc != '\n', we
     raise an "old mac?" error.  For any other np, including '\n' and
     'EOF', we stop reading as if it were '\n'. */

  /* Because there are many kinds of spaces, we will convert all of
     them to a single space.  Exception: \n will still be a line break. */

  /* We arbitrarily limit line size to some multiple of 4096 (less
     one, for the string terminator), rather than complicate this code
     with a dynamically-growing buffer. */

  pos = 0;  /* Always equal to length of string read so far. */
  while (1) {
    /* Get char, being careful about line termination. */
    c = pfgetc( f);
    if (c == '\r') {
      /* Read extra DOS/Windows newline. */
      c = pfgetc( f);
      if (c != '\n') {
        my_errorf( "mk_string_from_line: Found \\r followed by char %d='%c'.\n"
                   "We do not handle this case correctly.  If \\r is the\n"
                   "line terminator, please, convert the file to \\n or\n"
                   "\\r\\n line termination using a utility program, or\n"
                   "simply\n"
                   "\n"
                   "  cat fname | tr '\\r' '\\n' > newfname\n");
      }
    }

    /* Convert space-but-not-line-feed chars into single space. */
    if (c != '\n' && isspace(c)) c = ' ';

    /* Break anything non-printable char, now that crazy spaces and
       multi-char carriage returns are handled.  This will include
       EOF.
    */
    if (!isprint(c)) break;
    
    /* Store char and continue. */
    buf[pos] = c;
    pos += 1;
    if (pos == bufsize-1) {
      my_errorf( "mk_string_from_line: Line length exceeds %d, which is\n"
                 "the maximum supported length.  Edit this function to\n"
                 "increase this limit, and recompile.\n", pos+1);
    }
  }

  /* pos <= bufsize-1, and buf[pos] has not been written to in our loop. */

  /* Copy buffer into properly-sized string.  To retain some compatibility
     with the old mk_string_from_line, we return NULL on EOF. */
  if (pfeof(f)) s = NULL;
  else {
    buf[pos] = '\0';
    s = mk_copy_string( buf);
    AM_FREE_ARRAY( buf, char, bufsize);
  }

  return s;
}

void add_to_error_message(char **errmess, char *new_mess)
{
  if (!new_mess) return;
  else
  {
    if (*errmess)
    {
      int len = strlen(*errmess) + strlen(new_mess) + 2;
      char *buf = AM_MALLOC_ARRAY(char,len);
      sprintf(buf,"%s%s",*errmess,new_mess);
      am_free_string(*errmess);
      *errmess = mk_copy_string(buf);
      AM_FREE_ARRAY(buf,char,len);
    }
    else *errmess = mk_copy_string(new_mess);
  }
}

void prepend_error_message(char **errmess, char *new_mess)
{
  if (!new_mess) return;
  else
  {
    if (*errmess)
    {
      int len = strlen(*errmess) + strlen(new_mess) + 2;
      char *buf = AM_MALLOC_ARRAY(char,len);
      sprintf(buf,"%s%s",new_mess,*errmess);
      am_free_string(*errmess);
      *errmess = mk_copy_string(buf);
      AM_FREE_ARRAY(buf,char,len);
    }
    else *errmess = mk_copy_string(new_mess);
  }
}

int am_isspace(int c)
{ /* See http://bugzilla.redhat.com/bugzilla/show_bug.cgi?id=86465
     This is not an issue if we don't reference that symbol, which means
     avoiding locale stuff like the libc isspace(). */
  if(	(c == ' ')
||	(c == '\f')
||	(c == '\n')
||	(c == '\r')
||	(c == '\t')
||	(c == '\v') )
	return 1;
return 0;
}

char *mk_copy_string_trans_chars(char *string, char *chars, char replacement)
{
  int len = (int) strlen(string);
  int temp_size = len+1;
  char *temp = AM_MALLOC_ARRAY(char,temp_size);
  int i;
  int j = 0;
  char *result;

  for ( i = 0 ; i < len ; i++ )
    {
    char c = string[i];
    if ( !char_is_in(chars,c) )
      temp[j] = c;
    else
      temp[j] = replacement;
    j++;
    }

  temp[j] = '\0';

  /* Can't just return temp because it's probably the wrong length amount of
     memory */

  result = mk_copy_string(temp);

  AM_FREE_ARRAY(temp,char,temp_size);

  return result;
}

/* return a string that is the concatenation of s1 and s2.  
   If free_inputs is TRUE, frees s1 and s2 before returning.*/
char *mk_strcat(char *s1, char *s2, bool free_inputs)
{
    char *ret = NULL;
    if (s1 && s2)
	ret = mk_printf("%s%s", s1, s2);
    else if (s1)
	ret = mk_copy_string(s1);
    else if (s2)
	ret = mk_copy_string(s2);
	    
    if (free_inputs && s1) free_string(s1);
    if (free_inputs && s2) free_string(s2);
    return ret;
}

int string_num_bytes(char *s)
{
  int string_num_chars = (int)strlen(s);
  return string_num_chars + 1;
}
