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


/* 	File: am_string_array.h
 * 	Author(s): Andrew Moore
 * 	Date:
 * 	Purpose: This is the new home for most string_array functions.
 *      string_array is a dynamically resizable array of strings. Some
 *      of the functions and the declaration used to live in amiv.h.
 *
 */
#ifndef AM_STRING_ARRAY_H
#define AM_STRING_ARRAY_H

#include "standard.h"
#include "ambs.h"
#include "amiv.h"
#include "am_string.h"

typedef struct string_array_struct
{
  int string_array_code;
  int size;
  int sarr_size;
  char **sarr;
} string_array, *string_array_ptr;

string_array *mk_string_array(int size);
void free_string_array(string_array *sar);
string_array *mk_copy_string_array(const string_array *sa);

/* Warning: The compiler won't catch type errors here. Pass the values
   as integers, or behavior will be undefined (read: bad things will
   probably happen sometime) */
string_array *mk_string_array_x( int size, ...);

char **mk_array_from_string_array( string_array *sa);

/* Sort into lexicographic order */
void sort_string_array(string_array *sa);
string_array *mk_sort_string_array(string_array *sa);
string_array *mk_sort_string_array_remove_duplicates(string_array *sa);
string_array *mk_string_array_from_array(char **sarr,int size);

/*reverse string array order */
string_array *mk_reverse_string_array(string_array *sa);

void fprintf_string_array(FILE *s,
                          const char *m1,
                          const string_array *sar,
                          const char *m2);

void fprintf_string_array_contents(FILE *s,string_array *sar);
void fprintf_string_array_contents_on_one_line(FILE *s,string_array *sar);

int string_array_size(const string_array *sar);

bool string_array_equal(string_array *a,string_array *b);

char *string_array_ref(const string_array *sar, int i);
void string_array_set(string_array *sar,int i,char *value);

/* Extends the length of the string_array by one (efficiently, if
   amortized over time) and adds "string" as the size-1'th entry where
   size is the new size.

   COPIES in string. It's okay if string is NULL.  */
void add_to_string_array(string_array *sa,const char *string);

/* Only for people who really want to save time allocating memory. After 
   calling this, forget about the memory in string without freeing it.
*/
void add_to_string_array_no_copy(string_array *sa, char *string);

/* Following function adds an element at the pos location.  pos can be
   at the end of the array */
void insert_in_string_array(string_array *sa, int pos, char *string);

/* XXX Depreciated synonym for add_to_string_array */
void string_array_add(string_array *sa,char *string);

/* Returns a result of length equal to size of indices.
   result[i] = sa[indices[i]] */
string_array *mk_string_array_subset(string_array *sa,ivec *indices);

/* split/join -- would be cool to support regex at some point.
   note - delimiter is a string, not a set of characters
*/
string_array *mk_split_string(const char *str, const char *delimiter);
char *mk_join_string_array(const string_array *sa, const char *delimiter);


/*
Tokenizes a string_array given a separator:
For Example with the ";" token the string_array on the left becomes
that on the right:
[  0] = something  
[  1] = 1                
[  2] = 2                [  0] something 1 2
[  3] = ;            =>  [  1] 3 4
[  4] = 3                [  2] 5
[  5] = 4                 
[  6] = ;
[  7] = 5

Note: It is not necessary for the last element of the original 
string_array to be a separator
*/

string_array *parse_string_array(string_array *sa, char *sep);
string_array *parse_string_array_endpoints(string_array *sa, char *sep, int start, int end);

/* Stop characters are any whitespace character, or any
   character in the seppers string.

   If seppers is NULL it's ignored.

   Any stop character is removed and used as a separator breaking string
   into smaller chunks.
   Each word (chunk) is placed as a separate entry into the string array.

   If string is NULL or "" or has only sepper characters, a string_array of
   length 0 is returned.

   Modified 1-23-96 JS: A single string may include whitespace including
   a newline character if delimited by single or double quotes.  Note the
   delimiting characters are still left in the string, as are all backslashes
   used to produce the characters '"\ in the string.
   Modified 1-31-96 JS: The single quote is no longer considered a delimiter
*/
string_array *mk_broken_string_using_seppers(const char *string, const char *seppers);

/* Here stop characters are only those in seppers (not whitespace!) */
string_array *mk_broken_string_using_seppers_only(char *string,char *seppers);

/* whitespace is removed, and each word is placed as a separate entry
   int the string array.
   If string is NULL or "" or has only whitespace, a string_array of
   length 0 is returned.  */
string_array *mk_broken_string(const char *string);

/*  As in mk_broken_string except removes all double quote marks before doing the splitting. */
string_array *mk_broken_quoteless_string(char *string);

char *string_before_col(string_array *seps, int col);
char *rightmost_string(string_array *seps);

/* If n == 0 returns "",
   else returns rightmost "n" entries in sa, concatenated with...

     'sep' char between them (if use_sep TRUE)
     nothing    between them (if use_sep FALSE)

   PRE: n <= size */
char *mk_string_from_last_n_with_separator(string_array *sa,int n,
					   bool use_sep,char sep);

/* In the functions below, entries are separated by whitespace. If n == 0 returns "",
   else returns rightmost "n" entries in sa, concatenated with ' '
   between them.

   PRE: n <= size */
char *mk_string_from_last_n(string_array *sa,int n);

/* If n == 0 returns "",
   else returns all the entries in sa, concatenated with sep char between
   them. */
char *mk_string_from_string_array_with_separator(string_array *sa,char sep);

/* In below, entries are separated by whitespace. If sa has
   no entries returns "", else all entries in sa, concatenated with ' '
   between them. */
char *mk_string_from_string_array(string_array *sa);

/* In below, entries are separated by nothing. If sa has
   no entries returns "", else all entries in sa, concatenated directly,
   eg { "Andrew" , "Moore" } --> "AndrewMoore" */
char *mk_string_from_string_array_no_gaps(string_array *sa);


/* Reads to the end of line, or end of file, skipping white space
   and adding each new word into the string array. If nothing on
   the line, returns a string_array of size zero. If nothing
   on line and file ends immediately returns NULL */
string_array *mk_string_array_from_line(FILE *s);

void string_array_malloc_report(void);

/* Finds lowest index in sa with value string, or returns -1 if not found. */
int find_index_in_string_array(const string_array *sa, const char* string);
/* case-insensitive version of find_index_in_string_array */
int caseless_find_index_in_string_array(string_array *sa,char *string);
/* True if the value of string is a string in string_array */
bool string_array_member(string_array *sa,char *string);


/**** Removal functions on string_arrays ****/

/* Reduces the size of sa by one.
   string_array_ref(sa,index) disappears.
   Everything to the right of string_array_ref(sa,index) is copied one to the left.

   Formally: Let saold be the string_array value beore calling this function
             Let sanew be the string_array value after calling this function.

PRE: string_array_size(saold) > 0
     0 <= index < string_array_size(saold)

POST: string_array_size(sanew) = string_array_size(saold)-1
      for j = 0 , 1, 2 ... index-1  : 
         string_array_ref(sanew,j) == string_array_ref(saold,j)

      for j = index , index+1 , ... string_array_size(sanew)-1:
         string_array_ref(sanew,j) == string_array_ref(saold,j+1) */

void string_array_remove(string_array *sa,int idx); 
                                   /* Reduces size by 1, removes index'th 
                                      element. All elements to right of
                                      delete point copied one to left.
                                      See comments in amdm.c more details */

/* Shrinks sa by one element by removing the rightmost element. 
   Example:

  Before: sa == ( "I" "hate" "to" "see" "you" "leave" )
    string_array_remove_last_element(sa)
  Before: sa == ( "I" "hate" "to" "see" "you" ) */
void string_array_remove_last_element(string_array *sa); 

/* PRE: 0 <= start_index <= end_index <= size
   Returns a new string array of size end_index - start_index
   consisting of 
     { sa[start_index] , sa[start_index+1] ... sa[end_index-1] } */
string_array *mk_string_array_segment(string_array *sa,
                                      int start_index,int end_index);

string_array *mk_string_array_from_ivec(ivec *iv);

/* Returns a string array in which the i'th element is "<prefix><i>",
   eg mk_int_string_array("plop",3) would return

   { "plop0" , "plop1" , "plop2" } */
string_array *mk_int_string_array(char *prefix,int size);


/* *************************************************** */
/* Operations on Sorted String Arrays */
/* A sorted string array is an ordinary string array where
   the following property holds:
   FIXME: Insert description here

  You can tell if a string array is sorted by calling FIXME: do the obvious
  If you call sorted operations on a nonsorted string_array, undefined behavior will result */

/* SUGGESTION: If you're thinking of using sorted string arrays for
   sets of strings, let me recommend extra/namer.h instead -- ?? */

int find_index_in_sorted_string_array(string_array *sa,char *s);
/* Finds lowest index in sa with value >= string, between 0 and
   string_array_size(sa), inclusive.  NOTE that string_array_size(sa)
   is the index one past the end of the string_array! */
int find_closest_index_in_sorted_string_array(string_array *sa, char* string);
/* place s into string array whether or not a duplicate already exists */
void insert_in_sorted_string_array(string_array *sa,char *s);
/* place s into string array only if it does not already exist there */
void maybe_insert_in_sorted_string_array(string_array *sa,char *s);


/* If n == 0 returns "",
   else returns all the entries in sa, concatenated with...

     'sep' char between them (if use_sep TRUE)
     nothing    between them (if use_sep FALSE)

   PRE: n <= size */
char *mk_string_from_string_array_basic(string_array *sa,bool use_sep,char sep);

void check_string_array_access(const string_array *sar, int i, const char *name);
void assert_string_array_shape(string_array *sar,int size,char *name);
void swap_string_array(string_array *sa, int i, int j);
void qsort_string_array(string_array *sa, int left, int right);

/* Input is a string which may optionally have curly braces
   and commas in it. All square and curly braces and commas are turned into
   spaces and then a string_array is made from the remaining unique
   space-separated tokens */
string_array *mk_string_array_from_set_notation(char *setnot);

/* I wasn't able to find good documentation for this function. For now, read the
   source. Sorry. --Pat */
string_array *mk_string_array_from_parse(char *s,char separator);

/* Searches for key on the command line. If it doesn't find key returns
   a string_array of length zero.

   If it does find it then looks at the next token. If it's a single
   simple token, simply returns a string_array of length 1, containing
   that token.

   If this token begins and ends with "" or if it contains
   whitespace, then breaks it up by removing all whitespace and commas
   and returns a string_array with the remaining elements

   If this token begins with [, searches for the next token (including
   the start one) that ends with ] and uses all the items inside,
   again separating with whitespace and commas.

   EXAMPLES:
    Each of the command-lines would cause a string array of size three
    to be created with elements "william" "chye" and "lee-moore" if the
    caller had set key be "names".

      ./program names [william chye lee-moore]
      ./program names [ william chye lee-moore ]
      ./program names [william,chye,lee-moore]
      ./program names william,chye,lee-moore
      ./program names "william,chye,lee-moore"
      ./program names "william chye lee-moore"

    But this command line would result in an array of length 1 with only
    "william" in the string array:

      ./program names william chye lee-moore
*/
string_array *mk_string_array_from_args(char *key, int argc, char** argv);
/* Same as above, but returns NULL if not present. */
string_array *mk_maybe_string_array_from_args(char *key, int argc, char** argv);

/* Doesn't do that kind of complicated interpolation above, just drops off each
   pair into the returned value. */
string_array *mk_string_array_from_argc_argv(int argc,char *argv[]);


/* Parse a boolean from a string array. */
bool bool_from_argarray( char *key, string_array *argarray, bool defval);

/* Parse an int from a string array. */
int int_from_argarray( char *key, string_array *argarray, int defval);

/* Parse a double from a string array. */
double double_from_argarray( char *key, string_array *argarray,
                             double defval);

/* Parse a string from a string array. */
char *mk_string_from_argarray( char *key, string_array *argarray,
			       char *defval);

/* Split argument array into our options and extra options.  This
   allows one project to parse its options, then pass the unparsed
   options to another project's parse function. */
void mk_split_option_array( string_array *keys_with_args,
                            string_array *keys_without_args,
                            string_array *argarray,
                            string_array **my_opts,
                            string_array **extra_opts);


/********* DATFILE PARSING UTILITIES ************/
/* moved in from amdmex.h */

/* int linestyle formats (we recommend AUTO_FORMAT)....

   The following constants determine how lines are read from
   a datafile. COMMA_FORMAT expects commas between each line on
   the datafile. WHITESPACE expects one or more spaces between 
   each item. BAR_FORMAT expects vertical bars (|)

   And AUTO_FORMAT will accept commas, whitespace and bars
   as separators. */
#define COMMA_FORMAT      0
#define WHITESPACE_FORMAT 1
#define BAR_FORMAT        2
#define AUTO_FORMAT       3

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
  
  SPACE_STYLE parsing:

       All whitespace inside quotes are turned to dashes
       All other CONTIGUOUS blocks of whitespace are collapsed to one space
       Then the resulting tokens between whitespaces are used.
*/


/* Return TRUE if and only if all items in sa can be parsed as numbers */
bool are_numbers(string_array *sa);

void fprint_string_array_csv(FILE *s,string_array *x);
  
/* A line is interesting if its not all white space and
the leftmost non whitespace character isnt # */
bool line_string_is_interesting(char *line_string);

bool line_contains_only_character_c(char *line_string,char c);
bool line_contains_only_dashes(char *line_string);

/* Searches the file for the next line that isn't all whitespace and
   that doesn't have # as its first non-whitespace character. 

   If no-such line before file-end, returns NULL */
char *mk_next_interesting_line_string(PFILE *s,int *line_number);

/* As above excepts breaks resulting line into a string array of tokens... */
/* DON'T USE THIS!! It's old and inferior to...

string_array *mk_next_tokens(FILE *s,int *line_number,int line_format);

which is documented later in this file */
string_array *mk_next_interesting_line(PFILE *s,int *line_number);
string_array *mk_parse_data_line(char *string,int line_format);
string_array *mk_next_tokens(PFILE *s,int *line_number,int line_format);
string_array *mk_default_attribute_names(int num_atts);
bool contains_a_number(string_array *sa);
/* not used or defined?
bool line_has_unquoted_comma(char *string); */
void line_count_unquoted(char *string, int *r_num_unquoted_commas, int *r_num_unquoted_bars);
bool all_numeric(string_array *sa);

/* ******************************************* */
/* string_matrix operations. Maybe belongs in its own file */

typedef struct string_matrix
{
  int array_size;
  int rows;
  int cols;
  string_array **sas;
} string_matrix, *string_matrix_ptr;

string_matrix *mk_string_matrix(int rows,int cols);
int string_matrix_rows(string_matrix *sm);
int string_matrix_cols(string_matrix *sm);
char *string_matrix_ref(string_matrix *sm,int row,int col);

/* Makes a COPY of char *value ... identical to string_matrix_set */
void string_matrix_set(string_matrix *sm,int row,int col,char *value);

/* Alias to string_matrix_set */
void sm_set(string_matrix *sm,int row,int col,char *value);

string_array *mk_string_array_from_string_matrix_col(string_matrix *sm,
						     int col);

void fprintf_string_matrix(FILE *s,char *m1,string_matrix *sm,char *m2);
void free_string_matrix(string_matrix *sm);
void string_matrix_add_row(string_matrix *sm);
int max_string_length_in_column(string_matrix *sm,int col);

/* Returns a new string array in which the i'th entry is
   the i'th row of a plain text representation of the contents
   of the string matrix. The columns of the table are separated by
   the strings in "seps". seps[0] always appears before the left-hand
   column. seps[1] before the second column. seps[cols] always
   appears to the right of the final column. The number of
   entries in "seps" must be one more than the number of columns
   in sm. 

   Alternatively, seps may be NULL in which case one space is
   printed between each column. No spaces to left or to right. */
	/* FIXME: That description can probably be written more clearly -- Pat */
string_array *mk_tabular_string_array(string_matrix *sm,string_array *seps);

/* As above but the i'th char in seps in the separator */
string_array *mk_tabular_string_array_simple(string_matrix *sm,
					     char *seps);

void render_string_matrix(FILE *s,char *comment,string_matrix *sm);

/* See comments for mk_tabular_string_array_simple regaring the meaning
   of "seps". */
void render_string_matrix_with_seps(FILE *s,string_matrix *sm,char *seps);

/* See comments for mk_tabular_string_array_simple regaring the meaning
   of "seps". 

   If the first string on a line is a single dash, then this version
   draws a complete line of dashes for the full width of the string
   matrix output. 
*/
void render_string_matrix_with_seps_and_dashes(FILE *s,string_matrix *sm,
					       char *seps);

void string_matrix_real_set(string_matrix *sm,int row,int col,double value);

/* breaks row_string. Error unless the number of space-sparated substrings
   equals number of cols in string matrix. Sets sm(row,i) to i'th substring
   forall i */
void string_matrix_row_from_broken_string(string_matrix *sm, int row, char *row_string);

/* Makes and returns a string array of given size in which
   every entry contains a copy of the string stored in value */
string_array *mk_constant_string_array(int size,char *value);

/* Makes an 'LS' style string matrix. That means it takes the
   string array and puts each element into cells of a string
   matrix. The string matrix has "cols" columns. The order in
   which string_array elements are placed is

      sa[0]      sa[r+0]    ....    sa[(cols-1)r+0]
      sa[1]      sa[r+1]    ....    sa[(cols-1)r+1]
        :                                  :
        :                                  :
        :                                  :
      sa[r-1]    sa[2r-1]   ....    sa[(cols-1)r-1]

    where r is the least r such that r*cols >= string_array_size

   ...and some of the rightmost column might be filled with
      empty cells.
*/
string_matrix *mk_ls_style_string_matrix_given_cols(string_array *name,int cols);

/* Returns the max string length in sa */
int string_array_max_length(string_array *sa);

/* Makes an 'LS' style string matrix cleverly designed so that when
   printed it uses less than "max_chars" characters per line. (it auto-chooses
   and sizes the columns) */
string_matrix *mk_ls_style_string_matrix(string_array *names,int max_chars);

int sm_rows(string_matrix *sm);

void sm_set_string(string_matrix *sm,int row,int col,char *string);

void sm_set_double(string_matrix *sm,int row,int col,double x);

void sm_set_int(string_matrix *sm,int row,int col,double n);

/***************** sosarray ***************/

bool string_less(char *s1,char *s2);
bool string_greater(char *s1,char *s2);
bool string_leq(char *s1,char *s2);
bool string_geq(char *s1,char *s2);

string_array *mk_string_array_from_string(char *s);

/* A sosarray is a regular old string_array, except it is treated as a set of
   integers.

   An string_array is a legal sosarray if it is sorted in increasing order with no
   duplicates.

   The following set of functions consititute a reasonable simple
   package of set-theory operations.

   Note that everything is as efficient as possible for a set package 
   except for adding a single element and deleting a single element, 
   which (because of our representation by means of sorted string_arrays) could
   take time linear in set size. */

bool is_sosarray(string_array *sa);

/* Returns number of elements in sosarray */
int sosarray_size(string_array *sosarr);

/* If sosarr has 0 elements returns 0
   If value > string_array_max(sosarr) returns size
   If value <= string_array_min(sosarr) returns 0
   Else returns index such that
      value <= string_array_ref(sosarr,index) 
      string_array_ref(sosarr,index-1) < value
      
   It returns the value such that string_array_insert(sa,index,value)
   would represent the set with value added to sa (assuming value
   wasn't already in sa). */
int find_sosarray_insert_index(string_array *sosarr,char *string);

/* Adds the element while maintaining legal sosarraykiness.
   (If element already there, no change)
   Time cost: O(size) */
void add_to_sosarray(string_array *sosarr,char *string);

/* Returns -1 if the string does not exist in sosarr.
   Else returns index such that
      string == string_array_ref(sosarr,string) 
  Time cost: O(log(size)) */
int index_in_sosarray(string_array *sosarr,char *string);

/* Returns true iff sosarr contains string
   Time cost: O(log(size)) */
bool is_in_sosarray(string_array *sosarr,char *string);

void sosarray_remove_at_index(string_array *sosarr,int idx);

/* Does nothing if string is not in sosarr.
   If string is in sosarr, the sosarray is updated to
   represent sosarr \ { string } */
void sosarray_remove_string(string_array *sosarr,char *string);
  
/* Returns answer to A subset-of B?
   Returns true if and only if the set of integers in a is
   a subset of the set of integers in b */
bool sosarray_subset(string_array *sosarra,string_array *sosarrb);

bool equal_string_array(string_array *sa1,string_array *sa2);

bool sosarray_equal(string_array *sosarra,string_array *sosarrb);

/* Returns TRUE iff A is a subset of B and A != B */
bool sosarray_strict_subset(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_union(string_array *sosarra,string_array *sosarrb);

/* Returns A \ B.
   This is { x : x in A and x not in B } */
string_array *mk_sosarray_difference(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_intersection(string_array *sosarra,string_array *sosarrb);

/* Returns TRUE iff A intersect B is empty. O(size) time */
bool sosarray_disjoint(string_array *sosarra,string_array *sosarrb);

string_array *mk_sosarray_from_string_array(string_array *sa);

/* Turns a space separated string into a sosarray.
   Example: "3 1 4 1 5 9" ====> { 1 , 3 , 4 , 5 , 9 } */
string_array *mk_sosarray_from_string(char *s);
char *sosarray_first(string_array *sosarr);
char *sosarray_last(string_array *sosarr);

string_array *mk_append_string_arrays(string_array *a,string_array *b);

#endif /* AM_STRING_ARRAY_H */
