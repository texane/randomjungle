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
   File:        spardat.h
   Author:      Andrew W. Moore
   Created:     Fri Feb 15 14:06:41 EST 2002
   Description: Header for Sparse categorical dataset

   Copyright 2002, The Auton Lab, CMU
*/


#ifndef SPARDAT_H
#define SPARDAT_H

#include "amiv.h"
#include "amdyv.h"
#include "amdym.h"
#include "pfload.h"

/* A spardat is a representation of a machine learning dataset in which
   there is one output and many inputs.

   All attributes (inputs and output) are categorical and binary-valued.

   Most of the values in the dataset are zeroes. We want the memory size
   (and computation size) of things to be O(number of non-zero vales) and
   not O(num atts times num rows)
*/
typedef struct spardat
{
  string_array *attnum_to_name;  /* attnum_to_name[attnum] = name of attnum'th */
                                 /* input attribute. */

  ivec_array *attnum_to_rows;  /* attnum_to_rows[attnum] = sivec containing */
                               /* all rows in which attribute "attnum" takes */
                               /* a value of 1 */

  ivec_array *row_to_posatts; /* Reverse index:  */
                              /* row_to_posatts[row] = sivec containing all */
                              /* the attributes that take on a positive value */
                              /* in the given row */

  ivec *row_to_outval;        /* row_to_outval[row] = 1 if the row'th row is */
                              /*                        classed as active */
                              /* row_to_outval[row] = 0 if the row'th row is */
                              /*                        classed as inactive */

  ivec_array *outval_to_rows; /* An array of two sivecs. */
                              /* outval_to_rows[0] = sivec containing all */
                              /* rows with class == inactive */
                              /* outval_to_rows[1] = sivec containing all */
                              /* rows with class == active */
} spardat;

/* Glossary:

attnum: A small integer 0 <= attnum < num_atts denoting an attribute (a
"column" in the dataset)

row: A small integer 0 <= row < num_rows denoting an row (a
"record" in the dataset)

outval: A small integer that must be 0 or 1 denoting whether a row is
classed as active (1) or inactive (0)
*/

int spardat_num_atts(const spardat *sp);
int spardat_num_rows(const spardat *sp);

void fprintf_spardat(FILE *s,char *m1,spardat *x,char *m2);

/* As with ALL Auton function, you can tell this gives a pointer to
   something inside a data structure without returning anything that's
   a copy or that needs to be freed. You can tell this by the absence of
   "mk_" at the start of the function name. */
char *spardat_attnum_to_name(const spardat *sp,int attnum);

/* Returns a sivec of the rows in which the attribute has value "true" */
ivec *spardat_attnum_to_posrows(const spardat *sp,int attnum);

/* Returns the number of rows  in which the attribute has value "true" */
int spardat_attnum_to_num_posrows(const spardat *sp,int attnum);

/* Returns a sivec of the attributes that have the value "true" in the
   given row */
ivec *spardat_row_to_posatts(const spardat *sp,int row);

/* Returns the number of attributes that have the value "true" in the 
   given row */
int spardat_row_to_num_posatts(const spardat *sp,int row);

/* Returns the outval (0 or 1) associated with the given row */
int spardat_row_to_outval(const spardat *sp,int row);

int spardat_outval_to_num_rows(const spardat *sp,int outval);

/* Returns the set of rows in which the outval has the specified value */
ivec *spardat_outval_to_rows(const spardat *sp,int outval);

/* Returns the number of rows in "rows" in which attnum takes the
   value "True" */
int num_posrows_in_subset(const spardat *sp,ivec *rows,int attnum);

/* Returns the intersection of "rows" and the rows in which attnum
   is true. Note that the number of rows in the result equals
     num_posrows_in_subset(sp,rows,attnum)

   Remember that, following standard AUTON convention, this function MAKES
   something that will eventually need to be freed with free_ivec(...) */
ivec *mk_posrows_subset(const spardat *sp,ivec *rows,int attnum);

dyv *mk_spardat_times_dyv( const spardat *sp, const dyv *dv);
dyv *mk_spardat_transpose_times_dyv( const spardat *sp, const dyv *dv);

dym *mk_input_dym_from_spardat( const spardat *sp);
dyv *mk_output_dyv_from_spardat(const spardat *sp);

/* Makes a spardat. Following standard Auton conventions, you know that
   all arguments are COPIED in, so you can still access, alter or free
   the arguments afterwards and do so will not change the spardat data
   structure that was created. */
spardat *mk_spardat(string_array *attnum_to_name,
		    const ivec_array *row_to_posatts,
                    const ivec *row_to_outval);

/* Creates a copy of a spardat. */
spardat *mk_copy_spardat(const spardat *sp);

/* mk_spardat_from_filename: Makes a spardat from afc's strange file
     format.

   mk_filename_from_pfilename: Extracts the filename and parameters
     from a pfilename.  A pfilename is a filename "foo" plus a
     threshold and threshold direction: foo:5.7+ means file=foo,
     thresh=5.7, and output values >= 5.7 are "active" (1 in spardat
     outputs).

   Each Line of afc-format data corresponds to one row.

   Any row of the data-file begining with a # character is ignored.

   The first item on the line is the activation level of the row (a real
   value). The remining items are all of the form <attnum>:1

    Here's an example line:

     58.7 34:1 99:1 789:1

    That's a record in which three attnums (34 99 and 789) take positive
    values and the rest take 0 values. It has real-valud-activation 58.

   This function parses such a file and gets all the real-valued activation
   and then thresholds them.

   If high_mean_active TRUE then
     If real-valued-activation >= act_thresh the row is classed as active
     If real-valued-activation <  act_thresh the row is classed as inactive

   If high_mean_active is FALSE then
     If real-valued-activation >= act_thresh the row is classed as inactive
     If real-valued-activation <  act_thresh the row is classed as active
*/
char *mk_filename_from_pfilename(char *pfilename,double *r_act_thresh,
				 bool *r_high_means_active);
spardat *mk_spardat_from_filename(char *filename,double act_thresh,
				  bool high_means_active);

/* A pfilename is a single string that summarizes how to interpret the
   activation in a afc format datafile. For example

     test1.txt:50+ means use filename = test1.txt
                         use act_thresh = 50.0
                         use high_means_active == TRUE

     test1.txt:50- means use filename = test1.txt
                         use act_thresh = 50.0
                         use high_means_active == FALSE

  If pfilename ends in .csv then the spardat is made from a regular
  non-sparse cvs file, e.g. a-or-d.csv in this directory. In that case, you
  should specify "output attribute-name" on the command-line to specify
  which attribute is to be the output attribute that corresponds to "active"
*/
spardat *mk_spardat_from_pfilename(char *pfilename,int argc,char *argv[]);

precs *mk_precs_from_spardat( spardat *sd);
void write_spardat_to_filename( char *filename, spardat *sp);

/* Frees sp and all its subcomponents */
void free_spardat(spardat *sp);

/* Makes a subset of the original spardat in which only the rows
   mentioned in "rows" appear in the new spardat.

   Note the rows in the new spardat are renumbered to be contiguous between
   0 , 1 , ... spartdat_num_rows(new_spardat)-1

   Note spartdat_num_rows(new_spardat) == ivec_size(rows)

   Note that following standard Auton conventions, sp is unaffected by this */
spardat *mk_spardat_from_subset_of_rows(const spardat *sp,ivec *rows);

/* This is just like mk_spardat_from_subset_of_rows(), except it operates
   on columns. */
spardat *mk_spardat_from_subset_of_attnums( const spardat *sp, ivec *attnums);

/* Returns the number of active rows */
int spardat_num_active_rows(const spardat *sp);

/* Returns the number of rows mentioned in "rows" in which the
   output is active */
/* Takes time linear in rows and independent of number of active rows.
   This is more efficient than if we computed size of intersection
   of sp active rows and rows */
int spardat_num_active_in_rows(const spardat *sp,ivec *rows);

/* returns an array called posatts_array of size ivec_size(rows)

   posatts_array[q] = posatts from r'th row of sp, where r = rows[q]
*/
ivec_array *mk_posatts_array_from_spardat_rows(const spardat *sp,ivec *rows);



#endif /* #ifndef SPARDAT_H */
