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
   File:        spardat.c
   Author:      Andrew W. Moore
   Created:     Fri Feb 15 14:06:41 EST 2002
   Description: Sparse categorical dataset

   Copyright 2002, The Auton Lab, CMU
*/

#include "amiv.h"
#include "amdyv.h"

#include "spardat.h"

int spardat_num_atts(const spardat *sp)
{
  return string_array_size(sp->attnum_to_name);
}

int spardat_num_rows(const spardat *sp)
{
  return ivec_array_size(sp->row_to_posatts);
}

char *spardat_attnum_to_name(const spardat *sp,int attnum)
{
  return string_array_ref(sp->attnum_to_name,attnum);
}

ivec *spardat_attnum_to_posrows(const spardat *sp,int attnum)
{
  return ivec_array_ref(sp->attnum_to_rows,attnum);
}

int spardat_attnum_to_num_posrows(const spardat *sp,int attnum)
{
  return ivec_size(spardat_attnum_to_posrows(sp,attnum));
}

ivec *spardat_row_to_posatts(const spardat *sp,int row)
{
  return ivec_array_ref(sp->row_to_posatts,row);
}

int spardat_row_to_num_posatts(const spardat *sp,int row)
{
  return ivec_size(spardat_row_to_posatts(sp,row));
}

int spardat_row_to_outval(const spardat *sp,int row)
{
  return ivec_ref(sp->row_to_outval,row);
}

int spardat_outval_to_num_rows(const spardat *sp,int outval)
{
  return ivec_size(spardat_outval_to_rows(sp,outval));
}

ivec *spardat_outval_to_rows(const spardat *sp,int outval)
{
  my_assert(outval == 0 || outval == 1);
  return ivec_array_ref(sp->outval_to_rows,outval);
}

int num_posrows_in_subset(const spardat *sp,ivec *rows,int attnum)
{
  return size_of_sivec_intersection(spardat_attnum_to_posrows(sp,attnum),
				    rows);
}

ivec *mk_posrows_subset(const spardat *sp,ivec *rows,int attnum)
{
  return mk_sivec_intersection(spardat_attnum_to_posrows(sp,attnum),rows);
}

dyv *mk_spardat_times_dyv( const spardat *sp, const dyv *dv)
{
  int numrows, row;
  double sum;
  ivec *posatts;
  dyv *r_dv;

  numrows = spardat_num_rows( sp);
  r_dv = mk_dyv( numrows);

  for (row=0; row<numrows; ++row) {
    posatts = spardat_row_to_posatts( sp, row);
    sum = dyv_partial_sum( dv, posatts);
    dyv_set( r_dv, row, sum);
  }

  return r_dv;
}

dyv *mk_spardat_transpose_times_dyv( const spardat *sp, const dyv *dv)
{
  int numatts, att;
  double sum;
  ivec *posrows;
  dyv *r_dv;

  numatts = spardat_num_atts( sp);
  r_dv = mk_dyv( numatts);

  for (att=0; att<numatts; ++att) {
    posrows = spardat_attnum_to_posrows( sp, att);
    sum = dyv_partial_sum( dv, posrows);
    dyv_set( r_dv, att, sum);
  }

  return r_dv;
}

ivec_array *mk_row_to_posatts_from_precs(precs *ps)
{
  ivec_array *row_to_posatts = mk_empty_ivec_array();
  int row;
  for ( row = 0 ; row < precs_num_rows(ps) ; row++ )
  {
    ivec *posatts = precs_row_to_factors(ps,row);
    add_to_ivec_array(row_to_posatts,posatts);
  }
  return row_to_posatts;
}

ivec *mk_row_to_outval_from_precs(precs *ps,double act_thresh,
				  bool high_means_active)
{
  int num_rows = precs_num_rows(ps);
  ivec *row_to_outval = mk_ivec(num_rows);
  int row;
  for ( row = 0 ; row < num_rows ; row++ )
  {
    double activation = precs_row_to_activation(ps,row);
    bool active;
    int outval;

    if ( high_means_active )
      active = activation >= act_thresh;
    else
      active = activation <= act_thresh;

    outval = (active) ? 1 : 0;

    ivec_set(row_to_outval,row,outval);
  }

  return row_to_outval;
}

void fprintf_spardat(FILE *s,char *m1,spardat *x,char *m2)
{
  char *buff;

  buff = mk_printf("%s -> attnum_to_name",m1);
  fprintf_string_array(s,buff,x->attnum_to_name,m2);
  free_string(buff);

  buff = mk_printf("%s -> attnum_to_rows",m1);
  fprintf_ivec_array(s,buff,x->attnum_to_rows,m2);
  free_string(buff);

  buff = mk_printf("%s -> row_to_posatts",m1);
  fprintf_ivec_array(s,buff,x->row_to_posatts,m2);
  free_string(buff);

  buff = mk_printf("%s -> row_to_outval",m1);
  fprintf_ivec(s,buff,x->row_to_outval,m2);
  free_string(buff);

  buff = mk_printf("%s -> outval_to_rows",m1);
  fprintf_ivec_array(s,buff,x->outval_to_rows,m2);
  free_string(buff);
}

/* Returns the number of non-zero values in the sparse matrix of
   inputs */
int spardat_num_non_zero(spardat *sp)
{
  return sum_of_ivec_array_sizes(sp->row_to_posatts);
  }

spardat *mk_spardat(string_array *attnum_to_name,
		    const ivec_array *row_to_posatts,
                    const ivec *row_to_outval)
{
  spardat *sp = AM_MALLOC(spardat);
  int num_atts = string_array_size(attnum_to_name);
  int row;
  int num_rows = ivec_array_size(row_to_posatts);
  ivec *zero_length_ivec = mk_ivec(0);

  if ( ivec_size(row_to_outval) != ivec_array_size(row_to_posatts) )
  {
    my_errorf("The factor file contains information about %d compounds.\n"
	      "The activation file contains information about %d\n"
	      "activations. There should be the same number of records\n"
	      "in each.\n",ivec_size(row_to_outval),
	      ivec_array_size(row_to_posatts));
  }

  if (Verbosity >= 3) printf("Making spardat row_to_posatts...\n");
  sp->row_to_posatts = mk_copy_ivec_array(row_to_posatts);
  if (Verbosity >= 3) printf("Making spardat row_to_outval...\n");
  sp->row_to_outval = mk_copy_ivec(row_to_outval);
  if (Verbosity >= 3) printf("Making spardat attnum_to_rows...\n");
  sp->attnum_to_rows = mk_array_of_zero_length_ivecs(num_atts);

  for ( row = 0 ; row < num_rows ; row++ )
  {
    ivec *posatts = ivec_array_ref(row_to_posatts,row);
    int i;
    for ( i = 0 ; i < ivec_size(posatts) ; i++ )
    {
      int attnum = ivec_ref(posatts,i);
      ivec *rows_for_this_attnum;

      my_assert(ivec_array_size(sp->attnum_to_rows) > attnum);

      rows_for_this_attnum = ivec_array_ref(sp->attnum_to_rows,attnum);
      add_to_ivec(rows_for_this_attnum,row);
    }
  }

  if (Verbosity >= 3) printf("Making spardat attnum_to_names...\n");
  sp -> attnum_to_name = mk_copy_string_array(attnum_to_name);

  if (Verbosity >= 3) printf("Making spardat outval_to_rows...\n");
  sp -> outval_to_rows = mk_array_of_zero_length_ivecs(2);

  for ( row = 0 ; row < num_rows ; row++ )
  {
    int outval = ivec_ref(row_to_outval,row);
    ivec *rows_for_this_outval = ivec_array_ref(sp->outval_to_rows,outval);
    add_to_ivec(rows_for_this_outval,row);
  }

  free_ivec(zero_length_ivec);

  my_assert(string_array_size(attnum_to_name) ==
	    ivec_array_size(sp->attnum_to_rows));

  if (Verbosity >= 3) printf("...entire spardat constructed\n");

  if (Verbosity >= 3) {
    printf("spardat has %d rows, %d attributes, and %d non-zero "
	   "input values\n",
	   spardat_num_rows(sp),
	   spardat_num_atts(sp),
	   spardat_num_non_zero(sp));
  }

  return sp;
}

spardat *mk_copy_spardat(const spardat *sp)
{
  spardat *spcopy;
  spcopy = AM_MALLOC( spardat);
  spcopy->attnum_to_name = mk_copy_string_array( sp->attnum_to_name);
  spcopy->attnum_to_rows = mk_copy_ivec_array( sp->attnum_to_rows);
  spcopy->row_to_posatts = mk_copy_ivec_array( sp->row_to_posatts);
  spcopy->row_to_outval  = mk_copy_ivec( sp->row_to_outval);
  spcopy->outval_to_rows = mk_copy_ivec_array( sp->outval_to_rows);
  return spcopy;
}


spardat *mk_spardat_with_default_attnames(ivec_array *row_to_posatts,
					  ivec *row_to_outval)
{
  ivec *all_attnums = mk_identity_ivec(1+ivec_array_max_value(row_to_posatts));
  string_array *attnames = mk_string_array_from_ivec(all_attnums);
  spardat *sp = mk_spardat(attnames,row_to_posatts,row_to_outval);
  free_string_array(attnames);
  free_ivec(all_attnums);
  return sp;
}

spardat *mk_spardat_from_precs(precs *ps,double act_thresh,
			       bool high_means_active)
{
  ivec_array *row_to_posatts = mk_row_to_posatts_from_precs(ps);
  ivec *row_to_outval = mk_row_to_outval_from_precs(ps,act_thresh,
						    high_means_active);
  spardat *sp = mk_spardat_with_default_attnames(row_to_posatts,row_to_outval);
  free_ivec(row_to_outval);
  free_ivec_array(row_to_posatts);
  return sp;
}

/* Loads spardat using (first) afc pfile pformat (via precs). */
spardat *mk_spardat_from_filename(char *filename,double act_thresh,
				  bool high_means_active)
{
  precs *ps = mk_precs_from_filename(filename);
  spardat *sp = mk_spardat_from_precs(ps,act_thresh,high_means_active);
  free_precs(ps);
  return sp;
}

char *mk_filename_from_pfilename(char *pfilename,double *r_act_thresh,
				 bool *r_high_means_active)
{
  string_array *sa = mk_broken_string_using_seppers(pfilename,":");
  char *problem = NULL;
  char *filename = NULL;

  if ( string_array_size(sa) == 0) {
    problem = mk_copy_string( "The pfilename is an empty string.\n");
  }
  if ( string_array_size(sa) == 1 ) {
    problem=mk_copy_string("There's no colon in the middle of the pfilename");
  }
  else if ( string_array_size(sa) > 2 ) {
    problem=mk_copy_string("There's more than one colon");
  }

  if ( problem == NULL )
  {
    char *s1 = string_array_ref(sa,0);
    char *s2 = string_array_ref(sa,1);
    int s2len = (int)strlen(s2);
    if ( s2len < 2 )
      problem = mk_copy_string("Nothing useful after the colon");
    else
    {
      char c = s2[s2len-1];
      if ( c == '+' )
	*r_high_means_active = TRUE;
      else if ( c == '-' )
	*r_high_means_active = FALSE;
      else
	problem = mk_copy_string("pfilename should end in + or -");
      
      if ( problem == NULL )
      {
	s2[s2len-1] = '\0';
	if ( !is_a_number(s2) )
	  problem = mk_printf("%s is not a number",s2);
	else
	  *r_act_thresh = atof(s2);
	s2[s2len-1] = c;
      }
    }

    if ( problem == NULL )
      filename = mk_copy_string(s1);
  }
	
  if ( problem != NULL )
  {
    my_errorf("I was trying to parse \"%s\" as a pfilename, where a\n"
	      "pfilename should consist of a filename, followed by a\n"
	      "colon followed by a threshold follow by a + or -. For\n"
	      "example, \"data.txt:50+\" would have been good. But\n"
	      "with what you gave me there was the following problem:\n"
	      "%s\n",pfilename,problem);
  }

  free_string_array(sa);
  return filename;
}

ivec_array *mk_row_to_posatts_from_filename(char *fname)
{
  PFILE *s = safe_pfopen(fname,"r");
  bool finished = FALSE;
  int line = 0;
  ivec_array *row_to_posatts = mk_empty_ivec_array();
  int report_size = 1;

  printf("Beginning to load afc Factors from file %s...\n",fname);

  while ( !finished )
  {
    string_array *sa = mk_next_tokens(s,&line,AUTO_FORMAT);
    if ( sa == NULL )
      finished = TRUE;
    else
    {
      int size = string_array_size(sa);
      if ( size < 1 )
	my_errorf("On line %d of %s there are no tokens",line,fname);
      else
      {
	ivec *f = mk_factors_from_string_array(sa,line);
	add_to_ivec_array(row_to_posatts,f);

	if ( ivec_array_size(row_to_posatts) >= report_size )
	{
	  printf("Loaded %d line%s from %s\n",ivec_array_size(row_to_posatts),
		 (ivec_array_size(row_to_posatts)==1)?"":"s",fname);
	  report_size = 2 * report_size;
	}

	free_ivec(f);
      }

      free_string_array(sa);
    }
  }

  pfclose(s);

  printf("...Finished loading %d afc Factors from file %s\n",
	 ivec_array_size(row_to_posatts),fname);

  return row_to_posatts;
}

/*
We now have the ability to load spardats.

I suggest you do the following in your h/afc directory...

ln -s <revisit - where did these files go?>/jun02 jun02

(Maybe also do it in your Linux_2.4.17_x86_gcc.debug directory)

jun02 has the following files...

activations3.ssv  factors3.ssv  names3.ssv
activations4.ssv  factors4.ssv  names4.ssv

These correspond to two separate sets of high throughput runs called TEST3
and TEST4

The factors for test X (where X = 3 or 4) are in factorsX.ssv . The
only difference between this and the earlier spardat data files is that
factorsX.ssv does not contain the output activation level for the experiments
at the start of each line. It only contains the inputs (the posatts).

There are K sets of output activations for the X experiment, where
  K = 4 if X = 3
  K = 2 if X = 4

The j'th activation class for the i'th row of factorX is specified in
the i'th row of the Actj column in activationsX.ssv.

To load a spardat from this format, you must specify a wacky "jun02"
filename syntax of this form:

   <factorsfile>:<activationsfile>:<actcolumnname>

for example:

   jun02/factors3.ssv:jun02/activations3.ssv:Act1

The spardat (mk_spardat_from_pfilename) has been updated to cope
with this syntax as an additional option. It automatically figures
out when you are using this filename form as opposed to the other
spartdat loading options.
*/
spardat *mk_spardat_from_pfilename(char *pfilename,int argc,char *argv[])
{
  bool link;
  spardat *sp;
  string_array *sa;

  link = index_of_arg( "link", argc, argv) > 0;

  sp = NULL;
  sa = mk_broken_string_using_seppers(pfilename,":");

  if ( string_array_size(sa) == 2 ) {
    double act_thresh;
    bool high_means_active;
    
    char *filename = mk_filename_from_pfilename(pfilename,&act_thresh,
                                                &high_means_active);
    sp = mk_spardat_from_filename(filename,act_thresh,high_means_active);
    free_string(filename);
  }
  else {
    my_errorf( "mk_spardat_from_filename: filename='%s'\n"
               "****   needs exactly one colon in order to be loaded\n"
               "****   as a spardat", pfilename);
  }
  free_string_array(sa);

  if (Verbosity >= 1) {
    printf("spardat loaded from pfilename %s has...\n"
	   "  %d rows, %d attributes, and %d non-zero input values.\n"
	   "  It has %d positive output values.\n",
	   pfilename,spardat_num_rows(sp),
	   spardat_num_atts(sp),
	   spardat_num_non_zero(sp),spardat_outval_to_num_rows(sp,1));
  }

  return sp;
}

precs *mk_precs_from_spardat( spardat *sp)
{
  int numatts, numrows, i;
  double dval;
  ivec *factors;
  prec *p;
  precs *ps;

  numatts = spardat_num_atts( sp);
  numrows = spardat_num_rows( sp);

  /* Make precs. */
  ps = mk_empty_precs();
  for (i=0; i<numrows; ++i) {
    dval    = spardat_row_to_outval( sp, i);
    factors = spardat_row_to_posatts( sp, i);
    p = mk_prec( dval, factors);
    add_to_precs( ps, p);
    free_prec(p);
  }

  return ps;
}

void free_spardat(spardat *sp)
{
  free_string_array(sp->attnum_to_name);
  free_ivec_array(sp->attnum_to_rows);
  free_ivec_array(sp->row_to_posatts);
  free_ivec(sp->row_to_outval);
  free_ivec_array(sp->outval_to_rows);
  AM_FREE(sp,spardat);
}

spardat *mk_spardat_from_subset_of_rows(const spardat *sp,ivec *rows)
{
  ivec_array *sub_row_to_posatts=mk_ivec_array_subset(sp->row_to_posatts,rows);
  ivec *sub_row_to_outval = mk_ivec_subset(sp->row_to_outval,rows);
  spardat *sub = mk_spardat(sp->attnum_to_name,
			    sub_row_to_posatts,sub_row_to_outval);
  free_ivec_array(sub_row_to_posatts);
  free_ivec(sub_row_to_outval);
  return sub;
}  

spardat *mk_spardat_from_subset_of_attnums( const spardat *sp, ivec *attnums)
{
  int numrows, rowindex, attindex, att, newatt, numinv;
  ivec *attinv, *row, *newrow;
  string_array *sub_attnum_to_name;
  ivec_array *sub_row_to_posatts;
  spardat *sub;

  /*
    attinv is an ivec that maps an attnum to its position in attnums:
      attnums[ attinv[ attnum]] = attnum
  */
  if (ivec_size(attnums) == 0 ) attinv = mk_ivec(0);
  else attinv = mk_invert_nonneg_ivec( attnums);
  numinv = ivec_size( attinv);

  /* Create new attnames. */
  sub_attnum_to_name = mk_string_array_subset( sp->attnum_to_name, attnums);

  /* Build row array. */
  numrows = spardat_num_rows( sp);
  sub_row_to_posatts = mk_empty_ivec_array();
  for (rowindex=0; rowindex<numrows; ++rowindex) {
    row = spardat_row_to_posatts( sp, rowindex);
    newrow = mk_ivec(0);

    /* Renumber or drop atts of this row. */
    for (attindex=0; attindex < ivec_size( row); ++attindex) {
      att = ivec_ref( row, attindex);
      if (att >= numinv) continue;
      newatt = ivec_ref( attinv, att);
      if (newatt < 0) continue;
      add_to_ivec( newrow, newatt); /* Let's hope this isn't a bottleneck
                                       in practice. */
    }

    add_to_ivec_array( sub_row_to_posatts, newrow);
    free_ivec( newrow);
  }

  /* Make new spardat. */
  sub = mk_spardat( sub_attnum_to_name, sub_row_to_posatts, sp->row_to_outval);

  /* Clean up. */
  free_ivec( attinv);
  free_string_array( sub_attnum_to_name);
  free_ivec_array( sub_row_to_posatts);

  return sub;
}

int spardat_num_active_rows(const spardat *sp)
{
  return ivec_size(spardat_outval_to_rows(sp,1));
}

/* Returns the number of rows mentioned in "rows" in which the
   output is active */
/* Takes time linear in rows and independent of number of active rows.
   This is more efficient than if we computed size of intersection
   of sp active rows and rows */
int spardat_num_active_in_rows(const spardat *sp,ivec *rows)
{
  int sum = 0;
  int i;
  for ( i = 0 ; i < ivec_size(rows) ; i++ )
  {
    int row = ivec_ref(rows,i);
    if ( ivec_ref(sp->row_to_outval,row) == 1 )
      sum += 1;
  }
  return sum;
}
