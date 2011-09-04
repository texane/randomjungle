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
   File:        pfload.c
   Author:      Andrew W. Moore
   Created:     Mon Feb 11 10:50:18 EST 2002
   Description: Load in a sparse database

   Copyright 2002, The Auton Lab, CMU
*/

#include "pfload.h"

int facname_to_factor(char *facname,int line,int i)
{
  int len = (int)strlen(facname);
  int factor = -77;

  /*
    Removing :1 check 
    if ( len < 3 || facname[len-2] != ':' || facname[len-1] != '1' ) {
    my_errorf("factors should all be of the form <number>:1 but factor\n"
    "number %d on line %d is \"%s\"\n",
    i+1,line,facname);
    }
  */

  if (len == 0) {
    my_errorf( "facname_to_factor: empty facname on line %d pos %d, but how?",
	       line, i);
  }

  if (len > 2 && facname[len-2] == ':') {
    if (facname[len-1] != '1') {
      my_errorf("factor %d on line %d (%s) has : but is not followed by 1.",
		i+1, line, facname);
    }
    facname[len-2] = '\0';
  }

  if ( !is_a_number(facname) ) {
    my_errorf( "factors should all be of the form <number> or <number>:1 but\n"
	       " factor number %d on line %d has the following before the "
	       "\":1\":\n\"%s\"\n",
	      i+1,line,facname);
  }
  factor = atoi(facname);
  if ( factor < 0 )
    my_errorf("factor "
	      "number %d on line %d has value %d which is crazy\n",
	      i+1,line,factor);

  if (strlen(facname) != len) {
    if (strlen(facname) != len-2) {
      my_errorf( "facname_to_factor: internal error: len craziness on factor "
		 "%d, line %d, facname %s",
		 i+1, line, facname);
    }
    facname[len-2] = ':';
  }

  my_assert(factor >= 0);

  return factor;
}

ivec *mk_factors_from_string_array(string_array *sa,int line)
{
  int i, size;
  ivec *factors;

  size = string_array_size(sa);
  factors = mk_ivec( size);

  for ( i = 0 ; i < size ; i++ )
  {
    char *facname = string_array_ref(sa,i);
    int factor = facname_to_factor(facname,line,i);
    ivec_set(factors,i,factor);
  }

  /* We always want factors in increasing order. */
  ivec_sort( factors, factors);

  return factors;
}

precs *mk_precs_from_filename(char *fname)
{
  PFILE *s = safe_pfopen(fname,"r");
  bool finished = FALSE;
  int line = 0;
  precs *ps = mk_empty_precs();
  int report_size = 1;

  if (Verbosity >= 1) {
    printf("Beginning to load Records from file %s...\n",fname);
  }

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
	char *s0 = string_array_ref(sa,0);
        if ( !is_a_number(s0) )
	  my_errorf("First item on line %d of file %s is '%s', which is not a "
		    "number",line,fname,s0);
	else
	{
	  double activation = atof(s0);
	  ivec *factors;
	  prec *p;
          string_array_remove(sa,0);
	  factors = mk_factors_from_string_array(sa,line);
	  p = mk_prec(activation,factors);

	  add_pointer_to_precs(ps,p);
	  if ( precs_size(ps) >= report_size && Verbosity >= 1)
	  {
            if (Verbosity >= 4) {
              printf("Loaded %d line%s. Most recent line is:\n",
                     precs_size(ps),(precs_size(ps)==1)?"":"s");
              fprintf_prec(stdout,"p",precs_ref(ps,precs_size(ps)-1),"\n");
              report_size = 10 * report_size;
            }
	  }

	  free_ivec(factors);
	}
      }

      free_string_array(sa);
    }
  }

  pfclose(s);

  if (Verbosity >= 1) {
    printf("...Finished loading %d Records from file %s\n",
	   precs_size(ps),fname);
  }

  return ps;
}



int precs_num_rows(precs *ps)
{
  return precs_size(ps);
}

double precs_row_to_activation(precs *ps,int row)
{
  return prec_activation(precs_ref(ps,row));
}

ivec *precs_row_to_factors(precs *ps,int row)
{
  return prec_factors(precs_ref(ps,row));
}

int precs_num_factors(precs *ps)
{
  int result = 0;
  int i;
  ivec *factors;

  for ( i = 0 ; i < precs_size(ps) ; i++ ) {
    factors = precs_row_to_factors( ps, i);
    if (ivec_size( factors) == 0) continue;
    else result = int_max(result,1+sivec_max(factors));
  }

  return result;
}
