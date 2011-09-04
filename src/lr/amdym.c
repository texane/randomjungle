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


#include<errno.h>
#include "amma.h"
#include "ambs.h"
#include "amiv.h"
#include "amdyv.h"
#include "amdyv_array.h"
#include "am_string_array.h"
#include "amdym.h"


/*private header */
static int is_comment(const char* line);
static int find_number_of_columns(const char* line);
static int find_number_of_rows_and_columns(PFILE *fp , int *nrows, int *ncols);
static int assign_dym_values(PFILE* fp, dym* d);
static int assign_dym_values_for_csv(PFILE* fp, dym *factors, dyv *outputs);
static char *mk_readline(PFILE* fp );
static bool check_line_length(char *line, int lineno);

#define MAX_LINE_SIZE (10*1048576)
static char* delim = ",";



/** Allocate a matrix nrows X ncolumns */
dym* mk_dym(int nrows, int ncolumns)
{
  int i;
  dym *d;
  double **array;

  i = 0;
  d = (dym *) malloc(sizeof(dym));

  array = (double **) malloc(nrows * sizeof(double *));
  array[0] = (double *) malloc(nrows * ncolumns * sizeof(double));

  if( d && array && array[0] ) {
    for(i = 1; i < nrows; i++) array[i] = array[0] + (i * ncolumns);

    d->rows = nrows;
    d->cols = ncolumns;
    d->rows_allocated = nrows;
    d->tdarr = array;

  }
  else my_error( "mk_dym: Failed to allocate memory");

  return d;
}

dym *mk_copy_dym(const dym *d)
{
  return(mk_dym_scalar_mult(d,1.0));
}

void copy_dym(dym *d, dym *r_d)
{
  dym_scalar_mult(d,1.0,r_d);
}

void free_dym(dym* d)
{
  if(d != NULL) {
    if(d->tdarr != NULL) {
      if(d->tdarr[0]) free( d->tdarr[0]);
      free(d->tdarr);
    }
    free( d);
  }
  return;
}

int dym_rows(const dym *d)
{
  return(d->rows);
}

int dym_cols(const dym *d)
{
  return(d->cols);
}



/**
 *
 * 1. open the file to read
 * 2. Find the max rows and columns
 * 3. Allocate memory to hold the matrix
 * 4. Read the values from the file an assign them to the matrix
 *
 */
dym* mk_read_dym(const char* filename)
{
  int nrows =0 , ncols = 0;
  dym* d = NULL;
  PFILE* fp;

  fp = safe_pfopen(filename,"r");

  if( find_number_of_rows_and_columns(fp,&nrows,&ncols ) != -1 ) {
    d = mk_dym(nrows,ncols);
    if (d == NULL) {
      my_error( "mk_read_dym: unable to allocate memory for dym.");
    }

    prewind(fp); /* start from the beginning of the file */

    if( assign_dym_values(fp,d) == -1) {
      free_dym(d); /* error reading, free the allocated memory */
      my_error( "mk_read_dym: Unhandled error while reading.  Freeing dym.");
    }
  }

  pfclose(fp);
  return d;
}

/**
 *
 * 1. open the file to read
 * 2. Find the max rows and columns
 * 3. Allocate memory to hold the matrix
 * 4. Read the values from the file an assign them to the matrix
 *
 */
void mk_read_dym_for_csv(const char* filename, dym **factors, dyv **outputs)
{
  int nrows =0 , ncols = 0;
  PFILE* fp;

  if (factors == NULL) my_error("mk_read_dym_for_csv: Error: factors is NULL");
  if (outputs == NULL) my_error("mk_read_dym_for_csv: Error: outputs is NULL");

  fp = safe_pfopen(filename,"r");

  if( find_number_of_rows_and_columns(fp,&nrows,&ncols ) != -1 ) {
    *factors = mk_dym( nrows, ncols-1);
    if (*factors == NULL) {
      my_error( "mk_read_dym: unable to allocate memory for dym.");
    }
    *outputs = mk_dyv( nrows);

    prewind(fp); /* Jump back to the beginning of the file */

    if( assign_dym_values_for_csv( fp, *factors, *outputs) == -1) {
      free_dym( *factors); /* error reading, free the allocated memory */
      free_dyv( *outputs);
      my_error( "mk_read_dym: Unhandled error while reading.  Freeing dym.");
    }
  }

  pfclose(fp);
  return;
}

/** Returns 0 on success , -1 on failure */
int save_dym(PFILE* s, dym* d)
{
  int rc = 0 ;

  if( s && d ){
    int i;
    int j;

    for(i = 0 ; i < d->rows; i++){
      for(j =0 ; j < d->cols ; j++){
	rc = pfprintf(s,"%g",d->tdarr[i][j]);
	if( rc == -1) return rc;
	if( j == (d->cols - 1))/* last column */
	  rc = pfprintf(s,"\n");
	else
	  rc = pfprintf(s,delim);
	if( rc == -1) return rc;
      }
    }

  }
  pfflush(s);
  return rc;
}

void dym_scalar_mult(const dym *d, double alpha, dym *r_d)
{
  int i,j;

  for ( i = 0 ; i < r_d -> rows ; i++ )
    for ( j = 0 ; j < r_d -> cols ; j++ )
      r_d -> tdarr[i][j] = d->tdarr[i][j] * alpha;
}

dym *mk_dym_scalar_mult(const dym *d,double alpha)
{
  dym *result;

  result = mk_dym(d->rows,d->cols);
  dym_scalar_mult(d,alpha,result);
  return(result);
}

void dym_times_dyv(const dym *a, const dyv *b, dyv *result)
{
  int i;
  dyv *temp = mk_dyv(a->rows);
             /* We need a copy in case b and result share memory */

  if ( a->cols != b -> size )
    my_error("dym_times_dyv: sizes wrong");

  for ( i = 0 ; i < a->rows ; i++ )
  {
    double sum = 0.0;
    int j;
    for ( j = 0 ; j < a->cols ; j++ )
      sum += dym_ref( a, i, j) * dyv_ref( b, j);
    dyv_set( temp, i, sum);
  }

  copy_dyv(temp,result);
  free_dyv(temp);
}

dyv *mk_dym_times_dyv(const dym *a, const dyv *b)
{
  dyv *result = mk_dyv(a->rows);
  dym_times_dyv(a,b,result);
  return(result);
}

void dym_transpose( const dym *d, dym *r_d)
{
  dym *a = mk_dym(d->cols,d->rows);
             /* Note we have to first do the transpose to the result
                a, in case the routine was called with d's memory
                = r_d's memory */
  int i,j;

  for ( i = 0 ; i < d -> rows ; i++ )
    for ( j = 0 ; j < d -> cols ; j++ )
      a->tdarr[j][i] = d->tdarr[i][j];

  copy_dym(a,r_d);
  free_dym(a);
}

dym *mk_dym_transpose( const dym *a)
{
  dym *result = mk_dym(a->cols,a->rows);
  dym_transpose(a,result);
  return(result);
}

void dym_transpose_times_dyv(dym *p,dyv *v, dyv *r_dyv)
{
  int v_size, r_size, i, k;
  double factor, val;
  dyv *temp;

  v_size = dyv_size( v);
  r_size = dyv_size( r_dyv);

  temp = mk_zero_dyv( r_size);

  for (k=0; k<v_size; ++k) {
    factor = dyv_ref( v, k);
    for (i=0; i<r_size; ++i) {
      val = dym_ref( p, k, i);
      dyv_increment( temp, i, factor * val);
    }
  }

  copy_dyv(temp,r_dyv);
  free_dyv(temp);
  return;
}

dyv *mk_dym_transpose_times_dyv(dym *p,dyv *v)
/* Returns p^T v */
{
  dyv *result = mk_dyv(dym_cols(p));
  dym_transpose_times_dyv(p,v,result);
  return(result);
}

/*
  mk_dym_subset - builds a dym out of a subset of rows
  and columns from a matrix x such that if:
    r_given = rows(i)
    c_given = cols(j)
    val = x(r_given,c_given)
  then:
    new_dym(i,j) = val;
  Note: this can also be used to permute rows and columns.

  rows may be NULL meaning "use all rows"
  cols may be NULL meaning "use all columns"
*/
dym *mk_dym_subset( const dym *x, const ivec *rows, const ivec *cols)
{
  int Rn = (rows==NULL) ? dym_rows(x) : ivec_size(rows);
  int Cn = (cols==NULL) ? dym_cols(x) : ivec_size(cols);
  int r, c;
  dym* nu;

  nu = mk_dym(Rn,Cn);

  for(r=0; r<Rn; r++) {
    int row = (rows==NULL) ? r : ivec_ref(rows,r);
    for(c=0;c<Cn;c++) {
      int col = (cols==NULL) ? c : ivec_ref(cols,c);
      dym_set(nu,r,c,dym_ref(x,row,col));
    }
  }

  return nu;
}

dyv *mk_dyv_from_dym_column( const dym *dm, int col)
{
  int numrows, i;
  double val;
  dyv *dv;

  numrows = dym_rows( dm);
  dv = mk_dyv( numrows);

  for (i=0; i<numrows; ++i) {
    val = dym_ref( dm, i, col);
    dyv_set( dv, i, (int) val);
  }

  return dv;
}

dyv *mk_dyv_from_dym_row( const dym *dm, int row)
{
  int numcols;
  dyv *dv;
  numcols = dym_cols( dm);
  dv = mk_dyv_from_farr( dm->tdarr[row], numcols);
  return dv;
}

ivec *mk_ivec_from_binary_dym_column( const dym *dm, int col)
{
  int numrows, i;
  double val;
  ivec *iv;

  numrows = dym_rows( dm);
  iv = mk_ivec( numrows);

  for (i=0; i<numrows; ++i) {
    val = dym_ref( dm, i, col);

    if (val == 0.0) ivec_set( iv, i, 0);
    else if (val == 1.0) ivec_set( iv, i, 1);
    else {
      my_errorf( "mk_ivec_from_binary_dym_column: row=%d, col=%d, val=%f, "
                 "not binary.", i, col, val);
    }
  }

  return iv;
}

void fprintf_dym(FILE *s, const char *m1, const dym *d, const char *m2)
{
  if ( d == NULL ){
    fprintf(s,"%s = (dym *)NULL%s",m1,m2);
  }
  else if ( d->rows <= 0 || d->cols <= 0 )
    fprintf(s,"%s = <Dym with %d row%s and %d column%s>%s",
            m1,d->rows,(d->rows==-1)?"":"s",
            d->cols,(d->cols==-1)?"":"s",m2
           );
  else
  {
    int i;
    buftab bt;

    init_buftab(&bt,d->rows,d->cols + 4);

    for ( i = 0 ; i < d->rows ; i++ )
    {
      int j;
      set_buftab(&bt,i,0,(i == (d->rows-1)/2) ? m1 : "");
      set_buftab(&bt,i,1,(i == (d->rows-1)/2) ? "=" : "");
      set_buftab(&bt,i,2,"[");

      for ( j = 0 ; j < d -> cols ; j++ )
      {
        char buff[100];
        sprintf(buff," %g ",d->tdarr[i][j]);
        set_buftab(&bt,i,3+j,buff);
      }

      set_buftab(&bt,i,3+d->cols,"]");
    }

    fprint_buftab(s,&bt);
    free_buftab_contents(&bt);
  }
  fprintf(s,"\n");
}


 /** Read the values from the file and set them in dym 
 *  fp: File pointer at the beginning of the file
 *  d : allocated dym
 */
int assign_dym_values(PFILE* fp, dym* d)
{
  int rc, j, nrows, ncols, rowno, lineno;
  double v;
  char *line, *s, *endptr;
  string_array *sa;

  rc = 0;

  if (fp == NULL) my_error( "assign_dym_values: FILE *fp=NULL.");
  if (d == NULL) my_error( "assign_dym_values: dym *d=NULL.");

  nrows = dym_rows( d);
  ncols = dym_cols( d);
  lineno = 0;
  rowno = -1;

  while (!pfeof(fp)) {
    /* Get line. */
    line = mk_readline( fp);

    /* Check it is defined. */
    if (line==NULL) {
      if (pfeof(fp)) break; /* Read eof. */
      else my_error( "assign_dym_values: mk_readline() returned NULL.");
    }

    /* Line is happy. */
    lineno += 1;
    check_line_length(line, lineno);

    /* Skip comments. */
    if (is_comment(line)) {
      free( line);
      continue;
    }

    /* Line is not a comment. */
    rowno += 1;
    sa = mk_broken_string_using_seppers( line, (const char*) delim);
    if (string_array_size(sa) != ncols) {
      printf( "assign_dym_values: bad line:\n%s\n", line);
      my_errorf( "assign_dym_values: line %d has %d atts, but the data file\n"
                 "appears to have %d atts.",
                 lineno, string_array_size(sa), ncols);
    }

    /* Parse numbers from line. */
    for (j=0; j<ncols; ++j) {
      s = string_array_ref(sa,j);
      v = strtod( string_array_ref(sa,j), &endptr);
      if(errno == ERANGE) {
        printf( "assign_dym_values: bad line:\n%s\n", line);
        my_errorf( "assign_dym_values: value '%s' cause over- or underflow.\n"
                   "line=%d, att=%d\n", s, lineno, j+1);
      }
      if (endptr == s) {
        printf( "assign_dym_values: bad line:\n%s\n", line);
        my_errorf( "assign_dym_values: failed to convert '%s' to a number\n"
                   "line=%d, att=%d\n", s, lineno, j+1);
      }

      /* Set dym value. */
      v =d->tdarr[rowno][j] = v;
    }

    /* Done with line and tokens. */
    free_string_array( sa);
    free( line);
  }

  if (rowno+1 != nrows) {
    my_errorf( "assign_dym_values: read %d lines of data, but expected %d from"
               " first pass.\n",
               rowno+1, nrows);
  }

  return rc;
}

 /** Read the values from the file and set them in dym 
 *  fp: File pointer at the beginning of the file
 *  factors: allocated dym
 *  outputs: allocated dyv
 */
int assign_dym_values_for_csv( PFILE* fp, dym *factors, dyv *outputs)
{
  int rc, j, nrows, ncols, rowno, lineno;
  double v;
  char *line, *s, *endptr;
  string_array *sa;

  rc = 0;

  if (fp == NULL) my_error( "assign_dym_values_for_csv: FILE *fp=NULL.");
  if (factors == NULL) my_error( "assign_dym_values_for_csv: dym *d=NULL.");
  if (outputs == NULL) my_error( "assign_dym_values_for_csv: dym *d=NULL.");

  nrows = dym_rows( factors);
  ncols = dym_cols( factors);
  lineno = 0;
  rowno = -1;

  while (!pfeof(fp)) {
    /* Get line. */
    line = mk_readline( fp);

    /* Check it is defined. */
    if (line==NULL) {
      if (pfeof(fp)) break; /* Read eof. */
      else my_error( "assign_dym_values_for_csv: mk_readline() returned "
                     "NULL.");
    }

    /* Line is happy. */
    lineno += 1;
    check_line_length(line, lineno);

    /* Skip comments. */
    if (is_comment(line)) {
      free( line);
      continue;
    }

    /* Line is not a comment. */
    rowno += 1;
    sa = mk_broken_string_using_seppers( line, (const char*) delim);
    if (string_array_size(sa) != ncols+1) {
      printf( "assign_dym_values_for_csv: bad line:\n%s\n", line);
      my_errorf( "assign_dym_values_for_csv: line %d has %d values, but we\n"
                 "thought the data file had %d values per line.",
                 lineno, string_array_size(sa), ncols+1);
    }

    /* Parse numbers from line. */
    for (j=0; j < (ncols+1); ++j) {
      s = string_array_ref(sa,j);
      v = strtod( string_array_ref(sa,j), &endptr);
      if(errno == ERANGE) {
        printf( "assign_dym_values_for_csv: bad line:\n%s\n", line);
        my_errorf( "assign_dym_values_for_csv: value '%s' cause over- or "
                   "underflow.\n"
                   "line=%d, att=%d\n", s, lineno, j+1);
      }
      if (endptr == s) {
        printf( "assign_dym_values_for_csv: bad line:\n%s\n", line);
        my_errorf( "assign_dym_values_for_csv: failed to convert '%s' to a "
                   "number\n"
                   "line=%d, att=%d\n", s, lineno, j+1);
      }

      /* Set value. */
      if (j < ncols) dym_set( factors, rowno, j, v);
      else dyv_set( outputs, rowno, v);
    }

    /* Done with line and tokens. */
    free_string_array( sa);
    free( line);
  }

  if (rowno+1 != nrows) {
    my_errorf( "assign_dym_values_for_csv: read %d lines of data, but "
               "expected %d from first pass.\n",
               rowno+1, nrows);
  }

  return rc;
}

char *mk_readline(PFILE* fp )
{
   char *cptr, *in;

   in = (char *) malloc(MAX_LINE_SIZE * sizeof(char));

   if(in==NULL) {
     cptr = NULL;
     my_error( "mk_readline: Failed to allocate memory");
   }
   else {
     cptr = pfgets(in, MAX_LINE_SIZE, fp);
     if (cptr == NULL) free( in);
   }

   return cptr;
}

bool check_line_length(char *line, int lineno)
{
  if (line == NULL) {
    my_error("Internal error: check_line_length called for line==NULL.");
  }

  if (strlen(line) == MAX_LINE_SIZE && line[strlen(line)-1] != '\n') {
    fprintf(stderr, "Warning: check_line_length: line %d may be too long.\n",
            lineno);
    fprintf(stderr, "  If your file is okay, you may need to adjust\n");
    fprintf(stderr, "  utils/amdym.c:MAX_LINE_SIZE to be larger than\n");
    fprintf(stderr, "  its current size (%df KB)\n", MAX_LINE_SIZE);
    return 1;
  }
  return 0;
}

int find_number_of_columns(const char* line)
{
  int ncols = 0;

  if( line ){
    string_array *sa = mk_broken_string_using_seppers(line,(const char*)delim);
    ncols = string_array_size(sa);
    free_string_array(sa);
  }

  return ncols;

}

int find_number_of_rows_and_columns(PFILE* fp , int *nrows , int *ncols )
{
  int cols;
  char *line;
  *ncols = 0;
  *nrows = 0;

  while( ! pfeof(fp) ){
    line = mk_readline(fp);
    if(line) {
      check_line_length(line, 1+*nrows);

      if (!is_comment(line)) {
        cols = find_number_of_columns(line);
        if( *ncols < cols) *ncols = cols;
        ++(*nrows);
      }
      free(line);
    }
  }

  return 0;

}

int is_comment(const char* line)
{
  /* Check if first non-whitespace char is '#'. */
  size_t pos;

  /* First non-whitespace char with C99 function strspn. */
  pos = strspn( line, " \t");

  if (pos == strlen(line)-1) {
    /* First non-whitespace char is \n.  Treat as comment. */
    return 1;
  }
  else if (line[pos] == '#') {
    /* Comment. */
    return 1;
  }
  else {
    /* Not a comment. */
    return 0;
  }
}



/*************************************************************************/
/* TESTING                                                               */
/*************************************************************************/

#if 0
int assign_value(dym* d)
{

  if( d ){
    int i = 0;
    int j = 0;

    for(i =0 ; i < d->rows; i++){
      for(j=0 ; j < d->cols ; j++){
	d->tdarr[i][j] = 1.5;
      }
    }
  }

  return 0;
}

int display_dym(dym* d)
{

  if( d ){
    int i = 0;
    int j = 0;

    for(i =0 ; i < d->rows; i++){
      for(j=0 ; j < d->cols ; j++){
	printf("%g ",d->tdarr[i][j] );
      }
      printf("\n");
    }
  }

  return 0;
}



int display(PFILE *fp)
{
  int r, c, lineno;

  find_number_of_rows_and_columns(fp , &r ,&c);

  printf(" display %d %d\n",r,c);

  lineno = 0;
  while (! pfeof(fp) ){
    char* line = mk_readline(fp);
    if(line){
      lineno += 1;
      check_line_length(line, lineno);

      printf("line %s\n",line);
      printf("ncols %d\n",find_number_of_columns(line));
      free(line);
    }
  }

  return 0;
}

int main(int argc, char** argv)
{
  dym* d = mk_read_dym("dym.out");
  if( d ) {
    display_dym(d);
    free_dym(d);
  }

  return 0;

}

#endif
