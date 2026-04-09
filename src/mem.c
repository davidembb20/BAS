/* Modified to be compatible with R memory allocation */
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include "bas.h"
#include <gsl/gsl_vector.h> // GSL-ADD 
#include <gsl/gsl_matrix.h> // GSL-ADD 

// Allocates a vector of nr doubles.
double *vecalloc(int nr)
{
  double *x;
  x=(double *) R_alloc(nr,sizeof(double));
  return x;
}

// Allocates a vector of nr integers.
int *ivecalloc(int nr)
{
  int *x;
  x= (int *) R_alloc( nr, sizeof(int));
  return x;
}

/* Converts an R integer / numeric vector to a gsl_vector (for levels)  */
gsl_vector * sexp_to_gsl_vector (SEXP v) /* >>> GSL-ADD */
{
    if (TYPEOF(v) != INTSXP && TYPEOF(v) != REALSXP)
        error("`levels' must be a numeric/integer vector");

    const int len = LENGTH(v);
    gsl_vector *g = gsl_vector_alloc(len);

    if (TYPEOF(v) == INTSXP)
        for (int i = 0; i < len; ++i)
            gsl_vector_set(g, i, (double)INTEGER(v)[i]);
    else
        for (int i = 0; i < len; ++i)
            gsl_vector_set(g, i, REAL(v)[i]);

    return g;
}



/* ---------- helpers: turn SEXPs into GSL objects ---------------------  */
/* Converts an R logical / integer matrix (column-major) to a gsl_matrix
   in row-major order.  (for positions)          */
gsl_matrix * sexp_to_gsl_matrix(SEXP m)
{
    if (TYPEOF(m) != REALSXP && TYPEOF(m) != INTSXP)
        error("'positions' and 'var.costs' must be integer or numeric matrices");

    if (!isMatrix(m))
        error("'m' must be a matrix");

    int *dim = INTEGER(getAttrib(m, R_DimSymbol));
    const int nrow = dim[0], ncol = dim[1];

    gsl_matrix *M = gsl_matrix_alloc(nrow, ncol);

    if (TYPEOF(m) == INTSXP) {
        int *src = INTEGER(m);  // R matrix is column-major
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i)
                gsl_matrix_set(M, i, j, (double)src[i + nrow * j]);
    } else {  // REALSXP
        double *src = REAL(m);
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i)
                gsl_matrix_set(M, i, j, src[i + nrow * j]);
    }

    return M;
}



unsigned char  **cmatalloc(int nr, int nc)
{
int k;
unsigned char  **x;
x = (unsigned char **) R_alloc(nr, sizeof(unsigned char *));
for (k=0;k<nr;k++){
  x[k] = (unsigned char  *) R_alloc(nc, sizeof(unsigned char));
  memset(x[k], 0, nc*sizeof(unsigned char));
}
return x;
}



/* getListElement from Writing R Extensions */


     SEXP getListElement(SEXP list, char *str)
     {
       SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
       int i;

       for (i = 0; i < length(list); i++)
         if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
           elmt = VECTOR_ELT(list, i);
           break;
         }
       return elmt;
     }
