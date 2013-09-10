//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * A GSL sparse matrix interface for the UMFPACK library.
 *
 */

#ifndef GSL_EXT_SPARSE_H
#define GSL_EXT_SPARSE_H

#include "umfpack.h"

#define SPMAT_TYPE_TRIPLET    0
#define SPMAT_TYPE_COL        1

typedef struct _gsl_sparse_matrix_complex {

    int n, m;
    int type;
    
    //
    // Triplet format arrays:
    //
    UF_long nelem, nz;

    UF_long *row;
    UF_long *col;
    double  *val_real;
    double  *val_imag;
    

    //
    // Column format arrays:
    //
    UF_long *Ap, *Ai;
    double *Ax, *Az;
    
    double Info [UMFPACK_INFO], Control[UMFPACK_CONTROL];
    
    //
    // Temporary
    //
    gsl_matrix_complex *M_dense;

} gsl_sparse_matrix_complex;

//int        gsl_sparse_matrix_complex_resize(gsl_sparse_matrix_complex *spmat);

gsl_sparse_matrix_complex *    gsl_sparse_matrix_complex_alloc(int n, int m, int nelem);
void                gsl_sparse_matrix_complex_free(gsl_sparse_matrix_complex *spmat);

gsl_complex            gsl_sparse_matrix_complex_get(gsl_sparse_matrix_complex *spmat, long i, long j);
int                    gsl_sparse_matrix_complex_set(gsl_sparse_matrix_complex *spmat, long i, long j, gsl_complex z);

gsl_matrix_complex *        gsl_sparse_matrix_complex_convert_to_dense(gsl_sparse_matrix_complex *spmat);
int                     gsl_sparse_matrix_complex_convert_col(gsl_sparse_matrix_complex *spmat);


int                gsl_sparse_matrix_complex_print_col(gsl_sparse_matrix_complex *spmat);
int                gsl_sparse_matrix_complex_print_triplet(gsl_sparse_matrix_complex *spmat);



int                gsl_sparse_matrix_complex_LU_solve(gsl_sparse_matrix_complex *spmat, double *b_real, double *b_imag, double *x_real, double *x_imag);


#endif
