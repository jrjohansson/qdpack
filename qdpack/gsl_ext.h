//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * A collection of functions that extend the GSL.
 */

#ifndef GSL_EXT_H
#define GSL_EXT_H

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>

#ifdef RCSS
#include "gsl_ext_rscc.h"
#endif

int gsl_ext_matrix_convert(gsl_matrix_complex *dst, gsl_matrix *src);
gsl_matrix_complex *gsl_ext_matrix_convert_and_free(gsl_matrix *m);

int gsl_ext_expm_complex(gsl_matrix_complex *A, gsl_matrix_complex *eA);

int gsl_ext_eigen_zgeev(gsl_matrix_complex *A_gsl, gsl_matrix_complex *evec, gsl_vector_complex *eval);
    
#define GSL_EXT_EIGEN_SORT_ABS   0
#define GSL_EXT_EIGEN_SORT_REAL  1
#define GSL_EXT_EIGEN_SORT_IMAG  2
#define GSL_EXT_EIGEN_SORT_PHASE 3
int gsl_ext_eigen_sort(gsl_matrix_complex *evec, gsl_vector_complex *eval, int sort_order);

#endif
