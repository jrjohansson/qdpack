//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Description:
//
// A GSL based matrix backend for QDpack.
//

#ifndef _QDPACK_MATRIX
#define _QDPACK_MATRIX

#define GSL

#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <qdpack/gsl_ext.h>

//
// complex number abstraction layer: map to GSL functions
//

#define qdpack_complex gsl_complex
#define QDPACK_COMPLEX_ONE GSL_COMPLEX_ONE
#define QDPACK_COMPLEX_NEGONE GSL_COMPLEX_NEGONE
#define QDPACK_COMPLEX_ZERO GSL_COMPLEX_ZERO

#define QDPACK_REAL GSL_REAL
#define QDPACK_IMAG GSL_IMAG
#define QDPACK_SET_COMPLEX GSL_SET_COMPLEX

#define qdpack_complex_conjugate gsl_complex_conjugate

#define qdpack_complex_mul gsl_complex_mul
#define qdpack_complex_mul_real gsl_complex_mul_real
#define qdpack_complex_div gsl_complex_div

#define qdpack_complex_rect  gsl_complex_rect
#define qdpack_complex_polar gsl_complex_polar
#define qdpack_complex_exp   gsl_complex_exp

#define qdpack_complex_add   gsl_complex_add
#define qdpack_complex_sub   gsl_complex_sub

#define qdpack_complex_abs   gsl_complex_abs
#define qdpack_complex_abs2  gsl_complex_abs2
#define qdpack_complex_arg   gsl_complex_arg

#define qdpack_complex_negative gsl_complex_negative

#define QDPACK_TRANSPOSE_t CBLAS_TRANSPOSE_t
#define QDpackNoTrans      CblasNoTrans
#define QDpackTrans        CblasTrans
#define QDpackConjTrans    CblasConjTrans

//
// matrix representation
//
typedef struct _qdpack_matrix {

    int n, m;

    gsl_matrix_complex *data;

} qdpack_matrix_t;

// matrix
qdpack_matrix_t *qdpack_matrix_alloc(size_t m, size_t n);
void             qdpack_matrix_free(qdpack_matrix_t *mat);

qdpack_complex   qdpack_matrix_get(qdpack_matrix_t *mat, size_t i, size_t j);
void             qdpack_matrix_set(qdpack_matrix_t *mat, size_t i, size_t j, qdpack_complex value);
void             qdpack_matrix_set_all(qdpack_matrix_t *mat, qdpack_complex value);
void             qdpack_matrix_set_zero(qdpack_matrix_t *mat);
void             qdpack_matrix_set_identity(qdpack_matrix_t *mat);

void             qdpack_matrix_memcpy(qdpack_matrix_t *dst, qdpack_matrix_t *src);

void             qdpack_matrix_scale(qdpack_matrix_t *op, qdpack_complex z);
void             qdpack_matrix_add(qdpack_matrix_t *op1, qdpack_matrix_t *op2);
void             qdpack_matrix_sub(qdpack_matrix_t *op1, qdpack_matrix_t *op2);

void             qdpack_matrix_multiply(qdpack_matrix_t *A, qdpack_matrix_t *B, qdpack_matrix_t *C);

void             qdpack_matrix_blas_zgemm(QDPACK_TRANSPOSE_t transA, QDPACK_TRANSPOSE_t transB, 
                                          qdpack_complex alpha, qdpack_matrix_t *A, qdpack_matrix_t *B,
                                          qdpack_complex beta, qdpack_matrix_t *C);

qdpack_complex   qdpack_matrix_multiply_vmv(qdpack_matrix_t *mat, qdpack_matrix_t *v);
qdpack_complex   qdpack_matrix_dot(qdpack_matrix_t *a, qdpack_matrix_t *b); // multiply_vv

void             qdpack_matrix_exp(qdpack_matrix_t *m, qdpack_matrix_t *exp_m);

#define QDPACK_MATRIX_WRITE_FORMAT_REAL    0
#define QDPACK_MATRIX_WRITE_FORMAT_IMAG    1
#define QDPACK_MATRIX_WRITE_FORMAT_COMPLEX 2
#define QDPACK_MATRIX_WRITE_FORMAT_PYTHON  3
int              qdpack_matrix_write(FILE *f, qdpack_matrix_t *m, int format);
int              qdpack_matrix_python_write(FILE *f, qdpack_matrix_t *M);


// eigenstates
#define QDPACK_EIGEN_SORT_VAL_ASC   GSL_EIGEN_SORT_VAL_ASC
#define QDPACK_EIGEN_SORT_PHASE     GSL_EXT_EIGEN_SORT_PHASE

int qdpack_matrix_eigen_hermv(qdpack_matrix_t *m, qdpack_matrix_t *eval, qdpack_matrix_t *evec, int sort_order);
int qdpack_matrix_eigen_zgeev(qdpack_matrix_t *m, qdpack_matrix_t *eval, qdpack_matrix_t *evec, int sort_order);


#endif
