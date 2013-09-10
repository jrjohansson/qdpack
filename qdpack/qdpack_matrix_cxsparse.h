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

#include <suitesparse/cs.h>

#include <complex.h>

//
// complex number abstraction layer: map to C99 complex type/functions
//

#define qdpack_complex              double complex
#define QDPACK_COMPLEX_ONE          (1.0+0.0I)
#define QDPACK_COMPLEX_NEGONE       (-1.0+0.0I)
#define QDPACK_COMPLEX_ZERO         (0.0+0.0I)

#define QDPACK_REAL(x)              creal(x)
#define QDPACK_IMAG(x)              cimag(x)
#define QDPACK_SET_COMPLEX(z,x,y)   {*z=x+y*1.0I;}

#define qdpack_complex_conjugate    conj

#define qdpack_complex_mul(a,b)         (a*b)
#define qdpack_complex_mul_real(a,b)    (a*b)
#define qdpack_complex_div(a,b)         (a/b)

#define qdpack_complex_rect(x,y)        (x+y*1.0I)
#define qdpack_complex_polar(r,theta)   (r*cexp(theta))
#define qdpack_complex_exp              cexp

#define qdpack_complex_add(a,b)         (a+b)
#define qdpack_complex_sub(a,b)         (a-b)

#define qdpack_complex_abs          cabs
#define qdpack_complex_abs2(x)      cpow(cabs(x),2)
#define qdpack_complex_arg          carg

#define qdpack_complex_negative(x)  (-x)

#define QDPACK_TRANSPOSE_t int
#define QDpackNoTrans      0
#define QDpackTrans        1
#define QDpackConjTrans    2

//
// matrix representation
//
typedef struct _qdpack_matrix {

    int n, m;

   // gsl_matrix_complex *data;

    cs_cl *data;

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

#define QDPACK_MATRIX_WRITE_FORMAT_REAL    0
#define QDPACK_MATRIX_WRITE_FORMAT_IMAG    1
#define QDPACK_MATRIX_WRITE_FORMAT_COMPLEX 2
#define QDPACK_MATRIX_WRITE_FORMAT_PYTHON  3
int              qdpack_matrix_write(FILE *f, qdpack_matrix_t *m, int format);
int              qdpack_matrix_sparse_write(FILE *f, qdpack_matrix_t *M, int format);
int              qdpack_matrix_python_write(FILE *f, qdpack_matrix_t *M);

void             qdpack_matrix_exp(qdpack_matrix_t *m, qdpack_matrix_t *exp_m);


qdpack_complex * qdpack_vector_dense(qdpack_matrix_t *v);

//
void qdpack_matrix_triplet(qdpack_matrix_t *mat);
void qdpack_matrix_compressed(qdpack_matrix_t *mat);

// eigenstates
#define QDPACK_EIGEN_SORT_VAL_ASC   0
#define QDPACK_EIGEN_SORT_PHASE     1

int qdpack_matrix_eigen_hermv(qdpack_matrix_t *m, qdpack_matrix_t *eval, qdpack_matrix_t *evec, int sort_order);
int qdpack_matrix_eigen_zgeev(qdpack_matrix_t *m, qdpack_matrix_t *eval, qdpack_matrix_t *evec, int sort_order);


#endif
