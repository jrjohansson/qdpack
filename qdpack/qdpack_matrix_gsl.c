//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

// Description:
//
// A GSL based matrix backend for QDpack.
//
#include <stdio.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include "qdpack_matrix_gsl.h"
#include "gsl_ext_expm_complex.h"
#include "gsl_ext.h"

//==============================================================================
// thin abstraction layer for matrices
//==============================================================================

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
qdpack_matrix_t *
qdpack_matrix_alloc(size_t m, size_t n)
{
    qdpack_matrix_t *mat;

    if ((mat = (qdpack_matrix_t *)malloc(sizeof(qdpack_matrix_t))) == NULL)
    {
        fprintf(stderr, "%s: ERROR: Failed to allocate memory for matrix\n", __PRETTY_FUNCTION__);
        return NULL;
    }
   
    mat->data = gsl_matrix_complex_alloc(m, n);

    mat->m = m;
    mat->n = n;

    return mat;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void
qdpack_matrix_free(qdpack_matrix_t *mat)
{
    if (mat == NULL)
    {
        fprintf(stderr, "%s: ERROR: mat is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
   
    gsl_matrix_complex_free(mat->data);

    free(mat);
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
gsl_complex
qdpack_matrix_get(qdpack_matrix_t *mat, size_t i, size_t j)
{
    if (mat == NULL)
        return QDPACK_COMPLEX_ZERO;

    return gsl_matrix_complex_get(mat->data, i, j);
}

void
qdpack_matrix_set(qdpack_matrix_t *mat, size_t i, size_t j, qdpack_complex value)
{
    if (mat == NULL)
        return;

    return gsl_matrix_complex_set(mat->data, i, j, value);
}

void
qdpack_matrix_set_all(qdpack_matrix_t *mat, qdpack_complex value)
{
    if (mat == NULL)
        return;

    return gsl_matrix_complex_set_all(mat->data, value);
}

void
qdpack_matrix_set_zero(qdpack_matrix_t *mat)
{
    if (mat == NULL)
        return;

    return gsl_matrix_complex_set_zero(mat->data);
}

void
qdpack_matrix_set_identity(qdpack_matrix_t *mat)
{
    if (mat && mat->data)
    {    
        gsl_matrix_complex_set_identity(mat->data);
    }
}

void
qdpack_matrix_memcpy(qdpack_matrix_t *dst, qdpack_matrix_t *src)
{
    if (src == NULL || dst == NULL)
        return;

    dst->n = src->n;
    dst->m = src->m;
    
    gsl_matrix_complex_memcpy(dst->data, src->data);
}

void
qdpack_matrix_scale(qdpack_matrix_t *op, qdpack_complex z)
{
    if (op == NULL)
        return;

    gsl_matrix_complex_scale(op->data, z);
}

void
qdpack_matrix_add(qdpack_matrix_t *op1, qdpack_matrix_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    gsl_matrix_complex_add(op1->data, op2->data);
}

void
qdpack_matrix_sub(qdpack_matrix_t *op1, qdpack_matrix_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    gsl_matrix_complex_sub(op1->data, op2->data);
}

//
// BLAS style operator/matrix multiplication: C = alpha op(A) * op(B) + beta C
//
void
qdpack_matrix_blas_zgemm(QDPACK_TRANSPOSE_t transA,
                         QDPACK_TRANSPOSE_t transB, 
                         qdpack_complex alpha, 
                         qdpack_matrix_t *A, 
                         qdpack_matrix_t *B,
                         qdpack_complex beta, 
                         qdpack_matrix_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        gsl_blas_zgemm(transA, transB, alpha, A->data, B->data, beta, C->data);
    }
    else
    {
        fprintf(stderr, "%s: ERROR: NULL valued operator\n", __PRETTY_FUNCTION__);
    }
}

//
// plain operator/matrix multiplication: C = A * B
//
void
qdpack_matrix_multiply(qdpack_matrix_t *A, qdpack_matrix_t *B, qdpack_matrix_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        qdpack_matrix_blas_zgemm(CblasNoTrans, CblasNoTrans, QDPACK_COMPLEX_ONE, A, B, QDPACK_COMPLEX_ZERO, C);
    }
    else
    {
        fprintf(stderr, "%s: ERROR: NULL valued operator\n", __PRETTY_FUNCTION__);
    }
}

//
// XXX: implement with GSL matrices instead of QDPACK matrices
//
qdpack_complex
qdpack_matrix_multiply_vmv(qdpack_matrix_t *mat, qdpack_matrix_t *v)
{
    qdpack_matrix_t *prod, *expt;
    gsl_complex res;

    prod = qdpack_matrix_alloc(mat->m, 1);
    expt = qdpack_matrix_alloc(1, 1);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, GSL_COMPLEX_ONE, mat->data, v->data,  QDPACK_COMPLEX_ZERO, prod->data);
    gsl_blas_zgemm(CblasTrans,   CblasNoTrans, GSL_COMPLEX_ONE, v->data, prod->data, QDPACK_COMPLEX_ZERO, expt->data); // ConjTrans ??

    res = qdpack_matrix_get(expt, 0, 0);

    qdpack_matrix_free(expt);
    qdpack_matrix_free(prod);

    return res;
}

//
// Matrix exponentiation
//
void
qdpack_matrix_exp(qdpack_matrix_t *m, qdpack_matrix_t *exp_m)
{
    //gsl_ext_expm_complex(m->data, exp_m->data);
    gsl_linalg_complex_exponential_ss(m->data, exp_m->data, GSL_PREC_DOUBLE);
}


//
// dot product for vector-like matrices (should check that dimensions of both 
// are [m x 1])
//
qdpack_complex
qdpack_matrix_dot(qdpack_matrix_t *a, qdpack_matrix_t *b)
{
    // use gsl_blas_zdotc(a->data, b->data, &z); ?
    qdpack_matrix_t *prod;
    gsl_complex z;

    prod = qdpack_matrix_alloc(1, 1);

    gsl_blas_zgemm(CblasTrans, CblasNoTrans, GSL_COMPLEX_ONE, a->data, b->data, QDPACK_COMPLEX_ZERO, prod->data); // ConjTrans ??

    z = qdpack_matrix_get(prod, 0, 0);

    qdpack_matrix_free(prod);

    return z;
}

//==============================================================================
// I/O
//==============================================================================

//#define QDPACK_MATRIX_WRITE_FORMAT_REAL    0
//#define QDPACK_MATRIX_WRITE_FORMAT_IMAG    1
//#define QDPACK_MATRIX_WRITE_FORMAT_COMPLEX 2

int
qdpack_matrix_write(FILE *f, qdpack_matrix_t *M, int format)
{
    int i, j;

    if (f == NULL || M == NULL)
    {
        return -1;
    }
   
    fprintf(f, "# QDPACK [GSL] matrix: %d x %d \n", M->m, M->n);

    if (format == QDPACK_MATRIX_WRITE_FORMAT_PYTHON)
    {
        return qdpack_matrix_python_write(f, M);
    }
    
    for (i = 0; i < M->m; i++)
    {
        for (j = 0; j < M->n; j++)
        {
            if (format == QDPACK_MATRIX_WRITE_FORMAT_REAL)
            {
                fprintf(f, "%.6e%s", QDPACK_REAL(qdpack_matrix_get(M, i, j)), j == M->n-1 ? "" : " ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_IMAG)
            {
                fprintf(f, "%.6e%s", QDPACK_IMAG(qdpack_matrix_get(M, i, j)), j == M->n-1 ? "" : " ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_COMPLEX)
            {
                fprintf(f, "(%.6e,%.6e)%s", QDPACK_REAL(qdpack_matrix_get(M, i, j)), QDPACK_IMAG(qdpack_matrix_get(M, i, j)), j == M->n-1 ? "" : " ");
            }

        }
        fprintf(f, "\n");
    }

    return 0;
}

int
qdpack_matrix_python_write(FILE *f, qdpack_matrix_t *M)
{
    int i, j;

    if (f == NULL || M == NULL)
    {
        return -1;
    }
    
    fprintf(f, "array([\n");
    for (i = 0; i < M->m; i++)
    {
        fprintf(f, "[ ");
        for (j = 0; j < M->n; j++)
        {
            fprintf(f, "(%.6e+%.6ej)%s", QDPACK_REAL(qdpack_matrix_get(M, i, j)), QDPACK_IMAG(qdpack_matrix_get(M, i, j)), j == M->n-1 ? "" : ", ");
        }
        fprintf(f, "],\n");
    }
    fprintf(f, "])\n");

    return 0;
}


//==============================================================================
// eigenstates and eigenvalues
//==============================================================================

//
// Hermitian matrices
//
static gsl_eigen_hermv_workspace *hermv_workspace = NULL;
int
qdpack_matrix_eigen_hermv(qdpack_matrix_t *m,
                          qdpack_matrix_t *eval,
                          qdpack_matrix_t *evec, 
                          int sort_order)
{
    int i;
    gsl_vector *eval_vector;

    if (hermv_workspace == NULL)
    {
        hermv_workspace = gsl_eigen_hermv_alloc(m->m);
    }

    eval_vector = gsl_vector_alloc(eval->m);

    gsl_eigen_hermv(m->data, eval_vector, evec->data, hermv_workspace); 
    gsl_eigen_hermv_sort(eval_vector, evec->data, sort_order);

    for (i = 0; i < eval->m; i++)
    {   
        gsl_complex z;
        GSL_SET_COMPLEX(&z, gsl_vector_get(eval_vector, i), 0.0);
        qdpack_matrix_set(eval, i, 0, z);
    }

    return 0;
}

//
// General matrices
//
int
qdpack_matrix_eigen_zgeev(qdpack_matrix_t *m,
                          qdpack_matrix_t *eval,
                          qdpack_matrix_t *evec,
                          int sort_order)
{
    int i;
    gsl_vector_complex *eval_vector;

    eval_vector = gsl_vector_complex_alloc(eval->m);

    gsl_ext_eigen_zgeev(m->data, evec->data, eval_vector);
    
    gsl_ext_eigen_sort(evec->data, eval_vector, sort_order);

    for (i = 0; i < eval->m; i++)
    {   
        qdpack_matrix_set(eval, i, 0, gsl_vector_complex_get(eval_vector, i));
    }

    return 0;
}



