//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

// Description:
//
// A CXSparse based matrix backend for QDpack.
//
#include <string.h>
#include <suitesparse/cs.h>
#include "qdpack_matrix_cxsparse.h"

#define DEBUG 0

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

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if ((mat = (qdpack_matrix_t *)malloc(sizeof(qdpack_matrix_t))) == NULL)
    {
        fprintf(stderr, "%s: ERROR: Failed to allocate memory for matrix\n", __PRETTY_FUNCTION__);
        return NULL;
    }
   
    mat->data = cs_cl_spalloc(m, n, 2*n, 1, 1); // start with triplet form

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
    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
    {
        fprintf(stderr, "%s: ERROR: mat is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
   
    cs_cl_spfree(mat->data);

    free(mat);
}

//------------------------------------------------------------------------------
// convert from compressed format to triplet
//------------------------------------------------------------------------------
void
qdpack_matrix_triplet(qdpack_matrix_t *mat)
{
    int p_idx, i_idx;

    cs_cl *A;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (DEBUG)
        printf("%s: converting from compressed to triplet [%ld entries]\n", __PRETTY_FUNCTION__, mat->data->p[mat->data->n]);

    A = cs_cl_spalloc(mat->m, mat->n, mat->data->p[mat->data->n], 1, 1);

    A->nz = 0;
    for (p_idx = 0; p_idx < mat->data->n; p_idx++)
    {
       // printf("%s: j = %d -> i in [%d, %d] (%d)\n", __PRETTY_FUNCTION__, p_idx, mat->data->p[p_idx], mat->data->p[p_idx+1], A->nz);

        for (i_idx = mat->data->p[p_idx]; i_idx < mat->data->p[p_idx+1]; i_idx++)
        {
            A->p[A->nz] = p_idx;
            A->i[A->nz] = mat->data->i[i_idx];
            A->x[A->nz] = mat->data->x[i_idx];

            //printf("%s: M[%d, %d] (%d) = (%f, %f)\n", __PRETTY_FUNCTION__, mat->data->i[i_idx], p_idx, A->nz, QDPACK_REAL(mat->data->x[i_idx]), QDPACK_IMAG(mat->data->x[i_idx]));

            A->nz++;
        }
    }

    cs_cl_spfree(mat->data);

    mat->data = A;    
}

//------------------------------------------------------------------------------
// Convert from triplet to compressed column format
//------------------------------------------------------------------------------
void
qdpack_matrix_compressed(qdpack_matrix_t *mat)
{
    cs_cl *A;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (DEBUG)
        printf("%s: converting from triplet to compressed [%ld]\n", __PRETTY_FUNCTION__, mat->data->nz);

    A = cs_cl_compress(mat->data);

    cs_cl_spfree(mat->data);

    mat->data = A;    
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
qdpack_complex
qdpack_matrix_get(qdpack_matrix_t *mat, size_t i, size_t j)
{
    int n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
        return QDPACK_COMPLEX_ZERO;

    if (CS_TRIPLET(mat->data))
    {
        // triplet form

        for (n = 0; n < mat->data->nz; n++)
        {
            if (mat->data->i[n] == i && mat->data->p[n] == j)
            {
                return mat->data->x[n];
            }
        }
    }    
    else
    {
        // compressed column format

        for (n = mat->data->p[j]; n <= mat->data->p[j+1]; n++)
        {
            if (mat->data->i[n] == i)
            {
                return mat->data->x[n];
            }
        }
    }

    return QDPACK_COMPLEX_ZERO;
}

void
qdpack_matrix_set(qdpack_matrix_t *mat, size_t i, size_t j, qdpack_complex value)
{
    int n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
        return;

    if (CS_TRIPLET(mat->data))
    {
        // triplet form

        for (n = 0; n < mat->data->nz; n++)
        {
            if (mat->data->p[n] == j && mat->data->i[n] == i)
            {
                mat->data->x[n] = value;
                return;
            }
        }

        cs_cl_entry(mat->data, i, j, value);
    }    
    else
    {
        // compressed column format

        for (n = mat->data->p[j]; n <= mat->data->p[j+1]; n++)
        {
            if (mat->data->i[n] == i)
            {
                mat->data->x[n] = value;
                return;
            }
        }

        qdpack_matrix_compressed(mat);
        qdpack_matrix_set(mat, i, j, value);
    }
}

void
qdpack_matrix_set_all(qdpack_matrix_t *mat, qdpack_complex value)
{
    int m, n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
        return;

    qdpack_matrix_triplet(mat);

    mat->data->nz = 0;

    for (m = 0; m < mat->data->m; m++)
    {
        for (n = 0; n < mat->data->n; n++)
        {
            cs_cl_entry(mat->data, m, n, value);
        }
    }
}

void
qdpack_matrix_set_zero(qdpack_matrix_t *mat)
{
    int n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
        return;

    if (CS_TRIPLET(mat->data))
    {
        // triplet form
        mat->data->nz = 0;
    }    
    else
    {

        mat->data->nz = 0; // converts to triplet
        // compressed column format
        for (n = 0; n < mat->data->n + 1; n++)
        {
            mat->data->p[n] = 0;
        }
    }
}

void
qdpack_matrix_set_identity(qdpack_matrix_t *mat)
{
    int n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (!(mat && mat->data))
    {    
        // warning..
        return;
    }

    if (CS_TRIPLET(mat->data))
    {
        // triplet form

        for (n = 0; n < mat->data->n; n++)
        {
            mat->data->i[n] = n;
            mat->data->p[n] = n;
            mat->data->x[n] = QDPACK_COMPLEX_ONE;
        }
        mat->data->nz = n;
    }    
    else
    {
        // compressed column format

        for (n = 0; n < mat->data->n; n++)
        {
            mat->data->i[n] = n;
            mat->data->p[n] = n;
            mat->data->x[n] = QDPACK_COMPLEX_ONE;
        }
        mat->data->p[n] = n;
    }

}

void
qdpack_matrix_memcpy(qdpack_matrix_t *dst, qdpack_matrix_t *src)
{
    int size;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (src == NULL || dst == NULL)
        return;

    dst->n = src->n;
    dst->m = src->m;

    //
    // free up old
    //
    cs_free(dst->data->i);
    cs_free(dst->data->p);
    cs_free(dst->data->x);

    //
    // alloc and copy 
    //
    dst->data->i = (long *)malloc(src->data->nzmax * sizeof(long));
    memcpy((void *)dst->data->i, (void *)src->data->i, src->data->nzmax * sizeof(long));

    size = CS_TRIPLET(src->data) ? src->data->nzmax : src->data->n + 1;
    dst->data->p = (long *)malloc(size * sizeof(long));
    memcpy((void *)dst->data->p, (void *)src->data->p, size * sizeof(long)); 

    dst->data->x = (qdpack_complex *)malloc(src->data->nzmax * sizeof(qdpack_complex));
    memcpy((void *)dst->data->x, (void *)src->data->x, src->data->nzmax * sizeof(qdpack_complex));

    dst->data->n = src->data->n;
    dst->data->m = src->data->m;
    dst->data->nzmax = src->data->nzmax;
    dst->data->nz = src->data->nz;
}

void
qdpack_matrix_scale(qdpack_matrix_t *mat, qdpack_complex z)
{
    int n;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (mat == NULL)
        return;

    if (CS_TRIPLET(mat->data))
    {
        for (n = 0; n < mat->data->nz; n++)
        {
            mat->data->x[n] *= z;
        }
    }
    else
    {
        for (n = 0; n < mat->data->p[mat->data->n]; n++)
        {
            mat->data->x[n] *= z;
        }        
    }
}

void
qdpack_matrix_add(qdpack_matrix_t *op1, qdpack_matrix_t *op2)
{
    cs_cl *C;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
 
    if (!CS_CSC(op1->data))
    {
        qdpack_matrix_compressed(op1);
    }

    if (!CS_CSC(op2->data))
    {
        qdpack_matrix_compressed(op2);
    }   

    C = cs_cl_add(op1->data, op2->data, QDPACK_COMPLEX_ONE, QDPACK_COMPLEX_ONE);
      
    cs_cl_spfree(op1->data);

    op1->data = C;
}

void
qdpack_matrix_sub(qdpack_matrix_t *op1, qdpack_matrix_t *op2)
{
    cs_cl *C;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);


    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
 
    if (!CS_CSC(op1->data))
    {
        qdpack_matrix_compressed(op1);
    }

    if (!CS_CSC(op2->data))
    {
        qdpack_matrix_compressed(op2);
    }   

    C = cs_cl_add(op1->data, op2->data, QDPACK_COMPLEX_ONE, -QDPACK_COMPLEX_ONE);
      
    cs_cl_spfree(op1->data);

    op1->data = C;
}

//
// BLAS style operator/matrix multiplication: C = alpha op(A) * op(B) + beta C
//
// XXX: direct implementation
//
void
qdpack_matrix_blas_zgemm(QDPACK_TRANSPOSE_t op_A,
                         QDPACK_TRANSPOSE_t op_B, 
                         qdpack_complex alpha, 
                         qdpack_matrix_t *A, 
                         qdpack_matrix_t *B,
                         qdpack_complex beta, 
                         qdpack_matrix_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        cs_cl *A_op, *B_op, *X, *Y;

        if (!CS_CSC(A->data))
        {
            qdpack_matrix_compressed(A);
        }

        if (!CS_CSC(B->data))
        {
            qdpack_matrix_compressed(B);
        }   

        if (!CS_CSC(C->data))
        {
            qdpack_matrix_compressed(C);
        }   

        switch (op_A)
        {
            case QDpackNoTrans:
                A_op = A->data;
                break;

            case QDpackTrans:
                A_op = cs_cl_transpose(A->data, 1);
                break;

            case QDpackConjTrans:
                A_op = cs_cl_transpose(A->data, 1); // XXX
                break;

            default:
                fprintf(stderr, "%s: Error: illegal value of op_A.\n", __PRETTY_FUNCTION__);
                return;
        }

        switch (op_B)
        {
            case QDpackNoTrans:
                B_op = B->data;
                break;

            case QDpackTrans:
                B_op = cs_cl_transpose(B->data, 1);
                break;

            case QDpackConjTrans:
                B_op = cs_cl_transpose(B->data, 1);
                break;

            default:
                fprintf(stderr, "%s: Error: illegal value of op_A.\n", __PRETTY_FUNCTION__);
                return;
        }
        

        X = cs_cl_multiply(A_op, B_op);

        Y = cs_cl_add(X, C->data, alpha, beta);
        
        cs_cl_spfree(X);

        cs_cl_spfree(C->data);

        C->data = Y;

        if (op_A != QDpackNoTrans)
        {
            cs_cl_spfree(A_op);
        }

        if (op_B != QDpackNoTrans)
        {
            cs_cl_spfree(B_op);
        }
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
    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    if (A && A->data && B && B->data && C && C->data)
    {
        cs_cl *C_new;

        if (!CS_CSC(A->data))
        {
            qdpack_matrix_compressed(A);
        }

        if (!CS_CSC(B->data))
        {
            qdpack_matrix_compressed(B);
        }            

        C_new = cs_cl_multiply(A->data, B->data);
          
        cs_cl_spfree(C->data);

        C->data = C_new;
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
qdpack_matrix_multiply_vmv(qdpack_matrix_t *A, qdpack_matrix_t *v)
{
    int n;

    qdpack_complex *x, *y;

    qdpack_complex res = QDPACK_COMPLEX_ZERO;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    //
    // check inputs and allocate working memory
    //
    if (A->n != A->m || A->n != v->m)
    {
        return QDPACK_COMPLEX_ZERO;
    }

    if ((x = (qdpack_complex *)calloc(sizeof(qdpack_complex), v->m)) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate vector x\n", __PRETTY_FUNCTION__);
    }

    if ((y = (qdpack_complex *)calloc(sizeof(qdpack_complex), v->m)) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate vector y\n", __PRETTY_FUNCTION__);
    }

    //
    // convert sparse vector v to dense vector: v to x
    //
    if (CS_TRIPLET(v->data))
    {
        for (n = 0; n < v->data->nz; n++)
        {
            x[v->data->i[n]] += v->data->x[n];
        }
    }
    else
    {
        for (n = 0; n < v->data->p[1]; n++)
        {
            x[v->data->i[n]] += v->data->x[n];
        }
    }

    //
    // compute
    //
    cs_cl_gaxpy(A->data, x, y);

    // vector dot product
    for (n = 0; n < v->m; n++)
    {
        res += qdpack_complex_conjugate(x[n]) * y[n]; // conjugate ?
    }

    free(x);
    free(y);

    return res;
}

//
// Matrix exponentiation
//
void
qdpack_matrix_exp(qdpack_matrix_t *m, qdpack_matrix_t *exp_m)
{
    //gsl_ext_expm_complex(m->data, exp_m->data);
    //gsl_linalg_complex_exponential_ss(m->data, exp_m->data, GSL_PREC_DOUBLE);
    printf("%s: WARNING: Not Yet Implemented.\n", __PRETTY_FUNCTION__);

}


//
// dot product for vector-like matrices (should check that dimensions of both 
// are [m x 1])
//
qdpack_complex
qdpack_matrix_dot(qdpack_matrix_t *a, qdpack_matrix_t *b)
{
    int n;

    qdpack_complex *x, *y;

    qdpack_complex res = QDPACK_COMPLEX_ZERO;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    //
    // check inputs and allocate working memory
    //
    if (a->m != b->m || a->n != b->n || a->n != 1)
    {
        return QDPACK_COMPLEX_ZERO;
    }

    x = qdpack_vector_dense(a);
    y = qdpack_vector_dense(b);

    // vector dot product
    for (n = 0; n < a->m; n++)
    {
        res += qdpack_complex_conjugate(x[n]) * y[n];
    }

    free(x);
    free(y);

    return res;
}

qdpack_complex *
qdpack_vector_dense(qdpack_matrix_t *v)
{
    int n;
    qdpack_complex *x;

    if (DEBUG)
        printf("%s: DEBUG\n", __PRETTY_FUNCTION__);

    //
    // check inputs and allocate working memory
    //
    if (v->n != 1)
    {
        return NULL;
    }

    if ((x = (qdpack_complex *)calloc(sizeof(qdpack_complex), v->m)) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate vector x\n", __PRETTY_FUNCTION__);
    }

    //
    // convert sparse vector v to dense vector: v to x
    //
    if (CS_TRIPLET(v->data))
    {
        for (n = 0; n < v->data->nz; n++)
        {
            x[v->data->i[n]] += v->data->x[n];
        }
    }
    else
    {
        for (n = 0; n < v->data->p[1]; n++)
        {
            x[v->data->i[n]] += v->data->x[n];
        }
    }

    return x;
}

//==============================================================================
// I/O
//==============================================================================

int
qdpack_matrix_write(FILE *f, qdpack_matrix_t *M, int format)
{
    qdpack_complex z;
    int i, j;

    if (f == NULL || M == NULL)
    {
        return -1;
    }

    fprintf(f, "# QDPACK [CXSparse] [%s] matrix: %d x %d (density %.3f, n = %ld) \n", 
            CS_TRIPLET(M->data) ? "triplet" : "csc", M->m, M->n,
            CS_TRIPLET(M->data) ? (1.0 * M->data->nz) / (M->data->m * M->data->n) : (1.0 * M->data->p[M->data->n]) / (M->data->m * M->data->n),
            CS_TRIPLET(M->data) ? M->data->nz : M->data->p[M->data->n]);
    
    if (format == QDPACK_MATRIX_WRITE_FORMAT_PYTHON)
    {
        return qdpack_matrix_python_write(f, M);
    }
    
    for (i = 0; i < M->m; i++)
    {
        for (j = 0; j < M->n; j++)
        {
            z = qdpack_matrix_get(M, i, j);

            if (format == QDPACK_MATRIX_WRITE_FORMAT_REAL)
            {
                fprintf(f, "%.6e%s", QDPACK_REAL(z), j == M->n-1 ? "" : " ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_IMAG)
            {
                fprintf(f, "%.6e%s", QDPACK_IMAG(z), j == M->n-1 ? "" : " ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_COMPLEX)
            {
                fprintf(f, "(%.6e,%.6e)%s", QDPACK_REAL(z), QDPACK_IMAG(z), j == M->n-1 ? "" : " ");
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


int
qdpack_matrix_sparse_write(FILE *f, qdpack_matrix_t *M, int format)
{
    int i, j;

    if (f == NULL || M == NULL)
    {
        return -1;
    }
    
    fprintf(f, "# QDPACK [CXSparse] [%s] matrix: %d x %d \n", 
            CS_TRIPLET(M->data) ? "triplet" : "csc", M->m, M->n);
    

    if (CS_CSC(M->data))
    {
        fprintf(f, "p = [");
        for (j = 0; j < M->data->n; j++)
        {
            fprintf(f, "%ld%s", M->data->p[j], j == M->data->n-1 ? "" : ", ");
        }            
        fprintf(f, "]\n");


        fprintf(f, "i = [");
        for (i = 0; i < M->data->p[M->data->n]; i++)
        {
            fprintf(f, "%ld%s", M->data->i[i], i == M->data->p[M->data->n]-1 ? "" : ", ");
        }            
        fprintf(f, "]\n");


        fprintf(f, "x = [");
        for (i = 0; i < M->data->p[M->data->n]; i++)
        {
            if (format == QDPACK_MATRIX_WRITE_FORMAT_REAL)
            {
                fprintf(f, "%.6e%s", QDPACK_REAL(M->data->x[i]), i == M->data->p[M->data->n]-1 ? "" : ", ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_IMAG)
            {
                fprintf(f, "%.6e%s", QDPACK_IMAG(M->data->x[i]), i == M->data->p[M->data->n]-1 ? "" : ", ");
            }
            else //if (format == QDPACK_MATRIX_WRITE_FORMAT_COMPLEX)
            {
                fprintf(f, "(%.6e,%.6e)%s", QDPACK_REAL(M->data->x[i]), QDPACK_IMAG(M->data->x[i]), i == M->data->p[M->data->n]-1 ? "" : ", ");
            }

        }            
        fprintf(f, "]\n");
    }


    if (CS_TRIPLET(M->data))
    {
        fprintf(f, "p = [");
        for (j = 0; j < M->data->nz; j++)
        {
            fprintf(f, "%ld%s", M->data->p[j], j == M->data->nz-1 ? "" : ", ");
        }            
        fprintf(f, "]\n");

        fprintf(f, "i = [");
        for (j = 0; j < M->data->nz; j++)
        {
            fprintf(f, "%ld%s", M->data->i[j], j == M->data->nz-1 ? "" : ", ");
        }            
        fprintf(f, "]\n");

        fprintf(f, "x = [");
        for (i = 0; i < M->data->nz; i++)
        {
            if (format == QDPACK_MATRIX_WRITE_FORMAT_REAL)
            {
                fprintf(f, "%.6e%s", QDPACK_REAL(M->data->x[i]), i == M->data->nz-1 ? "" : ", ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_IMAG)
            {
                fprintf(f, "%.6e%s", QDPACK_IMAG(M->data->x[i]), i == M->data->nz-1 ? "" : ", ");
            }
            else if (format == QDPACK_MATRIX_WRITE_FORMAT_COMPLEX)
            {
                fprintf(f, "(%.6e,%.6e)%s", QDPACK_REAL(M->data->x[i]), QDPACK_IMAG(M->data->x[i]), i == M->data->nz-1 ? "" : ", ");
            }

        }            
        fprintf(f, "]\n");
    }

    return 0;
}

//==============================================================================
// eigenstates and eigenvalues
//==============================================================================

//
// Hermitian matrices
//
int
qdpack_matrix_eigen_hermv(qdpack_matrix_t *m,
                          qdpack_matrix_t *eval,
                          qdpack_matrix_t *evec, 
                          int sort_order)
{
    printf("%s: WARNING: Not Yet Implemented.\n", __PRETTY_FUNCTION__);
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
    printf("%s: WARNING: Not Yet Implemented.\n", __PRETTY_FUNCTION__);
    return 0;
}



