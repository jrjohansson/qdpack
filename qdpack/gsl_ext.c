//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#include <math.h>

#include "gsl_ext.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <f2c.h>


// XXX: find the correct header file to included instead of this declaration
int zgeev_(char *jobvl, char *jobvr, integer *n, 
	doublecomplex *a, integer *lda, doublecomplex *w, doublecomplex *vl, 
	integer *ldvl, doublecomplex *vr, integer *ldvr, doublecomplex *work, 
	integer *lwork, doublereal *rwork, integer *info);

/*
 * Convert a real valued GSL matrix to a complex valued GSL matrix.
 *
 */
int
gsl_ext_matrix_convert(gsl_matrix_complex *dst, gsl_matrix *src)
{
    unsigned int i, j;
    gsl_complex z;
    
    if (src->size1 != dst->size1 || src->size2 != dst->size2)
        return -1;

    for (i = 0; i < src->size1; i++)
    {
        for (j = 0; j < src->size1; j++)
        {
            GSL_SET_COMPLEX(&z, gsl_matrix_get(src, i, j), 0.0);
            gsl_matrix_complex_set(dst, i, j, z);
        }
    }
    
    return 0;
}

/* 
 * Allocate a complex matrix, convert the real matrix it complex form
 * and free the real matrix.
 *
 */
gsl_matrix_complex *
gsl_ext_matrix_convert_and_free(gsl_matrix *m)
{
    gsl_matrix_complex *zm;

    zm = gsl_matrix_complex_alloc(m->size1, m->size2);

    gsl_ext_matrix_convert(zm, m);

    gsl_matrix_free(m);
    
    return zm;
}

/*
 * Calculate the matrix exponent of A and store in eA.
 * Algorithm: Truncated Talyor series.
 *
 * WARNING: Large errors possible and it's slow.
 */
int 
gsl_ext_expm_complex(gsl_matrix_complex *A, gsl_matrix_complex *eA)
{
    int i;
    gsl_complex alpha, beta, z;
    gsl_matrix_complex *I, *T;
    
    I = gsl_matrix_complex_alloc(A->size1, A->size2);
    T = gsl_matrix_complex_alloc(A->size1, A->size2);
    
    GSL_SET_COMPLEX(&alpha, 1.0, 0.0);
    GSL_SET_COMPLEX(&beta,  0.0, 0.0);
    
    gsl_matrix_complex_set_identity(I);
    gsl_matrix_complex_set_identity(eA);
    
    for (i = 50; i > 0; i--)
    {
        GSL_SET_COMPLEX(&z, 1.0 / i, 0.0);
        gsl_matrix_complex_scale(eA, z);
        
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, eA, A, beta, T);    
        gsl_matrix_complex_add(T, I);
        gsl_matrix_complex_memcpy(eA, T);
    }
    
    return 0;
}

/*
 * Calculate eigenvectors and eigenvalues for a non-symmetric complex matrix
 * using the CLAPACK zgeev_ function.
 *
 */
int
gsl_ext_eigen_zgeev(gsl_matrix_complex *A_gsl, gsl_matrix_complex *evec, gsl_vector_complex *eval)
{
    //integer *pivot;
    integer n,i,j,info,lwork, ldvl, ldvr, lda;
    doublecomplex *A,*vr,*vl,*w,*work;
    doublereal *rwork;
    char jobvl,jobvr;

    n = A_gsl->size1;
    
//  pivot = (integer *)malloc((size_t)n * sizeof(int));
    A  = (doublecomplex *)malloc((size_t)n * n * sizeof(doublecomplex));
    w  = (doublecomplex *)malloc((size_t)n * sizeof(doublecomplex));
    vr = (doublecomplex *)malloc((size_t)n * n * sizeof(doublecomplex));
    vl = (doublecomplex *)malloc((size_t)n * n * sizeof(doublecomplex));
    lwork = 16 * n;
    work  = (doublecomplex *)malloc((size_t)lwork * sizeof(doublecomplex));
    rwork  = (doublereal *)malloc((size_t)lwork * sizeof(doublereal));

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            gsl_complex z;
            double re,im;
        
            z = gsl_matrix_complex_get(A_gsl, i, j);
            re = GSL_REAL(z);
            im = GSL_IMAG(z);

            A[j*n+i] = (doublecomplex){re,im};
        }
    }

    jobvl='N';
    jobvr='V';

    lda  = n;
    ldvr = n;
    ldvl = n;
    zgeev_(&jobvl, &jobvr, &n, A, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info);
    
    //ZGEEVX(BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, W, VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM, RCONDE, RCONDV, WORK, LWORK, RWORK, INFO ) 

    for (i = 0; i < n; i++)
    {
        gsl_complex z;
        GSL_SET_COMPLEX(&z, w[i].r, w[i].i);
        gsl_vector_complex_set(eval, i, z);
    }

    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            gsl_complex z;
            
            GSL_SET_COMPLEX(&z, vr[j*n+i].r, vr[j*n+i].i);
            
            gsl_matrix_complex_set(evec, i, j, z);
        }
    }

    if (info != 0)
    {
        printf("zgeev_: error: info = %d\n", (int)info);
    }
    
//    free(pivot);
    free(A);
    free(w);
    free(vr);
    free(vl);
    free(work);
    free(rwork);
    
    return 0;
}

/** ----------------------------------------------------------------------------
 * Sort eigenvectors
 */
int
gsl_ext_eigen_sort(gsl_matrix_complex *evec, gsl_vector_complex *eval, int sort_order)
{
    gsl_complex z;

      gsl_matrix_complex *evec_copy;
    gsl_vector_complex *eval_copy;
    
    int *idx_map, i, j, idx_temp;

    double p1, p2;

    if ((evec->size1 != evec->size2) || (evec->size1 != eval->size))
    {
        return -1;
    }

    evec_copy = gsl_matrix_complex_alloc(evec->size1, evec->size2);
    eval_copy = gsl_vector_complex_alloc(eval->size);

    idx_map  = (int *)malloc(sizeof(int) * eval->size);

    gsl_matrix_complex_memcpy(evec_copy, evec);
    gsl_vector_complex_memcpy(eval_copy, eval);

    // calculate new eigenvalue order

    for (i = 0; i < eval->size; i++)
    {
        idx_map[i] = i;
    }

    for (i = 0; i < eval->size - 1; i++)
    {
        for (j = i+1; j < eval->size; j++)
        {
            idx_temp = -1;

            if (sort_order == GSL_EXT_EIGEN_SORT_ABS)
            {
                p1 = gsl_complex_abs(gsl_vector_complex_get(eval, idx_map[i]));
                p2 = gsl_complex_abs(gsl_vector_complex_get(eval, idx_map[j]));
                if (p1 > p2)
                {
                    idx_temp = idx_map[i];
                }
            }

            if (sort_order == GSL_EXT_EIGEN_SORT_PHASE)
            {
                p1 = gsl_complex_arg(gsl_vector_complex_get(eval, idx_map[i]));
                p2 = gsl_complex_arg(gsl_vector_complex_get(eval, idx_map[j]));

                if (p1 >   M_PI) p1 -= 2*M_PI;        
                if (p1 <= -M_PI) p1 += 2*M_PI;        

                if (p2 >   M_PI) p2 -= 2*M_PI;        
                if (p2 <= -M_PI) p2 += 2*M_PI;        

                //if (((p2 < p1) && (p1 - p2 < M_PI)) || )
                if (p2 < p1)
                {
                    idx_temp = idx_map[i];
                }
            }

            if (idx_temp != -1)
            {
                // swap
                //idx_temp = idx_map[i];
                idx_map[i] = idx_map[j];
                idx_map[j] = idx_temp;
            }
        }
    } 

    // reshuffle the eigenvectors and eigenvalues
    for (i = 0; i < eval->size; i++)
    {   
        for (j = 0; j < eval->size; j++)
        {
            z = gsl_matrix_complex_get(evec_copy, idx_map[i], j);
            gsl_matrix_complex_set(evec, i, j, z);            
            //z = gsl_matrix_complex_get(evec_copy, i, idx_map[j]);
            //gsl_matrix_complex_set(evec, i, j, z);                    
        }

        z = gsl_vector_complex_get(eval_copy, idx_map[i]);
        gsl_vector_complex_set(eval, i, z);
    }    

    gsl_matrix_complex_free(evec_copy);
    gsl_vector_complex_free(eval_copy);

    free(idx_map);
 
    return 0;
}

/** ----------------------------------------------------------------------------
 * Normalize eigenvectors
 * /
int
gsl_ext_eigen_norm(gsl_matrix_complex *evec)
{
    gsl_complex z, w;

    if ((evec->size1 != evec->size2))
    {
        return -1;
    }

    // normalize each column
    for (i = 0; i < evec->size1; i++)
    {   
        GSL_SET_COMPLEX(&z, 0.0, 0.0);

        for (j = 0; j < evec->size2; j++)
        {
            z = gsl_complex_add(z, gsl_matrix_complex_get(evec, idx_i, j));
        }

        for (j = 0; j < evec->size2; j++)
        {
            z = gsl_complex_add(z, gsl_matrix_complex_get(evec, idx_i, j));
        }

            gsl_matrix_complex_set(evec, i, j, z);                    
        z = gsl_vector_complex_get(eval_copy, idx_map[i]);
        gsl_vector_complex_set(eval, i, z);
    }    
 
    return 0;
}
*/















