//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------


#include <stdio.h>
#include <math.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include <qdpack/qdpack.h>

/**
 * \breif    Calculate the von Neumann Entropy of a system with density matrix rho.
 *
 * This function takes a density matrix as input and calculates the von Neumann entropy.
 * Note that the density matrix is not restricted to be representing any particular
 * for of quantum system.
 *
 * @param    rho    The density matrix for which the von Neumann entropy is to be calculated.
 *
 */
double
qdpack_entanglement_neumann_entropy(qdpack_operator_t *rho)
{
    int i;
    double e, ne = 0.0;
    qdpack_matrix_t *rho_w;
    qdpack_matrix_t *eval, *evec;
    
    rho_w = qdpack_matrix_alloc(rho->m, rho->n);
    qdpack_matrix_memcpy(rho_w, rho->data);
        
    eval = qdpack_matrix_alloc(rho->m, 1);
    evec = qdpack_matrix_alloc(rho->m, rho->n);

    qdpack_matrix_eigen_hermv(rho_w, eval, evec, 0);

    for (i = 0; i < rho->m; i++)
    {
        e = QDPACK_REAL(qdpack_matrix_get(eval, i, 0));
        //ne += - e * log2(e);
        ne += - e * log(e);
    }

    /* clean-up */
    qdpack_matrix_free(rho_w);
    qdpack_matrix_free(eval);
    qdpack_matrix_free(evec);

    return ne;
}

/**
 * \breif    Calculate the log negativity for the density matrix of two TLS, rho.
 *
 * Calculate the negativity of a two-qubit system. 
 *
 * @param    rho    Density matrix for a two-qubit system.
 *
 */
double
qdpack_entanglement_log_neg(qdpack_operator_t *rho)
{
    int i;
    double ln, sum_eval;
    qdpack_matrix_t *rho_pt;

    qdpack_matrix_t *eval, *evec;

    /* -- calculate partial transpose of rho, and store in rho_pt 
     *
     * | x x A B |      | x x a b |
     * | x x C D |  ->  | x x c d |
     * | a b y y |      | A B y y |
     * | c d y y |      | C D y y |
     *
     */
    rho_pt = qdpack_matrix_alloc(4,4);
    qdpack_matrix_memcpy(rho_pt, rho->data); 
    // A,B,C,D -> a,b,c,d
    qdpack_matrix_set(rho_pt, 0, 2, qdpack_operator_get(rho, 2, 0));
    qdpack_matrix_set(rho_pt, 0, 3, qdpack_operator_get(rho, 2, 1));
    qdpack_matrix_set(rho_pt, 1, 2, qdpack_operator_get(rho, 3, 0));
    qdpack_matrix_set(rho_pt, 1, 3, qdpack_operator_get(rho, 3, 1));
    // a,b,c,d -> A,B,C,D
    qdpack_matrix_set(rho_pt, 2, 0, qdpack_operator_get(rho, 0, 2));
    qdpack_matrix_set(rho_pt, 2, 1, qdpack_operator_get(rho, 0, 3));
    qdpack_matrix_set(rho_pt, 3, 0, qdpack_operator_get(rho, 1, 2));
    qdpack_matrix_set(rho_pt, 3, 1, qdpack_operator_get(rho, 1, 3));

    /* --- diagonalize rho_pt --- */
    eval = qdpack_matrix_alloc(4, 1);
    evec = qdpack_matrix_alloc(4, 4);

    qdpack_matrix_eigen_hermv(rho_pt, eval, evec, 0);
    
    /* --- calculate the log negativity --- */
    sum_eval = 0;
    for (i = 0; i < 4; i++)
    {
        sum_eval += qdpack_complex_abs(qdpack_matrix_get(eval, i, 0));
    }
    ln = log2(sum_eval);

    /* clean-up */
    qdpack_matrix_free(rho_pt);
    qdpack_matrix_free(eval);
    qdpack_matrix_free(evec);

    return ln;
}


/**
 * \breif    Calculate the concurrence for the density matrix (two TLS) rho.
 *
 * Calculate the concurrence for the density matrix of a two-qubit system.
 *
 * @param    qs    Data structure that specifies the form of the quantum system.
 * @param    rho    The density matrix for a two-qubit system.
 *
 */
double
qdpack_entanglement_concurrence(qdpack_hilbert_space_t *qs, qdpack_operator_t *rho)
{
    double c = 0.0;
    qdpack_operator_t *y1, *y2, *y, *ws1, *ws2;
    qdpack_matrix_t *eval = qdpack_matrix_alloc(4, 1);
    qdpack_matrix_t *evec = qdpack_matrix_alloc(4, 4);

    y1  = qdpack_operator_alloc(qs);
    y2  = qdpack_operator_alloc(qs);
    
    operator_sigma_y(y1, qs, 0);
    operator_sigma_y(y2, qs, 1); 

    y   = qdpack_operator_alloc(qs);
    ws1 = qdpack_operator_alloc(qs);
    ws2 = qdpack_operator_alloc(qs);


    // calculate: rho * (y1 * y2) * rho' * (y1 * y2)
    //            rho *     y     * rho' *     y
    //            ---------------   -----------------
    //                  ws1                ws2
    //                   
    qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, y1, y2, QDPACK_COMPLEX_ZERO, y);
    
    qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, rho, y, QDPACK_COMPLEX_ZERO, ws1);
    qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, rho, y, QDPACK_COMPLEX_ZERO, ws2);
    
    qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, ws1, ws2, QDPACK_COMPLEX_ZERO, y);
    
    //printf("rr_real = \n"); qdpack_operator_print_real(y);
    //printf("rr_imag = \n"); qdpack_operator_print_imag(y);
    printf("rr = \n"); qdpack_operator_print(y);

    // calculate eigenvalues
    qdpack_matrix_eigen_zgeev(y->data, eval, evec, 0);

    {
        int i, j, n = 4;
        
        printf("Eigenvalues = \n");
        for (i = 0; i < n; i++)
        {
                qdpack_complex z;
                z = qdpack_matrix_get(eval, i, 0);
                printf("\t(%f, %f)", QDPACK_REAL(z), QDPACK_IMAG(z));
        }
        printf("\n");

        printf("evec = \n");
        for (i = 0; i < n; i++)
        {
                for (j = 0; j < n; j++)
                {
                        qdpack_complex z;
                        z = qdpack_matrix_get(evec, i, j);
                        printf("\t(%f, %f)", QDPACK_REAL(z), QDPACK_IMAG(z));
                }
                printf("\n");
        }
    }
    
    printf("c: eigvals: (%f, %f), (%f, %f), (%f, %f), (%f, %f)\n", 
           QDPACK_REAL(qdpack_matrix_get(eval, 0, 0)), QDPACK_IMAG(qdpack_matrix_get(eval, 0, 0)), 
           QDPACK_REAL(qdpack_matrix_get(eval, 1, 0)), QDPACK_IMAG(qdpack_matrix_get(eval, 1, 0)), 
           QDPACK_REAL(qdpack_matrix_get(eval, 2, 0)), QDPACK_IMAG(qdpack_matrix_get(eval, 2, 0)), 
           QDPACK_REAL(qdpack_matrix_get(eval, 3, 0)), QDPACK_IMAG(qdpack_matrix_get(eval, 3, 0)));

    c = + QDPACK_REAL(qdpack_matrix_get(eval, 0, 0)) - QDPACK_REAL(qdpack_matrix_get(eval, 1, 0))
        - QDPACK_REAL(qdpack_matrix_get(eval, 2, 0)) - QDPACK_REAL(qdpack_matrix_get(eval, 3, 0));

    c = (c < 0.0) ? 0.0 : c;

    qdpack_operator_free(y1);
    qdpack_operator_free(y2);
    qdpack_operator_free(y);
    qdpack_operator_free(ws1);
    qdpack_operator_free(ws2);
    qdpack_matrix_free(evec);
    qdpack_matrix_free(eval);
        
    return c;
}



