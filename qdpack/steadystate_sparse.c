//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Find the steady state of a quantum system with dissipation. 
 * Using sparse matrix algebra.
 */

#include <stdio.h>
#include <string.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <qdpack/qdpack.h>


/* ----------------------------------------------------------
 * Functions for calculating the Liouvillian.
 *
 */
static int
steady_state_L_add_op(gsl_sparse_matrix_complex *sop_L,
               qdpack_operator_t *op_left,
              qdpack_operator_t *op_right,
              qdpack_complex z, int N)
{
    qdpack_complex alpha, beta, gamma;

    int a, b, c, d, I, J;

    printf("DEBUG: steady_state_L_add_op: starting...");

    if (op_left != NULL && op_right == NULL)
    {
        /* 
         * Left operators:
         *
         * L_{ab}{cd} = O_{ac}KRON{b,d} 
         *
         * ->
         *
         * non-zero elements at L_{ab}{cb} = O_{ac}
             *
         */
        for (c = 0; c < N; c++)
        {
            for (a = 0; a < N; a++)
            {
                // alpha = op_left(a,c) * z
                                alpha = qdpack_operator_get(op_left, a, c);
                                alpha = qdpack_complex_mul(alpha, z);

                if (QDPACK_REAL(alpha) != 0.0 || QDPACK_IMAG(alpha) != 0.0)
                {
                    for (b = 0; b < N; b++)
                    {
                        I = ss_idx_major(N, a, b);
                        J = ss_idx_major(N, c, b);
    
                        // sop_L(I,J) = sop_L(I,J) + alpha
                        beta  = gsl_sparse_matrix_complex_get(sop_L, I, J);
                        gamma = qdpack_complex_add(alpha, beta);
                        gsl_sparse_matrix_complex_set(sop_L, I, J, gamma);
                    }
                }
            }
        }
    }
    
    
    /* 
     * right operators:
     *
     * L_{ab}{cd} = O_{bd}KRON{ca} 
     *
     * ->
     *
     * non-zero elements at L_{ab}{ad} = O_{bd}
         *
     */
    if (op_left == NULL && op_right != NULL)
    {
        for (d = 0; d < N; d++)
        {
            for (b = 0; b < N; b++)
            {
                // alpha = op_right(b,d) * z
                alpha = qdpack_operator_get(op_right, b, d);
                alpha = qdpack_complex_mul(alpha, z);
                
                if (QDPACK_REAL(alpha) != 0.0 || QDPACK_IMAG(alpha) != 0.0)
                {
                    for (a = 0; a < N; a++)
                    {
                        I = ss_idx_major(N, a, b);
                        J = ss_idx_major(N, a, d);
    
                        // sop_L(I,J) = sop_L(I,J) + alpha
                        beta  = gsl_sparse_matrix_complex_get(sop_L, I, J);
                        gamma = qdpack_complex_add(alpha, beta);
                        gsl_sparse_matrix_complex_set(sop_L, I, J, gamma);
                    }
                }
            }
        }
    }
    
    
    /* 
     * left and right operators simultaneously:
     *
     * L_{ab}{cd} = A_{ac}B_{db}
     *
     */
    if (op_left != NULL && op_right != NULL)
    {
        for (a = 0; a < N; a++)
        {
            for (c = 0; c < N; c++)
            {
                // alpha = op_left(a,c) * op_right(d,b) * z
                                alpha = qdpack_operator_get(op_left,  a, c);
                alpha = qdpack_complex_mul(alpha, z);
                
                if (QDPACK_REAL(alpha) != 0.0 || QDPACK_IMAG(alpha) != 0.0)
                {
                    for (b = 0; b < N; b++)
                    {
                        for (d = 0; d < N; d++)
                        {
                            I = ss_idx_major(N, a, b);
                            J = ss_idx_major(N, c, d);

                            beta  = qdpack_operator_get(op_right, d, b);

                            if (QDPACK_REAL(beta) != 0.0 || QDPACK_IMAG(beta) != 0.0)
                                            {
                                gamma = qdpack_complex_mul(alpha, beta);

                                                    // sop_L(I,J) = sop_L(I,J) + alpha
                                                    beta  = gsl_sparse_matrix_complex_get(sop_L, I, J);
                                                    gamma = qdpack_complex_add(gamma, beta);
                                                       gsl_sparse_matrix_complex_set(sop_L, I, J, gamma);
                            }
                        }
                    }
                }
            }
        }
    }
    
    printf("done.\n");
    
    return 0;
}

/*
 * Find the steady-state density matrix for a quantum system with dissipation by
 * calculate the system Liouville super operators, and solving L rho = 0, including
 * the density normalisation condition.
 *
 *
 */
int
qdpack_steadystate_dm_sparse(qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss, qdpack_simulation_t *param)
{
    gsl_sparse_matrix_complex *sop_L;
    double *rho_ss_vec_real, *rho_ss_vec_imag;
        double *ss_eq_rhs_real,  *ss_eq_rhs_imag;

    qdpack_complex z, alpha, beta;
    int ns, s, i;

    ns = qdpack_hilbert_space_nstates(qs);

    printf("DEBUG: solving for rho_ss: starting\n");

    /*
      * Allocate space for the L super-operator and
       * some working space vectors for solving the
     * LA problem.
      */
    sop_L             = gsl_sparse_matrix_complex_alloc(ns*ns, ns*ns, ns*ns);

    rho_ss_vec_real = (double *)malloc(sizeof(double) * ns*ns);
    rho_ss_vec_imag = (double *)malloc(sizeof(double) * ns*ns);
    ss_eq_rhs_real  = (double *)malloc(sizeof(double) * ns*ns);
    ss_eq_rhs_imag  = (double *)malloc(sizeof(double) * ns*ns);

    memset(rho_ss_vec_real, 0, sizeof(sizeof(double) * ns*ns));
    memset(rho_ss_vec_imag, 0, sizeof(sizeof(double) * ns*ns));
    memset(ss_eq_rhs_real, 0, sizeof(sizeof(double) * ns*ns));
    memset(ss_eq_rhs_imag, 0, sizeof(sizeof(double) * ns*ns));

        /* --- Precalculate lindblad operators --- */
    param->H0_func(param->H0, qs, param);

        for (i = 0; i < param->do_n; i++)
        {
                QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
                QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);

                param->do_aad[i] = NULL;
                param->do_ada[i] = NULL;

                if (param->do_g2[i] != 0.0)
                {
                        param->do_aad[i] = qdpack_operator_alloc(ns, ns);
                        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, param->do_a[i], param->do_a[i], beta, param->do_aad[i]);
                }
                if (param->do_g1[i] != 0.0)
                {
                        param->do_ada[i] = qdpack_operator_alloc(ns, ns);
                        qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, param->do_a[i], param->do_a[i], beta, param->do_ada[i]);
                }
        }
    
    /*
      * Add all operators to L:
      *
      */ 

    // L += -i H rho
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    steady_state_L_add_op(sop_L, param->H0, NULL, z, ns);

    // L +=  i rho H
    QDPACK_SET_COMPLEX(&z, 0.0, 1.0);
    steady_state_L_add_op(sop_L, NULL, param->H0, z, ns);

    // loop through all dissipation operators
    // drho_t += g(N+1)( a rho ad - 1/2ad a rho - 1/2rho ad a) + g N (ad rho a - 1/2a ad rho - 1/2rho a ad)
    for (i = 0; i < param->do_n; i++)
    {
        if (param->do_g1[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * param->do_g1[i], 0.0);
                steady_state_L_add_op(sop_L, NULL, param->do_ada[i], z, ns);
                steady_state_L_add_op(sop_L, param->do_ada[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * param->do_g1[i], 0.0);
                steady_state_L_add_op(sop_L, param->do_a[i], param->do_ad[i], z, ns);
        }
        if (param->do_g2[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * param->do_g2[i], 0.0);
                steady_state_L_add_op(sop_L, NULL, param->do_aad[i], z, ns);
                steady_state_L_add_op(sop_L, param->do_aad[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * param->do_g2[i], 0.0);
                steady_state_L_add_op(sop_L, param->do_ad[i], param->do_a[i], z, ns);
        }
    }

    /*
     * Replace first row with density matrix normalisation condition:
     *
     */
    QDPACK_SET_COMPLEX(&z, 1.0, 0.0);
    for (i = 0; i < ns; i++)
    {
        gsl_sparse_matrix_complex_set(sop_L, 0, ss_idx_major(ns, i,i), z);
    }
    // and create rhs: 
    //gsl_vector_complex_set_zero(ss_eq_rhs);
    //gsl_vector_complex_set(ss_eq_rhs, 0, z);
    ss_eq_rhs_real[0] = 1.0;    

    /*
      * Solve for rho_ss:
      *
      * L * rho_ss = [1, 0, 0, ...]'
      *
      */
    //gsl_sparse_matrix_complex_convert_to_dense(sop_L);
        
    //printf("Lm_real =\n");
        //qdpack_operator_print_real(sop_L->M_dense, ns*ns);
    //printf("Lm_imag =\n");
        //qdpack_operator_print_imag(sop_L->M_dense, ns*ns);

           //gsl_linalg_complex_LU_decomp(sop_L, p, &s);
           //gsl_linalg_complex_LU_solve (sop_L, p, ss_eq_rhs, rho_ss_vec);

    gsl_sparse_matrix_complex_convert_col(sop_L);

    gsl_sparse_matrix_complex_LU_solve(sop_L, ss_eq_rhs_real, ss_eq_rhs_imag, rho_ss_vec_real, rho_ss_vec_imag);

    //printf("rho_ss_vec = \n");
        //gsl_vector_complex_fprintf(stdout, rho_ss_vec, "%g");
    //for (i = 0; i < ns*ns; i++)
    //    printf("%f %f\n", rho_ss_vec_real[i], rho_ss_vec_imag[i]);
    //printf("ss_eq_rhs = \n");
    //for (i = 0; i < ns*ns; i++)
    //    printf("%f %f\n", ss_eq_rhs_real[i], ss_eq_rhs_imag[i]);
        //gsl_vector_complex_fprintf(stdout, ss_eq_rhs, "%g");


    /*
      * Copy rho_ss_vec to rho_ss:
      *
      */
    for (i = 0; i < ns*ns; i++)
    {
        //z = gsl_vector_complex_get(rho_ss_vec, i);
        QDPACK_SET_COMPLEX(&z, rho_ss_vec_real[i], rho_ss_vec_imag[i]);
        qdpack_operator_set(rho_ss, ss_idx_minor_ms(ns, i), ss_idx_minor_ls(ns, i), z);
    }

        //printf("rho_ss_real =\n");
        //qdpack_operator_print_real(rho_ss, ns);

    printf("DEBUG: solving for rho_ss: finished\n");
        
    return 0;
}




