//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Find the steady state of a quantum system with dissipation. 
 */

#include <stdio.h>


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>

#include <qdpack/qdpack.h>


#define QDPACK_TO_GSL_COMPLEX(z) (gsl_complex_rect(QDPACK_REAL(z),QDPACK_IMAG(z)))

/* ----------------------------------------------------------
 * Functions for calculating the Liouvillian.
 *
 */
static int
steady_state_L_add_op(qdpack_operator_t *sop_L,
                      qdpack_operator_t *op_left,
                      qdpack_operator_t *op_right,
                      qdpack_complex z, int N)
{
#ifdef GSL
    int a, b, c, d, i, j;
    qdpack_complex alpha, beta, gamma;

    //printf("DEBUG: steady_state_L_add_op: starting...");

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
                        i = ss_idx_major(N, a, b);
                        j = ss_idx_major(N, c, b);
    
                        // sop_L(I,J) = sop_L(I,J) + alpha
                        beta  = qdpack_operator_get(sop_L, i, j);
                        gamma = qdpack_complex_add(alpha, beta);
                        qdpack_operator_set(sop_L, i, j, gamma);
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
                        i = ss_idx_major(N, a, b);
                        j = ss_idx_major(N, a, d);
    
                        // sop_L(I,J) = sop_L(I,J) + alpha
                        beta  = qdpack_operator_get(sop_L, i, j);
                        gamma = qdpack_complex_add(alpha, beta);
                        qdpack_operator_set(sop_L, i, j, gamma);
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
                alpha = qdpack_operator_get(op_left, a, c);
                alpha = qdpack_complex_mul(alpha, z);
                
                if (QDPACK_REAL(alpha) != 0.0 || QDPACK_IMAG(alpha) != 0.0)
                {
                    for (b = 0; b < N; b++)
                    {
                        i = ss_idx_major(N, a, b);
                        
                        for (d = 0; d < N; d++)
                        {
                            j = ss_idx_major(N, c, d);

                            beta  = qdpack_operator_get(op_right, d, b);

                            if (QDPACK_REAL(beta) != 0.0 || QDPACK_IMAG(beta) != 0.0)
                            {
                                gamma = qdpack_complex_mul(alpha, beta);

                                // sop_L(I,J) = sop_L(I,J) + alpha
                                beta  = qdpack_operator_get(sop_L, i, j);
                                gamma = qdpack_complex_add(gamma, beta);
                                qdpack_operator_set(sop_L, i, j, gamma);
                            }
                        }
                    }
                }
            }
        }
    }
    
    //printf("done.\n");
#endif
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
qdpack_steadystate_dm(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss)
{
#ifdef GSL
    qdpack_operator_t *sop_L;
    gsl_vector_complex *rho_ss_vec, *ss_eq_rhs;
    gsl_permutation *p;

    qdpack_complex z, alpha, beta;
    int ns, s, i;

    qdpack_hilbert_space_t *qs_exp = qdpack_hilbert_space_new();

    ns = qdpack_hilbert_space_nstates(qs);

    qdpack_hilbert_space_add(qs_exp, ns*ns);


    printf("DEBUG: solving for rho_ss: starting building matrices\n");

    /*
      * Allocate space for the L super-operator and
       * some working space vectors for solving the
     * LA problem.
      */
    sop_L      = qdpack_operator_alloc(qs_exp);
    rho_ss_vec = gsl_vector_complex_alloc(ns*ns);
    ss_eq_rhs  = gsl_vector_complex_alloc(ns*ns);
    p        = gsl_permutation_alloc(ns*ns);

    qdpack_operator_set_zero(sop_L);

    /* --- Precalculate lindblad operators --- */
    sim->H0_func(sim, qs, sim->H0);

    for (i = 0; i < sim->do_n; i++)
    {
        QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
        QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);

        sim->do_aad[i] = NULL;
        sim->do_ada[i] = NULL;

        if (sim->do_g2[i] != 0.0)
        {
            sim->do_aad[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_aad[i]);
        }
        if (sim->do_g1[i] != 0.0)
        {
            sim->do_ada[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_ada[i]);
        }

        // debug
        //printf("a  =\n");
            //qdpack_operator_print_real(sim->do_a[i], qs->nstates[i]);
        //printf("ad =\n");
            //qdpack_operator_print_real(sim->do_ad[i], qs->nstates[i]);
        //printf("ada =\n");
            //qdpack_operator_print_real(sim->do_ada[i], qs->nstates[i]);
        //printf("aad =\n");
            //qdpack_operator_print_real(sim->do_aad[i], qs->nstates[i]);
    }    

    /*
      * Add all operators to L:
      *
      */ 

    // L += -i H rho
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    steady_state_L_add_op(sop_L, sim->H0, NULL, z, ns);

    // L +=  i rho H
    QDPACK_SET_COMPLEX(&z, 0.0, 1.0);
    steady_state_L_add_op(sop_L, NULL, sim->H0, z, ns);

    // loop through all dissipation operators
    // drho_t += g(N+1)( a rho ad - 1/2ad a rho - 1/2rho ad a) + g N (ad rho a - 1/2a ad rho - 1/2rho a ad)
    for (i = 0; i < sim->do_n; i++)
    {
        if (sim->do_g1[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * sim->do_g1[i], 0.0);
            steady_state_L_add_op(sop_L, NULL, sim->do_ada[i], z, ns);
            steady_state_L_add_op(sop_L, sim->do_ada[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * sim->do_g1[i], 0.0);
            steady_state_L_add_op(sop_L, sim->do_a[i], sim->do_ad[i], z, ns);
        }
        if (sim->do_g2[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * sim->do_g2[i], 0.0);
            steady_state_L_add_op(sop_L, NULL, sim->do_aad[i], z, ns);
            steady_state_L_add_op(sop_L, sim->do_aad[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * sim->do_g2[i], 0.0);
            steady_state_L_add_op(sop_L, sim->do_ad[i], sim->do_a[i], z, ns);
        }
    }

    /*
     * Replace first row with density matrix normalisation condition:
     *
     */
    QDPACK_SET_COMPLEX(&z, 1.0, 0.0);
    for (i = 0; i < ns; i++)
    {
        qdpack_operator_set(sop_L, 0, ss_idx_major(ns, i,i), z);
    }
    // and create rhs: 
    gsl_vector_complex_set_zero(ss_eq_rhs);
    gsl_vector_complex_set(ss_eq_rhs, 0, QDPACK_TO_GSL_COMPLEX(z));


    printf("DEBUG: solving for rho_ss: starting solving lin.eq. system\n");

    //
    // Solve for rho_ss:
    //
    // L * rho_ss = [1, 0, 0, ...]'
    //
    gsl_linalg_complex_LU_decomp(sop_L->data->data, p, &s);
    gsl_linalg_complex_LU_solve (sop_L->data->data, p, ss_eq_rhs, rho_ss_vec);

    //
    // Copy rho_ss_vec to rho_ss:
    //
    for (i = 0; i < ns*ns; i++)
    {
        z = gsl_vector_complex_get(rho_ss_vec, i);
        qdpack_operator_set(rho_ss, ss_idx_minor_ms(ns, i), ss_idx_minor_ls(ns, i), z);
    }

    printf("DEBUG: solving for rho_ss: finished\n");
        
#endif
    return 0;
}

/*
 * Find the steady-state density matrix for a quantum system with dissipation by
 * calculate the system Liouville super operators, and finding the eigenvector
 * with corresponding eigenvalue "0"
 */
int
qdpack_steadystate_dm_eig(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss)
{
    int ns, i = 0;
#ifdef GSL
    qdpack_operator_t *sop_L;
    gsl_vector_complex *rho_ss_vec, *ss_eq_rhs;

    qdpack_complex z, alpha, beta;

    qdpack_hilbert_space_t *qs_exp = qdpack_hilbert_space_new();

    ns = qdpack_hilbert_space_nstates(qs);

    qdpack_hilbert_space_add(qs_exp, ns*ns);

    printf("DEBUG: solving for rho_ss: starting building matrices\n");

    //
    // Allocate space for the L super-operator and
    // some working space vectors for solving the
    // LA problem.
    //
    sop_L      = qdpack_operator_alloc(qs_exp);
    rho_ss_vec = gsl_vector_complex_alloc(ns*ns);
    ss_eq_rhs  = gsl_vector_complex_alloc(ns*ns);

    qdpack_operator_set_zero(sop_L);

    /* --- Precalculate lindblad operators --- */
    sim->H0_func(sim, qs, sim->H0);

    for (i = 0; i < sim->do_n; i++)
    {
        QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
        QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);

        sim->do_aad[i] = NULL;
        sim->do_ada[i] = NULL;

        if (sim->do_g2[i] != 0.0)
        {
            sim->do_aad[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_aad[i]);
        }

        if (sim->do_g1[i] != 0.0)
        {
            sim->do_ada[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_ada[i]);
        }
    }    

    /*
     * Add all operators to L:
     *
     */ 

    // L += -i H rho
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    steady_state_L_add_op(sop_L, sim->H0, NULL, z, ns);

    // L +=  i rho H
    QDPACK_SET_COMPLEX(&z, 0.0, 1.0);
    steady_state_L_add_op(sop_L, NULL, sim->H0, z, ns);

    // loop through all dissipation operators
    // drho_t += g(N+1)( a rho ad - 1/2ad a rho - 1/2rho ad a) + g N (ad rho a - 1/2a ad rho - 1/2rho a ad)
    for (i = 0; i < sim->do_n; i++)
    {
        if (sim->do_g1[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * sim->do_g1[i], 0.0);
            steady_state_L_add_op(sop_L, NULL, sim->do_ada[i], z, ns);
            steady_state_L_add_op(sop_L, sim->do_ada[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * sim->do_g1[i], 0.0);
            steady_state_L_add_op(sop_L, sim->do_a[i], sim->do_ad[i], z, ns);
        }
        if (sim->do_g2[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&z, -0.5 * sim->do_g2[i], 0.0);
            steady_state_L_add_op(sop_L, NULL, sim->do_aad[i], z, ns);
            steady_state_L_add_op(sop_L, sim->do_aad[i], NULL, z, ns); // XXX: conjugate here...
            QDPACK_SET_COMPLEX(&z,  1.0 * sim->do_g2[i], 0.0);
            steady_state_L_add_op(sop_L, sim->do_ad[i], sim->do_a[i], z, ns);
        }
    }

    i = qdpack_steadystate_dm_propagator(sim, qs, sop_L, rho_ss);

    qdpack_operator_free(sop_L);
    gsl_vector_complex_free(rho_ss_vec);
    gsl_vector_complex_free(ss_eq_rhs);

    printf("DEBUG: solving for rho_ss: finished [%d]\n", i);
#endif        
    return i;
}

//
// Find the steady state given a (one-period) propagator for a peroidic system
//
int
qdpack_steadystate_dm_propagator(qdpack_simulation_t *sim, 
                                 qdpack_hilbert_space_t *qs, 
                                 qdpack_operator_t *U, 
                                 qdpack_operator_t *rho_ss)
{    
#ifdef GSL
    double dev_min = 1.0;
    int i, j = -1, qsn;
    qdpack_complex z_tr;

    gsl_matrix_complex *Uevec;
    gsl_vector_complex *Ueval;

    qsn = qdpack_hilbert_space_nstates(qs);

    Ueval = gsl_vector_complex_alloc(qsn*qsn);
    Uevec = gsl_matrix_complex_alloc(qsn*qsn, qsn*qsn);
    
    gsl_ext_eigen_zgeev(U->data->data, Uevec, Ueval); // xxx

    j = 0;
    for (i = 0; i < qsn*qsn; i++)
    {
        double dev;
        qdpack_complex z = gsl_vector_complex_get(Ueval, i);

        dev = qdpack_complex_abs(qdpack_complex_sub(z, QDPACK_COMPLEX_ONE));
        if (dev < dev_min)
        {
            dev_min = dev;
            j = i;
        }
    }

    if (dev_min < 0.0)
    {
        printf("%s: WARNING: dev_min < 0.0\n", __PRETTY_FUNCTION__);
    }

    if (j == -1)
    {
        printf("%s: failed to identify steady-state eigenvector.\n", __PRETTY_FUNCTION__);
        return -1;
    }

    for (i = 0; i < qsn*qsn; i++)
    {
        qdpack_operator_set(rho_ss, 
                            ss_idx_minor_ms(qsn, i), 
                            ss_idx_minor_ls(qsn, i),
                            gsl_matrix_complex_get(Uevec, i, j));
        //qdpack_complex_mul(qdpack_operator_get(Uevec, i, j), qdpack_operator_get(Uevec, i, j)));
    }

    // calculate the trace and rescale the density matrix:
    z_tr = operator_trace(rho_ss);
    z_tr = qdpack_complex_div(QDPACK_COMPLEX_ONE, z_tr);
    qdpack_operator_scale(rho_ss, z_tr);
 
    gsl_matrix_complex_free(Uevec);
    gsl_vector_complex_free(Ueval);
    
#endif
    return 0;
}



