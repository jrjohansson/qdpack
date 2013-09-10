//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file floquet.c
 *
 * @brief Functions for calculating floquet modes for (periodic) time-dependent
 * systems, and their evolution.
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <qdpack/qdpack.h>

/** ----------------------------------------------------------------------------
 * Floquet modes
 *
 */
int
qdpack_floquet_modes(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs)
{
    qdpack_operator_t *U; 
    qdpack_matrix_t *Uevec;
    qdpack_matrix_t *Ueval;
    
    int i, qsn = qdpack_hilbert_space_nstates(qs);

    U = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(U);

    sim->Ti = 0.0;
    sim->Tf = 2*M_PI/sim->h_td_w;
    sim->dT = 2*M_PI/sim->h_td_w;

    if (qdpack_propagator(sim, qs, U, NULL) == -1)
    {
        fprintf(stderr, "Evaluation of the propagator failed.\n");
        return -1;
    }

    Ueval = qdpack_matrix_alloc(qsn, 0);
    Uevec = qdpack_matrix_alloc(qsn, qsn);

    sim->floquet_quasienergies = qdpack_matrix_alloc(qsn, 0);
    sim->floquet_modes         = qdpack_operator_alloc(qs);

    qdpack_matrix_eigen_zgeev(U->data, Ueval, Uevec, QDPACK_EIGEN_SORT_PHASE);

    //
    // Uevec is the eigenvalues of the wave-function propagator, i.e.,
    // the Floquet modes.
    //    
    qdpack_matrix_memcpy(sim->floquet_modes->data, Uevec);

    //
    // calculate the quasi energies
    //
    for (i = 0; i < Ueval->m; i++)
    {   
        double arg = qdpack_complex_arg(qdpack_matrix_get(Ueval, i, 0));
        if (arg >   M_PI) arg -= 2*M_PI;        
        if (arg <= -M_PI) arg += 2*M_PI;    
                
        qdpack_matrix_set(sim->floquet_quasienergies, i, 0, qdpack_complex_rect(-arg/sim->Tf, 0.0));
    }

    qdpack_matrix_free(Uevec);
    qdpack_matrix_free(Ueval);
    qdpack_operator_free(U);
    
    return 0;
}

int
qdpack_floquet_modes_t(qdpack_simulation_t *sim,
                       qdpack_hilbert_space_t *qs, 
                       double t)
{
    qdpack_operator_t *U;

    int i, j;

    //
    // find the propagator for the time t
    //
    U = qdpack_operator_alloc(qs);

    sim->Ti = 0.0;
       sim->Tf = t;
    sim->dT = t;

    if (qdpack_propagator(sim, qs, U, NULL) == -1)
    {
        fprintf(stderr, "Evaluation of the propagator failed.\n");
        return -1;
    }

    //
    // Calculate the floquet mode at time t by propagating the mode at time 0:
    //
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, U, sim->floquet_modes, QDPACK_COMPLEX_ZERO, sim->floquet_states_t);

    for (i = 0; i < sim->floquet_modes_t->m; i++)
    {
        for (j = 0; j < sim->floquet_modes_t->n; j++)
        {
            qdpack_complex z1, z2;

            z1 = qdpack_operator_get(sim->floquet_states_t, i, j);

            z2 = qdpack_complex_mul(z1, qdpack_complex_polar(1.0, QDPACK_REAL(qdpack_matrix_get(sim->floquet_quasienergies, j, 0)) * t));
            qdpack_operator_set(sim->floquet_modes_t, i, j, z2);
        }
    }

    //
    // Uevec is the eigenvalues of the wave-function propagator, i.e.,
    // the Floquet states.
    //
    qdpack_operator_free(U);
    
    return 0;
}


/** ----------------------------------------------------------------------------
 * Calculates the dissipations rates to be used in a master equation in the 
 * floquet state basis. The rates are calculated based on the spectral density
 * callback function s_func_cb.
 *
 * @sim qs The quantum state structure.
 * @sim p  Pointer to solver simeter structure.
 * @sim s_func_cb Spectral density callback function.
 * @sim floquet_modes_store_cb Callback function for storing/processing the
 * floquet states.
 */
int
qdpack_floquet_dissipation_rates(qdpack_simulation_t *sim,
                                       qdpack_hilbert_space_t *qs, 
                                       operator_cb_func_t floquet_modes_store_cb,
                                       spectral_density_cb_t s_func_cb)
{
    qdpack_complex z;
    double t, delta_t, T;
    int i, i_limit = 100, j, qsn = qdpack_hilbert_space_nstates(qs);

    qdpack_operator_t *w1, *w2;

    int k, k_idx;

    //printf("%s: calculating the floquet dissipation rates... \n", __PRETTY_FUNCTION__);

    //
    // allocation
    //

    w1 = qdpack_operator_alloc(qs);
    w2 = qdpack_operator_alloc(qs);

    sim->floquet_k_max = 1;

    sim->floquet_gamma = malloc(sizeof(qdpack_operator_t *) * (2 * sim->floquet_k_max + 1));
    sim->floquet_delta = malloc(sizeof(qdpack_operator_t *) * (2 * sim->floquet_k_max + 1));
    sim->floquet_X     = malloc(sizeof(qdpack_operator_t *) * (2 * sim->floquet_k_max + 1));

    if (sim->floquet_X == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate floquet dissipation matrices.\n", __PRETTY_FUNCTION__);
        return -1;
    }

    for (k_idx = 0, k = - sim->floquet_k_max; k <= sim->floquet_k_max; k++, k_idx++)
    {
        sim->floquet_X[k_idx] = qdpack_operator_alloc(qs);
        qdpack_operator_set_zero(sim->floquet_X[k_idx]);

        sim->floquet_delta[k_idx] = qdpack_operator_alloc(qs);
        qdpack_operator_set_zero(sim->floquet_delta[k_idx]);

        sim->floquet_gamma[k_idx] = qdpack_operator_alloc(qs);
        qdpack_operator_set_zero(sim->floquet_gamma[k_idx]);
    }

    sim->floquet_modes_t  = qdpack_operator_alloc(qs);
    sim->floquet_states_t = qdpack_operator_alloc(qs);
    sim->floquet_A        = qdpack_operator_alloc(qs);

    if (sim->floquet_Gamma == NULL)
    {
        sim->floquet_Gamma = qdpack_matrix_alloc(qsn, qsn);
    }

    //
    //
    //
    floquet_modes_store_cb(sim, qs, sim->floquet_modes, 0.0);


    //
    // Evaluate the integral: 
    //
    // X_ijk = (1/T) int_0^T dt exp(-i k Omega t) <Phi_i(t) | \hat{X} | \Phi_j(t) >
    //
    //Omega   = sim->h_td_w / (2*M_PI) ;
    T       = 2*M_PI/sim->h_td_w;
    delta_t = T / i_limit;
    for (t = delta_t; t <= T+delta_t/2.0; t += delta_t)
    {
        // get the floquet states at time t
        qdpack_floquet_modes_t(sim, qs, t);

        // store/process the floquet states
        floquet_modes_store_cb(sim, qs, sim->floquet_modes_t, t);

        // loop though the k values -k_max ... k_max
        for (k_idx = 0, k = - sim->floquet_k_max; k <= sim->floquet_k_max; k++, k_idx++)
        {
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->floquet_modes_t, sim->floquet_operator, QDPACK_COMPLEX_ZERO, w1);
            qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, w1                  , sim->floquet_modes_t,  QDPACK_COMPLEX_ZERO, w2);

            z = qdpack_complex_mul_real(qdpack_complex_polar(1.0, - 2 * M_PI * k * t / T), delta_t / T);  

            qdpack_operator_scale(w2, z);

            qdpack_operator_add(sim->floquet_X[k_idx], w2);
        }   
    }

    //
    // calculate delta and gamma matrices. 
    //
    // XXX: these matrices could be real-valued.
    //
    for (k_idx = 0, k = - sim->floquet_k_max; k <= sim->floquet_k_max; k++, k_idx++)
    {
        printf("%s: floquet_X[%d]\n", __PRETTY_FUNCTION__, k_idx);
        qdpack_operator_print_real(sim->floquet_X[k_idx]);
        qdpack_operator_print_imag(sim->floquet_X[k_idx]);

        // calculate delta and gamma
        for (i = 0; i < qsn; i++)
        {
            for (j = 0; j < qsn; j++)
            {
                qdpack_complex z;
                double delta, gamma;

                // delta
                delta = QDPACK_REAL(qdpack_matrix_get(sim->floquet_quasienergies, i, 0))
                      - QDPACK_REAL(qdpack_matrix_get(sim->floquet_quasienergies, j, 0))
                      + 2 * M_PI * k / T;
                QDPACK_SET_COMPLEX(&z, delta, 0);
                qdpack_operator_set(sim->floquet_delta[k_idx], i, j, z);

                // gamma
                gamma = 2 * M_PI * heaviside_function(delta) * s_func_cb(sim, delta) 
                                 * qdpack_complex_abs2(qdpack_operator_get(sim->floquet_X[k_idx], i, j));
                QDPACK_SET_COMPLEX(&z, gamma, 0);
                qdpack_operator_set(sim->floquet_gamma[k_idx], i, j, z);
            }            
        }

        printf("%s: floquet_delta[%d]\n", __PRETTY_FUNCTION__, k_idx);
        qdpack_operator_print_real(sim->floquet_delta[k_idx]);
        printf("%s: floquet_gamma[%d]\n", __PRETTY_FUNCTION__, k_idx);
        qdpack_operator_print_real(sim->floquet_gamma[k_idx]);
    }


    // calculate A
    // todo: write as sum over k only, using delta and gamma in matrix form ?
    for (i = 0; i < qsn; i++)
    {
        for (j = 0; j < qsn; j++)
        {
            double a_ij_sum = 0.0;

            for (k = - sim->floquet_k_max; k <= sim->floquet_k_max; k++)
            {
                double g1, g2, d;
                int k1_idx, k2_idx;

                k1_idx =   k + sim->floquet_k_max;
                k2_idx = - k + sim->floquet_k_max;

                d  = fabs(QDPACK_REAL(qdpack_operator_get(sim->floquet_delta[k1_idx], i, j)));

                g1 = QDPACK_REAL(qdpack_operator_get(sim->floquet_gamma[k1_idx], i, j));
                g2 = QDPACK_REAL(qdpack_operator_get(sim->floquet_gamma[k2_idx], j, i));

                a_ij_sum += g1 + boltzmann_factor(d, sim->T_w) * (g1 + g2);

            }

            QDPACK_SET_COMPLEX(&z, a_ij_sum, 0);
            qdpack_operator_set(sim->floquet_A, i, j, z);
        }
    }

    printf("%s: floquet_A\n", __PRETTY_FUNCTION__);
    qdpack_operator_print_real(sim->floquet_A);
    
    // XXX: this is specific for 2-level systems...
    qdpack_matrix_set(sim->floquet_Gamma, 0, 1, 
                      qdpack_complex_rect(QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 0, 1)), 0.0));
    qdpack_matrix_set(sim->floquet_Gamma, 1, 0, 
                      qdpack_complex_rect(QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 1, 0)), 0.0));
    qdpack_matrix_set(sim->floquet_Gamma, 0, 0,
                      qdpack_complex_rect(
                        (QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 0, 0)) +
                         QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 0, 1)) +
                         QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 1, 0)) +
                         QDPACK_REAL(qdpack_operator_get(sim->floquet_A, 1, 1))) / 2.0, 0.0));
   
    qdpack_matrix_set(sim->floquet_Gamma, 1, 1, qdpack_matrix_get(sim->floquet_Gamma, 0, 0));

    if (sim->floquet_ss_prob == NULL)
    {   
        sim->floquet_ss_prob = qdpack_matrix_alloc(qsn, 0);
    }

    // XXX: specific to 2LS...
    // XXX: not working...
    qdpack_matrix_set(sim->floquet_ss_prob, 0, 0,
                      qdpack_complex_rect(
                          QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 1, 0)) / 
                          (QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 0, 1)) + 
                           QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 1, 0))), 0.0));
    qdpack_matrix_set(sim->floquet_ss_prob, 1, 0,
                      qdpack_complex_rect(
                          QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 0, 1)) / 
                          (QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 0, 1)) + 
                           QDPACK_REAL(qdpack_matrix_get(sim->floquet_Gamma, 1, 0))), 0.0));


    //
    // free
    //
    for (k_idx = 0, k = - sim->floquet_k_max; k < sim->floquet_k_max; k++, k_idx++)
    {
        //printf("%s: floquet_X[%d]\n", __PRETTY_FUNCTION__, k_idx);
        //qdpack_operator_print_real(sim->floquet_X[k_idx]);
        //qdpack_operator_print_imag(sim->floquet_X[k_idx]);

        //printf("%s: floquet_delta[%d]\n", __PRETTY_FUNCTION__, k_idx);
        //qdpack_operator_print_real(sim->floquet_delta[k_idx]);

        //printf("%s: floquet_gamma[%d]\n", __PRETTY_FUNCTION__, k_idx);
        //qdpack_operator_print_real(sim->floquet_gamma[k_idx]);

        qdpack_operator_free(sim->floquet_X[k_idx]);
        qdpack_operator_free(sim->floquet_delta[k_idx]);
        qdpack_operator_free(sim->floquet_gamma[k_idx]);
    }

    qdpack_operator_free(w1);
    qdpack_operator_free(w2);
    qdpack_operator_free(sim->floquet_A);

    //qdpack_operator_free(sim->floquet_Gamma);
    //qdpack_matrix_free(sim->floquet_ss_prob);

    free(sim->floquet_X);
    free(sim->floquet_delta);
    free(sim->floquet_gamma);
    
    return 0;
}

