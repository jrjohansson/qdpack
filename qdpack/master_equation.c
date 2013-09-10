//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file master_equation.c
 *
 * @brief Solvers for the time-evolution of Schrödinger equation and various
 * master equations.
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <qdpack/qdpack.h>


//==============================================================================
// WAVE FUNCTION EVOLUTION
//==============================================================================


/** ---------------------------------------------------------------------------- 
 * Calculate the derivative of the wave function 
 *
 */
static int
he_func_dpsi_dt(qdpack_simulation_t *sim,
                qdpack_hilbert_space_t *qs,
                qdpack_state_t *psi_t,
                qdpack_state_t *dpsi_dt, 
                double t)
{
    qdpack_complex z;

    if (sim->Ht_func)
    {
        sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, t);
    }
    else
    {
        sim->H_t = sim->H0;
    }
    
    // equation of motion: dpsi_dt = -i * H * psi
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->H_t, psi_t, QDPACK_COMPLEX_ZERO, dpsi_dt);
    qdpack_operator_state_multiply(sim->H_t, psi_t, dpsi_dt);
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    qdpack_state_scale(dpsi_dt, z);

    return 0;
}

int
qdpack_evolve_wf(qdpack_simulation_t *sim,
            qdpack_hilbert_space_t *qs, 
            qdpack_state_t *psi0,
            state_cb_func_t psi_cb_func)
{
    //printf("%s: starting evolution (Schrodinger equation)\n", __PRETTY_FUNCTION__);

    /** -- setup hamiltonians --- */
    if (sim->H0_func)
    {
        if (sim->H0 == NULL)
        {
            sim->H0 = qdpack_operator_alloc(qs);
        }
        sim->H0_func(sim, qs, sim->H0);
    }
    else
    {   
        fprintf(stderr, "%s: No hamiltonian.", __PRETTY_FUNCTION__);
        return -1;    
    }

    if (sim->H1_func)
    {
        if (sim->H1 == NULL)
        {
            sim->H1 = qdpack_operator_alloc(qs);
        }
        sim->H1_func(sim, qs, sim->H1);
    }

    if (sim->Ht_func)
    {
        if (sim->H_t == NULL)
        {
            sim->H_t = qdpack_operator_alloc(qs);
        }
    }

    /** -- Evolve the system with an ODE solver --- */
    odeint_state_rk_adaptive(sim, qs, psi0, he_func_dpsi_dt, psi_cb_func);
    
    return 0;
}



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// DENSITY MATRIX EVOLUTION
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------


/** 
 *
 * Evolve the density matrix rho0 assuming a non-dissipative dynamics
 * and a constant hamiltonian, using a very simple method. 
 * sim->Ti = start time
 * sim->Tf = end time
 *
 * XXX: not working: need to calculate exp_H
 */
int
qdpack_evolve_dm_unitary_const_simple(qdpack_simulation_t *sim,
                               qdpack_hilbert_space_t *qs,
                               qdpack_operator_t *rho0,
                               operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *rho_t, *ws1;
    double t;
    
    rho_t = qdpack_operator_alloc(qs);
    ws1   = qdpack_operator_alloc(qs);

    /** -- setup hamiltonians --- */
    if (sim->H0_func)
    {
        if (sim->H0 == NULL)
        {
            sim->H0 = qdpack_operator_alloc(qs);
        }
        sim->H0_func(sim, qs, sim->H0);
    }
    else
    {   
        fprintf(stderr, "%s: No hamiltonian.", __PRETTY_FUNCTION__);
        return -1;    
    }

    qdpack_operator_memcpy(rho_t, rho0);
    rho_cb_func(sim, qs, rho0, sim->Ti);

    for (t = sim->Ti; t < sim->Tf; t += sim->dT)
    {
        //qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, exp_H, rho_t, QDPACK_COMPLEX_ZERO, ws1);
        //qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, ws1,   exp_H, QDPACK_COMPLEX_ZERO, rho_t);
        
        rho_cb_func(sim, qs, rho_t, t);
    }

    qdpack_operator_free(rho_t);
    qdpack_operator_free(ws1);
    
    return 0;
}

/**
 * Evolve the density matrix rho0 assuming a non-dissipative dynamics
 * and a constant hamiltonian. 
 * sim->Ti = start time
 * sim->Tf = end time
 *
 *
 */
int
qdpack_evolve_dm_unitary_const(qdpack_simulation_t *sim,
                        qdpack_hilbert_space_t *qs,
                        qdpack_operator_t *rho0,
                        operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *eHr, *eHl, *rho_t, *H;
    double t;
    qdpack_complex alpha, beta, z;
    
    eHl = qdpack_operator_alloc(qs);
    eHr = qdpack_operator_alloc(qs);
    rho_t = qdpack_operator_alloc(qs);

    H = qdpack_operator_alloc(qs);
    sim->H0_func(sim, qs, H);
    
    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
        
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0 * sim->dT);
    qdpack_operator_scale(H, z);

    qdpack_matrix_exp(H->data, eHl->data);
    
    QDPACK_SET_COMPLEX(&z, -1.0, 0.0);
    qdpack_operator_scale(H, z);

    qdpack_matrix_exp(H->data, eHr->data);

    qdpack_operator_memcpy(rho_t, rho0);
    rho_cb_func(sim, qs, rho0, sim->Ti);

    for (t = sim->Ti; t < sim->Tf; t += sim->dT)
    {
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, eHr, beta, H);
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, eHl, H, beta, rho_t);
        
        rho_cb_func(sim, qs, rho_t, t);
    }

    qdpack_operator_free(eHl);
    qdpack_operator_free(eHr);
    qdpack_operator_free(rho_t);
    qdpack_operator_free(H);
    
    return 0;
}

/**
 * Evolve the density matrix rho0 assuming a non-dissipative dynamics
 * and a constant hamiltonian. 
 * sim->Ti = start time
 * sim->Tf = end time
 */
int
qdpack_evolve_dm_unitary_const_full_steps(qdpack_simulation_t *sim,
                                   qdpack_hilbert_space_t *qs,
                                   qdpack_operator_t *rho0,
                                   operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *eHr, *eHl, *rho_t, *H, *R;
    double t;
    qdpack_complex alpha, beta, z;

    eHl = qdpack_operator_alloc(qs);
    eHr = qdpack_operator_alloc(qs);
    R   = qdpack_operator_alloc(qs);
    rho_t = qdpack_operator_alloc(qs);

    H = qdpack_operator_alloc(qs);
    sim->H0_func(sim, qs, H);
    
    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
        
    rho_cb_func(sim, qs, rho0, sim->Ti);

    for (t = sim->Ti; t < sim->Tf; t += sim->dT)
    {
        qdpack_operator_memcpy(R, H);
        
        QDPACK_SET_COMPLEX(&z, 0.0, -1.0 * (t - sim->Ti));
        qdpack_operator_scale(R, z);
        qdpack_matrix_exp(R->data, eHl->data);
    
        QDPACK_SET_COMPLEX(&z, -1.0, 0.0);
        qdpack_operator_scale(R, z);
        qdpack_matrix_exp(R->data, eHr->data);

        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho0, eHr, beta, R);
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, eHl,  R,   beta, rho_t);
        
        rho_cb_func(sim, qs, rho_t, t);
    }

    qdpack_operator_free(H);
    qdpack_operator_free(R);
    qdpack_operator_free(eHl);
    qdpack_operator_free(eHr);
    return 0;
}



/** ----------------------------------------------------------------------------
 * Evolve the density matrix rho0 assuming a non-dissipative dynamics and
 * time-dependent hamiltonian. 
 *
 * sim->Ti = start time
 * sim->Tf = end time
 *
 *
 */
int
unitary_func_drho_dt(qdpack_simulation_t *sim,
                     qdpack_hilbert_space_t *qs,
                     qdpack_operator_t *rho_t, 
                     qdpack_operator_t *drho_dt, 
                     double t)
                     
{
    qdpack_complex z;

    if (sim->Ht_func)
    {
        sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, t);
    }
    else
    {
        sim->H_t = sim->H0;
    }

    // equation of motion:  drho_dt = -i ( H * rho - rho * H )
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, rho_t, sim->H_t, QDPACK_COMPLEX_ZERO, drho_dt);
    QDPACK_SET_COMPLEX(&z, -1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->H_t, rho_t, z, drho_dt);
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    qdpack_operator_scale(drho_dt, z);
    
    return 0;
}

int
qdpack_evolve_dm_unitary(qdpack_simulation_t *sim,
                    qdpack_hilbert_space_t *qs, 
                    qdpack_operator_t *rho0,
                    operator_cb_func_t rho_cb_func)
{
    printf("%s: starting unitary evolution of density matrix\n", __PRETTY_FUNCTION__);

    /** -- setup hamiltonians --- */
    if (sim->H0_func)
    {
        if (sim->H0 == NULL)
        {
            sim->H0 = qdpack_operator_alloc(qs);
        }
        sim->H0_func(sim, qs, sim->H0);
    }
    else
    {   
        fprintf(stderr, "%s: No hamiltonian.", __PRETTY_FUNCTION__);
        return -1;    
    }

    if (sim->H1_func)
    {
        if (sim->H1 == NULL)
        {
            sim->H1 = qdpack_operator_alloc(qs);
        }
        sim->H1_func(sim, qs, sim->H1);
    }

    if (sim->Ht_func)
    {
        if (sim->H_t == NULL)
        {
            sim->H_t = qdpack_operator_alloc(qs);
        }
    }

    //odeint_complex_rk(sim, qs, rho0, unitary_func_drho_dt, rho_cb_func);
    odeint_operator_rk_adaptive(sim, qs, rho0, unitary_func_drho_dt, rho_cb_func);

    return 0;
}

/** ----------------------------------------------------------------------------
 * Evolve the density matrix rho0 using a dissipative dynamics and time-dependent
 * hamiltonian. 
 * sim->Ti = start time
 * sim->Tf = end time
 *
 *
 */
int
lme_func_drho_dt(qdpack_simulation_t *sim,
                 qdpack_hilbert_space_t *qs,
                 qdpack_operator_t *rho_t,
                 qdpack_operator_t *drho_dt,
                 double t)
{
    int i;
    qdpack_complex z, alpha, beta;

    if (sim->Ht_func)
    {
        sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, t);
    }
    else
    {
        sim->H_t = sim->H0;
    }

    // drho_dt = -i ( H * rho - rho * H )
    QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, sim->H_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&beta,  -1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->H_t, rho_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    qdpack_operator_scale(drho_dt, z);

    // optional continuous basis transformation
    // drho_dt = ( dS/dt S^* * rho - rho * dS/dt S^* )
    if (sim->cont_basis_transform == 1)
    {
        qdpack_operator_t *S, *dS;
        
        printf("DEBUG: using continous basis transformation of master eqaution: %f\n", t);
        S  = basis_transform(sim, qs, 0, sim->ho_w_func, t);
        dS = basis_transform_derivative(sim, qs, 0, sim->ho_w_func, t);
        
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, QDPACK_COMPLEX_ONE, dS, S, QDPACK_COMPLEX_ZERO, sim->dSSt);
        
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE,    rho_t, sim->dSSt, QDPACK_COMPLEX_ONE, drho_dt);
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_NEGONE, sim->dSSt, rho_t, QDPACK_COMPLEX_ONE, drho_dt);

        qdpack_operator_free(S);
        qdpack_operator_free(dS);
    }
    
    // loop through all dissipation operators
    // drho_t += g/2 (N+1)( 2 a rho ad - ad a rho - rho ad a) + g/2 N ( 2 ad rho a - a ad rho - rho a ad)

    QDPACK_SET_COMPLEX(&alpha, -2.0, 0.0);

    for (i = 0; i < sim->do_n; i++)
    {
        if (sim->do_g1[i] != 0.0)
        {
            //printf("DEBUG: LME: operator %d, a. (%f, %f)\n", i, sim->do_g1[i], sim->do_g2[i]);
            //qdpack_operator_set_zero(sim->ws1);
            //qdpack_operator_set_zero(sim->ws2);
            
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, rho_t, sim->do_ada[i], QDPACK_COMPLEX_ZERO, sim->ws1);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->do_ada[i], rho_t, QDPACK_COMPLEX_ONE, sim->ws1);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->do_a[i], rho_t, QDPACK_COMPLEX_ZERO, sim->ws2);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->ws2, sim->do_a[i], QDPACK_COMPLEX_ONE, sim->ws1);
            
            QDPACK_SET_COMPLEX(&z, -sim->do_g1[i] / 2, 0.0);
            qdpack_operator_scale(sim->ws1, z);
            
            qdpack_operator_add(drho_dt, sim->ws1);
        }

        if (sim->do_g2[i] != 0.0)
        {
            //printf("DEBUG: LME: operator %d, ad. (%f, %f)\n", i, sim->do_g1[i], sim->do_g2[i]);
            //qdpack_operator_set_zero(sim->ws1);
            //qdpack_operator_set_zero(sim->ws2);
            
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, rho_t, sim->do_aad[i], QDPACK_COMPLEX_ZERO, sim->ws1);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->do_aad[i], rho_t, QDPACK_COMPLEX_ONE, sim->ws1);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->do_a[i], rho_t, QDPACK_COMPLEX_ZERO, sim->ws2);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->ws2, sim->do_a[i], QDPACK_COMPLEX_ONE, sim->ws1);
            
            QDPACK_SET_COMPLEX(&z, -sim->do_g2[i]/2, 0.0);
            qdpack_operator_scale(sim->ws1, z);
        
            qdpack_operator_add(drho_dt, sim->ws1);
        }        
    }

    return 0;
}

/** ----------------------------------------------------------------------------
 * Evolve the density matrix rho0 using a dissipative dynamics and time-dependent
 * hamiltonian, using the instantaneous basis for dissipation operators. 
 * 
 * sim->Ti = start time
 * sim->Tf = end time
 *
 */
int
lme_eb_func_drho_dt(qdpack_simulation_t *sim,
                    qdpack_hilbert_space_t *qs,
                    qdpack_operator_t *rho_t, 
                    qdpack_operator_t *drho_dt, 
                    double t)
{
    int i;
    qdpack_complex z, alpha, beta;

    if (sim->Ht_func)
    {
        sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, t);
    }
    else
    {
        sim->H_t = sim->H0;
    }

    printf("H(%.2f) =\n", t);
    qdpack_operator_print(sim->H_t);

    // drho_dt = -i ( H * rho - rho * H )
    QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, sim->H_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&beta,  -1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->H_t, rho_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    qdpack_operator_scale(drho_dt, z);

    //
    // recalculate the dissipation operators... transform to the instantaneous
    // eigenbasis

    qdpack_matrix_eigen_hermv(sim->H_t->data, sim->eval, sim->evec->data, QDPACK_EIGEN_SORT_VAL_ASC);

    for (i = 0; i < sim->do_n; i++)
    {    
        QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
        QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);

        // do_a_eb[i] = S do_a[i] S^{-1}
        qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->evec,   sim->do_a[i], QDPACK_COMPLEX_ZERO, sim->op_tmp);
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, QDPACK_COMPLEX_ONE, sim->op_tmp, sim->evec,    QDPACK_COMPLEX_ZERO, sim->do_a_eb[i]);

        if (sim->do_g2[i] != 0.0)
        {
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans,   QDPACK_COMPLEX_ONE, sim->evec,   sim->do_aad[i], QDPACK_COMPLEX_ZERO, sim->op_tmp);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, QDPACK_COMPLEX_ONE, sim->op_tmp, sim->evec,     QDPACK_COMPLEX_ZERO, sim->do_aad_eb[i]);
        }
        if (sim->do_g1[i] != 0.0)
        {
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans,   QDPACK_COMPLEX_ONE, sim->evec,   sim->do_ada[i], QDPACK_COMPLEX_ZERO, sim->op_tmp);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, QDPACK_COMPLEX_ONE, sim->op_tmp, sim->evec,     QDPACK_COMPLEX_ZERO, sim->do_ada_eb[i]);
        }
    }

    // optional continuous basis transformation
    // drho_dt = ( dS/dt S^* * rho - rho * dS/dt S^* )
    if (sim->cont_basis_transform == 1)
    {
        qdpack_operator_t *S, *dS;
        
        S  = basis_transform(sim, qs, 0, sim->ho_w_func, t);
        dS = basis_transform_derivative(sim, qs, 0, sim->ho_w_func, t);
        
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, QDPACK_COMPLEX_ONE, dS, S, QDPACK_COMPLEX_ZERO, sim->dSSt);
        
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE,    rho_t, sim->dSSt, QDPACK_COMPLEX_ONE, drho_dt);
        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_NEGONE, sim->dSSt, rho_t, QDPACK_COMPLEX_ONE, drho_dt);

        qdpack_operator_free(S);
        qdpack_operator_free(dS);
    }
    
    // loop through all dissipation operators
    // drho_t += g/2 (N+1)( 2 a rho ad - ad a rho - rho ad a) + g/2 N ( 2 ad rho a - a ad rho - rho a ad)
    for (i = 0; i < sim->do_n; i++)
    {
        if (sim->do_g1[i] != 0.0)
        {
            qdpack_operator_set_zero(sim->ws1);
            qdpack_operator_set_zero(sim->ws2);
            
            QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
            
            QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, sim->do_ada_eb[i], beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->do_ada_eb[i], rho_t, beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->do_a_eb[i], rho_t, beta, sim->ws2);
            QDPACK_SET_COMPLEX(&alpha, -2.0, 0.0);
            QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->ws2, sim->do_a_eb[i], beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&z, -sim->do_g1[i] / 2, 0.0);
            qdpack_operator_scale(sim->ws1, z);
            
            qdpack_operator_add(drho_dt, sim->ws1);
        }

        if (sim->do_g2[i] != 0.0)
        {
            qdpack_operator_set_zero(sim->ws1);
            qdpack_operator_set_zero(sim->ws2);
            
            QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
            
            QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, sim->do_aad_eb[i], beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->do_aad_eb[i], rho_t, beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, sim->do_a_eb[i], rho_t, beta, sim->ws2);
            QDPACK_SET_COMPLEX(&alpha, -2.0, 0.0);
            QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->ws2, sim->do_a_eb[i], beta, sim->ws1);
            
            QDPACK_SET_COMPLEX(&z, -sim->do_g2[i]/2, 0.0);
            qdpack_operator_scale(sim->ws1, z);
        
            qdpack_operator_add(drho_dt, sim->ws1);
        }        
    }

    return 0;
}

int
qdpack_evolve_dm_lme(qdpack_simulation_t *sim,
                     qdpack_hilbert_space_t *qs, 
                     qdpack_operator_t *rho0,
                     operator_cb_func_t rho_cb_func)
{
    qdpack_complex alpha, beta;
    int i;

    // work space matrices
    sim->ws1 = qdpack_operator_alloc(qs);
    sim->ws2 = qdpack_operator_alloc(qs);
    sim->op_tmp = qdpack_operator_alloc(qs);

    /** -- setup hamiltonians --- */
    if (sim->H0_func)
    {
        if (sim->H0 == NULL)
        {
            sim->H0 = qdpack_operator_alloc(qs);
        }
        sim->H0_func(sim, qs, sim->H0);
    }
    else
    {   
        fprintf(stderr, "%s: No hamiltonian.", __PRETTY_FUNCTION__);
        return -1;    
    }

    if (sim->H1_func)
    {
        if (sim->H1 == NULL)
        {
            sim->H1 = qdpack_operator_alloc(qs);
        }
        sim->H1_func(sim, qs, sim->H1);
    }

    if (sim->Ht_func)
    {
        if (sim->H_t == NULL)
        {
            sim->H_t = qdpack_operator_alloc(qs);
        }
    }

    if (sim->rho_f == NULL)
    {
        sim->rho_f = qdpack_operator_alloc(qs);
    }

    if (sim->cont_basis_transform == 1)
    {
        sim->dSSt = qdpack_operator_alloc(qs);
    }

    /** -- Precalculate lindblad operators --- */
    for (i = 0; i < sim->do_n; i++)
    {    
        QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
        QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);

        if (sim->option_me_eb)
        {
            //sim->do_ad_eb[i] = qdpack_operator_alloc(qs);
            sim->do_a_eb[i]  = qdpack_operator_alloc(qs);
        }

        sim->do_aad[i] = NULL;
        sim->do_ada[i] = NULL;

        //if (sim->do_g2[i] != 0.0)
        {
            sim->do_aad[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_aad[i]);

            if (sim->option_me_eb)
                sim->do_aad_eb[i] = qdpack_operator_alloc(qs);
        }

        //if (sim->do_g1[i] != 0.0)
        {
            sim->do_ada[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, sim->do_a[i], sim->do_a[i], beta, sim->do_ada[i]);

            if (sim->option_me_eb)
                sim->do_ada_eb[i] = qdpack_operator_alloc(qs);
        }
    }
    

    /** -- Evolve the system with an ODE solver --- */
    if (sim->option_adaptive)
    {
        if (sim->option_me_eb)
            odeint_operator_rk_adaptive(sim, qs, rho0, lme_eb_func_drho_dt, rho_cb_func);
        else
            odeint_operator_rk_adaptive(sim, qs, rho0, lme_func_drho_dt,    rho_cb_func);
    }
    else
    {
        if (sim->option_me_eb)
            odeint_operator_rk(sim, qs, rho0, lme_eb_func_drho_dt, rho_cb_func);
        else
            odeint_operator_rk(sim, qs, rho0, lme_func_drho_dt,    rho_cb_func);
    }
    
    /** -- free matrices --- */
    for (i = 0; i < sim->do_n; i++)
    {
        if (sim->option_me_eb)
        {
            //qdpack_operator_free(sim->do_ad_eb[i]);
            qdpack_operator_free(sim->do_a_eb[i]);
        }

        if (sim->do_aad[i] != NULL)
        {
            qdpack_operator_free(sim->do_aad[i]);
            if (sim->option_me_eb)
                qdpack_operator_free(sim->do_aad_eb[i]);
        }

        if (sim->do_ada[i] != NULL)
        {
            qdpack_operator_free(sim->do_ada[i]);
            if (sim->option_me_eb)    
                qdpack_operator_free(sim->do_ada_eb[i]);
        }
    }
    
    if (sim->cont_basis_transform == 1)
    {
        qdpack_operator_free(sim->dSSt);
    }

    qdpack_operator_free(sim->ws1);
    qdpack_operator_free(sim->ws2);
    qdpack_operator_free(sim->op_tmp);

    return 0;
}


/** ---------------------------------------------------------------------------- 
 * Evolve the density matrix rho0 using non-dissipative dynamics and
 * time-dependent hamiltonian. 
 *
 * sim->Ti = start time
 * sim->Tf = end time
 *
 */
int
he_func_drho_dt(qdpack_simulation_t *sim, 
                qdpack_hilbert_space_t *qs,
                qdpack_operator_t *rho_t, 
                qdpack_operator_t *drho_dt, double t)
{
    qdpack_complex z, alpha, beta;

    if (sim->Ht_func)
    {
        sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, t);
    }
    else
    {
        sim->H_t = sim->H0;
    }
    
    // equation of motion: drho_dt = -i ( H * rho - rho * H )
    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, rho_t, sim->H_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&beta, -1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, sim->H_t, rho_t, beta, drho_dt);
    QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
    qdpack_operator_scale(drho_dt, z);

    return 0;
}

int
qdpack_evolve_dm_he(qdpack_simulation_t *sim,
               qdpack_hilbert_space_t *qs,
               qdpack_operator_t *rho0,
               operator_cb_func_t rho_cb_func)
{
    /** -- calculate hamiltonians --- */
    
    sim->H0_func(sim, qs, sim->H0);
    if (sim->H1_func)
    {
        sim->H1_func(sim, qs, sim->H1);
    }
    else
    {
        sim->H1 = NULL;
    }

    /** -- Evolve the system with an ODE solver --- */
    odeint_operator_rk_adaptive(sim, qs, rho0, he_func_drho_dt, rho_cb_func);
    
    return 0;
}

//==============================================================================
// PROGAGATORS
//==============================================================================

/** ----------------------------------------------------------------------------
 * Calculate the evolution operator by evolving the density matrix according to
 * the master equation.
 */
int
qdpack_propagator_dm(qdpack_simulation_t *sim,
                     qdpack_hilbert_space_t *qs, 
                     qdpack_operator_t *U,
                     operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *rho0;
    
    int a, b, i, j, qsn = qdpack_hilbert_space_nstates(qs);

    rho0 = qdpack_operator_alloc(qs);
    qdpack_operator_set_zero(U);
    
    for (i = 0; i < qsn*qsn; i++)
    {
        //
        // setup (non-physical) initial density matrix
        //
        qdpack_operator_set_zero(rho0);
        qdpack_operator_set(rho0, ss_idx_minor_ms(qsn, i), ss_idx_minor_ls(qsn, i), QDPACK_COMPLEX_ONE);

        //
        // evolve
        //

        if (qdpack_evolve_dm_lme(sim, qs, rho0, rho_cb_func) != 0)
        //if (qdpack_evolve_dm_lme_ib(qs, rho0, sim, rho_cb_func) != 0) // relaxation acting in the eigenbasis.
        {
            fprintf(stderr, "%s: Evolution of the quantum system failed.\n", __PRETTY_FUNCTION__);
            return -1;
        }

        //
        // copy the resulting density matrix into U
        //        
        for (a = 0; a < qsn; a++)
        {
            for (b = 0; b < qsn; b++)
            {
                j = ss_idx_major(qsn, a, b);    
                qdpack_operator_set(U, j, i, qdpack_operator_get(sim->rho_f, a, b));    
            }
        }
    }    
    
    qdpack_operator_free(rho0);

    return 0;
}


/** ----------------------------------------------------------------------------
 * Calculate the evolution operator by evolving the wave function according to
 * the Hamiltonian.
 *
 */
int
qdpack_propagator(qdpack_simulation_t *sim,
                  qdpack_hilbert_space_t *qs,
                  qdpack_operator_t *U,
                  state_cb_func_t psi_cb_func)
{
    qdpack_state_t *psi0;
    
    int i, j, qsn = qdpack_hilbert_space_nstates(qs);

    psi0 = qdpack_state_alloc(qs);
    qdpack_operator_set_zero(U);

    for (i = 0; i < qsn; i++)
    {
        //
        // setup (non-physical) initial density matrix
        //
        qdpack_state_set_zero(psi0);
        qdpack_state_set(psi0, i, QDPACK_COMPLEX_ONE);

        //
        // evolve
        //

        if (qdpack_evolve_wf(sim, qs, psi0, psi_cb_func) != 0)
        {
            fprintf(stderr, "%s: Evolution of the quantum system failed.\n", __PRETTY_FUNCTION__);
            return -1;
        }

        //
        // copy the resulting density matrix into U
        //        
        for (j = 0; j < qsn; j++)
        {
            qdpack_operator_set(U, j, i, qdpack_operator_get(sim->rho_f, j, 0));    
        }
    }    
    
    qdpack_state_free(psi0);

    return 0;
}

/** ----------------------------------------------------------------------------
 * Calculates the evolution of the wave function using the monte carlo method.
 *
 */
int
qdpack_evolve_wfmc(qdpack_simulation_t *sim,
                   qdpack_hilbert_space_t *qs, 
                   qdpack_state_t *psi0,
                   state_cb_func_t psi_cb_func)
{
    int i, n, ntraj = 500;
    qdpack_state_t *psi, *psi1;
    qdpack_complex alpha;

    double r, t, dp = 0.0, dp_list[10], dp_cumsum[10]; // XXX: make dp_list size dynamic

    srand(time(NULL));

    /** -- setup hamiltonians --- */
    if (sim->H0_func)
    {
        if (sim->H0 == NULL)
        {
            sim->H0 = qdpack_operator_alloc(qs);
        }
        sim->H0_func(sim, qs, sim->H0);
    }
    else
    {   
        fprintf(stderr, "%s: No hamiltonian.", __PRETTY_FUNCTION__);
        return -1;    
    }

    if (sim->H1_func)
    {
        if (sim->H1 == NULL)
        {
            sim->H1 = qdpack_operator_alloc(qs);
        }
        sim->H1_func(sim, qs, sim->H1);
    }

    if (sim->Ht_func)
    {
        if (sim->H_t == NULL)
        {
            sim->H_t = qdpack_operator_alloc(qs);
        }
    }

    psi  = qdpack_state_alloc(qs);
    psi1 = qdpack_state_alloc(qs);

    /* -- calculate the effective non-hermitian hamiltonian --- */
    for (i = 0; i < sim->do_n; i++)
    {    
        if (sim->do_g1[i] != 0.0)
        {
            sim->do_ada[i] = qdpack_operator_alloc(qs);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->do_a[i], sim->do_a[i], QDPACK_COMPLEX_ZERO, sim->do_ada[i]);

            QDPACK_SET_COMPLEX(&alpha, 0.0, 0.5 * sim->do_g1[i]);
            qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, sim->do_a[i], sim->do_a[i], QDPACK_COMPLEX_ONE, sim->H0);
        }

        if (sim->do_g2[i] != 0.0)
        {
            QDPACK_SET_COMPLEX(&alpha, 0.0, 0.5 * sim->do_g2[i]);
            qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, sim->do_a[i], sim->do_a[i], QDPACK_COMPLEX_ONE, sim->H0);
        }
    }

    /* --- do monte carlo simulation --- */
    qdpack_state_memcpy(psi, psi0);
        
    for (t = 0; t <= sim->Tf; t += sim->dT)
    {
        psi_cb_func(sim, qs, psi, t);

        for (n = 0; n < ntraj; n++)
        {
            for (i = 0; i < sim->do_n; i++)
            {
                qdpack_complex z = qdpack_state_expectation_value(psi, sim->do_ada[i]);
                dp_list[i] = sim->do_g1[i] * QDPACK_REAL(z);
                dp = dp_list[i];
                dp_cumsum[i] = dp;
            }

            r = ((double)rand()) / RAND_MAX;
            
            if (dp < r)
            {   
                // normal evolution
                QDPACK_SET_COMPLEX(&alpha, 0.0, -0.5 * sim->dT);
                qdpack_operator_state_multiply(sim->H0, psi, psi1);
                qdpack_state_scale(psi1, alpha);
                QDPACK_SET_COMPLEX(&alpha, 1.0 / sqrt(1 - dp), 0.0);
                qdpack_state_scale(psi1, alpha);
                qdpack_state_add(psi, psi1);
            }
            else
            {
                // quantum jump

                r = ((double)rand()) / RAND_MAX;
                
                for (i = 0; dp_cumsum[i] < r; i++);

                printf("%s: performing quantum jump %d\n", __PRETTY_FUNCTION__, i);

                QDPACK_SET_COMPLEX(&alpha, sqrt(sim->do_g1[i]), 0.0);
                qdpack_operator_state_multiply(sim->do_a[i], psi, psi1);
                qdpack_state_scale(psi1, alpha);
                QDPACK_SET_COMPLEX(&alpha, 1.0 / sqrt(dp_list[i] / dp), 0.0);
                qdpack_state_scale(psi1, alpha);
                qdpack_state_add(psi, psi1);
            }
        }
    }

    return 0;
}

