//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file integrate.c
 *
 * @brief Functions for integrating systems of differential equations.  
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */


#include <math.h>
#include <stdio.h>

#include <qdpack/qdpack.h>

//
// operator integration
//

/**
 * \brief    Evolve a density matrix EQM using a standard Runge-Kutta.
 *
 * Standard Runge-Kutta method.
 *
 */
int 
odeint_operator_rk(qdpack_simulation_t *sim, 
                   qdpack_hilbert_space_t *qs, 
                   qdpack_operator_t *rho0, 
                   ode_operator_func_d_dt f_d_dt, 
                   operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *rho_t, *ws, *k1, *k2, *k3, *k4;
    qdpack_complex z;
    double t;

    rho_t  = qdpack_operator_alloc(qs);
    ws     = qdpack_operator_alloc(qs);
    k1     = qdpack_operator_alloc(qs);
    k2     = qdpack_operator_alloc(qs);
    k3     = qdpack_operator_alloc(qs);
    k4     = qdpack_operator_alloc(qs);

    qdpack_operator_memcpy(rho_t, rho0);
    
    if (rho_cb_func)
        rho_cb_func(sim, qs, rho_t, sim->Ti);
    
    for (t = sim->Ti; t <= sim->Tf; t += sim->dT)
    {
        qdpack_operator_memcpy(ws, rho_t);
        f_d_dt(sim, qs, ws, k1, t);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_operator_scale(k1, z);
        
        qdpack_operator_memcpy(ws, k1);
        QDPACK_SET_COMPLEX(&z, 0.5, 0);
        qdpack_operator_scale(ws, z);
        qdpack_operator_add(ws, rho_t);
        f_d_dt(sim, qs, ws, k2, t + sim->dT / 2.0);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_operator_scale(k2, z);
        
        qdpack_operator_memcpy(ws, k2);
        QDPACK_SET_COMPLEX(&z, 0.5, 0);
        qdpack_operator_scale(ws, z);
        qdpack_operator_add(ws, rho_t);
        f_d_dt(sim, qs, ws, k3, t + sim->dT / 2.0);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_operator_scale(k3, z);
        
        qdpack_operator_memcpy(ws, k3);
        qdpack_operator_add(ws, rho_t);
        f_d_dt(sim, qs, ws, k4, t + sim->dT);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_operator_scale(k4, z);
        
        QDPACK_SET_COMPLEX(&z, 2.0, 0);
        qdpack_operator_memcpy(ws, k1);
        qdpack_operator_scale(k2, z);
        qdpack_operator_add(ws, k2);    
        qdpack_operator_scale(k3, z);
        qdpack_operator_add(ws, k3);    
        qdpack_operator_add(ws, k4);    
        QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
        qdpack_operator_scale(ws, z);
        
        // rho_t(t+dt) = rho(t) + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
        qdpack_operator_add(rho_t, ws);    
        if (rho_cb_func)
            rho_cb_func(sim, qs, rho_t, t);
    }

    qdpack_operator_memcpy(sim->rho_f, rho_t);

    qdpack_operator_free(rho_t);
    qdpack_operator_free(ws);
    qdpack_operator_free(k1);
    qdpack_operator_free(k2);
    qdpack_operator_free(k3);
    qdpack_operator_free(k4);

    return 0;
}


/**
 * Two-step adaptive Runge-Kutta method with local error estimate
 * and step-length control.
 */
int 
odeint_operator_rk_adaptive(qdpack_simulation_t *sim, 
                            qdpack_hilbert_space_t *qs, 
                            qdpack_operator_t *rho0, 
                            ode_operator_func_d_dt f_d_dt, 
                            operator_cb_func_t rho_cb_func)
{
    qdpack_operator_t *rho1, *rho2;
    double t, e, emax = 0.0, etol = 0.001;

    rho1  = qdpack_operator_alloc(qs);
    rho2  = qdpack_operator_alloc(qs);
    qdpack_operator_memcpy(rho1, rho0);
    
    if (rho_cb_func)
        rho_cb_func(sim, qs, rho1, sim->Ti);
    
    for (t = sim->Ti; t <= sim->Tf - sim->dT; t += sim->dT)
    {
        //int rd;
        
        odeint_operator_rk_adaptive_rec(sim, qs, f_d_dt, rho1, rho2, 0, t, t+sim->dT, etol, &e);
        if (e > emax)
            emax = e;
        
        //printf("odeint_complex_rk_adaptive: Local error estimate %f (rec depth %d) at time t = %f\n", e, rd, t);
        
        qdpack_operator_memcpy(rho1, rho2);
        if (rho_cb_func)
            rho_cb_func(sim, qs, rho1, t + sim->dT);
    }
    
    qdpack_operator_memcpy(sim->rho_f, rho1);

    //printf("odeint_complex_rk_adaptive: Maximum (estimated) local error throughout simulation: %f\n", emax);
    
    qdpack_operator_free(rho1);
    qdpack_operator_free(rho2);

    return 0;
}

/**
 * Recursive method for calculating a two-step Runge-Kutta on shrinking interval.
 *
 */
int 
odeint_operator_rk_adaptive_rec(qdpack_simulation_t *sim, 
                                qdpack_hilbert_space_t *qs, 
                                ode_operator_func_d_dt f_d_dt,
                                qdpack_operator_t *rho0, 
                                qdpack_operator_t *rho_final,
                                int level, double ti, double tf, double etol, double *e)
{
    double e1, e2, le, h;
    qdpack_complex hz, z;
    qdpack_operator_t *rhot, *rho1, *rho2, *ws, *k1, *k2, *k3, *k4;
    unsigned int i, j; 
    int rd = level;

    rhot = qdpack_operator_alloc(qs);
    rho1 = qdpack_operator_alloc(qs);
    rho2 = qdpack_operator_alloc(qs);
    ws   = qdpack_operator_alloc(qs);
    k1   = qdpack_operator_alloc(qs); // can be optimized by using k and kacc
    k2   = qdpack_operator_alloc(qs);
    k3   = qdpack_operator_alloc(qs);
    k4   = qdpack_operator_alloc(qs);

    *e = 0.0;
    
    /*
     * calculate rho1 using a one-step RK 
     * 
     * rho1 = rho0 + 1/6 (k1 + 2 k2 + 2 k3 + k4)
     * alt: copy rho2 to rho1, if level != 0
     */
    h = tf - ti;
    QDPACK_SET_COMPLEX(&hz, h, 0.0);
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    
    // K1
    qdpack_operator_memcpy(ws, rho0);
    f_d_dt(sim, qs, ws, k1, ti);
    qdpack_operator_scale(k1, hz);
    
    // K2
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k2, ti + h / 2.0);
    qdpack_operator_scale(k2, hz);
    
    // K3
    qdpack_operator_memcpy(ws, k2);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k3, ti + h / 2.0);
    qdpack_operator_scale(k3, hz);
    
    // K4
    qdpack_operator_memcpy(ws, k3);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k4, ti + h);
    qdpack_operator_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(k2, z);
    qdpack_operator_add(ws, k2);    
    qdpack_operator_scale(k3, z);
    qdpack_operator_add(ws, k3);    
    qdpack_operator_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_operator_scale(ws, z);
        
    // rho_t(t+dt) = rho(t) + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_operator_memcpy(rho1, rho0);

    qdpack_operator_add(rho1, ws);    

    /* 
     * calculate rho2 using a two-step RK 
     *
     * 
     */
    h = (tf - ti) / 2.0;
    QDPACK_SET_COMPLEX(&hz, h, 0.0);
    
    /* --- STEP 1 --- */
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    // K1
    qdpack_operator_memcpy(ws, rho0);
    f_d_dt(sim, qs, ws, k1, ti);
    qdpack_operator_scale(k1, hz);
    
    // K2
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k2, ti + h / 2.0);
    qdpack_operator_scale(k2, hz);
    
    // K3
    qdpack_operator_memcpy(ws, k2);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k3, ti + h / 2.0);
    qdpack_operator_scale(k3, hz);
    
    // K4
    qdpack_operator_memcpy(ws, k3);
    qdpack_operator_add(ws, rho0);
    f_d_dt(sim, qs, ws, k4, ti + h);
    qdpack_operator_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(k2, z);
    qdpack_operator_add(ws, k2);    
    qdpack_operator_scale(k3, z);
    qdpack_operator_add(ws, k3);    
    qdpack_operator_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_operator_scale(ws, z);
        
    // rhot = rho0 + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_operator_memcpy(rhot, rho0);
    qdpack_operator_add(rhot, ws);    
    
    /* --- STEP 2 --- */
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    // K1
    qdpack_operator_memcpy(ws, rhot);
    f_d_dt(sim, qs, ws, k1, ti + h);
    qdpack_operator_scale(k1, hz);
    
    // K2
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rhot);
    f_d_dt(sim, qs, ws, k2, ti + 3 * h / 2.0);
    qdpack_operator_scale(k2, hz);
    
    // K3
    qdpack_operator_memcpy(ws, k2);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(ws, rhot);
    f_d_dt(sim, qs, ws, k3, ti + 3 * h / 2.0);
    qdpack_operator_scale(k3, hz);
    
    // K4
    qdpack_operator_memcpy(ws, k3);
    qdpack_operator_add(ws, rhot);
    f_d_dt(sim, qs, ws, k4, ti + 2 * h);
    qdpack_operator_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_operator_memcpy(ws, k1);
    qdpack_operator_scale(k2, z);
    qdpack_operator_add(ws, k2);    
    qdpack_operator_scale(k3, z);
    qdpack_operator_add(ws, k3);    
    qdpack_operator_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_operator_scale(ws, z);
        
    // rho2 = rhot + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_operator_memcpy(rho2, rhot);
    qdpack_operator_add(rho2, ws);    

    /* 
     * Combine rho1 and rho2 into a new estimate rho 
     *
     */
    qdpack_operator_memcpy(rho_final, rho2);
    QDPACK_SET_COMPLEX(&z, 16.0, 0.0);
    qdpack_operator_scale(rho_final, z);
    qdpack_operator_sub(rho_final, rho1);
    QDPACK_SET_COMPLEX(&z, 1/15.0, 0.0);
    qdpack_operator_scale(rho_final, z);

    /* 
     * Find maximum deviation 
     *
     */
    for (i = 0; i < rho_final->m; i++)    
    {
        for (j = 0; j < rho_final->n; j++)
        {
            qdpack_complex r_ij, r1_ij, r2_ij;
            
            r_ij  = qdpack_operator_get(rho_final,i,j);
        
            if (qdpack_complex_abs(r_ij) != 0.0)
            {
                r1_ij = qdpack_operator_get(rho1,i,j);
                r2_ij = qdpack_operator_get(rho2,i,j);
                
                z = qdpack_complex_sub(r1_ij, r2_ij);
                
                le = (qdpack_complex_abs(z)) / qdpack_complex_abs(r_ij);
                
                if (le > *e && qdpack_complex_abs(r_ij) > etol/10.0)
                {
                    *e = le;
                }
            }
        }
    }
    
//    printf("DEBUG: odeint_rk_adaptive_rec: local error e = %f,\t[ti, tf] = [%f, %f], level = %d\n", *e, ti, tf, level);

    /*
     * Recursively refine intervals if errors are too large
     *
     */
    if (*e > etol && level < RK_REC_MAX_LEVEL)
    {
        int rd1, rd2;
        rd1 = odeint_operator_rk_adaptive_rec(sim, qs, f_d_dt, rho0, rhot,      level + 1, ti, ti+(tf-ti)/2.0, etol, &e1);
        rd2 = odeint_operator_rk_adaptive_rec(sim, qs, f_d_dt, rhot, rho_final, level + 1, ti+(tf-ti)/2.0, tf, etol, &e2);
        *e = (e1 > e2) ? e1 : e2;
        rd = (rd1 > rd2) ? rd1 : rd2;
    }
    else
    {
        rd = level;
    }

    qdpack_operator_free(rhot);
    qdpack_operator_free(rho1);
    qdpack_operator_free(rho2);
    qdpack_operator_free(ws);
    qdpack_operator_free(k1);
    qdpack_operator_free(k2);
    qdpack_operator_free(k3);
    qdpack_operator_free(k4);
    
    return rd;
}


/**
 * Find maximum deviation between the two matrices m1 and m2.
 *
 */
double
operator_max_deviation(qdpack_operator_t *m1, qdpack_operator_t *m2, double etol)
{
    unsigned int i, j;
    double emax = 0.0, le;
    
    for (i = 0; i < m1->m; i++)    
    {
        for (j = 0; j < m2->n; j++)
        {
            qdpack_complex m1_ij, m2_ij;
            
            m2_ij  = qdpack_operator_get(m2,i,j);
        
            if (qdpack_complex_abs(m2_ij) != 0.0)
            {
                m1_ij = qdpack_operator_get(m1,i,j);
            
                le = qdpack_complex_abs(m1_ij) / qdpack_complex_abs(m2_ij);
            
                if (le > emax && qdpack_complex_abs(m2_ij) > etol/10.0)
                {
                    emax = le;
                }
            }
        }
    }

    return emax;
}


//
// state vector integration
//
/**
 * \brief    Evolve a density matrix EQM using a standard Runge-Kutta.
 *
 * Standard Runge-Kutta method.
 *
 */
int 
odeint_state_rk(qdpack_simulation_t *sim, 
                qdpack_hilbert_space_t *qs, 
                qdpack_state_t *psi0, 
                ode_state_func_d_dt f_d_dt, 
                state_cb_func_t psi_cb_func)
{
    qdpack_state_t *psi_t, *ws, *k1, *k2, *k3, *k4;
    qdpack_complex z;
    double t;

    psi_t  = qdpack_state_alloc(qs);
    ws     = qdpack_state_alloc(qs);
    k1     = qdpack_state_alloc(qs);
    k2     = qdpack_state_alloc(qs);
    k3     = qdpack_state_alloc(qs);
    k4     = qdpack_state_alloc(qs);

    qdpack_state_memcpy(psi_t, psi0);
    
    if (psi_cb_func)
        psi_cb_func(sim, qs, psi_t, sim->Ti);
    
    for (t = sim->Ti; t <= sim->Tf; t += sim->dT)
    {
        qdpack_state_memcpy(ws, psi_t);
        f_d_dt(sim, qs, ws, k1, t);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_state_scale(k1, z);
        
        qdpack_state_memcpy(ws, k1);
        QDPACK_SET_COMPLEX(&z, 0.5, 0);
        qdpack_state_scale(ws, z);
        qdpack_state_add(ws, psi_t);
        f_d_dt(sim, qs, ws, k2, t + sim->dT / 2.0);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_state_scale(k2, z);
        
        qdpack_state_memcpy(ws, k2);
        QDPACK_SET_COMPLEX(&z, 0.5, 0);
        qdpack_state_scale(ws, z);
        qdpack_state_add(ws, psi_t);
        f_d_dt(sim, qs, ws, k3, t + sim->dT / 2.0);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_state_scale(k3, z);
        
        qdpack_state_memcpy(ws, k3);
        qdpack_state_add(ws, psi_t);
        f_d_dt(sim, qs, ws, k4, t + sim->dT);
        QDPACK_SET_COMPLEX(&z, sim->dT, 0);
        qdpack_state_scale(k4, z);
        
        QDPACK_SET_COMPLEX(&z, 2.0, 0);
        qdpack_state_memcpy(ws, k1);
        qdpack_state_scale(k2, z);
        qdpack_state_add(ws, k2);    
        qdpack_state_scale(k3, z);
        qdpack_state_add(ws, k3);    
        qdpack_state_add(ws, k4);    
        QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
        qdpack_state_scale(ws, z);
        
        // psi_t(t+dt) = psi(t) + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
        qdpack_state_add(psi_t, ws);    
        if (psi_cb_func)
            psi_cb_func(sim, qs, psi_t, t);
    }

    qdpack_state_memcpy(sim->psi_f, psi_t);

    qdpack_state_free(psi_t);
    qdpack_state_free(ws);
    qdpack_state_free(k1);
    qdpack_state_free(k2);
    qdpack_state_free(k3);
    qdpack_state_free(k4);

    return 0;
}


/**
 * Two-step adaptive Runge-Kutta method with local error estimate
 * and step-length control.
 */
int 
odeint_state_rk_adaptive(qdpack_simulation_t *sim, 
                         qdpack_hilbert_space_t *qs, 
                         qdpack_state_t *psi0, 
                         ode_state_func_d_dt f_d_dt, 
                         state_cb_func_t psi_cb_func)
{
    qdpack_state_t *psi1, *psi2;
    double t, e, emax = 0.0, etol = 0.001;

    psi1  = qdpack_state_alloc(qs);
    psi2  = qdpack_state_alloc(qs);
    qdpack_state_memcpy(psi1, psi0);
    
    if (psi_cb_func)
        psi_cb_func(sim, qs, psi1, sim->Ti);
    
    for (t = sim->Ti; t <= sim->Tf - sim->dT; t += sim->dT)
    {
        //int rd;
        
        odeint_state_rk_adaptive_rec(sim, qs, f_d_dt, psi1, psi2, 0, t, t+sim->dT, etol, &e);
        if (e > emax)
            emax = e;
        
        //printf("odeint_complex_rk_adaptive: Local error estimate %f (rec depth %d) at time t = %f\n", e, rd, t);
        
        qdpack_state_memcpy(psi1, psi2);
        if (psi_cb_func)
            psi_cb_func(sim, qs, psi1, t + sim->dT);
    }
    
    qdpack_state_memcpy(sim->psi_f, psi1);

    //printf("odeint_complex_rk_adaptive: Maximum (estimated) local error throughout simulation: %f\n", emax);
    
    qdpack_state_free(psi1);
    qdpack_state_free(psi2);

    return 0;
}

/**
 * Recursive method for calculating a two-step Runge-Kutta on shrinking interval.
 *
 */
int 
odeint_state_rk_adaptive_rec(qdpack_simulation_t *sim, 
                             qdpack_hilbert_space_t *qs, 
                             ode_state_func_d_dt f_d_dt,
                             qdpack_state_t *psi0, 
                             qdpack_state_t *psi_final,
                             int level, double ti, double tf, double etol, double *e)
{
    double e1, e2, le, h;
    qdpack_complex hz, z;
    qdpack_state_t *psit, *psi1, *psi2, *ws, *k1, *k2, *k3, *k4;
    unsigned int i; 
    int rd = level;

    psit = qdpack_state_alloc(qs);
    psi1 = qdpack_state_alloc(qs);
    psi2 = qdpack_state_alloc(qs);
    ws   = qdpack_state_alloc(qs);
    k1   = qdpack_state_alloc(qs); // can be optimized by using k and kacc
    k2   = qdpack_state_alloc(qs);
    k3   = qdpack_state_alloc(qs);
    k4   = qdpack_state_alloc(qs);

    *e = 0.0;
    
    /*
     * calculate psi1 using a one-step RK 
     * 
     * psi1 = psi0 + 1/6 (k1 + 2 k2 + 2 k3 + k4)
     * alt: copy psi2 to psi1, if level != 0
     */
    h = tf - ti;
    QDPACK_SET_COMPLEX(&hz, h, 0.0);
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    
    // K1
    qdpack_state_memcpy(ws, psi0);
    f_d_dt(sim, qs, ws, k1, ti);
    qdpack_state_scale(k1, hz);
    
    // K2
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k2, ti + h / 2.0);
    qdpack_state_scale(k2, hz);
    
    // K3
    qdpack_state_memcpy(ws, k2);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k3, ti + h / 2.0);
    qdpack_state_scale(k3, hz);
    
    // K4
    qdpack_state_memcpy(ws, k3);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k4, ti + h);
    qdpack_state_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(k2, z);
    qdpack_state_add(ws, k2);    
    qdpack_state_scale(k3, z);
    qdpack_state_add(ws, k3);    
    qdpack_state_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_state_scale(ws, z);
        
    // psi_t(t+dt) = psi(t) + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_state_memcpy(psi1, psi0);

    qdpack_state_add(psi1, ws);    

    /* 
     * calculate psi2 using a two-step RK 
     *
     * 
     */
    h = (tf - ti) / 2.0;
    QDPACK_SET_COMPLEX(&hz, h, 0.0);
    
    /* --- STEP 1 --- */
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    // K1
    qdpack_state_memcpy(ws, psi0);
    f_d_dt(sim, qs, ws, k1, ti);
    qdpack_state_scale(k1, hz);
    
    // K2
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k2, ti + h / 2.0);
    qdpack_state_scale(k2, hz);
    
    // K3
    qdpack_state_memcpy(ws, k2);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k3, ti + h / 2.0);
    qdpack_state_scale(k3, hz);
    
    // K4
    qdpack_state_memcpy(ws, k3);
    qdpack_state_add(ws, psi0);
    f_d_dt(sim, qs, ws, k4, ti + h);
    qdpack_state_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(k2, z);
    qdpack_state_add(ws, k2);    
    qdpack_state_scale(k3, z);
    qdpack_state_add(ws, k3);    
    qdpack_state_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_state_scale(ws, z);
        
    // psit = psi0 + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_state_memcpy(psit, psi0);
    qdpack_state_add(psit, ws);    
    
    /* --- STEP 2 --- */
    QDPACK_SET_COMPLEX(&z, 0.5, 0);
    // K1
    qdpack_state_memcpy(ws, psit);
    f_d_dt(sim, qs, ws, k1, ti + h);
    qdpack_state_scale(k1, hz);
    
    // K2
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psit);
    f_d_dt(sim, qs, ws, k2, ti + 3 * h / 2.0);
    qdpack_state_scale(k2, hz);
    
    // K3
    qdpack_state_memcpy(ws, k2);
    qdpack_state_scale(ws, z);
    qdpack_state_add(ws, psit);
    f_d_dt(sim, qs, ws, k3, ti + 3 * h / 2.0);
    qdpack_state_scale(k3, hz);
    
    // K4
    qdpack_state_memcpy(ws, k3);
    qdpack_state_add(ws, psit);
    f_d_dt(sim, qs, ws, k4, ti + 2 * h);
    qdpack_state_scale(k4, hz);
        
    QDPACK_SET_COMPLEX(&z, 2.0, 0);
    qdpack_state_memcpy(ws, k1);
    qdpack_state_scale(k2, z);
    qdpack_state_add(ws, k2);    
    qdpack_state_scale(k3, z);
    qdpack_state_add(ws, k3);    
    qdpack_state_add(ws, k4);    
    QDPACK_SET_COMPLEX(&z, 1/6.0, 0);
    qdpack_state_scale(ws, z);
        
    // psi2 = psit + 1/6 * (k1 + 2 k2 + 2 k3 + k4)
    qdpack_state_memcpy(psi2, psit);
    qdpack_state_add(psi2, ws);    

    /* 
     * Combine psi1 and psi2 into a new estimate psi 
     *
     */
    qdpack_state_memcpy(psi_final, psi2);
    QDPACK_SET_COMPLEX(&z, 16.0, 0.0);
    qdpack_state_scale(psi_final, z);
    qdpack_state_sub(psi_final, psi1);
    QDPACK_SET_COMPLEX(&z, 1/15.0, 0.0);
    qdpack_state_scale(psi_final, z);

    /* 
     * Find maximum deviation 
     *
     */
    for (i = 0; i < psi_final->n; i++)    
    {
        qdpack_complex r_i, r1_i, r2_i;
        
        r_i  = qdpack_state_get(psi_final,i);
    
        if (qdpack_complex_abs(r_i) != 0.0)
        {
            r1_i = qdpack_state_get(psi1,i);
            r2_i = qdpack_state_get(psi2,i);
            
            z = qdpack_complex_sub(r1_i, r2_i);
            
            le = (qdpack_complex_abs(z)) / qdpack_complex_abs(r_i);
            
            if (le > *e && qdpack_complex_abs(r_i) > etol/10.0)
            {
                *e = le;
            }
        }
    }

    /*
     * Recursively refine intervals if errors are too large
     *
     */
    if (*e > etol && level < RK_REC_MAX_LEVEL)
    {
        int rd1, rd2;
        rd1 = odeint_state_rk_adaptive_rec(sim, qs, f_d_dt, psi0, psit,      level + 1, ti, ti+(tf-ti)/2.0, etol, &e1);
        rd2 = odeint_state_rk_adaptive_rec(sim, qs, f_d_dt, psit, psi_final, level + 1, ti+(tf-ti)/2.0, tf, etol, &e2);
        *e = (e1 > e2) ? e1 : e2;
        rd = (rd1 > rd2) ? rd1 : rd2;
    }
    else
    {
        rd = level;
    }

    qdpack_state_free(psit);
    qdpack_state_free(psi1);
    qdpack_state_free(psi2);
    qdpack_state_free(ws);
    qdpack_state_free(k1);
    qdpack_state_free(k2);
    qdpack_state_free(k3);
    qdpack_state_free(k4);
    
    return rd;
}


/**
 * Find maximum deviation between the two matrices m1 and m2.
 *
 */
double
state_max_deviation(qdpack_state_t *m1, qdpack_state_t *m2, double etol)
{
    unsigned int i;
    double emax = 0.0, le;
    
    for (i = 0; i < m1->n; i++)    
    {
        qdpack_complex m1_i, m2_i;
            
        m2_i  = qdpack_state_get(m2,i);
    
        if (qdpack_complex_abs(m2_i) != 0.0)
        {
            m1_i = qdpack_state_get(m1,i);
        
            le = qdpack_complex_abs(m1_i) / qdpack_complex_abs(m2_i);
        
            if (le > emax && qdpack_complex_abs(m2_i) > etol/10.0)
            {
                emax = le;
            }
        }
    }

    return emax;
}
