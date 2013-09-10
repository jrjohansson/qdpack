//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file hamiltonian.c
 *
 * @brief Functions for generating hamiltonians in matrix form. 
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */

#include <stdio.h>
#include <math.h>

#include <qdpack/qdpack.h>

//==============================================================================
// Generic time-dependence
//==============================================================================

/**
 *
 * For cases where the time-dependence of the Hamiltonian takes the form
 * H = H0 + H1 * f(t), this function is used to implement the specific
 * time-dependence in f(t).
 *
 * @param    H_t  Matrix where the resulting time-dependent Hamiltonian
 *                at time t should be stored (output from this function).
 * @param    H0   The constant part of the Hamiltonian.
 * @param    H1   The part of the Hamiltonian which is modulated by f(t).
 * @param    t    The instant of time for which we should calculate H(t).
 * @param    sim  Pointer to the data structure that contains all the 
 *                simulation parameters.
 *
 * @return On success 0, on failure 1.
 */
int
hamiltonian_t(qdpack_simulation_t *sim, qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t)
{
    if (sim->h_td_A == 0.0)
    {
        if (t == sim->Ti)
            qdpack_operator_memcpy(H_t, H0);
    }
    else 
    {
        // f(t) = A * sin(w * t)

        qdpack_complex z;

        qdpack_operator_memcpy(H_t, H1);
        QDPACK_SET_COMPLEX(&z, sim->h_td_A * sin(sim->h_td_w * t), 0.0);
        qdpack_operator_scale(H_t, z);
        qdpack_operator_add(H_t, H0);
    }

    return 0;
}

//==============================================================================
// coupled qubits
//==============================================================================

int
hamiltonian_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{    
    return operator_sigma_x(H, qs, 0);
}

int 
hamiltonian_z_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    return operator_sigma_z(H, qs, 0);
}

int
hamiltonian_x_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    return operator_sigma_x(H, qs, 0);
}

/** 
 * Generate the hamiltonian for a N-TLS system
 *
 * Two-level system (qubit) parameters: 
 *
 * H = sum_i(-0.5 espilon[i] sigma_z(i) - 0.5 delta[i] sigma_x(i)) + sum_ij(- 0.5 lambda[i][j] sigma_x[i] sigma_x[j]) 
 * 
 * The paramters espilon, delta and lambda are defined in h_qubits_parameters_t.
 */

int
hamiltonian_qubits(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    int i, j;
    qdpack_operator_t *s1, *s2, *s3;

    h_qubits_parameters_t *param = (h_qubits_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);

    s1 = qdpack_operator_alloc(qs);
    s2 = qdpack_operator_alloc(qs);
    s3 = qdpack_operator_alloc(qs);

    /*
     * the TLS energies
     */
    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_complex z;

        operator_sigma_z(s1, qs, i);
        QDPACK_SET_COMPLEX(&z, -0.5 * param->epsilon[i], 0.0);
        qdpack_operator_scale(s1, z);
        qdpack_operator_add(H, s1);

        operator_sigma_x(s1, qs, i);
        QDPACK_SET_COMPLEX(&z, -0.5 * param->delta[i], 0.0);
        qdpack_operator_scale(s1, z);
        qdpack_operator_add(H, s1);
    }

    /*
     * TLS interactions
     */
    for (i = 0; i < qs->nsubsys; i++)
    {
        for (j = 0; j < qs->nsubsys; j++)
        {
            if (param->lambda[i][j] != 0.0)
            {
                qdpack_complex z;

                operator_sigma_x(s1, qs, i);
                operator_sigma_x(s2, qs, j);

                QDPACK_SET_COMPLEX(&z, -0.5 * param->lambda[i][j], 0.0);
                qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, s1, s2, QDPACK_COMPLEX_ZERO, s3);
                qdpack_operator_add(H, s3);
            }
        }
    }

    qdpack_operator_free(s1);
    qdpack_operator_free(s2);
    qdpack_operator_free(s3);

    return 0;
}

//==============================================================================
// spin chain
//==============================================================================

/*
 * Generate the hamiltonian for a N-TLS system 
 *
 */
int 
hamiltonian_spinchain(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    int i, j;
    qdpack_operator_t *s1, *s2, *s3, *s1b, *s2b;
    
    h_qubits_parameters_t *param = (h_qubits_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);
    
    s1  = qdpack_operator_alloc(qs);
    s1b = qdpack_operator_alloc(qs);
    s2  = qdpack_operator_alloc(qs);
    s2b = qdpack_operator_alloc(qs);
    s3  = qdpack_operator_alloc(qs);

    /*
     * the TLS energies
     */
    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_complex z;
        operator_sigma_z(s1, qs, i);
        QDPACK_SET_COMPLEX(&z, -0.5 * param->epsilon[i], 0.0);
        qdpack_operator_scale(s1, z);
        qdpack_operator_add(H, s1);
                
        operator_sigma_x(s1, qs, i);
        QDPACK_SET_COMPLEX(&z, -0.5 * param->delta[i], 0.0);
        qdpack_operator_scale(s1, z);
        qdpack_operator_add(H, s1);
    }

    /*
     * TLS interactions
     */
    for (i = 0; i < qs->nsubsys; i++)
    {
        for (j = 0; j < qs->nsubsys; j++)
        {
            if (param->lambda[i][j] != 0.0)
            {
                qdpack_complex z;
                double theta;
                
                operator_sigma_x(s1,  qs, i);
                operator_sigma_z(s1b, qs, i);
                operator_sigma_x(s2,  qs, j);
                operator_sigma_z(s2b, qs, j);
                
                theta = tan(param->delta[i] / param->epsilon[i]);
                QDPACK_SET_COMPLEX(&z, sin(theta), 0.0);
                qdpack_operator_scale(s1, z);
                QDPACK_SET_COMPLEX(&z, cos(theta), 0.0);
                qdpack_operator_scale(s1b, z);
                qdpack_operator_add(s1, s1b);
                theta = tan(param->delta[j] / param->epsilon[j]);
                QDPACK_SET_COMPLEX(&z, sin(theta), 0.0);
                qdpack_operator_scale(s2, z);
                QDPACK_SET_COMPLEX(&z, cos(theta), 0.0);
                qdpack_operator_scale(s2b, z);
                qdpack_operator_add(s2, s2b);
                
                QDPACK_SET_COMPLEX(&z, -0.5 * param->lambda[i][j], 0.0);
                qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, s1, s2, QDPACK_COMPLEX_ZERO, s3);
                
                qdpack_operator_add(H, s3);                
            }
        }
    }

    qdpack_operator_free(s3);
    qdpack_operator_free(s1);
    qdpack_operator_free(s1b);
    qdpack_operator_free(s2);
    qdpack_operator_free(s2b);
    
    return 0;

}

// -----------------------------------------------------------------------------
/**
 * \brief    Generate the hamiltonian for two coupled harmonic oscillators
 *
 * This function allocate and initiates a matrix representation of a hamiltonian
 * of two coupled harmonic oscillators.
 *
 * @param    qs    Data structure that defines the quantum system.
 * @param    param    Data structure that contains all the parameters.
 *
 */
int 
hamiltonian_2ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *b, *bd, *ws, *Hint;
    qdpack_complex alpha, beta, z;
    
    h_ho_parameters_t *param = (h_ho_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);

    Hint = qdpack_operator_alloc(qs);

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    
    /* h_ w1 (ad a + 1/2) */
    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs,0, 0);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
    QDPACK_SET_COMPLEX(&beta, 0.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, beta, H);
    operator_unit(ws, qs, -1);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[0], 0.0);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(H, ws);
    
    /* h_ w (bd b + 1/2) */
    operator_ho_lowering(b, qs, 1, 0);
    operator_ho_raising(bd, qs, 1, 0);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[1], 0.0);

    QDPACK_SET_COMPLEX(&beta, 1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, bd, b, beta, H);
    operator_unit(ws, qs, -1);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[1], 0.0);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(H, ws);
    
    // g (a + ad) (b + bd)
    qdpack_operator_add(a, ad);
    qdpack_operator_add(b, bd);
    QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    QDPACK_SET_COMPLEX(&beta, 1.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, b, beta, H);

    // g (a + ad) (b + bd)^2
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    // g a * b^2
    //QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, b, b, beta, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, Hint, beta, H);
    //qdpack_operator_scale(Hint, );
    //qdpack_operator_add(H, Hint);


    // g (a + ad)^2 (b + bd)
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, a, a, QDPACK_COMPLEX_ZERO, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, Hint, b, QDPACK_COMPLEX_ONE, H);
    
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    qdpack_operator_free(ws);
    qdpack_operator_free(Hint);

    return 0;
}


// -----------------------------------------------------------------------------
/**
 * Generate the hamiltonian for the 2-HO driving field
 *
 */
int
hamiltonian_2ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad; //, *b, *bd;
    
    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);

    //qdpack_operator_set_zero(H);
    
    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs, 0, 0);
    //b  = operator_ho_lowering(qs, 1);
    //bd = operator_ho_raising(qs, 1);

    // a bd + ad b
    //QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, bd, beta, H);
    //QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, b, beta, H);
    
    // a + ad
    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, a);
    qdpack_operator_add(H, ad);
    
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    //qdpack_operator_free(b);
    //qdpack_operator_free(bd);
    
    return 0;
}

// -----------------------------------------------------------------------------
// OPTO-MECHANIAL SYSTEM: 2HO + nonlinear interaction
//
// -----------------------------------------------------------------------------

/**
 * \brief    Generate the hamiltonian for two coupled harmonic oscillators
 *
 * This function allocate and initiates a matrix representation of a hamiltonian
 * of two coupled harmonic oscillators.
 *
 * @param    qs    Data structure that defines the quantum system.
 * @param    param    Data structure that contains all the parameters.
 *
 */
int 
hamiltonian_2ho_optmech(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *na, *b, *bd, *nb, *ws, *Hint;
    qdpack_complex alpha;
    
    h_ho_parameters_t *param = (h_ho_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);

    Hint = qdpack_operator_alloc(qs);

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    na = qdpack_operator_alloc(qs);
    nb = qdpack_operator_alloc(qs);
    
    /* --- h_ w1 (ad a + 1/2) --- */
    printf("%s: ho_w[0] = %f [offset = %d]\n", __PRETTY_FUNCTION__, param->ho_w[0], param->n_offset);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);

    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs, 0, 0);
    operator_ho_N(na, qs, 0, param->n_offset);
    //printf("Na =\n"); qdpack_operator_print_real(na, ns);

    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ZERO, H);
//    qdpack_operator_scale(na, alpha);
//    qdpack_operator_add(H, na);

    /* --- h_ w (bd b + 1/2) --- */
    printf("%s: ho_w[1] = %f\n", __PRETTY_FUNCTION__, param->ho_w[1]);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[1], 0.0);

    operator_ho_lowering(b, qs, 1, 0);
    operator_ho_raising(bd, qs, 1, 0);
    operator_ho_N(nb, qs, 1, param->n_offset);
    //printf("Nb =\n"); qdpack_operator_print_real(nb, ns);

    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, bd, b, QDPACK_COMPLEX_ONE, H);
//    qdpack_operator_scale(nb, alpha);
//    qdpack_operator_add(H, nb);
    

    //operator_unit(ws, qs, -1);
    //QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[1], 0.0);
    //qdpack_operator_scale(ws, z);
    //qdpack_operator_add(H, ws);
    
    // g (a + ad) (b + bd)
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);

    // g (a + ad) (b + bd)^2
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    // g a * b^2
    //QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, b, b, beta, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, Hint, beta, H);
    //qdpack_operator_scale(Hint, );
    //qdpack_operator_add(H, Hint);


    // g (a + ad)^2 (b + bd)
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, a, a, QDPACK_COMPLEX_ZERO, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, Hint, b, QDPACK_COMPLEX_ONE, H);
    //qdpack_operator_free(Hint);

    printf("%s: ho_g = %f\n", __PRETTY_FUNCTION__, param->ho_g);
    // g (ad a) (b + bd) = g na (b + bd)
    //qdpack_operator_add(a, ad);
    qdpack_operator_add(b, bd);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, a, a, QDPACK_COMPLEX_ZERO, Hint);
    QDPACK_SET_COMPLEX(&alpha, param->ho_g/2.0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, na, b, QDPACK_COMPLEX_ONE, H);
    

    /* --- driving term in the rotating frame --- */
    printf("%s: ho_wd = h_td_A = %f\n", __PRETTY_FUNCTION__, sim->h_td_A);
    QDPACK_SET_COMPLEX(&alpha, sim->h_td_A, 0.0);
    qdpack_operator_add(a, ad);
    qdpack_operator_scale(a, alpha);
    qdpack_operator_add(H, a);

    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(na);
    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    qdpack_operator_free(nb);
    qdpack_operator_free(ws);
    qdpack_operator_free(Hint);

    return 0;
}

/**
 * \brief    Generate the hamiltonian for two coupled harmonic oscillators,
 * with linearized interaction (shifted operators)
 *
 * This function allocate and initiates a matrix representation of a hamiltonian
 * of two coupled harmonic oscillators (optical cavity coupled to a mechanical
 * resonator), with linearized interaction and in a rotating frame (s.t. the
 * cavity is nearly resonant with the mechanical oscillator).
 *
 * @param    qs    Data structure that defines the quantum system.
 * @param    param    Data structure that contains all the parameters.
 *
 */
int 
hamiltonian_2ho_optmech_linearized(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *na, *b, *bd;
    qdpack_complex alpha;
    
    h_ho_parameters_t *param = (h_ho_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    na = qdpack_operator_alloc(qs);

    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    
    /* --- h_ w1 (ad a + 1/2) --- */
    printf("%s: ho_w[0] = %f [offset = %d]\n", __PRETTY_FUNCTION__, param->ho_w[0], param->n_offset);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);

    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs, 0, 0);

    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ZERO, H);

    /* --- h_ w (bd b + 1/2) --- */
    printf("%s: ho_w[1] = %f\n", __PRETTY_FUNCTION__, param->ho_w[1]);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[1], 0.0);

    operator_ho_lowering(b, qs, 1, 0);
    operator_ho_raising(bd, qs, 1, 0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, bd, b, QDPACK_COMPLEX_ONE, H);
    

    /* -- couplings --- */
    qdpack_operator_add(a, ad);
    qdpack_operator_add(b, bd);

    /* --- original coupling --- * /
    // g (ad a) (b + bd) = g na (b + bd)
 
    printf("%s: ho_g = %f\n", __PRETTY_FUNCigned and manufactured by IBM, the ThinkPad brand was purchased by Lenovo in 2005.
﻿﻿TION__, param->ho_g);
    QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, na, b, QDPACK_COMPLEX_ONE, H);

    */

    /* --- effective coupling --- */
    // g = neta * ho_w[1]
    // Omega / Delta g (a + ad) (b + bd) 

    printf("%s: Omega ho_g / Delta = %f\n", __PRETTY_FUNCTION__, sim->h_td_A * param->ho_g / param->ho_w[0]);
    QDPACK_SET_COMPLEX(&alpha, (sim->h_td_A * param->ho_g / param->ho_w[0]), 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);

    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(na);

    qdpack_operator_free(b);
    qdpack_operator_free(bd);


    return 0;
}


// -----------------------------------------------------------------------------
// OPTO-MECHANIAL SYSTEM: 2HO + nonlinear interaction? : 
// for sideband cooling calculations
//
// -----------------------------------------------------------------------------

/**
 * \brief    Generate the hamiltonian for two coupled harmonic oscillators
 *
 * This function allocate and initiates a matrix representation of a hamiltonian
 * of two coupled harmonic oscillators.
 *
 * @param    qs    Data structure that defines the quantum system.
 * @param    param    Data structure that contains all the parameters.
 *
 */
int 
hamiltonian_2ho_cooling(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *na, *b, *bd, *nb, *ws, *Hint;
    qdpack_complex alpha;
    
    h_ho_parameters_t *param = (h_ho_parameters_t *)sim->parameters;

    qdpack_operator_set_zero(H);

    Hint = qdpack_operator_alloc(qs);

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    na = qdpack_operator_alloc(qs);
    nb = qdpack_operator_alloc(qs);
    
    /* --- h_ w1 (ad a + 1/2) --- */
    printf("%s: ho_w[0] = %f [offset = %d]\n", __PRETTY_FUNCTION__, param->ho_w[0], param->n_offset);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);

    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs, 0, 0);
    operator_ho_N(na, qs, 0, param->n_offset);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ZERO, H);

    /* --- h_ w2 (bd b + 1/2) --- */
    printf("%s: ho_w[1] = %f\n", __PRETTY_FUNCTION__, param->ho_w[1]);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[1], 0.0);

    operator_ho_lowering(b, qs, 1, 0);
    operator_ho_raising(bd, qs, 1, 0);
    operator_ho_N(nb, qs, 1, param->n_offset);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, bd, b, QDPACK_COMPLEX_ONE, H);


    /* --- interaction --- */

    // g (a + ad) (b + bd)
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);

    // g (a + ad) (b + bd)^2
    //qdpack_operator_add(a, ad);
    //qdpack_operator_add(b, bd);
    // g a * b^2
    //QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  0.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, b, b, QDPACK_COMPLEX_ZERO, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    //QDPACK_SET_COMPLEX(&beta,  1.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, Hint, beta, H);
    //qdpack_operator_scale(Hint, );
    //qdpack_operator_add(H, Hint);


    // g (a + ad)^2 (b + bd)
    qdpack_operator_add(a, ad);
    qdpack_operator_add(b, bd);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, a, a, QDPACK_COMPLEX_ZERO, Hint);
    QDPACK_SET_COMPLEX(&alpha, param->ho_g, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, Hint, b, QDPACK_COMPLEX_ONE, H);

    // g (ad a) (b + bd) = g na (b + bd)
    //qdpack_operator_add(b, bd);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, ad, a, QDPACK_COMPLEX_ZERO, Hint);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_g/2.0, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, Hint, b, QDPACK_COMPLEX_ONE, H);
    

    /* --- driving term in the rotating frame --- */
    printf("%s: ho_wd = %f\n", __PRETTY_FUNCTION__, param->ho_wd);
    QDPACK_SET_COMPLEX(&alpha, sim->h_td_A, 0.0);
    qdpack_operator_add(a, ad);
    qdpack_operator_scale(a, alpha);
    qdpack_operator_add(H, a);

    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(na);
    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    qdpack_operator_free(nb);
    qdpack_operator_free(ws);
    qdpack_operator_free(Hint);

    return 0;
}


//!-----------------------------------------------------------------------------
// @brief   Generate squeezing hamiltonian for a cavity mode 
//
// This function allocate and initiates a matrix representation of a hamiltonian
// for a parametric amplification: H = i eps [ (ad^\dag)^2 - a^2 ]
//
// @param    qs    Data structure that defines the quantum system.
// @param    param    Data structure that contains all the parameters.
//------------------------------------------------------------------------------
int 
hamiltonian_squeezing(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *aa, *adad;
    qdpack_complex z;
    
    h_ho_parameters_t *param = (h_ho_parameters_t *)sim->parameters;

    a    = qdpack_operator_alloc(qs);
    ad   = qdpack_operator_alloc(qs);
    aa   = qdpack_operator_alloc(qs);
    adad = qdpack_operator_alloc(qs);
    
    operator_ho_lowering(a, qs, 0, 0);
    operator_ho_raising(ad, qs, 0, 0);

    QDPACK_SET_COMPLEX(&z, 0.0, +param->xi);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, ad, ad, QDPACK_COMPLEX_ZERO, adad);
    QDPACK_SET_COMPLEX(&z, 0.0, -param->xi);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z,  a,  a, QDPACK_COMPLEX_ZERO, aa);

    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, adad);
    qdpack_operator_add(H,   aa);
    
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(aa);
    qdpack_operator_free(adad);

    return 0;
}


// -----------------------------------------------------------------------------
// Qubit + Cavity + Transmission line
// -----------------------------------------------------------------------------

/*
 * Generate the hamiltonian for the qubit + 2-HO driving field
 *
 */
int
hamiltonian_q2ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *b, *bd;
    
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);    

    operator_ho_lowering(b, qs, 2, 0);
    operator_ho_raising(bd, qs, 2, 0);

    // a + ad
    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, b);
    qdpack_operator_add(H, bd);

    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    
    return 0;
}


/*
 * Generate the hamiltonian for the Qubit + Cavity + Transmission line system/
 *
 */
int
hamiltonian_q2ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *b, *bd, *ws, *s1;
    qdpack_complex alpha, z;
    
    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    s1 = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    /* qubit */
    operator_sigma_z(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->epsilon[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);

    operator_sigma_x(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->delta[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);
    
    /* h_ w1 (ad a + 1/2) */
    operator_ho_lowering(a, qs, 1, 0);
    operator_ho_raising(ad, qs,1, 0);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ONE, H);
    //operator_unit(ws, qs, -1);
    //QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[0], 0.0);
    //qdpack_operator_scale(ws, z);
    qdpack_operator_add(H, ws);
    
    /* qubit - cavity interaction */
    operator_sigma_z(s1, qs, 0);
    qdpack_operator_memcpy(ws, a);
    qdpack_operator_add(ws, ad);
    QDPACK_SET_COMPLEX(&alpha, param->g_qc, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, s1, ws, QDPACK_COMPLEX_ONE, H);    

    /* h_ w (bd b + 1/2) */
    operator_ho_lowering(b, qs, 2, 0);
    operator_ho_raising(bd, qs, 2, 0);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[1], 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, bd, b, QDPACK_COMPLEX_ONE, H);
    operator_unit(ws, qs, -1);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[1], 0.0);
    qdpack_operator_scale(ws, z);
    qdpack_operator_add(H, ws);
        
    /* cavity - transmission line interaction : g_cr (a + ad) (b + bd) * /
    qdpack_operator_add(a, ad);
    qdpack_operator_add(b, bd);
    QDPACK_SET_COMPLEX(&alpha, param->g_cr, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);
    */
    
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    qdpack_operator_free(ws);
    qdpack_operator_free(s1);

    return 0;
}


// -----------------------------------------------------------------------------
// Qubit + Cavity
//
// -----------------------------------------------------------------------------

/*
 * Generate the hamiltonian for the qubit + HO + driving field
 *
 */
int
hamiltonian_jc_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *b, *bd;
    
    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;

    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    
    operator_ho_lowering(b, qs, 1, param->n_offset);
    operator_ho_raising(bd, qs, 1, param->n_offset);

    // a + ad
    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, b);
    qdpack_operator_add(H, bd);

    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    
    return 0;
}

/*
 * Generate the hamiltonian for the Qubit + Cavity
 *
 */
int
hamiltonian_jc_no_rwa(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *ws, *s1, *n;
    qdpack_complex z, alpha;
    
    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    s1 = qdpack_operator_alloc(qs);
    n  = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    /* qubit */
    operator_sigma_z(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->epsilon[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);

    operator_sigma_x(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->delta[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);
    
    /* h_ w1 ad a */
    // Not necessary
    //operator_unit(ws, qs, -1);
    //QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[0], 0.0);
    //qdpack_operator_scale(ws, z);
    
    operator_ho_N(n, qs, 1, param->n_offset);
    QDPACK_SET_COMPLEX(&z, param->ho_w[0], 0.0);
    qdpack_operator_scale(n, z);
    qdpack_operator_add(ws, n);
    
    qdpack_operator_add(H, ws);

    /* qubit - cavity interaction : without rotating-wave approximation */
    operator_ho_lowering(a, qs, 1, param->n_offset);
    operator_ho_raising(ad, qs, 1, param->n_offset);
    operator_sigma_x(s1, qs, 0);
    //qdpack_operator_memcpy(ws, a);
    //qdpack_operator_add(ws, ad);
    QDPACK_SET_COMPLEX(&alpha, param->g_qc, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, s1, ws, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, s1,  a, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, s1, ad, QDPACK_COMPLEX_ONE, H);

    qdpack_operator_free(ws);
    qdpack_operator_free(s1);    
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(n);

    return 0;
}

//------------------------------------------------------------------------------
// Generate the Jaynes-Cummings Hamiltonian for a Qubit + Cavity system
//
// H = ho_w * a^dag a + 0.5 eps sigma_z + g (a * sigma_+ + a^dag * sigma_-)
//
int
hamiltonian_jc(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *b, *n, *s;
    qdpack_complex alpha;

    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;
    
    a = qdpack_operator_alloc(qs);
    b = qdpack_operator_alloc(qs);
    n = qdpack_operator_alloc(qs);
    s = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    //
    // Qubit 
    //
    if (param->epsilon[0] != 0.0)
    {
        operator_sigma_z(s, qs, 0);
        QDPACK_SET_COMPLEX(&alpha, -0.5 * param->epsilon[0], 0.0);
        qdpack_operator_scale(s, alpha);
        qdpack_operator_add(H, s);
    }

    if (param->delta[0] != 0.0)
    {
        operator_sigma_x(s, qs, 0);
        QDPACK_SET_COMPLEX(&alpha, -0.5 * param->delta[0], 0.0);
        qdpack_operator_scale(s, alpha);
        qdpack_operator_add(H, s);
    }

    //
    // Cavity
    //
    operator_ho_N(n, qs, 1, param->n_offset);
    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
    qdpack_operator_scale(n, alpha);
    qdpack_operator_add(H, n);

    // 
    // Cavity-Qubit interaction
    //
    operator_ho_lowering(a, qs, 1, param->n_offset);
    operator_ho_lowering(b, qs, 0, 0);
    QDPACK_SET_COMPLEX(&alpha, param->g_qc, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackConjTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, alpha, a, b, QDPACK_COMPLEX_ONE, H);

    qdpack_operator_free(s);
    qdpack_operator_free(a);
    qdpack_operator_free(b);
    qdpack_operator_free(n);

    return 0;
}

// -----------------------------------------------------------------------------
// 3-level system + Cavity 
//
// -----------------------------------------------------------------------------

/*
 * Generate the hamiltonian for the 3-level atom + HO + driving field
 *
 */
int
hamiltonian_3ls_ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *b, *bd;

    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;
    
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    
    operator_ho_lowering(b, qs, 1, param->n_offset);
    operator_ho_raising(bd, qs, 1, param->n_offset);

    // a + ad
    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, b);
    qdpack_operator_add(H, bd);

    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    
    return 0;
}

/*
 * Generate the hamiltonian for the Qubit + Cavity
 *
 */
int
hamiltonian_3ls_ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *a, *ad, *ws, *s1, *n;
    qdpack_complex z, alpha;

    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;

    a  = qdpack_operator_alloc(qs);
    ad = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);
    s1 = qdpack_operator_alloc(qs);
    n  = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    /* qubit */
    operator_3ls_sigma_z(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->epsilon[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);
        
    operator_3ls_sigma_x(s1, qs, 0);
    QDPACK_SET_COMPLEX(&alpha, -0.5 * param->delta[0], 0.0);
    qdpack_operator_scale(s1, alpha);
    qdpack_operator_add(H, s1);

    /* h_ w1 (ad a + 1/2) */
    //a  = operator_ho_lowering(qs, 1, param->n_offset);
    //ad = operator_ho_raising(qs, 1, param->n_offset);
    //QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ONE, H);
    operator_unit(ws, qs, -1);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[0], 0.0);
    qdpack_operator_scale(ws, z);
    
    operator_ho_N(n, qs, 1, param->n_offset);
    QDPACK_SET_COMPLEX(&z, param->ho_w[0], 0.0);
    qdpack_operator_scale(n, z);
    qdpack_operator_add(ws, n);
    
    qdpack_operator_add(H, ws);

    /* qubit - cavity interaction */
    operator_ho_lowering(a, qs, 1, param->n_offset);
    operator_ho_raising(ad, qs, 1, param->n_offset);
    operator_3ls_sigma_z(s1, qs, 0);
    qdpack_operator_memcpy(ws, a);
    qdpack_operator_add(ws, ad);
    QDPACK_SET_COMPLEX(&alpha, param->g_qc, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, s1, ws, QDPACK_COMPLEX_ONE, H);

    qdpack_operator_free(ws);
    qdpack_operator_free(s1);
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(n);

    return 0;
}


// -----------------------------------------------------------------------------
// Resonator/cavity
//
// -----------------------------------------------------------------------------


/*
 * Generate the hamiltonian for the HO + driving field
 *
 */
int
hamiltonian_ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *b, *bd;
    
    b  = qdpack_operator_alloc(qs);
    bd = qdpack_operator_alloc(qs);
    
    operator_ho_lowering(b, qs, 0, 0);
    operator_ho_raising(bd, qs, 0, 0);

    // a + ad
    qdpack_operator_set_zero(H);
    qdpack_operator_add(H, b);
    qdpack_operator_add(H, bd);
    
    qdpack_operator_free(b);
    qdpack_operator_free(bd);
    
    return 0;
}

/*
 * Generate the hamiltonian for the harmonic oscillator
 *
 */
int
hamiltonian_ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *ws, *n; //*a, *ad, *ws;
    qdpack_complex z;

    h_jc_parameters_t *param = (h_jc_parameters_t *)sim->parameters;
    
    n  = qdpack_operator_alloc(qs);
    ws = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    /* h_ w1 (ad a + 1/2) */
//        a  = operator_ho_lowering(qs, 1, param->n_offset);
//        ad = operator_ho_raising(qs, 1, param->n_offset);
//        QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
//        qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ONE, H);
    operator_unit(ws, qs, -1);
    QDPACK_SET_COMPLEX(&z, 0.0 * param->ho_w[0], 0.0);
    qdpack_operator_scale(ws, z);

    operator_ho_N(n, qs, 0, param->n_offset);
    QDPACK_SET_COMPLEX(&z, param->ho_w[0], 0.0);
    qdpack_operator_scale(n, z);
    qdpack_operator_add(ws, n);

    qdpack_operator_add(H, ws);
        
    qdpack_operator_free(ws);
    qdpack_operator_free(n);

//    /* h_ w1 (ad a + 1/2) */
//    a  = operator_ho_lowering(qs, 0, 0);
//    ad = operator_ho_raising(qs, 0, 0);
//    QDPACK_SET_COMPLEX(&alpha, param->ho_w[0], 0.0);
//    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, alpha, ad, a, QDPACK_COMPLEX_ONE, H);
//    ws = operator_unit(qs, -1);
//    QDPACK_SET_COMPLEX(&z, 0.5 * param->ho_w[0], 0.0);
//    qdpack_operator_scale(ws, z);
//    qdpack_operator_add(H, ws);
//    qdpack_operator_free(ws);
//
//    qdpack_operator_free(a);
//    qdpack_operator_free(ad);

    return 0;
}

// -----------------------------------------------------------------------------
// 3LS + large spin + harmonic oscillator
// -----------------------------------------------------------------------------

/*
 * Generate the dicke model hamiltonian + 3LS 
 *
 */
int
hamiltonian_dicke(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    qdpack_operator_t *Hsub, *n, *sx, *Jp, *Jm, *J_sum, *a, *ad, *a_sum;
    qdpack_complex z;
    
    h_dicke_parameters_t *param = (h_dicke_parameters_t *)sim->parameters;

    Hsub  = qdpack_operator_alloc(qs);
    n     = qdpack_operator_alloc(qs);
    sx    = qdpack_operator_alloc(qs);
    Jp    = qdpack_operator_alloc(qs);
    Jm    = qdpack_operator_alloc(qs);
    J_sum = qdpack_operator_alloc(qs);
    a     = qdpack_operator_alloc(qs);
    ad    = qdpack_operator_alloc(qs);
    a_sum = qdpack_operator_alloc(qs);

    /* 
      * start by blanking the hamiltonian matrix
      */
    qdpack_operator_set_zero(H);

    /*
     * subsystem 1: 3LS artificial atom: 
     *
     * /
    Hsub = operator_sigma_z(qs, 0);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->epsilon[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);
    qdpack_operator_free(Hsub);

    Hsub = operator_sigma_x(qs, 0);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->delta[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);
    qdpack_operator_free(Hsub);
    */

    /* 
      * subsystem 2: N/2 spin:
      *
      */
    operator_Jz(Hsub, qs, 1-1, param->j);
    QDPACK_SET_COMPLEX(&z, param->J_w[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);

    /* 
     * subsystem 3: HO: h_ w1 (ad a + 1/2) 
     *
     */
    operator_ho_N(Hsub, qs, 2-1, param->n_offset);
    QDPACK_SET_COMPLEX(&z, -param->ho_w[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);

    /*
     * Interaction: 
     * 1) artificial atom <-> oscillator
     * 2) spin <-> oscillator
     *
     */
    operator_ho_lowering(a, qs, 2-1, param->n_offset);
    operator_ho_raising(ad, qs, 2-1, param->n_offset);
    qdpack_operator_memcpy(a_sum, a);
    qdpack_operator_add   (a_sum, ad);
        
    // 1)
    /*
    sx = operator_sigma_x(qs, 0);

    QDPACK_SET_COMPLEX(&z, param->g_ac, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sx, a_sum, QDPACK_COMPLEX_ONE, H);

    qdpack_operator_free(sx);
    */
    // 2)
    
    operator_J_plus (Jp, qs, 1-1, param->j);
    operator_J_minus(Jm, qs, 1-1, param->j);
    qdpack_operator_memcpy(J_sum, Jp);
    qdpack_operator_add   (J_sum, Jm);

    QDPACK_SET_COMPLEX(&z, param->g_sc, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, J_sum, a_sum, QDPACK_COMPLEX_ONE, H);

   
    qdpack_operator_free(Hsub);
    qdpack_operator_free(Jp);
    qdpack_operator_free(Jm);
    qdpack_operator_free(J_sum);
    qdpack_operator_free(a);
    qdpack_operator_free(ad);
    qdpack_operator_free(a_sum);
    qdpack_operator_free(sx);
    qdpack_operator_free(n);

    return 0;
}

//!-----------------------------------------------------------------------------
// @brief   
//
// 
// @param    qs    Data structure that defines the quantum system.
// @param    param    Data structure that contains all the parameters.
//------------------------------------------------------------------------------

/*
 * Generate the dicke model hamiltonian + 3LS 
 *
 */
int
hamiltonian_2ls_largespin(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    int qs_idx_spin, qs_idx_2ls;
    qdpack_operator_t *Hsub, *sx, *sz, *sm, *sp, *Jz, *Jp, *Jm;//, *n *J_sum, *a, *ad, *a_sum;
    qdpack_complex z;

    h_dicke_parameters_t *param = (h_dicke_parameters_t *)sim->parameters;

    printf("%s: DEBUG: Generating 2LS-largespin hamiltonian with total Hilbert space size %d x %d\n", 
           __PRETTY_FUNCTION__, qdpack_hilbert_space_nstates(qs), qdpack_hilbert_space_nstates(qs));
    
    Hsub  = qdpack_operator_alloc(qs);
    //n     = qdpack_operator_alloc(qs);
    sx    = qdpack_operator_alloc(qs);
    sz    = qdpack_operator_alloc(qs);
    sp    = qdpack_operator_alloc(qs);
    sm    = qdpack_operator_alloc(qs);

    Jz    = qdpack_operator_alloc(qs);
    Jp    = qdpack_operator_alloc(qs);
    Jm    = qdpack_operator_alloc(qs);
    //J_sum = qdpack_operator_alloc(qs);
    //a     = qdpack_operator_alloc(qs);
    //ad    = qdpack_operator_alloc(qs);
    //a_sum = qdpack_operator_alloc(qs);

    /* 
      * start by blanking the hamiltonian matrix
      */
    qdpack_operator_set_zero(H);

    /*
     * subsystem 1: 2LS quantum dot
     *
     */
    qs_idx_2ls = 0;

    operator_sigma_z(Hsub, qs, 0);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->epsilon[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);

    operator_sigma_x(Hsub, qs, 0);
    QDPACK_SET_COMPLEX(&z, 0.5 * param->delta[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);

    /* 
     * subsystem 2: N/2 spin:
     *
     */
    qs_idx_spin = 1;
    operator_Jz(Hsub, qs, qs_idx_spin, param->j);
    QDPACK_SET_COMPLEX(&z, param->J_w[0], 0.0);
    qdpack_operator_scale(Hsub, z);
    qdpack_operator_add(H, Hsub);

    /*
     * Interaction: 
     *
     */
       
    // 1)
    operator_sigma_z(sz, qs, qs_idx_2ls);    
    operator_sigma_x(sx, qs, qs_idx_2ls);    

    operator_sigma_plus (sm, qs, qs_idx_2ls);
    operator_sigma_minus(sp, qs, qs_idx_2ls);

    operator_Jz     (Jz, qs, qs_idx_spin, param->j);
    operator_J_plus (Jp, qs, qs_idx_spin, param->j);
    operator_J_minus(Jm, qs, qs_idx_spin, param->j);

    QDPACK_SET_COMPLEX(&z, param->g0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sz, Jz, QDPACK_COMPLEX_ONE, H);

    QDPACK_SET_COMPLEX(&z, param->g1, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sp, Jm, QDPACK_COMPLEX_ONE, H);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sm, Jp, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sx, Jm, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sx, Jp, QDPACK_COMPLEX_ONE, H);

    //qdpack_operator_memcpy(J_sum, Jp);
    //qdpack_operator_add   (J_sum, Jm);
    //QDPACK_SET_COMPLEX(&z, param->g1, 0.0);
    //qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, sx, J_sum, QDPACK_COMPLEX_ONE, H);

    // done   

    qdpack_operator_free(Hsub);
    qdpack_operator_free(Jz);
    qdpack_operator_free(Jp);
    qdpack_operator_free(Jm);
    //qdpack_operator_free(J_sum);
    //qdpack_operator_free(a);
    //qdpack_operator_free(ad);
    //qdpack_operator_free(a_sum);
    qdpack_operator_free(sx);
    qdpack_operator_free(sz);
    qdpack_operator_free(sp);
    qdpack_operator_free(sm);
    //qdpack_operator_free(n);

    return 0;
}



//------------------------------------------------------------------------------
// Generate a multilevel Hamiltonian for a quantum dot system
//
// State mapping:
//
// T(1,1) -> 0
// S(2,0) -> 1
// S(1,1) -> 2
//
//
// H = Delta_ |T(1,1)><T(1,1)| + 
//
//------------------------------------------------------------------------------
int
hamiltonian_mls_qdot(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    int n;
    qdpack_operator_t *op1, *op2;
    qdpack_complex z;
    
    h_dicke_parameters_t *param = (h_dicke_parameters_t *)sim->parameters;

    op1 = qdpack_operator_alloc(qs);
    op2 = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    /* first set the splittings between the states */
    n = 0; // |T(1,1)><T(1,1)| = |0><0|
    operator_project(op1, qs, 0, n, n);
    QDPACK_SET_COMPLEX(&z, param->delta[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);

    n = 1; // |S(1,1)><S(1,1)| = |0><0|
    operator_project(op1, qs, 0, n, n);
    QDPACK_SET_COMPLEX(&z, param->delta[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);

    n = 2; // |S(0,2)><S(0,2)| = |0><0|
    operator_project(op1, qs, 0, n, n);
    QDPACK_SET_COMPLEX(&z, param->delta[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
        
    //
    // setup coherent tunnelling interaction terms
    //
    
    // - |T(1,1)><S(1,1)| + hc = |0><1| + |1><0|
    operator_project(op1, qs, 0, 0, 1);
    operator_project(op2, qs, 0, 1, 0);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[0], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);

    // ELECTRON TUNNELLING WITHOUT SPIN FLIP |S(1,1)><S(0,2)| + hc = |1><2| + |2><1|
    operator_project(op1, qs, 0, 1, 2);
    operator_project(op2, qs, 0, 2, 1);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[1], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);

    // SPIN ORBIT |T(1,1)><S(0,2)| + hc = |0><1| + |1><0|
    operator_project(op1, qs, 0, 0, 2);
    operator_project(op2, qs, 0, 2, 0);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[2], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);



    qdpack_operator_free(op1);
    qdpack_operator_free(op2);

    return 0;
}

//------------------------------------------------------------------------------
//
//
//
//------------------------------------------------------------------------------
int
hamiltonian_mls_qdot_largespin(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H)
{
    int n;
    qdpack_operator_t *op1, *op2, *Jz, *Jm, *Jp;
    qdpack_complex z;
       
    h_dicke_parameters_t *param = (h_dicke_parameters_t *)sim->parameters;

    op1 = qdpack_operator_alloc(qs);
    op2 = qdpack_operator_alloc(qs);
    Jm = qdpack_operator_alloc(qs);
    Jp = qdpack_operator_alloc(qs);
    Jz = qdpack_operator_alloc(qs);

    qdpack_operator_set_zero(H);

    //
    // setup the quantum dot
    //
    //qs_idx_spin = 0;

    /* first set the splittings between the states */
    n = 0; // |T(1,1)><T(1,1)| = |0><0|
    operator_project(op1, qs, 0, 0, 0);
    QDPACK_SET_COMPLEX(&z, param->epsilon[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);

    n = 1; // |S(1,1)><S(1,1)| = |0><0|
    operator_project(op1, qs, 0, n, n);
    QDPACK_SET_COMPLEX(&z, param->epsilon[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);

    n = 2; // |S(0,2)><S(0,2)| = |0><0|
    operator_project(op1, qs, 0, n, n);
    QDPACK_SET_COMPLEX(&z, param->epsilon[n], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
        
    //
    // setup coherent tunnelling interaction terms
    //
    
    // - |T(1,1)><S(1,1)| + hc = |0><1| + |1><0|
    operator_project(op1, qs, 0, 0, 1);
    operator_project(op2, qs, 0, 1, 0);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[0], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);

    // ELECTRON TUNNELLING WITHOUT SPIN FLIP |S(1,1)><S(0,2)| + hc = |1><2| + |2><1|
    operator_project(op1, qs, 0, 1, 2);
    operator_project(op2, qs, 0, 2, 1);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[1], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);

    // SPIN ORBIT |T(1,1)><S(0,2)| + hc = |0><1| + |1><0|
    operator_project(op1, qs, 0, 0, 2);
    operator_project(op2, qs, 0, 2, 0);
    QDPACK_SET_COMPLEX(&z, param->g_interaction[2], 0.0);
    qdpack_operator_scale(op1, z);
    qdpack_operator_add(H, op1);
    qdpack_operator_scale(op2, z);
    qdpack_operator_add(H, op2);


    /* 
     * subsystem 2: N/2 spin:
     *
     */
    //qs_idx_spin = 1;
    operator_Jz(Jz, qs, 1, param->j);
    QDPACK_SET_COMPLEX(&z, param->J_w[0], 0.0);
    qdpack_operator_scale(Jz, z);
    qdpack_operator_add(H, Jz);

    /*
     * Interaction: 
     *
     */
       
    // 1)
    operator_project(op1, qs, 0, 0, 0);
    operator_Jz     (Jz,  qs, 1, param->j);
    QDPACK_SET_COMPLEX(&z, param->g0, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, op1, Jz, QDPACK_COMPLEX_ONE, H);

    operator_project(op1, qs, 0, 0, 1); // create
    operator_project(op2, qs, 0, 1, 0); // annihilate
    operator_J_plus (Jp,  qs, 1, param->j);
    operator_J_minus(Jm,  qs, 1, param->j);
    QDPACK_SET_COMPLEX(&z, param->g1, 0.0);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, op1, Jm, QDPACK_COMPLEX_ONE, H);
    qdpack_operator_blas_zgemm(QDpackNoTrans, QDpackNoTrans, z, op2, Jp, QDPACK_COMPLEX_ONE, H);

    // done   

    qdpack_operator_free(Jz);
    qdpack_operator_free(Jp);
    qdpack_operator_free(Jm);
    qdpack_operator_free(op1);
    qdpack_operator_free(op2);

    return 0;
}


