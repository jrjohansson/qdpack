//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @mainpage Quantum Dynamics Package
 *
 * A C library for time-evolution of quantum systems.
 * 
 * @author J Robert Johansson
 * 
 */

/**
 * @file simulation.h
 *
 * @brief Data structure that contains all the parameters and settings for a
 * simulation.
 * 
 * @author J Robert Johansson <robert@riken.jp>
 *
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include <math.h>
#include <string.h>

#include <qdpack/qdpack.h>

struct _qdpack_simulation;

typedef int(*hamiltonian_func_t)(struct _qdpack_simulation *, qdpack_hilbert_space_t *, qdpack_operator_t *);
typedef int(*hamiltonian_td_func_t)(struct _qdpack_simulation *, qdpack_operator_t *, qdpack_operator_t *, qdpack_operator_t *, double);
typedef double(*ho_w_func_t)(double, void*);

#define H_MAX_TERMS 20

typedef struct _qdpack_simulation {

    qdpack_hilbert_space_t *qs;
    
    /* -- time evolution -- */
    double Ti;
    double Tf;
    double dT;

    //double x, y;

    char simsig[1024];
    
    /* --- Temperature in frequency unit --- */
    //double wT;
    double T_w;

    qdpack_operator_t *H0;
    qdpack_operator_t *H1;
    qdpack_operator_t *H_t;
    qdpack_operator_t *H_a;
    qdpack_operator_t *H_b;


    /* -- dissipation -- */
    qdpack_operator_t *do_a[H_MAX_TERMS];
    qdpack_operator_t *do_ad[H_MAX_TERMS]; //XXX only used in steady state solver
    qdpack_operator_t *do_aad[H_MAX_TERMS];
    qdpack_operator_t *do_ada[H_MAX_TERMS];

    qdpack_operator_t *do_a_eb[H_MAX_TERMS];
    qdpack_operator_t *do_ad_eb[H_MAX_TERMS];
    qdpack_operator_t *do_aad_eb[H_MAX_TERMS];
    qdpack_operator_t *do_ada_eb[H_MAX_TERMS];

    qdpack_operator_t *ws1;
    qdpack_operator_t *ws2;
    qdpack_operator_t *op_tmp;
    
    double do_g1[H_MAX_TERMS];
    double do_g2[H_MAX_TERMS];

    int do_n;

    /* --- problem dependent parameters --- */

    void *parameters;

    /* --- continuous basis transform --- */
    qdpack_operator_t *dSSt;
    ho_w_func_t ho_w_func;

    /* --- options --- */
    int cont_basis_transform;
    int option_me_eb;
    int option_adaptive;
    
    /* --- Harmonic oscillator --- */

    double g0, g1, g2;

    qdpack_operator_t *N_op[10], *N2_op[10];
    qdpack_operator_t *a, *ad, *b, *bd, *b2, *bd2, *rho_f, *A, *B;
    qdpack_operator_t *n_op, *n_op_2, *x1, *x1_2, *x2, *x2_2;

    qdpack_state_t *psi_f;

    /* -- driving field -- */
    
    double h_td_A;
    double h_td_w;
    //double phi;

    //double Nexpt[2];
    
    /* --- hamiltonians --- */
    hamiltonian_func_t H0_func;
    hamiltonian_func_t H1_func;
    hamiltonian_td_func_t Ht_func;


    /* --- eigenvalues --- */
    qdpack_operator_t *evec;
    qdpack_matrix_t *eval;
    
    /* --- floquet state --- */
    
    qdpack_operator_t *floquet_states_t;
    qdpack_operator_t *floquet_modes, *floquet_modes_t;
    qdpack_matrix_t   *floquet_quasienergies;

    qdpack_operator_t **floquet_X, **floquet_gamma, **floquet_delta, *floquet_A, *floquet_operator;
    qdpack_matrix_t   *floquet_Gamma;
    qdpack_matrix_t   *floquet_ss_prob;

    int floquet_k_max;

} qdpack_simulation_t;


void qdpack_simulation_init(qdpack_simulation_t *sp);
void qdpack_simulation_free(qdpack_simulation_t *sp);

#endif
