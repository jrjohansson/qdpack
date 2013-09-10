//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file hamiltonian.h
 *
 * @brief Functions for generating hamiltonians in matrix form. 
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include <qdpack/qdpack.h>

#define H_MAX_TERMS 20

//------------------------------------------------------------------------------
// Generic time-dependence
//------------------------------------------------------------------------------

//typedef  int(*hamiltonian_func_real_t)(qdpack_hilbert_space_t *qs, qdpack_simulation_t *sim, qdpack_operator_t *H);

int hamiltonian_t(qdpack_simulation_t *sim, qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t);

//------------------------------------------------------------------------------
// coupled two-level systems
//------------------------------------------------------------------------------

typedef struct {

    /* -- Hamiltonian -- */
    int N;

    double epsilon[H_MAX_TERMS];
    double delta[H_MAX_TERMS];
    //double g_interaction[H_MAX_TERMS];
    double lambda[H_MAX_TERMS][H_MAX_TERMS];

} h_qubits_parameters_t;

int hamiltonian_qubits(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_z_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_x_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_spinchain(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// coupled oscillators
//------------------------------------------------------------------------------

typedef struct {

    int N, n_offset;

    double ho_w[H_MAX_TERMS];
    double ho_g;
    double ho_g2;
    double ho_wd;
    double g_qc;
    double g_cr;

    double xi;

} h_ho_parameters_t;

int hamiltonian_2ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_2ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_2ho_optmech(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_2ho_optmech_linearized(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_2ho_cooling(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// Jaynes-Cumming style hamiltonians (qubit coupled to oscillator mode) 
//------------------------------------------------------------------------------

typedef struct {

    int N, n_offset;

    double ho_w[H_MAX_TERMS];
    double ho_g;
    double ho_g2;
    double ho_wd;
    double g_qc;
    double g_cr;

    double epsilon[H_MAX_TERMS];
    double delta[H_MAX_TERMS];

} h_jc_parameters_t;

int hamiltonian_jc(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_jc_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_q2ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_q2ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_3ls_ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_3ls_ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_ho(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_ho_drive(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------

typedef struct {

    double j;
    double *J_w;
    double g_ac, g_sc;

    int N, n_offset;

    double ho_w[H_MAX_TERMS];
    double ho_g;
    double ho_g2;
    double ho_wd;
    double g_qc;
    double g_cr, g0, g1;

    //int N;

    double epsilon[H_MAX_TERMS];
    double delta[H_MAX_TERMS];
    double g_interaction[H_MAX_TERMS];
    //double lambda[H_MAX_TERMS][H_MAX_TERMS];

} h_dicke_parameters_t;

int hamiltonian_dicke(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

int hamiltonian_2ls_largespin(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
int hamiltonian_squeezing(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// 
//------------------------------------------------------------------------------
//int hamiltonian_jc(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

//------------------------------------------------------------------------------
// multi-level quantum model
//------------------------------------------------------------------------------
int hamiltonian_mls_qdot(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);
int hamiltonian_mls_qdot_largespin(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *H);

#endif

