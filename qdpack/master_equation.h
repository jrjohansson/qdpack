//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for generating hamiltonians in matrix form 
 */

#ifndef MASTER_EQUATION_H
#define MASTER_EQUATION_H

#include "hilbert_space.h"
#include "simulation.h"

// spectral density callback function
typedef double (*spectral_density_cb_t)(qdpack_simulation_t *sp, double f);

// density matrix solvers
typedef int (*operator_cb_func_t)(qdpack_simulation_t *sp, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_t, double t);

int    qdpack_evolve_dm_unitary_const       (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho0, operator_cb_func_t rho_cb_func);
int    qdpack_evolve_dm_unitary_const_simple(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho0, operator_cb_func_t rho_cb_func);
int    qdpack_evolve_dm_unitary             (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho0, operator_cb_func_t rho_cb_func);

int    qdpack_evolve_dm_lme                 (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho0, operator_cb_func_t rho_cb_func);
int    qdpack_evolve_dm_lme_ib              (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho0, operator_cb_func_t rho_cb_func);

int    qdpack_propagator_dm       (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *U, operator_cb_func_t rho_cb_func);

// wave function solvers
typedef int (*state_cb_func_t)(qdpack_simulation_t *sp, qdpack_hilbert_space_t *qs, qdpack_state_t *psi_t, double t);

int qdpack_evolve_wf    (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_state_t *psi0, state_cb_func_t psi_cb_func);

int qdpack_propagator   (qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *U, state_cb_func_t psi_cb_func);

// monte carlo
int qdpack_evolve_wfmc(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_state_t *psi0, state_cb_func_t psi_cb_func);

#endif

