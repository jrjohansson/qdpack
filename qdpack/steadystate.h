//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for calculating the steady state of a system with dissipation. 
 */

#ifndef STEADY_STATE_H
#define STEADY_STATE_H

#include <gsl/gsl_matrix.h>

#include "hilbert_space.h"
#include "simulation.h"

int qdpack_steadystate_dm(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss);
int qdpack_steadystate_dm_eig(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss);
int qdpack_steadystate_dm_propagator(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *U, qdpack_operator_t *rho_ss);

#endif

