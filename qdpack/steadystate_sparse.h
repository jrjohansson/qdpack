//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for calculating the steady state of a system with dissipation. 
 * Using sparse matrix algebra.
 */

#ifndef STEADY_STATE_SPARSE_H
#define STEADY_STATE_SPARSE_H

#include <gsl/gsl_matrix.h>

#include "hilbert_space.h"
#include "simulation.h"

int    qdpack_steadystate_dm_sparse(qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_ss, qdpack_simulation_t *param);
int    qdpack_steadystate_dm_propagator(qdpack_hilbert_space_t *qs, qdpack_simulation_t *param, qdpack_operator_t *rho_ss, qdpack_operator_t *U);

#endif

