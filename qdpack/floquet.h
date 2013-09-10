//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for calculating floquet modes and evolution.
 */

#ifndef MASTER_EQUATION_H
#define MASTER_EQUATION_H

#include <qdpack/qdpack.h>

int qdpack_floquet_modes(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs);
int qdpack_floquet_modes_t(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, double t);
int qdpack_floquet_dissipation_rates(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, operator_cb_func_t floquet_state_store_cb, spectral_density_cb_t s_func);

#endif

