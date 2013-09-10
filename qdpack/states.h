//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for combining density matrices and for tracing out sub-systems from
 * full-system density matrices.
 */

#ifndef STATES_H
#define STATES_H

#include <qdpack/qdpack.h>

// density matrix

qdpack_operator_t *qdpack_dm_pure_TLS(double p_ex);
qdpack_operator_t *qdpack_dm_fock_state(int n, int N);
qdpack_operator_t *qdpack_dm_coherent_state(double r, double theta, int N, int offset);
qdpack_operator_t *qdpack_dm_boson_thermal(double w, double w_th, int N);
qdpack_operator_t *qdpack_dm_uniform_superposition(int N);
qdpack_operator_t *qdpack_dm_uniform_mixture(int N);


// wave function

qdpack_state_t *qdpack_wf_fock_state(int n, int N);
qdpack_state_t *qdpack_wf_coherent_state(double r, double theta, int N, int offset);
qdpack_state_t *qdpack_wf_pure_TLS(double p_ex);

qdpack_state_t *qdpack_wf_superposition(qdpack_state_t *wf1, qdpack_state_t *wf2);

#endif

