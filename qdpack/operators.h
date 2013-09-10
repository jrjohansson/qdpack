//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for generating matrix representations of quantum operators.
 */

#ifndef OPERATORS_H
#define OPERATORS_H

#include <qdpack/qdpack.h>

// Qubit / 2LS operators
int operator_sigma_z(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_sigma_x(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_sigma_y(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_sigma_plus(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n);
int operator_sigma_minus(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n);

// 3LS operators
int operator_3ls_sigma_z(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_3ls_sigma_x(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);


// DQD 3-level operators
int operator_dqd_sigma_z(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_dqd_sigma_x(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);
int operator_dqd_sigma_y(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);


// large spin operators
int operator_J_minus(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, double jj);
int operator_J_plus (qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, double jj);
int operator_Jz     (qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, double jj);

// harmonic oscillator operators
int operator_ho_lowering(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, int offset);
int operator_ho_raising(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, int offset);
int operator_ho_N(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, int offset);

// generic operators
int operator_project(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, int si, int sj);
int operator_unit(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n);

qdpack_complex operator_trace(qdpack_operator_t *op);


#endif

