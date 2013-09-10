//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for calculating basis transformation matrices, etc.
 * 
 */

#ifndef BASIS_TRANSFORM_H
#define BASIS_TRANSFORM_H

#include "hilbert_space.h"
#include "simulation.h"

//typedef double (*ho_w_func_t)(double t, qdpack_simulation_t *sp);

qdpack_operator_t *basis_transform(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, int n, ho_w_func_t ho_w_cb, double t);
qdpack_operator_t *basis_transform_derivative(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, int n, ho_w_func_t ho_w_cb, double t);


#endif

