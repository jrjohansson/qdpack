//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for calculating quantifying entanglement. 
 */

#ifndef ENTANGLEMENT_H
#define ENTANGLEMENT_H

#include <gsl/gsl_matrix.h>

#include "qdpack.h"

double qdpack_entanglement_concurrence(qdpack_hilbert_space_t *qs, qdpack_operator_t *rho);
double qdpack_entanglement_log_neg(qdpack_operator_t *rho);
double qdpack_entanglement_neumann_entropy(qdpack_operator_t *rho);
    

#endif

