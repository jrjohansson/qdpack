//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Implementation of matrix-exponentiation for complex matrices.
 */

#ifndef GSL_EXPM_COMPLEX_H
#define GSL_EXPM_COMPLEX_H

#include <gsl/gsl_mode.h>
#include <gsl/gsl_matrix.h>

int gsl_linalg_complex_exponential_ss(const gsl_matrix_complex *A, gsl_matrix_complex *eA, gsl_mode_t mode);

#endif 
