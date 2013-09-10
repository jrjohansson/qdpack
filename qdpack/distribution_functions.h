//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------


#ifndef DIST_FUNCS_H
#define DIST_FUNCS_H

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include "qdpack.h"

int distribution_function_Q(qdpack_matrix_t *Q, qdpack_operator_t *rho, double alpha_max);


int distribution_function_characteristic_w(
            qdpack_matrix_t *phase_space, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex lambda_max);

double cb_integrand_wigner_cf(qdpack_complex lambda, qdpack_matrix_t *phase_space_cf, qdpack_complex alpha); 

int distribution_function_wigner_cf_quad(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max);

int distribution_function_wigner_cf_sum(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max);

int distribution_function_wigner_cf_fft(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max);

typedef struct _wigner_hermite_data {
    double p;
    double q;
    qdpack_operator_t *rho;
    int N;
} wigner_hermite_data_t;

typedef struct _wigner_data {
    double p;
    double q;
 
    qdpack_matrix_t *Cw;

    qdpack_complex lambda, alpha;
    qdpack_complex lambda_max, alpha_max;

    qdpack_operator_t *rho;

    gsl_function f_re, f_im;

    gsl_integration_workspace  *w_re, *w_im;
    gsl_integration_qawo_table *wf_re, *wf_im;

    int N;


} wigner_data_t;

double cb_integrand_wigner_hermite(double x, void *data);
int distribution_function_wigner_hermite_quad(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *rho,
            double x_limit);

#endif

