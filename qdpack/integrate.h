//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for integration of systems of differential equations. 
 */

#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <qdpack/qdpack.h>

#define RK_REC_MAX_LEVEL 15

//
// operator integration
//
typedef int (*ode_operator_func_d_dt)(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_operator_t *rho_t, qdpack_operator_t *drho_dt, double t);

int odeint_operator_rk(qdpack_simulation_t *sim, 
                       qdpack_hilbert_space_t *qs, 
                       qdpack_operator_t *rho0, 
                       ode_operator_func_d_dt f_d_dt, 
                       operator_cb_func_t rho_cb_func);

int odeint_operator_rk_adaptive(qdpack_simulation_t *sim, 
                                qdpack_hilbert_space_t *qs, 
                                qdpack_operator_t *rho0, 
                                ode_operator_func_d_dt f_d_dt, 
                                operator_cb_func_t rho_cb_func);

int odeint_operator_rk_adaptive_rec(qdpack_simulation_t *sim, 
                                    qdpack_hilbert_space_t *qs, 
                                    ode_operator_func_d_dt f_d_dt,
                                    qdpack_operator_t *rho0, 
                                    qdpack_operator_t *rho_final,
                                    int level, double ti, double tf, double etol, double *e);

double operator_max_deviation(qdpack_operator_t *m1, qdpack_operator_t *m2, double etol);
    

//
// state vector integration
//
typedef int (*ode_state_func_d_dt)(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, qdpack_state_t *psi_t, qdpack_state_t *dpsi_dt, double t);

int odeint_state_rk(qdpack_simulation_t *sim, 
                    qdpack_hilbert_space_t *qs,
                    qdpack_state_t *psi0, 
                    ode_state_func_d_dt f_d_dt, 
                    state_cb_func_t psi_cb_func);

int odeint_state_rk_adaptive(qdpack_simulation_t *sim,
                             qdpack_hilbert_space_t *qs,
                             qdpack_state_t *psi0, 
                             ode_state_func_d_dt f_d_dt, 
                             state_cb_func_t psi_cb_func);

int odeint_state_rk_adaptive_rec(qdpack_simulation_t *sim, 
                                 qdpack_hilbert_space_t *qs, 
                                 ode_state_func_d_dt f_d_dt,
                                 qdpack_state_t *psi0, 
                                 qdpack_state_t *psi_final,
                                 int level, double ti, double tf, double etol, double *e);

double state_max_deviation(qdpack_state_t *m1, qdpack_state_t *m2, double etol);
    
#endif

