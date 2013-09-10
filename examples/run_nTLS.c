//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Simulate the time evolution of a system of N two-level systems (TLS).
// 

#include <math.h>
#include <stdio.h> 
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <qdpack/qdpack.h>


/** 
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.). This is a callback function that the solver
 * call in each time-step.
 */
int
rho_store_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    double c;
    int i;
    char filename[128], row[512];
    qdpack_operator_t *dm_part;
    
        for (i = 0; i < qs->nsubsys; i++)
        {
        snprintf(filename, sizeof(filename), "data.qubits/tls_%d_dynamics.dat", i);
                
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f", t, 
                 QDPACK_REAL(qdpack_operator_get(dm_part, 0, 0)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 0, 1)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 1, 0)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1)));
        
        qdpack_file_row_add(filename, row);
                qdpack_operator_free(dm_part);
        }

    /* --- concurrence -- */
    snprintf(filename, sizeof(filename), "data.qubits/concurrence.dat");
    c = qdpack_dm_concurrence(qs, rho_t);
    snprintf(row, sizeof(row), "%f\t%f", t, c); 
    qdpack_file_row_add(filename, row);
    
    return 0;
}


/* ---
 * Two-level system parameters: 
 *
 * H = sum_i(-0.5 espilon[i] sigma_z(i) - 0.5 delta[i] sigma_x(i)) + sum_ij(- 0.5 lambda[i][j] sigma_x[i] sigma_x[j]) 
 * 
 */
double epsilon[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double *lambda[10];


/* ---
 * Program starts here.
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho1, *rho2;
    qdpack_operator_t *rho0;
    qdpack_operator_list_t dm_list;
    int qsn, i;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    /*
     * Setup solver parameters.
     *
     */

    // TLS hamiltonian parameters
    param.epsilon = epsilon;
    param.delta   = delta;
    for (i = 0; i < 10; i++)
    {
        param.epsilon[i] *= 2*M_PI;
        param.delta[i]   *= 2*M_PI;
        lambda[i] = (double *)calloc(10, sizeof(double));
    }
    lambda[0][1] = 0.25 * 2*M_PI;
    lambda[1][2] = 0.4 * 2*M_PI;
    lambda[2][3] = 0.3 * 2*M_PI;
    lambda[3][4] = 0.2 * 2*M_PI;
    lambda[4][5] = 0.1 * 2*M_PI;
    param.lambda  = lambda;
    
    
    // time range
    param.Ti =  0.0;
    param.Tf =  10.0;
    param.dT =  0.01;

    param.N = 6;

    // driving field (set to zero to disable)
    param.h_td_A = 0.0 * 2 * M_PI;
    param.h_td_w = 0.0 * 2 * M_PI;


    // function pointer to the functions generating the hamiltonian matrices
    param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    param.H1_func = (hamiltonian_func_t)hamiltonian_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_t;

    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();

    for (i = 0; i < param.N; i++)
        qdpack_hilbert_space_add(qs, 2);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);


    /*
     * Set up dissipation
     *
     */
    param.do_n = 0;
    param.do_a[param.do_n] = operator_ho_lowering(qs, 0); // qubit 0
    param.do_g1[param.do_n] = 0.50; // relaxation
    param.do_g2[param.do_n] = 0.05; // excitation
    param.do_n++;
    param.do_a[param.do_n] = operator_ho_lowering(qs, 1); // qubit 1
    param.do_g1[param.do_n] = 0.25; // relaxation
    param.do_g2[param.do_n] = 0.01; // excitation
    param.do_n++;
    param.do_a[param.do_n] = operator_sigma_z(qs, 0); // qubit 0
    param.do_g1[param.do_n] = 0.25; // dephasing
    param.do_n++;


    /*
     * debug
     *
     */
    printf("sigma_x(0) = \n"); qdpack_operator_print(operator_sigma_x(qs,0), 4);
    printf("sigma_x(0) = \n"); qdpack_operator_print(operator_sigma_x(qs,1), 4);
    printf("sigma_y(0) = \n"); qdpack_operator_print(operator_sigma_y(qs,0), 4);
    printf("sigma_y(0) = \n"); qdpack_operator_print(operator_sigma_y(qs,1), 4);

    
    /*
     * Setup initial state
     */
    rho1 = qdpack_dm_pure_TLS(0.8); // 0.8|0>+0.2|1>
    rho2 = qdpack_dm_pure_TLS(0.3); // 0.3|0>+0.7|1> 
    
        qdpack_operator_list_init(&dm_list);
        qdpack_operator_list_append(&dm_list, rho1);
        for (i = 1; i < param.N; i++)
        qdpack_operator_list_append(&dm_list, rho2);

    rho0 = qdpack_operator_tensor(qs, &dm_list);
    
    //printf("rho0 =\n");
    //qdpack_operator_print(rho0, qsn);


    /*
     * Evolve the quantum system.
     * 
     */
    //if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    return 0;
}



