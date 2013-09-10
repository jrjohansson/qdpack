//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Simulate an adiabatic quantum computation.
// 

#include <math.h>
#include <stdio.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <qdpack/qdpack.h>


#define DATADIR "data/adiabatic"

/* ---
 * This function defines the time-dependence of the hamiltonian.
 *
 */
int
hamiltonian_sweep_t(qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t, qdpack_simulation_t *p)
{
    double lambda = t;
        qdpack_complex z;

    /*
     *  H(t) = (lambda - 1) H0 + lambda H1
     */

    if (lambda == 0.0)
    {
            qdpack_operator_memcpy(H_t, H0);
    }
    else
    {
        /*
         * First calculate: H = (lambda-1)/lambda H0
         */
            qdpack_operator_memcpy(H_t, H0);
            QDPACK_SET_COMPLEX(&z, (1.0 - lambda)/lambda, 0.0);
            qdpack_operator_scale(H_t, z);
        
        /*
         * Next, add H1:    H = (lambda-1)/lambda H0 + H1
         */
        qdpack_operator_add(H_t, H1);

        /*
         * Finally, multiply with lambda:
         *
         * H = lambda ((1 - lambda)/lambda H0 + H1) = (1 - lambda) H0 + lambda H1 = H(t)
         */
            QDPACK_SET_COMPLEX(&z, lambda, 0.0);
            qdpack_operator_scale(H_t, z);
    }

    /* 
     * Diagonalize the hamiltonian and save eigenvalues to file:
     */
    {
        qdpack_operator_t *M, *evec;
        gsl_vector *eval;
        gsl_eigen_hermv_workspace *w;
        char filename[256], row[1024];
        int qsn, ri, j;

        qsn  = qdpack_hilbert_space_nstates(p->qs);
            M    = qdpack_operator_alloc(qsn, qsn);
            evec = qdpack_operator_alloc(qsn, qsn);
            eval = gsl_vector_alloc(qsn);
            w    = gsl_eigen_hermv_alloc(qsn);

            qdpack_operator_memcpy(M, H_t);
            gsl_eigen_hermv(M, eval, evec, w);
            gsl_eigen_hermv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

        ri = 0;
                snprintf(filename, sizeof(filename), "%s/H_evals%s.dat", DATADIR, p->simsig);
            ri = snprintf(row, sizeof(row), "%f", t);
            for (j = 0; j < qsn; j++)
            {
                  ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", gsl_vector_get(eval, j));
            }
            qdpack_file_row_add(filename, row);

        gsl_vector_free(eval);
        qdpack_operator_free(M);
        qdpack_operator_free(evec);

    }

        return 0;
}

/*
 * Reverse sigma_z coefficients in qubits hamiltonian
 *
 */
qdpack_operator_t *
hamiltonian_qubits_reverse(qdpack_hilbert_space_t *qs, qdpack_simulation_t *param)
{
    int i;
    qdpack_operator_t *H;

    /* reverse sigma_z coefficients */
    for (i = 0; i < 10; i++)
        param->epsilon[i] *= -1.0;        

    H = hamiltonian_qubits(qs, param);

    /* restore sigma_z coefficients */
    for (i = 0; i < 10; i++)
        param->epsilon[i] *= -1.0;        

    return H;
}

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
    char filename[256], row[512];
    qdpack_operator_t *dm_part;
    
        for (i = 0; i < qs->nsubsys; i++)
        {
        snprintf(filename, sizeof(filename), "%s/tls_%d_dynamics%s.dat", DATADIR, i, sp->simsig);
                
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f", t, 
                 QDPACK_REAL(qdpack_operator_get(dm_part, 0, 0)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 0, 1)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 1, 0)),
                 QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1)));
        
        qdpack_file_row_add(filename, row);
                qdpack_operator_free(dm_part);
        }

    return 0;
}


/* ---
 * Two-level system parameters: 
 *
 * H = sum_i(-0.5 espilon[i] sigma_z(i) - 0.5 delta[i] sigma_x(i)) + sum_ij(- 0.5 lambda[i][j] sigma_x[i] sigma_x[j]) 
 * 
 */
double epsilon[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double delta[]   = {100.0, 100.0, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1};
double *lambda[10];


/* ---
 * Program starts here.
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho1;
    qdpack_operator_t *rho0;
    qdpack_operator_list_t dm_list;
    int qsn, i;
    double Tf;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    /*
     * Setup solver parameters.
     *
     */

    Tf = 30 * delta[0] / epsilon[0];

    // TLS hamiltonian parameters
    param.epsilon = epsilon;
    param.delta   = delta;
    for (i = 0; i < 10; i++)
    {
        param.delta[i]   *= 2*M_PI;
        param.epsilon[i] *= 2*M_PI * Tf;
        lambda[i] = (double *)calloc(10, sizeof(double));
    }
    lambda[0][1] = 0.2 * 2*M_PI;
    lambda[1][2] = 0.2 * 2*M_PI;
    lambda[2][3] = 0.2 * 2*M_PI;
    lambda[3][4] = 0.2 * 2*M_PI;
    lambda[4][5] = 0.2 * 2*M_PI;
    param.lambda  = lambda;
    
    
    // time range
    param.Ti =  0.0;
    param.Tf =  1.0;
    param.dT =  0.005;

    param.N = 2;

    // driving field (set to zero to disable)
    param.h_td_A = 0.0 * 2 * M_PI;
    param.h_td_w = 0.0 * 2 * M_PI;


    // function pointer to the functions generating the hamiltonian matrices
    param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    param.H1_func = (hamiltonian_func_t)hamiltonian_qubits_reverse;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_sweep_t;

    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();

    for (i = 0; i < param.N; i++)
        qdpack_hilbert_space_add(qs, 2);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;

    /*
     * Set up dissipation
     *
     * /
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
    */

    /*
     * Setup initial state
     */
    rho1 = qdpack_dm_pure_TLS(0.0); // 0.8|0>+0.2|1>
    
        qdpack_operator_list_init(&dm_list);
        qdpack_operator_list_append(&dm_list, rho1);
        for (i = 1; i < param.N; i++)
        qdpack_operator_list_append(&dm_list, rho1);

    rho0 = qdpack_operator_tensor(qs, &dm_list);
    
    //printf("rho0 =\n");
    //qdpack_operator_print(rho0, qsn);

    /*
     * Create a simulation signature that will be appended as
     * an identifier on all output file names.
     */
    snprintf(param.simsig, sizeof(param.simsig), "_N_%d_V_%.3f_Delta_%.3f", param.N, param.epsilon[0] / (2*M_PI), param.delta[0] / (2*M_PI));

    /*
     * Evolve the quantum system.
     * 
     */
    //if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    return 0;
}



