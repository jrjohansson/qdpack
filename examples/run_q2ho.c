//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//
//------------------------------------------------------------------------------

#include <math.h>
#include <stdio.h> 
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <qdpack/qdpack.h>


/* ---
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.) 
 */
int
rho_store_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    int i, j;
    qdpack_complex N_expt; //, op_expt;
    char filename[128], row[1024];
    qdpack_operator_t *dm_part;
    
        for (i = 0; i < qs->nsubsys; i++)
        {
        int ri = 0;    
        
        snprintf(filename, sizeof(filename), "data.q2ho/rho_ho_%d_diag%s.dat", i, sp->simsig);
                
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        ri = snprintf(row, sizeof(row), "%f", t);

        for (j = 0; j < qs->nstates[i]; j++)
             ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", QDPACK_REAL(qdpack_operator_get(dm_part, j, j)));
    
        qdpack_file_row_add(filename, row);


        /* -- expectation values -- */
        N_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[i]);
        snprintf(filename, sizeof(filename), "data.q2ho/N_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
        qdpack_file_row_add(filename, row);        
        
/*
        op_expt = qdpack_operator_expectation_value(rho_t, sp->a);
        snprintf(filename, sizeof(filename), "data.q2ho/a_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->ad);
        snprintf(filename, sizeof(filename), "data.q2ho/ad_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b);
        snprintf(filename, sizeof(filename), "data.q2ho/b_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd);
        snprintf(filename, sizeof(filename), "data.q2ho/bd_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b2);
        snprintf(filename, sizeof(filename), "data.q2ho/b2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd2);
        snprintf(filename, sizeof(filename), "data.q2ho/bd2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        
*/
        
        qdpack_operator_free(dm_part);
        }
    
    return 0;
}


/* ---
 * Simulation parameters.
 * 
 */
double epsilon[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ho_w[]    = { 9.899, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double *lambda[10];

double Ej0 = 13.7;
double g0 = 80e-3;


/* ---
 * Program starts here.
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho_q, *rho_c, *rho_r;
    qdpack_operator_t *rho0;
    qdpack_operator_list_t dm_list;
    int qsn, i;
    double Navg, fw, fA, phi;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);


    if (argc != 1)
    {
        printf("usage: %s check code\n", argv[0]);
        exit(0);
    }
    
    phi     = strtod(argv[1], NULL);

    printf("Input arguments: phi = %f\n", phi);
    
    /*
     * Setup solver parameters.
     *
     */
    
    // the qubit
    epsilon[0] = 0;
    delta[0] = Ej0 * cos(fabs(M_PI * phi));

    printf("Qubit parameters: eps = %f, delta = %f\n", epsilon[0], delta[0]);
    
    param.epsilon = epsilon;
    param.delta   = delta;
    simulation.ho_w    = ho_w;
    for (i = 0; i < 10; i++)
    {
        param.epsilon[i] *= 2*M_PI;
        param.delta[i]   *= 2*M_PI;
        simulation.ho_w[i]    *= 2*M_PI;
        lambda[i] = (double *)calloc(10, sizeof(double));
    }
    lambda[0][1] = 0.0 * 2*M_PI;
    
    param.lambda  = lambda;

    param.g_qc = g0 * 2*M_PI; // 0.01
    param.g_cr = 0.0 * 2*M_PI; // 0.01
    
    param.Ti =  0.0;
    param.Tf =  2.0;
    param.dT =  0.02;

    param.N = 2;
    fA = 0.0;
    fw = 0.0;
    param.h_td_A = fA * 2 * M_PI;
    param.h_td_w = fw * ho_w[1];

    //param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    param.H0_func = (hamiltonian_func_t)hamiltonian_q2ho;
    param.H1_func = (hamiltonian_func_t)hamiltonian_q2ho_drive;
    
    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_add(qs, 20);
    qdpack_hilbert_space_add(qs, 1);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;


    /*
     * Set up dissipation
     *
     */
    param.wT = 0.0 * 2 * M_PI;
    param.do_n = 0;

    /* --- qubit --- * /
    param.do_a[param.do_n] = operator_ho_lowering(qs, 0); 
    Navg = 1/(exp(sqrt(param.epsilon[0]*param.epsilon[0] + param.delta[0]*param.delta[0])/param.wT) - 1);
    param.do_g1[param.do_n] = 0.075 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 0.075 * Navg;       // excitation
    param.do_n++;

    /* --- cavity --- * /
     param.do_a[param.do_n] = operator_ho_lowering(qs, 1); 
    Navg = 1/(exp(simulation.ho_w[0]/param.wT) - 1);
    param.do_g1[param.do_n] = 0.05 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 0.05 * Navg;       // excitation
    param.do_n++;
    
    /* --- transmission line --- * /
     param.do_a[param.do_n] = operator_ho_lowering(qs, 2); 
    Navg = 1/(exp(simulation.ho_w[1]/param.wT) - 1);
    param.do_g1[param.do_n] = 0.05 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 0.05 * Navg;       // excitation
    param.do_n++;
    */

    /* 
     * Pre-calc of operators to calculate expectation values for.
     */
    param.N_op[0] = operator_ho_N(qs, 0);
    param.N_op[1] = operator_ho_N(qs, 1);
    param.N_op[2] = operator_ho_N(qs, 2);
/*
    param.a  = operator_ho_lowering(qs, 0);
    param.ad = operator_ho_raising(qs, 0);
    param.b  = operator_ho_lowering(qs, 1);
    param.bd = operator_ho_raising(qs, 1);
    param.b2  = qdpack_operator_alloc(qsn, qsn);
    param.bd2  = qdpack_operator_alloc(qsn, qsn);
    {
        qdpack_complex alpha, beta;
        QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
        QDPACK_SET_COMPLEX(&beta, 0.0, 0.0);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, param.b,  param.b,  beta, param.b2);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, param.bd, param.bd, beta, param.bd2);
    }
*/
        
    snprintf(param.simsig, sizeof(param.simsig), "_NEC_Nq_%d_Nc_%d_Nr_%d_wq_%.2h_td_wc_%.2h_td_wr_%.2f_gqc_%.2f_gcr_%.2f",
            qs->nstates[0], qs->nstates[1], qs->nstates[2],
            sqrt(param.epsilon[0]*param.epsilon[0] + param.delta[0]*param.delta[0])/(2*M_PI), simulation.ho_w[0]/(2*M_PI), simulation.ho_w[1]/(2*M_PI),
            param.g_qc/(2*M_PI), param.g_cr/(2*M_PI));

    /*
     * Setup initial state
     */
    //rho1 = qdpack_dm_pure_TLS(0.8, 0.2);         // |> state
    rho_q = qdpack_dm_pure_TLS(0.9, 0.1);         // |> state
    rho_c = qdpack_dm_fock_state(0, qs->nstates[1]);     // |0> state
    rho_r = qdpack_dm_fock_state(0, qs->nstates[2]);     // |0> state
    
        qdpack_operator_list_init(&dm_list);
        qdpack_operator_list_append(&dm_list, rho_q);
        qdpack_operator_list_append(&dm_list, rho_c);
    qdpack_operator_list_append(&dm_list, rho_r);
    rho0 = qdpack_operator_tensor(qs, &dm_list);
    
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



