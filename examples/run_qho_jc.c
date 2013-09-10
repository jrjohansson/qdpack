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
 * Evaluate the correlation function <A(t+tau)A(t)> at tau
 *
 */
int
rho_corr_cb(qdpack_operator_t *op_tau, double tau, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
        int qsn;
        qdpack_complex corr1, corr2;
        qdpack_operator_t *c;
        char filename[1024], row[1024];

        qsn = qdpack_hilbert_space_nstates(qs);
        c = qdpack_operator_alloc(qsn, qsn);

        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, QDPACK_COMPLEX_ONE, sp->A,  op_tau,  QDPACK_COMPLEX_ZERO, c);

        corr1 = operator_trace(c);

        qdpack_operator_free(c);

        corr2 = qdpack_operator_expectation_value(sp->A, op_tau);

        snprintf(filename, sizeof(filename), "data.qho_es/C%s.dat", sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f", tau, QDPACK_REAL(corr1), QDPACK_IMAG(corr1), QDPACK_REAL(corr2), QDPACK_IMAG(corr2));
        qdpack_file_row_add(filename, row);

        return 0;
}

/* ---
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.) 
 */
int
rho_store_final_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    qdpack_complex N_expt; //, op_expt;
    char filename[128], row[1024];
    
    snprintf(filename, sizeof(filename), "data.qho_es/rho_f_expt_%d.dat", sp->N);
    
    /* -- expectation values -- */
    N_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[1]);
    snprintf(row, sizeof(row), "%f\t%f\t%f\t%f", sp->h_td_w/(2*M_PI), sp->h_td_A/(2*M_PI), QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
    qdpack_file_row_add(filename, row);        

    return 0;
}
    
/* ---
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.) 
 */
int
rho_store_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    int i, j, ri;
    qdpack_complex N_expt; //, op_expt;
    char filename[1024], row[16384];
    qdpack_operator_t *dm_part; //, *op;

    //ri = 0;    
    //snprintf(filename, sizeof(filename), "data.qho_es/rho_diag%s.dat", sp->simsig);
    //ri = snprintf(row, sizeof(row), "%f", t);
    //for (j = 0; j < rho_t->size1; j++)
    //{
    //     ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", (double)QDPACK_REAL(qdpack_operator_get(rho_t, j, j)));
    //    if (ri >= sizeof(row))
    //    {
    //        fprintf(stderr, "Error: Preparing row for %s overflowed.\n", sp->simsig);
    //        break;
    //    }
    //}
    //qdpack_file_row_add(filename, row);
    
    for (i = 0; i < qs->nsubsys; i++)
    {
        ri = 0;    

        snprintf(filename, sizeof(filename), "data.qho_es/rho_ho_%d_diag%s.dat", i, sp->simsig);
        
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        ri = snprintf(row, sizeof(row), "%f", t);

        for (j = 0; j < qs->nstates[i]; j++)
        {
             ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", (double)QDPACK_REAL(qdpack_operator_get(dm_part, j, j)));
        }
        qdpack_file_row_add(filename, row);
        
        /* -- expectation values -- */ 
        N_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[i]);
        snprintf(filename, sizeof(filename), "data.qho_es/N_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt), sp->do_g1[0]);
        qdpack_file_row_add(filename, row);        
        
        /* -- trace check values -- * /
        op = gsl_ext_matrix_convert_and_free(operator_unit(qs, i));
        N_expt = qdpack_operator_expectation_value(dm_part, op);
        snprintf(filename, sizeof(filename), "data.qho_es/Trace_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
        qdpack_file_row_add(filename, row);        
        qdpack_operator_free(op);
        */    
        
/*
        op_expt = qdpack_operator_expectation_value(rho_t, sp->a);
        snprintf(filename, sizeof(filename), "data.qho_corr/a_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->ad);
        snprintf(filename, sizeof(filename), "data.qho_corr/ad_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b);
        snprintf(filename, sizeof(filename), "data.qho_corr/b_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd);
        snprintf(filename, sizeof(filename), "data.qho_corr/bd_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b2);
        snprintf(filename, sizeof(filename), "data.qho_corr/b2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd2);
        snprintf(filename, sizeof(filename), "data.qho_corr/bd2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        
*/
        
        //qdpack_operator_free(dm_part);
        }
        
    /* -- expectation values -- * /
    N_expt = qdpack_operator_expectation_value(rho_t, sp->A);
    snprintf(filename, sizeof(filename), "data.qho_es/E_expt%s.dat", sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
    qdpack_file_row_add(filename, row);        
    */
        
    /* -- trace check values -- * /
    op = gsl_ext_matrix_convert_and_free(operator_unit(qs, -1));
    N_expt = qdpack_operator_expectation_value(rho_t, op);
    snprintf(filename, sizeof(filename), "data.qho_es/Trace_full%s.dat", sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
    qdpack_file_row_add(filename, row);        
    qdpack_operator_free(op);
    */
    
    if (t >= sp->Tf-sp->dT)
    {
        printf("Storing final rho: rho_f\n");
         qdpack_operator_memcpy(sp->rho_f, rho_t);
    }
    
    return 0;
}


/* ---
 * Simulation parameters.
 * 
 */
double epsilon[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ho_w[]    = {0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double *lambda[10];


/* ---
 * Program starts here.
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho_q, *rho_c;
    qdpack_operator_t *rho0, *rho_t;
    qdpack_operator_list_t dm_list;
    int qsn, i, noff, nsize;
    double Navg, fw, fA, phi, Gc, Gq, gqc;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    if (argc != 6)
    {
        printf("usage: %s check code\n", argv[0]);
        exit(0);
    }
    
    noff    = strtol(argv[1], NULL, 10);
    nsize   = strtol(argv[2], NULL, 10);
    //fw     = strtod(argv[2], NULL);
    fA      = strtod(argv[3], NULL);
    gqc    = strtod(argv[4], NULL);
    Gq    = strtod(argv[5], NULL);

    printf("Input arguments: phi = %f\n", phi);
    
    /*
     * Setup solver parameters.
     *
     */
    
    // the qubit
    epsilon[0] = 0.1;
    delta[0] = 0.0;// 98.99; //13.7 * cos(fabs(M_PI * phi));

    Gc  = 0.005;    // cavity decay rate
    //Gq  = 1.0;    // cavity decay rate
    gqc = 0.008;

    
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

    param.g_qc = gqc * 2*M_PI; // 0.01
    param.g_cr = 0.0 * 2*M_PI; // 0.01
    
    printf("Qubit-cavity coupling: g_qc = %f\n", param.g_qc);
    
    param.Ti = 0.0;
    param.Tf = 7500.0;
    param.dT = 1.0;

    param.N = nsize; 
    param.n_offset = noff;
    param.n_offset = (int)(200 * (Gq/2.0) * (1 - Gq/2.0) - nsize/2.0);
    if (param.n_offset < 0)
        param.n_offset = 0;
    
    //param.N = ((int)(Gq/(2*Gc)) + 15);
    //if (param.N > 65)
    //{
    //    param.N = 65;
    //}

    printf("Using cavity with %d states, and GND offset %d for Gq = %f.\n", param.N, param.n_offset, Gq);
    
    //fA = 0.25;
    //fw = 1.0;
    param.h_td_A = fA * 2 * M_PI;
    param.h_td_w = fw * ho_w[0];

    //param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    param.H0_func = (hamiltonian_func_t)hamiltonian_jc;
    param.H1_func = (hamiltonian_func_t)hamiltonian_jc_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_t;
    
    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_add(qs, param.N);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;

//        param.H0 = param.H0_func(param.qs, &param);
//        printf("H_real =\n");
//        qdpack_operator_print_real(param.H0, 2*param.N);
//        printf("H_imag =\n");
//        qdpack_operator_print_imag(param.H0, 2*param.N);
//        exit(0);

    /*
     * Set up dissipation
     *
     */
    param.wT = 0.0 * 2 * M_PI;
    param.do_n = 0;
    Navg = 0;

    /* --- qubit --- */
    param.do_a[param.do_n]  = qdpack_operator_alloc(qsn, qsn);
    operator_ho_lowering(param.do_a[param.do_n], qs, 0, 0); 
    param.do_ad[param.do_n] = qdpack_operator_alloc(qsn, qsn);
    operator_ho_raising(param.do_ad[param.do_n], qs, 0, 0); 
    Navg = 1/(exp(sqrt(param.epsilon[0]*param.epsilon[0] + param.delta[0]*param.delta[0])/param.wT) - 1);
    printf("qubit Navg = %f\n", Navg);
    param.do_g1[param.do_n] = Gq * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = Gq * Navg;       // excitation
    param.do_n++;
        
    /* --- cavity --- */
    param.do_a[param.do_n]  = qdpack_operator_alloc(qsn, qsn);
    operator_ho_lowering(param.do_a[param.do_n], qs, 1, param.n_offset); 
    param.do_ad[param.do_n] = qdpack_operator_alloc(qsn, qsn);
    operator_ho_raising(param.do_ad[param.do_n], qs, 1, param.n_offset); 
    Navg = 1/(exp(simulation.ho_w[0]/param.wT) - 1);
    param.do_g1[param.do_n] = Gc * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = Gc * Navg;       // excitation
    param.do_n++;
    
    /* --- transmission line --- * /
     param.do_a[param.do_n] = operator_ho_lowering(qs, 2, param.n_offset); 
    Navg = 1/(exp(simulation.ho_w[1]/param.wT) - 1);
    param.do_g1[param.do_n] = 0.05 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 0.05 * Navg;       // excitation
    param.do_n++;
    */

    /* 
     * Pre-calc of operators to calculate expectation values for.
     */
    param.N_op[0] = qdpack_operator_alloc(qsn, qsn);
    operator_ho_N(param.N_op[0], qs, 0, 0);
    param.N_op[1] = qdpack_operator_alloc(qsn, qsn);
    operator_ho_N(param.N_op[1], qs, 1, param.n_offset);

    param.a  = qdpack_operator_alloc(qsn, qsn);
    operator_ho_lowering(param.a, qs, 1, param.n_offset);
    param.ad = qdpack_operator_alloc(qsn, qsn);
    operator_ho_raising(param.ad, qs, 1, param.n_offset);

    param.A  = qdpack_operator_alloc(qsn, qsn);
    qdpack_operator_memcpy(param.A, param.a);
    qdpack_operator_add(param.A, param.ad);
        
    snprintf(param.simsig, sizeof(param.simsig),
            "_Nq_%d_Nc_%d_Noffset_%d_wq_%.3h_td_wc_%.3h_td_wr_%.3f_gqc_%.3f_Gq_%.3f_Gc_%.4f_qpump",
            qs->nstates[0], qs->nstates[1], param.n_offset,
            sqrt(param.epsilon[0]*param.epsilon[0] + param.delta[0]*param.delta[0])/(2*M_PI),
            simulation.ho_w[0]/(2*M_PI), simulation.ho_w[1]/(2*M_PI),
            param.g_qc/(2*M_PI), Gq, Gc);

    /*
     * Setup initial state
     */
    rho_q = qdpack_dm_pure_TLS(1.0);                     // |1> state
    rho_c = qdpack_dm_fock_state(0, qs->nstates[1]);     // |0> state
    //rho_c = qdpack_dm_coherent_state(sqrt(param.N), 0, qs->nstates[1]);     // |sqrt(param.N)> state
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho_q);
    qdpack_operator_list_append(&dm_list, rho_c);
    rho0 = qdpack_operator_alloc(qsn, qsn);
    qdpack_operator_tensor(rho0, qs, &dm_list);
        
    /* storage space for final rho */
    param.rho_f = qdpack_operator_alloc(qsn, qsn);
    
    /*
     * Evolve the quantum system until time t = T, where the steady state
     * is assumed to be reached.
     * 
     */
    printf("Start evolving system\n");
    if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_he_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const_full_steps(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const_simple(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    /* 
    * At this point the system (rho_f) is supposed to be in a steady-state,
     * such that C(t,tau) = C(tau).
     *
     */
    rho_store_final_cb(param.rho_f, param.Tf, qs, &param);
    
    printf("Start calculating the correlation function\n");
    // here, assume that rho0 = rho(t), proceed to calculate Corr(t,tau)

    param.Ti =  0.00;
    param.Tf =  5000.00;       // tau
    param.dT =  0.2;

    rho_t = qdpack_operator_alloc(qsn, qsn);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, QDPACK_COMPLEX_ONE, param.A,  param.rho_f,  QDPACK_COMPLEX_ZERO, rho_t);
    //gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, QDPACK_COMPLEX_ONE, param.A,  rho0,  QDPACK_COMPLEX_ZERO, rho_t);

    /*
     * Evolve the B rho(t)
     *
     */
    if (qdpack_evolve_dm_lme_t(qs, rho_t, &param, rho_corr_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho_t, &param, rho_corr_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the B rho failed.\n");
        return -1;
    }
    
    return 0;
}



