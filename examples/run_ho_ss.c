//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Calculate the steady state of the qubit + cavity.
// 

#include <math.h>
#include <stdio.h> 
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include <qdpack/qdpack.h>

#define DATADIR "/home/rob/qdpack-data/ho_ss"

/* ---
 * Evaluate the correlation function <A(t+tau)A(t)> at tau
 *
 */
int
rho_corr_cb(qdpack_operator_t *op_tau, double tau, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
        qdpack_complex corr;
        char filename[1024], row[1024];

        corr = qdpack_operator_expectation_value(sp->A, op_tau);

        snprintf(filename, sizeof(filename), "%s/C%s.dat", DATADIR, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", tau, QDPACK_REAL(corr), QDPACK_IMAG(corr));
        qdpack_file_row_add(filename, row);

        return 0;
}

/* ---
 * Store the final (steady-state) density matrix to file.
 */
int
rho_store_final_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    char filename[128];

        printf("rho_store_final_cb: executing final-time rho calculations: \n");
    
    snprintf(filename, sizeof(filename), "%s/rho_ss%s.dat", DATADIR, sp->simsig);

    qdpack_operator_write(filename, rho_t, rho_t->size1);
    

        {
                printf("rho_store_final_cb: calculating Q function\n");

                qdpack_operator_t *Q;
                int Npnts = 101;
                double alpha_max = 5.0;

                Q = qdpack_operator_alloc(Npnts, Npnts);

                distribution_function_Q(Q, rho_t, alpha_max);

//                printf("Q_real =\n");
//                qdpack_operator_print_real(Q, Npnts);
//                printf("Q_imag =\n");
//                qdpack_operator_print_imag(Q, Npnts);

                snprintf(filename, sizeof(filename), "%s/Q_real%s.dat", DATADIR, sp->simsig);
                qdpack_operator_write_double(filename, Q, Npnts, 1);
                snprintf(filename, sizeof(filename), "%s/Q_imag%s.dat", DATADIR, sp->simsig);
                qdpack_operator_write_double(filename, Q, Npnts, 0);
        }

        {
                printf("rho_store_final_cb: calculating wigner function\n");

                qdpack_operator_t *W;
                int Npnts = 201;
                qdpack_complex alpha_max;
                QDPACK_SET_COMPLEX(&alpha_max, 5.0, 5.0);

                W = qdpack_operator_alloc(Npnts, Npnts);

                //distribution_function_wigner(W, sp->a, sp->ad, rho_t, alpha_max);
//                distribution_function_characteristic_w(W, sp->a, sp->ad, rho_t, alpha_max);
////                distribution_function_wigner_cf_quad(W, sp->a, sp->ad, rho_t, alpha_max);
                distribution_function_wigner_cf_sum(W, sp->a, sp->ad, rho_t, alpha_max);

//                printf("Q_real =\n");
//                qdpack_operator_print_real(Q, Npnts);
//                printf("Q_imag =\n");
//                qdpack_operator_print_imag(Q, Npnts);

                snprintf(filename, sizeof(filename), "%s/W_real%s.dat", DATADIR, sp->simsig);
                qdpack_operator_write_double(filename, W, Npnts, 1);
                snprintf(filename, sizeof(filename), "%s/W_imag%s.dat", DATADIR, sp->simsig);
                qdpack_operator_write_double(filename, W, Npnts, 0);
        }

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
    qdpack_complex op_expt;
    char filename[1024], row[16384];
    qdpack_operator_t *dm_part;

    /*
     * Store the density matrix at time t:
      */ 
    //ri = 0;    
    //snprintf(filename, sizeof(filename), "%s/rho_diag%s.dat", DATADIR, sp->simsig);
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


    /*
     * For each sub-system:
     *
     */
    for (i = 0; i < qs->nsubsys; i++)
        {
        ri = 0;    
        
        snprintf(filename, sizeof(filename), "%s/rho_ho_%d_diag%s.dat", DATADIR, i, sp->simsig);
        
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        ri = snprintf(row, sizeof(row), "%f", t);

        for (j = 0; j < qs->nstates[i]; j++)
        {
             ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", (double)QDPACK_REAL(qdpack_operator_get(dm_part, j, j)));
        }
        qdpack_file_row_add(filename, row);
        
        /* -- expectation values -- */ 
        op_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[i]);
        snprintf(filename, sizeof(filename), "%s/N_expt_%d%s.dat", DATADIR, i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt), sp->do_g1[0]);
        qdpack_file_row_add(filename, row);        
        
        op_expt = qdpack_operator_expectation_value(rho_t, sp->a);
        snprintf(filename, sizeof(filename), "%s/a_expt_%d%s.dat", DATADIR, i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->ad);
        snprintf(filename, sizeof(filename), "%s/ad_expt_%d%s.dat", DATADIR, i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        
        //qdpack_operator_free(dm_part);
        }
        
    /* -- expectation values -- */
    op_expt = qdpack_operator_expectation_value(rho_t, sp->A);
    snprintf(filename, sizeof(filename), "%s/A_expt%s.dat", DATADIR, sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
    qdpack_file_row_add(filename, row);        
        
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
double ho_w[]    = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double *lambda[10];


/* ---
 * Program starts here.
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    gsl_vector_complex *wfs, *wf1, *wf2;
    qdpack_operator_t *rho_q, *rho_c;
    qdpack_operator_t *rho0, *rho_t, *rho1;
    qdpack_operator_list_t dm_list;
    int qsn, i, noff, nsize;
    double Navg, fw, fA, phi, Gc, Gq, gqc, wT;
    qdpack_complex z;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);


    if (argc != 6)
    {
        printf("usage: %s check code\n", argv[0]);
        exit(0);
    }
    
    nsize   = strtol(argv[1], NULL, 10);
    fw      = strtod(argv[2], NULL);
    fA      = strtod(argv[3], NULL);
    Gc    = strtod(argv[4], NULL);
    wT    = strtod(argv[5], NULL);

    /*
     * Setup solver parameters.
     *
     */
    simulation.ho_w    = ho_w;
    for (i = 0; i < 10; i++)
    {
        simulation.ho_w[i]    *= 2*M_PI;
    }
    
    param.Ti = 0.0;
    param.dT = 0.1;        // small compared to T = 1/(simulation.ho_w[0]/2*pi) = 1
    param.Tf = 1000.0;
    param.n_offset = 0; 
    param.N = nsize; 

    printf("Using cavity with %d states, GND offset %d, and decay rate Gc = %f and temperature wT = %f.\n", param.N, param.n_offset, Gc, wT);
    
    param.h_td_A = fA * 2 * M_PI;
    param.h_td_w = fw * ho_w[0];

    param.H0_func = (hamiltonian_func_t)hamiltonian_ho;
    param.H1_func = (hamiltonian_func_t)hamiltonian_ho_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_t;
    
    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, param.N);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;

    /*
     * Set up dissipation
     *
     */
    param.wT = wT * 2 * M_PI;
    param.do_n = 0;

    /* --- cavity --- */
    param.do_a[param.do_n]  = qdpack_operator_alloc(qsn, qsn);
    operator_ho_lowering(param.do_a[param.do_n], qs, 0, param.n_offset); 
    param.do_ad[param.do_n]    = qdpack_operator_alloc(qsn, qsn); 
    operator_ho_raising(param.do_ad[param.do_n], qs, 0, param.n_offset);

    if (param.wT > 0.0)
    {
        Navg = 1.0/(exp(simulation.ho_w[0]/param.wT) - 1.0);
        printf("Setting Navg = %f (param.wT = %f)\n", Navg, param.wT);
    }
    else
    {
        printf("setting Navg = 0 without calc. (param.wT = %f)\n", param.wT);
        Navg = 0;
    }
    param.do_g1[param.do_n] = Gc * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = Gc * Navg;       // excitation
    printf("Debug (Navg = %f): do_g1 = %f, go_g2 = %f\n", Navg, param.do_g1[param.do_n], param.do_g2[param.do_n]);
    param.do_n++;
    
    /* 
     * Pre-calc of operators to calculate expectation values for.
     */
    param.N_op[0] = qdpack_operator_alloc(qsn, qsn); 
    operator_ho_N(param.N_op[0], qs, 0, param.n_offset);

    param.a     = qdpack_operator_alloc(qsn, qsn); 
    operator_ho_lowering(param.a,  qs, 0, param.n_offset);
    param.ad    = qdpack_operator_alloc(qsn, qsn); 
    operator_ho_raising (param.ad, qs, 0, param.n_offset);

        param.A  = qdpack_operator_alloc(qsn, qsn);
        qdpack_operator_memcpy(param.A, param.a);
        qdpack_operator_add(param.A, param.ad);
        
    snprintf(param.simsig, sizeof(param.simsig), "_Nc_%d_Noffset_%d_wc_%.3f_Gc_%.4f_fA_%.4f_fw_%.3h_td_wT_%.3f",
        qs->nstates[0], param.n_offset, simulation.ho_w[0]/(2*M_PI), Gc, param.h_td_A / (2 * M_PI), param.h_td_w / (2 * M_PI), param.wT / (2 * M_PI));

        param.H_t = qdpack_operator_alloc(qsn, qsn);
        param.H0 = qdpack_operator_alloc(qsn, qsn);
        param.H1 = qdpack_operator_alloc(qsn, qsn);

    /*
     * Setup initial state
     */
    rho0 = qdpack_operator_alloc(qsn, qsn);
    rho1 = qdpack_operator_alloc(qsn, qsn);
    rho_c = qdpack_operator_alloc(qsn, qsn);
    wf1  = gsl_vector_complex_alloc(qsn);
    wf2  = gsl_vector_complex_alloc(qsn);
    wfs  = gsl_vector_complex_alloc(qsn);

//      rho_c = qdpack_dm_fock_state(1, qs->nstates[0]);

//    wf1 = qdpack_wf_fock_state(0, qs->nstates[0]); 
//    wf2 = qdpack_wf_fock_state(2, qs->nstates[0]);     
//    qdpack_wf_superposition(wfs, wf1, wf2);
//    qdpack_state_to_operator(rho_c, wfs);       

//    qdpack_wf_coherent_state(wf1, sqrt(4), 0 * M_PI, qs->nstates[0], 0);
//    qdpack_wf_coherent_state(wf2, sqrt(4), 1 * M_PI, qs->nstates[0], 0);
//    qdpack_wf_superposition(wfs, wf1, wf2);
//    qdpack_state_to_operator(rho_c, wfs);       
  
// core dump    
  rho0 = qdpack_dm_fock_state(3, qs->nstates[0]);     // |0> state
  rho1 = qdpack_dm_fock_state(12, qs->nstates[0]);     // |0> state
  qdpack_dm_mix(rho_c, rho0, rho1);

//  qdpack_wf_coherent_state(wfs, sqrt(4), M_PI/2.0, qs->nstates[0], 0);
//  qdpack_state_to_operator(rho_c, wfs);       

//    rho_c = qdpack_dm_coherent_state(sqrt(4), M_PI/2.0, qs->nstates[0], 0);

//  rho0 = qdpack_dm_coherent_state(sqrt(4), 0 * M_PI, qs->nstates[0], 0);
//  rho1 = qdpack_dm_coherent_state(sqrt(4), 1 * M_PI, qs->nstates[0], 0);
//  qdpack_dm_mix(rho_c, rho0, rho1);
    
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho_c);
    qdpack_operator_tensor(rho0, qs, &dm_list);

    z = qdpack_operator_trace(rho0);
    printf("initial dm trace check: real = %f, imag = %f\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    /* storage space for final rho */
    param.rho_f = qdpack_operator_alloc(qsn, qsn);
    
    /*
     * Evolve the quantum system until time t = T, where the steady state
         * is assumed to be reached.
     * 
     */
    //printf("Start evolving system\n");
    //if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_he_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const_full_steps(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const_simple(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    if (0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

//    printf("rho_dyn =\n");
//        qdpack_operator_print_real(param.rho_f, qsn);

//    qdpack_operator_set_zero(param.rho_f);
//    qdpack_steadystate_dm(qs, param.rho_f, &param);        
//    printf("rho_ss =\n");
//      qdpack_operator_print_real(param.rho_f, qsn);

//    qdpack_operator_set_zero(param.rho_f);
//    qdpack_steadystate_dm_sparse(qs, param.rho_f, &param);    
//    printf("rho_ss_sparse =\n");
//      qdpack_operator_print_real(param.rho_f, qsn);

    /* 
         * At this point the system (rho_f) is supposed to be in a steady-state,
     * such that C(t,tau) = C(tau).
     *
     */
//    rho_store_cb(param.rho_f, param.Tf, qs, &param);

    rho_store_final_cb(rho0, param.Tf, qs, &param);
//    rho_store_final_cb(param.rho_f, param.Tf, qs, &param);

    return 0;


    /*
      * Read steady-state density matrix form file:
      */
    //param.rho_f = qdpack_operator_alloc(qsn, qsn);
    //misc
    
    printf("Start calculating the correlation function\n");
        // here, assume that rho0 = rho(t), proceed to calculate Corr(t,tau)

        //param.Ti = 0.00;
        //param.Tf = 2000.00;       // tau
        //param.dT = 0.05;

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



