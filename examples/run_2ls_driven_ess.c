//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// The steady-state of a 2LS subject to a strong driving field: Find the steady
// state by evolving the density matrix for sufficiently long time.
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

static char DATADIR[256]; 
//#define DATADIR_TEMPLATE "/home/rob/qdpack-data/%s"
#define DATADIR_TEMPLATE "data/%s"

// excitation probability of the qubit, as calculated during the evolution.
static double pe_expt_sum = 0.0;
static int    pe_expt_no  = 0;

static qdpack_operator_t *rho_eb_t, *rho_tmp_t;

#define USE_EIGENBASIS 0

/**
 * Driving hamiltonian
 *
 */
int
hamiltonian_sd_t(qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t, qdpack_simulation_t *p)
{
    if (p->h_td_A == 0.0)
    {
        qdpack_operator_memcpy(H_t, H0);
    }
    else
    {
        qdpack_complex z;

        // H_t = H0 + f(t) H1
        // f(t) = A * sin(w * t)
 
        qdpack_operator_memcpy(H_t, H1);
        QDPACK_SET_COMPLEX(&z, p->h_td_A * sin(p->h_td_w * t) / 2, 0.0);
        qdpack_operator_scale(H_t, z);
        qdpack_operator_add(H_t, H0);
    }

    return 0;
}


/* ---
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.) 
 */
int
rho_store_final_cb(qdpack_operator_t *rho_t, double t, qdpack_hilbert_space_t *qs, qdpack_simulation_t *sp)
{
    qdpack_complex N_expt;
    char filename[1024], row[1024];
    
    snprintf(filename, sizeof(filename), "%s/p_e%s.dat", DATADIR, sp->simsig);
    
    //
    // store the average qubit excitation probability
    //
    snprintf(row, sizeof(row), "%f\t%f\t%f", sp->x, sp->y, pe_expt_sum / pe_expt_no);
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
    qdpack_complex N_expt, op1_expt, op2_expt;
    char filename[1024], row[16384];
    qdpack_operator_t *dm_part;
    
    if (USE_EIGENBASIS == 1)
    {
        gsl_eigen_hermv(sp->H_t, sp->eval, sp->evec, sp->w);
        gsl_eigen_hermv_sort(sp->eval, sp->evec, GSL_EIGEN_SORT_VAL_ASC); 

        gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, QDPACK_COMPLEX_ONE, sp->evec,  rho_t,    QDPACK_COMPLEX_ZERO, rho_tmp_t);
        gsl_blas_zgemm(CblasNoTrans,   CblasNoTrans, QDPACK_COMPLEX_ONE, rho_tmp_t, sp->evec, QDPACK_COMPLEX_ZERO, rho_eb_t);
    }
    else
    {
        qdpack_operator_memcpy(rho_eb_t, rho_t);
    }

    /* -- store expectation value of sigma z -- * / 
    N_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[i]);
    snprintf(filename, sizeof(filename), "%s/sigma_z_expt_%d%s_x_%.2f_y_%.2f.dat", DATADIR, i, sp->simsig, sp->x, sp->y);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
    qdpack_file_row_add(filename, row);        
    */

    //if (1) // if ( abs(sin(p->h_td_w * t)) < 0.1 )  // close to bias point
    if (t > sp->Tf/2.0)
    {
        pe_expt_sum += QDPACK_REAL(qdpack_operator_get(rho_eb_t, 1, 1));
        pe_expt_no++;
    }
    
    
    return 0;
}


/* ---
 * Simulation parameters.
 * 
 */
double epsilon[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ho_w[]    = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
    int qsn, i, noff, nsize, npeak;
    double Navg, fw, fA, phi, Gc, Gq, gqc, x, y, pDelta, G1, G2, T1, T2;
    double yi, yf, yd;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    if (argc != 5)
    {
        printf("usage: %s x y-init y-delta y-final\n", argv[0]);
        exit(0);
    }

    snprintf(DATADIR, sizeof(DATADIR), DATADIR_TEMPLATE, argv[0]);

    // ---- PARAMETERS START

    //
    // frequency in units of GHz and time in units of ns
    //

    //
    // x = coefficient of sigma_z in qubit hamiltonian, in units of driving frequency
    // x = epsilon / (fw * 2 * M_PI)
    //
    // y = the driving amplitude in units of driving frequency
    // y = h_td_A / (fw * 2 * M_PI)
    //
    
    // get values from command line parameters
    x      = strtod(argv[1], NULL);
    yi     = strtod(argv[2], NULL);
    yd     = strtod(argv[3], NULL);
    yf     = strtod(argv[4], NULL);
    
    fw     = 1.0;               // driving frequency
    pDelta  = fw/300.0;         // coefficient of sigma_x in qubit hamiltonain
    T1     = 1000.0 / fw;       // relaxation time
    T2     =    5.0 / fw;       // dephasing  time

    param.Ti = 0.0;
    param.Tf = T1 * 5.0;       // must be large enough to reach the steady state
    param.dT = (param.Tf-param.Ti) / 1000.0;

    // ---- PARAMETERS END

    //
    // setup parameter vectors
    //
    param.epsilon = epsilon;
    param.delta   = delta;
    simulation.ho_w    = ho_w;
    param.x = x;
    for (i = 0; i < 10; i++)
    {
        param.epsilon[i] *= 2*M_PI;
        param.delta[i]   *= 2*M_PI;
        simulation.ho_w[i]    *= 2*M_PI;
        lambda[i]         = (double *)calloc(10, sizeof(double));
    }
    param.lambda = lambda;
    param.cont_basis_transform = 0;
    param.h_td_w = fw * 2 * M_PI;

    G1 = 1.0/T1;
    G2 = 1.0/T2;

    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;

    // setup hamiltonian functions
    param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    param.H1_func = (hamiltonian_func_t)hamiltonian_z_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_sd_t;

    // preallocate matrices
    param.H_t = qdpack_operator_alloc(qsn, qsn);
    param.H0 = qdpack_operator_alloc(qsn, qsn);
    param.H1 = qdpack_operator_alloc(qsn, qsn);
      param.N_op[0] = qdpack_operator_alloc(qsn, qsn);
    operator_ho_N(param.N_op[0], qs, 0, 0);

    /*
     * Set up dissipation
     *
     */
    param.wT = 0.0 * 2 * M_PI;
    param.do_n = 0;
    Navg = 0;  // zero temperature

    // relaxation
    param.do_a[param.do_n]  = qdpack_operator_alloc(qsn, qsn);
    operator_ho_lowering(param.do_a[param.do_n], qs, 0, 0); 
    param.do_ad[param.do_n]    = qdpack_operator_alloc(qsn, qsn); 
    operator_ho_raising(param.do_ad[param.do_n], qs, 0, 0);
    param.do_g1[param.do_n] = G1 * (1 + Navg); // relaxation rate
    param.do_g2[param.do_n] = G1 * Navg;       // excitation rate
    param.do_n++;

    // dephasing
    param.do_a[param.do_n]  = qdpack_operator_alloc(qsn, qsn);
    operator_sigma_z(param.do_a[param.do_n], qs, 0);
    param.do_ad[param.do_n] = qdpack_operator_alloc(qsn, qsn);
    operator_sigma_z(param.do_ad[param.do_n], qs, 0); 
    param.do_g1[param.do_n] = G2 / 2; 
    param.do_g2[param.do_n] = G2 / 2; 
    param.do_n++;

    param.eval = gsl_vector_alloc(2);
    param.evec = qdpack_operator_alloc(2, 2);
    param.w    = gsl_eigen_hermv_alloc(2);

    snprintf(param.simsig, sizeof(param.simsig), "_dDelta_%.3f_G1_%.3f_G2_%.4f_fw_%.3f", pDelta, G1, G2, fw);

    /* storage space for final rho */
    param.rho_f    = qdpack_operator_alloc(qsn, qsn);
    rho_eb_t  = qdpack_operator_alloc(qsn, qsn);
    rho_tmp_t = qdpack_operator_alloc(qsn, qsn);

    /*
     * Setup initial state
     */
    rho_q = qdpack_dm_pure_TLS(0.0);
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho_q);
    rho0 = qdpack_operator_alloc(qsn, qsn);
    qdpack_operator_tensor(rho0, qs, &dm_list);

    delta[0]   = pDelta * (2*M_PI);
    epsilon[0] = x * fw * (2*M_PI);

    for (y = yi; y <= yf; y += yd)
    {    
        pe_expt_sum = 0.0;
        pe_expt_no  = 0;

        param.y = y;
        param.h_td_A  = y * fw * 2 * M_PI;

        printf("ITERATION: %s: x = %f, y = %f\n", param.simsig, param.x, param.y);

        if (USE_EIGENBASIS == 1)
        {
            param.H0_func(param.H0, qs, &param);

            gsl_eigen_hermv(param.H0, param.eval, param.evec, param.w);
            gsl_eigen_hermv_sort(param.eval, param.evec, GSL_EIGEN_SORT_VAL_ASC);

            // rho_0 = S rho_eb_0 S^{-1}
               gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   QDPACK_COMPLEX_ONE, param.evec, rho_q,      QDPACK_COMPLEX_ZERO, rho_tmp_t);
               gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, QDPACK_COMPLEX_ONE, rho_tmp_t,  param.evec, QDPACK_COMPLEX_ZERO, rho0);

            qdpack_operator_free(param.H0);
        }
            
        /*
         * Evolve the quantum system until time t = T, where the steady state
         * is assumed to be reached.
         * 
         */
        //if (qdpack_evolve_dm_lme_ib_t(qs, rho0, &param, rho_store_cb) != 0)
        if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
        {
            fprintf(stderr, "Evolution of the quantum system failed.\n");
            return -1;
        }
    
        rho_store_final_cb(param.rho_f, param.Tf, qs, &param);

    }

    return 0;

}



