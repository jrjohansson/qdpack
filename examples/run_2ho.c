//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Simulate the time evolution of two coupled oscillators. 
 *
 * Usage:
 *
 * 1) Run this program to generate data files (it outputs the entropy)
 *
 * data.2ho_cool/wT_xxxx.dat
 *
 * Format:
 * 
 *   t        wTeff[0]    wTeff[1]    wS[0]    wS[1]
 *
 *   wTeff = calculated from boltzman dist
 *   wS    = entropy
 *
 * Parameters:
 *
 *  time ./run_2ho 0.1 0.05 0.05 0.90 0.005
 *
 *  N = 8, 40        very good
 *  N = 5, 15        crap
 *
 * 2) Run MATLAB program that converts entropy to temperature.
 * 
 * 3) Plot the figure using gnuplot:
 */

#include <math.h>
#include <stdio.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#include <qdpack/qdpack.h>

/*
 * Generate the hamiltonian at time t.
 *
 */
int
hamiltonian_t_2ho(qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t, qdpack_simulation_t *p)
{
          static int field_turned_on = 0;
        qdpack_complex z;

        //qdpack_operator_memcpy(H_t, H0);
        //qdpack_operator_memcpy(H_t, H1);
        //QDPACK_SET_COMPLEX(&z, p->h_td_A * sin(p->h_td_w * t), 0.0);
        //qdpack_operator_scale(H_t, z);
        //qdpack_operator_add(H_t, H0);
        //return 0;

        if (t < 75.0)
        {
                qdpack_operator_memcpy(H_t, H0);
        }
//      else if (t < 75.0)
//      {
//              // turn off dissipation on the small oscillator
//              qdpack_operator_memcpy(H_t, H0);
//              p->do_n = 1;
//              p->do_g1[1] = 0;
//              p->do_g2[1] = 0;
//      }
        else
        {
                if (field_turned_on == 0)
                {
                        // restore original decay rates for the large oscillator
                        p->do_g1[0] = p->do_g1[0]/10;
                        p->do_g2[0] = p->do_g2[0]/10;
                        // turn off dissipation on the low-frequency oscillator
                        p->do_n = 1;
                        p->do_g1[1] = 0;
                        p->do_g2[1] = 0;
                        // Turn on coupling
                        p->ho_g = p->ho_g2; //0.01 * 2 * M_PI;
                        qdpack_operator_free(p->H0);
                        p->H0 = p->H0_func(p->qs, p);
                        field_turned_on++;
                }
                qdpack_operator_memcpy(H_t, H1);
                QDPACK_SET_COMPLEX(&z, p->h_td_A * sin(p->h_td_w * t), 0.0);
                qdpack_operator_scale(H_t, z);
                qdpack_operator_add(H_t, p->H0);
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
    int i, j;
    qdpack_complex N_expt; //, op_expt;
    char filename[128], row[1024];
    qdpack_operator_t *dm_part;
    double wTeff[2], wS[2];
    
        for (i = 0; i < qs->nsubsys; i++)
        {
        int ri = 0;    
        
        snprintf(filename, sizeof(filename), "data.2ho_cool/rho_ho_%d_diag%s.dat", i, sp->simsig);
                
        dm_part = qdpack_operator_traceout(qs, rho_t, i);

        wS[i] = qdpack_entanglement_neumann_entropy(dm_part);
        
        ri = snprintf(row, sizeof(row), "%f", t);

        for (j = 0; j < qs->nstates[i]; j++)
             ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", QDPACK_REAL(qdpack_operator_get(dm_part, j, j)));
    
        qdpack_file_row_add(filename, row);


        /* -- expectation values -- */
        N_expt = qdpack_operator_expectation_value(rho_t, sp->N_op[i]);
        snprintf(filename, sizeof(filename), "data.2ho_cool/N_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(N_expt), QDPACK_IMAG(N_expt));
        qdpack_file_row_add(filename, row);        
        
        sp->Nexpt[i] = QDPACK_REAL(N_expt);
/*
        op_expt = qdpack_operator_expectation_value(rho_t, sp->a);
        snprintf(filename, sizeof(filename), "data/a_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->ad);
        snprintf(filename, sizeof(filename), "data/ad_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b);
        snprintf(filename, sizeof(filename), "data/b_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd);
        snprintf(filename, sizeof(filename), "data/bd_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->b2);
        snprintf(filename, sizeof(filename), "data/b2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        

        op_expt = qdpack_operator_expectation_value(rho_t, sp->bd2);
        snprintf(filename, sizeof(filename), "data/bd2_expt_%d%s.dat", i, sp->simsig);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(op_expt), QDPACK_IMAG(op_expt));
        qdpack_file_row_add(filename, row);        
*/
        
        qdpack_operator_free(dm_part);
        }
    
    /* --- store effective temperatures calculated above --- */
    wTeff[0] = sp->ho_w[0] / log(1/sp->Nexpt[0] + 1);
    wTeff[1] = sp->ho_w[1] / log(1/sp->Nexpt[1] + 1);
    snprintf(filename, sizeof(filename), "data.2ho_cool/wT_eff%s.dat", sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f", t, wTeff[0], wTeff[1], wS[0], wS[1]);
    qdpack_file_row_add(filename, row);        

        
    return 0;
}


/* ---
 * Simulation parameters.
 * 
 */
double epsilon[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ho_w[]    = { 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
    double Navg, wTeff[2], g0, fw, fA, gho;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    if (argc != 6)
    {
        printf("usage: %s w1 fA g0 gh0\n", argv[0]);
    }
    
    ho_w[1] = strtod(argv[1], NULL);
    fA      = strtod(argv[2], NULL);
    g0      = strtod(argv[3], NULL);
    fw      = strtod(argv[4], NULL);
    gho     = strtod(argv[5], NULL);

    printf("Input arguments: ho_w[1] = %f, fA = %f, g0 = %f\n", ho_w[1], fA, g0);
    
    /*
     * Setup solver parameters.
     *
     */
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
    lambda[0][1] = 0.25 * 2*M_PI;
    lambda[1][2] = 0.40 * 2*M_PI;
    lambda[2][3] = 0.30 * 2*M_PI;
    lambda[3][4] = 0.20 * 2*M_PI;
    lambda[4][5] = 0.10 * 2*M_PI;
    
    param.lambda  = lambda;

    simulation.ho_g  = 0.0 * 2*M_PI; // 0.01
    simulation.ho_g2 = gho * 2*M_PI; // 0.005
    
    param.Ti = 0.0;
    param.Tf = 5000.0; // 5000
    param.dT = 0.25;

    param.N = 2;
    param.h_td_A = fA * 2 * M_PI;
    param.h_td_w = (ho_w[0] - fw * ho_w[1]); // XXX

    param.H0_func = (hamiltonian_func_t)hamiltonian_2ho;
    param.H1_func = (hamiltonian_func_t)hamiltonian_2ho_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_t_2ho;
    
    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    //for (i = 0; i < param.N; i++)
    //    qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_add(qs, 6);  //8
    qdpack_hilbert_space_add(qs, 25); //40
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;


    /*
     * Set up dissipation
     *
     */
    //param.wT = 0.5 * 2 * M_PI;
    //param.wT = 0.5 * 2 * M_PI;
    param.wT = 0.5 * 2 * M_PI;
    param.do_n = 0;
    param.do_a[param.do_n] = operator_ho_lowering(qs, 0, 0); // ho 0
    Navg = 1/(exp(simulation.ho_w[0]/param.wT) - 1);
    // use a factor of 10 increase to speed up the initial decay to ss, before coupling and driving is applied
    param.do_g1[param.do_n] = 10*g0 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 10*g0 * Navg;       // excitation
    param.do_n++;

     param.do_a[param.do_n] = operator_ho_lowering(qs, 1, 0); // ho 1
    Navg = 1/(exp(simulation.ho_w[1]/param.wT) - 1);
    param.do_g1[param.do_n] = 0.1 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = 0.1 * Navg; // excitation
    param.do_n++;

    /* 
     * Pre-calc of operators to calculate expectation values for.
     */
    param.N_op[0] = operator_ho_N(qs, 0, 0);
    param.N_op[1] = operator_ho_N(qs, 1, 0);
    //param.a  = operator_ho_lowering(qs, 0);
    //param.ad = operator_ho_raising(qs, 0);
    //param.b  = operator_ho_lowering(qs, 1);
    //param.bd = operator_ho_raising(qs, 1);
    //param.b2  = qdpack_operator_alloc(qsn, qsn);
    //param.bd2  = qdpack_operator_alloc(qsn, qsn);
    //{
    //    qdpack_complex alpha, beta;
    //    QDPACK_SET_COMPLEX(&alpha, 1.0, 0.0);
    //    QDPACK_SET_COMPLEX(&beta, 0.0, 0.0);
    //    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, param.b,  param.b,  beta, param.b2);
    //    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, param.bd, param.bd, beta, param.bd2);
    //}

    snprintf(param.simsig, sizeof(param.simsig), "_N_%d_M_%d_w0_%.2h_td_w1_%.2h_td_wT_%.2f_fA_%.2f_fw_%.2f_g0_%.4f_lambda_%.3f_C4", 
         qs->nstates[0], qs->nstates[1], simulation.ho_w[0]/(2*M_PI), simulation.ho_w[1]/(2*M_PI), param.wT/(2*M_PI),
         param.h_td_A/(2*M_PI), param.h_td_w/(2*M_PI), g0, gho);
    
    /*
     * Setup initial state
     */
    //rho1 = qdpack_dm_pure_TLS(0.8, 0.2);         // |> state
    //rho2 = qdpack_dm_pure_TLS(0.3, 0.7);         // |> state
    rho1 = qdpack_dm_fock_state(0, qs->nstates[0]);     // |0> state
    rho2 = qdpack_dm_fock_state(0, qs->nstates[1]);     // |0> state
    
        qdpack_operator_list_init(&dm_list);
        qdpack_operator_list_append(&dm_list, rho1);
        //for (i = 1; i < param.N; i++)
    qdpack_operator_list_append(&dm_list, rho2);

    rho0 = qdpack_operator_tensor(qs, &dm_list);

    /*
     * Evolve the quantum system.
     * 
     */
    if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    /*
     * Calculate final effective temperatures for HOs.
     *
     */
    {
        char filename[512], row[512];
        FILE *f;

        snprintf(filename, sizeof(filename), "data.2ho_cool/Teff%s.dat", param.simsig);
        f = fopen(filename, "w");
        
        wTeff[0] = simulation.ho_w[0] / log(1/param.Nexpt[0] + 1);
        wTeff[1] = simulation.ho_w[1] / log(1/param.Nexpt[1] + 1);

        fprintf(f, "wT = %f\n", param.wT);
        fprintf(f, "<N0(tf)> = %f\n", param.Nexpt[0]);
        fprintf(f, "<N1(tf)> = %f\n", param.Nexpt[1]);

        fprintf(f, "<N0_eq> = %f\n", 1/(exp(simulation.ho_w[0]/param.wT) - 1));
        fprintf(f, "<N1_eq> = %f\n", 1/(exp(simulation.ho_w[1]/param.wT) - 1));
    
        fprintf(f, "wTeff[0] = %f\n", wTeff[0]);
        fprintf(f, "wTeff[1] = %f\n", wTeff[1]);

        fprintf(f, "w[0] = %f\n", simulation.ho_w[0]);
        fprintf(f, "w[1] = %f\n", simulation.ho_w[1]);

        fprintf(f, "w[1] / w[0]\t= %f\n", simulation.ho_w[1] / simulation.ho_w[0]);
        fprintf(f, "wTeff[1] / wTeff[0]\t= %f\n", wTeff[1] / wTeff[0]);

        snprintf(filename, sizeof(filename), "data.2ho_cool/Teff.dat");
        snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", 
             simulation.ho_w[1]/(2*M_PI), fA, param.h_td_w / (2*M_PI), g0, wTeff[0], wTeff[1], 
             param.Nexpt[0], param.Nexpt[1], simulation.ho_w[1]/simulation.ho_w[0], wTeff[1]/wTeff[0],
             simulation.ho_w[1]/simulation.ho_w[0] - wTeff[1]/wTeff[0]);
        qdpack_file_row_add(filename, row);    
    }
    
    return 0;
}



