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

static char DATADIR[256]; 
static double pe_0_expt_sum = 0.0, pe_0_max = 0.0;
static int    pe_0_expt_no  = 0;
static double pe_1_expt_sum = 0.0;
static int    pe_1_expt_no  = 0;
#define DATADIR_TEMPLATE "/home/rob/qdpack-data/%s"

/**
 * Strong driving hamiltonian
 *
 */
int
hamiltonian_qubit_t(qdpack_operator_t *H_t, qdpack_operator_t *H0, qdpack_operator_t *H1, double t, qdpack_simulation_t *p)
{
        if (p->h_td_A == 0.0)
        {
                qdpack_operator_memcpy(H_t, H0);
                // p->H_t = p->H_0;
        }
        else
        {
                qdpack_complex z;

        // H_t = H0 + f(t) H1
                // f(t) = A * sin(w * t)

                //qdpack_operator_memcpy(H_t, H0);
                qdpack_operator_memcpy(H_t, H1);
                QDPACK_SET_COMPLEX(&z, p->h_td_A * cos(p->h_td_w * t) / 2, 0.0);
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
    char filename[128], row[1024];
    
    snprintf(filename, sizeof(filename), "%s/p_e%s.dat", DATADIR, sp->simsig);
    
    /* -- expectation values -- */
    snprintf(row, sizeof(row), "%f\t%f\t%f\t%f\t%f\t%f", sp->x, sp->y, sp->epsilon[0]/(2*M_PI),
         pe_0_expt_sum / pe_0_expt_no, pe_0_max, pe_1_expt_sum / pe_1_expt_no);
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
    qdpack_operator_t *dm_part; //, *op;

    /*
     * For each sub-system:
     *
     */
    for (i = 0; i < qs->nsubsys; i++)
        {
        ri = 0;    

        //printf("Procssing subsystem [%d] at time %f\n", i, t);
        
        //snprintf(filename, sizeof(filename), "%s/rho_%d_diag%s_x_%.3f_y_%.3f.dat", DATADIR, i, sp->simsig, sp->x, sp->y);
        
                dm_part = qdpack_operator_traceout(qs, rho_t, i);

        //ri = snprintf(row, sizeof(row), "%f", t);
        //
        //for (j = 0; j < qs->nstates[i]; j++)
        //{
        //     ri += snprintf(&row[ri], sizeof(row)-ri, "\t%f", (double)QDPACK_REAL(qdpack_operator_get(dm_part, j, j)));
        //}
        //qdpack_file_row_add(filename, row);
        
        if (i == 0)
        {
            pe_0_expt_sum += (double)QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1));
            pe_0_expt_no++;
            if ((double)QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1)) > pe_0_max)
            {
                pe_0_max = (double)QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1));
            }
        }
        if (i == 1)
        {
            pe_1_expt_sum += (double)QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1));
            pe_1_expt_no++;
        }

        qdpack_operator_free(dm_part);
    }    
    
    return 0;
}


/* ---
 * Simulation parameters.
 * 
 */
double epsilon[] = { 0.0, 8.8, 9.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double delta[]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double ho_w[]    = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
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
    double Navg, fw, fA, phi, Gc, Gq, gqc, x, y, pDelta, G1, G2, w0 = 27.3345, Ic = 11.659;
    double yi, yf, yd;

    qdpack_simulation_t param;
    qdpack_simulation_init(&param);

    if (argc != 5)
    {
        printf("usage: %s check code\n", argv[0]);
        exit(0);
    }

    snprintf(DATADIR, sizeof(DATADIR), DATADIR_TEMPLATE, argv[0]);
    
    x      = strtod(argv[1], NULL);
    //y      = strtod(argv[2], NULL);
    yi     = strtod(argv[2], NULL);
    yd     = strtod(argv[3], NULL);
    yf     = strtod(argv[4], NULL);
    
    param.x = x;

    pDelta = 0.0;
    fw = x;
    G1 = 0.0;
    G2 = 0.0;

    param.epsilon = epsilon;
    param.delta   = delta;
    simulation.ho_w    = ho_w;
    for (i = 0; i < 10; i++)
    {
        param.epsilon[i] *= 2*M_PI;
        param.delta[i]   *= 2*M_PI;
        simulation.ho_w[i]    *= 2*M_PI;
        lambda[i]         = (double *)calloc(10, sizeof(double));
    }
    param.lambda = lambda;
    param.lambda[0][1] = 0.000 * 2 * M_PI;
    param.lambda[0][2] = 0.000 * 2 * M_PI;

    param.Ti = 0.0;
    param.dT = 0.5;
    param.Tf = 500.0;
    param.cont_basis_transform = 0;

    param.h_td_A = 0.01 * 2 * M_PI;
    param.h_td_w =  fw * 2 * M_PI;

    param.H1_func = (hamiltonian_func_t)hamiltonian_x_drive;
    param.Ht_func = (hamiltonian_td_func_t)hamiltonian_qubit_t;

    /* 
     * create a new quantum system object
     */
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
//    qdpack_hilbert_space_add(qs, 2);
//    qdpack_hilbert_space_add(qs, 2);
    qsn = qdpack_hilbert_space_nstates(qs);
    qdpack_hilbert_space_print(qs);
    param.qs = qs;


    /*
     * Set up dissipation
     *
     */
    param.wT = 0.0 * 2 * M_PI;
    param.do_n = 0;
    Navg = 0;


    // relaxation
/*
    param.do_a[param.do_n]  = operator_ho_lowering(qs, 0, 0); 
    param.do_ad[param.do_n] = operator_ho_raising(qs, 0, 0);
    Navg = 0; //1/(exp(sqrt(param.epsilon[0]*param.epsilon[0] + param.delta[0]*param.delta[0])/param.wT) - 1);
    printf("qubit Navg = %f\n", Navg);
    param.do_g1[param.do_n] = G1 * (1 + Navg); // relaxation
    param.do_g2[param.do_n] = G1 * Navg;       // excitation
    param.do_n++;
*/
/*
    // dephasing
    param.do_a[param.do_n]  = operator_sigma_z(qs, 0);
        param.do_ad[param.do_n] = operator_sigma_z(qs, 0); 
        param.do_g1[param.do_n] = G2 / 2; 
        param.do_g2[param.do_n] = G2 / 2; 
        param.do_n++;
*/
    snprintf(param.simsig, sizeof(param.simsig), "_G1_%.3f_G2_%.3f_g1_%.3f_g2_%.3h_td_A_%.3f", 
         G1, G2, lambda[0][1]/(2*M_PI), lambda[0][2]/(2*M_PI), param.h_td_A/(2*M_PI));

    /* storage space for final rho */
    param.rho_f = qdpack_operator_alloc(qsn, qsn);

    /*
     * Setup initial state
     */
        qdpack_operator_list_init(&dm_list);
    rho_q = qdpack_dm_pure_TLS(0.0);          // |0> state
        qdpack_operator_list_append(&dm_list, rho_q); // qubit
  //      qdpack_operator_list_append(&dm_list, rho_q); // tls
   //     qdpack_operator_list_append(&dm_list, rho_q); // tls
    rho0 = qdpack_operator_tensor(qs, &dm_list);

    //delta[0]   = pDelta * fw * (2*M_PI);

    printf("ITERATION: %s: x = %f, y = ...\n", param.simsig, param.x);

    for (y = yi; y < yf; y += yd)
    {    
        pe_0_expt_sum = 0.0;
        pe_0_expt_no  = 0;
        pe_1_expt_sum = 0.0;
        pe_1_expt_no  = 0;
        pe_0_max = 0.0;

        param.y = y;
    
        //epsilon[0] = y * (2*M_PI);
        epsilon[0] = w0 * 2 * M_PI * sqrt( sqrt( 1.0 - (y / Ic) * (y / Ic) ) );
    
        if (abs(epsilon[0] - param.h_td_w) > 0.5*2*M_PI)
            continue;

        param.H0_func = (hamiltonian_func_t)hamiltonian_qubits;

            
        /*
         * Evolve the quantum system until time t = T, where the steady state
             * is assumed to be reached.
         * 
         */
        if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
        {
            fprintf(stderr, "Evolution of the quantum system failed.\n");
            return -1;
        }
    
        rho_store_final_cb(param.rho_f, param.Tf, qs, &param);

    }

    return 0;

}



