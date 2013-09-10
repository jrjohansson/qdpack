//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// The steady-state of a 2LS subject to a strong driving field: Find the steady
// state by finding the eigenstates of the propagator.
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

#define DATADIR "data/2ls.pss"

static double x, y;
static qdpack_operator_t *rho_eb_t, *rho_tmp_t;

/**
 * Strong driving hamiltonian :  H_t = H0 + f(t) H1, f(t) = A * sin(w * t)
 *
 */
int
hamiltonian_sd_t(qdpack_simulation_t *p,
                 qdpack_operator_t *H_t,
                 qdpack_operator_t *H0,
                 qdpack_operator_t *H1,
                 double t)
{
    if (p->h_td_A == 0.0)
    {
        qdpack_operator_memcpy(H_t, H0);
    }
    else
    {
        qdpack_operator_memcpy(H_t, H1);
        qdpack_operator_scale(H_t, qdpack_complex_rect(p->h_td_A * sin(p->h_td_w * t) / 2, 0.0));
        qdpack_operator_add(H_t, H0);
    }

    return 0;
}

//
// Process the density matrix for time t=T_final
//
int
process_rho_final(qdpack_simulation_t *sim,
                  qdpack_hilbert_space_t *qs,
                  qdpack_operator_t *rho_cb_t,
                  double t)
{
    char filename[128], row[1024];
    
    //
    // calculate the instantaneous eigenbasis
    //
    sim->Ht_func(sim, sim->H_t, sim->H0, sim->H1, sim->Tf);

    qdpack_matrix_eigen_hermv(sim->H_t->data, sim->eval, sim->evec->data, QDPACK_EIGEN_SORT_VAL_ASC);

    qdpack_operator_blas_zgemm(QDpackConjTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, sim->evec,  rho_cb_t, QDPACK_COMPLEX_ZERO, rho_tmp_t);
    qdpack_operator_blas_zgemm(QDpackNoTrans,   QDpackNoTrans, QDPACK_COMPLEX_ONE, rho_tmp_t, sim->evec, QDPACK_COMPLEX_ZERO, rho_eb_t);

    //
    // store the qubit excitation probability in both charge basis (rho_cb_t) 
    // and eigenbasis (rho_eb_t)
    //
    snprintf(filename, sizeof(filename), "%s/p_e_%s.dat", DATADIR, sim->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f\t%f",
             x, y, 
             QDPACK_REAL(qdpack_operator_get(rho_cb_t, 1, 1)),
             QDPACK_REAL(qdpack_operator_get(rho_eb_t, 1, 1)));

    qdpack_file_row_add(filename, row);        

    return 0;
}

//
// Program starts here.
// 
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs, *qs2;
    qdpack_operator_t *rho_ss, *U;
    int qsn;
    double Navg, fw, delta, G1, G2;
    double xi, xf, xd, yi, yf, yd;

    h_qubits_parameters_t qubit_params;

    // This data structure store all the parameters needed by the solver.
    qdpack_simulation_t sim;
    qdpack_simulation_init(&sim);

    sim.parameters = (void *)&qubit_params;
    memset((void *)&qubit_params, 0, sizeof(h_qubits_parameters_t));

    if (argc != 7)
    {
        printf("usage: %s x-init x-step x-final y-init y-step y-final\n", argv[0]);
        printf("\nexammple: $ %s -10.0 0.2 10.0 0.0 0.1 10.0\n", argv[0]);
        exit(0);
    }

    //
    // parameters
    //

    xi     = strtod(argv[1], NULL);
    xd     = strtod(argv[2], NULL);
    xf     = strtod(argv[3], NULL);

    yi     = strtod(argv[4], NULL);
    yd     = strtod(argv[5], NULL);
    yf     = strtod(argv[6], NULL);
    
    fw     = 2.0;   // driving frequency
    delta  = 0.1;   // coefficient of sigma_x in qubit hamiltonian
    Navg   = 0;     // zero temperature

    G1 = 0.00001;
    G2 = 0.001;

    //
    // hamiltonian
    //
    sim.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    sim.H1_func = (hamiltonian_func_t)hamiltonian_z_drive;
    sim.Ht_func = (hamiltonian_td_func_t)hamiltonian_sd_t;

    //
    // create a hilbert space
    //
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);        // a 2-level system
    qsn = qdpack_hilbert_space_nstates(qs);

    qs2 = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs2, 4);        // a 2-level system

    //
    // Set up dissipation
    //

    // relaxation
    if (G1 > 0.0)
    {
        sim.do_a[sim.do_n]  = qdpack_operator_alloc(qs);
        operator_ho_lowering(sim.do_a[sim.do_n], qs, 0, 0); 
        sim.do_g1[sim.do_n] = G1 * (1 + Navg); // relaxation
        sim.do_g2[sim.do_n] = G1 * Navg;       // excitation
        sim.do_n++;
    }

    // dephasing
    if (G2 > 0.0)
    {
        sim.do_a[sim.do_n]  = qdpack_operator_alloc(qs);
        operator_sigma_z(sim.do_a[sim.do_n], qs, 0);
        sim.do_g1[sim.do_n] = G2 / 2; 
        sim.do_g2[sim.do_n] = G2 / 2; 
        sim.do_n++;
    }

    // output data file name signature
    snprintf(sim.simsig, sizeof(sim.simsig), "dDelta_%.3f_G1_%.6f_G2_%.4f_fw_%.3f", delta, G1, G2, fw);

    //
    // pre-allocate memory
    //
    sim.eval  = qdpack_matrix_alloc(qsn, 1);
    sim.evec  = qdpack_operator_alloc(qs);
    U         = qdpack_operator_alloc(qs2);
    rho_eb_t  = qdpack_operator_alloc(qs);
    rho_tmp_t = qdpack_operator_alloc(qs);
    rho_ss    = qdpack_operator_alloc(qs);

    //
    // start sweep of parameter range
    //
    qubit_params.delta[0] = delta * 2 * M_PI; // qubit parameter: sigma_x coefficient
    sim.h_td_w            =    fw * 2 * M_PI; // driving frequency

    for (x = xi; x <= xf; x += xd)
    {    
        qubit_params.epsilon[0] = x * (2*M_PI); // qubit parameter: sigma_z coefficient

        printf("ITERATION: %s: x = %f, y = %.2f:%.2f:%.2f\n", sim.simsig, x, yi, yd, yf);

        for (y = yi; y <= yf; y += yd)
        {   
            sim.h_td_A = y * 2 * M_PI; // driving amplitude

            // calculate propagator for one period
            sim.Ti = 0.0; sim.Tf = 1/fw; sim.dT = 1/fw;
            if (qdpack_propagator_dm(&sim, qs, U, NULL) == -1)
            {
                fprintf(stderr, "Evaluation of the propagator failed.\n");
                return -1;
            }
    
            // Find the steady state from the evolution operator
            if (qdpack_steadystate_dm_propagator(&sim, qs, U, rho_ss) == -1)
            {
                fprintf(stderr, "Evaluation of the propagator failed.\n");
                return -1;
            }
        
            process_rho_final(&sim, qs, rho_ss, sim.Tf);
        }
    }

    return 0;
}



