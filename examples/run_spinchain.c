//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//
//------------------------------------------------------------------------------

//
// calculate the dynamics of a spin-1/2 chain
//

#include <math.h>
#include <stdio.h> 

#include <qdpack/qdpack.h>

/**
 * Directory where all the data files are to be stored,
 * relative to the directory from which the program is
 * started.
 */
#define DATADIR "data/spinchain/"

static qdpack_operator_t *N_op[10]; // XXX

/* ---
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.) 
 */
int
rho_store_cb(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs,
             qdpack_operator_t *rho_t, double t)
{
    double ent_ln, ent_ln_sum = 0.0, ent_ne, ent_ne_sum = 0.0;
    int i, j, row1_offset, row2_offset;
    char filename1[128], filename2[128], row1[4096], row2[4096];
    qdpack_operator_t *dm_part;
    qdpack_hilbert_space_t *qs_mask;


    //
    // calculate entanglement between qubits
    //

    qs_mask = qdpack_hilbert_space_copy(qs);
    for (i = 0; i < qs_mask->nsubsys; i++)
    {
        qs_mask->nstates[i] = 1;
    }

    snprintf(filename1, sizeof(filename1), "%s/ent_int_pairwise_%s.dat", DATADIR, sim->simsig);
    snprintf(filename2, sizeof(filename2), "%s/ent_ext_pairwise_%s.dat", DATADIR, sim->simsig);
    row1_offset = snprintf(row1, sizeof(row1), "%f", t);
    row2_offset = snprintf(row2, sizeof(row2), "%f", t);

    for (i = 0; i < qs->nsubsys-1; i++)
    {
        qs_mask->nstates[i] = 0;
        for (j = i+1; j < qs->nsubsys; j++)
        {
            qs_mask->nstates[j] = 0;
        
            dm_part = qdpack_operator_traceout_multi(qs, rho_t, qs_mask);
            ent_ln = qdpack_entanglement_log_neg(dm_part);        
            ent_ln_sum += ent_ln;    
            row1_offset += snprintf(&row1[row1_offset], sizeof(row1)-row1_offset, "\t%f", ent_ln); 
            ent_ne = qdpack_entanglement_neumann_entropy(dm_part);        
            ent_ne_sum += ent_ne;    
            row2_offset += snprintf(&row2[row2_offset], sizeof(row2)-row2_offset, "\t%f", ent_ne); 
        
            qdpack_operator_free(dm_part);
            qs_mask->nstates[j] = 1;
        }
        qs_mask->nstates[i] = 1;
    }
    
    snprintf(&row1[row1_offset], sizeof(row1)-row1_offset, "\t%f", ent_ln_sum); 
    qdpack_file_row_add(filename1, row1);
    snprintf(&row2[row1_offset], sizeof(row2)-row2_offset, "\t%f", ent_ln_sum); 
    qdpack_file_row_add(filename2, row2);

    /* --- expectation values -- */

    for (i = 0; i < ((h_qubits_parameters_t *)sim->parameters)->N; i++)
    {
        qdpack_complex expt;

        expt = qdpack_operator_expectation_value(rho_t, N_op[i]);
        snprintf(filename1, sizeof(filename1), "%s/qubit_occupation_%d_%s.dat", DATADIR, i, sim->simsig);
        snprintf(row1, sizeof(row1), "%f\t%f\t%f", t, QDPACK_REAL(expt), QDPACK_IMAG(expt));
        qdpack_file_row_add(filename1, row1);                
    }       
    
    return 0;
}

/* 
 * 
 * 
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho1, *rho2;
    qdpack_operator_t *rho0;
    qdpack_operator_list_t dm_list;
    int i;
    double G1, Navg;

    h_qubits_parameters_t qubit_params;

    // This data structure store all the parameters needed by the solver.
    qdpack_simulation_t sim;
    qdpack_simulation_init(&sim);

    sim.parameters = (void *)&qubit_params;
    memset((void *)&qubit_params, 0, sizeof(h_qubits_parameters_t));

    /*
     * Setup solver parameters.
     *
     */

    qubit_params.N = 6;

    /* TLS hamiltonian parameters, i.e. 
     * 
     *     - energy splitting:     qubit_params.epsilon[i]
     *     - tunneling energy:     qubit_params.delta[i]
     *     - inter-qubit coupling: qubit_params.lambda[i][j]
     *
     * where i and j is two qubit indices. The actual numerical
     * values are defined in the arrays above.
     */

    for (i = 0; i < qubit_params.N; i++)
    {
        qubit_params.epsilon[i] = 1.0  * 2*M_PI;
        qubit_params.delta[i]   = 0.0  * 2*M_PI;
    }

    for (i = 0; i < qubit_params.N-1; i++)
    {
        qubit_params.lambda[i][i+1] = 0.05 * 2*M_PI;
    }

   
    // Configure the time range for the evolution:
    //  - Ti = initial time
    //  - Tf = final time
    //  - dT = time-step between each call to the rho_store
    //         callback function. The solver itself uses
    //         an adaptive step length for accuracy control.
    //
    sim.Ti =  0.0;
    sim.Tf = 50.0;
    sim.dT =  0.01;

    // Select the pre-defined coupled-qubit hamiltonian
    //
    sim.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    
    // Create a new quantum system object: Here we configure
    // the quantum systems that make up the total system. We
    // call qdpack_hilbert_space_add for each subsystem, and the 
    // second simeter defines how many quantum state that
    // subsystem has.
    //
    qs = qdpack_hilbert_space_new();
    for (i = 0; i < qubit_params.N; i++)
        qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_print(qs);

    //
    // Simulation signature to be used in file name for output data files
    //
    snprintf(sim.simsig, sizeof(sim.simsig), "N_%d", qubit_params.N);

    //
    // Set up dissipation
    //
    //
    G1 = 0.01;
    for (i = 0; i < qubit_params.N; i++)
    {
        Navg = 0.0;

        if (G1 > 0.0)
        {
            sim.do_a[sim.do_n] = qdpack_operator_alloc(qs);
            operator_ho_lowering(sim.do_a[sim.do_n], qs, i, 0);
            sim.do_g1[sim.do_n] = (Navg + 1) * G1;  // relaxation rate
            sim.do_g2[sim.do_n] = Navg * G1;        // excitation rate 
            sim.do_n++;
        }
    }
    
    //
    // Setup initial state
    //
    rho1 = qdpack_dm_pure_TLS(1.0); 
    rho2 = qdpack_dm_pure_TLS(0.0); 
    
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho1);
    for (i = 1; i < qubit_params.N; i++)
        qdpack_operator_list_append(&dm_list, rho2);

    rho0 = qdpack_operator_alloc(qs);
    qdpack_operator_tensor(rho0, qs, &dm_list);

    //
    // pre-calculate operators for expectation values
    //
    for (i = 0; i < qubit_params.N; i++)
    {
        N_op[i] = qdpack_operator_alloc(qs);
        //operator_ho_N(N_op[i], qs, i, 0); // number operator for qubit i
        operator_sigma_x(N_op[i], qs, i);
    }

    printf("starting time evolution... \n");

    //
    // Evolve the quantum system.
    // 
    if (qdpack_evolve_dm_lme(&sim, qs, rho0, rho_store_cb) != 0) 
    //if (qdpack_evolve_dm_lme_t(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_const(qs, rho0, &param, rho_store_cb) != 0)
    //if (qdpack_evolve_dm_unitary_t(qs, rho0, &param, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    return 0;
}



