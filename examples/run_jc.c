//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//
//------------------------------------------------------------------------------

//
// Calculate the time evolution of a two-level system coupled to an oscillator
// (Jaynes-Cumming model)
//
// $ mkdir -p data/jc
// $ time ./run_jc

#include <math.h>
#include <stdio.h> 

#include <qdpack/qdpack.h>

// Directory where all the data files are to be stored.
#define DATADIR "data/jc/"

// pointers to operators to calculate expectation values for
static qdpack_operator_t *N_op[2];

// 
// Process the density matrix for time t: calculate expectation values for the
// precomputed operators N_op[]
//
int
rho_store_cb(qdpack_simulation_t *sp,
             qdpack_hilbert_space_t *qs,
             qdpack_operator_t *rho_t, double t)
{
    qdpack_complex expt;
    char filename[1024], row[16384];
      
    // qubit
    expt = qdpack_operator_expectation_value(rho_t, N_op[0]);
    snprintf(filename, sizeof(filename), "%s/qubit_excitation_prob_%s.dat", DATADIR, sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(expt), QDPACK_IMAG(expt));
    qdpack_file_row_add(filename, row);        

    // cavity
    expt = qdpack_operator_expectation_value(rho_t, N_op[1]);
    snprintf(filename, sizeof(filename), "%s/cavity_excitation_prob_%s.dat", DATADIR, sp->simsig);
    snprintf(row, sizeof(row), "%f\t%f\t%f", t, QDPACK_REAL(expt), QDPACK_IMAG(expt));
    qdpack_file_row_add(filename, row);                
       
    return 0;
}


// 
// Program starts here.
//
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho_q, *rho_c, *rho0;
    qdpack_operator_list_t dm_list;
    h_jc_parameters_t params;
    double Gamma_c, Gamma_q, Navg;

    // This data structure contains all the parameters needed by the solver.
    qdpack_simulation_t sim;
    qdpack_simulation_init(&sim);

    sim.parameters = (void *)&params;
    memset((void *)&params, 0, sizeof(h_jc_parameters_t));
   
    //
    // Setup system (hamiltonian) parameter.
    //
    
    params.epsilon[0] = 1.0  * 2 * M_PI; // qubit energy splitting
    params.delta[0]   = 0.0  * 2 * M_PI; // qubit tunnelling rate
    params.ho_w[0]    = 1.0  * 2 * M_PI; // oscillator energy
    params.g_qc       = 0.05 * 2 * M_PI; // coupling strength
    params.N          = 50;              // # fock states in the oscillator
  
    Gamma_c = 0.005;  // cavity decay rate
    Gamma_q = 0.05;   // qubit  decay rate
   
    //
    // Setup solver parameters
    //    
    sim.Ti = 0.0;
    sim.Tf = 25.0;
    sim.dT = sim.Tf / 100.0;

    // select a pre-defined hamiltonian for JC model
    sim.H0_func = (hamiltonian_func_t)hamiltonian_jc;
    
    // 
    // create a new hilbert space object
    //
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_add(qs, params.N);
    qdpack_hilbert_space_print(qs);

    //
    // create a "simulation signature" which is used in as a template for the
    // filenames for the output data
    //
    snprintf(sim.simsig, sizeof(sim.simsig),
            "Nq_%d_Nc_%d_wq_%.3f_who_%.3f_g_%.3f_Gq_%.3f_Gc_%.4f",
            qs->nstates[0], qs->nstates[1],
            sqrt(params.epsilon[0]*params.epsilon[0] + params.delta[0]*params.delta[0])/(2*M_PI),
            params.ho_w[0]/(2*M_PI),
            params.g_qc/(2*M_PI), Gamma_q, Gamma_c);

    //
    // Setup dissipation
    //
    Navg = 0;

    /* --- qubit --- */
    sim.do_a[sim.do_n]  = qdpack_operator_alloc(qs);
    operator_ho_lowering(sim.do_a[sim.do_n], qs, 0, 0); 
    sim.do_g1[sim.do_n] = Gamma_q * (1 + Navg); // relaxation
    sim.do_g2[sim.do_n] = Gamma_q * Navg;       // excitation
    sim.do_n++;
        
    /* --- cavity --- */
    sim.do_a[sim.do_n]  = qdpack_operator_alloc(qs);
    operator_ho_lowering(sim.do_a[sim.do_n], qs, 1, params.n_offset); 
    sim.do_g1[sim.do_n] = Gamma_c * (1 + Navg); // relaxation
    sim.do_g2[sim.do_n] = Gamma_c * Navg;       // excitation
    sim.do_n++;

    // 
    // Pre-calc of operators to calculate expectation values for.
    //
    N_op[0] = qdpack_operator_alloc(qs);
    operator_ho_N(N_op[0], qs, 0, 0);
    N_op[1] = qdpack_operator_alloc(qs);
    operator_ho_N(N_op[1], qs, 1, params.n_offset);
 
    //
    // Setup initial state
    //
    rho_q = qdpack_dm_pure_TLS(1.0);                 
    rho_c = qdpack_dm_fock_state(0, qs->nstates[1]); 

    // create a composite initial density matrix
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho_q);
    qdpack_operator_list_append(&dm_list, rho_c);
    rho0 = qdpack_operator_alloc(qs);
    qdpack_operator_tensor(rho0, qs, &dm_list);
        
    //
    // Evolve the quantum system until time t = T
    //
    printf("Start evolving system\n");
    if (qdpack_evolve_dm_lme(&sim, qs, rho0, rho_store_cb) != 0)
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }
  
    return 0;
}



