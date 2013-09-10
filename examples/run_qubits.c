//------------------------------------------------------------------------------ 
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Simulate the time evolution of a system of N two-level systems (TLS).
//
// $ mkdir -p data/qubit.dynamics
// $ 
//

#include <math.h>
#include <stdio.h> 

#include <qdpack/qdpack.h>

/**
 * Directory where all the data files are to be stored,
 * relative to the directory from which the program is
 * started.
 */
#define DATADIR "data/qubit.dynamics/"

/** 
 * Process the density matrix for time t (calc. expectations values, 
 * store to file, etc.). This is a callback function that the solver
 * call in each time-step.
 *
 * Here we can calculate expectation values etc. In this case we only
 * do simple ...
 */
int
rho_store_cb(qdpack_simulation_t *sp, qdpack_hilbert_space_t *qs,
             qdpack_operator_t *rho_t, double t)
{
    int i;
    char filename[128], row[512];
    qdpack_operator_t *dm_part;
    
    for (i = 0; i < qs->nsubsys; i++)
    {
        /*
         * For each subsystem's density matrix, extract the diagonal elements
         * that is the qubits occupation probabilities, and store to data
         * files.
         */               
        dm_part = qdpack_operator_traceout(qs, rho_t, i);

        snprintf(filename, sizeof(filename), "%s/qubit_population_%d.dat", DATADIR, i);
        snprintf(row, sizeof(row), "%f\t%f\t%f", t, 
                 QDPACK_REAL(qdpack_operator_get(dm_part, 0, 0)),    // rho^{qubit i}_{0,0}
                 QDPACK_REAL(qdpack_operator_get(dm_part, 1, 1)));   // rho^{qubit i}_{1,1}
        
        qdpack_file_row_add(filename, row);

        qdpack_operator_free(dm_part);
    }

    return 0;
}

/** 
 * Program starts here: Set up and start a calculation that simulates
 * the time-evolution of a coupled qubits.
 */
int
main(int argc, char **argv)
{
    qdpack_hilbert_space_t *qs;
    qdpack_operator_t *rho1, *rho2;
    qdpack_operator_t *rho0;
    qdpack_operator_list_t dm_list;
    int i;

    h_qubits_parameters_t qubit_params;

    // This data structure store all the parameters needed by the solver.
    qdpack_simulation_t sim;
    qdpack_simulation_init(&sim);

    sim.parameters = (void *)&qubit_params;
    memset((void *)&qubit_params, 0, sizeof(h_qubits_parameters_t));

    /*
     * Setup solver parameters: 
     *
     * Which simeters in qdpack_simulation_t that are used depends
     * on which hamiltonian generating function that is used
     * later on. Here we will use the hamiltonian for coupled
     * two-level systems, so we configure the following simeters
     * accordingly:
     *
     */

    /* Number of two-level systems (qubits) */
    qubit_params.N = 2;

    /* TLS hamiltonian parameters, i.e. 
     * 
     *     - energy splitting:     qubit_params.epsilon[i]
     *     - tunneling energy:     qubit_params.delta[i]
     *     - inter-qubit coupling: qubit_params.lambda[i][j]
     *
     * where i and j is two qubit indices. The actual numerical
     * values are defined in the arrays above.
     */

    
    qubit_params.epsilon[0] = 1.0 * 2*M_PI;
    qubit_params.delta[0]   = 0.0 * 2*M_PI;

    qubit_params.epsilon[1] = 1.0 * 2*M_PI;
    qubit_params.delta[1]   = 0.0 * 2*M_PI;

    /* qubits are coupled in a chain */
    qubit_params.lambda[0][1] = 0.100 * 2*M_PI;
    qubit_params.lambda[1][2] = 0.050 * 2*M_PI;
    
    
    /* Configure the time range for the evolution:
     *  - Ti = initial time
     *  - Tf = final time
     *  - dT = time-step between each call to the rho_store
     *         callback function. The solver itself uses
     *         an adaptive step length for accuracy control.
     */
    sim.Ti =  0.0;
    sim.Tf = 20.0;
    sim.dT =  0.01;

    /*
     * If we want to apply a driving field to the qubits (i.e.
     * for dynamic coupling, rabi oscillations, etc.), the following
     * two variables define the amplitude and frequency of the 
     * driving field. The operator that the driving field couple 
     * to is defined by the hamiltonian_drive function, see below.
     *
     * Set the amplitude h_td_A to zero to disable driving.
     */
    sim.h_td_A = 1.5  * 2 * M_PI;
    sim.h_td_w = 10.0 * 2 * M_PI;


    /* Function pointers to the functions that generate the Hamiltonian
     * matrices. Ht_func is called by the solver and expected to return
     * the instantanous value of the Hamiltonian (at time t). This function
     * can be customized for each simulation, but often it is sufficient
     * to use the standard from:
     *
     * H(t) = H_0 + H_1 * h_td_A * cos(h_td_w * t)
     *
     * which is implemented by the function hamiltonian_t from the 
     * QDpack's library of standard Hamiltonians.
     *
     */
    sim.H0_func = (hamiltonian_func_t)hamiltonian_qubits;
    sim.H1_func = (hamiltonian_func_t)hamiltonian_drive;
    sim.Ht_func = (hamiltonian_td_func_t)hamiltonian_t;


    /* Create a new quantum system object: Here we configure
     * the quantum systems that make up the total system. We
     * call qdpack_hilbert_space_add for each subsystem, and the 
     * second simeter defines how many quantum state that
     * subsystem has.
     *
     */
    qs = qdpack_hilbert_space_new();
    for (i = 0; i < qubit_params.N; i++)
        qdpack_hilbert_space_add(qs, 2);
    qdpack_hilbert_space_print(qs);


    /*
     * Set up dissipation: Here we can define operators and rates 
     * that appear in the master equation (relaxation rate, dephasing
     * rates, etc.) for each subsystem.
     *
     */
    sim.do_n = 0;

    // QUBIT 0:
    // no dissipation

    // QUBIT 1: 
    // Relaxation
    sim.do_a[sim.do_n] = qdpack_operator_alloc(qs);
    operator_ho_lowering(sim.do_a[sim.do_n], qs, 1, 0);     // qubit 1
    sim.do_g1[sim.do_n] = 0.075;  // relaxation rate
    sim.do_g2[sim.do_n] = 0.00;   // excitation rate (zero temperature)
    sim.do_n++;

    // Pure dephasing
    //sim.do_a[sim.do_n] = qdpack_operator_alloc(qs)
    // operator_sigma_z(sim.do_a[sim.do_n], qs, 1, 0);     // qubit 1
    // sim.do_g1[sim.do_n] = 0.25; // dephasing rate
    // sim.do_n++;
    
    // QUBIT 2:
    // Relaxation
    sim.do_a[sim.do_n] = qdpack_operator_alloc(qs);
    operator_ho_lowering(sim.do_a[sim.do_n] , qs, 2, 0);     // qubit 2
    sim.do_g1[sim.do_n] = 0.25;   // relaxation rate
    sim.do_g2[sim.do_n] = 0.05;   // excitation rate 
    sim.do_n++;
    

    /*
     * Setup initial state: Here we set the initial state for all qubits
     * to its ground state.
     */
    rho1 = qdpack_dm_pure_TLS(1.0); // = density matrix that correspond to |psi> = 1.0 |0> + 0.0 |1>
    rho2 = qdpack_dm_pure_TLS(0.0); // = density matrix that correspond to |psi> = 1.0 |0> + 0.0 |1> 
    
    qdpack_operator_list_init(&dm_list);
    qdpack_operator_list_append(&dm_list, rho1);
    for (i = 1; i < qubit_params.N; i++)
        qdpack_operator_list_append(&dm_list, rho2);

    rho0 = qdpack_operator_alloc(qs);
    qdpack_operator_tensor(rho0, qs, &dm_list);
    
    printf("rho0 =\n");
    qdpack_operator_print(rho0);


    /*
     * Evolve the quantum system: We can use different solvers
     * for different kind of evolution. Here we use the LME solver
     * (Lindblad Master Equation).
     *
     * The solver will call the callback function "rho_store_cb" for
     * each time-step. In the callback function we can calculate
     * expectation values, etc., and store to files.
     * 
     */
    if (qdpack_evolve_dm_lme(&sim, qs, rho0, rho_store_cb) != 0)           // dissipative evolution: Lindblad Master Equation
    //if (qdpack_evolve_dm_unitary_const(&sim, qs, rho0, rho_store_cb) != 0)    // unitary evolution with constant Hamiltonian
    //if (qdpack_evolve_dm_unitary_t(&sim, qs, rho0, rho_store_cb) != 0)        // unitary evolution with time-dependent Hamiltonian
    {
        fprintf(stderr, "Evolution of the quantum system failed.\n");
        return -1;
    }

    return 0;
}



