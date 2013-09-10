//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#include <stdio.h>

#include <qdpack/qdpack.h>

//------------------------------------------------------------------------------
// Density matrices
//------------------------------------------------------------------------------


/**
 * \breif    Initial density matrix for a two-level systems.
 *
 * Generate a pure-state density matrix for a TLS. (Used for constructing
 * initial states).
 *
 * p_ex = probability of being in the excited state
 *
 * @param    p_ex    Occupation probability of the two-level system excited state.
 */
qdpack_operator_t *
qdpack_dm_pure_TLS(double p_ex)
{
    double aa, bb;
    qdpack_operator_t *dm;
    qdpack_complex z;
    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);

    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate dm (Out of memory?).\n", __PRETTY_FUNCTION__);
        return NULL;
    }

    aa = sqrt(1-p_ex);
    bb = sqrt(p_ex);
    QDPACK_SET_COMPLEX(&z, aa * aa, 0.0); qdpack_operator_set(dm, 0, 0, z);
    QDPACK_SET_COMPLEX(&z, aa * bb, 0.0); qdpack_operator_set(dm, 0, 1, z);
    QDPACK_SET_COMPLEX(&z, bb * aa, 0.0); qdpack_operator_set(dm, 1, 0, z);
    QDPACK_SET_COMPLEX(&z, bb * bb, 0.0); qdpack_operator_set(dm, 1, 1, z);
    
    qdpack_hilbert_space_free(qs);

    return dm;
}


/**
 * \breif    Initial wave function for a two-level systems.
 *
 * Generate a pure-state wave function for a TLS. (Used for constructing
 * initial states).
 *
 * p_ex = probability of being in the excited state
 *
 * @param    p_ex    Occupation probability of the two-level system excited state.
 */
qdpack_state_t *
qdpack_wf_pure_TLS(double p_ex)
{
    double aa, bb;
    qdpack_state_t *wf;
    qdpack_complex z;
    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 2);
    
    if ((wf = qdpack_state_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_wf_pure_TLS: Failed to allocate wf (Out of memory?).\n");
        return NULL;
    }

    aa = sqrt(1-p_ex);
    bb = sqrt(p_ex);
    QDPACK_SET_COMPLEX(&z, aa, 0.0); qdpack_state_set(wf, 0, z);
    QDPACK_SET_COMPLEX(&z, bb, 0.0); qdpack_state_set(wf, 1, z);

    qdpack_hilbert_space_free(qs);
    
    return wf;
}

/**
 * \breif    Initial density matrix for a 3-level systems.
 *
 * Generate a pure-state density matrix for a 3LS. (Used for constructing
 * initial states).
 *
 * p_ex = probability of being in the excited state
 *
 * @param    p_ex    Occupation probability of the two-level system excited state.
 */
qdpack_operator_t *
qdpack_dm_pure_3LS(double p_ex)
{
    //double aa, bb;
    qdpack_operator_t *dm;
    qdpack_complex z;

    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, 3);    

    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_pure_3LS: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }

    qdpack_operator_set_zero(dm);
    
    QDPACK_SET_COMPLEX(&z, 1.0, 0.0);
    qdpack_operator_set(dm, 0, 0, z);

    qdpack_hilbert_space_free(qs); 

    return dm;
}    


/**
 * \breif    Initial density matrix for single-mode boson (thermal state)
 *
 * Generate a density matrix that represent a thermal boson state.
 * 
 * @param    w    Frequency (angular) of the boson field
 * @param    w_th    Temperature in frequency (angular) units.
 * @param    N    Number of fock states
 *
 */
qdpack_operator_t *
qdpack_dm_boson_thermal(double w, double w_th, int N)
{
    qdpack_operator_t *dm;
    qdpack_complex z;
    double Z, p;
    int n;

    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    Z = exp(-w / (2*w_th)) / (1 - exp(-w/w_th));

    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_boson_thermal: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }
    
    qdpack_operator_set_zero(dm);
    
    for (n = 0; n < N; n++)
    {
        p = exp(-w*(n + 0.5) / w_th) / Z;
        QDPACK_SET_COMPLEX(&z, p, 0.0);
        qdpack_operator_set(dm, n, n, z);
    }

    qdpack_hilbert_space_free(qs);
    
    return dm;
}


/**
 * \breif    Initial density matrix for a single-mode boson field (coherent state)
 *
 * Generate a density matrix that represents a coherent state:
 *
 * alpha = r exp(i theta)
 * @param       r        Coherent state amplitude.
 * @param       theta    Coherent state phase.
 * @param       N        Number of fock states.
 *
 */
qdpack_operator_t *
qdpack_dm_coherent_state(double r, double theta, int N, int offset)
{
    int m, n;
    qdpack_operator_t *dm;
    qdpack_complex z;

    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    //printf("DEBUG: coherent state: alpha = %f exp(i*%f), truncated to %d states.\n", r, theta, N);
    
    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_coherent_state: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }

    for (m = 0; m < N; m++)
    {
        for (n = 0; n < N; n++)
        {    
            double dmij = exp(-r * r) * (pow(r, (double)(m+offset)) / fact_sqrt(m+offset)) * (pow(r, (double)(n+offset)) / fact_sqrt(n+offset));
            QDPACK_SET_COMPLEX(&z, dmij * cos((m-n) * theta), dmij * sin((m-n) * theta));
            qdpack_operator_set(dm, m, n, z);
        }
    }

    qdpack_hilbert_space_free(qs);

    return dm;
}


/**
 * \breif    Initial density matrix, Fock state.
 *
 * Generate a density matrix that represent a fock state
 *
 * @param    n    Fock state excitation number.
 * @param    N    Number of fock states in the system.
 *
 */
qdpack_operator_t *
qdpack_dm_fock_state(int n, int N)
{
    qdpack_operator_t *dm;
    
    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    if (n >= N)
    {
                fprintf(stderr, "qdpack_dm_fock_state: n >= N.\n");
                return NULL;
    }
   
    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_fock_state: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }

    qdpack_operator_set_zero(dm);
    qdpack_operator_set(dm, n, n, QDPACK_COMPLEX_ONE);
    
    qdpack_hilbert_space_free(qs);

    return dm;
}

/**
 * \breif    Initial density matrix, uniform superposition
 *
 * Generate a density matrix that represent a uniform superposition of all states
 *
 * @param    N    Number of fock states in the system.
 *
 */
qdpack_operator_t *
qdpack_dm_uniform_superposition(int N)
{
    qdpack_complex z;
    qdpack_operator_t *dm;
    
    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    if (0 >= N)
    {
         fprintf(stderr, "qdpack_dm_fock_state: n >= N.\n");
         return NULL;
    }
    
    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_fock_state: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }

    QDPACK_SET_COMPLEX(&z, 1.0/sqrt(N), 0);
    qdpack_operator_set_all(dm, z);
    
    qdpack_hilbert_space_free(qs);

    return dm;
}

/**
 * \breif    Initial density matrix, uniform maximally mixed state
 *
 * Generate a density matrix that represent a uniform and maximally mixture of all states
 *
 * @param    N    Number of fock states in the system.
 *
 */
qdpack_operator_t *
qdpack_dm_uniform_mixture(int N)
{
    int i;
    qdpack_complex z;
    qdpack_operator_t *dm;    

    qdpack_hilbert_space_t *qs;
    
    if (0 >= N)
    {
         fprintf(stderr, "qdpack_dm_fock_state: n >= N.\n");
         return NULL;
    }

    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);
    
    if ((dm = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_dm_fock_state: Failed to allocate dm (Out of memory?).\n");
        return NULL;
    }

    QDPACK_SET_COMPLEX(&z, 1.0/N, 0);
    qdpack_operator_set_zero(dm);

    for (i = 0; i < dm->m; i++)
    {
        qdpack_operator_set(dm, i, i, z);
    }

    qdpack_hilbert_space_free(qs);
    
    return dm;
}

/*------------------------------------------------------------------------------------------------
 * Density matrix and wave function operations:
 *
 */
int
qdpack_dm_mix(qdpack_operator_t *rho_mix, qdpack_operator_t *rho1, qdpack_operator_t *rho2)
{
    qdpack_complex z;
    QDPACK_SET_COMPLEX(&z, 0.5, 0.0);

    qdpack_operator_memcpy(rho_mix, rho1);
    qdpack_operator_add(rho_mix, rho2);
    qdpack_operator_scale(rho_mix, z);

    return 0;
}

//------------------------------------------------------------------------------
// Wave functions (state vectors)
//------------------------------------------------------------------------------

/**
 * \breif    coherent state wave-function
 *
 * Generate a wave matrix that represents a coherent state:
 *
 * alpha = r exp(i theta)
 * @param       r        Coherent state amplitude.
 * @param       theta    Coherent state phase.
 * @param       N        Number of fock states.
 *
 */

qdpack_state_t *
qdpack_wf_coherent_state(double r, double theta, int N, int offset)
{
    int n;
    qdpack_complex z;
    qdpack_state_t *wf;
    qdpack_hilbert_space_t *qs;
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    if ((wf = qdpack_state_alloc(qs)) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate wf (Out of memory?).\n", __PRETTY_FUNCTION__);
        qdpack_hilbert_space_free(qs);
        return NULL;
    }

    qdpack_state_set_zero(wf);

    if (r == 0)
    {
        qdpack_state_set(wf, 0, QDPACK_COMPLEX_ONE);
    }
    else
    {
        for (n = 0; n < N; n++)
        {    
            double wfi = exp(-r * r / 2.0) * (pow(r, (double)(n+offset)) / fact_sqrt(n+offset));
            QDPACK_SET_COMPLEX(&z, wfi * cos(n * theta), wfi * sin(n * theta));
            qdpack_state_set(wf, n, z);
        }
    }

    qdpack_hilbert_space_free(qs);
    return wf;
}

qdpack_state_t *
qdpack_wf_superposition(qdpack_state_t *wf1, qdpack_state_t *wf2)
{
    qdpack_state_t *wf;
    qdpack_complex z;
    
    QDPACK_SET_COMPLEX(&z, sqrt(0.5), 0.0);
       
    if ((wf = qdpack_state_alloc(wf1->qs)) == NULL)
    {
        fprintf(stderr, "qdpack_wf_fock_state: Failed to allocate wf (Out of memory?).\n");
        wf = NULL;
    }
    else
    {
        qdpack_state_memcpy(wf, wf1);
        qdpack_state_add(wf, wf2);
        qdpack_state_scale(wf, z);
    }

    return wf;
}

/**
 * \breif    Fock state wave function.
 *
 * Generate a wave-function that represent a fock state
 *
 * @param    n    Fock state excitation number.
 * @param    N    Number of fock states in the system.
 *
 */
qdpack_state_t *
qdpack_wf_fock_state(int n, int N)
{
    qdpack_state_t *wf;
    qdpack_hilbert_space_t *qs;
    
    if (n >= N)
    {
        fprintf(stderr, "qdpack_wf_fock_state: n >= N.\n");
        return NULL;
    }
    
    qs = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs, N);

    if ((wf = qdpack_state_alloc(qs)) == NULL)
    {
        fprintf(stderr, "qdpack_wf_fock_state: Failed to allocate wf (Out of memory?).\n");
        wf = NULL;
    }
    else
    {
        qdpack_state_set_zero(wf);
        qdpack_state_set(wf, n, QDPACK_COMPLEX_ONE);
    }
    
    qdpack_hilbert_space_free(qs);
    return wf;
}

int
qdpack_state_to_operator(qdpack_operator_t *rho, qdpack_state_t *wf)
{
    qdpack_complex x, y;
    int i, j;

    if ((rho->m) != (wf->n) || (rho->n) != (wf->n))
    {
        fprintf(stderr, "qdpack_state_to_operator: vector-matrix size mismatch\n");
        return -1;
    }

    for (i = 0; i < wf->n; i++)
    {
        for (j = 0; j < wf->n; j++)
        {
            x = qdpack_state_get(wf, i);
            y = qdpack_state_get(wf, j);
            qdpack_operator_set(rho, i, j, qdpack_complex_mul(x, qdpack_complex_conjugate(y)));
        }
    }

    return 0;
}














