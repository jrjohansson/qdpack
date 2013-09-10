//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------


#include <stdio.h>

#include <gsl/gsl_sf_legendre.h>

#include <qdpack/qdpack.h>

/*
 * Generate basis transformation matrix for a harmonic oscillator with
 * frequency switch: S
 *
 * rho_t = S * rho_t0 * S'
 *
 */
qdpack_operator_t *
basis_transform(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, int n, ho_w_func_t ho_w_cb, double t)
{
    qdpack_complex z;
    qdpack_operator_t *S;
    int i, j, ns, mm, ll;
    double r, w0, wt, sij;
    qdpack_operator_list_t op_list;

    ns = qs->nstates[n];

    if ((S = qdpack_operator_alloc(qs)) == NULL)
    {
        fprintf(stderr, "basis_transform_ho: Failed to allocate basis tranformation matrix (Out of memory?).\n");
 
    }
    qdpack_operator_set_zero(S);
    
    w0 = ho_w_cb(0.0, sim);
    wt = ho_w_cb(t, sim);    

    r = 2 * sqrt(w0 * wt) / (wt + w0);
    if (r > 1.0)
        r = 1.0;
    else if (r < -1.0)
        r = -1.0;
    
    //printf("DEBUG: basis_transform: r = %f (%f, %f): %f\n", r, w0, wt, t);
    
    for (j = 0; j < ns; j++)
    {
        for (i = j; i < ns; i++)
        {
            if ((i == j) || ((i+j)%2 == 0))
            {
                ll = (i+j)/2;
                mm = abs(i-j)/2;
                sij = fact_sqrt(MIN(i,j)) / fact_sqrt(MAX(i,j)) * sqrt(r) * gsl_sf_legendre_Plm(ll, mm, r);
                //if (wt > w0)
                QDPACK_SET_COMPLEX(&z, sij, 0.0);
                //else    
                //    QDPACK_SET_COMPLEX(&z, -sij, 0.0);
                qdpack_operator_set(S, i, j, z);
                
                if (i != j)
                {
                    sij = pow(-1.0, (j-i)/2.0) * sij;
                    QDPACK_SET_COMPLEX(&z, sij, 0.0);
                    qdpack_operator_set(S, j, i, z);
                }
            }
        }
    }

    if (qs->nsubsys > 1)
    {
        qdpack_operator_list_init(&op_list);

        for (i = 0; i < qs->nsubsys; i++)
        {
            qdpack_operator_t *M;
            qdpack_hilbert_space_t *qs_sub = qdpack_hilbert_space_new();
            qdpack_hilbert_space_add(qs_sub, qs->nstates[i]);

            M = qdpack_operator_alloc(qs_sub);

            if (i != n)
            {
                operator_unit(M, qs, i);
            }
            else
            {
                M = S;
            }

            qdpack_operator_list_append(&op_list, M);

            // XXX: free qs, free M
        }
        qdpack_operator_tensor(S, qs, &op_list);
        qdpack_operator_list_free_entries(&op_list);
    }
    
    return S;
}

/*
 * Calculate dS/dt at time t. 
 * 
 * Method: dS/dt = (S(t+h)-S(t-h))/(2h)
 */
qdpack_operator_t *
basis_transform_derivative(qdpack_simulation_t *sim, qdpack_hilbert_space_t *qs, int n, ho_w_func_t ho_w_cb, double t)
{    
    double h;
    qdpack_complex z;
    qdpack_operator_t *S1, *S2;

    h = 1e-6;
    S1 = basis_transform(sim, qs, n, ho_w_cb, t-h);
    S2 = basis_transform(sim, qs, n, ho_w_cb, t+h);

    qdpack_operator_sub(S2, S1);
    
    if (ho_w_cb(t, sim) < ho_w_cb(0.0, sim))
    {
        QDPACK_SET_COMPLEX(&z, -1.0 / (2.0 * h), 0.0);
    }
    else
    {
        QDPACK_SET_COMPLEX(&z, 1.0 / (2.0 * h), 0.0);
    }
        
    qdpack_operator_scale(S2, z);

    qdpack_operator_free(S1);

    return S2;
}    



