//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file operators.c
 *
 * @brief Functions for generating matrix representations of quantum operators.  
 * 
 * @author Robert Johansson <robert@riken.jp>
 *
 */

#include <stdio.h>

#include <qdpack/qdpack.h>


/** ----------------------------------------------------------------------------------
 * Spin 1/2 (qubit) operators:
 *
 */

/**
 * Generate the unit-matrix for the sub-system k.
 *
 */
int 
operator_unit(qdpack_operator_t *u, qdpack_hilbert_space_t *qs, int k)
{
    qdpack_operator_set_identity(u);   
    return 0;
}

/**
 * Generate the sigma_z operator for the n:th subsystem.
 *
 */
int 
operator_sigma_z(qdpack_operator_t *sz, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i;

    if (qs == NULL)
    {
        fprintf(stderr, "%s: NULL-valued qs argument\n", __PRETTY_FUNCTION__);
        return -1;
    }

    if (qs->nstates[n] != 2)
    {
        fprintf(stderr, "%s: The sigma_z operator can only be generated for two-level systems.\n", __PRETTY_FUNCTION__);
        return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);

    qdpack_operator_set_zero(sz);

    for (i = 0; i < ns; i++) 
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv);
        if (qsv[n] == 0)
            qdpack_operator_set(sz, i, i, QDPACK_COMPLEX_ONE);
        else 
            qdpack_operator_set(sz, i, i, QDPACK_COMPLEX_NEGONE);
    }

    return 0;
}


/**
 * Generate the sigma_x operator for the n:th subsystem.
 *
 */
int
operator_sigma_x(qdpack_operator_t *sx, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i, j;

    if (qs == NULL)
    {
        fprintf(stderr, "operator_sigma_x called with NULL-valued qs argument\n");
        return -1;
    }

    if (qs->nstates[n] != 2)
    {
        fprintf(stderr, "The sigma_x operator cannot be generated for systems that are not TLSs.\n");
        return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);

    qdpack_operator_set_zero(sx);

    for (i = 0; i < ns; i++) 
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv);

        qsv[n] = (qsv[n] == 1) ? 0 : 1; 

        j = qdpack_hilbert_space_state_number(qs, qsv);
            
        qdpack_operator_set(sx, i, j, QDPACK_COMPLEX_ONE);
    }

    return 0;
}

/**
 * Generate the sigma_y operator for the n:th subsystem.
 *
 */
int 
operator_sigma_y(qdpack_operator_t *sy, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i, j;
    qdpack_complex z;

    if (qs == NULL)
    {
        fprintf(stderr, "operator_sigma_y called with NULL-valued qs argument\n");
        return -1;
    }

    if (qs->nstates[n] != 2)
    {
        fprintf(stderr, "The sigma_y operator cannot be generated for systems that are not TLSs.\n");
        return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);

    sy = qdpack_operator_alloc(qs);

    for (i = 0; i < ns; i++) 
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv);

        if (qsv[n] == 1)
        {
            qsv[n] = 0; 
            QDPACK_SET_COMPLEX(&z, 0.0, 1.0);
        }
        else
        {
            qsv[n] = 1; 
            QDPACK_SET_COMPLEX(&z, 0.0, -1.0);
        }
        
        j = qdpack_hilbert_space_state_number(qs, qsv);
            
        qdpack_operator_set(sy, i, j, z);
    }

    return 0;
}

int
operator_sigma_plus(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n)
{
    // equivalent to 
    // operator_project(param.H0, qs, n, 0, 1); 
    return operator_ho_raising(a, qs, n, 0);
}

int
operator_sigma_minus(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n)
{
    // equivalent to 
    // operator_project(param.H0, qs, n, 1, 0); 
    return operator_ho_lowering(a, qs, n, 0);
}


// ------------------------------------------------------------------------------
// Harmonic oscillator operators:
//
//

/**
 * Generate the lowering operator for an harmonic oscillator that is
 * sub-system n: a
 *
 */
int
operator_ho_lowering(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n, int offset)
{
    qdpack_operator_list_t dm_list;
    int i, j;    

    qdpack_operator_list_init(&dm_list);
    
    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_operator_t *M;
        qdpack_hilbert_space_t *qs_sub;

        qs_sub = qdpack_hilbert_subspace_get(qs, i);

        M = qdpack_operator_alloc(qs_sub);

        if (i != n)
        {
            operator_unit(M, qs_sub, 0);
        }
        else
        {
            qdpack_operator_set_zero(M);
            for (j = 1; j < qs->nstates[i]; j++)
            {
                qdpack_complex z;
                QDPACK_SET_COMPLEX(&z, sqrt(j+offset), 0);
                qdpack_operator_set(M, j-1, j, z);
            }
        }

        qdpack_operator_list_append(&dm_list, M);
        qdpack_hilbert_space_free(qs_sub);
    }
    qdpack_operator_tensor(a, qs, &dm_list);
    qdpack_operator_list_free_entries(&dm_list);
    
    return 0; 
}



/**
 * Generate the lowering operator for an harmonic oscillator that is
 * sub-system n: a dagger
 *
 */
int
operator_ho_raising(qdpack_operator_t *a, qdpack_hilbert_space_t *qs, int n, int offset)
{
    qdpack_operator_list_t dm_list;
    int i, j;    

    qdpack_operator_list_init(&dm_list);
    
    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_operator_t *M;
        qdpack_hilbert_space_t *qs_sub;

        qs_sub = qdpack_hilbert_subspace_get(qs, i);

        M = qdpack_operator_alloc(qs_sub);

        if (i != n)
        {
            operator_unit(M, qs_sub, 0);
        }
        else
        {
            qdpack_operator_set_zero(M);
            for (j = 1; j < qs->nstates[i]; j++)
            {
                qdpack_complex z;
                QDPACK_SET_COMPLEX(&z, sqrt(j+offset), 0);
                qdpack_operator_set(M, j, j-1, z);
            }
        }

        qdpack_operator_list_append(&dm_list, M);
        qdpack_hilbert_space_free(qs_sub);
    }
    qdpack_operator_tensor(a, qs, &dm_list);
    qdpack_operator_list_free_entries(&dm_list);
    
    return 0; 
}

/**
 * Calculate the N = ad a operator for the nth sub-system.
 *
 */
int
operator_ho_N(qdpack_operator_t *n_op, qdpack_hilbert_space_t *qs, int n, int offset)
{
    qdpack_operator_list_t dm_list;
    int i, j;

    qdpack_operator_list_init(&dm_list);

    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_operator_t *M;
        qdpack_hilbert_space_t *qs_sub;

        qs_sub = qdpack_hilbert_subspace_get(qs, i);

        M = qdpack_operator_alloc(qs_sub);


        if (i != n)
        {
            operator_unit(M, qs_sub, 0);
        }
        else
        {
            qdpack_operator_set_zero(M);
            for (j = 0; j < qs->nstates[i]; j++)
            {
                qdpack_complex z;
                QDPACK_SET_COMPLEX(&z, j+offset, 0);
                qdpack_operator_set(M, j, j, z);
            }
        }

        qdpack_operator_list_append(&dm_list, M);
        qdpack_hilbert_space_free(qs_sub);
    }
    qdpack_operator_tensor(n_op, qs, &dm_list);
    qdpack_operator_list_free_entries(&dm_list);

    return 0;
}

/**
 * Generic creation/annihilation operators:
 *
 */
int
operator_project(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n, int si, int sj)
{
    qdpack_operator_list_t dm_list;
    int i;

    qdpack_operator_list_init(&dm_list);

    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_operator_t *M;
        qdpack_hilbert_space_t *qs_sub;

        qs_sub = qdpack_hilbert_subspace_get(qs, i);

        M = qdpack_operator_alloc(qs_sub);

        if (i != n)
        {
            operator_unit(M, qs_sub, 0);
        }
        else
        {
            qdpack_operator_set_zero(M);
            qdpack_operator_set(M, sj, si, QDPACK_COMPLEX_ONE); // XXX: swapped indexes?
        }

        qdpack_operator_list_append(&dm_list, M);
        qdpack_hilbert_space_free(qs_sub);
    }
    qdpack_operator_tensor(op, qs, &dm_list);
    qdpack_operator_list_free_entries(&dm_list);

    return 0;
}

// -----------------------------------------------------------------------------
// 3-level system operators:
//
//
//

/**
 * Generate the 3LS sigma_z operator for the n:th subsystem.
 *
 */
int
operator_3ls_sigma_z(qdpack_operator_t *sz, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i;
    
    if (qs == NULL)
    {
        fprintf(stderr, "operator_sigma_z called with NULL-valued qs argument\n");
        return -1;
    }

    if (qs->nstates[n] != 3)
    {
        fprintf(stderr, "This sigma_z operator cannot be generated for systems that are not 3LS.\n");
        return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);

    sz = qdpack_operator_alloc(qs);

    for (i = 0; i < ns; i++) 
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv);
        if (qsv[n] == 0)
            qdpack_operator_set(sz, i, i, QDPACK_COMPLEX_ONE);
        else 
            qdpack_operator_set(sz, i, i, QDPACK_COMPLEX_NEGONE);
    }

    return 0;
}


/**
 * Generate the 3LS sigma_x operator for the n:th subsystem.
 *
 */
int
operator_3ls_sigma_x(qdpack_operator_t *sx, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i, j;

    if (qs == NULL)
    {
        fprintf(stderr, "operator_sigma_x called with NULL-valued qs argument\n");
        return -1;
    }

    if (qs->nstates[n] != 3)
    {
        fprintf(stderr, "This sigma_x operator cannot be generated for systems that are a 3LS.\n");
        return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);

    qdpack_operator_set_zero(sx);

    for (i = 0; i < ns; i++) 
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv);

        if (qsv[n] != 2)
        {
            if (qsv[n] == 1)
                qsv[n] = 0; 
            else 
                qsv[n] = 1; 
        
            j = qdpack_hilbert_space_state_number(qs, qsv);
            
            qdpack_operator_set(sx, i, j, QDPACK_COMPLEX_ONE);
        }
    }

    return 0;
}


// -----------------------------------------------------------------------------
// Angular momentum operators:
// -----------------------------------------------------------------------------

/**
 * Generate the sigma_z operator for the n:th subsystem.
 *
 */
int
operator_Jz(qdpack_operator_t *Jz, qdpack_hilbert_space_t *qs, int n, double j)
{
    qdpack_hilbert_space_state_vector_t qsv;
    int ns, i;
    double m;
    int N = (int)(2 * j + 1);

    if (qs == NULL)
    {
            fprintf(stderr, "operator_Jz called with NULL-valued qs argument\n");
            return -1;
    }

    if (qs->nstates[n] != N)
    {
            fprintf(stderr, "The sigma_Jz operator needs %d quantum states (given %d states).\n", N, qs->nstates[n]);
            return -1;
    }

    ns = qdpack_hilbert_space_nstates(qs);
    qdpack_operator_set_zero(Jz);

    for (i = 0; i < ns; i++)
    {
        qdpack_complex z;

        qdpack_hilbert_space_state_vector(qs, i, qsv);

        m = j - qsv[n];        

        QDPACK_SET_COMPLEX(&z, m, 0.0);
        qdpack_operator_set(Jz, i, i, z);
    }

    return 0;
}

/**
 * Generate matrix representation of the J+ operator
 * 
 */
int
operator_J_plus(qdpack_operator_t *Jp, qdpack_hilbert_space_t *qs, int n, double jj)
{
    qdpack_operator_list_t dm_list;
    int i, j, N = (int)(2 * jj + 1);
    double m;

    if (qs->nstates[n] != N)
    {
            fprintf(stderr, "The sigma_Jp operator needs %d quantum states (given %d states).\n", N, qs->nstates[n]);
            return -1;
    }

    qdpack_operator_list_init(&dm_list);

    for (i = 0; i < qs->nsubsys; i++)
    {
        qdpack_operator_t *M;
        qdpack_hilbert_space_t *qs_sub = qdpack_hilbert_space_new();
        qdpack_hilbert_space_add(qs_sub, qs->nstates[n]);

        M = qdpack_operator_alloc(qs_sub);

        if (i != n)
        {
             operator_unit(M, qs, i);
        }
        else
        {
           qdpack_operator_set_zero(M);

           for (j = 1; j < qs->nstates[n]; j++)
           {
               qdpack_complex z;

               m = jj - j;
               QDPACK_SET_COMPLEX(&z, sqrt(jj*(jj+1) - m*(m+1)), 0);

                qdpack_operator_set(M, j-1, j, z);
            }

            // XXX free qs_sub
        }

        qdpack_operator_list_append(&dm_list, M);
    }

    qdpack_operator_tensor(Jp, qs, &dm_list);
    qdpack_operator_list_free_entries(&dm_list);

    return 0;
}

/**
 * Generate matrix representation of the J- operator
 *
 */
int
operator_J_minus(qdpack_operator_t *Jm, qdpack_hilbert_space_t *qs, int n, double jj)
{
        qdpack_operator_list_t dm_list;
        int i, j, N = (int)(2 * jj + 1);
        double m;

        if (qs->nstates[n] != N)
        {
                fprintf(stderr, "The sigma_Jp operator needs %d quantum states (given %d states).\n", N, qs->nstates[n]);
                return -1;
        }

        qdpack_operator_list_init(&dm_list);

        for (i = 0; i < qs->nsubsys; i++)
        {
            qdpack_operator_t *M;
            qdpack_hilbert_space_t *qs_sub = qdpack_hilbert_space_new();
            qdpack_hilbert_space_add(qs_sub, qs->nstates[n]);

            M = qdpack_operator_alloc(qs_sub);

            if (i != n)
            {
                operator_unit(M, qs, i);
            }
            else
            {
                qdpack_operator_set_zero(M);

                for (j = 0; j < qs->nstates[n]-1; j++)
                {
                    qdpack_complex z;
                    m = jj - j;

                    QDPACK_SET_COMPLEX(&z, sqrt(jj*(jj+1) - m*(m-1)), 0);
                    //QDPACK_SET_COMPLEX(&z, m, 0);
                    qdpack_operator_set(M, j+1, j, z);
                 }
             }
             qdpack_operator_list_append(&dm_list, M);

            // XXX free qs_sub
        }
        qdpack_operator_tensor(Jm, qs, &dm_list);
        qdpack_operator_list_free_entries(&dm_list);

        return 0;
}


/* -----------------------------------------------------------------------------
 * Double quantum dot operators
 *
 */

int
operator_dqd_sigma_z(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_operator_t *op1, *op2;
    
    op1 = qdpack_operator_alloc(qs);
    op2 = qdpack_operator_alloc(qs);

    operator_project(op1, qs, n, 1, 1);
    operator_project(op2, qs, n, 2, 2);

    qdpack_operator_set_zero(op);

    qdpack_operator_scale(op2, QDPACK_COMPLEX_NEGONE);
    qdpack_operator_add(op, op1);
    qdpack_operator_add(op, op2);
    
    qdpack_operator_free(op1);
    qdpack_operator_free(op2);    
    
    return 0;
}

int
operator_dqd_sigma_y(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_complex z;
    qdpack_operator_t *op1, *op2;
    
    op1 = qdpack_operator_alloc(qs);
    op2 = qdpack_operator_alloc(qs);

    operator_project(op1, qs, n, 1, 2);
    operator_project(op2, qs, n, 2, 1);

    qdpack_operator_set_zero(op);

    qdpack_operator_scale(op2, QDPACK_COMPLEX_NEGONE);
    qdpack_operator_add(op, op1);
    qdpack_operator_add(op, op2);
    
    QDPACK_SET_COMPLEX(&z, 0.0, 1.0);
    qdpack_operator_scale(op, z);    
    
    qdpack_operator_free(op1);
    qdpack_operator_free(op2);    
    
    return 0;
}

int
operator_dqd_sigma_x(qdpack_operator_t *op, qdpack_hilbert_space_t *qs, int n)
{
    qdpack_operator_t *op1, *op2;

    op1 = qdpack_operator_alloc(qs);
    op2 = qdpack_operator_alloc(qs);

    operator_project(op1, qs, n, 1, 2);
    operator_project(op2, qs, n, 2, 1);

    qdpack_operator_set_zero(op);

    qdpack_operator_add(op, op1);
    qdpack_operator_add(op, op2);
    
    qdpack_operator_free(op1);
    qdpack_operator_free(op2);    
    
    return 0;

}




/* -----------------------------------------------------------------------------
 * Computational routines:
 *
 */

/**
 * Calculate the trace of a matrix (operator). Used for calculating expectation
 * values and two-time correlation functions.
 */
qdpack_complex
operator_trace(qdpack_operator_t *op) // XXX duplicate?
{
        qdpack_complex tr;
        int i;

        QDPACK_SET_COMPLEX(&tr, 0.0, 0.0);
        for (i = 0; i < op->m; i++)
                tr = qdpack_complex_add(tr, qdpack_operator_get(op, i, i));

        return tr;
}



