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
// density matrices
//------------------------------------------------------------------------------

/**
 * \breif Tensor productor of density matrices
 *
 * Combine a list of density matrices for the sub-systems into
 * a density matrix for the full system (product state), by taking
 * the tensor product.
 *
 * @param    qs        Data structure that defines the quantum system
 * @param    dm_list        A list of density matrices that are to be tensor multiplied
 *
 */
int 
qdpack_operator_tensor(qdpack_operator_t *dm_full, qdpack_hilbert_space_t *qs, qdpack_operator_list_t *dm_list)
{
    int ns, i, j, k;
    qdpack_complex dm_ij;
    qdpack_hilbert_space_state_vector_t qsv1, qsv2;

    ns = qdpack_hilbert_space_nstates(qs);

    if (qs == NULL || dm_list == NULL || qs->nsubsys != dm_list->len)
    {
        fprintf(stderr, "qdpack_operator_tensor: NULL arguments or system-size/dm_list mismatch.\n");
        return -1;
    }

    for (i = 0; i < ns; i++)
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv1);
        
        for (j = 0; j < ns; j++)
        {
            qdpack_hilbert_space_state_vector(qs, j, qsv2);
            
            QDPACK_SET_COMPLEX(&dm_ij, 1.0, 0.0);
            
            for (k = 0; k < qs->nsubsys; k++)
            {
                dm_ij = qdpack_complex_mul(dm_ij, qdpack_operator_get(dm_list->dm[k],qsv1[k],qsv2[k]));
            }
            qdpack_operator_set(dm_full, i, j, dm_ij);
        }
    }
    
    return 0;
}


/**
 * Initialize the data structure that describe a list of 
 * density matrices.
 *
 * @param    dm_list        List of density matrices to initiate.
 */
void
qdpack_operator_list_init(qdpack_operator_list_t *dm_list)
{
    int i;
    
    if (dm_list != NULL)
    {
        for (i = 0; i < NQS_MAX; i++)
            dm_list->dm[i] = NULL;
        dm_list->len = 0;
    }
}

/**
 * Append a density matrix to list of density matrices that describe 
 * sub-systems of the full quantum system.
 *
 * @param    dm_list        List of density matrices.
 * @param    dm        Density matrix of subsystem to add to the list.
 */
void
qdpack_operator_list_append(qdpack_operator_list_t *dm_list, qdpack_operator_t *dm)
{
    if (dm_list->len < NQS_MAX)
        dm_list->dm[dm_list->len++] = dm;
    else
        fprintf(stderr, "qdpack_operator_list_append: Reached max size.\n");
}

/**
 * Free density matrices in the list. 
 * 
 * @param    dm_list        List of density matrices.
 */
void
qdpack_operator_list_free_entries(qdpack_operator_list_t *dm_list)
{
    int i;

    for (i = 0; i < dm_list->len; i++)
        qdpack_operator_free(dm_list->dm[i]);
        
//        free(dm_list->dm[i]);
}

/**
 * \breif    Trace out all sussystem except the nth system
 *
 * Trace out all sub-system except sub-system n from the full-system
 * density matrix dm.
 *
 * @param    qs              Data structure that defines the quantum system
 * @param    dm        Full system density matrix
 * @param    n        The subsystem that is to be kept (all other subsystem traced out)
 *
 */
qdpack_operator_t *
qdpack_operator_traceout(qdpack_hilbert_space_t *qs, qdpack_operator_t *dm, int n)
{
    int i, j, ii, jj, k, ns_part, ns_sub;
    qdpack_operator_t *dm_part;
    qdpack_complex dm_ij;
    qdpack_hilbert_space_state_vector_t qsv1, qsv2;
    qdpack_hilbert_space_t *qs_sub, *qs_part;    

    qs_part = qdpack_hilbert_space_new();
    qdpack_hilbert_space_add(qs_part, qs->nstates[n]);

    ns_part = qs->nstates[n];
    if ((dm_part = qdpack_operator_alloc(qs_part)) == NULL)
    {
        fprintf(stderr, "qdpack_operator_traceout: Failed to allocate dm_part (Out of memory?).\n");
        return NULL;
    }

    qs_sub = qdpack_hilbert_space_copy(qs);
    qs_sub->nstates[n] = 1;
    ns_sub = qdpack_hilbert_space_nstates(qs_sub);
    
    for (i = 0; i < ns_part; i++)
    {
        for (j = 0; j < ns_part; j++)
        {
            QDPACK_SET_COMPLEX(&dm_ij, 0.0, 0.0);    
            for (k = 0; k < ns_sub; k++)
            {
                qdpack_hilbert_space_state_vector(qs_sub, k, qsv2);
                qdpack_hilbert_space_state_vector(qs_sub, k, qsv1);
            
                qsv1[n] = i;
                qsv2[n] = j;

                ii = qdpack_hilbert_space_state_number(qs, qsv1);
                jj = qdpack_hilbert_space_state_number(qs, qsv2);

                dm_ij = qdpack_complex_add(dm_ij, qdpack_operator_get(dm, ii, jj));
            }

            qdpack_operator_set(dm_part, i, j, dm_ij);
        }
    }
        
    return dm_part;
}

/**
 * \breif    Trace out a few of the subsystems from the full density matrix.
 *
 * Trace out subsystems marked in mask variable from the full-system density
 * matrix dm. 
 *
 * mask[n] = 1: -> trace out subsystem n
 * mask[n] = 0: -> keep subsystem n
 *
 * qs_sub:   part of system to be kept
 * qs_trace: part of system to be traced
 *
 * @param       qs        Data structure that defines the quantum system
 * @param       dm        Full system density matrix
 * @param       mask      Mask of subsystems to keep/trace out, see above.
 *
 */
qdpack_operator_t *
qdpack_operator_traceout_multi(qdpack_hilbert_space_t *qs, qdpack_operator_t *dm, qdpack_hilbert_space_t *mask)
{
    int i, j, ii, jj, k, ns_trace, ns_sub, n, m;
    qdpack_operator_t *dm_sub;
    qdpack_complex dm_ij;
    qdpack_hilbert_space_state_vector_t qsv_trace_i, qsv_trace_j, qsv_sub_i, qsv_sub_j;
    qdpack_hilbert_space_t *qs_sub, *qs_trace;    

    qs_trace = qdpack_hilbert_space_copy(qs);
    qs_sub   = qdpack_hilbert_space_new();
    for (n = 0; n < mask->nsubsys; n++)
    {
        if (mask->nstates[n] == 0)
        {    
            // If subsystem is masked "0", add it to qs_sub and remove from qs_trace
            qdpack_hilbert_space_add(qs_sub, qs->nstates[n]);
            qs_trace->nstates[n] = 1;
        }
    }
    ns_sub   = qdpack_hilbert_space_nstates(qs_sub);
    ns_trace = qdpack_hilbert_space_nstates(qs_trace);

    if ((dm_sub = qdpack_operator_alloc(qs_sub)) == NULL)
    {
        fprintf(stderr, "qdpack_operator_traceout: Failed to allocate dm_part (Out of memory?).\n");
        return NULL;
    }
    
    for (i = 0; i < ns_sub; i++)
    {
        qdpack_hilbert_space_state_vector(qs_sub, i, qsv_sub_i);
        for (j = 0; j < ns_sub; j++)
        {
            qdpack_hilbert_space_state_vector(qs_sub, j, qsv_sub_j);
            
            QDPACK_SET_COMPLEX(&dm_ij, 0.0, 0.0);    
            for (k = 0; k < ns_trace; k++)
            {
                qdpack_hilbert_space_state_vector(qs_trace, k, qsv_trace_i);
                qdpack_hilbert_space_state_vector(qs_trace, k, qsv_trace_j);
            
                m = 0;
                for (n = 0; n < mask->nsubsys; n++)
                {
                    if (mask->nstates[n] == 0)
                    {
                        qsv_trace_i[n] = qsv_sub_i[m];
                        qsv_trace_j[n] = qsv_sub_j[m];
                        m++;
                    }                
                }
                
                ii = qdpack_hilbert_space_state_number(qs, qsv_trace_i);
                jj = qdpack_hilbert_space_state_number(qs, qsv_trace_j);
            
                dm_ij = qdpack_complex_add(dm_ij, qdpack_operator_get(dm, ii, jj));
            }

            qdpack_operator_set(dm_sub, i, j, dm_ij);
        }
    }
        
    return dm_sub;
}



//------------------------------------------------------------------------------
// state vectors
//------------------------------------------------------------------------------

/**
 *
 */
int 
qdpack_wf_tensor(qdpack_state_t *wf_full, qdpack_hilbert_space_t *qs, qdpack_state_list_t *wf_list)
{
    int ns, i, k;
    qdpack_complex wf_i;
    qdpack_hilbert_space_state_vector_t qsv1;

    ns = qdpack_hilbert_space_nstates(qs);

    if (qs == NULL || wf_list == NULL || qs->nsubsys != wf_list->len)
    {
        fprintf(stderr, "qdpack_wf_tensor: NULL arguments or system-size/wf_list mismatch.\n");
        return -1;
    }

    for (i = 0; i < ns; i++)
    {
        qdpack_hilbert_space_state_vector(qs, i, qsv1);
        
        QDPACK_SET_COMPLEX(&wf_i, 1.0, 0.0);
        for (k = 0; k < qs->nsubsys; k++)
        {
            wf_i = qdpack_complex_mul(wf_i, qdpack_state_get(wf_list->wf[k], qsv1[k]));
        }
        qdpack_state_set(wf_full, i, wf_i);
    }
    
    return 0;
}

/**
 */
void
qdpack_wf_list_init(qdpack_state_list_t *wf_list)
{
    int i;
    
    if (wf_list != NULL)
    {
        for (i = 0; i < NQS_MAX; i++)
            wf_list->wf[i] = NULL;
        wf_list->len = 0;
    }
}

/**
 */
void
qdpack_wf_list_append(qdpack_state_list_t *wf_list, qdpack_state_t *wf)
{
    if (wf_list->len < NQS_MAX)
        wf_list->wf[wf_list->len++] = wf;
    else
        fprintf(stderr, "qdpack_wf_list_append: Reached max size.\n");
}

/**
 * Free density matrices in the list. 
 * 
 * @param    dm_list        List of density matrices.
 */
void
qdpack_wf_list_free_entries(qdpack_state_list_t *wf_list)
{
    int i;

    for (i = 0; i < wf_list->len; i++)
        qdpack_state_free(wf_list->wf[i]);
        
}


// -----------------------------------------------------------------------------
// Trace out a quantum system from the wave function, results in a density
// matrix. Can be optimized.
//
qdpack_operator_t *
qdpack_wf_traceout(qdpack_hilbert_space_t *qs, qdpack_state_t *wf, int n)
{
    qdpack_operator_t *dm, *dm_part;

    dm = qdpack_operator_alloc(qs);

    qdpack_state_to_operator(dm, wf);

    dm_part = qdpack_operator_traceout(qs, dm, n);

    qdpack_operator_free(dm);

    return dm_part;
}



