//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions related to definition of quantum systems (composed of subsystems
 * with different sizes) and enummeration of states.
 *
 */

#include <qdpack/hilbert_space.h>

#include <stdlib.h>
#include <stdio.h>

/* ----------------------------------------------------------
 * Index conversion functions:
 *
 */

/*
 * Calculate the superoperator index given two operator indices.
 * 
 * J = N*a+b.
 */
int
ss_idx_major(int N, int a, int b)
{
//    printf("idx_major(%d, %d)[%d] = %d\n", a, b, N, N*a+b);
    return N*a + b;
}

/*
 * Calculate the first (most significant: ms) operator index 
 * given one super-operator index (a).
 *
 * J = N*a + b,
 * Solve for a.
 */
int
ss_idx_minor_ms(int N, int J)
{
    return floor(((double)J)/N); // = a
}

/*
 * Calculate the second (least significant: ls) operator index 
 * given one super-operator index (b).
 *
 * J = N*a + b,
 * Solve for b.
 */
int
ss_idx_minor_ls(int N, int J)
{
    return J - N * floor(((double)J)/N); // = b
}

/* ----------------------------------------------------------
 * qdpack_hilbert_space_t object functions:
 *
 */

/** 
 * Create a new qdpack_hilbert_space_t object.
 *
 */
qdpack_hilbert_space_t *
qdpack_hilbert_space_new()
{
    qdpack_hilbert_space_t *qs;

    if ((qs = malloc(sizeof(qdpack_hilbert_space_t))) == NULL)
    {
        fprintf(stderr, "Failed ot allocate a new qdpack_hilbert_space_t object\n");
        return 0;
    }

    qs->nsubsys = 0;
    qs->nstates = NULL;

    return qs;
}

/*
 * Create a new qdpack_hilbert_space_t object using qs as a template.
 *
 *
 */
qdpack_hilbert_space_t *
qdpack_hilbert_space_copy(qdpack_hilbert_space_t *qs)
{
    int i;
    qdpack_hilbert_space_t *qs_new;

    qs_new = qdpack_hilbert_space_new();
    for (i = 0; i < qs->nsubsys; i++)
        qdpack_hilbert_space_add(qs_new, qs->nstates[i]);
    
    return qs_new;
}


/*
 *
 *
 *
 */
void
qdpack_hilbert_space_add(qdpack_hilbert_space_t *qs, int M)
{
    int *nss, i;

    nss = malloc(sizeof(int) * (qs->nsubsys + 1));

    for (i = 0; i < qs->nsubsys; i++)
        nss[i] = qs->nstates[i];
    nss[qs->nsubsys] = M;
    qs->nsubsys = qs->nsubsys + 1;
    free(qs->nstates);
    qs->nstates = nss;
}

/*
 *
 *
 */
void
qdpack_hilbert_space_print(qdpack_hilbert_space_t *qs)
{
    int i;

    printf("Quantum system with %d sub systems: ", qs->nsubsys);
    for (i = 0; i < qs->nsubsys; i++)
        printf("%d ", qs->nstates[i]);
    printf("\n");
}

/*
 *
 *
 */
void
qdpack_hilbert_space_free(qdpack_hilbert_space_t *qs)
{
    free(qs->nstates);
    free(qs);
}



/*
 * Calculate the number of quantum states in the full system.
 */
int
qdpack_hilbert_space_nstates(qdpack_hilbert_space_t *qs)
{
    int i, ns = 1;

    for (i = 0; i < qs->nsubsys; i++)
    {
        ns *= qs->nstates[i];
    }

    return ns;
}

/*
 * Calculate the quantum state vector that correspond to the quantum 
 * state number qsn.
 *
 */
void
qdpack_hilbert_space_state_vector(qdpack_hilbert_space_t *qs, int qsn, qdpack_hilbert_space_state_vector_t qsv)
{
    int i, acc;

    acc = qdpack_hilbert_space_nstates(qs);

    for (i = qs->nsubsys-1; i >= 0; i--)
    {
        acc = acc / qs->nstates[i];
        qsv[i] = floor(qsn/acc);
        qsn -= qsv[i] * acc;
    }
}

/*
 * State vectors
 * 
 */
int
qdpack_hilbert_space_state_vector_print(qdpack_hilbert_space_t *qs, qdpack_hilbert_space_state_vector_t qsv)
{
    int i;

    for (i = 0; i < qs->nsubsys; i++)
        printf("%d ", qsv[i]);

    return 0;
}

/*
 * Calculate the quantum state number that corresponds to the quantum 
 * state vector qsv.
 */
int
qdpack_hilbert_space_state_number(qdpack_hilbert_space_t *qs, qdpack_hilbert_space_state_vector_t qsv)
{
    int qsn, i, acc_nqs;
    
    qsn = 0;
    acc_nqs = 1;

    for (i = 0; i < qs->nsubsys; i++)
    {
        qsn += acc_nqs * qsv[i];
        acc_nqs *= qs->nstates[i];
    }

    return qsn;
}

//
// return a hilbert space that is a subspace of the given space
//
qdpack_hilbert_space_t *
qdpack_hilbert_subspace_get(qdpack_hilbert_space_t *qs, int n)
{
    qdpack_hilbert_space_t *subspace;

    if ((subspace = qdpack_hilbert_space_new()) == NULL)
    {
        // XXX: error message
        return NULL;
    }
    
    qdpack_hilbert_space_add(subspace, qs->nstates[n]);

    return subspace;
}



        
