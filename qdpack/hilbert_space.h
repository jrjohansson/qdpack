//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions related to definition of quantum systems (composed of subsystems
 * with different sizes) and enummeration of states
 */

#include <math.h>

#ifndef QUANTUM_SYSTEM_H
#define QUANTUM_SYSTEM_H

#define NQS_MAX 20

typedef int qdpack_hilbert_space_state_vector_t[NQS_MAX];

typedef struct {
    /* system 0: 2 levels
     * system 1: 2 levels
     * system 3: 10 levels
     * etc.
     */

    /* Number of subsystems in the quantum system */
    int nsubsys;

    /* Number of states of each subsystem */
    int *nstates;
} qdpack_hilbert_space_t;

int    ss_idx_major(int N, int a, int b);
int    ss_idx_minor_ms(int N, int J);
int    ss_idx_minor_ls(int N, int J);

qdpack_hilbert_space_t * qdpack_hilbert_space_new();
qdpack_hilbert_space_t * qdpack_hilbert_space_copy(qdpack_hilbert_space_t *qs);
void    qdpack_hilbert_space_print(qdpack_hilbert_space_t *qs);
void    qdpack_hilbert_space_free(qdpack_hilbert_space_t *qs);
void    qdpack_hilbert_space_add(qdpack_hilbert_space_t *qs, int M);
int     qdpack_hilbert_space_nstates(qdpack_hilbert_space_t *qs);
void    qdpack_hilbert_space_state_vector(qdpack_hilbert_space_t *qs, int qsn, qdpack_hilbert_space_state_vector_t qsv);
int     qdpack_hilbert_space_state_number(qdpack_hilbert_space_t *qs, qdpack_hilbert_space_state_vector_t qsv);

qdpack_hilbert_space_t *qdpack_hilbert_subspace_get(qdpack_hilbert_space_t *qs, int n);

#endif
