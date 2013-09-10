//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * Functions for combining density matrices and for tracing out sub-systems from
 * full-system density matrices.
 */

#ifndef COMPOSITE_H
#define COMPOSITE_H

#include <qdpack/qdpack.h>

typedef struct {
    qdpack_operator_t *dm[NQS_MAX];
    int len;
} qdpack_operator_list_t;

typedef struct {
    qdpack_state_t *wf[NQS_MAX];
    int len;
} qdpack_state_list_t;


// density matrix

void qdpack_operator_list_init(qdpack_operator_list_t *dm_list);
void qdpack_operator_list_append(qdpack_operator_list_t *dm_list, qdpack_operator_t *dm);
void qdpack_operator_list_free_entries(qdpack_operator_list_t *dm_list);
int  qdpack_operator_tensor(qdpack_operator_t *dm, qdpack_hilbert_space_t *qs, qdpack_operator_list_t *dm_list);

qdpack_operator_t *qdpack_operator_traceout(qdpack_hilbert_space_t *qs, qdpack_operator_t *dm, int n);
qdpack_operator_t *qdpack_operator_traceout_multi(qdpack_hilbert_space_t *qs, qdpack_operator_t *dm, qdpack_hilbert_space_t *mask);

// wave function

void qdpack_wf_list_init(qdpack_state_list_t *wf_list);
void qdpack_wf_list_append(qdpack_state_list_t *wf_list, qdpack_state_t *wf);
void qdpack_wf_list_free_entries(qdpack_state_list_t *wf_list);
int  qdpack_wf_tensor(qdpack_state_t *wf, qdpack_hilbert_space_t *qs, qdpack_state_list_t *wf_list);

qdpack_operator_t *qdpack_wf_traceout(qdpack_hilbert_space_t *qs, qdpack_state_t *wf, int n);


#endif

