//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// Description:
//
// QDpack objects for operators and state vectors. An abstraction layer between
// the quantum objects and the underlaying matrix representation, allowing for
// alternative matrix packages to be used (GSL, sparse, etc.)
//

#ifndef _QDPACK_OBJECT
#define _QDPACK_OBJECT

#include <qdpack/qdpack.h>

typedef struct _qdpack_operator_t {

    int n, m;

    qdpack_hilbert_space_t *qs;

    qdpack_matrix_t *data;

} qdpack_operator_t;

typedef struct _qdpack_state {

    // flag for bra or ket

    int n;

    qdpack_hilbert_space_t *qs;

    qdpack_matrix_t *data;

} qdpack_state_t;

// state vectors
qdpack_state_t *qdpack_state_alloc(qdpack_hilbert_space_t *qs);
void            qdpack_state_free(qdpack_state_t *op);

qdpack_complex  qdpack_state_get(qdpack_state_t *op, size_t i);
void            qdpack_state_set(qdpack_state_t *op, size_t i, qdpack_complex value);

void            qdpack_state_memcpy(qdpack_state_t *dst, qdpack_state_t *src);

void            qdpack_state_scale(qdpack_state_t *op, qdpack_complex z);
void            qdpack_state_add(qdpack_state_t *op1, qdpack_state_t *op2);
void            qdpack_state_sub(qdpack_state_t *op1, qdpack_state_t *op2);

void            qdpack_state_set_zero(qdpack_state_t *op);

qdpack_complex  qdpack_state_expectation_value(qdpack_state_t *wf, qdpack_operator_t *op);
int             qdpack_state_to_operator(qdpack_operator_t *rho, qdpack_state_t *wf);

qdpack_complex  qdpack_state_dot(qdpack_state_t *a, qdpack_state_t *b);


// operators
qdpack_operator_t *qdpack_operator_alloc(qdpack_hilbert_space_t *qs);
void               qdpack_operator_free(qdpack_operator_t *op);

qdpack_complex     qdpack_operator_get(qdpack_operator_t *op, size_t i, size_t j);
void               qdpack_operator_set(qdpack_operator_t *op, size_t i, size_t j, qdpack_complex value);
void               qdpack_operator_set_all(qdpack_operator_t *op, qdpack_complex value);
void               qdpack_operator_set_zero(qdpack_operator_t *op);
void               qdpack_operator_set_identity(qdpack_operator_t *op);

void               qdpack_operator_memcpy(qdpack_operator_t *dst, qdpack_operator_t *src);

void               qdpack_operator_scale(qdpack_operator_t *op, qdpack_complex z);
void               qdpack_operator_add(qdpack_operator_t *op1, qdpack_operator_t *op2);
void               qdpack_operator_sub(qdpack_operator_t *op1, qdpack_operator_t *op2);


void qdpack_operator_blas_zgemm(QDPACK_TRANSPOSE_t transA,
                                QDPACK_TRANSPOSE_t transB, 
                                qdpack_complex alpha, 
                                qdpack_operator_t *A, 
                                qdpack_operator_t *B,
                                qdpack_complex beta, 
                                qdpack_operator_t *C);

void qdpack_operator_multiply(qdpack_operator_t *A, qdpack_operator_t *B, qdpack_operator_t *C);
void qdpack_operator_state_multiply(qdpack_operator_t *A, qdpack_state_t *b, qdpack_state_t *c);

qdpack_complex qdpack_operator_trace(qdpack_operator_t *rho);
qdpack_complex qdpack_operator_expectation_value(qdpack_operator_t *rho, qdpack_operator_t *op);

#endif
