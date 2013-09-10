//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

// Description:
//
// objects and functions for representing and computing with quantum operators
// and state vectors 
//

#include <qdpack/qdpack.h>

//==============================================================================
// abstraction layer for quantum state vectors (ket and bras)
//==============================================================================

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------

qdpack_state_t *
qdpack_state_alloc(qdpack_hilbert_space_t *qs)
{
    qdpack_state_t *state;
    int n;

    if (qs == NULL)
    {
        return NULL;
    }

    n = qdpack_hilbert_space_nstates(qs);

    if ((state = (qdpack_state_t *)malloc(sizeof(qdpack_state_t))) == NULL)
    {
        return NULL;
    }

    state->n = n;

    state->qs = qdpack_hilbert_space_copy(qs);
   
    state->data = qdpack_matrix_alloc(n, 1);

    return state;
}

void
qdpack_state_free(qdpack_state_t *op)
{
    if (op == NULL)
    {
        return;
    }

    qdpack_hilbert_space_free(op->qs);   

    qdpack_matrix_free(op->data);

    free(op);
}

void
qdpack_state_memcpy(qdpack_state_t *dst, qdpack_state_t *src)
{
    if (src == NULL || dst == NULL)
        return;

    dst->n = src->n;
    if (dst->qs)
        qdpack_hilbert_space_free(dst->qs);
    dst->qs = qdpack_hilbert_space_copy(src->qs);
    
    qdpack_matrix_memcpy(dst->data, src->data);
}

qdpack_complex
qdpack_state_get(qdpack_state_t *op, size_t i)
{
    if (op == NULL)
        return QDPACK_COMPLEX_ZERO;

    return qdpack_matrix_get(op->data, i, 0);
}

void
qdpack_state_set(qdpack_state_t *state, size_t i, qdpack_complex value)
{
    if (state == NULL)
        return;

    qdpack_matrix_set(state->data, i, 0, value);
}

void
qdpack_state_set_zero(qdpack_state_t *op)
{
    if (op == NULL)
        return;

    return qdpack_matrix_set_zero(op->data);
}

void
qdpack_state_scale(qdpack_state_t *op, qdpack_complex z)
{
    if (op == NULL)
        return;

    qdpack_matrix_scale(op->data, z);
}

void
qdpack_state_add(qdpack_state_t *op1, qdpack_state_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    qdpack_matrix_add(op1->data, op2->data);
}

void
qdpack_state_sub(qdpack_state_t *op1, qdpack_state_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    qdpack_matrix_sub(op1->data, op2->data);
}

qdpack_complex
qdpack_state_dot(qdpack_state_t *a, qdpack_state_t *b)
{
    return qdpack_matrix_dot(a->data, b->data);
}


/**
 * \breif    Calculate the expectation value for a state vector
 *
 * Calculate the expectation value for the operator op for the system in the state
 * described by the wave function wf.
 * 
 * @param    wf    The state vector
 * @param    op    Operator for which the expectation value is to be calculated.
 * 
 */
qdpack_complex
qdpack_state_expectation_value(qdpack_state_t *wf, qdpack_operator_t *op)
{
    return qdpack_matrix_multiply_vmv(op->data, wf->data);
}

//==============================================================================
// abstraction layer for quantum operators
//==============================================================================

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
qdpack_operator_t *
qdpack_operator_alloc(qdpack_hilbert_space_t *qs)
{
    qdpack_operator_t *op;
    int n;

    if (qs == NULL)
    {
        return NULL;
    }

    n = qdpack_hilbert_space_nstates(qs);

    if ((op = (qdpack_operator_t *)malloc(sizeof(qdpack_operator_t))) == NULL)
    {
        return NULL;
    }

    op->n = n;
    op->m = n;

    op->qs = qdpack_hilbert_space_copy(qs);
   
    op->data = qdpack_matrix_alloc(n, n);

    return op;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void
qdpack_operator_free(qdpack_operator_t *op)
{
    if (op == NULL)
    {
        return;
    }
   
    qdpack_hilbert_space_free(op->qs);

    qdpack_matrix_free(op->data);

    free(op);
}

void
qdpack_operator_memcpy(qdpack_operator_t *dst, qdpack_operator_t *src)
{
    if (src == NULL || dst == NULL)
        return;

    dst->n = src->n;
    dst->m = src->m;
    if (dst->qs)
        qdpack_hilbert_space_free(dst->qs);
    dst->qs = qdpack_hilbert_space_copy(src->qs);
    
    qdpack_matrix_memcpy(dst->data, src->data);
}


qdpack_complex
qdpack_operator_get(qdpack_operator_t *op, size_t i, size_t j)
{
    if (op == NULL)
        return QDPACK_COMPLEX_ZERO;

    return qdpack_matrix_get(op->data, i, j);
}

void
qdpack_operator_set(qdpack_operator_t *op, size_t i, size_t j, qdpack_complex value)
{
    if (op == NULL)
        return;

    return qdpack_matrix_set(op->data, i, j, value);
}

void
qdpack_operator_set_all(qdpack_operator_t *op, qdpack_complex value)
{
    if (op == NULL)
        return;

    return qdpack_matrix_set_all(op->data, value);
}

void
qdpack_operator_set_zero(qdpack_operator_t *op)
{
    if (op == NULL)
        return;

    return qdpack_matrix_set_zero(op->data);
}

void
qdpack_operator_set_identity(qdpack_operator_t *op)
{
    if (op && op->data)
    {    
        qdpack_matrix_set_identity(op->data);
    }
}

void
qdpack_operator_scale(qdpack_operator_t *op, qdpack_complex z)
{
    if (op == NULL)
        return;

    qdpack_matrix_scale(op->data, z);
}

void
qdpack_operator_add(qdpack_operator_t *op1, qdpack_operator_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    qdpack_matrix_add(op1->data, op2->data);
}

void
qdpack_operator_sub(qdpack_operator_t *op1, qdpack_operator_t *op2)
{
    if (op1 == NULL || op2 == NULL)
    {
        fprintf(stderr, "%s: ERROR: op1 or op2 is NULL\n", __PRETTY_FUNCTION__);
        return;
    }
       
    qdpack_matrix_sub(op1->data, op2->data);
}

//
// BLAS style operator/matrix multiplication: C = alpha op(A) * op(B) + beta C
//
void
qdpack_operator_blas_zgemm(QDPACK_TRANSPOSE_t transA,
                           QDPACK_TRANSPOSE_t transB, 
                           qdpack_complex alpha, 
                           qdpack_operator_t *A, 
                           qdpack_operator_t *B,
                           qdpack_complex beta, 
                           qdpack_operator_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        qdpack_matrix_blas_zgemm(transA, transB, alpha, A->data, B->data, beta, C->data);
    }
    else
    {
        fprintf(stderr, "%s: ERROR: NULL valued operator\n", __PRETTY_FUNCTION__);
    }
}

//
// plain operator/matrix multiplication: C = A * B
//
void
qdpack_operator_multiply(qdpack_operator_t *A, qdpack_operator_t *B, qdpack_operator_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        qdpack_matrix_multiply(A->data, B->data, C->data);
    }
    else
    {
        fprintf(stderr, "%s: ERROR: NULL valued operator\n", __PRETTY_FUNCTION__);
    }
}

//
// plain operator - state multiplication: c = A * b
//
void
qdpack_operator_state_multiply(qdpack_operator_t *A, qdpack_state_t *B, qdpack_state_t *C)
{
    if (A && A->data && B && B->data && C && C->data)
    {
        qdpack_matrix_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, A->data, B->data, QDPACK_COMPLEX_ZERO, C->data);
    }
    else
    {
        fprintf(stderr, "%s: ERROR: NULL valued operator\n", __PRETTY_FUNCTION__);
    }
}


/**
 * \breif    Calculate the expectation value
 *
 * Calculate the expectation value for the operator op for the system in the state
 * described by the density matrix rho.
 * 
 * @param    rho    Full system density matrix.
 * @param    op    Operator for which the expectation value is to be calculated.
 * 
 */
qdpack_complex 
qdpack_operator_expectation_value(qdpack_operator_t *rho, qdpack_operator_t *op)
{
    qdpack_matrix_t *prod;
    qdpack_complex tr;
    int i;

    prod = qdpack_matrix_alloc(rho->m, rho->n);

    qdpack_matrix_multiply(op->data, rho->data, prod);
 
    QDPACK_SET_COMPLEX(&tr, 0.0, 0.0);

    for (i = 0; i < rho->m; i++)
    {
        tr = qdpack_complex_add(tr, qdpack_matrix_get(prod, i, i));
    }

    qdpack_matrix_free(prod);
    
    return tr;
}

//
//
//
qdpack_complex
qdpack_operator_trace(qdpack_operator_t *rho)
{
    int i;
    qdpack_complex tr;

    QDPACK_SET_COMPLEX(&tr, 0.0, 0.0);

    for (i = 0; i < rho->m; i++)
    {
        tr = qdpack_complex_add(tr, qdpack_operator_get(rho, i, i));
    }

    return tr;
}


