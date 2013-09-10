//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/**
 * @file simulation.c
 *
 * @brief Functions for managing the qdpack simulation object.  
 * 
 * @author J Robert Johansson <robert@riken.jp>
 *
 */

#include <stdio.h>

#include <qdpack/qdpack.h>

//
// initialize the default parameter values:
//
void
qdpack_simulation_init(qdpack_simulation_t *sp)
{
    if (sp)
    {
        memset((void *)sp, 0, sizeof(qdpack_simulation_t));

        sp->cont_basis_transform = 0;
        sp->option_me_eb         = 0;
        sp->option_adaptive      = 1;

        sp->H0  = NULL;
        sp->H1  = NULL;
        sp->H_t = NULL;

        sp->H0_func = NULL;
        sp->H1_func = NULL;
        sp->Ht_func = NULL;

        sp->h_td_A = 0.0;
        sp->h_td_w = 0.0;
    }

    return;
}

void
qdpack_simulation_free(qdpack_simulation_t *sp)
{
    if (sp)
    {
        if (sp->H0) qdpack_operator_free(sp->H0);
        if (sp->H1) qdpack_operator_free(sp->H1);
        //if (sp->H_t) qdpack_operator_free(sp->H_t);
    }

    return;
}
