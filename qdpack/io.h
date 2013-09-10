//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * A collection of various functions. Printing, saving to file, etc. 
 */

#ifndef MISC_H
#define MISC_H

#include <qdpack/qdpack.h>

// matrix print
// void qdpack_matrix_print(qdpack_matrix *m);

void qdpack_operator_print(qdpack_operator_t *m);
void qdpack_operator_print_real(qdpack_operator_t *m);
void qdpack_operator_print_imag(qdpack_operator_t *m);
void qdpack_operator_print_abs2(qdpack_operator_t *m);

// vector print
//void qdpack_state_print(qdpack_state_t *m);
//void qdpack_state_print(qdpack_state_t *v);
void qdpack_state_print_real(qdpack_state_t *v);
void qdpack_state_print_imag(qdpack_state_t *v);


//
// row based file io
//
int qdpack_file_row_add(const char *filename, const char *row);

//
// qdpack operator file io
//
int qdpack_operator_write_real(char *filename, qdpack_operator_t *m);
int qdpack_operator_write_imag(char *filename, qdpack_operator_t *m);
    
int qdpack_operator_write(char *filename, qdpack_operator_t *m);
int qdpack_operator_read(char *filename, qdpack_operator_t *m);

int qdpack_operator_write_double(char *filename, qdpack_operator_t *m, int real);

int qdpack_operator_diag_write(char *filename, qdpack_operator_t *m, char *msg);
int qdpack_operator_diag_read(char *filename, qdpack_operator_t *m);
int qdpack_operator_diag_real_write(char *filename, qdpack_operator_t *m, char *msg);
int qdpack_operator_diag_real_read(char *filename, qdpack_operator_t *m);


#endif
