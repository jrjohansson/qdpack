//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <qdpack/qdpack_matrix_gsl.h>

int
main(int argc, char **argv)
{
	int i, j, n = 5;

    qdpack_complex z;

	qdpack_matrix_t *m1, *m2, *m3;
	qdpack_matrix_t *v1, *v2;

	srand(time(NULL));

	m1 = qdpack_matrix_alloc(n,n);
	m2 = qdpack_matrix_alloc(n,n);
	m3 = qdpack_matrix_alloc(n,n);

	v1 = qdpack_matrix_alloc(n,1);
	v2 = qdpack_matrix_alloc(n,1);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			qdpack_matrix_set(m1, i, j, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
			qdpack_matrix_set(m2, i, j, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
		}

		qdpack_matrix_set(v1, i, 0, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
	}
	
	printf("\nm1 = \n");
	qdpack_matrix_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

	printf("\nm2 = \n");
	qdpack_matrix_write(stdout, m2, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

	printf("\nv1 = \n");
	qdpack_matrix_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    //
    // matrix multiplication
    //
    qdpack_matrix_multiply(m1, m2, m3);

	printf("\nm3 = \n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    //
    // matrix-vector multiplication
    //
    qdpack_matrix_multiply(m1, v1, v2);

	printf("\nv2 = \n");
	qdpack_matrix_write(stdout, v2, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);


    //
    // vector-matrix-vector multiplication
    //

    z = qdpack_matrix_multiply_vmv(m1, v1);
    printf("\nv1.T * m1 * v1 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    z = qdpack_matrix_multiply_vmv(m1, v2);
    printf("\nv2.T * m1 * v2 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

	return 0;
}
