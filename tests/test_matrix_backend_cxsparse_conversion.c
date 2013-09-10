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

#include <qdpack/qdpack_matrix_cxsparse.h>

int
main(int argc, char **argv)
{
	int i, j, n = 3;

	qdpack_matrix_t *m1;
	qdpack_matrix_t *v1;

	m1 = qdpack_matrix_alloc(n,n);
	v1 = qdpack_matrix_alloc(n,1);

	srand(time(NULL));

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			qdpack_matrix_set(m1, i, j, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
		}
		qdpack_matrix_set(v1, i, 0, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
	}

    //
    // matrix
    //
	
	printf("\nm1 = \n");
	qdpack_matrix_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

	qdpack_matrix_compressed(m1);

	printf("\nm1 = \n");
	qdpack_matrix_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

	qdpack_matrix_triplet(m1);

	printf("\nm1 = \n");
	qdpack_matrix_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

    //
    // vector
    //

	printf("\nv1 = \n");
	qdpack_matrix_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

	qdpack_matrix_compressed(v1);

	printf("\nv1 = \n");
	qdpack_matrix_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

	qdpack_matrix_triplet(m1);

	printf("\nv1 = \n");
	qdpack_matrix_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, v1, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);


	return 0;
}
