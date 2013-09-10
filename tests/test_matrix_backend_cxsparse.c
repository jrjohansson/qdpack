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
	int i, j, n = 4, k, K = 8;

    qdpack_complex z;

	qdpack_matrix_t *m1, *m2, *m3, *m4;
	qdpack_matrix_t *v1, *v2;

	srand(time(NULL));

	m1 = qdpack_matrix_alloc(n,n);
	m2 = qdpack_matrix_alloc(n,n);
	m3 = qdpack_matrix_alloc(n,n);
	m4 = qdpack_matrix_alloc(n,n);

	v1 = qdpack_matrix_alloc(n,1);
	v2 = qdpack_matrix_alloc(n,1);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			QDPACK_SET_COMPLEX(&z, rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX));
			qdpack_matrix_set(m1, i, j, z);
		}

		qdpack_matrix_set(v1, i, 0, qdpack_complex_rect(rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX)));
	}


	for (k = 0; k < K; k++)
	{
		qdpack_complex z;

		i = rand() % n;
		j = rand() % n;

	    QDPACK_SET_COMPLEX(&z, rand() / (0.1 * RAND_MAX), rand() / (0.1 * RAND_MAX));
	    qdpack_matrix_set(m2, i, j, z);
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
    printf("\n\n========== matrix-matrix multiplication ==========\n\n");

    qdpack_matrix_multiply(m1, m2, m3);

	printf("\nm3 = m1 * m2\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    //
    // matrix-vector multiplication
    //
    printf("\n\n========== matrix-vector multiplication ==========\n\n");

    qdpack_matrix_multiply(m1, v1, v2);

	printf("\nv2 = m1 * v1\n");
	qdpack_matrix_write(stdout, v2, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    //
    // vector-matrix-vector multiplication
    //
    printf("\n\n========== vector-matrix-vector multiplication ==========\n\n");

    z = qdpack_matrix_multiply_vmv(m1, v1);
    printf("\nv1.H * m1 * v1 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    z = qdpack_matrix_multiply_vmv(m2, v1);
    printf("\nv1.H * m2 * v1 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));


    //
    // vector-vector multiplication
    //
    printf("\n\n========== vector-vector multiplication ==========\n\n");

    z = qdpack_matrix_dot(v1, v1);
    printf("\nv1.H * v1 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    z = qdpack_matrix_dot(v2, v2);
    printf("\nv2.H * v2 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    z = qdpack_matrix_dot(v1, v2);
    printf("\nv1.H * v2 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));


    //
    // matrix/vector set/get
    //
    printf("\n\n========== matrix/vector set/get ==========\n\n");

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_set_all(m3, QDPACK_COMPLEX_ONE);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_set_zero(m3);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    //
    // matrix/vector copy
    //
    printf("\n\n========== matrix/vector copy ==========\n\n");

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

	printf("\nm2 =\n");
	qdpack_matrix_write(stdout, m2, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_memcpy(m3, m2);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);


    //
    // matrix/vector substraction
    //
    printf("\n\n========== matrix/vector sub ==========\n\n");


    qdpack_matrix_sub(m3, m2);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);
	qdpack_matrix_sparse_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);
    qdpack_matrix_triplet(m3);
	qdpack_matrix_sparse_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_COMPLEX);

    //
    // matrix/vector add
    //
    printf("\n\n========== matrix/vector add ==========\n\n");


    qdpack_matrix_add(m3, m2);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_add(m3, m2);

	printf("\nm3 =\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);


    //
    // matrix-matrix blas multiplication
    //
    printf("\n\n========== matrix-matrix blas multiplication ==========\n\n");

	printf("\nm1 = \n");
	qdpack_matrix_write(stdout, m1, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

	printf("\nm2 = \n");
	qdpack_matrix_write(stdout, m2, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, m1, m2, QDPACK_COMPLEX_ZERO, m3); 

	printf("\nm3 = m1 * m2\n");
	qdpack_matrix_write(stdout, m3, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_multiply(m1, m2, m4);

	printf("\nm4 = m1 * m2\n");
	qdpack_matrix_write(stdout, m4, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_sub(m4, m3);

	printf("\nm4 - m3 =\n");
	qdpack_matrix_write(stdout, m4, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_multiply(m1, m2, m4);

    qdpack_matrix_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, m1, m2, -QDPACK_COMPLEX_ONE, m4); 

	printf("\nm1 * m2 - m1 * m2\n");
	qdpack_matrix_write(stdout, m4, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

    qdpack_matrix_multiply(m1, m2, m4);

    qdpack_matrix_blas_zgemm(QDpackNoTrans, QDpackNoTrans, QDPACK_COMPLEX_ONE, m1, m2, -2*QDPACK_COMPLEX_ONE, m4); 

	printf("\nm1 * m2 - 2 * (m1 * m2)\n");
	qdpack_matrix_write(stdout, m4, QDPACK_MATRIX_WRITE_FORMAT_PYTHON);

	return 0;
}
