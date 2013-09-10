/*
 * Test matrix read and write to file operations.
 *
 *
 * JRJ, <robert@riken.jp>
 */

#include <stdio.h> 
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <qdpack/qdpack_matrix.h>

int
main(int argc, char **argv)
{
	int i, j, n = 10;

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
			QDPACK_SET_COMPLEX(&z, rand() / 128.0, rand() / 128.0);
			qdpack_matrix_set(m1, i, j, z);

			QDPACK_SET_COMPLEX(&z, rand() / 128.0, rand() / 128.0);
			qdpack_matrix_set(m2, i, j, z);
		}

		qdpack_matrix_set(v1, i, 0, qdpack_complex_rect(rand() / 128.0, rand() / 128.0));
	}
	
	printf("m1 = \n");
	qdpack_matrix_write(stdout, m1, 0);

	printf("m2 = \n");
	qdpack_matrix_write(stdout, m2, 0);

	printf("v1 = \n");
	qdpack_matrix_write(stdout, v1, 0);

    //
    // matrix multiplication
    //
    qdpack_matrix_multiply(m1, m2, m3);

	printf("m3 = \n");
	qdpack_matrix_write(stdout, m3, 0);

    //
    // matrix-vector multiplication
    //
    qdpack_matrix_multiply(m1, v1, v2);

	printf("v2 = \n");
	qdpack_matrix_write(stdout, v2, 0);


    //
    // vector-matrix-vector multiplication
    //

    z = qdpack_matrix_multiply_vmv(m1, v1);
    printf("v1.T * m1 * v1 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

    z = qdpack_matrix_multiply_vmv(m1, v2);
    printf("v2.T * m1 * v2 = (%.6e, %.6e)\n", QDPACK_REAL(z), QDPACK_IMAG(z));

	return 0;
}
