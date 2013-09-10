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

#include <qdpack/qdpack.h>

int
main(int argc, char **argv)
{
	//char *filename = "matrix.dat";
	int i, j, n = 10;
	qdpack_matrix_t *orig, *copy;

	srand(time(NULL));

	orig = qdpack_matrix_alloc(n,n);
	copy = qdpack_matrix_alloc(n,n);

	/*
	 * Create random matrix:
	 *
	 */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			qdpack_complex z;
			QDPACK_SET_COMPLEX(&z, rand() / 128.0, rand() / 128.0);
			qdpack_matrix_set(orig, i, j, z);
		}
	}
	
	printf("orig = \n");
	qdpack_matrix_write(stdout, orig, 0);
		
	/*
	 * Write matrix to file:
	 */ 
	//XXX misc_matrix_write(filename, orig, n);

	//XXX misc_matrix_read(filename, copy, n);

	printf("copy = \n");
	qdpack_matrix_write(stdout, copy, 0);

	/*
	 * Compare matrices:
	 *
	 */
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			qdpack_complex zo, zc;

			zo = qdpack_matrix_get(orig, i, j);
			zc = qdpack_matrix_get(copy, i, j);

			if (fabs(QDPACK_REAL(zo) - QDPACK_REAL(zc)) > 0.0001 ||
			    fabs(QDPACK_IMAG(zo) - QDPACK_IMAG(zc)) > 0.0001)
			{
				printf("Error: readback missmatch at location [%d, %d]\n", i, j);
			}
		}
	}

	/*
	 * Write matrix to file:
	 * / 
	misc_matrix_diag_write(filename, orig, n, "");

	misc_matrix_diag_read(filename, copy, n);

	printf("copy = \n");
	misc_matrix_complex_print_real(copy, n);

	/ *
	 * Write matrix to file:
	 * / 
	misc_matrix_diag_real_write(filename, orig, n, "");

	misc_matrix_diag_real_read(filename, copy, n);

	printf("copy = \n");
	misc_matrix_complex_print_real(copy, n);


	return 0;
	
	
	/ *
	 * Compare matrices:
	 *
	 * /
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			gsl_complex zo, zc;

			zo = gsl_matrix_complex_get(orig, i, j);
			zc = gsl_matrix_complex_get(copy, i, j);

			if ( fabs(GSL_REAL(zo) - GSL_REAL(zc)) > 0.0001 ||
			     fabs(GSL_IMAG(zo) - GSL_IMAG(zc)) > 0.0001)
			{
				printf("Error: readback missmatch at location [%d, %d]\n", i, j);
			}
		}
	}

	printf("Test finished.\n");
	*/
	return 0;
}
