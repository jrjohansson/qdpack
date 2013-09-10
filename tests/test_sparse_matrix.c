/*
 * Test sparse matrix extention.
 *
 * JRJ, <robert@riken.jp>
 */
#include <stdio.h> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include "misc.h"
#include "gsl_ext_sparse.h"

int
main(int argc, char **argv)
{
	int i, j, n = 10, k, K = 30;
	gsl_matrix_complex *M_dense;
	gsl_sparse_matrix_complex *M_sparse;

	srand(time(NULL));

	M_sparse = gsl_sparse_matrix_complex_alloc(n,n,n);

	/*
	 * Create random matrix:
	 *
	 */
	for (k = 0; k < K; k++)
	{
		gsl_complex z;

		i = rand() % n;
		j = rand() % n;

		z = gsl_sparse_matrix_complex_get(M_sparse, i, j);

		GSL_SET_COMPLEX(&z, 1.0 + GSL_REAL(z), 1.0 + GSL_IMAG(z));

		gsl_sparse_matrix_complex_set(M_sparse, i, j, z);
	} 

	// print sparse in column form
	gsl_sparse_matrix_complex_print_triplet(M_sparse);

	// convert from triplet to column format
	gsl_sparse_matrix_complex_convert_col(M_sparse);
	
	// print sparse in column form
	gsl_sparse_matrix_complex_print_col(M_sparse);
	
	// convert to dense
	M_dense = gsl_sparse_matrix_complex_convert_to_dense(M_sparse);

	printf("real M = \n");
	misc_matrix_complex_print_real(M_dense, n);
	printf("imag M = \n");
	misc_matrix_complex_print_real(M_dense, n);
		
	return 0;
}
