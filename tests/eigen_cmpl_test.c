/*
 * Test code for computing eigenvalues of a general (non-symmetric) complex matrix.
 *
 * Example intput:

3
(0.816417425405234,0.144827013835311) (0.324151253793389,0.523759318049997) (0.303578688297421,0.976534355431795)
(0.181262716185302,0.117307899985462) (0.261452937033027,0.0836455137468874) (0.374216732569039,0.696820253506303)
(0.538486659061164,0.929988969117403) (0.865246564149857,0.289147571194917) (0.855454510543495,0.186055930331349)
 
 *
 * Complile:
 
gcc eigen_cmpl_test.c gsl_ext.c -o eigen_cmpl_test  -L/usr/local/lib/lapack/ -llapack -lblas -lm -lF77 -lI77 -I/usr/local/include/lapack -lgsl -lcblas
 
  
 * Robert Johansson, <robert@riken.jp>
 */
#include <stdio.h> 
#include <stdlib.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>

#include "gsl_ext.h"

int
main(int argc, char **argv)
{
	int n, i, j;
	gsl_matrix_complex *A, *evec;
	gsl_vector_complex *eval;

	// Read a complex matrix from stdin and allocate memory
	scanf("%d",&n);
	printf("\nReading complex matrix of size %dx%d from stdin:\n",n, n);
	
	A    = gsl_matrix_complex_alloc(n,n);
	evec = gsl_matrix_complex_alloc(n,n);
	eval = gsl_vector_complex_alloc(n);

	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{	
			gsl_complex z;
			double re, im;
			scanf(" (%lf,%lf)", &re, &im);
			GSL_SET_COMPLEX(&z, re, im);
			gsl_matrix_complex_set(A, i, j, z);
		}
	}

	printf("A = \n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{	
			gsl_complex z;
			z = gsl_matrix_complex_get(A, i, j);
			printf(" (%.15g,%.15g)", GSL_REAL(z), GSL_IMAG(z));
		}
		printf("\n");
	}

	// Start calculating eigen values and eigen vectors
	gsl_ext_eigen_zgeev(A, evec, eval);

	printf("Eigenvalues = \n");
	for (i = 0; i < n; i++)
	{
		gsl_complex z;
		z = gsl_vector_complex_get(eval, i);
		printf("\t(%f, %f)", GSL_REAL(z), GSL_IMAG(z));
	}
	printf("\n");

	printf("evec = \n");
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
		{
			gsl_complex z;
			z = gsl_matrix_complex_get(evec, i, j);
			printf("\t(%f, %f)", GSL_REAL(z), GSL_IMAG(z));
		}
		printf("\n");
	}
	
	/*	
	for (j = 0; j < n; j++)
	{
		printf("Eigenvector%d = \n", j);
		for (i = 0; i < n; i++)
		{
			gsl_complex z;
			z = gsl_matrix_complex_get(evec, i, j);
			printf("\t(%f, %f)", GSL_REAL(z), GSL_IMAG(z));
		}
		printf("\n");
	}
	*/
	
	return 0;
}
