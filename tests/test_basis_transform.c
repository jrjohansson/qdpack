/*
 * Test basis transform matrix generation.
 *
 * J Robert Johansson, <robert@riken.jp>
 */

#include <quantum_system.h>
#include <basis_transform.h>
#include <misc.h>
#include <stdio.h> 


#define wm 1.00
#define wp 0.85
#define g  1.0
#define T  10.0

static double
ho_w_func(double t, void *data)
{
	double w = sqrt(wm*wm + (wp*wp-wm*wm) / (1 + exp(-2.0 * g * (t-T))));

	printf("DEBUG: ho_w: w(t=%f) = %f\n", t, w);

	return w;
}

int
main(int argc, char **argv)
{
	gsl_matrix_complex *S, *dS;
	quantum_system *qs;
	int qsn, i;

	qs = quantum_system_new();
	quantum_system_print(qs);
	quantum_system_add(qs, 10);
//	quantum_system_add(qs, 5);
	quantum_system_print(qs);
	qsn = quantum_system_nstates(qs);
	
	S = basis_transform(qs, 0, (ho_w_func_t)ho_w_func, 8.5);  
	printf("S =\n");
	misc_matrix_complex_print_real(S, qsn);
	dS = basis_transform_derivative(qs, 0, (ho_w_func_t)ho_w_func, 8.5);  
	printf("dS =\n");
	misc_matrix_complex_print_real(dS, qsn);
	
	return 0;
}
