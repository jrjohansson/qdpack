/*
 * Test the angular momentum operators:
 *
 */
#include "quantum_system.h"
#include "operators.h"
#include "misc.h"

#include <stdio.h> 
#include <gsl/gsl_matrix.h>


int
main(int argc, char **argv)
{
	gsl_matrix_complex *Jz, *Jp, *Jm;
	quantum_system *qs;
	quantum_system_state_vector_t qsv;
	int qsn, n, N = 3;
	double j, jmax = 5.0/2.0;

	for (j = 0; j <= jmax; j += 0.5)
	{
		qs = quantum_system_new();
		quantum_system_add(qs, 2 * j + 1);
		quantum_system_print(qs);

		qsn = quantum_system_nstates(qs);

		printf("======================================================================\n");
		printf("=== Calculating angular momentum operators for a spin = %.3f particle\n", j);

		Jz = operator_Jz(qs, 0, j);
		printf("Jz =\n");
		misc_matrix_complex_print_real(Jz, qsn);

		Jp = operator_J_plus(qs, 0, j);
		printf("Jp =\n");
		misc_matrix_complex_print_real(Jp, qsn);

		Jm = operator_J_minus(qs, 0, j);
		printf("Jm =\n");
		misc_matrix_complex_print_real(Jm, qsn);

	}
	
	return 0;
}
