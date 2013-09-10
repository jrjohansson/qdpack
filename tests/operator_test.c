/*
 *
 *
 *
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
	gsl_matrix_complex *sx1, *sx2, *sz1, *sz2;
	quantum_system *qs;
	quantum_system_state_vector_t qsv;
	int qsn, i;

	qs = quantum_system_new();
	
	quantum_system_print(qs);
	quantum_system_add(qs, 2);
	//quantum_system_add(qs, 4);
	quantum_system_print(qs);

	qsn = quantum_system_nstates(qs);


	sx1 = operator_sigma_x(qs, 0);
	printf("sx1 =\n");
	misc_matrix_complex_print(sx1, qsn);
	
	sz1 = operator_sigma_z(qs, 0);
	printf("sz1 =\n");
	misc_matrix_complex_print(sz1, qsn);
	
	sx2 = operator_ho_raising(qs, 1);
	printf("sx2 =\n");
	misc_matrix_complex_print(sx2, qsn);
	
	sz2 = operator_ho_lowering(qs, 1);
	printf("sz1 =\n");
	misc_matrix_complex_print(sz2, qsn);
	

	

	
	return 0;
}
