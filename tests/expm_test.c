/*
 *
 *
 *
 *
 */
#include "quantum_system.h"
#include "density_matrix.h"
#include "operators.h"
#include "misc.h"
#include "math.h"

#include <stdio.h> 
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>


int
main(int argc, char **argv)
{
	gsl_matrix *sx1, *sx2, *sz1, *sz2, *u, *rho1, *rho2;
	gsl_matrix_complex *dm_full, *dm_part, *expm_test;
	quantum_system *qs;
	quantum_system_state_vector_t qsv;
	density_matrix_list_t dm_list;
	int qsn, i;
	gsl_complex z;
	
	qs = quantum_system_new();
	quantum_system_add(qs, 2);
	
	qsn = quantum_system_nstates(qs);
	sx1 = operator_sigma_x(qs, 0);
	u = operator_unit(qs,0);
	
	rho1 = density_matrix_pure_TLS(0.74, 0.26);
	rho2 = density_matrix_pure_TLS(0.5, 0.5);
	
	quantum_system_add(qs, 2);
	qsn = quantum_system_nstates(qs);
	quantum_system_print(qs);
	
        density_matrix_list_init(&dm_list);
        density_matrix_list_append(&dm_list, rho1);
        density_matrix_list_append(&dm_list, rho2);

	dm_full = density_matrix_combine(qs, &dm_list);
	
	printf("dm_full =\n");
	misc_matrix_complex_print(dm_full, qsn);

	expm_test = gsl_matrix_complex_alloc(qsn, qsn);
	
	math_expm_complex(dm_full, expm_test);
	printf("expm_test =\n");
	misc_matrix_complex_print(expm_test, qsn);

	GSL_SET_COMPLEX(&z, 0.0, 1.0);
	gsl_matrix_complex_scale(dm_full, z);

	math_expm_complex(dm_full, expm_test);
	printf("expm_test_real =\n");
	misc_matrix_complex_print(expm_test, qsn);
	printf("expm_test_imag =\n");
	misc_matrix_complex_print_imag(expm_test, qsn);
	
	return 0;
}
