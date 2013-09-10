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

#include <stdio.h> 
#include <gsl/gsl_matrix.h>


int
main(int argc, char **argv)
{
	gsl_matrix *sx1, *sx2, *sz1, *sz2, *u, *rho1, *rho2;
	gsl_matrix_complex *dm_full, *dm_part;
	quantum_system *qs;
	quantum_system_state_vector_t qsv;
	density_matrix_list_t dm_list;
	int qsn, i;
	
	qs = quantum_system_new();
	quantum_system_add(qs, 2);
	
	qsn = quantum_system_nstates(qs);
	sx1 = operator_sigma_x(qs, 0);
	u = operator_unit(qs,0);
	
	rho1 = density_matrix_pure_TLS(0.74, 0.26);
	rho2 = density_matrix_pure_TLS(0.5, 0.5);
	
	quantum_system_add(qs, 2);
//	quantum_system_add(qs, 2);
//	quantum_system_add(qs, 2);
	qsn = quantum_system_nstates(qs);
	quantum_system_print(qs);
	
        density_matrix_list_init(&dm_list);
        density_matrix_list_append(&dm_list, rho1);
        density_matrix_list_append(&dm_list, rho2);
//        density_matrix_list_append(&dm_list, u);
//        density_matrix_list_append(&dm_list, u);

	dm_full = density_matrix_combine(qs, &dm_list);
	
	printf("dm_full =\n");
	misc_matrix_complex_print(dm_full, qsn);

	for (i = 0; i < qs->nsubsys; i++)
	{
		dm_part = density_matrix_traceout(qs, dm_full, i);
		printf("dm_part[%d] =\n", i);
		misc_matrix_complex_print(dm_part, qs->nstates[i]);
		gsl_matrix_complex_free(dm_part);
	}

	
	return 0;
}
