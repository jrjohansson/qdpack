/*
 *
 *
 *
 *
 */
#include "quantum_system.h"
#include <stdio.h> 

int
main(int argc, char **argv)
{
	quantum_system *qs;
	quantum_system_state_vector_t qsv;
	int qsn, i;

	qs = quantum_system_new();
	
	quantum_system_print(qs);

	quantum_system_add(qs, 5);
	quantum_system_add(qs, 5);

	quantum_system_print(qs);


	qsn = quantum_system_nstates(qs);

	for (i = 0; i < qsn; i++)
	{
		quantum_system_state_vector(qs, i, qsv);
		
		printf("qsn = %d\t->\t", i);
		quantum_system_state_vector_print(qs, qsv);
		printf("\t->\tqsn = %d\n", quantum_system_state_number(qs, qsv));
		
	}
	
	return 0;
}
