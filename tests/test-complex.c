//
// Test complex number abstraction
//
//

#include <stdio.h>
#include "qdpack_matrix_cxsparse.h"

int
main(int argc, char *argv)
{
    qdpack_complex z1, z2, z3;

    QDPACK_SET_COMPLEX(&z1, 1.3, 0.3);
    z2 = qdpack_complex_rect(4.3, 1.2);

    printf("z1 = %.3f + %.3fi\n", QDPACK_REAL(z1), QDPACK_IMAG(z1));
    printf("z2 = %.3f + %.3fi\n", QDPACK_REAL(z2), QDPACK_IMAG(z2));

    z3 = qdpack_complex_add(z1, z2);

    printf("z1+z2 = %.3f + %.3fi\n", QDPACK_REAL(z3), QDPACK_IMAG(z3));

    z3 = qdpack_complex_sub(z1, z2);

    printf("z1-z2 = %.3f + %.3fi\n", QDPACK_REAL(z3), QDPACK_IMAG(z3));

    return 0;
}
