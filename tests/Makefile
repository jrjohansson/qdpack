#------------------------------------------------------------------------------- 
# Copyright (C) 2012, Robert Johansson <robert@riken.jp>
# All rights reserved.
#
# This file is part of QDpack, and licensed under the LGPL.
# http://dml.riken.jp/~rob/qdpack.html
#-------------------------------------------------------------------------------

PROGRAMS=test_matrix_fops test_matrix_backend 

all: 
	for i in ${PROGRAMS}; do \
                echo Building  $$i; \
		$(CC) -o $$i $$i.c -I.. $(CFLAGS) ../qdpack/libqdpack.a $(LDFLAGS) ; \
	done
	$(CC) -o test_matrix_backend_gsl test_matrix_backend_gsl.c ../qdpack/qdpack_matrix_gsl.c ../qdpack/gsl_ext_expm_complex.c ../qdpack/gsl_ext.c -I.. $(CFLAGS) $(LDFLAGS)
	$(CC) -o test_matrix_backend_cxsparse test_matrix_backend_cxsparse.c ../qdpack/qdpack_matrix_cxsparse.c -I.. $(CFLAGS) $(LDFLAGS)
	$(CC) -o test_matrix_backend_cxsparse_conversion test_matrix_backend_cxsparse_conversion.c ../qdpack/qdpack_matrix_cxsparse.c -I.. $(CFLAGS) $(LDFLAGS)

install:
	echo "nothing..."

clean:
	rm -rf *.o *~ $(PROGRAMS)

