#------------------------------------------------------------------------------- 
# Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
# All rights reserved.
#
# This file is part of QDpack, and licensed under the LGPL.
# http://dml.riken.jp/~rob/qdpack.html
#-------------------------------------------------------------------------------

QDLIB=mlib.o integrate.o entanglement.o hamiltonian.o hilbert_space.o states.o \
      composite.o operators.o master_equation.o basis_transform.o steadystate.o  \
      distribution_functions.o qdpack_object.o simulation.o io.o floquet.o \
      qdpack_matrix.o gsl_ext_expm_complex.o gsl_ext_expm_complex.o gsl_ext.o

#gsl_ext_sparse.o steady_sparse.o

lib: $(QDLIB)
	ar $(ARFLAGS) libqdpack.a $(QDLIB)

install: lib
	cp libqdpack.a /usr/local/lib/libqdpack.a
	mkdir -p /usr/local/include/qdpack
	cp *.h /usr/local/include/qdpack/

clean:
	rm -rf *.o *~ libqdpack.a
