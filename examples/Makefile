#------------------------------------------------------------------------------- 
# Copyright (C) 2012, Robert Johansson <robert@riken.jp>
# All rights reserved.
#
# This file is part of QDpack, and licensed under the LGPL.
# http://dml.riken.jp/~rob/qdpack.html
#-------------------------------------------------------------------------------

PROGRAMS=run_qubits run_jc run_spinchain run_2ls_driven_pss

all: 
	for i in ${PROGRAMS}; do \
                echo Building  $$i; \
		$(CC) -o $$i $$i.c -I.. $(CFLAGS) ../qdpack/libqdpack.a $(LDFLAGS) ; \
	done

install:
	echo "nothing..."

clean:
	rm -rf *.o *~ $(PROGRAMS)
