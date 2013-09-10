#------------------------------------------------------------------------------- 
# Copyright (C) 2012, Robert Johansson <robert@riken.jp>
# All rights reserved.
#
# This file is part of QDpack, and licensed under the LGPL.
# http://dml.riken.jp/~rob/qdpack.html
#-------------------------------------------------------------------------------

# ubuntu desktop
CFLAGS= -Wall -I/usr/include/suitesparse -O3 -std=gnu99 -g -Iqdpack -I..
LDFLAGS=-lm -lumfpack -lamd -lgsl -lf77blas -llapack -lcblas -latlas -lcxsparse
#CC=clang
#CC=gcc

# bento and tofu clusters
#CFLAGS=-I/beowulf/include/ -I/beowulf/include/amd -Ilib -O3
#LDFLAGS=-L/beowulf/lib -L/beowulf/lib/atlas -llapack -lcblaswr -lm -lF77 -lI77 -lf77blas -lgsl -lcblas -latlas -lumfpack -lamd

# RSCC
#CFLAGS=-DRSCC
#LDFLAGS=-llapack -latlas -lcblas -lblas -lgsl -lgslcblas -lm
#ARFLAGS=-pc
#CC=cc -gnu2 -static

all: 
	(cd qdpack   && CC='$(CC)' CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)' make)
	(cd examples && CC='$(CC)' CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)' make)
	(cd tests && CC='$(CC)' CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)' make)

install:
	(cd qdpack   && CC='$(CC)' CFLAGS='$(CFLAGS)' LDFLAGS='$(LDFLAGS)' make install)


doc:
	doxygen config.dox

clean:
	(cd qdpack && make clean)
	(cd examples && make clean)
