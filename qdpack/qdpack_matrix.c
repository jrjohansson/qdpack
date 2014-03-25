//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#include <qdpack/qdpack.h>

// select matrix backend
#ifdef USE_DENSE
#include "qdpack_matrix_gsl.c"
#endif

#ifdef USE_SPARSE
#include "qdpack_matrix_cxsparse.c"
#endif
