//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

// select matrix backend
#include <qdpack/qdpack.h>

#ifndef USE_DENSE
#ifndef USE_SPARSE
#error Must compile with USE_DENSE or USE_SPARSE
#endif
#endif

#ifdef USE_DENSE
#include <qdpack/qdpack_matrix_gsl.h>
#endif

#ifdef USE_SPARSE
#include <qdpack/qdpack_matrix_cxsparse.h>
#endif
