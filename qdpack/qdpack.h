//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//
//------------------------------------------------------------------------------

#ifndef _QDPACK
#define _QDPACK

#define USE_DENSE

//
// @mainpage Quantum Dynamics Package
//
// A C library for the time evolution of quantum systems.
// 
// @author J Robert Johansson <robert@riken.jp>
// 
//

#include <qdpack/qdpack_matrix.h>

#include <qdpack/hilbert_space.h>

#include <qdpack/qdpack_object.h>

#include <qdpack/composite.h>
#include <qdpack/operators.h>
#include <qdpack/states.h>

#include <qdpack/simulation.h>
#include <qdpack/hamiltonian.h>

#include <qdpack/master_equation.h>
#include <qdpack/integrate.h>
#include <qdpack/steadystate.h>

#include <qdpack/entanglement.h>

#include <qdpack/basis_transform.h>
#include <qdpack/distribution_functions.h>

#include <qdpack/io.h>
#include <qdpack/mlib.h>

#endif
