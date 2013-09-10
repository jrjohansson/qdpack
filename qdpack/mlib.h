//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

/*
 * A stand-alone library of math functions.
 */
#ifndef MLIB_H
#define MLIB_H

#include <math.h>
#include <ctype.h>
#include <complex.h>

/*
 * Polynomials. Hermite polynomials.
 */
struct _poly {
        int n;
        long double *coeff;
} typedef poly;

int poly_print(poly *p);
poly *poly_hermite(int N);
poly *poly_new(int n);
int poly_free(poly *p);
long double poly_eval(poly *p, long double z);

/*
 * Factorial and product sequences
 */
unsigned long long fact_l(unsigned long long n);
long double fact_d(long double n);
long double fact_sqrt(long double n);
long double fact_seq(long double s, long double f);

/*
 * various helper functions
 */
#define MAX(a,b)    ((a<b)?b:a)
#define MIN(a,b)    ((a<b)?a:b)

double heaviside_function(double x);
double boltzmann_factor(double x, double x_env);

/* 
 * Quandrature functions
 */
#define QUAD_TOL    1e-7

long double quad_adaptive_real(long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param);
complex long double quad_adaptive_complex(complex long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param);

long double quad_adaptive_real_inf(long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param);
complex long double quad_adaptive_complex_inf(complex long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param);
    
#endif /* MLIB_H */
