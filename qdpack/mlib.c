//------------------------------------------------------------------------------
// Copyright (C) 2012, J Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

#include <stdlib.h>
#include <strings.h>
#include <stdio.h>

#include <ctype.h>
#include <math.h>

#include <qdpack/mlib.h>

/*---------------------------------------------------------------------------------------- 
 * Polynomial and Hermite polynomial functions 
 *
 */

/**
 * poly_free: Free memory associated with the polynomial p.
 *
 * @param p    data structure that represent a polynomial.
 */
int
poly_free(poly *p)
{
    if (p)
    {
        if (p->coeff)
        {
            free(p->coeff);
        }
        free(p);
    }
    
    return 0;
}

/**
 * Create a polynomial of degree n.
 *
 * @param n     degree of the polynomial.
 */
poly *
poly_new(int n)
{
    poly *pc;

    if ((pc = (poly *)malloc((size_t)sizeof(poly))) == NULL)
    {
        return NULL;
    }
    if ((pc->coeff = (long double *)malloc((n+1)*sizeof(long double))) == NULL)
    {
        free(pc);
        return NULL; 
    }
    pc->n = n;
    bzero(pc->coeff, sizeof(long double) * n);
    return pc;
}

/*
 * poly_eval: Evaluate the polynomial
 *
 * Parameters: pc - the polynomial data structure
 * Parameters:  z - value for which the polynomial shoud be evaluated
 */
long double
poly_eval(poly *pc, long double z)
{
    int i;
    long double res;

    res = pc->coeff[0];
    for (i = 1; i <= pc->n; i++)
    {
        res += pc->coeff[i] * pow(z, i);
    }

    return res;
}

/*
 * poly_hermite: Generate a polynomial data structure that describe a Hermitian
 *               polynomial of degree N.
 * 
 */
poly *
poly_hermite(int N)
{
    poly *k, *km1, *km2;
    int n, i, p;
    
    if (N == 0)
    {
        k = poly_new(0);
        k->coeff[0] = 1.0;
    }
    else if (N == 1)
    {
        k = poly_new(1);
        k->coeff[0] = 0.0;
        k->coeff[1] = 2.0;
    }
    else
    {
        k   = poly_new(N+1);
        km1 = poly_new(N+1);
        km2 = poly_new(N+1);

        km1->coeff[0] = 1.0;
        km2->coeff[0] = 0.0;
        km2->coeff[1] = 2.0;

        for (n = 1; n <= N+1; n++)
        {
            k->coeff[0] = -2.0 * n * km1->coeff[0];
            for (i = 1; i <= N; i++)
            {
                k->coeff[i] = 2 * (km2->coeff[i-1] - n * km1->coeff[i]);
            }
    

            for (p = 0; p <= N+1; p++)
            {
                long double tmp = km1->coeff[p];
                km1->coeff[p] = km2->coeff[p];
                km2->coeff[p] = k->coeff[p];
                k->coeff[p] = tmp;
            }
        }
    
        poly_free(km1);
        poly_free(km2);
    }
    
    k->n = N;
    
    return k;
}

/*
 * Print the polynomial described by pc to the standard output.
 * 
 */
int
poly_print(poly *pc)
{
    int i;

    printf("H_%d(x) = ", pc->n);
    for (i = 0; i <= pc->n; i++)
    {
        printf("+ %d x^%d ", (int)(pc->coeff[i]), i);
    }
    printf("\n");

    return 0;
}


/* -----------------------------------------------------------------------------------------
 * Factorial and product sequences 
 *
 */

/*
 * Factorial of a long integer
 */
unsigned long long
fact_l(unsigned long long n)
{
    if (n == 1 || n == 0)
        return 1;
    else
        return n * fact_l(n-1); 
}

/*
 * Factorial in double precision.
 */
long double
fact_d(long double n)
{
    if (n <= 1.0)
        return 1.0;
    else if (n > 20.0)
    {
        printf("WARNING: calculating factorial of a too large number\n");
        return sqrt(2*M_PI) * pow(n, n+1/2) * exp(-n);
    }
    else
        return n * fact_d(n-1.0); 
}

/*
 * Calculate the square root of a factorial: sqrt(n!)
 */
long double
fact_sqrt(long double n)
{
    if (n <= 1.0)
        return 1.0;
    else
        return sqrt(n) * fact_sqrt(n-1.0); 
}

/*
 * Calculate product sequence: s*(s+1)*(s+2)*...(f-2)*(f-1)*f
 */
long double
fact_seq(long double s, long double f)
{
    if (f - s > 20)
        printf("fact_seq: warning: too large sequence.\n");

    if (s > f)
        return 1;
    else if (s == f)
        return f;
    else
        return s * fact_seq(s+1, f); 
}


/*
 *
 */
double
heaviside_function(double x)
{
    if (x > 0.0)
    {
        return 1.0;
    }
    
    return 0.0;
}

double
boltzmann_factor(double x, double x_env)
{       
    if (x > 0.0 && x_env > 0.0)
    {
        return 1.0 / (exp(x / x_env) - 1);
    }
    else
    {
        return 0.0;
    }
}

/* ------------------------------------------------------------------------------------------------------
 * Quadrature functions. Evaluates integrals: one set of function for real valued function, and one
 * set for complex valued functions.
 *
 */


/*
 * Internal recursive function for adaptive evaluation of complex valued integral.
 *
 */
complex long double
quad_adaptive_complex_rec(complex long double (*f)(long double, void *), void *param, 
              long double a, long double b, complex long double fa, 
              complex long double fc, complex long double fb, long double tol)
{
    static int level = 0;
    complex long double Q1, Q2, Q;
    
    complex long double fd, fe;
    long double h, c;

    h = b - a;
    c = (a+b)/2.0;

    fd = f((a+c)/2.0, param);
    fe = f((c+b)/2.0, param);

    Q1 = h/6.0  * (fa + 4*fc + fb);
    Q2 = h/12.0 * (fa + 4*fd + 2*fc + 4*fe + fb);

    if (fabsl(creal(Q1-Q2)) <= tol && fabsl(cimag(Q1-Q2)) <= tol)
    {
        Q = Q2 + (Q2 - Q1)/15.0;
    }
    else
    {
        level++;
        //printf("quad_adaptive_complex_rec: Refining interval (%f, %f : %g) : [%g, %g] - [%g, %g] = [%g, %g].\n", 
        //       a, b, b-a, creal(Q2), cimag(Q2), creal(Q1), cimag(Q1), creal(Q2-Q1), cimag(Q2-Q1));
        Q1 = quad_adaptive_complex_rec(f, param, a, c, fa, fd, fc, tol);
        Q2 = quad_adaptive_complex_rec(f, param, c, b, fc, fe, fb, tol);
        Q = Q1 + Q2;
        level--;
    }
    
    return Q;
}


/*
 * Evaluate complex valued integral over the interval [x_min, x_max].
 *
 */
complex long double
quad_adaptive_complex(complex long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param)
    
{
    complex long double sum = 0.0, dsum = 0.0;

    complex long double fa, fm, fb;
    long double a, b;
    
    int iter = 0;

    /* evaluate integral on finite interval */
    a = x_min;
    b = x_min+dx;

    fa = f(a, param);
        
    do {
        /* simpsons rule (Ref. Heath, p 241) applied to sub-interval in left and right tail */
        
        fb = f(b, param);
        fm = f((a+b)/2, param);
        
        //dsum = (b-a)/6.0 * (fa + 4*fm + fb);

        dsum = quad_adaptive_complex_rec(f, param, a, b, fa, fm, fb, QUAD_TOL);
        
        fa = fb;
        a  = b;
        b += dx;
        
        sum += dsum;
        iter++;
    } while (b <= x_max);

    return sum;
}


/*
 * Evaluate complex valued integral from -inf to inf, where the function is monotonically decaying 
 * outside the region [xmin, xmax].
 * 
 */
complex long double
quad_adaptive_complex_inf(complex long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param)
    
{
    complex long double sum = 0.0, sum_trunc = 0.0, dsum = 0.0;

    complex long double faL, fmL, fbL;
    long double aL, aR, bL, bR;

    int iter = 0;

    /* evaluate integral on finite interval */
    sum = quad_adaptive_complex(f, x_min, x_max, dx, param);
    
    /* evaluate correction due to truncation */
    aL = x_min - dx;
    bL = x_min;
    aR = x_max;
    bR = x_max + dx;

    fbL = f(bL, param);
    //faR = f(aR, param);
    
    iter = 0;
    
    do {
        /* simpsons rule (Ref. Heath, p 241) applied to sub-interval in left and right tail */
        
        fmL = f((aL+bL)/2.0, param);
        faL = f(aL, param);
        
        //fmR = f((aR+bR)/2.0, param);
        //fbR = f(bR, param);
        
        dsum = (bL-aL)/6.0 * (faL + 4*fmL + fbL) + (bL-aL)/6.0 * (faL + 4*fmL + fbL);    
    
        aR = bR;
        bR = bR + dx;
        //faR = fbR;

        bL = aL;
        aL = aL - dx;
        fbL = faL;
        
        sum_trunc += dsum;
    
        printf("quad_adpt_inf: correcting for pre-mature truncation: (%Lf, %Lf) (%Lf, %Lf) %Lf [%f, %f]\n", aL, bL, aR, bR, dx, creal(dsum), cimag(dsum));
        
        iter++;
    } while (fabsl(creal(dsum)) >= QUAD_TOL*10 || fabsl(cimag(dsum)) >= QUAD_TOL*10);

    printf("quad_adpt_inf: add correction due to truncation error: (%f, %f) [%Lf, %d]\n", creal(sum_trunc), cimag(sum_trunc), bR, iter);
    
    return sum+sum_trunc;
}


/*
 * Real valued functions. Internal recursive function for adaptive evaluation of an integral.
 *
 */
long double
quad_adaptive_real_rec(long double (*f)(long double, void *), void *param, long double a, long double b, long double fa, long double fc, long double fb, long double tol)
{
    static int level = 0;
    long double Q1, Q2, Q;
    
    long double fd, fe, h, c;

    if (level > 30)
    {
        printf("Exceedingly high recursion level = %d:\t (%Le, %Le) = (%Lf, %Lf)\n", level, a, b, fa, fb);
    }
    
    h = b - a;
    c = (a+b)/2.0;

    fd = f((a+c)/2.0, param);
    fe = f((c+b)/2.0, param);

    Q1 = h/6.0  * (fa + 4*fc + fb);
    Q2 = h/12.0 * (fa + 4*fd + 2*fc + 4*fe + fb);

    if (fabsl(Q1-Q2) <= tol)
    {
        Q = Q2 + (Q2 - Q1)/15.0;
    }
    else
    {
        level++;
        //if (level > 16)
        //printf("quad_adaptive_complex_rec: [%d] Refining interval (%Lf, %Lf : %Lg) (%Le, %Le): [%g] - [%g] = [%g].\n", 
        //       level, a, b, b-a, fa, fb, Q2, Q1, fabsl(Q2-Q1));
        Q1 = quad_adaptive_real_rec(f, param, a, c, fa, fd, fc, tol);
        Q2 = quad_adaptive_real_rec(f, param, c, b, fc, fe, fb, tol);
        Q = Q1 + Q2;
        level--;
    }
    
    return Q;
}


/*
 * Evaluate an real-valued integral from xmin to xmax.
 * 
 */
long double
quad_adaptive_real(long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param)
{
    long double sum = 0.0, dsum = 0.0;
    long double fa, fm, fb, a, b;

    int iter = 0;

    /* evaluate integral on finite interval */
    a = x_min;
    b = x_min+dx;

    fa = f(a, param);
        
    do {
        /* simpsons rule (Ref. Heath, p 241) applied to sub-interval in left and right tail */
        
        fb = f(b, param);
        fm = f((a+b)/2, param);
        
        //dsum = (b-a)/6.0 * (fa + 4*fm + fb);

        dsum = quad_adaptive_real_rec(f, param, a, b, fa, fm, fb, QUAD_TOL);
        
        fa = fb;
        a  = b;
        b += dx;
        
        sum += dsum;
        iter++;
    } while (b <= x_max);

    return sum;
}    


/*
 * Evaluate a real-valued integral from -inf to inf, which is monotonically decaying outside the region [xmin, xmax].
 * 
 */
long double
quad_adaptive_real_inf(long double (*f)(long double,void *), long double x_min, long double x_max, long double dx, void *param)
{
    long double sum = 0.0, sum_trunc = 0.0, dsum = 0.0;
    long double faL, fmL, fbL;//, fbR; //, faR, fmR
    long double aL, aR, bL, bR;

    int iter = 0;

    /* evaluate integral on finite interval */
    sum = quad_adaptive_real(f, x_min, x_max, dx, param);
    
    /* evaluate correction due to truncation */
    aL = x_min - dx;
    bL = x_min;
    aR = x_max;
    bR = x_max + dx;

    fbL = f(bL, param);
    //faR = f(aR, param);
    
    iter = 0;
    
    do {
        /* simpsons rule (Ref. Heath, p 241) applied to sub-interval in left and right tail */
        
        fmL = f((aL+bL)/2.0, param);
        faL = f(aL, param);
        
        //fmR = f((aR+bR)/2.0, param);
        //fbR = f(bR, param);
        
        dsum = (bL-aL)/6.0 * (faL + 4*fmL + fbL) + (bL-aL)/6.0 * (faL + 4*fmL + fbL);    
    
        aR = bR;
        bR = bR + dx;
        //faR = fbR;

        bL = aL;
        aL = aL - dx;
        fbL = faL;
        
        sum_trunc += dsum;
    
        printf("quad_adpt_inf: correcting for pre-mature truncation: (%Lf, %Lf) (%Lf, %Lf) %Lf [%Lf]\n", aL, bL, aR, bR, dx, dsum);
        
        iter++;
    } while (fabsl(dsum) >= QUAD_TOL*10);

    printf("quad_adpt_inf: add correction due to truncation error: (%Lf, %Lf) [%d]\n", sum_trunc, bR, iter);
    
    return sum+sum_trunc;
}    


