/* 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002 Gerard Jungman, Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */
/* Author:  R. Johansson */

/* Calculate the matrix exponential, following
 * Moler + Van Loan, SIAM Rev. 20, 801 (1978).
 */

#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_mode.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_linalg.h>
#include "gsl_ext_expm_complex.h"

/* store one of the suggested choices for the
 * Taylor series / square  method from Moler + VanLoan
 */
struct moler_vanloan_optimal_suggestion
{
  int k;
  int j;
};
typedef  struct moler_vanloan_optimal_suggestion  mvl_suggestion_t;


/* table from Moler and Van Loan
 * mvl_tab[gsl_mode_t][matrix_norm_group]
 */
static mvl_suggestion_t mvl_tab[3][6] =
{
  /* double precision */
  {
    { 5, 1 }, { 5, 4 }, { 7, 5 }, { 9, 7 }, { 10, 10 }, { 8, 14 }
  },

  /* single precision */
  {
    { 2, 1 }, { 4, 0 }, { 7, 1 }, { 6, 5 }, { 5, 9 }, { 7, 11 }
  },

  /* approx precision */
  {
    { 1, 0 }, { 3, 0 }, { 5, 1 }, { 4, 5 }, { 4, 8 }, { 2, 11 }
  }
};


inline
static double
sup_norm(const gsl_matrix_complex * A)
{
    unsigned int i, j;
      double max = 0.0;
    
    for (i = 0; i < A->size1; i++)
    {
        for (j = 0; j < A->size2; j++)
        {
            gsl_complex z = gsl_matrix_complex_get(A, i, j);

            if (gsl_complex_abs(z) > max)
                max = gsl_complex_abs(z);
        }
    }
  
      return max;
}


static
mvl_suggestion_t
obtain_suggestion(const gsl_matrix_complex * A, gsl_mode_t mode)
{
  const unsigned int mode_prec = GSL_MODE_PREC(mode);
  const double norm_A = sup_norm(A);
  if(norm_A < 0.01) return mvl_tab[mode_prec][0];
  else if(norm_A < 0.1) return mvl_tab[mode_prec][1];
  else if(norm_A < 1.0) return mvl_tab[mode_prec][2];
  else if(norm_A < 10.0) return mvl_tab[mode_prec][3];
  else if(norm_A < 100.0) return mvl_tab[mode_prec][4];
  else if(norm_A < 1000.0) return mvl_tab[mode_prec][5];
  else
  {
    /* outside the table we simply increase the number
     * of squarings, bringing the reduced matrix into
     * the range of the table; this is obviously suboptimal,
     * but that is the price paid for not having those extra
     * table entries
     */
    const double extra = log(1.01*norm_A/1000.0) / M_LN2;
    const int extra_i = (unsigned int) ceil(extra);
    mvl_suggestion_t s = mvl_tab[mode][5];
    s.j += extra_i;
    return s;
  }
}


/* use series representation to calculate matrix exponential;
 * this is used for small matrices; we use the sup_norm
 * to measure the size of the terms in the expansion
 */
static void
matrix_complex_exp_series(
  const gsl_matrix_complex * B,
  gsl_matrix_complex * eB,
  int number_of_terms
  )
{
  int count;
  gsl_matrix_complex *temp = gsl_matrix_complex_calloc(B->size1, B->size2);
  gsl_complex z;
  gsl_complex alpha, beta;

  GSL_SET_COMPLEX(&alpha, 1.0, 0.0);
  GSL_SET_COMPLEX(&beta, 0.0, 0.0);
  
  /* init the Horner polynomial evaluation,
   * eB = 1 + B/number_of_terms; we use
   * eB to collect the partial results
   */  
  gsl_matrix_complex_memcpy(eB, B);
  GSL_SET_COMPLEX(&z, 1.0/number_of_terms, 0.0);
  gsl_matrix_complex_scale(eB, z);
  GSL_SET_COMPLEX(&z, 1.0, 0.0);
  gsl_matrix_complex_add_diagonal(eB, z);
  for(count = number_of_terms-1; count >= 1; --count)
  {
    /*  mult_temp = 1 + B eB / count  */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, B, eB, beta, temp);
    GSL_SET_COMPLEX(&z, 1.0/count, 0.0);
    gsl_matrix_complex_scale(temp, z);
    GSL_SET_COMPLEX(&z, 1.0, 0.0);
    gsl_matrix_complex_add_diagonal(temp, z);

    /*  transfer partial result out of temp */
    gsl_matrix_complex_memcpy(eB, temp);
  }

  /* now eB holds the full result; we're done */
  gsl_matrix_complex_free(temp);
}


int
gsl_linalg_complex_exponential_ss(
  const gsl_matrix_complex * A,
  gsl_matrix_complex * eA,
  gsl_mode_t mode
  )
{
  if(A->size1 != A->size2)
  {
    GSL_ERROR("cannot exponentiate a non-square matrix", GSL_ENOTSQR);
  }
  else if(A->size1 != eA->size1 || A->size2 != eA->size2)
  {
    GSL_ERROR("exponential of matrix must have same dimension as matrix", GSL_EBADLEN);
  }
  else
  {
    int i;
    const mvl_suggestion_t sugg = obtain_suggestion(A, mode);
    const double divisor = exp(M_LN2 * sugg.j);
    gsl_complex z_divisor;
    gsl_complex alpha, beta;

    gsl_matrix_complex * reduced_A = gsl_matrix_complex_alloc(A->size1, A->size2);

    GSL_SET_COMPLEX(&alpha, 1.0, 0.0);
    GSL_SET_COMPLEX(&beta,  0.0, 0.0);
    
    /*  decrease A by the calculated divisor  */
    gsl_matrix_complex_memcpy(reduced_A, A);
    GSL_SET_COMPLEX(&z_divisor, 1.0/divisor, 0.0);
    gsl_matrix_complex_scale(reduced_A, z_divisor);

    /*  calculate exp of reduced matrix; store in eA as temp  */
    matrix_complex_exp_series(reduced_A, eA, sugg.k);

    /*  square repeatedly; use reduced_A for scratch */
    for(i = 0; i < sugg.j; ++i)
    {
      gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, eA, eA, beta, reduced_A);
      gsl_matrix_complex_memcpy(eA, reduced_A);
    }

    gsl_matrix_complex_free(reduced_A);

    return GSL_SUCCESS;
  }
}

