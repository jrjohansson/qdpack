//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------



#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>



int 
gsl_matrix_complex_add_diagonal (gsl_matrix_complex * a, const gsl_complex x)
{
  const size_t M = a->size1;
  const size_t N = a->size2;
  const size_t tda = a->tda;
  const size_t loop_lim = (M < N ? M : N);
  size_t i;
  for (i = 0; i < loop_lim; i++)
    {
      a->data[2 * (i * tda + i)] += GSL_REAL (x);
      a->data[2 * (i * tda + i) + 1] += GSL_IMAG (x);
    }

  return GSL_SUCCESS;
}


int 
gsl_matrix_complex_scale (gsl_matrix_complex * a, const gsl_complex x)
{
  const size_t M = a->size1;
  const size_t N = a->size2;
  const size_t tda = a->tda;

  size_t i, j;

  double xr = GSL_REAL(x);
  double xi = GSL_IMAG(x);

  for (i = 0; i < M; i++)
    {
      for (j = 0; j < N; j++)
        {
          const size_t aij = 2 * (i * tda + j);

          double ar = a->data[aij];
          double ai = a->data[aij + 1];
          
          a->data[aij] = ar * xr - ai * xi;
          a->data[aij + 1] = ar * xi + ai * xr;
        }
    }

  return GSL_SUCCESS;
}


int 
gsl_matrix_complex_add (gsl_matrix_complex * a, const gsl_matrix_complex * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR ("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t tda_a = a->tda;
      const size_t tda_b = b->tda;

      size_t i, j;

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < N; j++)
            {
              const size_t aij = 2 * (i * tda_a + j);
              const size_t bij = 2 * (i * tda_b + j);

              a->data[aij] += b->data[bij];
              a->data[aij + 1] += b->data[bij + 1];
            }
        }

      return GSL_SUCCESS;
    }
}

int 
gsl_matrix_complex_sub (gsl_matrix_complex * a, const gsl_matrix_complex * b)
{
  const size_t M = a->size1;
  const size_t N = a->size2;

  if (b->size1 != M || b->size2 != N)
    {
      GSL_ERROR ("matrices must have same dimensions", GSL_EBADLEN);
    }
  else
    {
      const size_t tda_a = a->tda;
      const size_t tda_b = b->tda;

      size_t i, j;

      for (i = 0; i < M; i++)
        {
          for (j = 0; j < N; j++)
            {
              const size_t aij = 2 * (i * tda_a + j);
              const size_t bij = 2 * (i * tda_b + j);

              a->data[aij] -= b->data[bij];
              a->data[aij + 1] -= b->data[bij + 1];
            }
        }

      return GSL_SUCCESS;
    }
}

