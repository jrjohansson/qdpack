//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------

//
// A GSL sparse matrix interface for the UMFPACK library.
//

#include "math.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>

#include "gsl_ext_sparse.h"

#include "umfpack.h"


/* use a cheap approximate absolute value for complex numbers: */
#define ABS(x,z) ((x) >= 0 ? (x) : -(x)) + ((z) >= 0 ? (z) : -(z))

#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif


/* -----------------------------------------------------------------------------------
 * look for the element with index (i,j) in the triplet vectors. if found, return
 * the index of the element.
 *
 */
static int
triplet_find_element(gsl_sparse_matrix_complex *spmat, long i, long j)
{
    int n;

    for (n = 0; n < spmat->nz; n++)
    {
        if (spmat->row[n] == i && spmat->col[n] == j)
        {
            return n;
        }
    }

    return -1;
}

static double resid(UF_long transpose, UF_long Ap[ ], UF_long Ai[ ], double Ax[ ], double Az[ ], double *x, double *xz, double *b, double *bz, int n)
{
    UF_long i, j, p ;
    double norm ;

    static double *r;
    static double *rz;

    r  = (double *)malloc(sizeof(double) * n);
    rz = (double *)malloc(sizeof(double) * n);


    for (i = 0 ; i < n ; i++)
    {
        r [i] = -b [i] ;
        rz[i] = -bz[i] ;
    }
    if (transpose)
    {
        for (j = 0 ; j < n ; j++)
        {
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                i = Ai [p] ;
                /* complex: r(j) += conj (Aij) * x (i) */
                r [j] += Ax [p] * x [i] ;
                r [j] += Az [p] * xz[i] ;
                rz[j] -= Az [p] * x [i] ;
                rz[j] += Ax [p] * xz[i] ;
            }
        }
    }
    else
    {
        for (j = 0 ; j < n ; j++)
        {
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
                i = Ai [p] ;
                r [i] += Ax [p] * x [j] ;
                r [i] -= Az [p] * xz[j] ;
                rz[i] += Az [p] * x [j] ;
                rz[i] += Ax [p] * xz[j] ;
            }
        }
    }
    norm = 0. ;
    for (i = 0 ; i < n ; i++)
    {
        norm = MAX (ABS (r [i], rz [i]), norm);
    }
    return (norm);
}




/* -----------------------------------------------------------------------------------
 * Allocate a sparse matrix structure for a [n x m] with nelem non-zero elements.
 *
 *
 */
gsl_sparse_matrix_complex *
gsl_sparse_matrix_complex_alloc(int n, int m, int nelem)
{
    gsl_sparse_matrix_complex *spmat;

    if ((spmat = malloc(sizeof(gsl_sparse_matrix_complex))) == NULL)
    {
        fprintf(stderr, "%s: Failed to allocate memory for the sparse matrix data structure.\n", __PRETTY_FUNCTION__);
        return NULL;
    }

    spmat->M_dense = NULL;
        
    // real sparse stuff
    spmat->n = n;
    spmat->m = m;
    spmat->nz = 0;
    spmat->nelem = nelem;
    
    spmat->row = (UF_long *)malloc((nelem) * sizeof(UF_long));
    spmat->col = (UF_long *)malloc((nelem) * sizeof(UF_long));

    spmat->val_real = (double *)malloc((nelem) * sizeof(double));
    spmat->val_imag = (double *)malloc((nelem) * sizeof(double));


    //
    // type of data storage
    //
    spmat->type = SPMAT_TYPE_TRIPLET;    

    umfpack_zl_defaults(spmat->Control);
    spmat->Control[UMFPACK_PRL] = 6;
        
    return spmat;    
}

/* -----------------------------------------------------------------------------------
 * Allocate a sparse matrix structure for a [n x m] with nelem non-zero elements.
 *
 *
 */
void
gsl_sparse_matrix_complex_free(gsl_sparse_matrix_complex *spmat)
{
    if (spmat == NULL)
    {
        return;
    }
    
    if (spmat->M_dense != NULL)
    {
        gsl_matrix_complex_free(spmat->M_dense);
    }
    
    free(spmat);
}


/* -----------------------------------------------------------------------------------
 * Retreive an element from the sparse matrix.
 *
 *
 */
gsl_complex
gsl_sparse_matrix_complex_get(gsl_sparse_matrix_complex *spmat, long i, long j)
{
    gsl_complex z;
    int n;
    
    if (spmat == NULL)
    {
        return z;
    }
    
    //
    // Search for element in triplet form,
    // return 0 if not found:
    //
        if ((n = triplet_find_element(spmat, i, j)) != -1)
    {
        GSL_SET_COMPLEX(&z, spmat->val_real[n], spmat->val_imag[n]);
    }
    else
    {
        GSL_SET_COMPLEX(&z, 0.0, 0.0);
    }

    return z;
}


/* -----------------------------------------------------------------------------------
 * Set an element in the sparse matrix.
 *
 *
 */
int
gsl_sparse_matrix_complex_set(gsl_sparse_matrix_complex *spmat, long i, long j, gsl_complex z)
{
    int n;

    if (spmat == NULL)
    {
        return -1;
    }

    //
    // Add the element to the sparse matrix in triplet form:
    //
    
    if (spmat->nz >= spmat->nelem)
    {
        printf("%s: INFO: reallocating the sparse triplet form vectors.\n", __PRETTY_FUNCTION__);
        spmat->nelem *= 2;

        spmat->row      = (UF_long *)realloc(spmat->row,      (spmat->nelem) * sizeof(UF_long));
        spmat->col      = (UF_long *)realloc(spmat->col,      (spmat->nelem) * sizeof(UF_long));
        spmat->val_real =  (double *)realloc(spmat->val_real, (spmat->nelem) * sizeof(double));
        spmat->val_imag =  (double *)realloc(spmat->val_imag, (spmat->nelem) * sizeof(double));
        
    }
    
    //
    // Check if the i,j element is already present in the vectors...
    //
    n = triplet_find_element(spmat, i, j);

    if (n == -1)
    {
        if (GSL_REAL(z) == 0.0 && GSL_IMAG(z) == 0.0)
        {
            return 0;
        }

        //
        // the (i,j) element does not yet exist in the sparse
        // matrix, so add it now:
        //
        spmat->val_real[spmat->nz] = GSL_REAL(z);
        spmat->val_imag[spmat->nz] = GSL_IMAG(z);
    
        spmat->row[spmat->nz] = i;
        spmat->col[spmat->nz] = j;

        spmat->nz++;
    }
    else
    {
        // 
        // the (i,j) element is already added to the sparse
        // matrix, so we replace it.
        //
        spmat->val_real[n] = GSL_REAL(z);
        spmat->val_imag[n] = GSL_IMAG(z);
    }

    return 0;
}

/* -----------------------------------------------------------------------------------
 * Convert a sparse matrix to dense format (standard GSL matrix). 
 *
 */
gsl_matrix_complex *
gsl_sparse_matrix_complex_convert_to_dense(gsl_sparse_matrix_complex *spmat)
{
    int i;

    if (spmat == NULL)
    {
        return NULL;
    }

    spmat->M_dense = gsl_matrix_complex_alloc(spmat->n, spmat->m);
    gsl_matrix_complex_set_zero(spmat->M_dense);

    //
    // triplet form to dense form:
    //
    for (i = 0; i < spmat->nz; i++)
    {
        gsl_complex z;
        GSL_SET_COMPLEX(&z, spmat->val_real[i], spmat->val_imag[i]);
        gsl_matrix_complex_set(spmat->M_dense, spmat->row[i], spmat->col[i], z);
    }
    
    return spmat->M_dense;
}

/* -----------------------------------------------------------------------------------
 * Convert a sparse matrix from triplet to column format.
 *
 */
int
gsl_sparse_matrix_complex_convert_col(gsl_sparse_matrix_complex *spmat)
{
    static char *__PRETTY_FUNCTION__ = "gsl_sparse_matrix_complex_convert_col";
    
    UF_long status;

    spmat->Ap = (UF_long *)malloc((spmat->n+1) * sizeof(UF_long));
    spmat->Ai = (UF_long *)malloc(spmat->nz    * sizeof(UF_long));
    spmat->Ax = (double *) malloc(spmat->nz    * sizeof(double));
    spmat->Az = (double *) malloc(spmat->nz    * sizeof(double));

    if (!spmat->Ap || !spmat->Ai || !spmat->Ax || !spmat->Az)
    {
        fprintf(stderr, "%s: Out of memory", __PRETTY_FUNCTION__);
        return -1;
    }

    status = umfpack_zl_triplet_to_col(spmat->n, spmat->m, spmat->nz,
                   spmat->row, spmat->col, spmat->val_real, spmat->val_imag,
                   spmat->Ap,  spmat->Ai,  spmat->Ax,       spmat->Az,
                   (UF_long *)NULL);

    if (status < 0)
    {
        umfpack_zl_report_status(spmat->Control, status);
        fprintf(stderr, "gsl_sparse_matrix_complex_convert_col: umfpack_zl_triplet_to_col failed\n");
        return -1;
    }

    //
    // Now that we have converted the triplet form matrix to column form
    // we might want to free up the vectors used for storing the matrix
    // in the triplet form.
    //

    return 0;
}


/* -----------------------------------------------------------------------------------
 * Print column-form sparse matrix to stdout.
 *
 */
int
gsl_sparse_matrix_complex_print_col(gsl_sparse_matrix_complex *spmat)
{
    /* print the column-form of A */
    printf("A (column [%dx%d] %d): ", spmat->n, spmat->m, (int)spmat->nz);
    umfpack_zl_report_matrix(spmat->n, spmat->m, spmat->Ap, spmat->Ai, spmat->Ax, spmat->Az, 1, spmat->Control);

    return 0;    
}



/* -----------------------------------------------------------------------------------
 * Print triplet-form sparse matrix to stdout.
 *
 */
int
gsl_sparse_matrix_complex_print_triplet(gsl_sparse_matrix_complex *spmat)
{
    /* print the triplet-form of A */
    printf("A (triplet [%dx%d] %d): ", spmat->n, spmat->m, (int)spmat->nz);
    umfpack_zl_report_triplet(spmat->n, spmat->n, spmat->nz, spmat->row, spmat->col, spmat->val_real, spmat->val_imag, spmat->Control);

    return 0;    
}


/* -----------------------------------------------------------------------------------
 * Solve: M x = b
 *
 */
int
gsl_sparse_matrix_complex_LU_solve(gsl_sparse_matrix_complex *spmat, double *b_real, double *b_imag, double *x_real, double *x_imag)
{
    void *Symbolic, *Numeric;
    int status;

    //gsl_sparse_matrix_complex_print_col(spmat);

    /* --- symbolic factorization --- */

    status = umfpack_zl_symbolic(spmat->n, spmat->n, spmat->Ap, spmat->Ai, spmat->Ax, spmat->Az, &Symbolic, spmat->Control, spmat->Info);
   
    if (status < 0)
        {
            umfpack_zl_report_info(spmat->Control, spmat->Info);
            //umfpack_zl_report_status(spmat->Control, status);
            fprintf(stderr, "%s: umfpack_zl_symbolic failed\n", __PRETTY_FUNCTION__);
        return -1;
        }

        //printf("%s: Symbolic factorization of A: ", __PRETTY_FUNCTION__);
    //umfpack_zl_report_symbolic(Symbolic, spmat->Control);

    /* --- numeric factorization --- */

    status = umfpack_zl_numeric(spmat->Ap, spmat->Ai, spmat->Ax, spmat->Az, Symbolic, &Numeric, spmat->Control, spmat->Info);

    if (status < 0)
    {
            umfpack_zl_report_info(spmat->Control, spmat->Info);
            umfpack_zl_report_status(spmat->Control, status);
            fprintf(stderr, "%s: umfpack_zl_numeric failed", __PRETTY_FUNCTION__);
        return -1;
        }

        //printf("%s: Numeric factorization of A: ", __PRETTY_FUNCTION__);
    //umfpack_zl_report_numeric(Numeric, spmat->Control);    

    /* --- Solve M x = b --- */
    
    status = umfpack_zl_solve(UMFPACK_A, spmat->Ap, spmat->Ai, spmat->Ax, spmat->Az, x_real, x_imag, b_real, b_imag, Numeric, spmat->Control, spmat->Info);

    //umfpack_zl_report_info(spmat->Control, spmat->Info);
        //umfpack_zl_report_status(spmat->Control, status);

        if (status < 0)
        {
            fprintf(stderr, "%s: umfpack_zl_solve failed\n", __PRETTY_FUNCTION__);
        }
       //printf("%s: x (solution of Ax=b): ") ;
    //umfpack_zl_report_vector(spmat->n, x_real, x_imag, spmat->Control);
    
    {
        double rnorm = resid(FALSE, spmat->Ap, spmat->Ai, spmat->Ax, spmat->Az, x_real, x_imag, b_real, b_imag, spmat->n);
           printf ("maxnorm of residual: %g\n\n", rnorm) ;
    }
    
    return 0;
}



