//------------------------------------------------------------------------------
// Copyright (C) 2012, Robert Johansson <robert@riken.jp>
// All rights reserved.
//
// This file is part of QDpack, and licensed under the LGPL.
// http://dml.riken.jp/~rob/qdpack.html
//------------------------------------------------------------------------------


#include <stdio.h>
#include <math.h>

#include <gsl/gsl_integration.h>

#include <qdpack/qdpack.h>

static const size_t qawo_n=100;
static const size_t qawo_limit=1001;
static double abserr = 1e-6;
static double relerr = 1e-3;

/*
 * Coordinate to index and index to coordinate conversion functions:
 * Status: Done 
 */
static double
coordinate_from_index(double max, int i, int N)
{
    return (i - N/2.0) * max / (N/2.0);
}

static int
coordinate_to_index(double val, double max, int N)
{
    int n = (int)(floor((val + max) / (2*max) * (N)));

//    printf("debug: val+max = %f, 2*max = %f, (val + max) / (2*max) * N = %f, floor(...) = %f\n",
//            (val + max), (2*max), (val + max) / (2*max) * N, floor((val + max) / (2*max) * N));

    if (n < 0)
        n = 0;

    if (n >= N)
        n = N - 1;

//    printf("debug: coordinate_to_index: val = %f => n = %d: Val should be contained in [%f, %f]\n", 
//           val, n, coordinate_from_index(max, n, N), coordinate_from_index(max, n+1, N));

    return n;
}

/*
 * Interpolate an arbitrary point [between grid points] given a grid distribution.
 * Status: Done 
 */
static double
ch_td_wigner_interpolate(qdpack_complex z, void *rawdata)
{
    wigner_data_t *data = (wigner_data_t *)rawdata;
      
    qdpack_complex z00, z10, z01, f_z00, f_z01, f_z10; //, f_z;

    int z0_re_idx, z0_im_idx, z1_re_idx, z1_im_idx;

    double f_z_re; //, f_z_im;

//    printf("debug: ch_td_wigner_interpolate: lambda_max = %f + i %f\n", QDPACK_REAL(data->lambda_max), QDPACK_IMAG(data->lambda_max));
//    printf("debug: ch_td_wigner_interpolate: lambda = %f + i %f [%d]\n", QDPACK_REAL(lambda), QDPACK_IMAG(lambda), data->N);

    z0_re_idx = coordinate_to_index(QDPACK_REAL(z), QDPACK_REAL(data->lambda_max), data->N);
    z0_im_idx = coordinate_to_index(QDPACK_IMAG(z), QDPACK_IMAG(data->lambda_max), data->N);

    z0_re_idx = (z0_re_idx < 0) ? 0 : z0_re_idx;
    z0_im_idx = (z0_im_idx < 0) ? 0 : z0_im_idx;

    z1_re_idx = (z0_re_idx < (data->N-1)) ? z0_re_idx + 1 : z0_re_idx;
    z1_im_idx = (z0_im_idx < (data->N-1)) ? z0_im_idx + 1 : z0_im_idx;

 //   printf("debug: ch_td_wigner_interpolate idx_re(1,2) = (%d -> %d), idx_im(1,2) = (%d -> %d)\n", 
 //           z0_re_idx, z1_re_idx, z0_im_idx, z1_im_idx);
  
    QDPACK_SET_COMPLEX(&z00, 
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z0_re_idx, data->N),
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z0_im_idx, data->N)); 

    QDPACK_SET_COMPLEX(&z01, 
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z0_re_idx, data->N),
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z1_im_idx, data->N)); 

    QDPACK_SET_COMPLEX(&z10, 
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z1_re_idx, data->N),
        coordinate_from_index(QDPACK_REAL(data->lambda_max), z0_im_idx, data->N)); 

    f_z00 = qdpack_matrix_get(data->Cw, z0_re_idx, z0_im_idx);
    f_z01 = qdpack_matrix_get(data->Cw, z0_re_idx, z1_im_idx);
    f_z10 = qdpack_matrix_get(data->Cw, z1_re_idx, z0_im_idx);

    printf("debug: ch_td_wigner_interpolate: z           = %f + i %f. f(z) = ?\n", 
           QDPACK_REAL(z), QDPACK_IMAG(z));
    printf("debug: ch_td_wigner_interpolate: z00[%d, %d] = %f + i %f. f(z00) = %f + i %f\n", 
           z0_re_idx, z0_im_idx, QDPACK_REAL(z00), QDPACK_IMAG(z00), QDPACK_REAL(f_z00), QDPACK_IMAG(f_z00));
    printf("debug: ch_td_wigner_interpolate: z01[%d, %d] = %f + i %f. f(z01) = %f + i %f\n", 
           z0_re_idx, z1_im_idx, QDPACK_REAL(z01), QDPACK_IMAG(z01), QDPACK_REAL(f_z01), QDPACK_IMAG(f_z01));
    printf("debug: ch_td_wigner_interpolate: z10[%d, %d] = %f + i %f. f(z10) = %f + i %f\n", 
           z1_re_idx, z0_im_idx, QDPACK_REAL(z10), QDPACK_IMAG(z10), QDPACK_REAL(f_z10), QDPACK_IMAG(f_z10));


    printf("debug: ch_td_wigner_interpolate: dx = %f [%f], dy = %f [%f]\n", 
           QDPACK_REAL(z)   - QDPACK_REAL(z00),
           QDPACK_REAL(z10) - QDPACK_REAL(z00),
           QDPACK_IMAG(z)   - QDPACK_IMAG(z00),
           QDPACK_IMAG(z01) - QDPACK_IMAG(z00));

    f_z_re = QDPACK_REAL(f_z00);
    if ((z1_re_idx != z0_re_idx) && QDPACK_REAL(z10) != QDPACK_REAL(z00))
    {
        f_z_re += (QDPACK_REAL(f_z10) - QDPACK_REAL(f_z00)) * 
                  (QDPACK_REAL(z) - QDPACK_REAL(z00)) / (QDPACK_REAL(z10) - QDPACK_REAL(z00));
    }
    if ((z1_im_idx != z0_im_idx) && (QDPACK_IMAG(z01) != QDPACK_IMAG(z00)))
    {
        f_z_re += (QDPACK_REAL(f_z01) - QDPACK_REAL(f_z00)) * 
                  (QDPACK_IMAG(z) - QDPACK_IMAG(z00)) / (QDPACK_IMAG(z01) - QDPACK_IMAG(z00));
    }

/*
    f_z_im = QDPACK_IMAG(f_z00);
    if(0) //if (QDPACK_REAL(z10) != QDPACK_REAL(z00))
    {
        f_z_im += (QDPACK_IMAG(f_z10) - QDPACK_IMAG(f_z00)) * 
                  (QDPACK_REAL(z) - QDPACK_REAL(z00)) / (QDPACK_REAL(z10) - QDPACK_REAL(z00));
    }
    if(0) //if (QDPACK_IMAG(z01) != QDPACK_IMAG(z00))
    {
        f_z_im += (QDPACK_IMAG(f_z01) - QDPACK_IMAG(f_z00)) * 
                  (QDPACK_IMAG(z) - QDPACK_IMAG(z00)) / (QDPACK_IMAG(z01) - QDPACK_IMAG(z00));
    }
*/
    if (isnan(f_z_re))
    {
        printf("ch_td_wigner_interpolate: z_re was NaN... \n");
        f_z_re = 0.0;
    }

/*
    if (isnan(f_z_im))
    {
        printf("ch_td_wigner_interpolate: z_im was NaN... \n");
        f_z_im = 0.0;
    }
*/

//    QDPACK_SET_COMPLEX(&f_z, f_z_re, 0.0);
//    QDPACK_SET_COMPLEX(&f_z, exp(-qdpack_complex_abs(z)), 0.0);
//    QDPACK_SET_COMPLEX(&f_z, 0.0, 0.0);

//    return f_z;

   // f_z_re = exp(-qdpack_complex_abs(z));

    return f_z_re;
//    return ;
}


/**
 * \breif    Calculate the Q-function corresponding to a system with density matrix rho.
 *
 */
int
distribution_function_Q(qdpack_matrix_t *Q, qdpack_operator_t *RHO, double alpha_max)
{
    qdpack_state_t *wf;
    qdpack_state_t *prod;
    qdpack_complex z;
    int i, j;

    double a_re, a_im, r, theta;

    if (Q == NULL || RHO == NULL)
    {
        printf("distribution_function_Q: parameter errors\n");
        return -1;
    }

    prod = qdpack_state_alloc(RHO->qs);     
    qdpack_matrix_set_zero(Q);
        
    for (i = 0; i < Q->m; i++)
    {
        a_re = coordinate_from_index(alpha_max, i, Q->m);
        for (j = 0; j < Q->m; j++)
        {
            a_im = coordinate_from_index(alpha_max, j, Q->m);
            
            r     = sqrt(a_re * a_re + a_im * a_im);
            theta = atan(a_im / a_re);
            if (a_re < 0)
            {
                theta = theta + M_PI;
            }

            wf = qdpack_wf_coherent_state(r, theta, RHO->m, 0);

            qdpack_operator_state_multiply(RHO, wf, prod);
            z = qdpack_state_dot(wf, prod);

            // use
            // z = qdpack_matrix_multiply_vmv(RHO->data, wf->data)

            qdpack_matrix_set(Q, i, j, z);

            qdpack_state_free(wf);
        }
    }
    qdpack_state_free(prod);
    
    return 0;       
}


/*
 * Compute the characteristic function corresponding to an harmonic
 * oscillator density matrix:
 */
int
distribution_function_characteristic_w(
            qdpack_matrix_t *phase_space, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex lambda_max)
{
    qdpack_operator_t *a_exp, *ad_exp, *ws;

    int i, j, qsn, N;

    qdpack_complex lambda, z;

    if (!phase_space || !a || !ad || !rho)
    {
        fprintf(stderr, "%s: invalid arguments (null pointer).\n", __PRETTY_FUNCTION__);
        return 0;
    }
    else
    {
        printf("%s: calculating... \n", __PRETTY_FUNCTION__);
    }
    
    N   = phase_space->m;
    qsn = rho->m;

    if ( (rho->n != qsn) ||
         (a->m   != qsn) ||
         (a->n   != qsn) ||
         (ad->m  != qsn) ||
         (ad->n  != qsn) )
    {
        fprintf(stderr, "%s: invalid arguments (size mismatch).\n", __PRETTY_FUNCTION__);
        return 0;
    }

    ws     = qdpack_operator_alloc(rho->qs);
    a_exp  = qdpack_operator_alloc(rho->qs);
    ad_exp = qdpack_operator_alloc(rho->qs);

    // i - real
    // j - imag
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            // calculate the exponents of a and ad:
            QDPACK_SET_COMPLEX(&lambda, 
                coordinate_from_index(QDPACK_REAL(lambda_max), i, N),
                coordinate_from_index(QDPACK_IMAG(lambda_max), j, N));

            qdpack_operator_memcpy(ws, a);
            qdpack_operator_scale(ws, qdpack_complex_negative(qdpack_complex_conjugate(lambda)));
            qdpack_matrix_exp(ws->data, a_exp->data);

            qdpack_operator_memcpy(ws, ad);
            qdpack_operator_scale(ws, lambda);
            qdpack_matrix_exp(ws->data, ad_exp->data);

            qdpack_operator_multiply(ad_exp, a_exp, ws);

            z = qdpack_operator_expectation_value(rho, ws);

            qdpack_matrix_set(phase_space, i, j, qdpack_complex_mul_real(z, exp(-0.5 * qdpack_complex_abs2(lambda))));
        }
    }    

    qdpack_operator_free(ws);
    qdpack_operator_free(a_exp);
    qdpack_operator_free(ad_exp);

    return 0;
}


double
cb_integrand_wigner_cf(qdpack_complex lambda, qdpack_matrix_t *phase_space_cf, qdpack_complex alpha)
{
    return 0.0;    
}

static double
cb_qawo_wigner_f_im(double lambda_im, void *rawdata)
{
//    qdpack_complex z;
    double z;
    wigner_data_t *data = (wigner_data_t *)rawdata;

    // set the real and imaginary part
    QDPACK_SET_COMPLEX(&(data->lambda), QDPACK_REAL(data->lambda), lambda_im);

    z = ch_td_wigner_interpolate(data->lambda, data);

//    printf("debug: cb_wigner_f_im [lambda = %f + i %f]: z = %f + i %f\n", 
//           QDPACK_REAL(data->lambda), QDPACK_IMAG(data->lambda), z, 0.0); //QDPACK_REAL(z), QDPACK_IMAG(z));

//    return QDPACK_REAL(z);
    return z;
}

static double
cb_qawo_wigner_f_re(double lambda_re, void *rawdata)
{
    double result = 1.0, abserr;
    int status;

    wigner_data_t *data = (wigner_data_t *)rawdata;

//    printf("debug [f_re(%f)]: calculating...\n", lambda_re);

    // set the real part
    QDPACK_SET_COMPLEX(&(data->lambda), lambda_re, 0);

//    data->wf_im = gsl_integration_qawo_table_alloc(
//                      2 * QDPACK_REAL(data->alpha), 
//                      2 * QDPACK_IMAG(data->lambda_max), 
//                      GSL_INTEG_COSINE, qawo_n); 
//    gsl_integration_qawo_table_set_length(data->wf_im, 2 * QDPACK_IMAG(data->lambda_max));

    gsl_integration_qawo_table_set(data->wf_im, 2 * QDPACK_REAL(data->alpha), 2 * QDPACK_IMAG(data->lambda_max), GSL_INTEG_COSINE);


    if ((status = gsl_integration_qawo(&(data->f_im), -QDPACK_IMAG(data->lambda_max), abserr, relerr, qawo_limit, data->w_im, data->wf_im, &result, &abserr)) != 0)
    {
        printf("Warning: gsl_integration_qawo failed [status = %d]...\n", status);
    }
    else
    {
    //    printf("debug [f_re(%f)]: gsl_integration_qawo converged to %f with abserr %f\n", lambda_re, result, abserr);
    }        

//    gsl_integration_qawo_table_free(data->wf_im);

    return result;
}


/*
 * Compute the Wigner function corresponding to an harmonic
 * oscillator density matrix, by 2d numberical integration
 * of the characteristic function.
 */
int
distribution_function_wigner_cf_qawo(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max)
{
    //qdpack_matrix_t *phase_space_cfw;
    int i, j, N, status;
    double result, abserr;
    qdpack_complex z;

    wigner_data_t data;

    N   = phase_space_wigner->m;
    data.N = N;
    data.alpha_max  = alpha_max;
    data.lambda_max = alpha_max;

//    phase_space_cfw = qdpack_operator_alloc(phase_space_wigner->m, phase_space_wigner->n);
    data.Cw = qdpack_matrix_alloc(phase_space_wigner->m, phase_space_wigner->n);

    // first calculate the characteristic function
    // lambda_max = 2*pi * alpha_max ???
    distribution_function_characteristic_w(data.Cw, a, ad, rho, data.lambda_max); 
    
  
//    data.Cw = phase_space_cfw;

    // set-up integration spaces
    data.w_re  = gsl_integration_workspace_alloc(qawo_limit);
    data.w_im  = gsl_integration_workspace_alloc(qawo_limit);

    data.f_re.function = &cb_qawo_wigner_f_re;
    data.f_re.params   = &data;
    data.f_im.function = &cb_qawo_wigner_f_im; 
    data.f_im.params   = &data;

    data.wf_im = gsl_integration_qawo_table_alloc(2 * QDPACK_IMAG(data.alpha), 2 * QDPACK_IMAG(data.lambda_max), GSL_INTEG_COSINE, qawo_n); 
    gsl_integration_qawo_table_set_length(data.wf_im, 2 * QDPACK_IMAG(data.lambda_max));
    data.wf_re = gsl_integration_qawo_table_alloc(2 * QDPACK_IMAG(data.alpha), 2 * QDPACK_REAL(data.lambda_max), GSL_INTEG_COSINE, qawo_n); 
    gsl_integration_qawo_table_set_length(data.wf_re, 2 * QDPACK_REAL(data.lambda_max));

    // i - real
    // j - imag
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            QDPACK_SET_COMPLEX(&(data.alpha), 
                coordinate_from_index(QDPACK_REAL(alpha_max), i, N),
                coordinate_from_index(QDPACK_IMAG(alpha_max), j, N));
       
//            data.wf_re = gsl_integration_qawo_table_alloc(
//                             2 * QDPACK_IMAG(alpha), 
//                             2 * QDPACK_REAL(lambda_max), 
//                             GSL_INTEG_COSINE, qawo_n); 

              gsl_integration_qawo_table_set(data.wf_re, 2 * QDPACK_IMAG(data.alpha), 2 * QDPACK_REAL(data.lambda_max), GSL_INTEG_COSINE);
//            gsl_integration_qawo_table_set_length(data.wf_re, 2 * QDPACK_REAL(lambda_max));

            if ((status = gsl_integration_qawo(&(data.f_re), -QDPACK_REAL(data.lambda_max), abserr, relerr, qawo_limit, data.w_re, data.wf_re, &result, &abserr)) != 0)
            {
                printf("Warning: gsl_integration_qawo failed [status = %d]...\n", status);
            }
            else
            {
                printf("debug: alpha = [%f, %f]: gsl_integration_qawo converged to %f with abserr %f\n", 
                        QDPACK_REAL(data.alpha), QDPACK_IMAG(data.alpha), result, abserr);
            }
            
            //result = ch_td_wigner_interpolate(data.alpha, (void *)&data);

            QDPACK_SET_COMPLEX(&z, result, 0.0);
            //z = qdpack_operator_get(data.Cw, i, j);

            qdpack_matrix_set(phase_space_wigner, i, j, z);
            printf("debug: W(alpha = [%f, %f]) = [%f, %f]\n", 
                   QDPACK_REAL(data.alpha), QDPACK_IMAG(data.alpha), QDPACK_REAL(z), QDPACK_IMAG(z));

//            gsl_integration_qawo_table_free(data.wf_re);
        }
    }    

    gsl_integration_workspace_free(data.w_re);
    gsl_integration_workspace_free(data.w_im);

    return 0;
}

/*
 * Compute the Wigner function corresponding to an harmonic
 * oscillator density matrix, by 2d FFT of the characteristic function.
 */
int
distribution_function_wigner_cf_fft(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max)
{
    qdpack_matrix_t *phase_space_cf;
    qdpack_complex lambda_max;

    phase_space_cf = qdpack_matrix_alloc(phase_space_wigner->m, phase_space_wigner->n);

    // first calculate the characteristic function
    // lambda_max = 2*pi * alpha_max ???
    distribution_function_characteristic_w(phase_space_cf, a, ad, rho, lambda_max); 
    
    //gsl_ifft(phase_space_wigner, phase_space_cf);

    return 0;
}

/*
 * Compute the Wigner function corresponding to an harmonic
 * oscillator density matrix, by 2d numberical integration
 * of the characteristic function.
 */
static double
cb_qags_wigner_f_im(double lambda_im, void *rawdata)
{
    double z;
    wigner_data_t *data = (wigner_data_t *)rawdata;

    QDPACK_SET_COMPLEX(&(data->lambda), QDPACK_REAL(data->lambda), lambda_im);
    //z = exp(-qdpack_complex_abs2(data->alpha)*qdpack_complex_abs2(data->lambda));
    z = ch_td_wigner_interpolate(data->lambda, data);
    return z;
}

static double
cb_qags_wigner_f_re(double lambda_re, void *rawdata)
{
    double result, error;
    int status;

    wigner_data_t *data = (wigner_data_t *)rawdata;

    QDPACK_SET_COMPLEX(&(data->lambda), lambda_re, 0);

    if ((status = gsl_integration_qags(
                    &(data->f_im), 
                    -QDPACK_IMAG(data->lambda_max),
                    +QDPACK_IMAG(data->lambda_max), 
                    abserr, relerr, 
                    10000, data->w_im, 
                    &result, &error)) != 0)
    {
        printf("Warning: gsl_integration_qawo failed [status = %d]...\n", status);
    }
    else
    {
     //  printf("debug [f_re(%f)]: gsl_integration_qawo converged to %f with abserr %f [interval = %d]\n", lambda_re, result, error, (int)(data->w_im->size));
    }        

    return result;
//    return exp(-qdpack_complex_abs2(data->alpha)*qdpack_complex_abs2(data->lambda));
}

int
distribution_function_wigner_cf_quad(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max)
{
    //qdpack_operator_t *phase_space_cfw;
    int i, j, status;
    double result, error;
    qdpack_complex z;

    wigner_data_t data;

    data.N = phase_space_wigner->m;
    data.alpha_max  = alpha_max;
    data.lambda_max = alpha_max;

    // first calculate the characteristic function
    data.Cw = qdpack_matrix_alloc(phase_space_wigner->m, phase_space_wigner->n);
    distribution_function_characteristic_w(data.Cw, a, ad, rho, data.lambda_max); 

    // set-up integration spaces
    data.w_re = gsl_integration_workspace_alloc(10000);
    data.w_im = gsl_integration_workspace_alloc(10000);
       
    data.f_re.function = &cb_qags_wigner_f_re;
    data.f_re.params   = &data;
    data.f_im.function = &cb_qags_wigner_f_im; 
    data.f_im.params   = &data;

    // i - real
    // j - imag
    for (i = 0; i < data.N; i++)
    {
        for (j = 0; j < data.N; j++)
        {
            QDPACK_SET_COMPLEX(&(data.alpha), 
                coordinate_from_index(QDPACK_REAL(alpha_max), i, data.N),
                coordinate_from_index(QDPACK_IMAG(alpha_max), j, data.N));
       
            if ((status = gsl_integration_qags(
                                &(data.f_re), 
                                -QDPACK_REAL(data.lambda_max),
                                +QDPACK_REAL(data.lambda_max), 
                                abserr, relerr, 
                                10000, data.w_re, 
                                &result, &error)) != 0)
            {
                printf("Warning: gsl_integration_qawo failed [status = %d]...\n", status);
            }
            else
            {
                printf("debug: alpha = [%f, %f]: gsl_integration_qawo converged to %f with abserr %f [interval = %d]\n", 
                        QDPACK_REAL(data.alpha), QDPACK_IMAG(data.alpha), result, error, (int)(data.w_re->size));
            }
            
            QDPACK_SET_COMPLEX(&z, result, 0.0);
            qdpack_matrix_set(phase_space_wigner, i, j, z);
        }
    }    

    gsl_integration_workspace_free(data.w_re);
    gsl_integration_workspace_free(data.w_im);

    return 0;
}




/*
 * Integrand for quadrature:
 * XXX: unfinished
 */
double
cb_integrand_wigner_hermite(double x, void *data)
{
    //wigner_data_t *integrand_data =  (wigner_data_t *)data;
    
//    double x = 0;
    int n, N = 0;

    for (n = 0; n < N; n++)
    {
//        x += QDPACK_REAL(qdpack_operator_get(integrand_data->rho, n, n)) *
//             hermite_poly_eval(n, x + integrand_data->q / 2) *
//             hermite_poly_eval(n, x - integrand_data->q / 2);
    }

    return x;
}

/*
 * XXX unfinished
 */
int
distribution_function_wigner_hermite_quad(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *rho,
            double x_limit)
{
    //qdpack_matrix_t *phase_space_cf;
    int i, j, N, M;
    qdpack_complex z;
    //double q, p; 

    //wigner_data_t integrand_data;
   
    N = phase_space_wigner->m;
    M = phase_space_wigner->n;

    // i - real
    // j - imag
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < M; j++)
        {
            //integrand_data.q = (i - N/2) * x_limit / (N/2);
            //integrand_data.p = (j - M/2) * x_limit / (M/2);
            
    //        z = gsl_quad_1d(-x_limit, x_limit,
    //                        cb_integrand_wigner_hermite,
    //                        (void *)&integrand_data);

            qdpack_matrix_set(phase_space_wigner, i, j, z);
        }
    }    

    return 0;
}



int
distribution_function_wigner_cf_sum(
            qdpack_matrix_t *phase_space_wigner, 
            qdpack_operator_t *a,
            qdpack_operator_t *ad,
            qdpack_operator_t *rho,
            qdpack_complex alpha_max)
{
    //qdpack_matrix_t *phase_space_cfw;
    int i, j, m, n;
    double dlambda_im, dlambda_re;
    qdpack_complex z, ez, sum_outer, sum_inner;

    wigner_data_t data;

    data.N = phase_space_wigner->m;
    data.alpha_max  = alpha_max;
    data.lambda_max = alpha_max;

    printf("rho = [%d, %d], a = [%d, %d], ad = [%d, %d]\n",
            rho->m, rho->n,
            a->m, a->n,
            ad->m, ad->n);


    // first calculate the characteristic function
    data.Cw = qdpack_matrix_alloc(phase_space_wigner->m, phase_space_wigner->n);
    distribution_function_characteristic_w(data.Cw, a, ad, rho, data.lambda_max); 

    dlambda_re = 2 * QDPACK_REAL(data.lambda_max) / data.N;
    dlambda_im = 2 * QDPACK_IMAG(data.lambda_max) / data.N;

    // i - real
    // j - imag
    for (i = 0; i < data.N; i++)
    {
        for (j = 0; j < data.N; j++)
        {
            QDPACK_SET_COMPLEX(&(data.alpha), 
                coordinate_from_index(QDPACK_REAL(alpha_max), i, data.N),
                coordinate_from_index(QDPACK_IMAG(alpha_max), j, data.N));

            //result = 0.0;
            QDPACK_SET_COMPLEX(&sum_outer, 0.0, 0.0);

            for (m = 0; m < data.N; m++)
            {
                //result_inner = 0.0;
                QDPACK_SET_COMPLEX(&sum_inner, 0.0, 0.0);

                for (n = 0; n < data.N; n++)
                {
                    QDPACK_SET_COMPLEX(&(data.lambda), 
                        coordinate_from_index(QDPACK_REAL(data.lambda_max), m, data.N),
                        coordinate_from_index(QDPACK_IMAG(data.lambda_max), n, data.N));

                    // value of the characteristic function
                    z = qdpack_matrix_get(data.Cw, m, n);

//                    result_inner += cos(2*QDPACK_IMAG(data.lambda)*QDPACK_REAL(data.alpha)) * QDPACK_REAL(z) * dlambda_im;

                    ez = qdpack_complex_mul_real(
                            qdpack_complex_exp(
                               qdpack_complex_sub(
                                    qdpack_complex_mul(qdpack_complex_conjugate(data.lambda), data.alpha),
                                    qdpack_complex_mul(data.lambda, qdpack_complex_conjugate(data.alpha))
                                              )
                                           ), dlambda_im);

//                    ez = qdpack_complex_exp(la);
//                    QDPACK_SET_COMPLEX(&ez, cos(QDPACK_REAL(la)), sin(QDPACK_IMAG(la));
                   
                 
                    sum_inner = qdpack_complex_add(sum_inner, qdpack_complex_mul(ez, z));                

                    // including sin * sin
//                    result_inner += (cos(2*QDPACK_REAL(data.lambda)*QDPACK_IMAG(data.alpha)) * cos(2*QDPACK_IMAG(data.lambda)*QDPACK_REAL(data.alpha)) + 
  //                                   sin(2*QDPACK_REAL(data.lambda)*QDPACK_IMAG(data.alpha)) * sin(2*QDPACK_IMAG(data.lambda)*QDPACK_REAL(data.alpha))
    //                                ) * QDPACK_REAL(z) * dlambda_im;
//

                }
//                result += cos(2*QDPACK_REAL(data.lambda)*QDPACK_IMAG(data.alpha)) * result_inner  * dlambda_re;
//                sum_inner = qdpack_complex_mul_real(result_inner, dlambda_re);

                sum_outer = qdpack_complex_add(sum_outer, qdpack_complex_mul_real(sum_inner, dlambda_re));
//                result += result_inner  * dlambda_re;
            }
            //result = result / (M_PI * M_PI);
            sum_outer = qdpack_complex_mul_real(sum_outer, 1.0/(M_PI * M_PI));

//            QDPACK_SET_COMPLEX(&z, result, 0.0);
            qdpack_matrix_set(phase_space_wigner, i, j, sum_outer);
        }
    }    

    return 0;
}


