#include <clover.h>

#include <math.h>
#if defined(HAVE_LAPACK)
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#else
#  error "no linear algebra library"
#endif

int
q(df_solve_in_eigenspace)(
        struct Q(State)         *s,
        struct Q(Deflator)      *d,
        struct Fermion          *x,
        struct Fermion          *b)
{
    assert(NULL != s &&
            NULL != d &&
            NULL != x &&
            NULL != b);
    latvec_c lv_x   = latvec_c_view(d->dim, x);
    if (d->usize <= 0) {
        latvec_c_zero(lv_x);
        return 0;
    }    
    latvec_c lv_b   = latvec_c_view(d->dim, b);
    latmat_c cur_U  = latmat_c_submat_col(d->U, 0, d->usize);
    lat_lmH_dot_lv(d->usize, cur_U, lv_b, d->zwork);

#if defined(HAVE_LAPACK)
    long int usize  = d->usize,
             ONE    = 1,
             umax   = d->umax,
             info   = 0;
    char cU = 'U';
    zpotrs_(&cU, &usize, &ONE, d->C, &umax, d->zwork, &umax, &info, 1);
    assert(0 == info);
#elif defined(HAVE_GSL)

    gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            d->C, d->usize, d->usize, d->umax);
    gsl_vector_complex_view gsl_pb = gsl_vector_complex_view_array(
            d->zwork, d->usize);
    gsl_vector_complex_view gsl_px = gsl_vector_complex_view_array(
            d->zwork2, d->usize);
    CHECK_GSL_STATUS(gsl_linalg_complex_cholesky_solve(
            &gsl_C.matrix, 
            &gsl_pb.vector, 
            &gsl_px.vector));
#else
#  error "no linear algebra library"
#endif

    lat_lm_dot_zv(d->usize, cur_U, d->zwork, lv_x);
    return 0;
}

int
q(df_preamble)(
        struct Q(State)         *s,
        struct Q(Deflator)      *d,
        struct Fermion          *x,
        struct Fermion          *b)
{
    assert(NULL != s &&
            NULL != d &&
            NULL != x &&
            NULL != b);

    latvec_c lv_x   = latvec_c_view(d->dim, x);
    latvec_c lv_b   = latvec_c_view(d->dim, b);

#define cur_r       (d->work_c_1)
#define cur_r_aux   (d->work_c_2)
    if (d->usize <= 0) {
        latvec_c_zero(lv_x);
        latvec_c_copy(lv_b, cur_r);
    } else {
        if (q(df_solve_in_eigenspace)(s, d, x, b))
            return 1;
        /* compute residual */
        latvec_c_linop(cur_r, lv_x, cur_r_aux);
        /* FIXME optimize the code below with a special primitive;
           cur_r <- lv_b - cur_r */
        lat_cc_axpy_d(-1., lv_b, cur_r);
        lat_c_scal_d(-1., cur_r);   
    }
    double rnorm = sqrt(lat_c_nrm2(cur_r));
    lat_c_scal_d(1. / rnorm, cur_r);
    
    /* save normalized residual as the first vector */
    latmat_c_insert_col(d->V, 0, cur_r);
    
    /* init eigenvalue search: vsize, T, V */
    d->vsize = 0;
    memset(d->T, 0, d->vmax * d->vmax * sizeof(d->T[0]));

    /* init stopping threshold */
    d->resid_norm_sq_min = d->eps * d->eps * rnorm * rnorm;
    
    return 0;
}
