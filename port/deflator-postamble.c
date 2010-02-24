#include <clover.h>

#if defined(HAVE_LAPACK)
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#else
#  error "no linear algebra library"
#endif

int
q(df_postamble)(
        struct Q(State)         *s,
        struct Q(Deflator)      *d)
{
    assert(NULL != s &&
            NULL != d);
    if (NULL == s ||
            NULL == d ||
            d->frozen || 
            d->umax <= d->usize ||
            d->vsize < d->nev)
        return 0;


    double v_norm2_min = eps_reortho * eps_reortho;
    int unew = 0,
        i_v = 0;
    long int usize_old = d->usize;
    while ((d->usize < d->umax) && (i_v < d->nev)) {
        /* reorthogonalize n_reortho times */
/* macros to reuse workspace */
#define cur_v       (d->work_c_1)
#define cur_Av      (d->work_c_2)
#define cur_z_v     (d->work_z_1)
#define cur_z_Av    (d->work_z_2)
#define cur_z_aux   (d->work_z_3)
        latmat_c_get_col(d->V, i_v, cur_v);
        
        if (0 < d->usize) {
            for (int i_reortho = n_reortho; i_reortho--; ) {
                latmat_c cur_U = latmat_c_submat_col(d->U, 0, d->usize);
                lat_lmH_dot_lv(d->usize, 
                               cur_U, 
                               cur_v, 
                               d->zwork);
                lat_lm_dot_zv(d->usize, 
                              cur_U, 
                              d->zwork, 
                              cur_Av);
                /* FIXME cur_v <- cur_v - cur_Av surely may be optimized */
                doublecomplex tmone = {-1.0, 0.0 };
                lat_cc_axpy(tmone, cur_Av, cur_v);
            }
        }
        double v_norm2 = lat_c_nrm2(cur_v);
        if (v_norm2 < v_norm2_min) {    /* skip */
            i_v++;
            continue;
        }
        lat_c_scal_d(1 / sqrt(v_norm2), cur_v);
        latmat_c_insert_col(d->U, d->usize, cur_v);

        /* compute A U and store in V ; 
           XXX requires(?) double precision operator */
        latvec_cz_copy(cur_v, cur_z_v);
        /* TODO apply double-prec operator to cur_z_v: */
        latvec_z_linop(s, d, cur_z_Av, cur_z_v, cur_z_aux);
        /* FIXME need to implement latmat_c^H (dot) latvec_z to make efficient */
        for (int i = 0; i < d->usize; i++) {
            d->H[i + d->usize*d->umax] = lat_cz_dotu(d->U.store[i], cur_z_Av);
            d->H[d->usize + i*d->umax].r = d->H[i + d->usize*d->umax].r;
            d->H[d->usize + i*d->umax].i =-d->H[i + d->usize*d->umax].i;
        }
        d->H[d->usize * (d->umax + 1)].r = lat_zz_dotu(cur_z_v, cur_z_Av);
        d->H[d->usize * (d->umax + 1)].i = 0.0;

        d->usize ++;
        unew ++;
        i_v ++;
#undef cur_v
#undef cur_Av
#undef cur_z_v
#undef cur_z_Av
#undef cur_z_aux
    }
    assert(usize_old + unew == d->usize);

    /* compute Cholesky decomposition */
    memcpy(d->C, d->H, d->usize * d->umax * sizeof(d->C[0]));

#if HAVE_LAPACK
    long int usize  = d->usize,
             umax   = d->umax,
             info   = 0;
    char cU = 'U';
    zpotrf_(&cU, &usize, d->C, &umax, &info, 1);
    assert(0 == info);
#elif HAVE_GSL
    gsl_matrix_complex_view gsl_C = gsl_matrix_complex_view_array_with_tda(
            d->C, usize, usize, umax);
    CHECK_GSL_STATUS(gsl_matrix_transpose(&gsl_C.matrix));
    CHECK_GSL_STATUS(gsl_linalg_complex_cholesky_decomp(&gsl_C.matrix));
#else
#  error "no linear algebra library"
#endif

    return unew;
}
