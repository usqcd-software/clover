#include <clover.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#if defined(HAVE_LAPACK)
#  include <blas.h>
#  include <lapack.h>
#elif defined(HAVE_GSL)
#  include <gsl/gsl_linalg.h>
#  include <gsl/gsl_eigen.h>
#  include <gsl/gsl_blas.h>
#else
#  error "no linear algebra library"
#endif

int
q(df_update1)(
        struct Q(State)         *s,
        struct Q(Deflator)      *d,
        double alpha, 
        double beta, 
        double alpha_prev, 
        double beta_prev, 
        double resid_norm_sq, 
        struct Fermion *resid,
        struct Fermion *A_resid
        )
{
    if (NULL == d || 
            NULL == s ||
            d->frozen)
        return 0;       /* skip call */

    if (d->vsize < d->vmax) 
        return -1;      /* should call update0 instead */

    assert(d->vmax == d->vsize && 
            1 < d->vsize);

    /* [vsize-1, vsize-1] elem */
    doublecomplex *pT = d->T + (d->vsize - 1 ) * (1 + d->vmax);
    pT->r = 1. / alpha + beta_prev / alpha_prev;
    pT->i = 0.0;


    /* do restart */
    d->vsize = 2 * d->nev;
    long int vmax = d->vmax;
    long int vsize = d->vsize;
#if defined(HAVE_LAPACK)
    char cV = 'V',
         cU = 'U',
         cL = 'L',
         cR = 'R',
         cC = 'C',
         cN = 'N';
    long int info = 0;
    long int tmp_i;


    /* diagonalize T:vmax*vmax matrix 
       requires zwork size = lwork >= 2*vmax-1, rwork size >=3*vmax-2 */
    memcpy(d->hevecs1, d->T, vmax * vmax * sizeof(d->T[0]));
    zheev_(&cV, &cU, &vmax, d->hevecs1, &vmax, 
           d->hevals, d->zwork, &(d->lwork), 
           d->rwork, &info, 1, 1);
    PRINT_STATUS(info, 0);
    assert(0 == info);

    /* diagonalize T:(vmax-1)*(vmax-1) matrix
       requires zwork size = lwork >= 2*vmax-1 */
    memcpy(d->hevecs2, d->T, vmax * vmax * sizeof(d->T[0]));
    tmp_i = vmax - 1;
    zheev_(&cV, &cU, &tmp_i, d->hevecs2, &vmax,
           d->hevals, d->zwork, &(d->lwork),
           d->rwork, &info, 1, 1);
    PRINT_STATUS(info, 0);
    assert(0 == info);

    /* select first nev vectors from both spaces */
    memcpy(d->hevecs1 + d->nev * vmax, d->hevecs2, 
           d->nev * vmax * sizeof(d->hevecs1[0]));
    /* fill [vmax-1, nev:2*nev] with zeros */
    for (int i = d->nev; i < 2 * d->nev; i++) {
        doublecomplex *p = d->hevecs1 + (i + 1) * vmax - 1;
        p->r = 0.0;
        p->i = 0.0;
    }

    /* QR factorization/orthogonalization of 2*nev columns
       requires zwork size = lwork >= 2*vmax, tau size >= 2*nev */
    zgeqrf_(&vmax, &vsize, d->hevecs1, &vmax, 
            d->tau, d->zwork, &(d->lwork),
            &info);
    assert(0 == info);

    /* compute hevecs2 <- Q^H T Q */
    memcpy(d->hevecs2, d->T, vmax * vmax * sizeof(d->T[0]));
    /* requires size of zwork = lwork >= vmax ; optimal-?*/
    
    zunmqr_(&cL, &cC, &vmax, &vmax, &vsize, 
            d->hevecs1, &vmax, d->tau,
            d->hevecs2, &vmax, d->zwork,
            &(d->lwork), &info, 1, 1);
    assert(0 == info);
    zunmqr_(&cR, &cN, &vsize, &vmax, &vsize, 
            d->hevecs1, &vmax, d->tau,
            d->hevecs2, &vmax, d->zwork,
            &(d->lwork), &info, 1, 1);
    assert(0 == info);

    /* compute eigenpairs of Q^H T Q = Z M Z^H, hevecs2 <- Z */
    zheev_(&cV, &cU, &vsize, d->hevecs2, &vmax, 
           d->hevals, d->zwork, &(d->lwork),
           d->rwork, &info, 1, 1);
    assert(0 == info);

    /* fill Z[2nev:vmax, 0:2nev] with zeros */
    for (int j = 0; j < vsize; j++) {
        doublecomplex *p = d->hevecs2 + j * vmax;
        for (int i = vsize ; i < vmax; i++)
            p[i].r = p[i].i = 0.0;
    }
    /* compute hevecs2 <- Q Z */
    zunmqr_(&cL, &cN, &vmax, &vsize, &vsize,
            d->hevecs1, &vmax, d->tau,
            d->hevecs2, &vmax, d->zwork,
            &(d->lwork), &info, 1, 1);
    assert(0 == info);
    

#elif defined(HAVE_GSL)

    gsl_complex gsl_tzero = gsl_complex_rect(0., 0.);
    gsl_complex gsl_tpone = gsl_complex_rect(1., 0.);
    gsl_complex gsl_tmone = gsl_complex_rect(-1., 0.);

    gsl_matrix_complex_view gsl_T = gsl_matrix_complex_view_array(
            d->T, 
            vmax, vmax);
    gsl_matrix_complex_transpose(&gsl_T.matrix);    /* d->T uses FORTRAN matrix indexing */
    /* eigenpairs of T */
    gsl_matrix_complex_memcpy(d->gsl_T_full, &gsl_T.matrix);
    CHECK_GSL_STATUS(gsl_eigen_hermv(
            d->gsl_T_full, 
            d->gsl_hevals1, 
            d->gsl_hevecs1,
            d->gsl_wkspace1))
    CHECK_GSL_STATUS(gsl_sort_smallest_index(d->hevals_select1, d->nev, 
            gsl_vector_const_ptr(d->gsl_hevals1, 0), 1, vmax))

    /* eigenpairs of T[:-1, :-1] */
    gsl_matrix_complex_view gsl_T_m1 = gsl_matrix_complex_submatrix(
            &gsl_T.matrix, 
            0, 0, vmax-1, vmax-1);
    CHECK_GSL_STATUS(gsl_matrix_complex_memcpy(
                d->gsl_T_m1, 
                &gsl_T_m1.matrix));

    CHECK_GSL_STATUS(gsl_eigen_hermv(
            &gsl_T_m1.matrix, 
            d->gsl_hevals2,
            d->gsl_hevecs2,
            d->gsl_wkspace2))
    CHECK_GSL_STATUS(gsl_sort_smallest_index(d->hevals_select2, d->nev,
            gsl_vector_const_ptr(d->gsl_hevals2, 0), 1, vmax-1));

    /* construct Q = (Y[:nev], Y1[:nev]) */
    for (int j = 0; j < d->nev; j++) {
        int j1 = d->hevals_select1[j];
        for (int i = 0; i < vmax; i++)
            gsl_matrix_complex_set(d->gsl_QR, i, j,
                    gsl_matrix_complex_get(d->gsl_hevecs1, i, j1));
    }
    for (int j = 0; j < d->nev; j++) {
        int j2 = d->hevals_select2[j];
        for (int i = 0; i < vmax-1; i++)
            gsl_matrix_complex_set(d->gsl_QR, i, d->nev + j,
                    gsl_matrix_complex_get(d->gsl_hevecs2, i, j2));
        gsl_matrix_complex_set(d->gsl_QR, vmax-1, d->nev + j,
                gsl_complex_rect(0., 0.));
    }

    /* QR decomp */
    CHECK_GSL_STATUS(gsl_linalg_complex_QR_decomp(
            d->gsl_QR, d->gsl_tau));
    CHECK_GSL_STATUS(gsl_linalg_complex_QR_unpack(
            d->gsl_QR, d->gsl_tau, 
            d->gsl_Q_unpack, d->gsl_tmp_MxS));
    gsl_matrix_complex_view gsl_Q = gsl_matrix_complex_submatrix(
            d->gsl_Q_unpack,
            0, 0, vsize, vmax);

    /* projecting matrix T on Q: Q^H . T . Q */
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
            gsl_tpone, &gsl_T.matrix, gsl_Q.matrix,
            gsl_tzero, d->gsl_tmp_MxS));
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasConjTrans, CblasNoTrans,
            gsl_tpone, &gsl_Q.matrix, d->gsl_tmp_MxS,
            gsl_tzero, d->gsl_T_proj));

    /* eigenpairs of Q^H . T . Q */
    gsl_vector_complex_view gsl_hevals = gsl_vector_complex_view_array(
            d->hevals, vsize);
    CHECK_GSL_STATUS(gsl_eigen_hermv(
            d->gsl_T_proj,
            &gsl_hevals.vector,
            d->gsl_hevecs3,
            d->gsl_wkspace3));
    CHECK_GSL_STATUS(gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,
            gsl_tpone, &gsl_Q.matrix, d->gsl_hevecs3,
            gsl_tzero, d->gsl_tmp_MxS));

    gsl_matrix_complex_view gsl_QZ_transp = gsl_matrix_complex_view_array_with_tda(
            d->hevecs2, vsize, vmax, vmax);
    CHECK_GSL_STATUS(gsl_matrix_complex_transpose_memcpy( /* back to FORTRAN */
                &gsl_QZ_transp.matrix, d->gsl_tmp_MxS));


    /* TODO thorough review of the code above */NOT_IMPLEMENTED;
#else
#  error "no linear algebra library"
#endif


    /* rotate V[:, 0:vmax] space with (Q Z) */
    latmat_c tmp_V = latmat_c_submat_col(d->tmp_V, 0, vsize);
    lat_lm_dot_zm(vsize, vmax, 
                  d->V,
                  d->hevecs2, vmax, 
                  tmp_V);
    latmat_c_copy(tmp_V, latmat_c_submat_col(d->V, 0, vsize));
    
    /* check eig convergence */
    if (resid_norm_sq < d->resid_norm_sq_min)
        return 3;   /* eig converged */
    double resid_norm = sqrt(resid_norm_sq);
    
    /* compute new T */
    memset(d->T, 0, vmax * vmax * sizeof(d->T[0]));
    lat_lmH_dot_lv(vsize, 
                   tmp_V, 
                   latvec_c_view(d->dim, Ar), 
                   d->T + vsize * vmax);
    for (int i = 0 ; i < vsize ; i++) {
        d->T[i * (vmax + 1)].r      = d->hevals[i];
        d->T[i * (vmax + 1)].i      = 0.0;

        d->T[i + vsize * vmax].r    /= resid_norm;
        d->T[i + vsize * vmax].i    /= resid_norm;

        d->T[vsize + i * vmax].r    =  d->T[i + vsize * vmax].r;
        d->T[vsize + i * vmax].i    = -d->T[i + vsize * vmax].i;
    }


    /* remember the vector ||resid|| */
#define cur_r   (d->work_c_1)
    latvec_c_copy(latvec_c_view(d->dim, resid), 
                  cur_r);
    lat_c_scal_d(1. / resid_norm, cur_r);
    latmat_c_insert_col(d->V, d->vsize, cur_r)

    d->vsize += 1;
    
    return 0; /* normal */
}
