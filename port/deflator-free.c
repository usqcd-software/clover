#include <clover.h>

void
q(df_free)(struct Q(Deflator) **deflator_ptr)
{
    struct Q(State) *s;

    if (deflator_ptr == 0 || *deflator_ptr == 0)
        return;
    s = (*deflator_ptr)->state;
    BEGIN_TIMING(s);
    
    /* XXX free other components of the deflator */
#define guarded_free(v, cmd) { if (NULL == v) cmd; }
    if (!latmat_c_is_null(&(d->V))) latmat_c_free(s, d->V);
    guarded_free(d->T,          q(free)(s, d->T));

    if (!latmat_c_is_null(&(d->U))) latmat_c_free(s, d->U);
    guarded_free(d->H,          q(free)(s, d->H));
    guarded_free(d->C,          q(free)(s, d->C));

    guarded_free(d->hevecs2,    q(free)(s, d->hevecs2));
    guarded_free(d->hevals,     q(free)(s, d->hevals));
    guarded_free(d->zwork,      q(free)(s, d->zwork));

#if defined(HAVE_LAPACK)
    guarded_free(d->hevecs1,    q(free)(s, d->hevecs1));
    guarded_free(d->tau,        q(free)(s, d->tau));
    guarded_free(d->rwork,      q(free)(s, d->rwork));
#elif defined(HAVE_GSL)
    guarded_free(d->zwork2,     q(free)(s, d->zwork2));
    guarded_free(d->gsl_T_full,     gsl_matrix_complex_free(d->gsl_T_full));
    guarded_free(d->gsl_hevecs1,    gsl_matrix_complex_free(d->gsl_hevecs1));
    guarded_free(d->gsl_hevals1,    gsl_vector_free(d->gsl_hevals1));
    guarded_free(d->gsl_wkspace1,   gsl_eigen_herm_free(d->gsl_wkspace1));
    guarded_free(d->gsl_T_m1,       gsl_matrix_complex_free(d->gsl_T_m1));
    guarded_free(d->gsl_hevecs2,    gsl_matrix_complex_free(d->gsl_hevecs2);
    guarded_free(d->gsl_hevals2,    gsl_vector_free(d->gsl_hevals2));
    guarded_free(d->gsl_wkspace2,   gsl_eigen_herm_free(d->gsl_wkspace2));
    guarded_free(d->gsl_T_proj,     gsl_matrix_complex_free(d->gsl_T_proj));
    guarded_free(d->gsl_hevecs3,    gsl_matrix_complex_free(d->gsl_hevecs3));
    guarded_free(d->gsl_wkspace3,   gsl_eigen_herm_free(d->gsl_wkspace3));
    guarded_free(d->gsl_QR,         gsl_matrix_complex_free(d->gsl_QR));
    guarded_free(d->gsl_Q_unpack,   gsl_matrix_complex_free(d->gsl_Q_unpack));
    guarded_free(d->gsl_tmp_MxS,    gsl_matrix_complex_free(d->gsl_tmp_MxS));
    guarded_free(d->gsl_tau,        gsl_vector_complex_free(d->gsl_tau));
    guarded_free(d->hevals_select1, q(free)(s, d->hevals_select1);
    guarded_free(d->hevals_select2, q(free)(s, d->hevals_select2));
#  error "no linear algebra library"
#endif

    if (!latmat_c_is_null(&(d->tmp_V)))     latmat_c_free(s, tmp_V);
    if (!latvec_z_is_null(&(d->work_z_1)))  latvec_z_free(s, d->work_z_1);
    if (!latvec_z_is_null(&(d->work_z_2)))  latvec_z_free(s, d->work_z_2);
    if (!latvec_z_is_null(&(d->work_z_3)))  latvec_z_free(s, d->work_z_3);
    if (!latvec_c_is_null(&(d->work_c_1)))  latvec_c_free(s, d->work_c_1);
    if (!latvec_c_is_null(&(d->work_c_2)))  latvec_c_free(s, d->work_c_2));


    END_TIMING(s, 0, 0, 0);
    *deflator_ptr = 0;
}
