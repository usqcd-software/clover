#include <clover.h>

int
Q(create_deflator)(struct Q(Deflator) **deflator_ptr,
                   struct Q(State) *s,
                   int vmax, int nev, 
                   double eps, int umax)
{
    struct Q(Deflator) *d;
    int dim;

    if (s == NULL || s->error_latched)
        return 1;

    if (deflator_ptr == NULL)
        return q(set_error)(s, 0, "create_deflator(): NULL pointer");
  
    dim = s->even.full_size;
    *deflator_ptr = NULL;
    d = q(malloc)(s, sizeof (struct Q(Deflator)));
    if (d == 0)
        return q(set_error)(s, 0, "allocate_deflator(): not enough memory");

    /* check data types */
#if defined(HAVE_LAPACK)
    {
        doublecomplex dc;
        assert( &(dc.r) == (double *)(&dc) &&
                &(dc.i) == (double *)(&dc) + 1 &&
                sizeof(dc) == 2 * sizeof(double));
    }
#elif defined(HAVE_GSL)
    {
        doublecomplex dc;
        assert( &(dc.r) == (double *)(&dc) &&
                &(dc.i) == (double *)(&dc) + 1 &&
                sizeof(dc) == 2 * sizeof(double));
        gsl_complex *gc = (gsl_complex *)(&dc);
        assert( &(dc.r) == &(gc->dat[0]) &&
                &(dc.i) == &(gc->dat[1]) &&
                sizeof(dc) == sizeof(*gc));
    }
#else 
#  error "no linear algebra library"
#endif


    /* first, set all to null */
    q(latmat_c_null)(&(d->V));
    d->T                = NULL;

    q(latmat_c_null)(&(d->U));
    d->H                = NULL;
    d->C                = NULL;
    
    d->zwork            = NULL;
    d->hevals           = NULL;
    d->hevecs2          = NULL;

#if defined(HAVE_LAPACK)
    d->hevecs1          = NULL;
    d->tau              = NULL;
    d->rwork            = NULL;
#elif defined(HAVE_GSL)
    d->zwork2           = NULL;
    d->gsl_T_full       = NULL;
    d->gsl_hevecs1      = NULL;
    d->gsl_hevals1      = NULL;
    d->gsl_wkspace1     = NULL;
    d->gsl_T_m1         = NULL;
    d->gsl_hevecs2      = NULL;
    d->gsl_hevals2      = NULL;
    d->gsl_wkspace2     = NULL;
    d->gsl_T_proj       = NULL;
    d->gsl_hevecs3      = NULL;
    d->gsl_wkspace3     = NULL;
    d->gsl_QR           = NULL;
    d->gsl_Q_unpack     = NULL;
    d->gsl_tmp_MxS      = NULL;
    d->gsl_tau          = NULL;
    d->hevals_select1   = NULL;
    d->hevals_select2   = NULL;
#else
#  error "no linear algebra library"
#endif


    q(latmat_c_null)(&(d->tmp_V));
    q(latvec_c_null)(&(d->work_c_1));
    q(latvec_c_null)(&(d->work_c_2));
    q(latvec_c_null)(&(d->work_c_3));



    /* allocate */
#define ds      sizeof(double)
#define zs      sizeof(doublecomplex)
    d->V                = q(latmat_c_alloc)(s, dim, vmax);
    d->T                = q(malloc)(s, vmax * vmax * zs);

    d->U                = q(latmat_c_alloc)(s, dim, umax);
    d->H                = q(malloc)(s, umax * umax * zs);
    d->C                = q(malloc)(s, umax * umax * zs);

    d->hevecs2          = q(malloc)(s, vmax * vmax * zs);
    d->hevals           = q(malloc)(s, vmax * ds);
    /* reuse zwork in the  outer eigcg update */
    d->lwork            = (2*vmax > umax ? 2*vmax : umax);
    d->zwork            = q(malloc)(s, d->lwork * zs);

#if defined(HAVE_LAPACK)
    d->hevecs1          = q(malloc)(s, vmax * vmax * zs);
    d->tau              = q(malloc)(s, vmax * zs);
    d->rwork            = q(malloc)(s, 3 * vmax * ds);
#elif defined(HAVE_GSL)
    d->zwork2           = q(malloc)(s, umax * zs);
    d->gsl_T_full       = gsl_matrix_complex_alloc(vmax, vmax);
    d->gsl_hevecs1      = gsl_matrix_complex_alloc(vmax, vmax);
    d->gsl_hevals1      = gsl_vector_alloc(vmax);
    d->gsl_wkspace1     = gsl_eigen_hermv_alloc(vmax);
    d->gsl_T_m1         = gsl_matrix_complex_alloc(vmax-1, vmax-1);
    d->gsl_hevecs2      = gsl_matrix_complex_alloc(vmax-1, vmax-1);
    d->gsl_hevals2      = gsl_vector_alloc(vmax-1);
    d->gsl_wkspace2     = gsl_eigen_hermv_alloc(vmax-1);
    d->gsl_T_proj       = gsl_matrix_complex_alloc(2*nev, 2*nev);
    d->gsl_hevecs3      = gsl_matrix_complex_alloc(2*nev, 2*nev);
    d->gsl_wkspace3     = gsl_eigen_hermv_alloc(2*nev);
    d->gsl_QR           = gsl_matrix_complex_alloc(vmax, 2*nev);
    d->gsl_Q_unpack     = gsl_matrix_complex_alloc(vmax, vmax);
    d->gsl_tmp_MxS      = gsl_matrix_complex_alloc(vmax, 2*nev);
    d->gsl_tau          = gsl_vector_complex_alloc(2*nev);
    d->hevals_select1   = q(malloc)(s, vmax * sizeof(d->hevals_select1[0]));
    d->hevals_select2   = q(malloc)(s, vmax * sizeof(d->hevals_select2[0]));
#else
#  error "no linear algebra library"
#endif

    d->tmp_V            = q(latmat_c_alloc)(s, dim, 2*nev);
    d->work_c_1         = q(latvec_c_alloc)(s, dim);
    d->work_c_2         = q(latvec_c_alloc)(s, dim);
    d->work_c_3         = q(latvec_c_alloc)(s, dim);

    /* check allocation */
    if (
            q(latmat_c_is_null)(&(d->V))           ||
            NULL == d->T                        ||

            q(latmat_c_is_null)(&(d->U))           ||
            NULL == d->H                        ||
            NULL == d->C                        ||
                    
            NULL == d->zwork                    ||
            NULL == d->hevals                   ||
            NULL == d->hevecs2                  ||
                                               
#if defined(HAVE_LAPACK)                       
            NULL == d->hevecs1                  ||
            NULL == d->tau                      ||
            NULL == d->rwork                    ||
#elif defined(HAVE_GSL)                        
            NULL == d->zwork2                   ||
            NULL == d->gsl_T_full               ||
            NULL == d->gsl_hevecs1              ||
            NULL == d->gsl_hevals1              ||
            NULL == d->gsl_wkspace1             ||
            NULL == d->gsl_T_m1                 ||
            NULL == d->gsl_hevecs2              ||
            NULL == d->gsl_hevals2              ||
            NULL == d->gsl_wkspace2             ||  
            NULL == d->gsl_T_proj               ||
            NULL == d->gsl_hevecs3              ||
            NULL == d->gsl_wkspace3             ||
            NULL == d->gsl_QR                   ||
            NULL == d->gsl_Q_unpack             ||
            NULL == d->gsl_tmp_MxS              ||
            NULL == d->gsl_tau                  ||
            NULL == d->hevals_select1           ||
            NULL == d->hevals_select2           ||
#else
#  error "no linear algebra library"
#endif
            q(latmat_c_is_null)(&(d->tmp_V))       ||
            q(latvec_c_is_null)(&(d->work_c_1))    ||
            q(latvec_c_is_null)(&(d->work_c_2))    ||
            q(latvec_c_is_null)(&(d->work_c_3))
    ) {
        Q(free_deflator)(&d);
        return q(set_error)(s, 0, "allocate_deflator(): not enough memory");
    }

    BEGIN_TIMING(s);
    d->state = s;

    d->dim   = dim;
    d->vmax  = vmax;
    d->vsize  = 0;
    d->nev   = nev;
    d->eps   = eps;
    d->umax  = umax;
    d->usize  = 0;
    d->frozen = 0;
    

    *deflator_ptr = d;
    END_TIMING(s, 0, 0, 0);
    /* TODO check that everything is allocated; otherwise, call 'free' */

  return 0;
}
