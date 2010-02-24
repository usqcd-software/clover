#include <assert.h>
#include <clover.h>


/* allocate & free */
latvec_c
latvec_c_alloc(struct Q(State) *state, int dim)
{
    latvec_c res;
    res.dim     = dim;
    res.f       = NULL;
    /* TODO allocate fermion with all the icing needed */NOT_IMPLEMENTED;

    return res;
}
latvec_c 
latvec_c_view(int dim, struct FermionF *f)
{
    latvec_c res;
    res.dim     = dim;
    res.f       = f;
    return res;
}
void 
latvec_c_copy(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));
    /* TODO copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void
latvec_c_free(struct Q(State) *state, latvec_c *v)
{
    if (NULL != v && NULL != v->f) {
        /* TODO free fermion */NOT_IMPLEMENTED;
        v->f = NULL;
    }
}
latvec_z
latvec_z_alloc(struct Q(State) *state, int dim)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = NULL;
    /* TODO allocate f */NOT_IMPLEMENTED;
    return res;
}
latvec_z
latvec_z_view(int dim, struct FermionD *f)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = f;
    return res;
}
void 
latvec_z_copy(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void
latvec_z_free(struct Q(State) *state, latvec_z *v)
{
    if (NULL != v && NULL != v->f) {
        /* free fermion */NOT_IMPLEMENTED;
        v->f = NULL;
    }
}

void 
latvec_zc_copy(latvec_z x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_c_is_null(&y));
    /* TODO copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void 
latvec_cz_copy(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}


    
    
latmat_c 
latmat_c_alloc(struct Q(State) *state, int dim, int ncol)
{
    latmat_c res;
    res.dim     = dim;
    res.size    = ncol;
    res.begin   = 0;
    res.len     = ncol;
    res.fv      = NULL;
    /* TODO allocate fermion "matrix" */NOT_IMPLEMENTED;
    return res;
}
latmat_c 
latmat_c_view(int dim, int size, struct vFermion *fv)
{
    latmat_c res;
    res.dim     = dim;
    res.size    = size;
    res.begin   = 0;
    res.len     = size;
    res.fv       = fv;
    return res;
}
void
latmat_c_free(struct Q(State) *state, latmat_c *m) 
{
    if (NULL != m && NULL != m->fv) {
        assert(0 == m->begin &&     
            m->size == m->len); /* at least some chance that this is not a sub-matrix */
        /* TODO free matrix */NOT_IMPLEMENTED;
        m->fv = NULL;
    }
}
void 
latmat_c_copy(latmat_c m1, latmat_c m2)
{
    assert(m1.len == m2.len && m1.dim == m2.dim);
    assert(!latmat_c_is_null(&m1) && 
            !latmat_c_is_null(&m2));
    
    qx(fv_copy)(m1.dim, m1.len, 
                m2.fv, m2.size, m2.begin,
                m1.fv, m1.size, m1.begin);
}
latmat_c
latmat_c_submat_col(latmat_c m, int col, int ncol)
{
    assert(0 <= col && 
            0 < ncol &&
            col < m.len && 
            col + ncol <= m.len);
    assert(!latmat_c_is_null(&(m)));
    latmat_c res;
    res.dim     = m.dim;
    res.size    = m.size;
    res.begin   = m.begin + col;
    res.len     = ncol;
    res.fv      = m.fv;

    return res;
}
void
latmat_c_insert_col(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len && v.dim == m.dim);
    assert(!latmat_c_is_null(&m) &&
            !latvec_c_is_null(&v));

    qx(fv_put)(m.dim,
               m.fv, m.size, col,
               v.f);
}
void
latmat_c_get_col(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len && v.dim == m.dim);
    assert(!latmat_c_is_null(&m) &&
            !latvec_c_is_null(&v));

    qx(fv_get)(m.dim,
               v.f,
               m.fv, m.size, col);
}

/* linear algebra */
/* dot products */
doublecomplex 
lat_cc_dotu(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
doublecomplex 
lat_cz_dotu(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
doublecomplex 
lat_zz_dotu(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
/* norm */
double 
lat_c_nrm2(latvec_c x) 
{
    assert(!latvec_c_is_null(&x));
    
    double res = 0.;
    /* TODO check that all is correct */NOT_IMPLEMENTED;
    qop_f3_clover_f_norm(&res, x.dim, x.f);
    res = res * res;
    
    return res;
}
double 
lat_z_nrm2(latvec_z x) 
{
    assert(!latvec_z_is_null(&x));
    double res = 0.;
    
    /* TODO check that all is correct */NOT_IMPLEMENTED;
    qop_d3_clover_f_norm(&res, x.dim, x.f);
    res = res * res;
    
    return res;
}

/* vector ops */
#if 0
void 
lat_c_scal(doublecomplex alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    /* TODO implement x <- alpha * x */NOT_IMPLEMENTED;
}
void 
lat_cc_axpy(doublecomplex alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, alpha, x.f);
}
void 
lat_cz_axpy(doublecomplex alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
lat_zz_axpy(doublecomplex alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
#endif
void 
lat_c_scal_d(double alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    /* TODO implement x <- alpha * x */NOT_IMPLEMENTED;
}

void 
lat_cc_axpy_d(double alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, alpha, x.f);
}
void 
lat_cz_axpy_d(double alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
lat_zz_axpy_d(double alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO implement y <- y + alpha * x */NOT_IMPLEMENTED;
}


/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
void 
lat_lmH_dot_lm(int m, int n,
               latmat_c a, 
               latmat_c b,
               doublecomplex *c, int ldc)
{
    assert(a.dim == b.dim &&
            a.len == m &&
            b.len == n &&
            m <= ldc);
    assert(!latmat_c_is_null(&a) &&
            !latmat_c_is_null(&b));
    
    qx(fvH_dot_fv)(a.dim,
                   (double *)c, ldc,
                   a.fv, a.size, a.begin, a.len,
                   b.fv, b.size, b.begin, b.len);

}
/* y <- A^\dag x, A:lat*m, x:lat, y:m */
void 
lat_lmH_dot_lv(int m, 
               latmat_c a, 
               latvec_c x,  
               doublecomplex *y)
{
    assert(a.dim == x.dim &&
            a.len == m);
    assert(!latmat_c_is_null(&a) &&
            !latvec_c_is_null(&x));

    qx(fvH_dot_f)(a.dim,
                  (double *)y, 
                  a.fv, a.size, a.begin, a.len,
                  x.f);
}
/* C <- A * B, A:lat*k, B:k*n, C:lat*n */
void 
lat_lm_dot_zm(int n, int k, 
              latmat_c a, 
              doublecomplex *b, int ldb, 
              latmat_c c)
{
    assert(a.dim == c.dim &&
            a.len == k &&
            c.len == n &&
            k <= ldb);
    assert(!latmat_c_is_null(&a) && 
            !latmat_c_is_null(&c));

    qx(fv_dot_zm)(a.dim,
                  c.fv, c.size, c.begin, c.len,
                  a.fv, a.size, a.begin, a.len,
                  (double *)b, ldb);
}
/* y <- A * x, A:lat*n, x:n, y:lat */
void 
lat_lm_dot_zv(int n, 
              latmat_c a, 
              doublecomplex *x,
              latvec_c y)
{
    assert(a.dim == y.dim &&
            a.len == n);
    assert(!latmat_c_is_null(&a) &&
            !latvec_c_is_null(&y));

    qx(fv_dot_zv)(a.dim,
                  y.f,
                  a.fv, a.size, a.begin, a.len,
                  (double *)x);
}
