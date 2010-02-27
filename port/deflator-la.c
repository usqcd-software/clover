#include <assert.h>
#define QOP_CLOVER_DEFAULT_PRECISION 'F'
#include <clover.h>
#include <qmp.h>


/* allocate & free */
latvec_c
q(latvec_c_alloc)(struct Q(State) *state, int dim)
{
    latvec_c res;
    res.dim     = dim;
    res.mem_ptr = qx(allocate_eo)(state, &res.mem_size,
                                  (void *)&res.f, 0, 1, 0);
    if (res.mem_ptr == NULL)
        res.f = NULL;
    return res;
}
latvec_c 
q(latvec_c_view)(int dim, struct FermionF *f)
{
    latvec_c res;
    res.dim     = dim;
    res.f       = f;
    res.mem_ptr = NULL;
    res.mem_size = 0;
    return res;
}
void 
q(latvec_c_copy)(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));
    qx(f_copy)(y.f, y.dim, x.f);
}
void
q(latvec_c_zero)(latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    qx(f_zero)(x.f, x.dim);
}
void
q(latvec_c_free)(struct Q(State) *state, latvec_c *v)
{
    if (NULL != v && NULL != v->mem_ptr && NULL != v->f) {
        q(free)(state, v->mem_ptr, v->mem_size);
        v->mem_ptr = NULL;
        v->f = NULL;
    }
}


doublecomplex 
q(lat_c_dotu)(latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    doublecomplex res = { 0., 0. };
    double s[2];

    qx(f_dot)(&s[0], &s[1], x.dim, x.f, y.f);
    QMP_sum_double_array(s, 2);
    res.r = s[0];
    res.i = s[1];
    
    return res;
}
void 
q(lat_c_scal_d)(double alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    qx(f_rmul1)(x.f, x.dim, alpha);
}

void 
q(lat_c_axpy_d)(double alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, alpha, x.f);
}
double 
q(lat_c_nrm2)(latvec_c x) 
{
    assert(!latvec_c_is_null(&x));
    
    double res = 0.;
    qx(f_norm)(&res, x.dim, x.f);
    QMP_sum_double(&res);
    
    return res;
}

#if 0
latvec_z
q(latvec_z_alloc)(struct Q(State) *state, int dim)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = NULL;
    /* TODO LATER allocate f */NOT_IMPLEMENTED;
    return res;
}
latvec_z
q(latvec_z_view)(int dim, struct FermionD *f)
{
    latvec_z res;
    res.dim     = dim;
    res.f       = f;
    return res;
}
void 
q(latvec_z_copy)(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void
q(latvec_z_free)(struct Q(State) *state, latvec_z *v)
{
    if (NULL != v && NULL != v->f) {
        /* TODO LATER free fermion */NOT_IMPLEMENTED;
        v->f = NULL;
    }
}

void 
q(latvec_zc_copy)(latvec_z x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_c_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}
void 
q(latvec_cz_copy)(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER copy fermion y <- x (note order!) */NOT_IMPLEMENTED;
}


doublecomplex 
q(lat_cz_dotu)(latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO LATER return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
doublecomplex 
q(lat_z_dotu)(latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));

    doublecomplex res = { 0., 0. };
    /* TODO LATER return x^H . y */NOT_IMPLEMENTED;
    
    return res;
}
/* norm */
double 
q(lat_z_nrm2)(latvec_z x) 
{
    assert(!latvec_z_is_null(&x));
    double res = 0.;
    
    /* TODO LATER check that all is correct */NOT_IMPLEMENTED;
    qop_d3_clover_f_norm(&res, x.dim, x.f);
    res = res * res;
    
    return res;
}


void 
q(lat_c_scal)(doublecomplex alpha, latvec_c x)
{
    assert(!latvec_c_is_null(&x));
    /* TODO LATER implement x <- alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_c_axpy)(doublecomplex alpha, latvec_c x, latvec_c y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_c_is_null(&y));

    qx(f_add2)(y.f, y.dim, alpha, x.f);
}
void 
q(lat_cz_axpy)(doublecomplex alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_z_axpy)(doublecomplex alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_cz_axpy_d)(double alpha, latvec_c x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_c_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}
void 
q(lat_z_axpy_d)(double alpha, latvec_z x, latvec_z y)
{
    assert(x.dim == y.dim);
    assert(!latvec_z_is_null(&x) &&
            !latvec_z_is_null(&y));
    /* TODO LATER implement y <- y + alpha * x */NOT_IMPLEMENTED;
}

#endif
    
    
latmat_c 
q(latmat_c_alloc)(struct Q(State) *state, int dim, int ncol)
{
    latmat_c res;
    res.dim     = dim;
    res.size    = ncol;
    res.begin   = 0;
    res.len     = ncol;
    res.fv      = NULL;
    res.mem_ptr = q(allocate_aligned)(state, &res.mem_size, (void *)&res.fv, 0,
                                      qx(sizeof_vfermion(dim, ncol)));
    if (res.mem_ptr == NULL)
        res.fv = NULL;

    return res;
}
latmat_c 
q(latmat_c_view)(int dim, int size, struct vFermion *fv)
{
    latmat_c res;
    res.dim     = dim;
    res.size    = size;
    res.begin   = 0;
    res.len     = size;
    res.fv       = fv;
    res.mem_ptr = NULL;
    res.mem_size = 0;
    return res;
}
void
q(latmat_c_free)(struct Q(State) *state, latmat_c *m) 
{
    if (NULL != m && NULL != m->mem_ptr && NULL != m->fv) {
        assert(0 == m->begin &&     
            m->size == m->len);
        /* at least some chance that this is not a sub-matrix */
        q(free)(state, m->mem_ptr, m->mem_size);
        m->mem_ptr = 0;
        m->fv = NULL;
    }
}
void 
q(latmat_c_copy)(latmat_c m1, latmat_c m2)
{
    assert(m1.len == m2.len && m1.dim == m2.dim);
    assert(!latmat_c_is_null(&m1) && 
            !latmat_c_is_null(&m2));
    
    qx(fv_copy)(m1.dim, m1.len, 
                m2.fv, m2.size, m2.begin,
                m1.fv, m1.size, m1.begin);
}
latmat_c
q(latmat_c_submat_col)(latmat_c m, int col, int ncol)
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
    res.mem_ptr = NULL;
    res.mem_size = 0;

    return res;
}
void
q(latmat_c_insert_col)(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len && v.dim == m.dim);
    assert(!latmat_c_is_null(&m) &&
            !latvec_c_is_null(&v));

    qx(fv_put)(m.dim,
               m.fv, m.size, col,
               v.f);
}
void
q(latmat_c_get_col)(latmat_c m, int col, latvec_c v)
{
    assert(col < m.len && v.dim == m.dim);
    assert(!latmat_c_is_null(&m) &&
            !latvec_c_is_null(&v));

    qx(fv_get)(m.dim,
               v.f,
               m.fv, m.size, col);
}


/* C <- A^\dag * B, A:lat*m, B:lat*n, C:m*n */
void 
q(lat_lmH_dot_lm)(int m, int n,
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
q(lat_lmH_dot_lv)(int m, 
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
q(lat_lm_dot_zm)(int n, int k, 
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
q(lat_lm_dot_zv)(int n, 
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

void
q(latvec_c_linop)(struct Q(State)         *s,
                  latvec_c                 y,
                  latvec_c                 x,
                  const struct Q(Gauge)   *g,
                  struct FermionF         *tmp_e,    
                  struct FermionF         *tmp_o)
{
    long long flops = 0;
    long long sent = 0;
    long long received = 0;

    /* XXX keep track of flops */
    qx(cg_operator)(s, y.f, g, x.f, tmp_e, tmp_o,
                    &flops, &sent, &received);
}
