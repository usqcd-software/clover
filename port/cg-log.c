#include <clover.h>

static double
qx(cg_true_residual)(const struct Fermion *xi_e,
                     struct Q(State) *state,
                     const struct QX(Gauge) *gauge,
                     const struct Fermion *chi_e,
                     long long *flops,
                     long long *sent,
                     long long *received,
                     struct Fermion *t0_e,
                     struct Fermion *t1_e,
                     struct Fermion *t0_o)
{
    double norm;

    qx(op_even_M)(t0_e, state, gauge, xi_e, flops, sent, received,
                  t0_o);
    qx(op_even_Mx)(t1_e, state, gauge, t0_e, flops, sent, received,
                   t0_o);
    *flops += qx(f_diff_norm)(&norm, state->even.full_size, t1_e, chi_e);
    QMP_sum_double(&norm);

    return norm;
}

void
qx(cg_log)(double cg_res, const char *source, int iter,
           const struct Fermion *xi_e,
           struct Q(State) *state,
           const struct QX(Gauge) *gauge,
           const struct Fermion *chi_e,
           long long *flops,
           long long *sent,
           long long *received,
           unsigned int options,
           struct Fermion *t0_e,
           struct Fermion *t1_e,
           struct Fermion *t0_o,
           struct Fermion *t1_o)
{
    double true_res = 0.0;

    if (options & Q(LOG_TRUE_RESIDUAL)) {
        true_res = qx(cg_true_residual)(xi_e, state, gauge, chi_e,
                                        flops, sent, received,
                                        t0_e, t1_e, t0_o);
    }
#define ITER_LOG (Q(LOG_CG_RESIDUAL) |    \
                  Q(LOG_TRUE_RESIDUAL))
    switch (options & ITER_LOG) {
#undef ITER_LOG
    default:
        break;
    case Q(LOG_CG_RESIDUAL):
        qx(zprint)(state, source,
                   "CG step %5d"
                   "  CG residual %11.4e",
                   iter, cg_res);
        break;
    case Q(LOG_TRUE_RESIDUAL):
        qx(zprint)(state, source,
                   "CG step %5d"
                   "  true residual %11.4e",
                   iter, true_res);
        break;
    case Q(LOG_CG_RESIDUAL) | Q(LOG_TRUE_RESIDUAL):
        qx(zprint)(state, source,
                   "CG step %5d"
                   "  CG residual %11.4e"
                   "  true residual %11.4e",
                   iter, cg_res, true_res);
        break;
    }
}
