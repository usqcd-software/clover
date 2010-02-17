#include <clover.h>

int
qx(cg_solver)(struct Fermion *psi_e, const char *source,
              int *out_iter,
              double *out_epsilon,
              struct Q(State) *state,
              const struct QX(Gauge) *gauge,
              const struct Fermion *chi_e,
              int max_iter,
              double epsilon,
              unsigned options,
              long long *flops,
              long long *sent,
              long long *received,
              struct Fermion *rho_e,
              struct Fermion *pi_e,
              struct Fermion *zeta_e,
              struct Fermion *t0_e,
              struct Fermion *t1_e,
              struct Fermion *t0_o,
              struct Fermion *t1_o)
{
    int e_size = state->even.full_size;
    double a, b, g, r, norm_omega;
    int i;

    qx(f_copy)(rho_e, e_size, chi_e);
    *flops += qx(f_norm)(&r, e_size, rho_e);
    QMP_sum_double(&r);
    qx(f_copy)(pi_e, e_size, rho_e);
    qx(f_zero)(psi_e, e_size);
    if (r < epsilon) {
        i = 0;
        goto end;
    }
    for (i = 0; i < max_iter; i++) {
        qx(op_even_Mn)(t0_e, &norm_omega, state, gauge, pi_e,
                       flops, sent, received,
                       t0_o);
        qx(op_even_Mx)(zeta_e, state, gauge, t0_e,
                       flops, sent, received,
                       t0_o);
        if (norm_omega == 0.0) {
            *out_iter = i;
            *out_epsilon = r;
            q(set_error)(state, 0, "cg_solver() hit zero mode");
            return 2;
        }
        a = r / norm_omega;
        *flops += qx(f_add2_norm)(rho_e, &g, e_size, -a, zeta_e);
        QMP_sum_double(&g);
        if (g < epsilon) {
            *flops += qx(f_add2)(psi_e, e_size, a, pi_e);
            r = g;
            break;
        }
        b = g / r;
        r = g;
        qx(cg_xp)(psi_e, pi_e, e_size, a, b, rho_e);
        if (options)
            qx(cg_log)(r, source,
                       i, psi_e, state, gauge, chi_e,
                       flops, sent, received,
                       options,
                       t0_e, t1_e, t0_o, t1_o);
    }
end:
    *out_iter = i;
    *out_epsilon = r;
    if (i == max_iter)
        return 1;
    return 0;
}
