#include <clover.h>

#if QOP_CLOVER_DEFAULT_PRECISION == 'F'
#define DF_PREAMBLE(psi_e, chi_e) do { \
    if (NULL == deflator) qx(f_zero)(psi_e, e_size); \
    else { \
        if (q(df_preamble)(state, deflator, psi_e, chi_e)) { \
            q(set_error)(state, 0, "cg_solver() not enough memory"); \
            return 3; \
        } } } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)           \
    q(df_update0)(state, deflator, a1, b1, a0, b0, r, rho)
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho) \
    q(df_update1)(state, deflator, a1, b1, a0, b0, r, rho, A_rho)
#define DF_POSTAMBLE() \
    do { q(df_postamble)(state, deflator); } while (0)
#else
#define DF_PREAMBLE(psi_e, chi_e) do { \
        qx(f_zero)(psi_e, e_size);\
    } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)  0
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho)  0
#define DF_POSTAMBLE()  do {} while (0)
#endif

int
qx(cg_solver)(struct Fermion            *psi_e,
              const char                *name,
              int                       *out_iter,
              double                    *out_epsilon,
              struct Q(State)           *state,
              const struct QX(Gauge)    *gauge,
              const struct Fermion      *chi_e,
              struct Q(Deflator)        *deflator,
              int                        max_iter,
              double                     epsilon,
              unsigned                   options,
              long long                 *flops,
              long long                 *sent,
              long long                 *received,
              struct Fermion            *rho_e,
              struct Fermion            *pi_e,
              struct Fermion            *zeta_e,
              struct Fermion            *t0_e,
              struct Fermion            *t1_e,
              struct Fermion            *t0_o,
              struct Fermion            *t1_o)
{
    double a0 = 1, b0 = 0;
    int df_status;
    int e_size = state->even.full_size;
    double a, b, g, r, norm_omega;
    int i;

    DF_PREAMBLE(psi_e, (struct Fermion *) chi_e);
    qx(op_even_M)(t0_e, state, gauge, psi_e,
                  flops, sent, received,
                  t0_o);
    qx(op_even_Mx)(zeta_e, state, gauge, t0_e,
                   flops, sent, received,
                   t0_o);
    qx(f_add3)(rho_e, e_size, chi_e, -1, zeta_e);
    *flops += qx(f_norm)(&r, e_size, rho_e);
    QMP_sum_double(&r);
    qx(f_copy)(pi_e, e_size, rho_e);
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
            DF_POSTAMBLE();
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
        df_status = DF_UPDATE0(a, b, a0, b0, g, rho_e);
        if (-1 == df_status) {
            qx(op_even_M)(t0_e, state, gauge, rho_e,
                          flops, sent, received,
                          t0_o);
            qx(op_even_Mx)(zeta_e, state, gauge, t0_e,
                           flops, sent, received,
                           t0_o);
            df_status = DF_UPDATE1(a, b, a0, b0, g, rho_e, zeta_e);
        } 
        if (3 == df_status) {
            /* TODO: restart CG with new deflation */NOT_IMPLEMENTED;
        }
        a0 = a;
        b0 = b;
        if (options)
            qx(cg_log)(r, name,
                       i, psi_e, state, gauge, chi_e,
                       flops, sent, received,
                       options,
                       t0_e, t1_e, t0_o, t1_o);
    }
end:
    *out_iter = i;
    *out_epsilon = r;
    DF_POSTAMBLE();
    if (i == max_iter)
        return 1;
    return 0;
}
