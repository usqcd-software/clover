#include <clover.h>

#define XXX_DEBUG

#if QOP_CLOVER_DEFAULT_PRECISION == 'F'
#define DF_PREAMBLE(psi_e, rho_e, r, chi_e) do {                        \
        if (q(df_preamble)(state, deflator, psi_e, rho_e, r, chi_e,     \
                           &ws)) {                                      \
            q(set_error)(state, 0, "cg_solver() not enough memory");    \
            return CG_NOEMEM;                                           \
        } } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)                           \
    q(df_update0)(state, deflator, a1, b1, a0, b0, r, rho)
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho)                             \
    q(df_update1)(state, deflator, a1, b1, a0, b0, r, rho, A_rho)
#define DF_POSTAMBLE() \
    do { q(df_postamble)(state, deflator, &ws); } while (0)
#else
#define DF_PREAMBLE(psi_e, rho_e, r, chi_e) do {        \
        qx(f_zero)(psi_e, e_size);                      \
        qx(f_copy)(rho_e, e_size, chi_e);               \
        qx(f_norm)(r, e_size, rho_e);                   \
        QMP_sum_double(r);                              \
    } while (0)
#define DF_UPDATE0(a1,b1,a0,b0,r,rho)  0
#define DF_UPDATE1(a1,b1,a0,b0,r,rho,A_rho)  0
#define DF_POSTAMBLE()  do {} while (0)
#endif

void
qx(cg_operator)(struct Fermion           *res_e,
                const struct Fermion     *psi_e,
                struct MxM_workspace     *ws)
{
    qx(op_even_M)(ws->tmp_e, ws->state, ws->gauge, psi_e,
                  ws->flops, ws->sent, ws->received,
                  ws->tmp_o);
    qx(op_even_Mx)(res_e, ws->state, ws->gauge, ws->tmp_e,
                   ws->flops, ws->sent, ws->received,
                   ws->tmp_o);
}

CG_STATUS
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
    struct MxM_workspace  ws;

    ws.state = state;
    ws.gauge = gauge;
    ws.tmp_e = t0_e;
    ws.tmp_o = t0_o;
    ws.flops = flops;
    ws.sent = sent;
    ws.received = received;

    DF_PREAMBLE(psi_e, rho_e, &r, (struct Fermion *) chi_e);
#ifdef XXX_DEBUG
    {
        double r_chi, r_rho, r_psi;

        qx(f_norm)(&r_chi, e_size, chi_e);
        qx(f_norm)(&r_rho, e_size, rho_e);
        qx(f_norm)(&r_psi, e_size, psi_e);

        QMP_sum_double(&r_chi);
        QMP_sum_double(&r_rho);
        QMP_sum_double(&r_psi);

        printf("\nCG entry: chi %15.7e  rho %15.7e  psi %15.7e\n",
               r_chi, r_rho, r_psi);
        printf("  preamble r2 %15.7e\n", r);
    }
#endif
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
#ifdef XXX_DEBUG
            printf("Zmode exit\n");
#endif      
            return CG_ZEROMODE;
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
            qx(cg_operator)(zeta_e, rho_e, &ws);
#ifdef XXX_DEBUG
            {
                double zr[2];
                double zp[2];
                qx(f_dot)(&zr[0], &zr[1], e_size, zeta_e, rho_e);
                QMP_sum_double_array(zr, 2);
                qx(f_dot)(&zp[0], &zp[1], e_size, zeta_e, pi_e);
                QMP_sum_double_array(zp, 2);
                printf("cg %5d: update1  %15.7e %15.7e  %15.7e %15.7e\n",
                       i, zr[0], zr[1], zp[0], zp[1]);
            }
#endif
            df_status = DF_UPDATE1(a, b, a0, b0, g, rho_e, zeta_e);
        } 
        if (3 == df_status) {
            *out_iter = i;
            *out_epsilon = r;
            DF_POSTAMBLE();
#ifdef XXX_DEBUG
            printf("EC exit\n");
#endif
            return CG_EIGCONV;
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
    if (i == max_iter) {
#ifdef XXX_DEBUG
        printf("EC exit\n");
#endif
        return CG_MAXITER;
    }
#ifdef XXX_DEBUG
    printf("OK exit\n");
#endif
    return CG_SUCCESS;
}
