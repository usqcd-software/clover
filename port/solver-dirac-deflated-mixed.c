#define QOP_CLOVER_DEFAULT_PRECISION 'F'
#include <clover.h>
#undef QOP_CLOVER_DEFAULT_PRECISION
#define QOP_CLOVER_DEFAULT_PRECISION 'D'
#include <clover.h>

/* Solve
 *   D_cl psi = eta
 *
 * with psi_0 as an initial guess
 *
 */ 

int
Q(deflated_mixed_D_CG)(struct QD(Fermion)          *psi,
                       int                         *out_iterations,
                       double                      *out_epsilon,
                       const struct QD(Fermion)    *psi_0,
                       const struct QD(Gauge)      *gauge,
                       const struct QD(Fermion)    *eta,
                       struct Q(Deflator)          *deflator,
                       int                          f_iter,
                       double                       f_epsilon,
                       int                          max_iterations,
                       double                       min_epsilon,
                       unsigned int                 options)
{
    DECLARE_STATE;

    /* check arguments */
    CHECK_ARG0(psi);
    CHECK_ARGn(psi_0, "eigD_CG");
    CHECK_ARGn(gauge, "eigD_CG");
    CHECK_ARGn(eta, "eigD_CG");
    CHECK_ARGn(deflator, "eigD_CG");

    return q(mixed_cg)(state, "eigD_CG",
                       psi, out_iterations, out_epsilon,
                       psi_0, gauge, eta, deflator,
                       f_iter, f_epsilon, max_iterations, min_epsilon,
                       options);
}
