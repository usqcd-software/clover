#include "clover-test.h"

/****
 **       t[1]
 **  ^-----------X
 **  |           |
 **  | nu        | t[0]
 **  |     mu    |
 **  o----------->
 **
 **
 **     C
 **   D   B
 **     A
 */

void
clover(QDP_ColorMatrix *cl[], QDP_ColorMatrix *U[])
{
    int mu, nu, i;
    QDP_ColorMatrix *t[5];

    create_Mvector(t, NELEMS(t));

    for (i = mu = 0; mu < NDIM; mu++) {
        for (nu = mu + 1; nu < NDIM; nu++, i++) {
            /* [0]: B */
            QDP_M_eq_sM(t[0], U[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
            /* [1]: C */
            QDP_M_eq_sM(t[1], U[mu], QDP_neighbor[nu], QDP_forward, QDP_all);
            /* [2]: BC */
            QDP_M_eq_M_times_Ma(t[2], t[0], t[1], QDP_all);
            /* [3]: BCD */
            QDP_M_eq_M_times_Ma(t[3], t[2], U[nu], QDP_all);
            /* [4]: ABCD */
            QDP_M_eq_M_times_M(t[4], U[mu], t[3], QDP_all); /* cl */
            /* [5]: BCDA */
            QDP_M_eq_M_times_M(t[5], t[3], U[mu], QDP_all);
            /* [6]: BCDA|mu */
            QDP_M_eq_sM(t[6], t[5], QDP_neighbor[mu], QDP_backward, QDP_all);
            /* [4]: ABCD+BCDA|mu */
            QDP_M_peq_M(t[4], t[6], QDP_all); /* cl */
            /* [3]: DA */
            QDP_M_eq_Ma_times_M(t[3], U[nu], U[mu], QDP_all);
            /* [5]: DABC */
            QDP_M_eq_M_times_M(t[5], t[3], t[2], QDP_all);
            /* [6]: DABC|nu */
            QDP_M_eq_sM(t[6], t[5], QDP_neighbor[nu], QDP_backward, QDP_all);
            /* [4]: ABCD+BCDA|mu+DABC|nu */
            QDP_M_peq_M(t[4], t[6], QDP_all); /* cl */
            /* [2]: CDA */
            QDP_M_eq_Ma_times_M(t[2], t[1], t[3], QDP_all);
            /* [5]: CDAB */
            QDP_M_eq_M_times_M(t[5], t[2], t[0], QDP_all);
            /* [6]: CDAB|mu */
            QDP_M_eq_sM(t[6], t[5], QDP_neighbor[mu], QDP_backward, QDP_all);
            /* [3]: CDAB|mu|nu */
            QDP_M_eq_sM(t[3], t[6], QDP_neighbor[nu], QDP_backward, QDP_all);
            /* [4]: ABCD+BCDA|mu+DABC|nu+CDAB|mu|nu */
            QDP_M_peq_M(t[4], t[6], QDP_all);
            /* cl[i]: [4]-[4]^+ */
            QDP_M_eq_M(cl[i], t[4], QDP_all);
            QDP_M_meq_conj_M(cl[i], t[4], QDP_all);
        }
    }

    destroy_Mvector(t, NELEMS(t));
    return;
}
