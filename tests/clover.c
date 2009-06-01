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
 **     C               3 0
 **   D   B             2 1
 **     A
 **
 ** clover indices:   mu nu   gamma           wilson mu
 **  0                 0  1     -i [3]         0          [1]
 **  1                 0  2     -i [5]         1          [2]
 **  2                 0  3     -i [9]         2          [4]
 **  3                 1  2     -i [6]         3          [8]
 **  4                 1  3     -i [10]
 **  5                 2  3     -i [12]
 */

void
clover(QDP_ColorMatrix *cl[], QDP_ColorMatrix *U[])
{
    int mu, nu, i;
    QDP_ColorMatrix *t[6];

    create_Mvector(t, NELEMS(t));

    for (i = mu = 0; mu < NDIM; mu++) {
        for (nu = mu + 1; nu < NDIM; nu++, i++) {
            /* [0]: B */
            QDP_M_eq_sM(t[0], U[nu], QDP_neighbor[mu], QDP_forward, QDP_all);
            /* [1]: DA */
            QDP_M_eq_Ma_times_M(t[1], U[nu], U[mu], QDP_all);
            /* [2]: DAB */
            QDP_M_eq_M_times_M(t[2], t[1], t[0], QDP_all);
            /* [3]: DAB|nu */
            QDP_M_eq_sM(t[3], t[2], QDP_neighbor[nu], QDP_backward, QDP_all);
            /* [4]: DABC|nu = p1 */
            QDP_M_eq_M_times_Ma(t[4], t[3], U[mu], QDP_all);
            /* [1]: C */
            QDP_M_eq_sM(t[1], U[mu], QDP_neighbor[nu], QDP_forward, QDP_all);
            /* [5]: CDAB|nu */
            QDP_M_eq_Ma_times_M(t[5], U[mu], t[3], QDP_all);
            /* [2]: BC */
            QDP_M_eq_M_times_Ma(t[2], t[0], t[1], QDP_all);
            /* [3]: BCD */
            QDP_M_eq_M_times_Ma(t[3], t[2], U[nu], QDP_all);
            /* [4]: ABCD + DABC|nu */
            QDP_M_peq_M_times_M(t[4], U[mu], t[3], QDP_all);
            /* [5]: BCDA + CDAB|nu */
            QDP_M_peq_M_times_M(t[5], t[3], U[mu], QDP_all);
            /* [2]: BCDA|mu + CDAB|nu|mu */
            QDP_M_eq_sM(t[2], t[5], QDP_neighbor[mu], QDP_backward, QDP_all);
            /* [4]: clover */
            QDP_M_peq_M(t[4], t[2], QDP_all);
            
            QDP_M_eq_M(cl[i], t[4], QDP_all);
            QDP_M_meq_conj_M(cl[i], t[4], QDP_all);
        }
    }

    destroy_Mvector(t, NELEMS(t));
    return;
}
