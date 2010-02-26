#include <clover.h>

/*
 * c[i,j] = herm(fv[fv_begin + i]) . g[gv_begin+j] 
 *      for all i = (0 .. fv_len-1), 
 *              j = (0 .. gv_len-1),
 * c is a complex matrix as [re:0/im:1 + 2 * (i + ldc * j)]
 */
unsigned int
qx(fvH_dot_fv)(int size,
               double *c, int ldc,
               const struct vFermion *fv,
               int fv_size, int fv_begin, int fv_len,
               const struct vFermion *gv,
               int gv_size, int gv_begin, int gv_len)
{
    int j;
    unsigned int res;
    
    if (fv_len == ldc) {
        memset(c, 0, 2 * gv_len * fv_len * sizeof (double));
    } else {
        for (j = 0; j < gv_len; j++)
            memset(c + 2 * ldc * j, 0, 2 * fv_len * sizeof (double));
    }

    res = qx(do_fvH_dot_fv)(size, c, ldc,
                            fv, fv_size, fv_begin, fv_len,
                            gv, gv_size, gv_begin, gv_len);
    if (fv_len == ldc) {
        QMP_sum_double_array(c, 2 * gv_len * fv_len);
    } else {
        for (j = 0; j < gv_len; j++)
            QMP_sum_double_array(c + 2 * ldc * j, 2 * fv_len);
    }

    return res;
}
