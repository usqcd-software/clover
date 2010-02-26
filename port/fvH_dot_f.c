#include <clover.h>
#include <string.h>

/*
 *  c[i] = herm(fv[fv_begin+i]) * g 
 *      for all i = (0 .. fv_len-1)
 *  c is complex vector as [re:0/im:1 + 2 * i]
 */
unsigned int
qx(fvH_dot_f)(int size,
              double *c,
              const struct vFermion *fv, int fv_size, int fv_begin, int fv_len,
              const struct Fermion *g)
{
    unsigned int res;

    memset(c, 0, 2 * fv_len * sizeof (double));
    res = qx(do_fvH_dot_f)(size, c, fv, fv_size, fv_begin, fv_len, g);
    QMP_sum_double_array(c, 2 * fv_len);

    return res;
}
