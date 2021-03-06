#include <stdio.h>
#include <stdlib.h>
#include <clover.h>
#include <math.h>
#include "vfermion-test.h"

static void
test_dot_vf(const char *name,
            int esize, 
            struct vFermion *v,
            int v_w, int v_b, int v_l,
            struct vFermion *w,
            int w_w, int w_b, int w_l,
            struct Fermion *t0,
            struct Fermion *t1,
            int ldc)
{
    int m_size = 2 * w_l * ldc * sizeof (double);
    int i, j;
    double *m_ab = malloc(m_size);
    double n_f, n_t;
    double ff_r, ff_i;
    double delta;
    double m_e = 0;

    memset(m_ab, 0, m_size);
    qx(do_fvH_dot_fv)(esize, m_ab, ldc,
                      v, v_w, v_b, v_l,
                      w, w_w, w_b, w_l);
    for (i = 0; i < v_l; i++) {
        n_f = 0;
        qx(fv_get)(esize, t0, v, v_w, i + v_b);
        qx(f_norm)(&n_f, esize, t0);
        for (j = 0; j < w_l; j++) {
            qx(fv_get)(esize, t1, w, w_w, j + w_b);
            ff_r = ff_i = 0;
            qx(f_dot)(&ff_r, &ff_i, esize, t0, t1);
            n_t = 0;
            qx(f_norm)(&n_t, esize, t1);
            ff_r -= m_ab[0 + 2 * (i + ldc * j)];
            ff_i -= m_ab[1 + 2 * (i + ldc * j)];
            delta = (ff_r * ff_r + ff_i * ff_i) / sqrt(n_f * n_t);
            printf("%s %3d %3d %15.7e\n", name, i, j, delta);
            if (delta > m_e)
                m_e = delta;
        }
    }
    printf("Max error: %s %15.7e\n", name, m_e);
}

int
main(int argc, char *argv[])
{
    int esize;
    int width;
    int v_begin, v_len;
    int w_begin, w_len;
    int ldc;
    struct vFermion *v0;
    struct vFermion *v1;
    struct Fermion *t0;
    struct Fermion *t1;

    if (argc != 8) {
        fprintf(stderr, "Usage: %s esize width v_b v_l w_b w_ln ldc\n",
                argv[0]);
        return 1;
    }
    esize = atoi(argv[1]);
    width = atoi(argv[2]);
    v_begin = atoi(argv[3]);
    v_len = atoi(argv[4]);
    w_begin = atoi(argv[5]);
    w_len = atoi(argv[6]);
    ldc = atoi(argv[7]);

    t0 = new_fermion(esize);
    t1 = new_fermion(esize);
    v0 = new_vfermion(esize, width);
    v1 = new_vfermion(esize, 2 * width);

    construct_vf(esize, width, v0, t0, 12.45);
    construct_vf(esize, 2 * width, v1, t0, -7.4565);

    test_dot_vf("test A", esize,
                v0, width, v_begin, v_len,
                v1, 2 * width, w_begin, w_len,
                t0, t1, ldc);

#if 0
    show_vfermion("V0", v0, esize, width, tmp);
    show_vfermion("V1", v1, esize, 2 * width, tmp);
#endif

    return 0;
}
