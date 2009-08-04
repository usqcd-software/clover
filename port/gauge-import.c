#include <clover.h>

/* F memory order
 *  [0][1][0][0].re
 *  [0][1][0][0].im
 *  ...
 *  [0][1][0][C-1].im
 *  [0][1][1][0].re
 *  ...
 *  [0][1][C-1][C-1].im
 *  ... [0][2][0][0].re
 *  ... [0][3][0][0].re
 *  ... [1][2]...
 *  ... [1][3]...
 *  ... [2][3]...
 */
 
#if Q(PROJECTED_FERMION_DIM) == 2
#define Fxx(m,s,r,a,b)   ((m)[(r)+((a)*Q(COLORS)+(b))*2         \
                              +(s)*2*Q(COLORS)*Q(COLORS)])
#define F01_re(m,a,b)  Fxx(m, 0, 0, a, b)
#define F01_im(m,a,b)  Fxx(m, 0, 1, a, b)
#define F02_re(m,a,b)  Fxx(m, 1, 0, a, b)
#define F02_im(m,a,b)  Fxx(m, 1, 1, a, b)
#define F03_re(m,a,b)  Fxx(m, 2, 0, a, b)
#define F03_im(m,a,b)  Fxx(m, 2, 1, a, b)
#define F12_re(m,a,b)  Fxx(m, 3, 0, a, b)
#define F12_im(m,a,b)  Fxx(m, 3, 1, a, b)
#define F13_re(m,a,b)  Fxx(m, 4, 0, a, b)
#define F13_im(m,a,b)  Fxx(m, 4, 1, a, b)
#define F23_re(m,a,b)  Fxx(m, 5, 0, a, b)
#define F23_im(m,a,b)  Fxx(m, 5, 1, a, b)

#define X(m,a,i,b,j) ((m)                                               \
                      + 2 * ((j) + Q(PROJECTED_FERMION_DIM)             \
                             * ((b) + Q(COLORS)                         \
                                * ((i) + Q(PROJECTED_FERMION_DIM)       \
                                   * (a)))))
static void
qx(mk_cllo)(double v[], const double r[])
{
    int a, b;
    double f01r, f01i, f23r, f23i, f0123r, f0123i;
    double f12r, f12i, f03r, f03i, f1203r, f1203i;
    double f02r, f02i, f13r, f13i, f0213r, f0213i;
    double one;

    for (a = 0; a < Q(COLORS); a++) {
        for (b = 0; b < Q(COLORS); b++) {
            f01r = F01_re(r, a, b); f01i = F01_im(r, a, b);
            f02r = F02_re(r, a, b); f02i = F02_im(r, a, b);
            f03r = F03_re(r, a, b); f03i = F03_im(r, a, b);
            f12r = F12_re(r, a, b); f12i = F12_im(r, a, b);
            f13r = F13_re(r, a, b); f13i = F13_im(r, a, b);
            f23r = F23_re(r, a, b); f23i = F23_im(r, a, b);

            f0123r = f01r - f23r; f0123i = f01i - f23i;
            f1203r = f12r - f03r; f1203i = f12i - f03i;
            f0213r = f02r + f13r; f0213i = f02i + f13i;

            one = (a == b)? 1.0: 0.0;
            X(v, a, 0, b, 0)[0] = one - f0123r;
            X(v, a, 0, b, 0)[1] = - f0123i;
            X(v, a, 0, b, 1)[0] = -(f1203r + f0213i);
            X(v, a, 0, b, 1)[1] = -(f1203i - f0213r);
            X(v, a, 1, b, 0)[0] = -(f1203r - f0213i);
            X(v, a, 1, b, 0)[1] = -(f1203i + f0213r);
            X(v, a, 1, b, 1)[0] = one + f0123r;
            X(v, a, 1, b, 1)[1] = f0123i;
        }
    }
}

static void
qx(mk_clhi)(double v[], const double r[])
{
    int a, b;
    double f01r, f01i, f23r, f23i, f0123r, f0123i;
    double f12r, f12i, f03r, f03i, f1203r, f1203i;
    double f02r, f02i, f13r, f13i, f0213r, f0213i;
    double one;

    for (a = 0; a < Q(COLORS); a++) {
        for (b = 0; b < Q(COLORS); b++) {
            f01r = F01_re(r, a, b); f01i = F01_im(r, a, b);
            f02r = F02_re(r, a, b); f02i = F02_im(r, a, b);
            f03r = F03_re(r, a, b); f03i = F03_im(r, a, b);
            f12r = F12_re(r, a, b); f12i = F12_im(r, a, b);
            f13r = F13_re(r, a, b); f13i = F13_im(r, a, b);
            f23r = F23_re(r, a, b); f23i = F23_im(r, a, b);

            f0123r = f01r + f23r; f0123i = f01i + f23i;
            f1203r = f12r + f03r; f1203i = f12i + f03i;
            f0213r = f02r - f13r; f0213i = f02i - f13i;

            one = (a == b)? 1.0: 0.0;
            X(v, a, 0, b, 0)[0] = one - f0123r;
            X(v, a, 0, b, 0)[1] = - f0123i;
            X(v, a, 0, b, 1)[0] = -(f1203r + f0213i);
            X(v, a, 0, b, 1)[1] = -(f1203i - f0213r);
            X(v, a, 1, b, 0)[0] = -(f1203r - f0213i);
            X(v, a, 1, b, 0)[1] = -(f1203i + f0213r);
            X(v, a, 1, b, 1)[0] = one + f0123r;
            X(v, a, 1, b, 1)[1] = f0123i;
        }
    }
}

void qx(cl_invert)(double v[], const double r[]);
/* XXXX */

#endif /* Q(PROJECTED_FERMION_DIM) == 2 */

int
QX(import_gauge)(struct QX(Gauge) **gauge_ptr,
                 struct Q(State) *state,
                 double kappa,
                 double c_sw,
                 double (*U_reader)(int dir,
                                    const int pos[Q(DIM)],
                                    int a,
                                    int b,
                                    int re_im,
                                    void *env),
                 double (*F_reader)(int mu,
                                    int nu,
                                    const int pos[Q(DIM)],
                                    int a,
                                    int b,
                                    int re_im,
                                    void *env),
                 void *env)
{
  struct QX(Gauge) *gauge;
  void *ptr;
  size_t size;
  double r[Q(DIM) * Q(COLORS) * Q(COLORS) * 2];
  double cl[Q(DIM) * (Q(DIM) - 1) * Q(COLORS) * Q(COLORS)]; /* /2 *2 */
  double rr[Q(COLORS) * Q(COLORS)
            * Q(PROJECTED_FERMION_DIM) * Q(PROJECTED_FERMION_DIM) * 2];
  double rx[Q(COLORS) * Q(COLORS)
            * Q(PROJECTED_FERMION_DIM) * Q(PROJECTED_FERMION_DIM) * 2];
  double *v;
  double kc = kappa * c_sw;
  int p, d, a, b, mu, nu;
  int x[Q(DIM)];
  int u_s, ce_s, co_s;

  if (state == NULL || state->error_latched)
    return 1;

  if (gauge_ptr == NULL)
    return q(set_error)(state, 0, "import_gauge(): NULL pointer");

  *gauge_ptr = NULL;
  u_s = qx(sizeof_gauge)(state->volume);
  ce_s = qx(sizeof_clover)(state->even.full_size);
  co_s = qx(sizeof_clover)(state->odd.full_size);
  gauge = q(allocate_aligned)(state, &size, &ptr,
                              sizeof (struct QX(Gauge)), u_s + ce_s + 2 * co_s);
  if (gauge == 0)
    return q(set_error)(state, 0, "import_gauge(): not enough memory");

  BEGIN_TIMING(state);
  gauge->state = state;
  gauge->size = size;
  gauge->g_data = ptr;
  gauge->ce_data = (void *)(((char *)(gauge->g_data)) + u_s);
  gauge->co_data = (void *)(((char *)(gauge->ce_data)) + ce_s);
  gauge->cox_data = (void *)(((char *)(gauge->co_data)) + co_s);
  *gauge_ptr = gauge;

  for (p = 0; p < state->volume; p++) {
      q(l2v)(x, &state->local, state->lx2v[p]);
      for (v = r, d = 0; d < Q(DIM); d++) {
          for (a = 0; a < Q(COLORS); a++) {
              for (b = 0; b < Q(COLORS); b++) {
                  *v++ = kappa * U_reader(d, x, a, b, 0, env);
                  *v++ = kappa * U_reader(d, x, a, b, 1, env);
              }
          }
      }
      qx(put_gauge)(gauge->g_data, p, r);
  }

  for (p = 0; p < state->even.full_size; p++) {
      q(l2v)(x, state->even.local, state->even.lx2v[p]);
      for (v = cl, mu = 0; mu < Q(DIM); mu++) {
          for (nu = mu + 1; nu < Q(DIM); nu++) {
              for (a = 0; a < Q(COLORS); a++) {
                  for (b = 0; b < Q(COLORS); b++) {
                      *v++ = kc * F_reader(mu, nu, x, a, b, 0, env);
                      *v++ = kc * F_reader(mu, nu, x, a, b, 1, env);
                  }
              }
          }
      }
      qx(mk_cllo)(rr, cl);
      qx(put_clover_lo)(gauge->ce_data, p, rr);
      qx(mk_clhi)(rr, cl);
      qx(put_clover_hi)(gauge->ce_data, p, rr);
  }

  for (p = 0; p < state->odd.full_size; p++) {
      q(l2v)(x, state->odd.local, state->odd.lx2v[p]);
      for (v = cl, mu = 0; mu < Q(DIM); mu++) {
          for (nu = mu + 1; nu < Q(DIM); nu++) {
              for (a = 0; a < Q(COLORS); a++) {
                  for (b = 0; b < Q(COLORS); b++) {
                      *v++ = kc * F_reader(mu, nu, x, a, b, 0, env);
                      *v++ = kc * F_reader(mu, nu, x, a, b, 1, env);
                  }
              }
          }
      }
      qx(mk_cllo)(rr, cl);
      qx(put_clover_lo)(gauge->co_data, p, rr);
      qx(cl_invert)(rx, rr);
      qx(put_clover_lo)(gauge->cox_data, p, rx);
      qx(mk_clhi)(rr, cl);
      qx(put_clover_hi)(gauge->co_data, p, rr);
      qx(cl_invert)(rx, rr);
      qx(put_clover_hi)(gauge->cox_data, p, rx);
  }

  END_TIMING(state, 0, 0, 0);
  return 0;
}
