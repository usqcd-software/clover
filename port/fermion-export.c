#include <clover.h>

int
QX(export_fermion)(void (*writer)(const int pos[4],
                                  int color,
                                  int dirac,
                                  int re_im,
                                  double value,
                                  void *env),
                   void *env,
                   const struct QX(Fermion) *fermion)
{
  struct Q(State) *state;
  double m[Q(FERMION_DIM) * Q(COLORS) * 2 * sizeof (double)];

  if (fermion == 0)
    return 1;

  state = fermion->state;

  BEGIN_TIMING(state);
  qx(x_export)(&state->even, m, fermion->even, writer, env);
  qx(x_export)(&state->odd, m, fermion->odd, writer, env);
  END_TIMING(state, 0, 0, 0);

  return 0;
}
