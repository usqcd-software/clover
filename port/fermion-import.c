#include <clover.h>

int
QX(import_fermion)(struct QX(Fermion) **fermion_ptr,
                   struct Q(State) *state,
                   double (*reader)(const int pos[4],
                                    int color,
                                    int dirac,
                                    int re_im,
                                    void *env),
                   void *env)
{
  double m[Q(FERMION_DIM) * Q(COLORS) * 2 * sizeof (double)];

  if (QX(allocate_fermion)(fermion_ptr, state) != 0)
    return 1;

  BEGIN_TIMING(state);
  qx(x_import)(&state->even, m, (*fermion_ptr)->even, reader, env);
  qx(x_import)(&state->odd, m, (*fermion_ptr)->odd, reader, env);
  END_TIMING(state, 0, 0, 0);

  return 0;
}
