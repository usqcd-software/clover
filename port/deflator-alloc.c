#include <clover.h>

int
Q(create_deflator)(struct Q(Deflator) **deflator_ptr,
                   struct Q(State) *state,
                   int nev, int vsize, double eps)
{
    struct Q(Deflator) *deflator;

    if (state == NULL || state->error_latched)
        return 1;

    if (deflator_ptr == NULL)
        return q(set_error)(state, 0, "create_deflator(): NULL pointer");
  
    *deflator_ptr = NULL;
    deflator = q(malloc)(state, sizeof (struct Q(Deflator)));
    if (deflator == 0)
        return q(set_error)(state, 0, "allocate_deflator(): not enough memory");

    /* XXX allocate other deflator components */

    BEGIN_TIMING(state);
    deflator->state = state;
    deflator->nev = nev;
    deflator->vsize = vsize;
    deflator->eps = eps;
    /* XXX initialize other deflator components */
    *deflator_ptr = deflator;
    END_TIMING(state, 0, 0, 0);

  return 0;
}
