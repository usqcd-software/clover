#include <clover.h>

void
Q(free_deflator)(struct Q(Deflator) **deflator_ptr)
{
    struct Q(State) *state;

    if (deflator_ptr == 0 || *deflator_ptr == 0)
        return;
    state = (*deflator_ptr)->state;
    BEGIN_TIMING(state);
    /* XXX free other components of the deflator */
    q(free)(state, *deflator_ptr, sizeof (struct Q(Deflator)));
    END_TIMING(state, 0, 0, 0);
    *deflator_ptr = 0;
}
