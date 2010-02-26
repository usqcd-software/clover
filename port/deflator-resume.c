#include <clover.h>

void
Q(deflator_resume)(struct Q(Deflator) *deflator_ptr)
{
    if (NULL == d)
        return;
    if (d->usize < d->umax)
        d->frozen = 0;
}
