#include <clover.h>

void
Q(deflator_stop)(struct Q(Deflator) *deflator_ptr)
{
    if (NULL == d)
        return;
    d->frozen = 1;
}
