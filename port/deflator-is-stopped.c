#include <clover.h>

int
Q(deflator_is_stopped)(struct Q(Deflator) *d)
{
    if (d == NULL)
        return 1;
    return d->frozen;
}
