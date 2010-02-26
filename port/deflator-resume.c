#include <clover.h>

void
q(df_resume)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    if (d->usize < d->umax)
        d->frozen = 0;
}
