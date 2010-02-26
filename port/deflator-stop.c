#include <clover.h>

void
q(df_stop)(struct Q(Deflator) *d)
{
    if (NULL == d)
        return;
    d->frozen = 1;
}
