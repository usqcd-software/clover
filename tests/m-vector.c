#include "clover-test.h"

void
create_Mvector(QDP_ColorMatrix *U[], int size)
{
    int i;

    for (i = 0; i < size; i++)
        U[i] = QDP_create_M();
}

void
destroy_Mvector(QDP_ColorMatrix *U[], int size)
{
    int i;
    
    for (i = 0; i < size; i++) {
        QDP_destroy_M(U[i]); 
        U[i] = 0;
    }
}
